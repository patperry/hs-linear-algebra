{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Internal
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Internal (
    -- * Banded matrix data types
    BMatrix(..),
    Banded,
    IOBanded,

    module BLAS.Matrix.Base,
    module BLAS.Tensor,

    -- * Converting to and from foreign pointers
    toForeignPtr,
    fromForeignPtr,
    
    -- * To and from the underlying storage matrix
    toRawMatrix,
    fromRawMatrix,
    
    -- * Bandwith properties
    bandwidth,
    numLower,
    numUpper,
    
    -- * Creating new matrices
    -- ** Pure
    banded,
    listsBanded,
    -- ** Impure
    newBanded_,
    newBanded,
    newListsBanded,
    
    -- * Getting rows and columns
    row,
    col,
    rows,
    cols,
    getRow,
    getCol,
    toLists,
    
    -- * Vector views
    diag,
    rowView,
    colView,
    
    -- * Casting matrices
    coerceBanded,
    
    -- * Unsafe operations
    unsafeBanded,
    unsafeNewBanded,
    unsafeFreeze,
    unsafeThaw,
    unsafeWithElemPtr,
    unsafeWithBasePtr,
    unsafeDiag,
    unsafeGetRow,
    unsafeGetCol,
    unsafeRow,
    unsafeCol,
    unsafeRowView,
    unsafeColView,
    ) where

import Control.Arrow ( second )
import Control.Monad ( zipWithM_ )
import Data.Ix ( inRange, range )
import Data.List ( foldl' )
import Data.Maybe ( fromJust )
import Foreign
import System.IO.Unsafe
import Unsafe.Coerce          

import BLAS.Access
import BLAS.Elem ( Elem, BLAS1 )
import qualified BLAS.Elem as E
import BLAS.Internal ( checkedRow, checkedCol, checkedDiag, diagStart, 
    diagLen, clearArray, inlinePerformIO )
                  
import BLAS.Matrix.Base hiding ( Matrix )
import qualified BLAS.Matrix.Base as C
import BLAS.Tensor

import Data.AEq

import Data.Matrix.Dense.Internal ( DMatrix )
import qualified Data.Matrix.Dense.Internal as M
                                
import Data.Vector.Dense.Internal ( DVector, Vector, conj, dim, newListVector )
import qualified Data.Vector.Dense.Internal as V

        
data BMatrix t mn e 
    = BM { fptrOf   :: {-# UNPACK #-} !(ForeignPtr e)
         , offsetOf :: {-# UNPACK #-} !Int
         , size1    :: {-# UNPACK #-} !Int
         , size2    :: {-# UNPACK #-} !Int
         , lowBW    :: {-# UNPACK #-} !Int
         , upBW     :: {-# UNPACK #-} !Int
         , ldaOf    :: {-# UNPACK #-} !Int
         , isHerm   :: {-# UNPACK #-} !Bool
         }

type Banded = BMatrix Imm
type IOBanded = BMatrix Mut


toLists :: (BLAS1 e) => Banded (m,n) e -> ((Int,Int), (Int,Int),[[e]])
toLists a = ( (m,n)
            , (kl,ku)
            , map paddedDiag [(-kl)..ku]
            )
  where
    (m,n)   = shape a
    (kl,ku) = (numLower a, numUpper a)
    
    padBegin i = replicate (max (-i) 0)    0
    padEnd   i = replicate (max (m-n+i) 0) 0
    paddedDiag i = (padBegin i) ++ (elems $ diag a i) ++ (padEnd i)

                



banded :: (BLAS1 e) => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> Banded (m,n) e
banded mn kl ijes = unsafePerformIO $ newBanded mn kl ijes
{-# NOINLINE banded #-}

unsafeBanded :: (BLAS1 e) => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> Banded (m,n) e
unsafeBanded mn kl ijes = unsafePerformIO $ unsafeNewBanded mn kl ijes
{-# NOINLINE unsafeBanded #-}

newBanded :: (BLAS1 e) => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> IO (BMatrix t (m,n) e)
newBanded = newBandedHelp writeElem

unsafeNewBanded :: (BLAS1 e) => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> IO (BMatrix t (m,n) e)
unsafeNewBanded = newBandedHelp unsafeWriteElem

newBandedHelp :: (BLAS1 e) => 
       (IOBanded (m,n) e -> (Int,Int) -> e -> IO ()) 
    -> (Int,Int) -> (Int,Int) -> [((Int,Int),e)] -> IO (BMatrix t (m,n) e)
newBandedHelp set (m,n) (kl,ku) ijes = do
    x <- newBanded_ (m,n) (kl,ku)
    withForeignPtr (fptrOf x) $ flip clearArray ((kl+1+ku)*n)
    mapM_ (uncurry $ set $ unsafeThaw x) ijes
    return x

listsBanded :: (BLAS1 e) => (Int,Int) -> (Int,Int) -> [[e]] -> Banded (m,n) e
listsBanded mn kl xs = unsafePerformIO $ newListsBanded mn kl xs
{-# NOINLINE listsBanded #-}

newListsBanded :: (BLAS1 e) => (Int,Int) -> (Int,Int) -> [[e]] -> IO (BMatrix t (m,n) e)
newListsBanded (m,n) (kl,ku) xs = do
    a <- newBanded_ (m,n) (kl,ku)
    zipWithM_ (writeDiagElems (unsafeThaw a)) [(negate kl)..ku] xs
    return a
  where
    writeDiagElems a i es =
        let d   = diag a i
            nb  = max 0 (negate i)
            es' = drop nb es
        in zipWithM_ (unsafeWriteElem d) [0..(dim d - 1)] es'

unsafeDiag :: (Elem e) => BMatrix t (m,n) e -> Int -> DVector t k e
unsafeDiag a d
    | isHerm a = conj $ unsafeDiag (herm a) (negate d)
    | otherwise =
        let f      = fptrOf a
            off    = indexOf a (diagStart d)
            len    = diagLen (shape a) d
            stride = ldaOf a
            c      = False
        in V.fromForeignPtr f off len stride c
        
diag :: (Elem e) => BMatrix t (m,n) e -> Int -> DVector t k e
diag a = checkedDiag (shape a) (unsafeDiag a) 



              
            
unsafeWithElemPtr :: (Elem e) => BMatrix t (m,n) e -> (Int,Int) -> (Ptr e -> IO a) -> IO a
unsafeWithElemPtr a (i,j) f
    | isHerm a  = unsafeWithElemPtr (herm a) (j,i) f
    | otherwise = withForeignPtr (fptrOf a) $ \ptr ->
                      f $ ptr `advancePtr` (indexOf a (i,j))

unsafeWithBasePtr :: (Elem e) => BMatrix t (m,n) e -> (Ptr e -> IO a) -> IO a
unsafeWithBasePtr a f =
    withForeignPtr (fptrOf a) $ \ptr ->
        f $ ptr `advancePtr` (offsetOf a)

row :: (BLAS1 e) => Banded (m,n) e -> Int -> Vector n e
row a = checkedRow (shape a) (unsafeRow a)

unsafeRow :: (BLAS1 e) => Banded (m,n) e -> Int -> Vector n e
unsafeRow a i = unsafePerformIO $ getRow a i
{-# NOINLINE unsafeRow #-}

getRow :: (BLAS1 e) => BMatrix t (m,n) e -> Int -> IO (DVector r n e)
getRow a = checkedRow (shape a) (unsafeGetRow a)

unsafeGetRow :: (BLAS1 e) => BMatrix t (m,n) e -> Int -> IO (DVector r n e)
unsafeGetRow a i = 
    let (nb,x,na) = unsafeRowView a i
        n = numCols a
    in do
        es <- getElems x
        newListVector n $ (replicate nb 0) ++ es ++ (replicate na 0)

col :: (BLAS1 e) => Banded (m,n) e -> Int -> Vector m e
col a = checkedCol (shape a) (unsafeCol a)

unsafeCol :: (BLAS1 e) => Banded (m,n) e -> Int -> Vector m e
unsafeCol a i = unsafePerformIO $ getCol a i
{-# NOINLINE unsafeCol #-}

getCol :: (BLAS1 e) => BMatrix t (m,n) e -> Int -> IO (DVector r m e)
getCol a = checkedCol (shape a) (unsafeGetCol a)

unsafeGetCol :: (BLAS1 e) => BMatrix t (m,n) e -> Int -> IO (DVector r m e)
unsafeGetCol a j = unsafeGetRow (herm a) j >>= return . conj

rows :: (BLAS1 e) => Banded (m,n) e -> [Vector n e]
rows a = [ unsafeRow a i | i <- [0..(numRows a - 1)] ]

cols :: (BLAS1 e) => Banded (m,n) e -> [Vector m e]
cols a = [ unsafeCol a i | i <- [0..(numCols a - 1)] ]

unsafeColView :: (Elem e) => BMatrix t (m,n) e -> Int -> (Int, DVector t k e, Int)
unsafeColView a@(BM f off m _ kl ku ld _) j
    | isHerm a = 
        case unsafeRowView (herm a) j of (nb, v, na) -> (nb, conj v, na)
    | otherwise =
        let nb     = max (j - ku)         0
            na     = max (m - 1 - j - kl) 0
            r      = max (ku - j) 0 
            c      = j 
            off'   = off + r + c * ld
            stride = 1
            len    = m - (nb + na)
        in if len >= 0
            then (nb, V.fromForeignPtr f off' len stride False, na)
            else (m , V.fromForeignPtr f off' 0   stride False,  0)


unsafeRowView :: (Elem e) => BMatrix t (m,n) e -> Int -> (Int, DVector t k e, Int)
unsafeRowView a@(BM f off _ n kl ku ld _) i
    | isHerm a =
        case unsafeColView (herm a) i of (nb, v, na) -> (nb, conj v, na)        
    | otherwise =
        let nb     = max (i - kl)         0
            na     = max (n - 1 - i - ku) 0
            r      = min (ku + i)         (kl + ku)
            c      = max (i - kl)         0 
            off'   = off + r + c * ld
            stride = ld - 1
            len    = n - (nb + na)
        in if len >= 0 
            then (nb, V.fromForeignPtr f off' len stride False, na)
            else (n , V.fromForeignPtr f off' 0   stride False,  0)


rowView :: (Elem e) => BMatrix t (m,n) e -> Int -> (Int, DVector t k e, Int)
rowView a = checkedRow (shape a) (unsafeRowView a) 

colView :: (Elem e) => BMatrix t (m,n) e -> Int -> (Int, DVector t k e, Int)
colView a = checkedCol (shape a) (unsafeColView a)


instance C.Matrix (BMatrix t) where
    numRows a | isHerm a  = size2 a
              | otherwise = size1 a
              
    numCols a | isHerm a  = size1 a
              | otherwise = size2 a

    herm a = let h' = (not . isHerm) a
             in coerceBanded $ a{ isHerm=h' }

    
instance (BLAS1 e) => ITensor (BMatrix Imm (m,n)) (Int,Int) e where
    size = inlinePerformIO . getSize
    
    unsafeAt a = inlinePerformIO . (unsafeReadElem a)
    
    indices = inlinePerformIO . getIndices
    elems   = inlinePerformIO . getElems
    assocs  = inlinePerformIO . getAssocs

    (//)          = replaceHelp writeElem
    unsafeReplace = replaceHelp unsafeWriteElem

    amap f a = banded (shape a) (numLower a, numUpper a) ies
      where
        ies = map (second f) (assocs a)

replaceHelp :: (BLAS1 e) => 
       (IOBanded (m,n) e -> (Int,Int) -> e -> IO ())
    -> Banded (m,n) e -> [((Int,Int), e)] -> Banded (m,n) e
replaceHelp set x ies =
    unsafeFreeze $ unsafePerformIO $ do
        y  <- newCopy (unsafeThaw x)
        mapM_ (uncurry $ set y) ies
        return y
{-# NOINLINE replaceHelp #-}
    
    
    
instance (BLAS1 e) => RTensor (BMatrix t (m,n)) (Int,Int) e IO where
    newCopy b = 
        let (mn,kl,a,h) = toRawMatrix b
        in do
            a' <- newCopy a
            return $ fromJust $ fromRawMatrix mn kl a' h
    


instance (BLAS1 e) => Show (BMatrix Imm (m,n) e) where
    show a 
        | isHerm a = 
           "herm (" ++ show (herm a) ++ ")"
        | otherwise = 
             let (mn,kl,es) = toLists a 
             in "listsBanded " ++ show mn ++ " " ++ show kl ++ " " ++ show es
       
       
compareHelp :: (BLAS1 e) => 
    (e -> e -> Bool) -> Banded (m,n) e -> Banded (m,n) e -> Bool
compareHelp cmp x y
    | isHerm x && isHerm y =
        compareHelp cmp (herm x) (herm y)
compareHelp cmp x y =
    (shape x == shape y) && (and $ zipWith cmp (elems x) (elems y))

instance (BLAS1 e, Eq e) => Eq (BMatrix Imm (m,n) e) where
    (==) = compareHelp (==)

instance (BLAS1 e, AEq e) => AEq (BMatrix Imm (m,n) e) where
    (===) = compareHelp (===)
    (~==) = compareHelp (~==)       
             