{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
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
    = BM { fptrOf   :: !(ForeignPtr e)
         , offsetOf :: !Int
         , size1    :: !Int
         , size2    :: !Int
         , lowBW    :: !Int
         , upBW     :: !Int
         , ldaOf    :: !Int
         , isHerm   :: !Bool
         }

type Banded = BMatrix Imm
type IOBanded = BMatrix Mut

fromForeignPtr :: ForeignPtr e -> Int -> (Int,Int) -> (Int,Int) -> Int -> Bool
    -> BMatrix t (m,n) e
fromForeignPtr f o (m,n) (kl,ku) l h = BM f o m n kl ku l h

toForeignPtr :: BMatrix t (m,n) e -> (ForeignPtr e, Int, (Int,Int), (Int,Int), Int, Bool)
toForeignPtr (BM f o m n kl ku l h) = (f, o, (m,n), (kl,ku), l, h)

unsafeFreeze :: BMatrix t mn e -> Banded mn e
unsafeFreeze = unsafeCoerce

unsafeThaw :: BMatrix t mn e -> IOBanded mn e
unsafeThaw = unsafeCoerce

-- | Coerce the phantom shape type from one type to another.
coerceBanded :: BMatrix t mn e -> BMatrix t kl e
coerceBanded = unsafeCoerce

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

toRawMatrix :: (Elem e) => BMatrix t (m,n) e -> ((Int,Int), (Int,Int), DMatrix t (m',n') e, Bool)
toRawMatrix (BM f o m n kl ku ld h) = 
    ((m,n), (kl,ku), M.fromForeignPtr f o (kl+1+ku,n) ld, h)

fromRawMatrix :: (Elem e) => (Int,Int) -> (Int,Int) -> DMatrix t (m,n) e -> Bool -> Maybe (BMatrix t (m',n') e)
fromRawMatrix (m,n) (kl,ku) a h = 
    if M.isHerm a 
        then Nothing
        else let (f,o,(m',n'),ld) = M.toForeignPtr a
             in case undefined of
                 _ | m' /= kl+1+ku -> 
                     error $ "fromMatrix: number of rows must be equal to number of diagonals"
                 _ | n' /= n ->
                     error $ "fromMatrix: numbers of columns must be equal"
                 _ ->
                     Just $ BM f o m n kl ku ld h
                

bandwidth :: BMatrix t (m,n) e -> (Int,Int)
bandwidth a = 
    let (kl,ku) = (numLower a, numUpper a)
    in (negate kl, ku)

numLower :: BMatrix t (m,n) e -> Int
numLower a | isHerm a  = upBW a
           | otherwise = lowBW a

numUpper :: BMatrix t (m,n) e -> Int
numUpper a | isHerm a  = lowBW a
           | otherwise = upBW a


newBanded_ :: (Elem e) => (Int,Int) -> (Int,Int) -> IO (BMatrix t (m,n) e)
newBanded_ (m,n) (kl,ku)
    | m < 0 || n < 0 =
        err "dimensions must be non-negative."
    | kl < 0 =
        err "lower bandwdth must be non-negative."
    | m /= 0 && kl >= m =
        err "lower bandwidth must be less than m."
    | ku < 0 =
        err "upper bandwidth must be non-negative."
    | n /= 0 && ku >= n =
        err "upper bandwidth must be less than n."
    | otherwise =
        let off = 0
            m'  = kl + 1 + ku
            l   = m'
            h   = False
        in do    
            ptr <- mallocForeignPtrArray (m' * n)
            return $ fromForeignPtr ptr off (m,n) (kl,ku) l h
    where
      err s = ioError $ userError $ 
                  "newBanded_ " ++ show (m,n) ++ " " ++ show (kl,ku) ++ ": " ++ s

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


indexOf :: BMatrix t (m,n) e -> (Int,Int) -> Int
indexOf (BM _ off _ _ _ ku ld h) (i,j) =
    let (i',j') = if h then (j,i) else (i,j)
    in off + ku + (i' - j') + j' * ld
    --off + i' * tda + (j' - i' + kl)
           

hasStorage :: BMatrix t (m,n) e -> (Int,Int) -> Bool
hasStorage (BM _ _ m n kl ku _ h) (i,j) =
    let (i',j') = if h then (j,i) else (i,j)
    in (  inRange (0,m-1) i'
       && inRange (0,n-1) j'
       && inRange (max 0 (j'-ku), min (m-1) (j'+kl)) i'
       )
              
            
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
    numRows = fst . shape
    numCols = snd . shape

    herm a = let h' = (not . isHerm) a
             in coerceBanded $ a{ isHerm=h' }
        
instance Tensor (BMatrix t (m,n)) (Int,Int) e where
    shape a | isHerm a  = (size2 a, size1 a)
            | otherwise = (size1 a, size2 a)

    bounds a = let (m,n) = shape a in ((0,0), (m-1,n-1))
    
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
    
    getSize (BM _ _ m n kl ku _ _) = 
        return $ foldl' (+) 0 $ 
             map (diagLen (m,n)) [(-kl)..ku]
    
    unsafeReadElem a (i,j)
        | isHerm a = 
            unsafeReadElem (herm a) (j,i) >>= return . E.conj
        | hasStorage a (i,j) =
            withForeignPtr (fptrOf a) $ \ptr ->
                peekElemOff ptr (indexOf a (i,j))
        | otherwise =
            return 0

    getIndices a =
        return $ filter (\ij -> inlinePerformIO $ canModifyElem (unsafeThaw a) ij)
                        (range $ bounds a)

    getElems a = getAssocs a >>= return . (map snd)

    getAssocs a = do
        is <- unsafeInterleaveIO $ getIndices a
        mapM (\i -> unsafeReadElem a i >>= \e -> return (i,e)) is

instance (BLAS1 e) => MTensor (BMatrix Mut (m,n)) (Int,Int) e IO where
    setZero a = case toRawMatrix a of (_,_,a',_) -> setZero a'
    
    setConstant e a = 
        let (_,_,a',h) = toRawMatrix a
            e' = if h then E.conj e else e
        in setConstant e' a'
        
    canModifyElem a ij = return $ hasStorage a ij

    unsafeWriteElem a (i,j) e
        | isHerm a = 
            unsafeWriteElem (herm a) (j,i) (E.conj e)
        | otherwise =
            withForeignPtr (fptrOf a) $ \ptr ->
                pokeElemOff ptr (indexOf a (i,j)) e

    modifyWith f a = 
        let (_,_,a',h) = toRawMatrix a
        in if h then modifyWith f a'
                else modifyWith f (herm a)

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
             