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
    ldaOf,
    isHerm,
    
    -- * To and from the underlying storage matrix
    toMatrix,
    fromMatrix,
    
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
    getRow,
    getCol,
    
    -- * Vector views
    diag,
    rowView,
    colView,
    
    -- * Casting matrices
    coerceMatrix,
    
    -- * Unsafe operations
    unsafeBanded,
    unsafeNewBanded,
    unsafeFreeze,
    unsafeThaw,
    unsafeWithElemPtr,
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

import Data.Matrix.Dense.Internal ( DMatrix )
import qualified Data.Matrix.Dense.Internal as M
                                
import Data.Vector.Dense.Internal ( DVector, Vector, conj, dim, newListVector )
import qualified Data.Vector.Dense.Internal as V

        
data BMatrix t mn e 
    = BM { fptr   :: !(ForeignPtr e)
         , offset :: !Int
         , size1  :: !Int
         , size2  :: !Int
         , lowBW  :: !Int
         , upBW   :: !Int
         , lda    :: !Int
         }
    | H !(BMatrix t mn e) 

type Banded = BMatrix Imm
type IOBanded = BMatrix Mut

fromForeignPtr :: ForeignPtr e -> Int -> (Int,Int) -> (Int,Int) -> Int 
    -> BMatrix t (m,n) e
fromForeignPtr f o (m,n) (kl,ku) l = BM f o m n kl ku l

toForeignPtr :: BMatrix t (m,n) e -> (ForeignPtr e, Int, (Int,Int), (Int,Int), Int)
toForeignPtr (H a) = toForeignPtr a
toForeignPtr (BM f o m n kl ku l) = (f, o, (m,n), (kl,ku), l)

ldaOf :: BMatrix t (m,n) e -> Int
ldaOf (H a) = ldaOf a
ldaOf a     = lda a

isHerm :: BMatrix t (m,n) e -> Bool
isHerm (H a) = not (isHerm a)
isHerm _     = False

unsafeFreeze :: BMatrix t mn e -> Banded mn e
unsafeFreeze = unsafeCoerce

unsafeThaw :: BMatrix t mn e -> IOBanded mn e
unsafeThaw = unsafeCoerce

-- | Coerce the phantom shape type from one type to another.
coerceMatrix :: BMatrix t mn e -> BMatrix t kl e
coerceMatrix = unsafeCoerce

toMatrix :: (Elem e) => BMatrix t (m,n) e -> (DMatrix t (m',n') e, (Int,Int), (Int,Int))
toMatrix (H a) = 
    case toMatrix (herm a) of (b, (m,n), (kl,ku)) -> (herm b, (n,m), (ku,kl))
toMatrix (BM f o m n kl ku ld) = 
    (M.fromForeignPtr f o (kl+1+ku,n) ld, (m,n), (kl,ku))

fromMatrix :: (Elem e) => DMatrix t (m,n) e -> (Int,Int) -> (Int,Int) -> BMatrix t (m',n') e
fromMatrix a (m,n) (kl,ku) = case a of
    (M.H a') -> herm (fromMatrix a' (n,m) (ku,kl))
    _        -> 
        let (f,o,(m',n'),ld) = M.toForeignPtr a
        in case undefined of
            _ | m' /= kl+1+ku -> 
                error $ "fromMatrix: number of rows must be equal to number of diagonals"
            _ | n' /= n ->
                error $ "fromMatrix: numbers of columns must be equal"
            _ ->
                BM f o m n kl ku ld
                

bandwidth :: BMatrix t (m,n) e -> (Int,Int)
bandwidth a = 
    let (kl,ku) = (numLower a, numUpper a)
    in (negate kl, ku)

numLower :: BMatrix t (m,n) e -> Int
numLower (H a) = numUpper a
numLower a     = lowBW a

numUpper :: BMatrix t (m,n) e -> Int
numUpper (H a) = numLower a
numUpper a     = upBW a


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
        in do    
            ptr <- mallocForeignPtrArray (m' * n)
            return $ fromForeignPtr ptr off (m,n) (kl,ku) l
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
    withForeignPtr (fptr x) $ flip clearArray ((kl+1+ku)*n)
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
unsafeDiag (H a) d = 
    conj $ unsafeDiag a (negate d)
unsafeDiag a d =
    let f      = fptr a
        off    = indexOf a (diagStart d)
        len    = diagLen (shape a) d
        stride = lda a
    in V.fromForeignPtr f off len stride
        
diag :: (Elem e) => BMatrix t (m,n) e -> Int -> DVector t k e
diag a = checkedDiag (shape a) (unsafeDiag a) 


indexOf :: BMatrix t mn e -> (Int,Int) -> Int
indexOf (H a) (i,j) = indexOf a (j,i)
indexOf (BM _ off _ _ _ ku ld) (i,j) =
    off + ku + (i - j) + j * ld
    --off + i * tda + (j - i + kl)
           

            
unsafeWithElemPtr :: (Elem e) => BMatrix t (m,n) e -> (Int,Int) -> (Ptr e -> IO a) -> IO a
unsafeWithElemPtr a (i,j) f = case a of
    (H a') -> unsafeWithElemPtr a' (j,i) f
    _      -> withForeignPtr (fptr a) $ \ptr ->
                  f $ ptr `advancePtr` (indexOf a (i,j))

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


unsafeColView :: (Elem e) => BMatrix t (m,n) e -> Int -> (Int, DVector t k e, Int)
unsafeColView (BM f off m _ kl ku ld) j =
    let nb     = max (j - ku)         0
        na     = max (m - 1 - j - kl) 0
        r      = max (ku - j) 0 
        c      = j 
        off'   = off + r + c * ld
        stride = 1
        len    = m - (nb + na)
    in if len >= 0
        then (nb, V.fromForeignPtr f off' len stride, na)
        else (m , V.fromForeignPtr f off' 0   stride,  0)
      
unsafeColView a j = 
    case unsafeRowView (herm a) j of (nb, v, na) -> (nb, conj v, na)

unsafeRowView :: (Elem e) => BMatrix t (m,n) e -> Int -> (Int, DVector t k e, Int)
unsafeRowView (BM f off _ n kl ku ld) i =
    let nb     = max (i - kl)         0
        na     = max (n - 1 - i - ku) 0
        r      = min (ku + i)         (kl + ku)
        c      = max (i - kl)         0 
        off'   = off + r + c * ld
        stride = ld - 1
        len    = n - (nb + na)
    in if len >= 0 
        then (nb, V.fromForeignPtr f off' len stride, na)
        else (n , V.fromForeignPtr f off' 0   stride,  0)

unsafeRowView a i =
    case unsafeColView (herm a) i of (nb, v, na) -> (nb, conj v, na)


rowView :: (Elem e) => BMatrix t (m,n) e -> Int -> (Int, DVector t k e, Int)
rowView a = checkedRow (shape a) (unsafeRowView a) 

colView :: (Elem e) => BMatrix t (m,n) e -> Int -> (Int, DVector t k e, Int)
colView a = checkedCol (shape a) (unsafeColView a)


instance C.Matrix (BMatrix t) where
    numRows = fst . shape
    numCols = snd . shape

    herm a = case a of
        (H a')   -> coerceMatrix a'
        _        -> H (coerceMatrix a)
        
instance Tensor (BMatrix t (m,n)) (Int,Int) e where
    shape a = case a of
        (H a')   -> case shape a' of (m,n) -> (n,m)
        _        -> (size1 a, size2 a)

    bounds a = let (m,n) = shape a in ((0,0), (m-1,n-1))
    
instance (BLAS1 e) => ITensor (BMatrix Imm (m,n)) (Int,Int) e where
    size = inlinePerformIO . getSize
    
    unsafeAt a = inlinePerformIO . (unsafeReadElem a)
    
    indices = inlinePerformIO . getIndices
    elems   = inlinePerformIO . getElems
    assocs  = inlinePerformIO . getAssocs

    (//)          = replaceHelp writeElem
    unsafeReplace = replaceHelp unsafeWriteElem

    amap f a = banded (shape a) (bandwidth a) ies
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
        let (a,mn,kl) = toMatrix b
        in do
            a' <- newCopy a
            return $ fromMatrix a' mn kl
    
    getSize a = case a of
        (H a')               -> getSize a'
        (BM _ _ m n kl ku _) -> return $ foldl' (+) 0 $ 
                                         map (diagLen (m,n)) [(-kl)..ku]
    
    unsafeReadElem a (i,j) = case a of
        (H a') -> unsafeReadElem a' (j,i) >>= return . E.conj
        _      -> withForeignPtr (fptr a) $ \ptr ->
                      peekElemOff ptr (indexOf a (i,j))

    getIndices a =
        return $ filter (\ij -> inlinePerformIO $ canModifyElem (unsafeThaw a) ij)
                        (range $ bounds a)

    getElems a = getAssocs a >>= return . (map snd)

    getAssocs a = do
        is <- unsafeInterleaveIO $ getIndices a
        mapM (\i -> unsafeReadElem a i >>= \e -> return (i,e)) is

instance (BLAS1 e) => MTensor (BMatrix Mut (m,n)) (Int,Int) e IO where
    setZero a = case toMatrix a of (a',_,_) -> setZero a'
    setConstant e a = case toMatrix a of (a',_,_) -> setConstant e a'
        
    canModifyElem a (i,j) = case a of
        (H a')               -> canModifyElem a' (j,i)
        (BM _ _ m n kl ku _) -> return $ inRange (0,m-1) i && 
                                         inRange (0,n-1) j &&
                                         inRange (max 0 (j-ku), min (m-1) (j+kl)) i
    
    unsafeWriteElem a (i,j) e = case a of
        (H a') -> unsafeWriteElem a' (j,i) (E.conj e)
        _      -> withForeignPtr (fptr a) $ \ptr ->
                      pokeElemOff ptr (indexOf a (i,j)) e

    modifyWith f a = case toMatrix a of (a',_,_) -> modifyWith f a'
