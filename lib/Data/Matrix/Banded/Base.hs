{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances,
        TypeFamilies, Rank2Types #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Base
    where

import Control.Monad
import Control.Monad.ST
import Data.AEq
import Data.Ix
import System.IO.Unsafe
import Unsafe.Coerce

import BLAS.Internal( clearArray, checkedRow, checkedCol, checkedDiag,
    diagLen, inlinePerformIO )

import Data.Elem.BLAS( Elem, BLAS1, BLAS2, BLAS3 )

import Data.Tensor.Class
import Data.Tensor.Class.ITensor
import Data.Tensor.Class.MTensor

import Data.Matrix.Class
import Data.Matrix.Class.IMatrixBase
import Data.Matrix.Class.MMatrixBase
import Data.Matrix.Class.ISolveBase
import Data.Matrix.Class.MSolveBase

import Data.Matrix.Herm
import Data.Matrix.Tri

import Data.Vector.Dense.ST( runSTVector )
import Data.Vector.Dense.Base( BaseVector, ReadVector, WriteVector, Vector(..), 
    dim, unsafeVectorToIOVector, unsafePerformIOWithVector,
    unsafeConvertIOVector, newZeroVector, newCopyVector )
import Data.Matrix.Dense.ST( runSTMatrix )
import Data.Matrix.Dense.Base( BaseMatrix, ReadMatrix, WriteMatrix, Matrix(..),
    unsafeMatrixToIOMatrix, unsafePerformIOWithMatrix )

import Data.Matrix.Banded.IOBase( IOBanded )
import qualified Data.Matrix.Banded.IOBase as IO


-- | Immutable banded matrices.  The type arguments are as follows:
--
--     * @np@: a phantom type for the shape of the matrix.  Most functions
--       will demand that this be specified as a pair.  When writing a function
--       signature, you should always prefer @Banded (n,p) e@ to
--       @Banded np e@.
--
--     * @e@: the element type of the matrix.  Only certain element types
--       are supported.
--
newtype Banded np e = Banded (IOBanded np e)

freezeIOBanded :: (BLAS1 e) => IOBanded np e -> IO (Banded np e)
freezeIOBanded x = do
    y <- IO.newCopyIOBanded x
    return (Banded y)

thawIOBanded :: (BLAS1 e) => Banded np e -> IO (IOBanded np e)
thawIOBanded (Banded x) =
    IO.newCopyIOBanded x

unsafeFreezeIOBanded :: IOBanded np e -> IO (Banded np e)
unsafeFreezeIOBanded = return . Banded

unsafeThawIOBanded :: Banded np e -> IO (IOBanded np e)
unsafeThawIOBanded (Banded x) = return x

-- | Common functionality for all banded matrix types.
class ( MatrixShaped a e, HasVectorView a, HasMatrixStorage a, Elem e
      , BaseVector (VectorView a) e, BaseMatrix (MatrixStorage a) e
      ) => BaseBanded a e where
    
    -- | Get the number of lower diagonals in the banded matrix.
    numLower :: a (n,p) e -> Int
    
    -- | Get the number of upper diagonals in the banded matrix
    numUpper :: a (n,p) e -> Int
        
    -- | Get the range of valid diagonals in the banded matrix.
    -- @bandwidthds a@ is equal to @(numLower a, numUpper a)@.
    bandwidths :: a (n,p) e -> (Int,Int)
          
    -- | Get the leading dimension of the underlying storage of the
    -- banded matrix.
    ldaBanded :: a (n,p) e -> Int
    
    -- | Indicate whether or not the banded matrix is transposed and
    -- conjugated.
    isHermBanded :: a (n,p) e -> Bool

    -- | Cast the shape type of the banded matrix.
    coerceBanded :: a np e -> a np' e
    coerceBanded = unsafeCoerce
    {-# INLINE coerceBanded #-}

    -- | Get a matrix with the underlying storage of the banded matrix.
    matrixBanded :: a (n,p) e -> MatrixStorage a (k,l) e

    -- | View a vector as a banded matrix of the given shape.  The vector
    -- must have length equal to one of the specified dimensions.
    viewVectorAsBanded :: (Int,Int) -> VectorView a k e -> a (n,p) e
    
    -- | View a vector as a diagonal banded matrix.
    viewVectorAsDiagBanded :: VectorView a n e -> a (n,n) e
    viewVectorAsDiagBanded x = let
        n = dim x
        in viewVectorAsBanded (n,n) x
    {-# INLINE viewVectorAsBanded #-}
    
    -- | If the banded matrix has only a single diagonal, return a view
    -- into that diagonal.  Otherwise, return @Nothing@.
    maybeViewBandedAsVector :: a (n,p) e -> Maybe (VectorView a k e)

    -- | Given a shape and bandwidths, possibly view the elements stored
    -- in a dense matrix as a banded matrix.  This will fail unless
    -- the dense matrix has @isHermMatrix@ to be false, has the same
    -- number of columns as the desired banded matrix, and has number of
    -- rows equal to the desired number of diagonals.
    maybeBandedFromMatrix :: (Int,Int)
                          -> (Int,Int)
                          -> MatrixStorage a (k,p) e
                          -> Maybe (a (n,p) e)

    unsafeDiagViewBanded :: a (n,p) e -> Int -> VectorView a k e
    unsafeRowViewBanded :: a (n,p) e -> Int -> (Int, VectorView a k e, Int)
    unsafeColViewBanded :: a (n,p) e -> Int -> (Int, VectorView a k e, Int)
    
    -- | Unsafe cast from a matrix to an 'IOBanded'.
    unsafeBandedToIOBanded :: a (n,p) e -> IOBanded (n,p) e
    unsafeIOBandedToBanded :: IOBanded (n,p) e -> a (n,p) e


-- | Banded matrices that can be read in a monad.
class ( BaseBanded a e, BLAS2 e, ReadTensor a (Int,Int) e m
      , MMatrix a e m, MMatrix (Herm a) e m, MMatrix (Tri a) e m
      , MSolve (Tri a) e m
      , ReadVector (VectorView a) e m, ReadMatrix (MatrixStorage a) e m
      ) => ReadBanded a e m where

    -- | Cast the banded matrix to an 'IOBanded', perform an @IO@ action, and
    -- convert the @IO@ action to an action in the monad @m@.  This
    -- operation is /very/ unsafe.
    unsafePerformIOWithBanded :: a (n,p) e -> (IOBanded (n,p) e -> IO r) -> m r

    -- | Convert a mutable banded matrix to an immutable one by taking a 
    -- complete copy of it.
    freezeBanded :: a (n,p) e -> m (Banded (n,p) e)
    unsafeFreezeBanded :: a (n,p) e -> m (Banded (n,p) e)

-- | Banded matrices that can be created or modified in a monad.
class ( ReadBanded a e m, WriteTensor a (Int,Int) e m
      , WriteVector (VectorView a) e m
      , WriteMatrix (MatrixStorage a) e m
      ) => WriteBanded a e m where

    -- | Creates a new banded matrix of the given shape and bandwidths.  
    -- The elements will be uninitialized.
    newBanded_ :: (Int,Int) -> (Int,Int) -> m (a (n,p) e)

    -- | Unsafely convert an 'IO' action that creates an 'IOBanded' into
    -- an action in @m@ that creates a matrix.
    unsafeConvertIOBanded :: IO (IOBanded (n,p) e) -> m (a (n,p) e)

    -- | Convert an immutable banded matrix to a mutable one by taking a 
    -- complete copy of it.
    thawBanded :: Banded (n,p) e -> m (a (n,p) e)
    unsafeThawBanded :: Banded (n,p) e -> m (a (n,p) e)

-- | Create a banded matrix with the given shape, bandwidths, and 
-- associations.  The indices in the associations list must all fall
-- in the bandwidth of the matrix.  Unspecified elements will be set
-- to zero.
newBanded :: (WriteBanded a e m) => 
    (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> m (a (n,p) e)
newBanded = newBandedHelp writeElem
{-# INLINE newBanded #-}

unsafeNewBanded :: (WriteBanded a e m) => 
    (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> m (a (n,p) e)
unsafeNewBanded = newBandedHelp unsafeWriteElem
{-# INLINE unsafeNewBanded #-}

newBandedHelp :: (WriteBanded a e m) => 
       (IOBanded (n,p) e -> (Int,Int) -> e -> IO ())
    -> (Int,Int) -> (Int,Int) -> [((Int,Int),e)] -> m (a (n,p) e)
newBandedHelp set (m,n) (kl,ku) ijes = 
    unsafeConvertIOBanded $ do
        x <- newBanded_ (m,n) (kl,ku)
        IO.withIOBanded x $ flip clearArray ((kl+1+ku)*n)
        mapM_ (uncurry $ set x) ijes
        return x
{-# INLINE newBandedHelp #-}

-- | Create a banded matrix of the given shape and bandwidths by specifying
-- its diagonal elements.  The lists must all have the same length, equal
-- to the number of elements in the main diagonal of the matrix.  The 
-- sub-diagonals are specified first, then the super-diagonals.  In 
-- subdiagonal @i@, the first @i@ elements of the list are ignored.
newListsBanded :: (WriteBanded a e m) => 
    (Int,Int) -> (Int,Int) -> [[e]] -> m (a (n,p) e)
newListsBanded (m,n) (kl,ku) xs = do
    a <- newBanded_ (m,n) (kl,ku)
    zipWithM_ (writeDiagElems a) [(negate kl)..ku] xs
    return a
  where
    writeDiagElems :: (WriteBanded a e m) => a (n,p) e -> Int -> [e] -> m ()
    writeDiagElems a i es =
        let d   = unsafeDiagViewBanded a i
            nb  = max 0 (negate i)
            es' = drop nb es
        in zipWithM_ (unsafeWriteElem d) [0..(dim d - 1)] es'
{-# INLINE newListsBanded #-}

-- | Create a zero banded matrix with the specified shape and bandwidths.
newZeroBanded :: (WriteBanded a e m) => (Int,Int) -> (Int,Int) -> m (a (n,p) e)
newZeroBanded mn bw = unsafeConvertIOBanded $
    IO.newZeroIOBanded mn bw
{-# INLINE newZeroBanded #-}
 
-- | Create a constant banded matrix of the specified shape and bandwidths.
newConstantBanded :: (WriteBanded a e m) 
                  => (Int,Int) -> (Int,Int) -> e -> m (a (n,p) e)
newConstantBanded mn bw e = unsafeConvertIOBanded $
    IO.newConstantIOBanded mn bw e
{-# INLINE newConstantBanded #-}
 
-- | Set every element of a banded matrix to zero.
setZeroBanded :: (WriteBanded a e m) => a (n,p) e -> m ()
setZeroBanded a =
    unsafePerformIOWithBanded a $ IO.setZeroIOBanded
{-# INLINE setZeroBanded #-}
 
-- | Set every element of a banded matrix to a constant.
setConstantBanded :: (WriteBanded a e m) => e -> a (n,p) e -> m ()
setConstantBanded e a =
    unsafePerformIOWithBanded a $ IO.setConstantIOBanded e
{-# INLINE setConstantBanded #-}

-- | Create a new banded matrix by taking a copy of another one.
newCopyBanded :: (ReadBanded a e m, WriteBanded b e m)
              => a (n,p) e -> m (b (n,p) e)
newCopyBanded a = unsafeConvertIOBanded $
    IO.newCopyIOBanded (unsafeBandedToIOBanded a)
{-# INLINE newCopyBanded #-}

-- | Copy the elements of one banded matrix into another.  The two matrices
-- must have the same shape and badwidths.
copyBanded :: (WriteBanded b e m, ReadBanded a e m) =>
    b (n,p) e -> a (n,p) e -> m ()
copyBanded dst src
    | shape dst /= shape src =
        error "Shape mismatch in copyBanded."
    | bandwidths dst /= bandwidths src =
        error "Bandwidth mismatch in copyBanded."
    | otherwise =
        unsafeCopyBanded dst src
{-# INLINE copyBanded #-}

unsafeCopyBanded :: (WriteBanded b e m, ReadBanded a e m)
                 => b (n,p) e -> a (n,p) e -> m ()
unsafeCopyBanded dst src =
    unsafePerformIOWithBanded dst $ \dst' ->
        IO.unsafeCopyIOBanded dst' (unsafeBandedToIOBanded src)
{-# INLINE unsafeCopyBanded #-}

-- | Get a view of a diagonal of the banded matrix.  This will fail if
-- the index is outside of the bandwidth.
diagViewBanded :: (BaseBanded a e)
               => a (n,p) e -> Int -> VectorView a k e
diagViewBanded a i
    | i < -(numLower a) || i > numUpper a =
        error $ "Tried to get a diagonal view outside of the bandwidth."
    | otherwise =
        unsafeDiagViewBanded a i
{-# INLINE diagViewBanded #-}

-- | Get a view into the partial row of the banded matrix, along with the
-- number of zeros to pad before and after the view.
rowViewBanded :: (BaseBanded a e) => 
    a (n,p) e -> Int -> (Int, VectorView a k e, Int)
rowViewBanded a = checkedRow (shape a) (unsafeRowViewBanded a) 
{-# INLINE rowViewBanded #-}

-- | Get a view into the partial column of the banded matrix, along with the
-- number of zeros to pad before and after the view.
colViewBanded :: (BaseBanded a e) => 
    a (n,p) e -> Int -> (Int, VectorView a k e, Int)
colViewBanded a = checkedCol (shape a) (unsafeColViewBanded a)
{-# INLINE colViewBanded #-}

-- | Get a copy of the given diagonal of a banded matrix.
getDiagBanded :: (ReadBanded a e m, WriteVector y e m) =>
    a (n,p) e -> Int -> m (y k e)
getDiagBanded a i | i >= -kl && i <= ku =
                       newCopyVector $ diagViewBanded a i
                  | otherwise =
                       newZeroVector $ diagLen (m,n) i
  where
    (m,n)   = shape a
    (kl,ku) = bandwidths a
{-# INLINE getDiagBanded #-}
 
unsafeGetDiagBanded :: (ReadBanded a e m, WriteVector y e m) =>
    a (n,p) e -> Int -> m (y k e)
unsafeGetDiagBanded a i = 
    newCopyVector $ unsafeDiagViewBanded a i
{-# INLINE unsafeGetDiagBanded #-}
 
unsafeGetRowBanded :: (ReadBanded a e m, WriteVector y e m) =>
    a (n,p) e -> Int -> m (y p e)
unsafeGetRowBanded a i = unsafeConvertIOVector $
    IO.unsafeGetRowIOBanded (unsafeBandedToIOBanded a) i
{-# INLINE unsafeGetRowBanded #-}
 
unsafeGetColBanded :: (ReadBanded a e m, WriteVector y e m) =>
    a (n,p) e -> Int -> m (y n e)
unsafeGetColBanded a i = unsafeConvertIOVector $
    IO.unsafeGetColIOBanded (unsafeBandedToIOBanded a) i
{-# INLINE unsafeGetColBanded #-}

gbmv :: (ReadBanded a e m, ReadVector x e m, WriteVector y e m) =>
    e -> a (k,l) e -> x l e -> e -> y k e -> m ()
gbmv alpha a x beta y =
    unsafePerformIOWithVector y $
        IO.gbmv alpha (unsafeBandedToIOBanded a) (unsafeVectorToIOVector x) beta
{-# INLINE gbmv #-}

gbmm :: (ReadBanded a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    e -> a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
gbmm alpha a b beta c =
    unsafePerformIOWithMatrix c $
        IO.gbmm alpha (unsafeBandedToIOBanded a) (unsafeMatrixToIOMatrix b) beta
{-# INLINE gbmm #-}

hbmv :: (ReadBanded a e m, ReadVector x e m, WriteVector y e m) =>
    e -> Herm a (k,l) e -> x l e -> e -> y k e -> m ()
hbmv alpha a x beta y =
    unsafePerformIOWithVector y $
        IO.hbmv alpha (mapHerm unsafeBandedToIOBanded a) (unsafeVectorToIOVector x) beta
{-# INLINE hbmv #-}

hbmm :: (ReadBanded a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    e -> Herm a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
hbmm alpha a b beta c =
    unsafePerformIOWithMatrix c $
        IO.hbmm alpha (mapHerm unsafeBandedToIOBanded a) (unsafeMatrixToIOMatrix b) beta
{-# INLINE hbmm #-}

tbmv :: (ReadBanded a e m, WriteVector y e m) =>
    e -> Tri a (k,k) e -> y k e -> m ()
tbmv alpha a x =
    unsafePerformIOWithVector x $
        IO.tbmv alpha (mapTri unsafeBandedToIOBanded a)
{-# INLINE tbmv #-}

tbmm :: (ReadBanded a e m, WriteMatrix b e m) =>
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
tbmm alpha a b =
    unsafePerformIOWithMatrix b $
        IO.tbmm alpha (mapTri unsafeBandedToIOBanded a)
{-# INLINE tbmm #-}

tbmv' :: (ReadBanded a e m, ReadVector x e m, WriteVector y e m) =>
    e -> Tri a (k,l) e -> x l e -> e -> y k e -> m ()
tbmv' alpha a x beta y =
    unsafePerformIOWithVector y $
        IO.tbmv' alpha (mapTri unsafeBandedToIOBanded a) (unsafeVectorToIOVector x) beta
{-# INLINE tbmv' #-}

tbmm' :: (ReadBanded a e m, ReadMatrix b e m, WriteMatrix c e m) =>
    e -> Tri a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
tbmm' alpha a b beta c =
    unsafePerformIOWithMatrix c $
        IO.tbmm' alpha (mapTri unsafeBandedToIOBanded a) (unsafeMatrixToIOMatrix b) beta
{-# INLINE tbmm' #-}

tbsv :: (ReadBanded a e m, WriteVector y e m) =>
    e -> Tri a (k,k) e -> y k e -> m ()
tbsv alpha a x =
    unsafePerformIOWithVector x $
        IO.tbmv alpha (mapTri unsafeBandedToIOBanded a)
{-# INLINE tbsv #-}

tbsm :: (ReadBanded a e m, WriteMatrix b e m) =>
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
tbsm alpha a b =
    unsafePerformIOWithMatrix b $
        IO.tbsm alpha (mapTri unsafeBandedToIOBanded a)
{-# INLINE tbsm #-}

tbsv' :: (ReadBanded a e m, ReadVector y e m, WriteVector x e m)
      => e -> Tri a (k,l) e -> y k e -> x l e -> m ()
tbsv' alpha a y x = 
    unsafePerformIOWithVector x $
        IO.tbsv' alpha (mapTri unsafeBandedToIOBanded a) (unsafeVectorToIOVector y)
{-# INLINE tbsv' #-}

tbsm' :: (ReadBanded a e m, ReadMatrix c e m, WriteMatrix b e m) 
      => e -> Tri a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()
tbsm' alpha a c b =
    unsafePerformIOWithMatrix b $
        IO.tbsm' alpha (mapTri unsafeBandedToIOBanded a) (unsafeMatrixToIOMatrix c)
{-# INLINE tbsm' #-}

instance (Elem e) => BaseBanded IOBanded e where
    numLower = IO.numLowerIOBanded
    {-# INLINE numLower #-}
    numUpper = IO.numUpperIOBanded
    {-# INLINE numUpper #-}
    bandwidths = IO.bandwidthsIOBanded
    {-# INLINE bandwidths #-}
    ldaBanded = IO.ldaIOBanded
    {-# INLINE ldaBanded #-}
    isHermBanded = IO.isHermIOBanded
    {-# INLINE isHermBanded #-}
    matrixBanded = IO.matrixIOBanded
    {-# INLINE matrixBanded #-}
    maybeBandedFromMatrix = IO.maybeBandedFromIOMatrix
    {-# INLINE maybeBandedFromMatrix #-}
    viewVectorAsBanded = IO.viewVectorAsIOBanded
    {-# INLINE viewVectorAsBanded #-}
    maybeViewBandedAsVector = IO.maybeViewIOBandedAsVector
    {-# INLINE maybeViewBandedAsVector #-}
    unsafeDiagViewBanded = IO.unsafeDiagViewIOBanded
    {-# INLINE unsafeDiagViewBanded #-}
    unsafeRowViewBanded = IO.unsafeRowViewIOBanded
    {-# INLINE unsafeRowViewBanded #-}
    unsafeColViewBanded = IO.unsafeColViewIOBanded
    {-# INLINE unsafeColViewBanded #-}
    unsafeIOBandedToBanded = id
    {-# INLINE unsafeIOBandedToBanded #-}
    unsafeBandedToIOBanded = id
    {-# INLINE unsafeBandedToIOBanded #-}

instance (BLAS3 e) => ReadBanded IOBanded e IO where
    unsafePerformIOWithBanded a f = f a
    {-# INLINE unsafePerformIOWithBanded #-}
    freezeBanded = freezeIOBanded
    {-# INLINE freezeBanded #-}
    unsafeFreezeBanded = unsafeFreezeIOBanded
    {-# INLINE unsafeFreezeBanded #-}
    
instance (BLAS3 e) => WriteBanded IOBanded e IO where
    newBanded_ = IO.newIOBanded_
    {-# INLINE newBanded_ #-}
    unsafeConvertIOBanded = id
    {-# INLINE unsafeConvertIOBanded #-}
    thawBanded = thawIOBanded
    {-# INLINE thawBanded #-}
    unsafeThawBanded = unsafeThawIOBanded
    {-# INLINE unsafeThawBanded #-}

-- | Create a banded matrix with the given shape, bandwidths, and 
-- associations.  The indices in the associations list must all fall
-- in the bandwidth of the matrix.  Unspecified elements will be set
-- to zero.
banded :: (BLAS3 e) => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> Banded (n,p) e
banded mn kl ies = unsafePerformIO $
    unsafeFreezeIOBanded =<< newBanded mn kl ies
{-# NOINLINE banded #-}

unsafeBanded :: (BLAS3 e) => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> Banded (n,p) e
unsafeBanded mn kl ies = unsafePerformIO $
    unsafeFreezeIOBanded =<< unsafeNewBanded mn kl ies
{-# NOINLINE unsafeBanded #-}

-- | Create a banded matrix of the given shape and bandwidths by specifying
-- its diagonal elements.  The lists must all have the same length, equal
-- to the number of elements in the main diagonal of the matrix.  The 
-- sub-diagonals are specified first, then the super-diagonals.  In 
-- subdiagonal @i@, the first @i@ elements of the list are ignored.
listsBanded :: (BLAS3 e) => (Int,Int) -> (Int,Int) -> [[e]] -> Banded (n,p) e
listsBanded mn kl xs = unsafePerformIO $
    unsafeFreezeIOBanded =<< newListsBanded mn kl xs
{-# NOINLINE listsBanded #-}

-- | Create a zero banded matrix with the specified shape and bandwidths.
zeroBanded :: (BLAS3 e) => (Int,Int) -> (Int,Int) -> Banded (n,p) e
zeroBanded mn kl = unsafePerformIO $
    unsafeFreezeIOBanded =<< newZeroBanded mn kl
{-# NOINLINE zeroBanded #-}

-- | Create a constant banded matrix of the specified shape and bandwidths.
constantBanded :: (BLAS3 e) => (Int,Int) -> (Int,Int) -> e -> Banded (n,p) e
constantBanded mn kl e = unsafePerformIO $
    unsafeFreezeIOBanded =<< newConstantBanded mn kl e
{-# INLINE constantBanded #-}

-- | Create a banded matrix from a vector.  The vector must have length
-- equal to one of the specified dimension sizes.
bandedFromVector :: (Elem e) => (Int,Int) -> Vector k e -> Banded (n,p) e
bandedFromVector = viewVectorAsBanded
{-# INLINE bandedFromVector #-}

-- | Create a diagonal banded matrix from a vector.
diagBandedFromVector :: (Elem e) => Vector n e -> Banded (n,n) e
diagBandedFromVector = viewVectorAsDiagBanded
{-# INLINE diagBandedFromVector #-}

-- | Convert a diagonal banded matrix to a vector.  Fail if the banded
-- matrix has more than one diagonal
maybeVectorFromBanded :: (Elem e) => Banded (n,p) e -> Maybe (Vector k e)
maybeVectorFromBanded = maybeViewBandedAsVector
{-# INLINE maybeVectorFromBanded #-}

-- | Get a the given diagonal in a banded matrix.  Negative indices correspond 
-- to sub-diagonals.
diagBanded :: (BLAS1 e) => Banded (n,p) e -> Int -> Vector k e
diagBanded a = checkedDiag (shape a) (unsafeDiagBanded a)
{-# INLINE diagBanded #-}

unsafeDiagBanded :: (BLAS1 e) => Banded (n,p) e -> Int -> Vector k e
unsafeDiagBanded a i | i >= -kl && i <= ku = unsafeDiagViewBanded a i
               | otherwise = runSTVector $ 
    newZeroVector $ diagLen (shape a) i
  where
    (kl,ku) = bandwidths a
{-# INLINE unsafeDiagBanded #-}

listsFromBanded :: (BLAS1 e) => Banded np e -> ((Int,Int), (Int,Int),[[e]])
listsFromBanded a = ( (m,n)
            , (kl,ku)
            , map paddedDiag [(-kl)..ku]
            )
  where
    (m,n)   = shape a
    (kl,ku) = bandwidths (coerceBanded a)
    
    padBegin i   = replicate (max (-i) 0)    0
    padEnd   i   = replicate (max (m-n+i) 0) 0
    paddedDiag i = (  padBegin i
                   ++ elems (unsafeDiagViewBanded (coerceBanded a) i) 
                   ++ padEnd i 
                   )

instance (BLAS3 e) => Show (Banded (n,p) e) where
    show a 
        | isHermBanded a = 
           "herm (" ++ show (herm a) ++ ")"
        | otherwise = 
             let (mn,kl,es) = listsFromBanded a 
             in "listsBanded " ++ show mn ++ " " ++ show kl ++ " " ++ show es

compareBandedHelp :: (BLAS3 e) => 
    (e -> e -> Bool) -> Banded (n,p) e -> Banded (n,p) e -> Bool
compareBandedHelp cmp a b
    | shape a /= shape b =
        False
    | isHermBanded a == isHermBanded b && bandwidths a == bandwidths b =
        let elems' = if isHermBanded a then elems . herm
                                       else elems
        in
            and $ zipWith cmp (elems' a) (elems' b)
    | otherwise =
        let l = max (numLower a) (numLower b)
            u = max (numUpper a) (numUpper b)
        in
            and $ zipWith cmp (diagElems (-l,u) a) (diagElems (-l,u) b)
  where
    diagElems bw c = concatMap elems [ diagBanded c i | i <- range bw ]

instance (BLAS3 e, Eq e) => Eq (Banded (n,p) e) where
    (==) = compareBandedHelp (==)

instance (BLAS3 e, AEq e) => AEq (Banded (n,p) e) where
    (===) = compareBandedHelp (===)
    (~==) = compareBandedHelp (~==)


replaceBandedHelp :: (BLAS3 e) => 
       (forall n. IOBanded n e -> (Int,Int) -> e -> IO ())
    -> Banded mn e -> [((Int,Int), e)] -> Banded mn e
replaceBandedHelp set x ies = unsafePerformIO $ do
    y  <- newCopyBanded =<< unsafeThawIOBanded (coerceBanded x)
    mapM_ (uncurry $ set y) ies
    unsafeFreezeIOBanded (coerceBanded y)
{-# NOINLINE replaceBandedHelp #-}

instance (BLAS3 e) => ITensor Banded (Int,Int) e where
    (//)          = replaceBandedHelp writeElem
    {-# INLINE (//) #-}
    unsafeReplace = replaceBandedHelp unsafeWriteElem
    {-# INLINE unsafeReplace #-}
    unsafeAt (Banded a) i  = inlinePerformIO (unsafeReadElem a i)
    {-# INLINE unsafeAt #-}
    size (Banded a) = IO.sizeIOBanded a
    {-# INLINE size #-}
    elems (Banded a) = inlinePerformIO $ getElems a
    {-# INLINE elems #-}
    indices (Banded a) = IO.indicesIOBanded a
    {-# INLINE indices #-}
    assocs (Banded a) = inlinePerformIO $ getAssocs a
    {-# INLINE assocs #-}
    tmap f a      = coerceBanded $ listsBanded mn bw (map (map f) es)
      where (mn,bw,es) = listsFromBanded a

instance (BLAS3 e, Monad m) => ReadTensor Banded (Int,Int) e m where
    getSize = return . size
    {-# INLINE getSize #-}
    getAssocs = return . assocs
    {-# INLINE getAssocs #-}
    getIndices = return . indices
    {-# INLINE getIndices #-}
    getElems = return . elems
    {-# INLINE getElems #-}
    getAssocs' = return . assocs
    {-# INLINE getAssocs' #-}
    getIndices' = return . indices
    {-# INLINE getIndices' #-}
    getElems' = return . elems
    {-# INLINE getElems' #-}
    unsafeReadElem x i = return (unsafeAt x i)
    {-# INLINE unsafeReadElem #-}


instance HasVectorView Banded where
    type VectorView Banded = Vector

instance HasMatrixStorage Banded where
    type MatrixStorage Banded = Matrix

instance (Elem e) => Shaped Banded (Int,Int) e where
    shape (Banded a) = IO.shapeIOBanded a
    {-# INLINE shape #-}
    bounds (Banded a) = IO.boundsIOBanded a
    {-# INLINE bounds #-}

instance (Elem e) => MatrixShaped Banded e where
    herm (Banded a) = Banded $ IO.hermIOBanded a
    {-# INLINE herm #-}
    
instance (Elem e) => BaseBanded Banded e where
    numLower (Banded a) = IO.numLowerIOBanded a
    {-# INLINE numLower #-}
    numUpper (Banded a) = IO.numUpperIOBanded a
    {-# INLINE numUpper #-}
    bandwidths (Banded a) = IO.bandwidthsIOBanded a
    {-# INLINE bandwidths #-}
    ldaBanded (Banded a) = IO.ldaIOBanded a
    {-# INLINE ldaBanded #-}
    isHermBanded (Banded a) = IO.isHermIOBanded a
    {-# INLINE isHermBanded #-}
    matrixBanded (Banded a) = Matrix $ IO.matrixIOBanded a
    {-# INLINE matrixBanded #-}
    maybeBandedFromMatrix mn kl (Matrix a) = 
        liftM Banded $ IO.maybeBandedFromIOMatrix mn kl a
    {-# INLINE maybeBandedFromMatrix #-}
    viewVectorAsBanded mn (Vector x) = Banded $ IO.viewVectorAsIOBanded mn x
    {-# INLINE viewVectorAsBanded #-}
    maybeViewBandedAsVector (Banded a) = 
        liftM Vector $ IO.maybeViewIOBandedAsVector a
    {-# INLINE maybeViewBandedAsVector #-}    
    unsafeDiagViewBanded (Banded a) i = Vector $ IO.unsafeDiagViewIOBanded a i
    {-# INLINE unsafeDiagViewBanded #-}
    unsafeRowViewBanded (Banded a) i = 
        case IO.unsafeRowViewIOBanded a i of (nb,x,na) -> (nb, Vector x, na)
    {-# INLINE unsafeRowViewBanded #-}
    unsafeColViewBanded (Banded a) j = 
        case IO.unsafeColViewIOBanded a j of (nb,x,na) -> (nb, Vector x, na)
    {-# INLINE unsafeColViewBanded #-}
    unsafeIOBandedToBanded = Banded
    {-# INLINE unsafeIOBandedToBanded #-}
    unsafeBandedToIOBanded (Banded a) = a
    {-# INLINE unsafeBandedToIOBanded #-}

instance (BLAS3 e) => ReadBanded Banded e IO where
    unsafePerformIOWithBanded (Banded a) f = f a
    {-# INLINE unsafePerformIOWithBanded #-}
    freezeBanded (Banded a) = freezeIOBanded a
    {-# INLINE freezeBanded #-}
    unsafeFreezeBanded (Banded a) = unsafeFreezeIOBanded a
    {-# INLINE unsafeFreezeBanded #-}

instance (BLAS3 e) => MMatrix Banded e IO where
    unsafeDoSApplyAdd    = gbmv
    {-# INLINE unsafeDoSApplyAdd #-}
    unsafeDoSApplyAddMat = gbmm
    {-# INLINE unsafeDoSApplyAddMat #-}
    unsafeGetRow         = unsafeGetRowBanded
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol         = unsafeGetColBanded
    {-# INLINE unsafeGetCol #-}
    getRows = getRowsIO
    {-# INLINE getRows #-}
    getCols = getColsIO
    {-# INLINE getCols #-}

instance (BLAS3 e) => MMatrix (Herm Banded) e IO where
    unsafeDoSApplyAdd    = hbmv
    {-# INLINE unsafeDoSApplyAdd #-}    
    unsafeDoSApplyAddMat = hbmm
    {-# INLINE unsafeDoSApplyAddMat #-}    
    getRows = getRowsIO
    {-# INLINE getRows #-}
    getCols = getColsIO
    {-# INLINE getCols #-}

instance (BLAS3 e) => MMatrix (Tri Banded) e IO where
    unsafeDoSApply_      = tbmv
    {-# INLINE unsafeDoSApply_ #-}        
    unsafeDoSApplyMat_   = tbmm
    {-# INLINE unsafeDoSApplyMat_ #-}    
    unsafeDoSApplyAdd    = tbmv'
    {-# INLINE unsafeDoSApplyAdd #-}    
    unsafeDoSApplyAddMat = tbmm'
    {-# INLINE unsafeDoSApplyAddMat #-}    
    getRows = getRowsIO
    {-# INLINE getRows #-}
    getCols = getColsIO
    {-# INLINE getCols #-}

instance (BLAS3 e) => MSolve (Tri Banded) e IO where
    unsafeDoSSolve_    = tbsv
    {-# INLINE unsafeDoSSolve_ #-}    
    unsafeDoSSolveMat_ = tbsm
    {-# INLINE unsafeDoSSolveMat_ #-}    
    unsafeDoSSolve     = tbsv'
    {-# INLINE unsafeDoSSolve #-}
    unsafeDoSSolveMat  = tbsm'
    {-# INLINE unsafeDoSSolveMat #-}    

instance (BLAS3 e) => ReadBanded Banded e (ST s) where
    unsafePerformIOWithBanded (Banded a) f = unsafeIOToST $ f a
    {-# INLINE unsafePerformIOWithBanded #-}
    freezeBanded (Banded a) = unsafeIOToST $ freezeIOBanded a
    {-# INLINE freezeBanded #-}
    unsafeFreezeBanded (Banded a) = unsafeIOToST $ unsafeFreezeIOBanded a
    {-# INLINE unsafeFreezeBanded #-}

instance (BLAS3 e) => MMatrix Banded e (ST s) where
    unsafeDoSApplyAdd    = gbmv
    {-# INLINE unsafeDoSApplyAdd #-}
    unsafeDoSApplyAddMat = gbmm
    {-# INLINE unsafeDoSApplyAddMat #-}
    unsafeGetRow         = unsafeGetRowBanded
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol         = unsafeGetColBanded
    {-# INLINE unsafeGetCol #-}
    getRows = getRowsST
    {-# INLINE getRows #-}
    getCols = getColsST
    {-# INLINE getCols #-}

instance (BLAS3 e) => MMatrix (Herm Banded) e (ST s) where
    unsafeDoSApplyAdd    = hbmv
    {-# INLINE unsafeDoSApplyAdd #-}    
    unsafeDoSApplyAddMat = hbmm
    {-# INLINE unsafeDoSApplyAddMat #-}    
    getRows = getRowsST
    {-# INLINE getRows #-}
    getCols = getColsST
    {-# INLINE getCols #-}

instance (BLAS3 e) => MMatrix (Tri Banded) e (ST s) where
    unsafeDoSApply_      = tbmv
    {-# INLINE unsafeDoSApply_ #-}        
    unsafeDoSApplyMat_   = tbmm
    {-# INLINE unsafeDoSApplyMat_ #-}    
    unsafeDoSApplyAdd    = tbmv'
    {-# INLINE unsafeDoSApplyAdd #-}    
    unsafeDoSApplyAddMat = tbmm'
    {-# INLINE unsafeDoSApplyAddMat #-}    
    getRows = getRowsST
    {-# INLINE getRows #-}
    getCols = getColsST
    {-# INLINE getCols #-}

instance (BLAS3 e) => MSolve (Tri Banded) e (ST s) where
    unsafeDoSSolve_    = tbsv
    {-# INLINE unsafeDoSSolve_ #-}    
    unsafeDoSSolveMat_ = tbsm
    {-# INLINE unsafeDoSSolveMat_ #-}    
    unsafeDoSSolve     = tbsv'
    {-# INLINE unsafeDoSSolve #-}
    unsafeDoSSolveMat  = tbsm'
    {-# INLINE unsafeDoSSolveMat #-}    

instance (BLAS3 e) => IMatrix Banded e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    
    unsafeRow a i             = runSTVector $ unsafeGetRow a i
    unsafeCol a j             = runSTVector $ unsafeGetCol a j

instance (BLAS3 e) => IMatrix (Herm Banded) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    

instance (BLAS3 e) => IMatrix (Tri Banded) e where
    unsafeSApply alpha a x    = runSTVector $ unsafeGetSApply    alpha a x
    unsafeSApplyMat alpha a b = runSTMatrix $ unsafeGetSApplyMat alpha a b    

instance (BLAS3 e) => ISolve (Tri Banded) e where
    unsafeSSolve    alpha a y = runSTVector $ unsafeGetSSolve    alpha a y
    unsafeSSolveMat alpha a c = runSTMatrix $ unsafeGetSSolveMat alpha a c
    