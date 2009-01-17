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
import Data.AEq
import Data.Maybe
import Data.Ix
import System.IO.Unsafe
import Unsafe.Coerce

import BLAS.Internal( clearArray, checkedRow, checkedCol, checkedDiag,
    diagLen, inlinePerformIO )

import Data.Elem.BLAS

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
    dim, unsafeIOVectorToVector, unsafeVectorToIOVector,
    unsafePerformIOWithVector, unsafeConvertIOVector, newZeroVector,
    newCopyVector )
import Data.Matrix.Dense.ST( runSTMatrix )
import Data.Matrix.Dense.Base( BaseMatrix, ReadMatrix, WriteMatrix, Matrix(..),
    unsafeMatrixToIOMatrix, unsafePerformIOWithMatrix )

import Data.Matrix.Banded.IOBase( IOBanded(..) )
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

freezeIOBanded :: IOBanded np e -> IO (Banded np e)
freezeIOBanded x = do
    y <- IO.newCopyIOBanded x
    return (Banded y)

thawIOBanded :: Banded np e -> IO (IOBanded np e)
thawIOBanded (Banded x) =
    IO.newCopyIOBanded x

unsafeFreezeIOBanded :: IOBanded np e -> IO (Banded np e)
unsafeFreezeIOBanded = return . Banded

unsafeThawIOBanded :: Banded np e -> IO (IOBanded np e)
unsafeThawIOBanded (Banded x) = return x

-- | Common functionality for all banded matrix types.
class ( MatrixShaped a, 
        HasVectorView a, 
        HasMatrixStorage a,
        BaseVector (VectorView a), 
        BaseMatrix (MatrixStorage a) ) => BaseBanded a where
    
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
    
    -- | Get the storage type of the banded matrix.
    transEnumBanded :: a (n,p) e -> TransEnum

    -- | Indicate whether or not the banded matrix storage is 
    -- transposed and conjugated.
    isHermBanded :: a (n,p) e -> Bool
    isHermBanded = (ConjTrans ==) . transEnumBanded
    {-# INLINE isHermBanded #-}

    -- | Cast the shape type of the banded matrix.
    coerceBanded :: a np e -> a np' e
    coerceBanded = unsafeCoerce
    {-# INLINE coerceBanded #-}

    -- | Get a matrix with the underlying storage of the banded matrix.
    -- This will fail if the banded matrix is hermed.
    maybeMatrixStorageFromBanded :: a (n,p) e -> Maybe (MatrixStorage a (k,p) e)

    -- | Given a shape and bandwidths, possibly view the elements stored
    -- in a dense matrix as a banded matrix.  This will if the matrix
    -- storage is hermed.  An error will be called if the number of rows
    -- in the matrix does not equal the desired number of diagonals or
    -- if the number of columns in the matrix does not equal the desired
    -- number of columns.
    maybeBandedFromMatrixStorage :: (Int,Int)
                                 -> (Int,Int)
                                 -> MatrixStorage a (k,p) e
                                 -> Maybe (a (n,p) e)

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

    unsafeDiagViewBanded :: a (n,p) e -> Int -> VectorView a k e
    unsafeRowViewBanded :: a (n,p) e -> Int -> (Int, VectorView a k e, Int)
    unsafeColViewBanded :: a (n,p) e -> Int -> (Int, VectorView a k e, Int)
    
    -- | Unsafe cast from a matrix to an 'IOBanded'.
    unsafeBandedToIOBanded :: a (n,p) e -> IOBanded (n,p) e
    unsafeIOBandedToBanded :: IOBanded (n,p) e -> a (n,p) e


-- | Banded matrices that can be read in a monad.
class ( BaseBanded a, 
        ReadTensor a (Int,Int) m,
        MMatrix a m, 
        MMatrix (Herm a) m, 
        MMatrix (Tri a) m, MSolve (Tri a) m,
        ReadVector (VectorView a) m,
        ReadMatrix (MatrixStorage a) m ) => ReadBanded a m where

    -- | Cast the banded matrix to an 'IOBanded', perform an @IO@ action, and
    -- convert the @IO@ action to an action in the monad @m@.  This
    -- operation is /very/ unsafe.
    unsafePerformIOWithBanded :: a (n,p) e -> (IOBanded (n,p) e -> IO r) -> m r

    -- | Convert a mutable banded matrix to an immutable one by taking a 
    -- complete copy of it.
    freezeBanded :: a (n,p) e -> m (Banded (n,p) e)
    unsafeFreezeBanded :: a (n,p) e -> m (Banded (n,p) e)

-- | Banded matrices that can be created or modified in a monad.
class ( ReadBanded a m, 
        WriteTensor a (Int,Int) m,
        WriteVector (VectorView a) m,
        WriteMatrix (MatrixStorage a) m ) => WriteBanded a m where

    -- | Creates a new banded matrix of the given shape and bandwidths.  
    -- The elements will be uninitialized.
    newBanded_ :: (Elem e) => (Int,Int) -> (Int,Int) -> m (a (n,p) e)

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
newBanded :: (WriteBanded a m, Elem e) => 
    (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> m (a (n,p) e)
newBanded = newBandedHelp writeElem
{-# INLINE newBanded #-}

unsafeNewBanded :: (WriteBanded a m, Elem e) => 
    (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> m (a (n,p) e)
unsafeNewBanded = newBandedHelp unsafeWriteElem
{-# INLINE unsafeNewBanded #-}

newBandedHelp :: (WriteBanded a m, Elem e) => 
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
newListsBanded :: (WriteBanded a m, Elem e) => 
    (Int,Int) -> (Int,Int) -> [[e]] -> m (a (n,p) e)
newListsBanded (m,n) (kl,ku) xs = do
    a <- newBanded_ (m,n) (kl,ku)
    zipWithM_ (writeDiagElems a) [(negate kl)..ku] xs
    return a
  where
    writeDiagElems :: (WriteBanded a m, Elem e) 
                   => a (n,p) e -> Int -> [e] -> m ()
    writeDiagElems a i es =
        let d   = unsafeDiagViewBanded a i
            nb  = max 0 (negate i)
            es' = drop nb es
        in zipWithM_ (unsafeWriteElem d) [0..(dim d - 1)] es'
{-# INLINE newListsBanded #-}

-- | Create a zero banded matrix with the specified shape and bandwidths.
newZeroBanded :: (WriteBanded a m, Elem e) 
              => (Int,Int) -> (Int,Int) -> m (a (n,p) e)
newZeroBanded mn bw = unsafeConvertIOBanded $
    IO.newZeroIOBanded mn bw
{-# INLINE newZeroBanded #-}
 
-- | Create a constant banded matrix of the specified shape and bandwidths.
newConstantBanded :: (WriteBanded a m, Elem e) 
                  => (Int,Int) -> (Int,Int) -> e -> m (a (n,p) e)
newConstantBanded mn bw e = unsafeConvertIOBanded $
    IO.newConstantIOBanded mn bw e
{-# INLINE newConstantBanded #-}
 
-- | Set every element of a banded matrix to zero.
setZeroBanded :: (WriteBanded a m) => a (n,p) e -> m ()
setZeroBanded a =
    unsafePerformIOWithBanded a $ IO.setZeroIOBanded
{-# INLINE setZeroBanded #-}
 
-- | Set every element of a banded matrix to a constant.
setConstantBanded :: (WriteBanded a m) => e -> a (n,p) e -> m ()
setConstantBanded e a =
    unsafePerformIOWithBanded a $ IO.setConstantIOBanded e
{-# INLINE setConstantBanded #-}

-- | Create a new banded matrix by taking a copy of another one.
newCopyBanded :: (ReadBanded a m, WriteBanded b m)
              => a (n,p) e -> m (b (n,p) e)
newCopyBanded a = unsafeConvertIOBanded $
    IO.newCopyIOBanded (unsafeBandedToIOBanded a)
{-# INLINE newCopyBanded #-}

-- | Copy the elements of one banded matrix into another.  The two matrices
-- must have the same shape and badwidths.
copyBanded :: (WriteBanded b m, ReadBanded a m) =>
    b (n,p) e -> a (n,p) e -> m ()
copyBanded dst src
    | shape dst /= shape src =
        error "Shape mismatch in copyBanded."
    | bandwidths dst /= bandwidths src =
        error "Bandwidth mismatch in copyBanded."
    | otherwise =
        unsafeCopyBanded dst src
{-# INLINE copyBanded #-}

unsafeCopyBanded :: (WriteBanded b m, ReadBanded a m)
                 => b (n,p) e -> a (n,p) e -> m ()
unsafeCopyBanded dst src =
    unsafePerformIOWithBanded dst $ \dst' ->
        IO.unsafeCopyIOBanded dst' (unsafeBandedToIOBanded src)
{-# INLINE unsafeCopyBanded #-}

-- | Get a view of a diagonal of the banded matrix.  This will fail if
-- the index is outside of the bandwidth.
diagViewBanded :: (BaseBanded a)
               => a (n,p) e -> Int -> VectorView a k e
diagViewBanded a i
    | i < -(numLower a) || i > numUpper a =
        error $ "Tried to get a diagonal view outside of the bandwidth."
    | otherwise =
        unsafeDiagViewBanded a i
{-# INLINE diagViewBanded #-}

-- | Get a view into the partial row of the banded matrix, along with the
-- number of zeros to pad before and after the view.
rowViewBanded :: (BaseBanded a) => 
    a (n,p) e -> Int -> (Int, VectorView a k e, Int)
rowViewBanded a = checkedRow (shape a) (unsafeRowViewBanded a) 
{-# INLINE rowViewBanded #-}

-- | Get a view into the partial column of the banded matrix, along with the
-- number of zeros to pad before and after the view.
colViewBanded :: (BaseBanded a) => 
    a (n,p) e -> Int -> (Int, VectorView a k e, Int)
colViewBanded a = checkedCol (shape a) (unsafeColViewBanded a)
{-# INLINE colViewBanded #-}

-- | Get a copy of the given diagonal of a banded matrix.
getDiagBanded :: (ReadBanded a m, WriteVector y m, Elem e) =>
    a (n,p) e -> Int -> m (y k e)
getDiagBanded a i | i >= -kl && i <= ku =
                       newCopyVector $ diagViewBanded a i
                  | otherwise =
                       newZeroVector $ diagLen (m,n) i
  where
    (m,n)   = shape a
    (kl,ku) = bandwidths a
{-# INLINE getDiagBanded #-}
 
unsafeGetDiagBanded :: (ReadBanded a m, WriteVector y m, Elem e) =>
    a (n,p) e -> Int -> m (y k e)
unsafeGetDiagBanded a i = 
    newCopyVector $ unsafeDiagViewBanded a i
{-# INLINE unsafeGetDiagBanded #-}
 
unsafeGetRowBanded :: (ReadBanded a m, WriteVector y m, Elem e) =>
    a (n,p) e -> Int -> m (y p e)
unsafeGetRowBanded a i = unsafeConvertIOVector $
    IO.unsafeGetRowIOBanded (unsafeBandedToIOBanded a) i
{-# INLINE unsafeGetRowBanded #-}
 
unsafeGetColBanded :: (ReadBanded a m, WriteVector y m, Elem e) =>
    a (n,p) e -> Int -> m (y n e)
unsafeGetColBanded a i = unsafeConvertIOVector $
    IO.unsafeGetColIOBanded (unsafeBandedToIOBanded a) i
{-# INLINE unsafeGetColBanded #-}

gbmv :: (ReadBanded a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    e -> a (k,l) e -> x l e -> e -> y k e -> m ()
gbmv alpha a x beta y =
    unsafePerformIOWithVector y $
        IO.gbmv alpha (unsafeBandedToIOBanded a) (unsafeVectorToIOVector x) beta
{-# INLINE gbmv #-}

gbmm :: (ReadBanded a m, ReadMatrix b m, WriteMatrix c m, BLAS2 e) =>
    e -> a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
gbmm alpha a b beta c =
    unsafePerformIOWithMatrix c $
        IO.gbmm alpha (unsafeBandedToIOBanded a) (unsafeMatrixToIOMatrix b) beta
{-# INLINE gbmm #-}

unsafeGetColHermBanded :: (ReadBanded a m, WriteVector x m, Elem e)
                      => Herm a (n,p) e -> Int -> m (x n e)
unsafeGetColHermBanded a i = unsafeConvertIOVector $
    IO.unsafeGetColHermIOBanded (mapHerm unsafeBandedToIOBanded a) i
{-# INLINE unsafeGetColHermBanded #-}

unsafeGetRowHermBanded :: (ReadBanded a m, WriteVector x m, Elem e)
                      => Herm a (n,p) e -> Int -> m (x p e)
unsafeGetRowHermBanded a i = unsafeConvertIOVector $
    IO.unsafeGetRowHermIOBanded (mapHerm unsafeBandedToIOBanded a) i
{-# INLINE unsafeGetRowHermBanded #-}

hbmv :: (ReadBanded a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    e -> Herm a (k,l) e -> x l e -> e -> y k e -> m ()
hbmv alpha a x beta y =
    unsafePerformIOWithVector y $
        IO.hbmv alpha (mapHerm unsafeBandedToIOBanded a) (unsafeVectorToIOVector x) beta
{-# INLINE hbmv #-}

hbmm :: (ReadBanded a m, ReadMatrix b m, WriteMatrix c m, BLAS2 e) =>
    e -> Herm a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
hbmm alpha a b beta c =
    unsafePerformIOWithMatrix c $
        IO.hbmm alpha (mapHerm unsafeBandedToIOBanded a) (unsafeMatrixToIOMatrix b) beta
{-# INLINE hbmm #-}

unsafeGetColTriBanded :: (ReadBanded a m, WriteVector x m, Elem e)
                      => Tri a (n,p) e -> Int -> m (x n e)
unsafeGetColTriBanded a i = unsafeConvertIOVector $
    IO.unsafeGetColTriIOBanded (mapTri unsafeBandedToIOBanded a) i
{-# INLINE unsafeGetColTriBanded #-}

unsafeGetRowTriBanded :: (ReadBanded a m, WriteVector x m, Elem e)
                      => Tri a (n,p) e -> Int -> m (x p e)
unsafeGetRowTriBanded a i = unsafeConvertIOVector $
    IO.unsafeGetRowTriIOBanded (mapTri unsafeBandedToIOBanded a) i
{-# INLINE unsafeGetRowTriBanded #-}

tbmv :: (ReadBanded a m, WriteVector y m, BLAS2 e) =>
    e -> Tri a (k,k) e -> y k e -> m ()
tbmv alpha a x =
    unsafePerformIOWithVector x $
        IO.tbmv alpha (mapTri unsafeBandedToIOBanded a)
{-# INLINE tbmv #-}

tbmm :: (ReadBanded a m, WriteMatrix b m, BLAS2 e) =>
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
tbmm alpha a b =
    unsafePerformIOWithMatrix b $
        IO.tbmm alpha (mapTri unsafeBandedToIOBanded a)
{-# INLINE tbmm #-}

tbmv' :: (ReadBanded a m, ReadVector x m, WriteVector y m, BLAS2 e) =>
    e -> Tri a (k,l) e -> x l e -> e -> y k e -> m ()
tbmv' alpha a x beta y =
    unsafePerformIOWithVector y $
        IO.tbmv' alpha (mapTri unsafeBandedToIOBanded a) (unsafeVectorToIOVector x) beta
{-# INLINE tbmv' #-}

tbmm' :: (ReadBanded a m, ReadMatrix b m, WriteMatrix c m, BLAS2 e) =>
    e -> Tri a (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
tbmm' alpha a b beta c =
    unsafePerformIOWithMatrix c $
        IO.tbmm' alpha (mapTri unsafeBandedToIOBanded a) (unsafeMatrixToIOMatrix b) beta
{-# INLINE tbmm' #-}

tbsv :: (ReadBanded a m, WriteVector y m, BLAS2 e) =>
    e -> Tri a (k,k) e -> y k e -> m ()
tbsv alpha a x =
    unsafePerformIOWithVector x $
        IO.tbmv alpha (mapTri unsafeBandedToIOBanded a)
{-# INLINE tbsv #-}

tbsm :: (ReadBanded a m, WriteMatrix b m, BLAS2 e) =>
    e -> Tri a (k,k) e -> b (k,l) e -> m ()
tbsm alpha a b =
    unsafePerformIOWithMatrix b $
        IO.tbsm alpha (mapTri unsafeBandedToIOBanded a)
{-# INLINE tbsm #-}

tbsv' :: (ReadBanded a m, ReadVector y m, WriteVector x m, BLAS2 e)
      => e -> Tri a (k,l) e -> y k e -> x l e -> m ()
tbsv' alpha a y x = 
    unsafePerformIOWithVector x $
        IO.tbsv' alpha (mapTri unsafeBandedToIOBanded a) (unsafeVectorToIOVector y)
{-# INLINE tbsv' #-}

tbsm' :: (ReadBanded a m, ReadMatrix c m, WriteMatrix b m, BLAS2 e) 
      => e -> Tri a (r,s) e -> c (r,t) e -> b (s,t) e -> m ()
tbsm' alpha a c b =
    unsafePerformIOWithMatrix b $
        IO.tbsm' alpha (mapTri unsafeBandedToIOBanded a) (unsafeMatrixToIOMatrix c)
{-# INLINE tbsm' #-}

instance BaseBanded IOBanded where
    numLower = IO.numLowerIOBanded
    {-# INLINE numLower #-}
    numUpper = IO.numUpperIOBanded
    {-# INLINE numUpper #-}
    bandwidths = IO.bandwidthsIOBanded
    {-# INLINE bandwidths #-}
    ldaBanded = IO.ldaIOBanded
    {-# INLINE ldaBanded #-}
    transEnumBanded = IO.transEnumIOBanded
    {-# INLINE transEnumBanded #-}
    maybeMatrixStorageFromBanded = IO.maybeMatrixStorageFromIOBanded
    {-# INLINE maybeMatrixStorageFromBanded #-}
    maybeBandedFromMatrixStorage = IO.maybeIOBandedFromMatrixStorage
    {-# INLINE maybeBandedFromMatrixStorage #-}
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

instance ReadBanded IOBanded IO where
    unsafePerformIOWithBanded a f = f a
    {-# INLINE unsafePerformIOWithBanded #-}
    freezeBanded = freezeIOBanded
    {-# INLINE freezeBanded #-}
    unsafeFreezeBanded = unsafeFreezeIOBanded
    {-# INLINE unsafeFreezeBanded #-}
    
instance WriteBanded IOBanded IO where
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
banded :: (Elem e) 
       => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> Banded (n,p) e
banded mn kl ies = unsafePerformIO $
    unsafeFreezeIOBanded =<< IO.newIOBanded mn kl ies

unsafeBanded :: (Elem e) 
             => (Int,Int) -> (Int,Int) -> [((Int,Int), e)] -> Banded (n,p) e
unsafeBanded mn kl ies = unsafePerformIO $
    unsafeFreezeIOBanded =<< IO.unsafeNewIOBanded mn kl ies

-- | Create a banded matrix of the given shape and bandwidths by specifying
-- its diagonal elements.  The lists must all have the same length, equal
-- to the number of elements in the main diagonal of the matrix.  The 
-- sub-diagonals are specified first, then the super-diagonals.  In 
-- subdiagonal @i@, the first @i@ elements of the list are ignored.
listsBanded :: (Elem e) 
            => (Int,Int) -> (Int,Int) -> [[e]] -> Banded (n,p) e
listsBanded mn kl xs = unsafePerformIO $
    unsafeFreezeIOBanded =<< IO.newListsIOBanded mn kl xs

-- | Create a zero banded matrix with the specified shape and bandwidths.
zeroBanded :: (Elem e)
           => (Int,Int) -> (Int,Int) -> Banded (n,p) e
zeroBanded mn kl = unsafePerformIO $
    unsafeFreezeIOBanded =<< IO.newZeroIOBanded mn kl

-- | Create a constant banded matrix of the specified shape and bandwidths.
constantBanded :: (Elem e) 
               => (Int,Int) -> (Int,Int) -> e -> Banded (n,p) e
constantBanded mn kl e = unsafePerformIO $
    unsafeFreezeIOBanded =<< IO.newConstantIOBanded mn kl e

-- | Create a banded matrix from a vector.  The vector must have length
-- equal to one of the specified dimension sizes.
bandedFromVector :: (Int,Int) -> Vector k e -> Banded (n,p) e
bandedFromVector = viewVectorAsBanded
{-# INLINE bandedFromVector #-}

-- | Create a diagonal banded matrix from a vector.
diagBandedFromVector :: Vector n e -> Banded (n,n) e
diagBandedFromVector = viewVectorAsDiagBanded
{-# INLINE diagBandedFromVector #-}

-- | Convert a diagonal banded matrix to a vector.  Fail if the banded
-- matrix has more than one diagonal
maybeVectorFromBanded :: Banded (n,p) e -> Maybe (Vector k e)
maybeVectorFromBanded = maybeViewBandedAsVector
{-# INLINE maybeVectorFromBanded #-}

-- | Get a the given diagonal in a banded matrix.  Negative indices correspond 
-- to sub-diagonals.
diagBanded :: Banded (n,p) e -> Int -> Vector k e
diagBanded a = checkedDiag (shape a) (unsafeDiagBanded a)
{-# INLINE diagBanded #-}

unsafeDiagBanded :: Banded (n,p) e -> Int -> Vector k e
unsafeDiagBanded a@(Banded (IOBanded _ _ _ _ _ _ _ _)) i 
    | i >= -kl && i <= ku = unsafeDiagViewBanded a i
    | otherwise           = runSTVector $ newZeroVector $ diagLen (shape a) i
  where
    (kl,ku) = bandwidths a
{-# INLINE unsafeDiagBanded #-}

listsFromBanded :: Banded np e -> ((Int,Int), (Int,Int),[[e]])
listsFromBanded a@(Banded (IOBanded _ _ _ _ _ _ _ _)) = ( (m,n)
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

instance (Show e) => Show (Banded (n,p) e) where
    show a 
        | isHermBanded a = 
           "herm (" ++ show (herm a) ++ ")"
        | otherwise = 
             let (mn,kl,es) = listsFromBanded a 
             in "listsBanded " ++ show mn ++ " " ++ show kl ++ " " ++ show es

compareBandedHelp :: (e -> e -> Bool) 
                  -> Banded (n,p) e -> Banded (n,p) e -> Bool
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

instance (Eq e) => Eq (Banded (n,p) e) where
    (==) = compareBandedHelp (==)

instance (AEq e) => AEq (Banded (n,p) e) where
    (===) = compareBandedHelp (===)
    (~==) = compareBandedHelp (~==)


replaceBandedHelp :: (forall n. IOBanded n e -> (Int,Int) -> e -> IO ())
                  -> Banded mn e -> [((Int,Int), e)] -> Banded mn e
replaceBandedHelp set x ies = unsafePerformIO $ do
    y  <- IO.newCopyIOBanded =<< unsafeThawIOBanded (coerceBanded x)
    mapM_ (uncurry $ set y) ies
    unsafeFreezeIOBanded (coerceBanded y)
{-# NOINLINE replaceBandedHelp #-}

tmapBanded :: (e -> e) -> Banded (n,p) e -> Banded (n,p) e
tmapBanded f a@(Banded (IOBanded _ _ _ _ _ _ _ _)) =
    case maybeMatrixStorageFromBanded a of
        Just a' -> fromJust $
            maybeBandedFromMatrixStorage (shape a) (bandwidths a) (tmap f a')
        Nothing ->
            let f' = conjugate . f . conjugate
            in herm (tmapBanded f' (herm a))


instance ITensor Banded (Int,Int) where
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
    tmap f = coerceBanded . tmapBanded f . coerceBanded
    {-# INLINE tmap #-}

instance (Monad m) => ReadTensor Banded (Int,Int) m where
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

instance Shaped Banded (Int,Int) where
    shape (Banded a) = IO.shapeIOBanded a
    {-# INLINE shape #-}
    bounds (Banded a) = IO.boundsIOBanded a
    {-# INLINE bounds #-}

instance MatrixShaped Banded where
    herm (Banded a) = Banded $ IO.hermIOBanded a
    {-# INLINE herm #-}
    
instance BaseBanded Banded where
    numLower (Banded a) = IO.numLowerIOBanded a
    {-# INLINE numLower #-}
    numUpper (Banded a) = IO.numUpperIOBanded a
    {-# INLINE numUpper #-}
    bandwidths (Banded a) = IO.bandwidthsIOBanded a
    {-# INLINE bandwidths #-}
    ldaBanded (Banded a) = IO.ldaIOBanded a
    {-# INLINE ldaBanded #-}
    transEnumBanded (Banded a) = IO.transEnumIOBanded a
    {-# INLINE transEnumBanded #-}
    maybeMatrixStorageFromBanded (Banded a) = liftM Matrix $ IO.maybeMatrixStorageFromIOBanded a
    {-# INLINE maybeMatrixStorageFromBanded #-}
    maybeBandedFromMatrixStorage mn kl (Matrix a) = 
        liftM Banded $ IO.maybeIOBandedFromMatrixStorage mn kl a
    {-# INLINE maybeBandedFromMatrixStorage #-}
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

-- The NOINLINE pragmas and the strictness annotations here are *really*
-- important.  Otherwise, the compiler might think that certain actions
-- don't need to be run.
instance (Monad m) => ReadBanded Banded m where
    freezeBanded (Banded a) = return $! unsafePerformIO $ freezeIOBanded a
    {-# NOINLINE freezeBanded #-}
    unsafeFreezeBanded = return . id
    {-# INLINE unsafeFreezeBanded #-}
    unsafePerformIOWithBanded (Banded a) f = return $! unsafePerformIO $ f a
    {-# NOINLINE unsafePerformIOWithBanded #-}

instance (Monad m) => MMatrix Banded m where
    unsafeDoSApplyAddVector    = gbmv
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = gbmm
    {-# INLINE unsafeDoSApplyAddMatrix #-}
    unsafeGetRow         = unsafeGetRowBanded
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol         = unsafeGetColBanded
    {-# INLINE unsafeGetCol #-}
    getRows =
        return . map (unsafeIOVectorToVector . unsafeVectorToIOVector ) . rows    
    {-# INLINE getRows #-}
    getCols = 
        return . map (unsafeIOVectorToVector . unsafeVectorToIOVector ) . cols
    {-# INLINE getCols #-}

instance (Monad m) => MMatrix (Herm Banded) m where
    unsafeDoSApplyAddVector    = hbmv
    {-# INLINE unsafeDoSApplyAddVector #-}    
    unsafeDoSApplyAddMatrix = hbmm
    {-# INLINE unsafeDoSApplyAddMatrix #-}    
    getRows =
        return . map (unsafeIOVectorToVector . unsafeVectorToIOVector ) . rows    
    {-# INLINE getRows #-}
    getCols = 
        return . map (unsafeIOVectorToVector . unsafeVectorToIOVector ) . cols
    {-# INLINE getCols #-}
    unsafeGetRow = unsafeGetRowHermBanded
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol = unsafeGetColHermBanded
    {-# INLINE unsafeGetCol #-}

instance (Monad m) => MMatrix (Tri Banded) m where
    unsafeDoSApplyVector_       = tbmv
    {-# INLINE unsafeDoSApplyVector_  #-}        
    unsafeDoSApplyMatrix_   = tbmm
    {-# INLINE unsafeDoSApplyMatrix_ #-}    
    unsafeDoSApplyAddVector    = tbmv'
    {-# INLINE unsafeDoSApplyAddVector #-}    
    unsafeDoSApplyAddMatrix = tbmm'
    {-# INLINE unsafeDoSApplyAddMatrix #-}    
    getRows =
        return . map (unsafeIOVectorToVector . unsafeVectorToIOVector ) . rows    
    {-# INLINE getRows #-}
    getCols = 
        return . map (unsafeIOVectorToVector . unsafeVectorToIOVector ) . cols
    {-# INLINE getCols #-}
    unsafeGetRow = unsafeGetRowTriBanded
    {-# INLINE unsafeGetRow #-}
    unsafeGetCol = unsafeGetColTriBanded
    {-# INLINE unsafeGetCol #-}

instance (Monad m) => MSolve (Tri Banded) m where
    unsafeDoSSolveVector_    = tbsv
    {-# INLINE unsafeDoSSolveVector_ #-}    
    unsafeDoSSolveMatrix_ = tbsm
    {-# INLINE unsafeDoSSolveMatrix_ #-}    
    unsafeDoSSolveVector     = tbsv'
    {-# INLINE unsafeDoSSolveVector #-}
    unsafeDoSSolveMatrix  = tbsm'
    {-# INLINE unsafeDoSSolveMatrix #-}    

instance IMatrix Banded where
    unsafeSApplyVector alpha a x = runSTVector $ unsafeGetSApplyVector    alpha a x
    {-# INLINE unsafeSApplyVector #-}
    unsafeSApplyMatrix alpha a b = runSTMatrix $ unsafeGetSApplyMatrix alpha a b    
    {-# INLINE unsafeSApplyMatrix #-}
    unsafeRow a i = runSTVector $ unsafeGetRow a i
    {-# INLINE unsafeRow #-}
    unsafeCol a j = runSTVector $ unsafeGetCol a j
    {-# INLINE unsafeCol #-}

instance IMatrix (Herm Banded) where
    unsafeSApplyVector alpha a x = runSTVector $ unsafeGetSApplyVector    alpha a x
    {-# INLINE unsafeSApplyVector #-}
    unsafeSApplyMatrix alpha a b = runSTMatrix $ unsafeGetSApplyMatrix alpha a b    
    {-# INLINE unsafeSApplyMatrix #-}
    unsafeRow a i = runSTVector $ unsafeGetRow a i
    {-# INLINE unsafeRow #-}
    unsafeCol a j = runSTVector $ unsafeGetCol a j
    {-# INLINE unsafeCol #-}

instance IMatrix (Tri Banded) where
    unsafeSApplyVector alpha a x = runSTVector $ unsafeGetSApplyVector    alpha a x
    {-# INLINE unsafeSApplyVector #-}
    unsafeSApplyMatrix alpha a b = runSTMatrix $ unsafeGetSApplyMatrix alpha a b    
    {-# INLINE unsafeSApplyMatrix #-}
    unsafeRow a i = runSTVector $ unsafeGetRow a i
    {-# INLINE unsafeRow #-}
    unsafeCol a j = runSTVector $ unsafeGetCol a j
    {-# INLINE unsafeCol #-}

instance ISolve (Tri Banded) where
    unsafeSSolveVector    alpha a y = runSTVector $ unsafeGetSSolveVector    alpha a y
    {-# NOINLINE unsafeSSolveVector #-}
    unsafeSSolveMatrix alpha a c = runSTMatrix $ unsafeGetSSolveMatrix alpha a c
    {-# NOINLINE unsafeSSolveMatrix #-}
