{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances,
        RankNTypes #-}
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
import Unsafe.Coerce

import BLAS.Internal( clearArray, checkedRow, checkedCol )

import Data.Elem.BLAS( Elem, BLAS1, BLAS2, BLAS3 )

import Data.Matrix.Class
import Data.Matrix.Class.IMatrix
import Data.Matrix.Class.MMatrix
import Data.Matrix.Class.ISolve
import Data.Matrix.Class.MSolve

import Data.Matrix.Herm
import Data.Matrix.Tri

import Data.Matrix.Banded.IOBase
import Data.Vector.Dense.Class
import Data.Matrix.Dense.Class


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
    y <- newCopyIOBanded x
    return (Banded y)

thawIOBanded :: (BLAS1 e) => Banded np e -> IO (IOBanded np e)
thawIOBanded (Banded x) =
    newCopyIOBanded x

unsafeFreezeIOBanded :: IOBanded np e -> IO (Banded np e)
unsafeFreezeIOBanded = return . Banded

unsafeThawIOBanded :: Banded np e -> IO (IOBanded np e)
unsafeThawIOBanded (Banded x) = return x

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

    -- | Unsafe cast from an 'IOBanded' to a banded matrix.
    unsafeIOBandedToBanded :: IOBanded (n,p) e -> a (n,p) e
    
    -- | Unsafe cast from a matrix to an 'IOBanded'.
    unsafeBandedToIOBanded :: a (n,p) e -> IOBanded (n,p) e

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
        withIOBanded x $ flip clearArray ((kl+1+ku)*n)
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
    newZeroIOBanded mn bw
{-# INLINE newZeroBanded #-}
 
-- | Create a constant banded matrix of the specified shape and bandwidths.
newConstantBanded :: (WriteBanded a e m) 
                  => (Int,Int) -> (Int,Int) -> e -> m (a (n,p) e)
newConstantBanded mn bw e = unsafeConvertIOBanded $
    newConstantIOBanded mn bw e
{-# INLINE newConstantBanded #-}
 
-- | Set every element of a banded matrix to zero.
setZeroBanded :: (WriteBanded a e m) => a (n,p) e -> m ()
setZeroBanded a =
    unsafePerformIOWithBanded a $ setZeroIOBanded
{-# INLINE setZeroBanded #-}
 
-- | Set every element of a banded matrix to a constant.
setConstantBanded :: (WriteBanded a e m) => e -> a (n,p) e -> m ()
setConstantBanded e a =
    unsafePerformIOWithBanded a $ setConstantBanded e
{-# INLINE setConstantBanded #-}

newCopyBanded :: (ReadBanded a e m, WriteBanded b e m)
              => a (n,p) e -> m (b (n,p) e)
newCopyBanded a = unsafeConvertIOBanded $
    newCopyIOBanded (unsafeBandedToIOBanded a)
{-# INLINE newCopyBanded #-}

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
        unsafeCopyIOBanded dst' (unsafeBandedToIOBanded src)
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


instance (Elem e) => BaseBanded IOBanded e where
    numLower = numLowerIOBanded
    {-# INLINE numLower #-}
    numUpper = numUpperIOBanded
    {-# INLINE numUpper #-}
    bandwidths = bandwidthsIOBanded
    {-# INLINE bandwidths #-}
    ldaBanded = ldaIOBanded
    {-# INLINE ldaBanded #-}
    isHermBanded = isHermIOBanded
    {-# INLINE isHermBanded #-}
    matrixBanded = matrixIOBanded
    {-# INLINE matrixBanded #-}
    maybeBandedFromMatrix = maybeBandedFromIOMatrix
    {-# INLINE maybeBandedFromMatrix #-}
    unsafeDiagViewBanded = unsafeDiagViewIOBanded
    {-# INLINE unsafeDiagViewBanded #-}
    unsafeRowViewBanded = unsafeRowViewIOBanded
    {-# INLINE unsafeRowViewBanded #-}
    unsafeColViewBanded = unsafeColViewBanded
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
    newBanded_ = newIOBanded_
    {-# INLINE newBanded_ #-}
    unsafeConvertIOBanded = id
    {-# INLINE unsafeConvertIOBanded #-}
    thawBanded = thawIOBanded
    {-# INLINE thawBanded #-}
    unsafeThawBanded = unsafeThawIOBanded
    {-# INLINE unsafeThawBanded #-}
    