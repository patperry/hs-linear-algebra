{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances #-}
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

import Unsafe.Coerce

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
    unsafeRowViewBanded :: a (n,p) e -> Int -> (Int, VectorView a p e, Int)
    unsafeColViewBanded :: a (n,p) e -> Int -> (Int, VectorView a n e, Int)

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
    