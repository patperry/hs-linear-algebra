{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.HermBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.HermBase
    where

import Unsafe.Coerce

import BLAS.Internal( checkSquare )
import Data.Matrix.Class
import Data.Tensor.Class
import BLAS.Types ( UpLoEnum(..) )

-- | A hermitian view of an underlying matrix.  The view can either be
-- of the upper or lower triangular part of the matrix.  The type arguments
-- are as follows:
--
--     * @a@: the underlyting matrix type.
--
--     * @nn@: a phantom type for the shape of the matrix.
--
--     * @e@: the element type of the matrix.
--
data Herm a nn e = Herm UpLoEnum (a nn e)

-- | Cast the phantom shape type.
coerceHerm :: Herm a mn e -> Herm a mn' e
coerceHerm = unsafeCoerce

-- | ApplyVector a function to the unerlying matrix.
mapHerm :: (a nn e -> b nn' e) -> Herm a nn e -> Herm b nn' e
mapHerm f (Herm u a) = Herm u $ f a

-- | Convert from a base matrix type to a Herm matrix type.
hermFromBase :: UpLoEnum -> a (n,n) e -> Herm a (n,n) e
hermFromBase = Herm
        
-- | Convert from a Herm matrix type to a base matrix type.        
hermToBase :: Herm a (n,n) e -> (UpLoEnum, a (n,n) e)
hermToBase (Herm u a) = (u,a)

-- | Construct a lower-triangular hermitian view into a matrix.  This also
-- checks to see if the base matrix is square.
hermL :: (MatrixShaped a) => a (n,n) e -> Herm a (n,n) e
hermL a = checkSquare "hermL" (shape a) (Herm Lower a)

-- | Construct an upper-triangular hermitian view into a matrix.  This also
-- checks to see if the base matrix is square.
hermU :: (MatrixShaped a) => a (n,n) e -> Herm a (n,n) e
hermU a = checkSquare "hermU" (shape a) (Herm Upper a)
      
instance (MatrixShaped a) => Shaped (Herm a) (Int,Int) where
    shape  (Herm _ a) = (n,n)             where n = min (numRows a) (numCols a)
    {-# INLINE shape #-}
    bounds (Herm _ a) = ((0,0),(n-1,n-1)) where n = min (numRows a) (numCols a)
    {-# INLINE bounds #-}
      
instance (MatrixShaped a) => MatrixShaped (Herm a) where

instance (HasHerm a) => HasHerm (Herm a) where
    herm = coerceHerm
    {-# INLINE herm #-}
    
instance Show (a (n,n) e) => Show (Herm a (n,n) e) where
    show (Herm u a) = constructor ++ " (" ++ show a ++ ")"
      where
        constructor = case u of
            Lower -> "hermL"
            Upper -> "hermU"
