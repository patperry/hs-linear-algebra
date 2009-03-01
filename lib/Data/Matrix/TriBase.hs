{-# LANGUAGE MultiParamTypeClasses, FlexibleContexts, FlexibleInstances #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.TriBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.TriBase
    where

import Unsafe.Coerce

import Data.Matrix.Class
import Data.Tensor.Class

-- | A triangular or trapezoidal view of an underlying matrix.  The view 
-- can either be of the upper or lower triangular part of the matrix, and
-- can optionally include or exclude the diagonal.  If the diagonal enum 
-- is @Unit@, the diagonal entries of the underlying matrix are not
-- referenced, but are assumed to be @1@.  The type arguments are as follows:
--
--     * @a@: the underlyting matrix type.
--
--     * @e@: the element type of the matrix.
--
data Tri a e = Tri UpLoEnum DiagEnum (a e)

-- | Cast the phantom shape type.
coerceTri :: Tri a e -> Tri a e
coerceTri = unsafeCoerce

-- | ApplyVector a function to the base matrix.
mapTri :: (a e -> b e) -> Tri a e -> Tri b e
mapTri f (Tri u d a) = Tri u d $ f a

-- | Convert from a base matrix type to a triangular view.
triFromBase :: UpLoEnum -> DiagEnum -> a e -> Tri a e
triFromBase = Tri
        
-- | Convert from a triangular view to the base matrix.
triToBase :: Tri a e -> (UpLoEnum, DiagEnum, a e)
triToBase (Tri u d a) = (u,d,a)

-- | Get a lower triangular view of a matrix.
lower :: (MatrixShaped a) => a e -> Tri a e
lower = Tri Lower NonUnit

-- | Get a lower triangular view of a matrix, with unit diagonal.
lowerU :: (MatrixShaped a) => a e -> Tri a e
lowerU = Tri Lower Unit

-- | Get an upper triangular view of a matrix.
upper :: (MatrixShaped a) => a e -> Tri a e
upper = Tri Upper NonUnit

-- | Get an upper triangular view of a matrix, with unit diagonal.
upperU :: (MatrixShaped a) => a e -> Tri a e
upperU = Tri Upper Unit

      
instance (MatrixShaped a) => Shaped (Tri a) (Int,Int) where
    shape (Tri _ _ a) = (numRows a, numCols a)
    {-# INLINE shape #-}    
    bounds a = ((0,0),(m-1,n-1)) where (m,n) = shape a
    {-# INLINE bounds #-}
    
instance (MatrixShaped a) => MatrixShaped (Tri a) where

instance (HasHerm a) => HasHerm (Tri a) where
    herm (Tri u d a) = Tri (flipUpLo u) d (herm a)
    {-# INLINE herm #-}

instance (Show (a e), MatrixShaped a) => Show (Tri a e) where
    show (Tri u d a) =
        constructor ++ " (" ++ show a ++ ")"
        where
          constructor = case (u,d) of
              (Lower, NonUnit) -> "lower"
              (Lower, Unit   ) -> "lowerU"
              (Upper, NonUnit) -> "upper"
              (Upper, Unit   ) -> "upperU"
