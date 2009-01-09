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
import BLAS.Types ( UpLoEnum(..), DiagEnum(..), flipUpLo )

-- | A triangular or trapezoidal view of an underlying matrix.  The view 
-- can either be of the upper or lower triangular part of the matrix, and
-- can optionally include or exclude the diagonal.  If the diagonal enum 
-- is @Unit@, the diagonal entries of the underlying matrix are not
-- referenced, but are assumed to be @1@.  The type arguments are as follows:
--
--     * @a@: the underlyting matrix type.
--
--     * @np@: a phantom type for the shape of the view.
--
--     * @e@: the element type of the matrix.
--
data Tri a np e = Tri UpLoEnum DiagEnum (a np e)

-- | Cast the phantom shape type.
coerceTri :: Tri a np e -> Tri a np' e
coerceTri = unsafeCoerce

-- | Apply a function to the base matrix.
mapTri :: (a np e -> b np' e) -> Tri a np e -> Tri b np' e
mapTri f (Tri u d a) = Tri u d $ f a

-- | Convert from a base matrix type to a triangular view.
triFromBase :: UpLoEnum -> DiagEnum -> a (n,p) e -> Tri a (n,p) e
triFromBase = Tri
        
-- | Convert from a triangular view to the base matrix.
triToBase :: Tri a (n,p) e -> (UpLoEnum, DiagEnum, a (n,p) e)
triToBase (Tri u d a) = (u,d,a)

-- | Get a lower triangular view of a matrix.
lower :: (MatrixShaped a e) => a (n,p) e -> Tri a (n,p) e
lower = Tri Lower NonUnit

-- | Get a lower triangular view of a matrix, with unit diagonal.
lowerU :: (MatrixShaped a e) => a (n,p) e -> Tri a (n,p) e
lowerU = Tri Lower Unit

-- | Get an upper triangular view of a matrix.
upper :: (MatrixShaped a e) => a (n,p) e -> Tri a (n,p) e
upper = Tri Upper NonUnit

-- | Get an upper triangular view of a matrix, with unit diagonal.
upperU :: (MatrixShaped a e) => a (n,p) e -> Tri a (n,p) e
upperU = Tri Upper Unit

      
instance MatrixShaped a e => Shaped (Tri a) (Int,Int) e where
    shape (Tri Lower _ a) = (numRows a, min (numRows a) (numCols a))
    shape (Tri Upper _ a) = (min (numRows a) (numCols a), numCols a)
    {-# INLINE shape #-}    
    bounds a = ((0,0),(m-1,n-1)) where (m,n) = shape a
    {-# INLINE bounds #-}
    
instance MatrixShaped a e => MatrixShaped (Tri a) e where
    herm (Tri u d a) = Tri (flipUpLo u) d (herm a)
    {-# INLINE herm #-}

instance (Show (a (n,p) e), MatrixShaped a e) => Show (Tri a (n,p) e) where
    show (Tri u d a) =
        constructor ++ suffix ++ " (" ++ show a ++ ")"
        where
          constructor = case (u,d) of
              (Lower, NonUnit) -> "lower"
              (Lower, Unit   ) -> "lowerU"
              (Upper, NonUnit) -> "upper"
              (Upper, Unit   ) -> "upperU"

          suffix = case undefined of
                       _ | isSquare a -> ""
                       _ | isFat a    -> "Fat"
                       _              -> "Tall"
