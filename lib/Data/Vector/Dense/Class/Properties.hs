-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Class.Properties
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Class.Properties (
    -- * Vector Properties
    getSumAbs,
    getNorm2,
    getWhichMaxAbs,
    getDot,
    unsafeGetDot,

    ) where

import BLAS.C( BLAS1 )
import Data.Elem.Conj
import BLAS.Internal (  checkVecVecOp )
import qualified BLAS.C as BLAS

import Data.Tensor.Class.MTensor( unsafeReadElem )
import Data.Vector.Dense.Class.Internal


-- | Gets the sum of the absolute values of the vector entries.
getSumAbs :: (ReadVector x e m) => x n e -> m Double
getSumAbs = vectorCall BLAS.asum
    
-- | Gets the 2-norm of a vector.
getNorm2 :: (ReadVector x e m) => x n e -> m Double
getNorm2 = vectorCall BLAS.nrm2

-- | Gets the index and norm of the element with maximum magnitude.  This is 
-- undefined if any of the elements are @NaN@.  It will throw an exception if 
-- the dimension of the vector is 0.
getWhichMaxAbs :: (ReadVector x e m) => x n e -> m (Int, e)
getWhichMaxAbs x =
    case (dim x) of
        0 -> fail $ "getWhichMaxAbs of an empty vector"
        _ -> do
            i <- vectorCall BLAS.iamax x
            e <- unsafeReadElem x i
            return (i,e)


-- | Computes the dot product of two vectors.
getDot :: (ReadVector x e m, ReadVector y e m) => 
    x n e -> y n e -> m e
getDot x y = checkVecVecOp "getDot" (dim x) (dim y) $ unsafeGetDot x y
{-# INLINE getDot #-}

unsafeGetDot :: (ReadVector x e m, ReadVector y e m) => 
    x n e -> y n e -> m e
unsafeGetDot x y =
    case (isConj x, isConj y) of
        (False, False) -> vectorCall2 BLAS.dotc x y
        (True , False) -> vectorCall2 BLAS.dotu x y
        (False, True ) -> vectorCall2 BLAS.dotu x y >>= return . conj
        (True , True)  -> vectorCall2 BLAS.dotc x y >>= return . conj
{-# INLINE unsafeGetDot #-}
