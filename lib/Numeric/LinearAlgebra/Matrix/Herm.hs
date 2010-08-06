-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix.Herm
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Hermitian views of matrices.
--

module Numeric.LinearAlgebra.Matrix.Herm (
    -- * Hermitian views of matrices
    Herm(..),
    withHerm,
    uploHerm,
    
    -- * Herm Matrix operations
    copyToHermMatrix,
    
    -- ** Immutable Matrix-Vector
    mulHermMatrixVector,
    mulHermMatrixVectorWithScale,
    mulHermMatrixAddVector,
    mulHermMatrixAddVectorWithScales,
    
    -- ** Immutable Matrix-Matrix
    mulHermMatrixMatrix,
    mulHermMatrixMatrixWithScale,
    mulHermMatrixAddMatrix,
    mulHermMatrixAddMatrixWithScales,

    -- ** Mutable Matrix-Vector
    mulHermMatrixToVector,
    mulHermMatrixToVectorWithScale,
    mulHermMatrixAddToVector,
    mulHermMatrixAddToVectorWithScales,
    
    -- ** Mutable Matrix-Matrix
    mulHermMatrixToMatrix,
    mulHermMatrixToMatrixWithScale,
    mulHermMatrixAddToMatrix,
    mulHermMatrixAddToMatrixWithScales,
    ) where

import Control.Monad( when )
import Control.Monad.ST( ST, unsafeIOToST )
import Text.Printf( printf )

import Numeric.LinearAlgebra.Vector.Base
import Numeric.LinearAlgebra.Vector.STBase
import Numeric.LinearAlgebra.Matrix.Base
import Numeric.LinearAlgebra.Matrix.STBase
import Numeric.LinearAlgebra.Types
import qualified Numeric.LinearAlgebra.Types.BLAS as BLAS

-- | A hermitian view of an underlying matrix.  The view can either be
-- of the upper or lower triangular part of the matrix.  The type arguments
-- are as follows:
--
--     * @m@: the underlyting matrix type.
--
--     * @e@: the element type of the matrix.
--
data Herm m e = Herm Uplo (m e) deriving (Show)

-- | Apply a function to the unerlying matrix.
withHerm :: (m e -> a) -> Herm m e -> a
withHerm f (Herm _ m) = f m

-- | Returns the @Uplo@ enum of the herm.
uploHerm :: Herm m e -> Uplo
uploHerm (Herm u _) = u

copyToHermMatrix :: (Storable e, RMatrix m)
                 => Herm m e
                 -> Herm (STMatrix s) e
                 -> ST s ()
copyToHermMatrix = undefined

-- | @mulHermMatrixVector a x@ returns @a * x@.
mulHermMatrixVector :: (BLAS2 e)
                    => Herm Matrix e
                    -> Vector e
                    -> Vector e
mulHermMatrixVector a x =
    runVector $ do
        y <- newVector_ (dimVector x)
        mulHermMatrixToVector a x y
        return y

-- | @mulHermMatrixVectorWithScale alpha a x@ retunrs @alpha * a * x@.
mulHermMatrixVectorWithScale :: (BLAS2 e)
                             => e
                             -> Herm Matrix e
                             -> Vector e
                             -> Vector e
mulHermMatrixVectorWithScale alpha a x =
    runVector $ do
        y <- newVector_ (dimVector x)
        mulHermMatrixToVectorWithScale alpha a x y
        return y
                       
-- | @mulHermMatrixAddVector a x y@ returns @a * x + y@.
mulHermMatrixAddVector :: (BLAS2 e)
                       => Herm Matrix e
                       -> Vector e
                       -> Vector e
                       -> Vector e
mulHermMatrixAddVector a x y =
    runVector $ do
        y' <- newVector_ (dimVector y)
        mulHermMatrixAddToVector a x y y'
        return y'

-- | @mulHermMatrixAddVectorWithScales alpha a x y@
-- returns @alpha * a * x + beta * y@.
mulHermMatrixAddVectorWithScales :: (BLAS2 e)
                                 => e
                                 -> Herm Matrix e
                                 -> Vector e
                                 -> e
                                 -> Vector e
                                 -> Vector e
mulHermMatrixAddVectorWithScales alpha a x beta y =
    runVector $ do
        y' <- newVector_ (dimVector y)
        mulHermMatrixAddToVectorWithScales alpha a x beta y y'
        return y'

-- | @mulHermMatrixMatrix side a b@
-- returns @alpha * a * b@ when @side@ is @LeftSide@ and
-- @alpha * b * a@ when @side@ is @RightSide@.
mulHermMatrixMatrix :: (BLAS3 e)
                    => Side -> Herm Matrix e
                    -> Matrix e
                    -> Matrix e
mulHermMatrixMatrix side a b = 
    runMatrix $ do
        c <- newMatrix_ (dimMatrix b)
        mulHermMatrixToMatrix side a b c
        return c

-- | @mulHermMatrixMatrixWithScale alpha side a b@
-- returns @alpha * a * b@ when @side@ is @LeftSide@ and
-- @alpha * b * a@ when @side@ is @RightSide@.
mulHermMatrixMatrixWithScale :: (BLAS3 e)
                             => e
                             -> Side -> Herm Matrix e
                             -> Matrix e
                             -> Matrix e
mulHermMatrixMatrixWithScale alpha side a b =
    runMatrix $ do
        c <- newMatrix_ (dimMatrix b)
        mulHermMatrixToMatrixWithScale alpha side a b c
        return c

-- | @mulHermMatrixAddMatrix transa a transb b c@
-- returns @a * b + c@ when @side@ is @LeftSide@ and
-- @b * a + c@ when @side@ is @RightSide@.
mulHermMatrixAddMatrix :: (BLAS3 e)
                       => Side -> Herm Matrix e
                       -> Matrix e
                       -> Matrix e
                       -> Matrix e
mulHermMatrixAddMatrix side a b c =
    runMatrix $ do
        c' <- newMatrix_ (dimMatrix c)
        mulHermMatrixAddToMatrix side a b c c'
        return c'

-- | @mulHermMatrixAddMatrixWithScales alpha side a b beta c@
-- returns @alpha * a * b + beta * c@ when @side@ is @LeftSide@ and
-- @alpha * b * a + beta * c@ when @side@ is @RightSide@.
mulHermMatrixAddMatrixWithScales :: (BLAS3 e)
                                 => e
                                 -> Side -> Herm Matrix e
                                 -> Matrix e
                                 -> e
                                 -> Matrix e
                                 -> Matrix e
mulHermMatrixAddMatrixWithScales alpha side a b beta c = 
    runMatrix $ do
        c' <- newMatrix_ (dimMatrix c)
        mulHermMatrixAddToMatrixWithScales alpha side a b beta c c'
        return c'

-- | @mulHermMatrixToVector transa a x y@
-- sets @y := a * x@.
mulHermMatrixToVector :: (RMatrix m, RVector v, BLAS2 e)
                      => Herm m e
                      -> v e
                      -> STVector s e
                      -> ST s ()
mulHermMatrixToVector = mulHermMatrixToVectorWithScale 1

-- | @mulHermMatrixToVectorWithScale alpha transa a x y@
-- sets @y := alpha * a * x@.
mulHermMatrixToVectorWithScale :: (RMatrix m, RVector v, BLAS2 e)
                               => e
                               -> Herm m e
                               -> v e
                               -> STVector s e
                               -> ST s ()
mulHermMatrixToVectorWithScale alpha a x y =
    mulHermMatrixAddToVectorWithScales alpha a x 0 y y

-- | @mulHermMatrixAddToVector a x y y'@
-- sets @y' := a * x + y@.
mulHermMatrixAddToVector :: (RMatrix m, RVector v1, RVector v2, BLAS2 e)
                         => Herm m e
                         -> v1 e
                         -> v2 e
                         -> STVector s e
                         -> ST s ()
mulHermMatrixAddToVector a x y y' =
    mulHermMatrixAddToVectorWithScales 1 a x 1 y y'

-- | @mulHermMatrixAddToVectorWithScales alpha a x beta y y'@
-- sets @y' := alpha * a * x + beta * y@.
mulHermMatrixAddToVectorWithScales :: (RMatrix m, RVector v1, RVector v2, BLAS2 e)
                                   => e
                                   -> Herm m e
                                   -> v1 e
                                   -> e
                                   -> v2 e
                                   -> STVector s e
                                   -> ST s ()
mulHermMatrixAddToVectorWithScales alpha (Herm uplo a) x beta y y'
    | ma /= na = error $
        printf ("mulHermMatrixAddToVectorWithScales _"
                ++ " (Herm %s <matrix with dim (%d,%d)>)"
                ++ " %s <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>: Herm matrix is not square")
               (show uplo) ma na
               nx ny ny'
               
    | (not . and) [ (ma,na) == (n,n)
                  , nx == n
                  , ny == n
                  , ny' == n
                  ] = error $
        printf ("mulHermMatrixAddToVectorWithScales _"
                ++ " (Herm %s <matrix with dim (%d,%d)>)"
                ++ " %s <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>: dimension mismatch")
               (show uplo) ma na
               nx ny ny'

    | otherwise = do
        when (beta /= 0) $ unsafeCopyToVector y y'
        unsafeIOToST $
            unsafeWithMatrix a $ \pa lda ->
            unsafeWithVector x $ \px ->
            unsafeWithVector y' $ \py ->
                BLAS.hemv uplo n alpha pa lda px 1 beta py 1
  where
    (ma,na) = dimMatrix a
    nx = dimVector x
    ny = dimVector y
    ny' = dimVector y'
    n = ny'

-- | @mulHermMatrixToMatrix side a b c@
-- sets @c := a * b@ when @side@ is @LeftSide@ and
-- @c := b * a@ when @side@ is @RightSide@.
mulHermMatrixToMatrix :: (RMatrix m1, RMatrix m2, BLAS3 e)
                      => Side -> Herm m1 e
                      -> m2 e
                      -> STMatrix s e
                      -> ST s ()
mulHermMatrixToMatrix = mulHermMatrixToMatrixWithScale 1

-- | @mulHermMatrixToMatrixWithScale alpha side a b c@
-- sets @c := alpha * a * b@ when @side@ is @LeftSide@ and
-- @c := alpha * b * a@ when @side@ is @RightSide@.
mulHermMatrixToMatrixWithScale :: (RMatrix m1, RMatrix m2, BLAS3 e)
                               => e
                               -> Side -> Herm m1 e
                               -> m2 e
                               -> STMatrix s e
                               -> ST s ()
mulHermMatrixToMatrixWithScale alpha side a b c =
    mulHermMatrixAddToMatrixWithScales alpha side a b 0 c c

-- | @mulHermMatrixAddToMatrix side a b c c'@
-- sets @c' := a * b + c@ when @side@ is @LeftSide@ and
-- @c' := b * a + c@ when @side@ is @RightSide@.
mulHermMatrixAddToMatrix :: (RMatrix m1, RMatrix m2, RMatrix m3, BLAS3 e)
                         => Side -> Herm m1 e
                         -> m2 e
                         -> m3 e
                         -> STMatrix s e
                         -> ST s ()
mulHermMatrixAddToMatrix side a b c c' =
    mulHermMatrixAddToMatrixWithScales 1 side a b 1 c c'

-- | @mulHermMatrixAddToMatrixWithScales alpha side a b beta c c'@
-- sets @c' := alpha * a * b + beta * c@ when @side@ is @LeftSide@ and
-- @c' := alpha * b * a + beta * c@ when @side@ is @RightSide@.
mulHermMatrixAddToMatrixWithScales :: (RMatrix m1, RMatrix m2, RMatrix m3, BLAS3 e)
                                   => e
                                   -> Side -> Herm m1 e
                                   -> m2 e
                                   -> e
                                   -> m3 e
                                   -> STMatrix s e
                                   -> ST s ()
mulHermMatrixAddToMatrixWithScales alpha side (Herm uplo a) b beta c c'
    | ma /= na = error $
        printf ("mulHermMatrixAddToMatrixWithScales _"
                ++ " %s (Herm %s <matrix with dim (%d,%d)>)" 
                ++ " <matrix with dim (%d,%d)>"
                ++ " _"
                ++ " <matrix with dim (%d,%d)>"
                ++ " <matrix with dim (%d,%d)>: Herm matrix is not square")
               (show side) (show uplo) ma na
               mb nb
               mc nc
               mc' nc'
    | (not . and) [ case side of LeftSide  -> (ma,na) == (m,m)
                                 RightSide -> (ma,na) == (n,n)
                  , (mb, nb ) == (m,n)
                  , (mc, nc ) == (m,n)
                  , (mc',nc') == (m,n)
                  ] = error $
        printf ("mulHermMatrixAddToMatrixWithScales _"
                ++ " %s (Herm %s <matrix with dim (%d,%d)>)" 
                ++ " <matrix with dim (%d,%d)>"
                ++ " _"
                ++ " <matrix with dim (%d,%d)>"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
               (show side) (show uplo) ma na
               mb nb
               mc nc
               mc' nc'
    | otherwise = do
        when (beta /= 0) $ unsafeCopyToMatrix c c'
        unsafeIOToST $
            unsafeWithMatrix a $ \pa lda ->
            unsafeWithMatrix b $ \pb ldb ->
            unsafeWithMatrix c' $ \pc ldc ->
                BLAS.hemm side uplo m n alpha pa lda pb ldb beta pc ldc
  where
    (ma,na) = dimMatrix a
    (mb,nb) = dimMatrix b
    (mc,nc) = dimMatrix c
    (mc',nc') = dimMatrix c'
    (m,n) = dimMatrix c'
