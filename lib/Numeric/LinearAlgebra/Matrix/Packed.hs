{-# LANGUAGE DeriveDataTypeable, FlexibleContexts, TypeFamilies #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix.Packed
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Packed matrices.
--
module Numeric.LinearAlgebra.Matrix.Packed (
    -- * Packed matrix datatypes
    Packed,
    STPacked,
    IOPacked,

    RPacked(..),
    toPacked,
    
    -- * Conversions between mutable and immutable packed matrices
    runPacked,
    freezePacked,
    unsafeFreezePacked,
    
    -- * Immutable interface

    -- ** Herm Packed-Vector multiplication
    mulHermPackedVector,
    mulHermPackedVectorWithScale,
    mulHermPackedAddVectorWithScales,
    
    -- ** Herm Packed updates
    rank1UpdateHermPacked,
    rank2UpdateHermPacked,

    -- * Mutable interface
    newCopyPacked,
    runHermPacked,

    -- ** Herm Packed-Vector multiplication
    mulHermPackedToVector,
    mulHermPackedToVectorWithScale,
    mulHermPackedAddToVectorWithScales,
    
    -- ** Herm Packed updates
    rank1UpdateToHermPacked,
    rank2UpdateToHermPacked,
    
    ) where

import Control.Monad.ST( ST, RealWorld, runST, unsafeIOToST )
import Data.Typeable( Typeable )
import Foreign( Storable, Ptr )
import Text.Printf( printf )

import Numeric.LinearAlgebra.Types( HasVectorView(..) )
import Numeric.LinearAlgebra.Matrix.Herm( Herm(..) )
import Numeric.LinearAlgebra.Vector( Vector, dimVector   )
import Numeric.LinearAlgebra.Vector.ST( STVector, RVector, freezeVector,
    newVector_, newCopyVector, runVector, unsafeFreezeVector,
    unsafeWithVector )
import Foreign.BLAS( BLAS2 )
import qualified Foreign.BLAS as BLAS


-- | Immutable packed matrices, stored in column-major order.
data Packed e = Packed !Int !(Vector e)
    deriving (Typeable)
    
-- | Mutable packed matrices in the 'ST' monad.
data STPacked s e = STPacked !Int !(STVector s e)
    deriving (Typeable)
    
-- | Mutable packed matrices in the 'IO' monad.
type IOPacked = STPacked RealWorld


-- | Read-only packed matrices.
class (HasVectorView p, RVector (VectorView p)) => RPacked p where
    -- | Create a packed matrix view of a vector, wihtout checking
    -- the dimension of the vector.
    unsafeToPacked :: (Storable e) => Int -> VectorView p e -> p e
    
    -- | Returns the dimension and underlying vector storage of a
    -- packed matrix.
    fromPacked :: (Storable e) => p e -> (Int, VectorView p e)
    
    -- | Returns the dimension of the packed matrix.
    dimPacked :: (Storable e) => p e -> Int

    -- | Returns the underlying storage of the packed matrix.
    unPacked :: (Storable e) => p e -> VectorView p e
    
    -- | Perform an IO action with a pointer to the first element of
    -- the packed matrix.
    unsafeWithPacked :: (Storable e) => p e -> (Ptr e -> IO a) -> IO a

instance HasVectorView Packed where
    type VectorView Packed = Vector

instance HasVectorView (STPacked s) where
    type VectorView (STPacked s) = STVector s

instance RPacked Packed where
    unsafeToPacked = Packed
    {-# INLINE unsafeToPacked #-}
    fromPacked (Packed n v) = (n,v)
    {-# INLINE fromPacked #-}
    dimPacked (Packed n _) = n
    {-# INLINE dimPacked #-}
    unPacked (Packed _ v) = v
    {-# INLINE unPacked #-}
    unsafeWithPacked (Packed _ v) = unsafeWithVector v
    {-# INLINE unsafeWithPacked #-}

instance RPacked (STPacked s) where
    unsafeToPacked = STPacked
    {-# INLINE unsafeToPacked #-}
    fromPacked (STPacked n v) = (n,v)
    {-# INLINE fromPacked #-}
    dimPacked (STPacked n _) = n
    {-# INLINE dimPacked #-}
    unPacked (STPacked _ v) = v
    {-# INLINE unPacked #-}
    unsafeWithPacked (STPacked _ v) = unsafeWithVector v
    {-# INLINE unsafeWithPacked #-}

-- | Create a packed matrix view of a vector, ensurint that the
-- vector has dimension @n * (n+1)/2@, where @n@ is the desired dimension.
toPacked :: (Storable e, RPacked p) => Int -> VectorView p e -> p e
toPacked n x
    | not $ 2 * nx == n * (n+1) = error $
        printf ("toPacked %d <vector with dim %d>: dimension mismatch")
               n nx
    | otherwise =
        unsafeToPacked n x
  where
    nx = dimVector x


-- | Create a new copy of a packed matrix.
newCopyPacked :: (Storable e, RPacked p)
              => p e -> ST s (STPacked s e)
newCopyPacked p = let
    (n,v) = fromPacked p
    in unsafeToPacked n `fmap` newCopyVector v
{-# INLINE newCopyPacked #-}

-- | Converts a mutable packed matrix to an immutable one by taking a complete
-- copy of it.
freezePacked :: (Storable e) => STPacked s e -> ST s (Packed e)
freezePacked (STPacked n mp) = do
    p <- freezeVector mp
    return $ Packed n p

-- | Converts a mutable matrix into an immutable matrix. This simply casts
-- the matrix from one type to the other without copying the matrix. Note
-- that because the matrix is possibly not copied, any subsequent
-- modifications made to the mutable version of the matrix may be shared with
-- the immutable version. It is safe to use, therefore, if the mutable
-- version is never modified after the freeze operation.
unsafeFreezePacked :: (Storable e) => STPacked s e -> ST s (Packed e)
unsafeFreezePacked (STPacked n mp) = do
    p <- unsafeFreezeVector mp
    return $ Packed n p

-- | A safe way to create and work with a mutable Packed before returning 
-- an immutable one for later perusal.
runPacked :: (Storable e)
          => (forall s. ST s ((STPacked s) e))
          -> Packed e
runPacked stmp = runST $ do
    mp <- stmp
    unsafeFreezePacked mp


-- | A safe way to create and work with a mutable Herm Packed before returning 
-- an immutable one for later perusal.
runHermPacked :: (Storable e)
              => (forall s. ST s (Herm (STPacked s) e))
              -> Herm Packed e
runHermPacked stmh = runST $ do
    (Herm u mp) <- stmh
    p <- unsafeFreezePacked mp
    return $ Herm u p

-- | @rank1UpdateHermPacked alpha x a@ returns
-- @alpha * x * x^H + a@.
rank1UpdateHermPacked :: (BLAS2 e)
                      => Double -> Vector e -> Herm Packed e -> Herm Packed e
rank1UpdateHermPacked alpha x (Herm uplo ap) = runHermPacked $ do
    hp' <- Herm uplo `fmap` newCopyPacked ap
    rank1UpdateToHermPacked alpha x hp'
    return hp'

-- | @rank2UpdateHermPacked alpha x y a@ returns
-- @alpha * x * y^H + conj(alpha) * y * x^H + a@.
rank2UpdateHermPacked :: (BLAS2 e)
                      => e -> Vector e -> Vector e -> Herm Packed e
                      -> Herm Packed e
rank2UpdateHermPacked alpha x y (Herm uplo ap) = runHermPacked $ do
    hp' <- Herm uplo `fmap` newCopyPacked ap
    rank2UpdateToHermPacked alpha x y hp'
    return hp'

-- | @rank1UpdateToHermPacked alpha x a@ sets
-- @a := alpha * x * x^H + a@.
rank1UpdateToHermPacked :: (RVector v, BLAS2 e)
                        => Double -> v e -> Herm (STPacked s) e -> ST s ()
rank1UpdateToHermPacked alpha x (Herm uplo a)
    | (not . and) [ nx == n, na == n ] = error $
        printf ("rank1UpdateToHermPacked _ <vector with dim %d>"
                 ++ " (Herm _ <packed matrix with dim %d>):"
                 ++ " invalid dimensions") nx na
    | otherwise =
        unsafeIOToST $
        unsafeWithVector x $ \px ->
        unsafeWithPacked a $ \pa ->
            BLAS.hpr uplo n alpha px 1 pa
  where
    nx = dimVector x
    na = dimPacked a
    n = nx

-- | @rank2UpdateToHermPacked alpha x y a@ sets
-- @a := alpha * x * y^H + conj(alpha) * y * x^H + a@.
rank2UpdateToHermPacked :: (RVector v1, RVector v2, BLAS2 e)
                        => e -> v1 e -> v2 e -> Herm (STPacked s) e -> ST s ()
rank2UpdateToHermPacked alpha x y (Herm uplo a)
    | (not . and) [ nx == n, ny == n, na == n ] = error $
        printf ("rank2UpdateToHermPacked _ <vector with dim %d>"
                 ++ " <vector with dim %d>"
                 ++ " (Herm _ <packed matrix with dim %d>):"
                 ++ " invalid dimensions") nx ny na
    | otherwise =
        unsafeIOToST $
        unsafeWithVector x $ \px ->
        unsafeWithVector y $ \py ->        
        unsafeWithPacked a $ \pa ->
            BLAS.hpr2 uplo n alpha px 1 py 1 pa
  where
    nx = dimVector x
    ny = dimVector y
    na = dimPacked a
    n = nx

-- | @mulHermPackedVector a x@ returns @a * x@.
mulHermPackedVector :: (BLAS2 e)
                    => Herm Packed e
                    -> Vector e
                    -> Vector e
mulHermPackedVector a x =
    runVector $ do
        y <- newVector_ (dimVector x)
        mulHermPackedToVector a x y
        return y

-- | @mulHermPackedVectorWithScale alpha a x@ retunrs @alpha * a * x@.
mulHermPackedVectorWithScale :: (BLAS2 e)
                             => e
                             -> Herm Packed e
                             -> Vector e
                             -> Vector e
mulHermPackedVectorWithScale alpha a x =
    runVector $ do
        y <- newVector_ (dimVector x)
        mulHermPackedToVectorWithScale alpha a x y
        return y
                       
-- | @mulHermPackedAddVectorWithScales alpha a x y@
-- returns @alpha * a * x + beta * y@.
mulHermPackedAddVectorWithScales :: (BLAS2 e)
                                 => e
                                 -> Herm Packed e
                                 -> Vector e
                                 -> e
                                 -> Vector e
                                 -> Vector e
mulHermPackedAddVectorWithScales alpha a x beta y =
    runVector $ do
        y' <- newCopyVector y
        mulHermPackedAddToVectorWithScales alpha a x beta y'
        return y'

-- | @mulHermPackedToVector a x y@ sets @y := a * x@.
mulHermPackedToVector :: (RPacked p, RVector v, BLAS2 e)
                      => Herm p e
                      -> v e
                      -> STVector s e
                      -> ST s ()
mulHermPackedToVector = mulHermPackedToVectorWithScale 1

-- | @mulHermPackedToVectorWithScale alpha a x y@
-- sets @y := alpha * a * x@.
mulHermPackedToVectorWithScale :: (RPacked p, RVector v, BLAS2 e)
                               => e
                               -> Herm p e
                               -> v e
                               -> STVector s e
                               -> ST s ()
mulHermPackedToVectorWithScale alpha a x y =
    mulHermPackedAddToVectorWithScales alpha a x 0 y

-- | @mulHermPackedAddToVectorWithScales alpha a x beta y@
-- sets @y := alpha * a * x + beta * y@.
mulHermPackedAddToVectorWithScales :: (RPacked p, RVector v, BLAS2 e)
                                   => e
                                   -> Herm p e
                                   -> v e
                                   -> e
                                   -> STVector s e
                                   -> ST s ()
mulHermPackedAddToVectorWithScales alpha (Herm uplo a) x beta y
    | (not . and) [ na == n
                  , nx == n
                  , ny == n
                  ] = error $
        printf ("mulHermPackedAddToVectorWithScales _"
                ++ " (Herm %s <packed matrix with dim %d>)"
                ++ " %s <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>: dimension mismatch")
               (show uplo) na
               nx ny

    | otherwise =
        unsafeIOToST $
            unsafeWithPacked a $ \pa ->
            unsafeWithVector x $ \px ->
            unsafeWithVector y $ \py ->
                BLAS.hpmv uplo n alpha pa px 1 beta py 1
  where
    na = dimPacked a
    nx = dimVector x
    ny = dimVector y
    n = ny
