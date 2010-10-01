{-# LANGUAGE DeriveDataTypeable, FlexibleContexts, Rank2Types #-}
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
    fromVector,
    unsafeFromVector,
    toVector,
    viewFromSTVector,
    unsafeViewFromSTVector,
    toSTVectorView,
    withSTVectorView,
    
    -- * Conversions between mutable and immutable packed matrices
    create,
    freeze,
    unsafeFreeze,
    
    -- * Immutable interface

    -- ** Herm Packed-Vector multiplication
    hermMulVector,
    hermMulVectorWithScale,
    hermMulAddVectorWithScales,
    
    -- ** Herm Packed updates
    hermRank1Update,
    hermRank2Update,

    -- * Mutable interface
    newCopy,
    hermCreate,

    -- ** Herm Packed-Vector multiplication
    hermMulToVector,
    hermMulToVectorWithScale,
    hermMulAddToVectorWithScales,
    
    -- ** Herm Packed updates
    hermRank1UpdateTo,
    hermRank2UpdateTo,
    
    ) where

import Control.Monad.ST( ST, RealWorld, runST, unsafeIOToST )
import Data.Typeable( Typeable )
import Foreign( Storable, Ptr )
import Text.Printf( printf )

import Numeric.LinearAlgebra.Matrix.Herm( Herm(..) )
import Numeric.LinearAlgebra.Vector( Vector )
import qualified Numeric.LinearAlgebra.Vector as V
import Numeric.LinearAlgebra.Vector.ST( STVector, RVector )
import qualified Numeric.LinearAlgebra.Vector.ST as V
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

-- | Create a packed matrix view of a vector, ensurint that the
-- vector has dimension @n * (n+1)/2@, where @n@ is the desired dimension.
fromVector :: (Storable e) => Int -> Vector e -> Packed e
fromVector n x
    | not $ 2 * nx == n * (n+1) = error $
        printf ("fromVector %d <vector with dim %d>: dimension mismatch")
               n nx
    | otherwise =
        unsafeFromVector n x
  where
    nx = V.dim x
{-# INLINE fromVector #-}

-- | Create a packed matrix view of a vector, wihtout checking
-- the dimension of the vector.
unsafeFromVector :: (Storable e) => Int -> Vector e -> Packed e
unsafeFromVector = Packed
{-# INLINE unsafeFromVector #-}

-- | Returns the dimension and underlying vector storage of a
-- packed matrix.
toVector :: (Storable e) => Packed e -> (Int, Vector e)
toVector (Packed n v) = (n,v)
{-# INLINE toVector #-}

-- | Create a packed matrix view of a vector, ensurint that the
-- vector has dimension @n * (n+1)/2@, where @n@ is the desired dimension.
viewFromSTVector :: (Storable e) => Int -> STVector s e -> STPacked s e
viewFromSTVector n x
    | not $ 2 * nx == n * (n+1) = error $
        printf ("fromVectorST %d <vector with dim %d>: dimension mismatch")
               n nx
    | otherwise =
        STPacked n x
  where
    nx = V.dim x
{-# INLINE viewFromSTVector #-}

-- | Create a packed matrix view of a vector, wihtout checking
-- the dimension of the vector.
unsafeViewFromSTVector :: (Storable e) => Int -> STVector s e -> STPacked s e
unsafeViewFromSTVector = STPacked
{-# INLINE unsafeViewFromSTVector #-}

-- | Returns the dimension and underlying vector storage of a
-- packed matrix.
toSTVectorView :: (Storable e) => STPacked s e -> (Int, STVector s e)
toSTVectorView (STPacked n v) = (n,v)
{-# INLINE toSTVectorView #-}


-- | Read-only packed matrices.
class RPacked p where
    -- | Returns the dimension of the packed matrix.
    dim :: (Storable e) => p e -> Int

    -- | Perform an action with the underlying vector storage of
    -- the packed matrix.
    withVectorView :: (Storable e)
                   => p e
                   -> (forall v . RVector v => v e -> ST s a)
                   -> ST s a
    
    -- | Perform an IO action with a pointer to the first element of
    -- the packed matrix.
    unsafeWith :: (Storable e) => p e -> (Ptr e -> IO a) -> IO a

-- | Perform an action with the underlying vector storage of
-- the mutable packed matrix.
withSTVectorView :: (Storable e)
             => STPacked s e
             -> (STVector s e -> ST s a)
             -> ST s a
withSTVectorView (STPacked _ v) f = f v
{-# INLINE withSTVectorView #-}

instance RPacked Packed where
    dim (Packed n _) = n
    {-# INLINE dim #-}
    withVectorView (Packed _ v) f = f v
    {-# INLINE withVectorView #-}
    unsafeWith (Packed _ v) = V.unsafeWith v
    {-# INLINE unsafeWith #-}

instance RPacked (STPacked s) where
    dim (STPacked n _) = n
    {-# INLINE dim #-}
    withVectorView (STPacked _ v) f = f v
    {-# INLINE withVectorView #-}
    unsafeWith (STPacked _ v) = V.unsafeWith v
    {-# INLINE unsafeWith #-}


-- | Create a new copy of a packed matrix.
newCopy :: (Storable e, RPacked p)
              => p e -> ST s (STPacked s e)
newCopy p =
    withVectorView p $ \x ->
        unsafeViewFromSTVector (dim p) `fmap` V.newCopy x
{-# INLINE newCopy #-}

-- | Converts a mutable packed matrix to an immutable one by taking a complete
-- copy of it.
freeze :: (Storable e) => STPacked s e -> ST s (Packed e)
freeze (STPacked n mp) = do
    p <- V.freeze mp
    return $ Packed n p

-- | Converts a mutable matrix into an immutable matrix. This simply casts
-- the matrix from one type to the other without copying the matrix. Note
-- that because the matrix is possibly not copied, any subsequent
-- modifications made to the mutable version of the matrix may be shared with
-- the immutable version. It is safe to use, therefore, if the mutable
-- version is never modified after the freeze operation.
unsafeFreeze :: (Storable e) => STPacked s e -> ST s (Packed e)
unsafeFreeze (STPacked n mp) = do
    p <- V.unsafeFreeze mp
    return $ Packed n p

-- | A safe way to V.create and work with a mutable Packed before returning 
-- an immutable one for later perusal.
create :: (Storable e)
       => (forall s. ST s ((STPacked s) e))
       -> Packed e
create stmp = runST $ do
    mp <- stmp
    unsafeFreeze mp


-- | A safe way to V.create and work with a mutable Herm Packed before returning 
-- an immutable one for later perusal.
hermCreate :: (Storable e)
           => (forall s. ST s (Herm (STPacked s) e))
           -> Herm Packed e
hermCreate stmh = runST $ do
    (Herm u mp) <- stmh
    p <- unsafeFreeze mp
    return $ Herm u p

-- | @hermRank1Update alpha x a@ returns
-- @alpha * x * x^H + a@.
hermRank1Update :: (BLAS2 e)
                => Double -> Vector e -> Herm Packed e -> Herm Packed e
hermRank1Update alpha x (Herm uplo ap) = hermCreate $ do
    hp' <- Herm uplo `fmap` newCopy ap
    hermRank1UpdateTo alpha x hp'
    return hp'

-- | @hermRank2Update alpha x y a@ returns
-- @alpha * x * y^H + conj(alpha) * y * x^H + a@.
hermRank2Update :: (BLAS2 e)
                => e -> Vector e -> Vector e -> Herm Packed e
                -> Herm Packed e
hermRank2Update alpha x y (Herm uplo ap) = hermCreate $ do
    hp' <- Herm uplo `fmap` newCopy ap
    hermRank2UpdateTo alpha x y hp'
    return hp'

-- | @hermRank1UpdateTo alpha x a@ sets
-- @a := alpha * x * x^H + a@.
hermRank1UpdateTo :: (RVector v, BLAS2 e)
                  => Double -> v e -> Herm (STPacked s) e -> ST s ()
hermRank1UpdateTo alpha x (Herm uplo a)
    | (not . and) [ nx == n, na == n ] = error $
        printf ("hermRank1UpdateTo _ <vector with dim %d>"
                 ++ " (Herm _ <packed matrix with dim %d>):"
                 ++ " invalid dimensions") nx na
    | otherwise =
        unsafeIOToST $
        V.unsafeWith x $ \px ->
        unsafeWith a $ \pa ->
            BLAS.hpr uplo n alpha px 1 pa
  where
    nx = V.dim x
    na = dim a
    n = nx

-- | @hermRank2UpdateTo alpha x y a@ sets
-- @a := alpha * x * y^H + conj(alpha) * y * x^H + a@.
hermRank2UpdateTo :: (RVector v1, RVector v2, BLAS2 e)
                  => e -> v1 e -> v2 e -> Herm (STPacked s) e -> ST s ()
hermRank2UpdateTo alpha x y (Herm uplo a)
    | (not . and) [ nx == n, ny == n, na == n ] = error $
        printf ("hermRank2UpdateTo _ <vector with dim %d>"
                 ++ " <vector with dim %d>"
                 ++ " (Herm _ <packed matrix with dim %d>):"
                 ++ " invalid dimensions") nx ny na
    | otherwise =
        unsafeIOToST $
        V.unsafeWith x $ \px ->
        V.unsafeWith y $ \py ->        
        unsafeWith a $ \pa ->
            BLAS.hpr2 uplo n alpha px 1 py 1 pa
  where
    nx = V.dim x
    ny = V.dim y
    na = dim a
    n = nx

-- | @hermMulVector a x@ returns @a * x@.
hermMulVector :: (BLAS2 e)
                    => Herm Packed e
                    -> Vector e
                    -> Vector e
hermMulVector a x =
    V.create $ do
        y <- V.new_ (V.dim x)
        hermMulToVector a x y
        return y

-- | @hermMulVectorWithScale alpha a x@ retunrs @alpha * a * x@.
hermMulVectorWithScale :: (BLAS2 e)
                       => e
                       -> Herm Packed e
                       -> Vector e
                       -> Vector e
hermMulVectorWithScale alpha a x =
    V.create $ do
        y <- V.new_ (V.dim x)
        hermMulToVectorWithScale alpha a x y
        return y
                       
-- | @hermMulAddVectorWithScales alpha a x y@
-- returns @alpha * a * x + beta * y@.
hermMulAddVectorWithScales :: (BLAS2 e)
                           => e
                           -> Herm Packed e
                           -> Vector e
                           -> e
                           -> Vector e
                           -> Vector e
hermMulAddVectorWithScales alpha a x beta y =
    V.create $ do
        y' <- V.newCopy y
        hermMulAddToVectorWithScales alpha a x beta y'
        return y'

-- | @hermMulToVector a x y@ sets @y := a * x@.
hermMulToVector :: (RPacked p, RVector v, BLAS2 e)
                => Herm p e
                -> v e
                -> STVector s e
                -> ST s ()
hermMulToVector = hermMulToVectorWithScale 1

-- | @hermMulToVectorWithScale alpha a x y@
-- sets @y := alpha * a * x@.
hermMulToVectorWithScale :: (RPacked p, RVector v, BLAS2 e)
                         => e
                         -> Herm p e
                         -> v e
                         -> STVector s e
                         -> ST s ()
hermMulToVectorWithScale alpha a x y =
    hermMulAddToVectorWithScales alpha a x 0 y

-- | @hermMulAddToVectorWithScales alpha a x beta y@
-- sets @y := alpha * a * x + beta * y@.
hermMulAddToVectorWithScales :: (RPacked p, RVector v, BLAS2 e)
                             => e
                             -> Herm p e
                             -> v e
                             -> e
                             -> STVector s e
                             -> ST s ()
hermMulAddToVectorWithScales alpha (Herm uplo a) x beta y
    | (not . and) [ na == n
                  , nx == n
                  , ny == n
                  ] = error $
        printf ("hermMulAddToVectorWithScales _"
                ++ " (Herm %s <packed matrix with dim %d>)"
                ++ " %s <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>: dimension mismatch")
               (show uplo) na
               nx ny

    | otherwise =
        unsafeIOToST $
            unsafeWith a $ \pa ->
            V.unsafeWith x $ \px ->
            V.unsafeWith y $ \py ->
                BLAS.hpmv uplo n alpha pa px 1 beta py 1
  where
    na = dim a
    nx = V.dim x
    ny = V.dim y
    n = ny
