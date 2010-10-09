{-# LANGUAGE Rank2Types #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Immutable dense matrices.

module Numeric.LinearAlgebra.Matrix (
    -- * Immutable matrices
    Matrix,
    dim,

    -- * Matrix construction
    fromList,
    fromCols,
    fromRows,
    constant,
    zero,

    -- * Accessing matrices
    at,
    indices,
    elems,
    assocs,

    -- * Incremental matrix updates
    update,
    unsafeUpdate,
    accum,
    unsafeAccum,

    -- * Derived matrices
    map,
    zipWith,

    -- * Matrix views
    slice,
    takeRows,
    dropRows,
    splitRowsAt,
    takeCols,
    dropCols,
    splitColsAt,

    -- * Matrix rows and columns
    col,
    cols,
    row,
    rows,
    
    -- * Matrix diagonals
    diag,
    
    -- * Conversions to vectors
    toVector,
    
    -- * Conversions from vectors
    fromVector,
    fromCol,
    fromRow,

    -- * Matrix math operations
    shiftDiagBy,
    shiftDiagByWithScale,    
    add,
    addWithScale,
    sub,
    scaleBy,
    scaleRowsBy,
    scaleColsBy,
    negate,
    conjugate,

    -- * Linear algebra
    trans,
    conjTrans,
    rank1Update,
    
    -- ** Matrix-Vector multiplication
    mulVector,
    mulVectorWithScale,
    addMulVectorWithScales,
    
    -- ** Matrix-Matrix multiplication
    mulMatrix,
    mulMatrixWithScale,
    addMulMatrixWithScales,

    -- * Conversions between foreign pointers
    unsafeFromForeignPtr,
    unsafeToForeignPtr,

    -- * Mutable interface
    module Numeric.LinearAlgebra.Matrix.ST,
    
    -- * Hermitian views
    module Numeric.LinearAlgebra.Matrix.Herm,
    
    -- * Cholesky factorizations
    module Numeric.LinearAlgebra.Matrix.Cholesky,
    
    -- * Eigenvalues and eigenvectors
    module Numeric.LinearAlgebra.Matrix.Eigen,

    -- * Basic multivariate statistics
    module Numeric.LinearAlgebra.Matrix.Statistics,

    ) where

import Prelude hiding ( map, zipWith, negate, )
import Control.Monad( zipWithM_ )
import Control.Monad.ST( ST )
import Foreign( Storable )
import Text.Printf( printf )

import Foreign.BLAS( BLAS1, BLAS2, BLAS3, Trans(..) )
import Foreign.VMath( VNum )


import Numeric.LinearAlgebra.Matrix.Base
import Numeric.LinearAlgebra.Matrix.Herm
import Numeric.LinearAlgebra.Matrix.ST
import Numeric.LinearAlgebra.Matrix.Cholesky
import Numeric.LinearAlgebra.Matrix.Eigen
import Numeric.LinearAlgebra.Matrix.Statistics
import Numeric.LinearAlgebra.Vector( Vector )
import qualified Numeric.LinearAlgebra.Vector as V


infixl 7 `scaleBy`, `scaleRowsBy`, `scaleColsBy`
infixl 6 `add`, `shiftDiagBy`, `sub`


-- | Create a matrix of the given dimension with the given vectors as
-- columns.
fromCols :: (Storable e) => (Int,Int) -> [Vector e] -> Matrix e
fromCols mn cs = create $ do
    a <- new_ mn
    withColsM a $ \cs' -> zipWithM_ V.copyTo cs' cs
    return a

-- | Create a matrix of the given dimension with the given vectors as
-- rows.
fromRows :: (Storable e) => (Int,Int) -> [Vector e] -> Matrix e
fromRows (m,n) rs = create $ do
    a <- new_ (m,n)
    sequence_ [ setRow a i r | (i,r) <- zip [ 0..m-1 ] rs ]
    return a

-- | Get the given row of the matrix.
row :: (BLAS1 e) => Matrix e -> Int -> Vector e
row a i
    | i < 0 || i >= m = error $
        printf ("row <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n i
    | otherwise =
        unsafeRow a i
  where
    (m,n) = dim a
{-# INLINE row #-}

-- | Version of 'row' that doesn't range-check indices.
unsafeRow :: (BLAS1 e) => Matrix e -> Int -> Vector e
unsafeRow a@(Matrix v _ n lda) i
    | lda == 1 = v
    | otherwise = V.create $ do
        r <- V.new_ n
        unsafeRowTo r a i
        return r

-- | Get a list of the rows in a matrix
rows :: (BLAS1 e) => Matrix e -> [Vector e]
rows a = [ unsafeRow a i | i <- [ 0..m-1 ] ]
  where
    (m,_) = dim a


-- | Get the diagonal of the matrix.
diag :: (Storable e) => Matrix e -> Vector e
diag a = V.create $ do
    x <- V.new_ mn
    diagTo x a
    return x
  where
    (m,n) = dim a
    mn = min m n

-- | @shiftDiagBy d a@ returns @diag(d) + a@.
shiftDiagBy :: (BLAS1 e) => Vector e -> Matrix e -> Matrix e
shiftDiagBy s a = create $ do
    a' <- newCopy a
    shiftDiagByM_ s a'
    return a'

-- | @shiftDiagByWithScale alpha d a@ returns @alpha * diag(d) + a@.
shiftDiagByWithScale :: (BLAS1 e) => e -> Vector e -> Matrix e -> Matrix e
shiftDiagByWithScale e s a = create $ do
    a' <- newCopy a
    shiftDiagByWithScaleM_ e s a'
    return a'

-- | @add a b@ returns @a + b@.
add :: (VNum e) => Matrix e -> Matrix e -> Matrix e
add = result2 addTo

-- | @sub a b@ returns @a - b@.
sub :: (VNum e) => Matrix e -> Matrix e -> Matrix e
sub = result2 subTo

-- | @scaleBy k a@ returns @k * a@.
scaleBy :: (BLAS1 e) => e -> Matrix e -> Matrix e
scaleBy k a = create $ do
    a' <- newCopy a
    scaleByM_ k a'
    return a'

-- | @addWithScale alpha x y@ returns @alpha * x + y@.
addWithScale :: (BLAS1 e) => e -> Matrix e -> Matrix e -> Matrix e
addWithScale alpha x y = create $ do
    y' <- newCopy y
    addWithScaleM_ alpha x y'
    return y'

-- | @scaleRowsBy s a@ returns @diag(s) * a@.
scaleRowsBy :: (BLAS1 e) => Vector e -> Matrix e -> Matrix e
scaleRowsBy s a = create $ do
    a' <- newCopy a
    scaleRowsByM_ s a'
    return a'

-- | @scaleColsBy s a@ returns @a * diag(s)@.
scaleColsBy :: (BLAS1 e) => Vector e -> Matrix e -> Matrix e
scaleColsBy s a = create $ do
    a' <- newCopy a
    scaleColsByM_ s a'
    return a'

-- | @negate a@ returns @-a@.
negate :: (VNum e) => Matrix e -> Matrix e
negate = result negateTo

-- | @conjugate a@ returns @conjugate(a)@.
conjugate :: (VNum e) => Matrix e -> Matrix e
conjugate = result conjugateTo

-- | @trans a@ retunrs @trans(a)@.
trans :: (BLAS1 e)
      => Matrix e
      -> Matrix e
trans a = let
    (m,n) = dim a
    in create $ do
        a' <- new_ (n,m)
        transTo a' a
        return a'

-- | @conjTrans a@ retunrs @conj(trans(a))@.
conjTrans :: (BLAS1 e)
          => Matrix e
          -> Matrix e
conjTrans a = let
    (m,n) = dim a
    in create $ do
        a' <- new_ (n,m)
        conjTransTo a' a
        return a'

-- | @rank1Update alpha x y a@ returns @alpha * x * y^H + a@.
rank1Update :: (BLAS2 e)
            => e
            -> Vector e
            -> Vector e
            -> Matrix e
            -> Matrix e
rank1Update alpha x y a =
    create $ do
        a' <- newCopy a
        rank1UpdateM_ alpha x y a'
        return a'

-- | @mulVector transa a x@
-- returns @op(a) * x@, where @op(a)@ is determined by @transa@.                   
mulVector :: (BLAS2 e)
          => Trans -> Matrix e
          -> Vector e
          -> Vector e
mulVector transa a x = let
    m = case transa of NoTrans -> (fst . dim) a
                       _       -> (snd . dim) a
    in V.create $ do
        y <- V.new_ m
        mulVectorTo y transa a x
        return y

-- | @mulVectorWithScale alpha transa a x@
-- retunrs @alpha * op(a) * x@, where @op(a)@ is determined by @transa@.                   
mulVectorWithScale :: (BLAS2 e)
                   => e
                   -> Trans -> Matrix e
                   -> Vector e
                   -> Vector e
mulVectorWithScale alpha transa a x = let
    m = case transa of NoTrans -> (fst . dim) a
                       _       -> (snd . dim) a
    in V.create $ do
        y <- V.new_ m
        mulVectorWithScaleTo y alpha transa a x
        return y
                       
-- | @addMulVectorWithScales alpha transa a x beta y@
-- returns @alpha * op(a) * x + beta * y@, where @op(a)@ is
-- determined by @transa@.
addMulVectorWithScales :: (BLAS2 e)
                       => e
                       -> Trans -> Matrix e
                       -> Vector e
                       -> e
                       -> Vector e
                       -> Vector e
addMulVectorWithScales alpha transa a x beta y =
    V.create $ do
        y' <- V.newCopy y
        addMulVectorWithScalesM_ alpha transa a x beta y'
        return y'

-- | @mulMatrix transa a transb b@
-- returns @op(a) * op(b)@, where @op(a)@ and @op(b)@ are determined
-- by @transa@ and @transb@.                   
mulMatrix :: (BLAS3 e)
          => Trans -> Matrix e
          -> Trans -> Matrix e
          -> Matrix e
mulMatrix transa a transb b = let
    m = case transa of NoTrans -> (fst . dim) a
                       _       -> (snd . dim) a
    n = case transb of NoTrans -> (snd . dim) b
                       _       -> (fst . dim) b
    in create $ do
        c <- new_ (m,n)
        mulMatrixTo c transa a transb b
        return c

-- | @mulMatrixWithScale alpha transa a transb b@
-- returns @alpha * op(a) * op(b)@, where @op(a)@ and @op(b)@ are determined
-- by @transa@ and @transb@.                   
mulMatrixWithScale :: (BLAS3 e)
                   => e
                   -> Trans -> Matrix e
                   -> Trans -> Matrix e
                   -> Matrix e
mulMatrixWithScale alpha transa a transb b = let
    m = case transa of NoTrans -> (fst . dim) a
                       _       -> (snd . dim) a
    n = case transb of NoTrans -> (snd . dim) b
                       _       -> (fst . dim) b
    in create $ do
        c <- new_ (m,n)
        mulMatrixWithScaleTo c alpha transa a transb b
        return c

-- | @addMulMatrixWithScales alpha transa a transb b beta c@
-- returns @alpha * op(a) * op(b) + beta * c@, where @op(a)@ and
-- @op(b)@ are determined by @transa@ and @transb@.
addMulMatrixWithScales :: (BLAS3 e)
                       => e
                       -> Trans -> Matrix e
                       -> Trans -> Matrix e
                       -> e
                       -> Matrix e
                       -> Matrix e
addMulMatrixWithScales alpha transa a transb b beta c = 
    create $ do
        c' <- newCopy c
        addMulMatrixWithScalesM_ alpha transa a transb b beta c'
        return c'

result :: (Storable e, Storable f)
       => (forall s . STMatrix s f -> Matrix e -> ST s a)
       -> Matrix e
       -> Matrix f
result f a = create $ newResult f a
{-# INLINE result #-}

result2 :: (Storable e, Storable f, Storable g)
        => (forall s . STMatrix s g -> Matrix e -> Matrix f -> ST s a)
        -> Matrix e
        -> Matrix f
        -> Matrix g
result2 f a1 a2 = create $ newResult2 f a1 a2
{-# INLINE result2 #-}

newResult :: (RMatrix m, Storable e, Storable f)
          => (STMatrix s f -> m e -> ST s a)
          -> m e
          -> ST s (STMatrix s f)
newResult f a = do
    mn <- getDim a
    c <- new_ mn
    _ <- f c a
    return c
{-# INLINE newResult #-}

newResult2 :: (RMatrix m1, RMatrix m2, Storable e, Storable f, Storable g)
           => (STMatrix s g -> m1 e -> m2 f -> ST s a)
           -> m1 e
           -> m2 f
           -> ST s (STMatrix s g)
newResult2 f a1 a2 = do
    mn <- getDim a1
    c <- new_ mn
    _ <- f c a1 a2
    return c
{-# INLINE newResult2 #-}
