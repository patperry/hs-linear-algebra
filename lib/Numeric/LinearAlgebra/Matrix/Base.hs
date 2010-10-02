{-# LANGUAGE DeriveDataTypeable, GeneralizedNewtypeDeriving, Rank2Types #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix.Base
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Numeric.LinearAlgebra.Matrix.Base
    where

import Prelude hiding ( read, map, zipWith )
import qualified Prelude as P

import Control.Monad( forM_, zipWithM_ )
import Control.Monad.ST
import Data.AEq( AEq(..) )
import Data.Typeable( Typeable )
import Foreign( peekElemOff )
import Text.Printf( printf )
import Unsafe.Coerce( unsafeCoerce )

import Numeric.LinearAlgebra.Internal( inlinePerformIO )

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Vector.Base( Vector, unVector, unSTVector )
import qualified Numeric.LinearAlgebra.Vector as V

import Numeric.LinearAlgebra.Matrix.STBase


infixl 7 `scale`, `scaleRows`, `scaleCols`
infixl 6 `add`, `shift`, `shiftDiag`, `sub`


-- | Immutable dense matrices. The type arguments are as follows:
--
--     * @e@: the element type of the matrix.
--
newtype Matrix e = Matrix { unMatrix :: STMatrix RealWorld e }
    deriving (RMatrix, Typeable)

unSTMatrix :: STMatrix s e -> Matrix e
unSTMatrix = Matrix . unsafeCoerce
{-# INLINE unSTMatrix #-}

-- | A safe way to create and work with a mutable matrix before returning 
-- an immutable matrix for later perusal. This function avoids copying
-- the matrix before returning it - it uses 'unsafeFreeze' internally,
-- but this wrapper is a safe interface to that function. 
create :: (Storable e) => (forall s . ST s (STMatrix s e)) -> Matrix e
create mx = runST $ mx >>= unsafeFreeze
{-# INLINE create #-}

-- | Converts a mutable matrix to an immutable one by taking a complete
-- copy of it.
freeze :: (Storable e) => STMatrix s e -> ST s (Matrix e)
freeze = fmap unSTMatrix . newCopy
{-# INLINE freeze #-}

-- | Converts a mutable matrix into an immutable matrix. This simply casts
-- the matrix from one type to the other without copying the matrix.
-- Note that because the matrix is possibly not copied, any subsequent
-- modifications made to the mutable version of the matrix may be shared with
-- the immutable version. It is safe to use, therefore, if the mutable
-- version is never modified after the freeze operation.
unsafeFreeze :: (Storable e) => STMatrix s e -> ST s (Matrix e)
unsafeFreeze = return . unSTMatrix
{-# INLINE unsafeFreeze #-}


-- | Create a matrix with the given dimension and elements.  The elements
-- given in the association list must all have unique indices, otherwise
-- the result is undefined.
--
-- Not every index within the bounds of the matrix need appear in the
-- association list, but the values associated with indices that do not
-- appear will be undefined.
fromAssocs :: (Storable e) => (Int,Int) -> [((Int,Int), e)] -> Matrix e
fromAssocs mn ies = create $ do
    a <- new_ mn
    setAssocs a ies
    return a
{-# INLINE fromAssocs #-}

-- | Same as 'fromAssocs', but does not range-check the indices.
unsafeFromAssocs :: (Storable e) => (Int,Int) -> [((Int,Int), e)] -> Matrix e
unsafeFromAssocs mn ies = create $ do
    a <- new_ mn
    unsafeSetAssocs a ies
    return a
{-# INLINE unsafeFromAssocs #-}

-- | Create a matrix of the given dimension with elements initialized
-- to the values from the list, in column major order.
fromList :: (Storable e) => (Int,Int) -> [e] -> Matrix e
fromList mn es = create $ do
    a <- new_ mn
    setElems a es
    return a
{-# INLINE fromList #-}

-- | Create a matrix of the given dimension with the given vectors as
-- columns.
fromCols :: (Storable e) => (Int,Int) -> [Vector e] -> Matrix e
fromCols mn cs = create $ do
    a <- new_ mn
    withSTColViews a $ zipWithM_ V.copyTo cs
    return a

-- | Create a matrix of the given dimension with the given vectors as
-- rows.
fromRows :: (Storable e) => (Int,Int) -> [Vector e] -> Matrix e
fromRows (m,n) rs = create $ do
    a <- new_ (m,n)
    sequence_ [ setRow a i r | (i,r) <- zip [ 0..m-1 ] rs ]
    return a

-- | Create a matrix of the given dimension with all elements initialized
-- to the given value
constant :: (Storable e) => (Int,Int) -> e -> Matrix e
constant mn e = create $ new mn e
{-# INLINE constant #-}

-- | Returns the element of a matrix at the specified index.
at :: (Storable e) => Matrix e -> (Int,Int) -> e
at a ij@(i,j)
    | i < 0 || i >= m || j < 0 || j >= n = error $
        printf ("at <matrix with dim (%d,%d)> (%d,%d):"
                ++ " invalid index") m n i j
    | otherwise =
        unsafeAt a ij
  where
      (m,n) = dim a
{-# INLINE at #-}

unsafeAt :: (Storable e) => Matrix e -> (Int,Int) -> e
unsafeAt (Matrix (STMatrix a _ _ lda)) (i,j) = inlinePerformIO $ 
    V.unsafeWith a $ \p ->
        peekElemOff p (i + j * lda)
{-# INLINE unsafeAt #-}

-- | Returns a list of the elements of a matrix, in the same order as their
-- indices.
elems :: (Storable e) => Matrix e -> [e]
elems a = concatMap V.elems (cols a)
{-# INLINE elems #-}

-- | Returns the contents of a matrix as a list of associations.
assocs :: (Storable e) => Matrix e -> [((Int,Int),e)]
assocs x = zip (indices x) (elems x)
{-# INLINE assocs #-}

unsafeReplace :: (Storable e) => Matrix e -> [((Int,Int),e)] -> Matrix e
unsafeReplace a ies = create $ do
    a' <- newCopy a
    unsafeSetAssocs a' ies
    return a'

-- | Create a new matrix by replacing the values at the specified indices.
replace :: (Storable e) => Matrix e -> [((Int,Int),e)] -> Matrix e
replace a ies = create $ do
    a' <- newCopy a
    setAssocs a' ies
    return a'

-- | @accum f@ takes a matrix and an association list and accumulates
-- pairs from the list into the matrix with the accumulating function @f@.
accum :: (Storable e)
            => (e -> e' -> e) 
            -> Matrix e
            -> [((Int,Int), e')]
            -> Matrix e
accum f a ies = create $ do
    a' <- newCopy a
    forM_ ies $ \(i,e') -> do
        e <- read a' i
        unsafeWrite a' i (f e e') -- index checked on prev. line
    return a'

-- | Same as 'accum' but does not range-check indices.
unsafeAccum :: (Storable e)
                  => (e -> e' -> e)
                  -> Matrix e
                  -> [((Int,Int), e')]
                  -> Matrix e
unsafeAccum f a ies = create $ do
    a' <- newCopy a
    forM_ ies $ \(i,e') -> do
        e <- unsafeRead a' i
        unsafeWrite a' i (f e e')
    return a'

-- | Construct a new matrix by applying a function to every element of
-- a matrix.
map :: (Storable e, Storable e')
    => (e -> e')
    -> Matrix e
    -> Matrix e'
map f a = create $ do
    a' <- new_ (dim a)
    unsafeMapTo f a a'
    return a'
{-# INLINE map #-}

-- | Construct a new matrix by applying a function to every pair of elements
-- of two matrices.  The two matrices must have identical dimensions.
zipWith :: (Storable e, Storable e', Storable f)
        => (e -> e' -> f)
        -> Matrix e
        -> Matrix e'
        -> Matrix f
zipWith f a a'
    | m /= m' || n /= n' = error $
        printf ("zipWith <function> <matrix with dim (%d,%d)>"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
                m n m' n'
    | otherwise =
        unsafeZipWith f a a'
  where
    (m,n) = dim a
    (m',n') = dim a'
{-# INLINE zipWith #-}

-- | Version of 'zipWith' that does not check if the input matrices
-- have the same dimensions.
unsafeZipWith :: (Storable e, Storable e', Storable f)
              => (e -> e' -> f)
              -> Matrix e
              -> Matrix e'
              -> Matrix f
unsafeZipWith f a a' =
    fromList (dim a') $ P.zipWith f (elems a) (elems a')
{-# INLINE unsafeZipWith #-}

-- | Get the given column of the matrix.
col :: (Storable e) => Matrix e -> Int -> Vector e
col a j
    | j < 0 || j >= n = error $
        printf ("col <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n j
    | otherwise =
        unsafeCol a j
  where
    (m,n) = dim a
{-# INLINE col #-}

unsafeCol :: (Storable e) => Matrix e -> Int -> Vector e
unsafeCol a j = unSTVector $ unsafeSTColView (unMatrix a) j
{-# INLINE unsafeCol #-}

-- | Get a list of the columns of the matrix.
cols :: (Storable e) => Matrix e -> [Vector e]
cols a = P.map unSTVector $ stCols (unMatrix a)
{-# INLINE cols #-}

-- | Get the given row of the matrix.
row :: (Storable e) => Matrix e -> Int -> Vector e
row a i = V.create $ do
    x <- V.new_ n
    getRow a i x
    return x
  where
    (_,n) = dim a

unsafeRow :: (Storable e) => Matrix e -> Int -> Vector e
unsafeRow a i = V.create $ do
    x <- V.new_ n
    unsafeGetRow a i x
    return x
  where
    (_,n) = dim a

-- | Get the diagonal of the matrix.
diag :: (Storable e) => Matrix e -> Vector e
diag a = V.create $ do
    x <- V.new_ mn
    getDiag a x
    return x
  where
    (m,n) = dim a
    mn = min m n

-- | Get a list of the rows of the matrix.
rows :: (Storable e) => Matrix e -> [Vector e]
rows a = [ unsafeRow a i | i <- [ 0..m-1 ] ]
  where
    (m,_) = dim a
{-# INLINE rows #-}

-- | Convert a matrix to a vector by stacking its columns.
toVector :: (Storable e)
         => Matrix e
         -> Vector e
toVector a = case maybeSTVectorView (unMatrix a) of
    Just v  -> unSTVector v
    Nothing -> V.concat $ cols a

-- | Cast a vector to a matrix of the given shape.
fromVector :: (Storable e)
           => (Int,Int)
           -> Vector e
           -> Matrix e
fromVector mn v = unSTMatrix $ fromSTVector mn (unVector v)

-- | Cast a vector to a matrix with one column.
fromCol :: (Storable e)
              => Vector e
              -> Matrix e
fromCol v = unSTMatrix $ fromSTCol (unVector v)

-- | Cast a vector to a matrix with one row.
fromRow :: (Storable e)
              => Vector e
              -> Matrix e
fromRow v = unSTMatrix $ fromSTRow (unVector v)

instance (Storable e, Show e) => Show (Matrix e) where
    show x = "fromList " ++ show (dim x) ++ " " ++ show (elems x)
    {-# INLINE show #-}

instance (Storable e, Eq e) => Eq (Matrix e) where
    (==) = compareWith (==)
    {-# INLINE (==) #-}

instance (Storable e, AEq e) => AEq (Matrix e) where
    (===) = compareWith (===)
    {-# INLINE (===) #-}
    (~==) = compareWith (~==)
    {-# INLINE (~==) #-}

-- | @shift k a@ returns @k + a@.
shift :: (VNum e) => e -> Matrix e -> Matrix e
shift k = result $ shiftTo k

-- | @shiftDiag d a@ returns @diag(d) + a@.
shiftDiag :: (BLAS1 e) => Vector e -> Matrix e -> Matrix e
shiftDiag s = result $ shiftDiagTo s

-- | @shiftDiagWithScale alpha d a@ returns @alpha * diag(d) + a@.
shiftDiagWithScale :: (BLAS1 e) => e -> Vector e -> Matrix e -> Matrix e
shiftDiagWithScale e s = result $ shiftDiagToWithScale e s

-- | @add a b@ returns @a + b@.
add :: (VNum e) => Matrix e -> Matrix e -> Matrix e
add = result2 addTo

-- | @addWithScales alpha a beta b@ returns @alpha*a + beta*b@.
addWithScales :: (VNum e) => e -> Matrix e -> e -> Matrix e -> Matrix e
addWithScales alpha a beta b =
    (result2 $ \a' b' -> addToWithScales alpha a' beta b') a b

-- | @sub a b@ returns @a - b@.
sub :: (VNum e) => Matrix e -> Matrix e -> Matrix e
sub = result2 subTo

-- | @scale k a@ returns @k * a@.
scale :: (VNum e) => e -> Matrix e -> Matrix e
scale k = result $ scaleTo k

-- | @scaleRows s a@ returns @diag(s) * a@.
scaleRows :: (VNum e) => Vector e -> Matrix e -> Matrix e
scaleRows s = result $ scaleRowsTo s

-- | @scaleCols s a@ returns @a * diag(s)@.
scaleCols :: (VNum e) => Vector e -> Matrix e -> Matrix e
scaleCols s = result $ scaleColsTo s

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
        transTo a a'
        return a'

-- | @conjTrans a@ retunrs @conj(trans(a))@.
conjTrans :: (BLAS1 e)
          => Matrix e
          -> Matrix e
conjTrans a = let
    (m,n) = dim a
    in create $ do
        a' <- new_ (n,m)
        conjTransTo a a'
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
        rank1UpdateTo alpha x y a'
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
        mulToVector transa a x y
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
        mulToVectorWithScale alpha transa a x y
        return y
                       
-- | @mulAddVectorWithScales alpha transa a x beta y@
-- returns @alpha * op(a) * x + beta * y@, where @op(a)@ is
-- determined by @transa@.
mulAddVectorWithScales :: (BLAS2 e)
                       => e
                       -> Trans -> Matrix e
                       -> Vector e
                       -> e
                       -> Vector e
                       -> Vector e
mulAddVectorWithScales alpha transa a x beta y =
    V.create $ do
        y' <- V.newCopy y
        mulAddToVectorWithScales alpha transa a x beta y'
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
        mulToMatrix transa a transb b c
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
        mulToMatrixWithScale alpha transa a transb b c
        return c

-- | @mulAddMatrixWithScales alpha transa a transb b beta c@
-- returns @alpha * op(a) * op(b) + beta * c@, where @op(a)@ and
-- @op(b)@ are determined by @transa@ and @transb@.
mulAddMatrixWithScales :: (BLAS3 e)
                       => e
                       -> Trans -> Matrix e
                       -> Trans -> Matrix e
                       -> e
                       -> Matrix e
                       -> Matrix e
mulAddMatrixWithScales alpha transa a transb b beta c = 
    create $ do
        c' <- newCopy c
        mulAddToMatrixWithScales alpha transa a transb b beta c'
        return c'

compareWith :: (Storable e, Storable e')
            => (e -> e' -> Bool)
            -> Matrix e
            -> Matrix e'
            -> Bool
compareWith cmp a a' =
    dim a == dim a'
    && and (P.zipWith cmp (elems a) (elems a'))
{-# INLINE compareWith #-}

result :: (Storable e, Storable f)
       => (forall s . Matrix e -> STMatrix s f -> ST s a)
       -> Matrix e
       -> Matrix f
result f a = create $ newResult f a
{-# INLINE result #-}

result2 :: (Storable e, Storable f, Storable g)
        => (forall s . Matrix e -> Matrix f -> STMatrix s g -> ST s a)
        -> Matrix e
        -> Matrix f
        -> Matrix g
result2 f a1 a2 = create $ newResult2 f a1 a2
{-# INLINE result2 #-}
