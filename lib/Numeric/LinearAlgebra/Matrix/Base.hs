{-# LANGUAGE DeriveDataTypeable, GeneralizedNewtypeDeriving, Rank2Types,
        TypeFamilies #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix.Base
-- Copyright  : Copyright (c) , Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Numeric.LinearAlgebra.Matrix.Base
    where

import Control.Monad( forM_, zipWithM_ )
import Control.Monad.ST
import Data.AEq( AEq(..) )
import Data.Typeable( Typeable )
import Foreign( peekElemOff )
import Text.Printf( printf )
import Unsafe.Coerce( unsafeCoerce )

import Numeric.LinearAlgebra.Internal( inlinePerformIO )

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Vector.Base
import Numeric.LinearAlgebra.Vector.STBase

import Numeric.LinearAlgebra.Matrix.STBase


infixl 7 `scaleMatrix`, `scaleRowsMatrix`, `scaleColsMatrix`
infixl 6 `addMatrix`, `shiftMatrix`, `shiftDiagMatrix`, `subMatrix`


-- | Immutable dense matrices. The type arguments are as follows:
--
--     * @e@: the element type of the matrix.
--
newtype Matrix e = Matrix { unMatrix :: STMatrix RealWorld e }
    deriving (RMatrix, Typeable)

instance HasVectorView Matrix where
    type VectorView Matrix = Vector

-- | A safe way to create and work with a mutable matrix before returning 
-- an immutable matrix for later perusal. This function avoids copying
-- the matrix before returning it - it uses 'unsafeFreezeMatrix' internally,
-- but this wrapper is a safe interface to that function. 
runMatrix :: (Storable e) => (forall s . ST s (STMatrix s e)) -> Matrix e
runMatrix mx = runST $ mx >>= unsafeFreezeMatrix
{-# INLINE runMatrix #-}

-- | Converts a mutable matrix to an immutable one by taking a complete
-- copy of it.
freezeMatrix :: (Storable e) => STMatrix s e -> ST s (Matrix e)
freezeMatrix = fmap Matrix . unsafeCoerce . newCopyMatrix
{-# INLINE freezeMatrix #-}

-- | Converts a mutable matrix into an immutable matrix. This simply casts
-- the matrix from one type to the other without copying the matrix.
-- Note that because the matrix is possibly not copied, any subsequent
-- modifications made to the mutable version of the matrix may be shared with
-- the immutable version. It is safe to use, therefore, if the mutable
-- version is never modified after the freeze operation.
unsafeFreezeMatrix :: (Storable e) => STMatrix s e -> ST s (Matrix e)
unsafeFreezeMatrix = return . Matrix . unsafeCoerce
{-# INLINE unsafeFreezeMatrix #-}


-- | Create a matrix with the given dimension and elements.  The elements
-- given in the association list must all have unique indices, otherwise
-- the result is undefined.
--
-- Not every index within the bounds of the matrix need appear in the
-- association list, but the values associated with indices that do not
-- appear will be undefined.
matrix :: (Storable e) => (Int,Int) -> [((Int,Int), e)] -> Matrix e
matrix mn ies = runMatrix $ do
    a <- newMatrix_ mn
    setAssocsMatrix a ies
    return a
{-# INLINE matrix #-}

-- | Same as 'matrix', but does not range-check the indices.
unsafeMatrix :: (Storable e) => (Int,Int) -> [((Int,Int), e)] -> Matrix e
unsafeMatrix mn ies = runMatrix $ do
    a <- newMatrix_ mn
    unsafeSetAssocsMatrix a ies
    return a
{-# INLINE unsafeMatrix #-}

-- | Create a matrix of the given dimension with elements initialized
-- to the values from the list, in column major order.
listMatrix :: (Storable e) => (Int,Int) -> [e] -> Matrix e
listMatrix mn es = runMatrix $ do
    a <- newMatrix_ mn
    setElemsMatrix a es
    return a
{-# INLINE listMatrix #-}

-- | Create a matrix of the given dimension with the given vectors as
-- columns.
colListMatrix :: (Storable e) => (Int,Int) -> [Vector e] -> Matrix e
colListMatrix mn cs = runMatrix $ do
    a <- newMatrix_ mn
    zipWithM_ copyToVector cs (colsMatrix a)
    return a

-- | Create a matrix of the given dimension with the given vectors as
-- rows.
rowListMatrix :: (Storable e) => (Int,Int) -> [Vector e] -> Matrix e
rowListMatrix (m,n) rs = runMatrix $ do
    a <- newMatrix_ (m,n)
    sequence_ [ setRowMatrix a i r | (i,r) <- zip [ 0..m-1 ] rs ]
    return a

-- | Create a matrix of the given dimension with all elements initialized
-- to the given value
constantMatrix :: (Storable e) => (Int,Int) -> e -> Matrix e
constantMatrix mn e = runMatrix $ newMatrix mn e
{-# INLINE constantMatrix #-}

-- | Returns the element of a matrix at the specified index.
atMatrix :: (Storable e) => Matrix e -> (Int,Int) -> e
atMatrix a ij@(i,j)
    | i < 0 || i >= m || j < 0 || j >= n = error $
        printf ("atMatrix <matrix with dim (%d,%d)> (%d,%d):"
                ++ " invalid index") m n i j
    | otherwise =
        unsafeAtMatrix a ij
  where
      (m,n) = dimMatrix a
{-# INLINE atMatrix #-}

unsafeAtMatrix :: (Storable e) => Matrix e -> (Int,Int) -> e
unsafeAtMatrix (Matrix (STMatrix a _ _ lda)) (i,j) = inlinePerformIO $ 
    unsafeWithVector a $ \p ->
        peekElemOff p (i + j * lda)
{-# INLINE unsafeAtMatrix #-}

-- | Returns a list of the elements of a matrix, in the same order as their
-- indices.
elemsMatrix :: (Storable e) => Matrix e -> [e]
elemsMatrix a = concatMap elemsVector (colsMatrix a)
{-# INLINE elemsMatrix #-}

-- | Returns the contents of a matrix as a list of associations.
assocsMatrix :: (Storable e) => Matrix e -> [((Int,Int),e)]
assocsMatrix x = zip (indicesMatrix x) (elemsMatrix x)
{-# INLINE assocsMatrix #-}

unsafeReplaceMatrix :: (Storable e) => Matrix e -> [((Int,Int),e)] -> Matrix e
unsafeReplaceMatrix a ies = runMatrix $ do
    a' <- newCopyMatrix a
    unsafeSetAssocsMatrix a' ies
    return a'

-- | Create a new matrix by replacing the values at the specified indices.
replaceMatrix :: (Storable e) => Matrix e -> [((Int,Int),e)] -> Matrix e
replaceMatrix a ies = runMatrix $ do
    a' <- newCopyMatrix a
    setAssocsMatrix a' ies
    return a'

-- | @accumMatrix f@ takes a matrix and an association list and accumulates
-- pairs from the list into the matrix with the accumulating function @f@.
accumMatrix :: (Storable e)
            => (e -> e' -> e) 
            -> Matrix e
            -> [((Int,Int), e')]
            -> Matrix e
accumMatrix f a ies = runMatrix $ do
    a' <- newCopyMatrix a
    forM_ ies $ \(i,new) -> do
        old <- readMatrix a' i
        unsafeWriteMatrix a' i (f old new) -- index checked on prev. line
    return a'

-- | Same as 'accumMatrix' but does not range-check indices.
unsafeAccumMatrix :: (Storable e)
                  => (e -> e' -> e)
                  -> Matrix e
                  -> [((Int,Int), e')]
                  -> Matrix e
unsafeAccumMatrix f a ies = runMatrix $ do
    a' <- newCopyMatrix a
    forM_ ies $ \(i,new) -> do
        old <- unsafeReadMatrix a' i
        unsafeWriteMatrix a' i (f old new)
    return a'

-- | Construct a new matrix by applying a function to every element of
-- a matrix.
mapMatrix :: (Storable e, Storable e')
          => (e -> e')
          -> Matrix e
          -> Matrix e'
mapMatrix f a = runMatrix $ do
    a' <- newMatrix_ (dimMatrix a)
    unsafeMapToMatrix f a a'
    return a'
{-# INLINE mapMatrix #-}

-- | Construct a new matrix by applying a function to every pair of elements
-- of two matrices.  The two matrices must have identical dimensions.
zipWithMatrix :: (Storable e, Storable e', Storable f)
              => (e -> e' -> f)
              -> Matrix e
              -> Matrix e'
              -> Matrix f
zipWithMatrix f a a'
    | m /= m' || n /= n' = error $
        printf ("zipWithMatrix <function> <matrix with dim (%d,%d)>"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
                m n m' n'
    | otherwise =
        unsafeZipWithMatrix f a a'
  where
    (m,n) = dimMatrix a
    (m',n') = dimMatrix a'
{-# INLINE zipWithMatrix #-}

-- | Version of 'zipWithMatrix' that does not check if the input matrices
-- have the same dimensions.
unsafeZipWithMatrix :: (Storable e, Storable e', Storable f)
                    => (e -> e' -> f)
                    -> Matrix e
                    -> Matrix e'
                    -> Matrix f
unsafeZipWithMatrix f a a' =
    listMatrix (dimMatrix a') $ zipWith f (elemsMatrix a) (elemsMatrix a')
{-# INLINE unsafeZipWithMatrix #-}

-- | Get the given row of the matrix.
rowMatrix :: (Storable e) => Matrix e -> Int -> Vector e
rowMatrix a i = runVector $ do
    x <- newVector_ n
    getRowMatrix a i x
    return x
  where
    (_,n) = dimMatrix a

unsafeRowMatrix :: (Storable e) => Matrix e -> Int -> Vector e
unsafeRowMatrix a i = runVector $ do
    x <- newVector_ n
    unsafeGetRowMatrix a i x
    return x
  where
    (_,n) = dimMatrix a

-- | Get a list of the rows of the matrix.
rowsMatrix :: (Storable e) => Matrix e -> [Vector e]
rowsMatrix a = [ unsafeRowMatrix a i | i <- [ 0..m-1 ] ]
  where
    (m,_) = dimMatrix a

instance (Storable e, Show e) => Show (Matrix e) where
    show x = "listMatrix " ++ show (dimMatrix x) ++ " " ++ show (elemsMatrix x)
    {-# INLINE show #-}

instance (Storable e, Eq e) => Eq (Matrix e) where
    (==) = compareMatrixWith (==)
    {-# INLINE (==) #-}

instance (Storable e, AEq e) => AEq (Matrix e) where
    (===) = compareMatrixWith (===)
    {-# INLINE (===) #-}
    (~==) = compareMatrixWith (~==)
    {-# INLINE (~==) #-}

-- | @shiftMatrix k a@ returns @k + a@.
shiftMatrix :: (VNum e) => e -> Matrix e -> Matrix e
shiftMatrix k = resultMatrix $ shiftToMatrix k

-- | @shiftDiagMatrix d a@ returns @diag(d) + a@.
shiftDiagMatrix :: (BLAS1 e) => Vector e -> Matrix e -> Matrix e
shiftDiagMatrix s = resultMatrix $ shiftDiagToMatrix s

-- | @shiftDiagMatrixWithScale alpha d a@ returns @alpha * diag(d) + a@.
shiftDiagMatrixWithScale :: (BLAS1 e) => e -> Vector e -> Matrix e -> Matrix e
shiftDiagMatrixWithScale e s = resultMatrix $ shiftDiagToMatrixWithScale e s

-- | @addMatrix a b@ returns @a + b@.
addMatrix :: (VNum e) => Matrix e -> Matrix e -> Matrix e
addMatrix = resultMatrix2 addToMatrix

-- | @addMatrixWithScales alpha a beta b@ returns @alpha*a + beta*b@.
addMatrixWithScales :: (VNum e) => e -> Matrix e -> e -> Matrix e -> Matrix e
addMatrixWithScales alpha a beta b =
    (resultMatrix2 $ \a' b' -> addToMatrixWithScales alpha a' beta b') a b

-- | @subMatrix a b@ returns @a - b@.
subMatrix :: (VNum e) => Matrix e -> Matrix e -> Matrix e
subMatrix = resultMatrix2 subToMatrix

-- | @scaleMatrix k a@ returns @k * a@.
scaleMatrix :: (VNum e) => e -> Matrix e -> Matrix e
scaleMatrix k = resultMatrix $ scaleToMatrix k

-- | @scaleRowsMatrix s a@ returns @diag(s) * a@.
scaleRowsMatrix :: (VNum e) => Vector e -> Matrix e -> Matrix e
scaleRowsMatrix s = resultMatrix $ scaleRowsToMatrix s

-- | @scaleColsMatrix s a@ returns @a * diag(s)@.
scaleColsMatrix :: (VNum e) => Vector e -> Matrix e -> Matrix e
scaleColsMatrix s = resultMatrix $ scaleColsToMatrix s

-- | @negateMatrix a@ returns @-a@.
negateMatrix :: (VNum e) => Matrix e -> Matrix e
negateMatrix = resultMatrix negateToMatrix

-- | @conjMatrix a@ returns @conj(a)@.
conjMatrix :: (VNum e) => Matrix e -> Matrix e
conjMatrix = resultMatrix conjToMatrix

-- | @transMatrix a@ retunrs @trans(a)@.
transMatrix :: (BLAS1 e)
            => Matrix e
            -> Matrix e
transMatrix a = let
    (m,n) = dimMatrix a
    in runMatrix $ do
        a' <- newMatrix_ (n,m)
        transToMatrix a a'
        return a'

-- | @conjTransMatrix a@ retunrs @conj(trans(a))@.
conjTransMatrix :: (BLAS1 e, VNum e)
                => Matrix e
                -> Matrix e
conjTransMatrix a = let
    (m,n) = dimMatrix a
    in runMatrix $ do
        a' <- newMatrix_ (n,m)
        conjTransToMatrix a a'
        return a'

-- | @rank1UpdateMatrix alpha x y a@ returns @alpha * x * y^H + a@.
rank1UpdateMatrix :: (BLAS2 e)
                  => e
                  -> Vector e
                  -> Vector e
                  -> Matrix e
                  -> Matrix e
rank1UpdateMatrix alpha x y a =
    runMatrix $ do
        a' <- newMatrix_ (dimMatrix a)
        rank1UpdateToMatrix alpha x y a a'
        return a'

-- | @mulMatrixVector transa a x@
-- returns @op(a) * x@, where @op(a)@ is determined by @transa@.                   
mulMatrixVector :: (BLAS2 e)
                => Trans -> Matrix e
                -> Vector e
                -> Vector e
mulMatrixVector transa a x = let
    m = case transa of NoTrans -> (fst . dimMatrix) a
                       _       -> (snd . dimMatrix) a
    in runVector $ do
        y <- newVector_ m
        mulMatrixToVector transa a x y
        return y

-- | @mulMatrixVectorWithScale alpha transa a x@
-- retunrs @alpha * op(a) * x@, where @op(a)@ is determined by @transa@.                   
mulMatrixVectorWithScale :: (BLAS2 e)
                         => e
                         -> Trans -> Matrix e
                         -> Vector e
                         -> Vector e
mulMatrixVectorWithScale alpha transa a x = let
    m = case transa of NoTrans -> (fst . dimMatrix) a
                       _       -> (snd . dimMatrix) a
    in runVector $ do
        y <- newVector_ m
        mulMatrixToVectorWithScale alpha transa a x y
        return y
                       
-- | @mulMatrixAddVector transa a x y@
-- returns @op(a) * x + y@, where @op(a)@ is determined by @transa@.                   
mulMatrixAddVector :: (BLAS2 e)
                   => Trans -> Matrix e
                   -> Vector e
                   -> Vector e
                   -> Vector e
mulMatrixAddVector transa a x y =
    runVector $ do
        y' <- newVector_ (dimVector y)
        mulMatrixAddToVector transa a x y y'
        return y'

-- | @mulMatrixAddVectorWithScales alpha transa a x beta y@
-- returns @alpha * op(a) * x + beta * y@, where @op(a)@ is
-- determined by @transa@.
mulMatrixAddVectorWithScales :: (BLAS2 e)
                             => e
                             -> Trans -> Matrix e
                             -> Vector e
                             -> e
                             -> Vector e
                             -> Vector e
mulMatrixAddVectorWithScales alpha transa a x beta y =
    runVector $ do
        y' <- newVector_ (dimVector y)
        mulMatrixAddToVectorWithScales alpha transa a x beta y y'
        return y'

-- | @mulMatrixMatrix transa a transb b@
-- returns @op(a) * op(b)@, where @op(a)@ and @op(b)@ are determined
-- by @transa@ and @transb@.                   
mulMatrixMatrix :: (BLAS3 e)
                => Trans -> Matrix e
                -> Trans -> Matrix e
                -> Matrix e
mulMatrixMatrix transa a transb b = let
    m = case transa of NoTrans -> (fst . dimMatrix) a
                       _       -> (snd . dimMatrix) a
    n = case transb of NoTrans -> (snd . dimMatrix) b
                       _       -> (fst . dimMatrix) b
    in runMatrix $ do
        c <- newMatrix_ (m,n)
        mulMatrixToMatrix transa a transb b c
        return c

-- | @mulMatrixMatrixWithScale alpha transa a transb b@
-- returns @alpha * op(a) * op(b)@, where @op(a)@ and @op(b)@ are determined
-- by @transa@ and @transb@.                   
mulMatrixMatrixWithScale :: (BLAS3 e)
                         => e
                         -> Trans -> Matrix e
                         -> Trans -> Matrix e
                         -> Matrix e
mulMatrixMatrixWithScale alpha transa a transb b = let
    m = case transa of NoTrans -> (fst . dimMatrix) a
                       _       -> (snd . dimMatrix) a
    n = case transb of NoTrans -> (snd . dimMatrix) b
                       _       -> (fst . dimMatrix) b
    in runMatrix $ do
        c <- newMatrix_ (m,n)
        mulMatrixToMatrixWithScale alpha transa a transb b c
        return c

-- | @mulMatrixAddMatrix transa a transb b c@
-- returns @op(a) * op(b) + c@, where @op(a)@ and @op(b)@ are determined
-- by @transa@ and @transb@.                   
mulMatrixAddMatrix :: (BLAS3 e)
                   => Trans -> Matrix e
                   -> Trans -> Matrix e
                   -> Matrix e
                   -> Matrix e
mulMatrixAddMatrix transa a transb b c =
    runMatrix $ do
        c' <- newMatrix_ (dimMatrix c)
        mulMatrixAddToMatrix transa a transb b c c'
        return c'

-- | @mulMatrixAddMatrixWithScales alpha transa a transb b beta c@
-- returns @alpha * op(a) * op(b) + beta * c@, where @op(a)@ and
-- @op(b)@ are determined by @transa@ and @transb@.
mulMatrixAddMatrixWithScales :: (BLAS3 e)
                             => e
                             -> Trans -> Matrix e
                             -> Trans -> Matrix e
                             -> e
                             -> Matrix e
                             -> Matrix e
mulMatrixAddMatrixWithScales alpha transa a transb b beta c = 
    runMatrix $ do
        c' <- newMatrix_ (dimMatrix c)
        mulMatrixAddToMatrixWithScales alpha transa a transb b beta c c'
        return c'

compareMatrixWith :: (Storable e, Storable e')
                  => (e -> e' -> Bool)
                  -> Matrix e
                  -> Matrix e'
                  -> Bool
compareMatrixWith cmp a a' =
    dimMatrix a == dimMatrix a'
    && and (zipWith cmp (elemsMatrix a) (elemsMatrix a'))
{-# INLINE compareMatrixWith #-}

resultMatrix :: (Storable e, Storable f)
             => (forall s . Matrix e -> STMatrix s f -> ST s a)
             -> Matrix e
             -> Matrix f
resultMatrix f a = runMatrix $ newResultMatrix f a
{-# INLINE resultMatrix #-}

resultMatrix2 :: (Storable e, Storable f, Storable g)
              => (forall s . Matrix e -> Matrix f -> STMatrix s g -> ST s a)
              -> Matrix e
              -> Matrix f
              -> Matrix g
resultMatrix2 f a1 a2 = runMatrix $ newResultMatrix2 f a1 a2
{-# INLINE resultMatrix2 #-}
