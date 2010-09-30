{-# LANGUAGE DeriveDataTypeable, FlexibleContexts, Rank2Types #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix.STBase
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : asinerimental
--

module Numeric.LinearAlgebra.Matrix.STBase
    where
      
import Control.Monad( forM_, zipWithM_ )
import Control.Monad.ST( ST, RealWorld, unsafeIOToST )
import Data.Maybe( fromMaybe )
import Data.Typeable( Typeable )
import Foreign( ForeignPtr, Ptr, advancePtr, peekElemOff, pokeElemOff,
    mallocForeignPtrArray )
import Text.Printf( printf )
import Unsafe.Coerce( unsafeCoerce )

import Numeric.LinearAlgebra.Types
import qualified Foreign.BLAS as BLAS
import Numeric.LinearAlgebra.Vector.STBase( STVector, RVector )
import qualified Numeric.LinearAlgebra.Vector.STBase as V

-- | Dense matrices in the 'ST' monad.  The type arguments are as follows:
--
--     * @s@: the state variable argument for the 'ST' type
--
--     * @e@: the element type of the matrix.
--
data STMatrix s e =
      STMatrix {-# UNPACK #-} !(STVector s e) -- matrix data
               {-# UNPACK #-} !Int            -- row dimension
               {-# UNPACK #-} !Int            -- column dimension
               {-# UNPACK #-} !Int            -- leading dimension
  deriving (Typeable)

-- | Dense matrices in the 'IO' monad.
type IOMatrix = STMatrix RealWorld

-- | Read-only matrices
class RMatrix m where
    -- | The dimensions of the matrix (number of rows and columns).
    dimMatrix :: (Storable e) => m e -> (Int,Int)
    
    unsafeSliceMatrix :: (Storable e)
                      => (Int,Int) -> (Int,Int) -> m e -> m e

    unsafeWithColMatrix :: (Storable e)
                        => m e
                        -> Int 
                        -> (forall v . RVector v => v e -> ST s a)
                        -> ST s a

    -- | Perform an action with a list of views of the matrix columns.
    withColsMatrix :: (Storable e)
                   => m e
                   -> (forall v . RVector v => [v e] -> ST s a)
                   -> ST s a

    -- | Possibly view a matrix as a vector and perform an action on the
    -- view.  This only succeeds if the matrix is stored contiguously in
    -- memory, i.e. if the matrix contains a single column or the \"lda\"
    -- of the matrix is equal to the number of rows.
    maybeWithVectorViewMatrix :: (Storable e)
                              => m e
                              -> (forall v . RVector v => v e -> ST s a)
                              -> Maybe (ST s a)

    -- | Execute an 'IO' action with a pointer to the first element in the
    -- matrix and the leading dimension (lda).
    unsafeWithMatrix :: (Storable e) => m e -> (Ptr e -> Int -> IO a) -> IO a

    -- | Convert a vector to a @ForeignPtr@, and offset, dimensions,
    -- and leading dimension (lda).
    unsafeMatrixToForeignPtr :: (Storable e)
                             => m e -> (ForeignPtr e, Int, (Int,Int), Int)
                             
    -- | Cast a @ForeignPtr@ to a matrix.
    unsafeMatrixFromForeignPtr :: (Storable e)
                               => ForeignPtr e -- ^ the pointer
                               -> Int          -- ^ the offset
                               -> (Int,Int)    -- ^ the dimensions
                               -> Int          -- ^ the leading dimension
                               -> m e


instance RMatrix (STMatrix s) where
    dimMatrix (STMatrix _ m n _) = (m,n)
    {-# INLINE dimMatrix #-}

    unsafeWithColMatrix a j f = cast $ unsafeWithColMatrixST a j (cast . f)
      where
        cast :: ST s' a -> ST s a
        cast = unsafeCoerce
    {-# INLINE unsafeWithColMatrix #-}

    withColsMatrix a f = cast $ withColsMatrixST a (cast . f)
      where
        cast :: ST s' a -> ST s a
        cast = unsafeCoerce
    {-# INLINE withColsMatrix #-}

    unsafeSliceMatrix (i,j) (m',n') (STMatrix a _ _ lda) = let
        o = i + j*lda
        l = if m' == 0 || n' == 0
                then 0
                else lda * (n' - 1) + m'
        a' = V.unsafeSlice o l a
        in STMatrix a' m' n' lda
    {-# INLINE unsafeSliceMatrix #-}
    
    maybeWithVectorViewMatrix a f =
        cast `fmap` maybeWithVectorViewMatrixST a (cast . f)
      where
        cast :: ST s' a -> ST s a
        cast = unsafeCoerce
    {-# INLINE maybeWithVectorViewMatrix #-}

    unsafeWithMatrix (STMatrix a _ _ lda) f =
        V.unsafeWith a $ \p -> f p lda
    {-# INLINE unsafeWithMatrix #-}

    unsafeMatrixToForeignPtr (STMatrix a m n lda) = let
        (f,o,_) = V.unsafeToForeignPtr a
        in (f, o, (m,n), lda)
    {-# INLINE unsafeMatrixToForeignPtr #-}

    unsafeMatrixFromForeignPtr f o (m,n) lda = let
        d = if m == 0 then 0 else lda * n
        a = V.unsafeFromForeignPtr f o d
        in (STMatrix a m n lda)
    {-# INLINE unsafeMatrixFromForeignPtr #-}


unsafeColMatrixST :: (Storable e)
                  => STMatrix s e
                  -> Int
                  -> STVector s e
unsafeColMatrixST (STMatrix a m _ lda) j = let
    o = j * lda
    in V.unsafeSlice o m a
{-# INLINE unsafeColMatrixST #-}

-- | Perform an action with a view of a matrix column (no index checking).
unsafeWithColMatrixST :: (Storable e)
                      => STMatrix s e
                      -> Int
                      -> (STVector s e -> ST s a)
                      -> ST s a
unsafeWithColMatrixST a j f = f $ unsafeColMatrixST a j
{-# INLINE unsafeWithColMatrixST #-}

-- | Get a list of views of the matrix columns.
colsMatrixST :: (Storable e)
             => STMatrix s e
             -> [STVector s e]
colsMatrixST (STMatrix a m n lda) = let
    os = take n $ [ 0, lda .. ]
    xs = [ V.unsafeSlice o m a | o <- os ]
    in xs
{-# INLINE colsMatrixST #-}

-- | Perform an action with a list of views of the matrix columns.
withColsMatrixST :: (Storable e)
                 => STMatrix s e
                 -> ([STVector s e] -> ST s a)
                 -> ST s a
withColsMatrixST a f = f $ colsMatrixST a
{-# INLINE withColsMatrixST #-}

-- | Possibly view a matrix as a vector.  This only succeeds if the
-- matrix is stored contiguously in memory, i.e. if the matrix contains
-- a single column or the \"lda\" of the matrix is equal to the number
-- of rows.
maybeVectorViewMatrixST :: (Storable e)
                        => STMatrix s e
                        -> Maybe (STVector s e)
maybeVectorViewMatrixST (STMatrix a m n lda) | lda == max 1 m = Just a
                                             | n   == 1       = Just a
                                             | otherwise      = Nothing
{-# INLINE maybeVectorViewMatrixST #-}

-- | Possibly view a matrix as a vector and perform an action on the
-- view.  This only succeeds when 'maybeVectorViewMatrixST' does.
maybeWithVectorViewMatrixST :: (Storable e)
                            => STMatrix s e
                            -> (STVector s e -> ST s a)
                            -> Maybe (ST s a)
maybeWithVectorViewMatrixST a f = f `fmap` maybeVectorViewMatrixST a
{-# INLINE maybeWithVectorViewMatrixST #-}

-- | Cast a vector to a matrix of the given shape.  See also
-- 'withMatrixViewVector'.
matrixViewVectorST :: (Storable e)
                   => (Int,Int)
                   -> STVector s e
                   -> STMatrix s e
matrixViewVectorST mn@(m,n) v
    | V.dim v /= m*n = error $
        printf ("matrixViewVectorST (%d,%d) <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim v)
    | otherwise =
        cast v
  where
    cast x = let
        (fptr,o,_) = V.unsafeToForeignPtr x
        lda = max 1 m
        in unsafeMatrixFromForeignPtr fptr o mn lda
{-# INLINE matrixViewVectorST #-}

-- | Cast a vector to a matrix with one column.  See also
-- 'withMatrixViewColVector'.
matrixViewColVectorST :: (Storable e)
                      => STVector s e
                      -> STMatrix s e
matrixViewColVectorST v = matrixViewVectorST (V.dim v, 1) v
{-# INLINE matrixViewColVectorST #-}

-- | Cast a vector to a matrix with one row.  See also
-- 'withMatrixViewRowVector'.
matrixViewRowVectorST :: (Storable e)
                      => STVector s e
                      -> STMatrix s e
matrixViewRowVectorST v = matrixViewVectorST (V.dim v, 1) v
{-# INLINE matrixViewRowVectorST #-}

-- | Cast a vector to a matrix of the given shape and pass it to
-- the specified function.
withMatrixViewVector :: (RVector v, Storable e)
                     => (Int,Int)
                     -> v e
                     -> (forall m . RMatrix m => m e -> a)
                     -> a
withMatrixViewVector mn@(m,n) v f
    | V.dim v /= m*n = error $
        printf ("withMatrixViewVector (%d,%d) <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim v)
    | otherwise =
        f (cast v)
  where
    cast :: (RVector v, Storable e) => v e -> STMatrix s e
    cast x = let
        (fptr,o,_) = V.unsafeToForeignPtr x
        lda = max 1 m
        in unsafeMatrixFromForeignPtr fptr o mn lda
{-# INLINE withMatrixViewVector #-}

-- | Cast a vector to a matrix with one column and pass it to
-- the specified function.
withMatrixViewColVector :: (RVector v, Storable e)
                        => v e
                        -> (forall m . RMatrix m => m e -> a)
                        -> a
withMatrixViewColVector v = withMatrixViewVector (V.dim v, 1) v
{-# INLINE withMatrixViewColVector #-}

-- | Cast a vector to a matrix with one row and pass it to
-- the specified function.
withMatrixViewRowVector :: (RVector v, Storable e)
                        => v e
                        -> (forall m . RMatrix m => m e -> a)
                        -> a
withMatrixViewRowVector v = withMatrixViewVector (1, V.dim v) v
{-# INLINE withMatrixViewRowVector #-}

-- | Perform an action with a view of a matrix column.
withColMatrix :: (RMatrix m, Storable e)
              => m e
              -> Int
              -> (forall v . RVector v => v e -> ST s a)
              -> ST s a
withColMatrix a j f
    | j < 0 || j >= n = error $
        printf ("withColMatrix <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n j
    | otherwise =
        unsafeWithColMatrix a j f
  where
    (m,n) = dimMatrix a
{-# INLINE withColMatrix #-}

-- | Perform an action with a view of a mutable matrix column.
withColMatrixST :: (Storable e)
                => STMatrix s e
                -> Int
                -> (STVector s e -> ST s a)
                -> ST s a
withColMatrixST a j f
    | j < 0 || j >= n = error $
        printf ("withColMatrixST <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n j
    | otherwise =
        unsafeWithColMatrixST a j f
  where
    (m,n) = dimMatrix a
{-# INLINE withColMatrixST #-}

-- | @sliceMatrix (i,j) (m,n) a@ creates a submatrix view of @a@ starting at
-- element @(i,j)@ and having dimensions @(m,n)@.
sliceMatrix :: (RMatrix m, Storable e)
            => (Int,Int)
            -> (Int,Int)
            -> m e
            -> m e
sliceMatrix (i,j) (m',n') a
    | (i < 0 || m' < 0 || i + m' > m 
       || j < 0 || n' < 0 || j + n' > n) = error $
        printf ("sliceMatrix (%d,%d) (%d,%d) <matrix with dim (%d,%d)>:"
                ++ " index out of range") i j m' n' m n
    | otherwise =
        unsafeSliceMatrix (i,j) (m',n') a
  where
    (m,n) = dimMatrix a
{-# INLINE sliceMatrix #-}

-- | Create a new matrix of given shape, but do not initialize the elements.
newMatrix_ :: (Storable e) => (Int,Int) -> ST s (STMatrix s e)
newMatrix_ (m,n) 
    | m < 0 || n < 0 = error $
        printf "newMatrix_ (%d,%d): invalid dimensions" m n
    | otherwise = unsafeIOToST $ do
        f <- mallocForeignPtrArray (m*n)
        return $ unsafeMatrixFromForeignPtr f 0 (m,n) (max 1 m)

-- | Create a matrix with every element initialized to the same value.
newMatrix :: (Storable e) => (Int,Int) -> e -> ST s (STMatrix s e)
newMatrix (m,n) e = do
    a <- newMatrix_ (m,n)
    setElemsMatrix a $ replicate (m*n) e
    return a

-- | Creates a new matrix by copying another one.    
newCopyMatrix :: (RMatrix m, Storable e) => m e -> ST s (STMatrix s e)
newCopyMatrix a = do
    b <- newMatrix_ (dimMatrix a)
    unsafeCopyToMatrix a b
    return b

-- | @copyToMatrix src dst@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyToMatrix :: (RMatrix m, Storable e) => m e -> STMatrix s e -> ST s ()
copyToMatrix = checkMatrixOp2 "copyToMatrix" unsafeCopyToMatrix
{-# INLINE copyToMatrix #-}

unsafeCopyToMatrix :: (RMatrix m, Storable e) => m e -> STMatrix s e -> ST s ()
unsafeCopyToMatrix = matrixVectorOp2 V.unsafeCopyTo
{-# INLINE unsafeCopyToMatrix #-}

-- | Create a view of a matrix by taking the initial rows.
takeRowsMatrix :: (RMatrix m, Storable e) => Int -> m e -> m e
takeRowsMatrix i a = sliceMatrix (0,0) (i,n) a
  where
    (_,n) = dimMatrix a

-- | Create a view of a matrix by droping the initial rows.
dropRowsMatrix :: (RMatrix m, Storable e) => Int -> m e -> m e
dropRowsMatrix i a = sliceMatrix (i,0) (m-i,n) a
  where
    (m,n) = dimMatrix a

-- | Split a matrix into two blocks and returns views into the blocks.  If
-- @(a1, a2) = splitRowsMatrixAt i a@, then
-- @a1 = sliceMatrix (0,0) (i,n) a@ and
-- @a2 = sliceMatrix (i,0) (m-i,n) a@, where @(m,n)@ is the dimension of @a@.
splitRowsMatrixAt :: (RMatrix m, Storable e) => Int -> m e -> (m e, m e)
splitRowsMatrixAt i a
    | i < 0 || i > m = error $
        printf ("splitRowsMatrixAt %d <matrix with dim (%d,%d)>:"
                ++ " invalid index") i m n
    | otherwise = let
        a1 = unsafeSliceMatrix (0,0) (i,n)   a
        a2 = unsafeSliceMatrix (i,0) (m-i,n) a
    in (a1,a2)
  where
    (m,n) = dimMatrix a
{-# INLINE splitRowsMatrixAt #-}

-- | Create a view of a matrix by taking the initial columns.
takeColsMatrix :: (RMatrix m, Storable e) => Int -> m e -> m e
takeColsMatrix j a = sliceMatrix (0,0) (m,j) a
  where
    (m,_) = dimMatrix a

-- | Create a view of a matrix by droping the initial columns.
dropColsMatrix :: (RMatrix m, Storable e) => Int -> m e -> m e
dropColsMatrix j a = sliceMatrix (0,j) (m,n-j) a
  where
    (m,n) = dimMatrix a

-- | Split a matrix into two blocks and returns views into the blocks.  If
-- @(a1, a2) = splitColsMatrixAt j a@, then
-- @a1 = sliceMatrix (0,0) (m,j) a@ and
-- @a2 = sliceMatrix (0,j) (m,n-j) a@, where @(m,n)@ is the dimension of @a@.
splitColsMatrixAt :: (RMatrix m, Storable e) => Int -> m e -> (m e, m e)
splitColsMatrixAt j a
    | j < 0 || j > n = error $
        printf ("splitColsMatrixAt %d <matrix with dim (%d,%d)>:"
                ++ " invalid index") j m n
    | otherwise = let
        a1 = unsafeSliceMatrix (0,0) (m,j)   a
        a2 = unsafeSliceMatrix (0,j) (m,n-j) a
    in (a1,a2)
  where
    (m,n) = dimMatrix a
{-# INLINE splitColsMatrixAt #-}

-- | Get the indices of the elements in the matrix, in column-major order.
indicesMatrix :: (RMatrix m, Storable e) => m e -> [(Int,Int)]
indicesMatrix a = [ (i,j) | j <- [ 0..n-1 ], i <- [ 0..m-1 ] ]
  where (m,n) = dimMatrix a
  
-- | Lazily get the elements of the matrix, in column-major order.  
getElemsMatrix :: (RMatrix m, Storable e) => m e -> ST s [e]
getElemsMatrix a = case maybeWithVectorViewMatrix a V.getElems of
    Just es -> es
    Nothing -> withColsMatrix a $ \xs ->
                   concat `fmap` mapM V.getElems xs

-- | Get the elements of the matrix, in column-major order.
getElemsMatrix' :: (RMatrix m, Storable e) => m e -> ST s [e]
getElemsMatrix' a = case maybeWithVectorViewMatrix a V.getElems' of
    Just es -> es
    Nothing -> withColsMatrix a $ \xs ->
                   concat `fmap` mapM V.getElems' xs

-- | Lazily get the association list of the matrix, in column-major order.
getAssocsMatrix :: (RMatrix m, Storable e) => m e -> ST s [((Int,Int),e)]
getAssocsMatrix a = do
    es <- getElemsMatrix a
    return $ zip (indicesMatrix a) es

-- | Get the association list of the matrix, in column-major order.
getAssocsMatrix' :: (RMatrix m, Storable e) => m e -> ST s [((Int,Int),e)]
getAssocsMatrix' a = do
    es <- getElemsMatrix' a
    return $ zip (indicesMatrix a) es

-- | Set all of the values of the matrix from the elements in the list,
-- in column-major order.
setElemsMatrix :: (Storable e) => STMatrix s e -> [e] -> ST s ()
setElemsMatrix a es =
    case maybeWithVectorViewMatrixST a (`V.setElems` es) of
        Just st  -> st
        Nothing -> go 0 es
  where
    (m,n) = dimMatrix a
    go j [] | j == n = return ()
    go j [] | j < n = error $ 
        printf ("setElemsMatrix <matrix with dim (%d,%d>"
                ++ "<list with length %d>: not enough elements)") m n (j*m)
    go j es' =
        let (es1', es2') = splitAt m es'
        in do
            withColMatrixST a j (`V.setElems` es1')
            go (j+1) es2'

-- | Set the given values in the matrix.  If an index is repeated twice,
-- the value is implementation-defined.
setAssocsMatrix :: (Storable e) => STMatrix s e -> [((Int,Int),e)] -> ST s ()
setAssocsMatrix a ies =
    sequence_ [ writeMatrix a i e | (i,e) <- ies ]

unsafeSetAssocsMatrix :: (Storable e) => STMatrix s e -> [((Int,Int),e)] -> ST s ()
unsafeSetAssocsMatrix a ies =
    sequence_ [ unsafeWriteMatrix a i e | (i,e) <- ies ]

-- | Set the specified row of the matrix to the given vector.
setRowMatrix :: (RVector v, Storable e)
             => STMatrix s e -> Int -> v e -> ST s ()
setRowMatrix a i x
    | i < 0 || i >= m = error $
        printf ("setRowMatrix <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n i
    | V.dim x /= n = error $
        printf ("setRowMatrix <matrix with dim (%d,%d)> _"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim x)
    | otherwise =
        unsafeSetRowMatrix a i x
  where
    (m,n) = dimMatrix a
{-# INLINE setRowMatrix #-}

unsafeSetRowMatrix :: (RVector v, Storable e)
                   => STMatrix s e -> Int -> v e -> ST s ()
unsafeSetRowMatrix a i x = do
    jes <- V.getAssocs x
    sequence_ [ unsafeWriteMatrix a (i,j) e | (j,e) <- jes ]
{-# INLINE unsafeSetRowMatrix #-}

-- | Copy the specified row of the matrix to the vector.
getRowMatrix :: (RMatrix m, Storable e)
             => m e -> Int -> STVector s e -> ST s ()
getRowMatrix a i x
    | i < 0 || i >= m = error $
        printf ("getRowMatrix <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n i
    | V.dim x /= n = error $
        printf ("getRowMatrix <matrix with dim (%d,%d)> _"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim x)
    | otherwise =
        unsafeGetRowMatrix a i x
  where
    (m,n) = dimMatrix a
{-# INLINE getRowMatrix #-}

unsafeGetRowMatrix :: (RMatrix m, Storable e)
                   => m e -> Int -> STVector s e -> ST s ()
unsafeGetRowMatrix a i x = 
    forM_ [ 0..n-1 ] $ \j -> do
        e <- unsafeReadMatrix a (i,j)
        V.unsafeWrite x j e
  where
    (_,n) = dimMatrix a
{-# INLINE unsafeGetRowMatrix #-}

-- | Set the diagonal of the matrix to the given vector.
setDiagMatrix :: (RVector v, Storable e)
              => STMatrix s e -> v e -> ST s ()
setDiagMatrix a x
    | V.dim x /= mn = error $
        printf ("setRowMatrix <matrix with dim (%d,%d)>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim x)
    | otherwise =
        unsafeSetDiagMatrix a x
  where
    (m,n) = dimMatrix a
    mn = min m n
{-# INLINE setDiagMatrix #-}

unsafeSetDiagMatrix :: (RVector v, Storable e)
                    => STMatrix s e -> v e -> ST s ()
unsafeSetDiagMatrix a x = do
    ies <- V.getAssocs x
    sequence_ [ unsafeWriteMatrix a (i,i) e | (i,e) <- ies ]
{-# INLINE unsafeSetDiagMatrix #-}

-- | Copy the diagonal of the matrix to the vector.
getDiagMatrix :: (RMatrix m, Storable e)
             => m e -> STVector s e -> ST s ()
getDiagMatrix a x
    | V.dim x /= mn = error $
        printf ("getDiagMatrix <matrix with dim (%d,%d)>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim x)
    | otherwise =
        unsafeGetDiagMatrix a x
  where
    (m,n) = dimMatrix a
    mn = min m n
{-# INLINE getDiagMatrix #-}

unsafeGetDiagMatrix :: (RMatrix m, Storable e)
                    => m e -> STVector s e -> ST s ()
unsafeGetDiagMatrix a x = 
    forM_ [ 0..mn-1 ] $ \i -> do
        e <- unsafeReadMatrix a (i,i)
        V.unsafeWrite x i e
  where
    (m,n) = dimMatrix a
    mn = min m n
{-# INLINE unsafeGetDiagMatrix #-}

-- | Get the element stored at the given index.
readMatrix :: (RMatrix m, Storable e) => m e -> (Int,Int) -> ST s e
readMatrix a (i,j)
    | i < 0 || i >= m || j < 0 || j >= n = error $
        printf ("readMatrix <matrix with dim (%d,%d)> (%d,%d):"
                ++ " index out of range") m n i j
    | otherwise =
        unsafeReadMatrix a (i,j)
  where
    (m,n) = dimMatrix a
{-# INLINE readMatrix #-}

unsafeReadMatrix :: (RMatrix m, Storable e) => m e -> (Int,Int) -> ST s e
unsafeReadMatrix a (i,j) = unsafeIOToST $
    unsafeWithMatrix a $ \p lda ->
        peekElemOff p (i + j * lda)
{-# INLINE unsafeReadMatrix #-}

-- | Set the element stored at the given index.
writeMatrix :: (Storable e)
            => STMatrix s e -> (Int,Int) -> e -> ST s ()
writeMatrix a (i,j)
    | i < 0 || i >= m || j < 0 || j >= n = error $
        printf ("writeMatrix <matrix with dim (%d,%d)> (%d,%d):"
                ++ " index out of range") m n i j
    | otherwise =
        unsafeWriteMatrix a (i,j)
  where
    (m,n) = dimMatrix a
{-# INLINE writeMatrix #-}

unsafeWriteMatrix :: (Storable e)
                  => STMatrix s e -> (Int,Int) -> e -> ST s ()
unsafeWriteMatrix a (i,j) e = unsafeIOToST $
    unsafeWithMatrix a $ \p lda ->
        pokeElemOff p (i + j * lda) e
{-# INLINE unsafeWriteMatrix #-}

-- | Modify the element stored at the given index.
modifyMatrix :: (Storable e)
             => STMatrix s e -> (Int,Int) -> (e -> e) -> ST s ()
modifyMatrix a (i,j)
    | i < 0 || i >= m || j < 0 || j >= n = error $
        printf ("modifyMatrix <matrix with dim (%d,%d)> (%d,%d):"
                ++ " index out of range") m n i j
    | otherwise =
        unsafeModifyMatrix a (i,j)
  where
    (m,n) = dimMatrix a
{-# INLINE modifyMatrix #-}

unsafeModifyMatrix :: (Storable e)
                   => STMatrix s e -> (Int,Int) -> (e -> e) -> ST s ()
unsafeModifyMatrix a (i,j) f = unsafeIOToST $
    unsafeWithMatrix a $ \p lda -> 
        let o = i + j * lda
        in do
            e <- peekElemOff p o
            pokeElemOff p o $ f e
{-# INLINE unsafeModifyMatrix #-}

-- | @mapToMatrix f a c@ replaces @c@ elementwise with @f(a)@.
mapToMatrix :: (RMatrix m, Storable e, Storable f)
            => (e -> f)
            -> m e
            -> STMatrix s f
            -> ST s ()
mapToMatrix f = checkMatrixOp2 "mapToMatrix _" $ unsafeMapToMatrix f
{-# INLINE mapToMatrix #-}
                            
unsafeMapToMatrix :: (RMatrix m, Storable e, Storable f)
                  => (e -> f)
                  -> m e
                  -> STMatrix s f
                  -> ST s ()
unsafeMapToMatrix f a c =
    fromMaybe colwise $ maybeWithVectorViewMatrix a $ \x ->
    fromMaybe colwise $ maybeWithVectorViewMatrixST c $ \z ->
        V.unsafeMapTo f x z
  where
    colwise = withColsMatrix   a $ \xs ->
              withColsMatrixST c $ \zs ->
                  sequence_ [ V.unsafeMapTo f x z
                            | (x,z) <- zip xs zs
                            ]

-- | @zipWithToMatrix f a b c@ replaces @c@ elementwise with @f(a, b)@.
zipWithToMatrix :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
                => (e1 -> e2 -> f)
                -> m1 e1
                -> m2 e2
                -> STMatrix s f
                -> ST s ()
zipWithToMatrix f = checkMatrixOp3 "zipWithToMatrix _" $
    unsafeZipWithToMatrix f
{-# INLINE zipWithToMatrix #-}

unsafeZipWithToMatrix :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
                      => (e1 -> e2 -> f)
                      -> m1 e1
                      -> m2 e2
                      -> STMatrix s f
                      -> ST s ()
unsafeZipWithToMatrix f a b c =
    fromMaybe colwise $ maybeWithVectorViewMatrix a $ \x ->
    fromMaybe colwise $ maybeWithVectorViewMatrix b $ \y ->
    fromMaybe colwise $ maybeWithVectorViewMatrixST c $ \z ->
        V.unsafeZipWithTo f x y z
  where
    colwise = withColsMatrix   a $ \xs ->
              withColsMatrix   b $ \ys ->
              withColsMatrixST c $ \zs ->
                  sequence_ [ V.unsafeZipWithTo f x y z
                            | (x,y,z) <- zip3 xs ys zs
                            ]

-- | Set every element in the matrix to a default value.  For
-- standard numeric types (including 'Double', 'Complex Double', and 'Int'),
-- the default value is '0'.
clearMatrix :: (Storable e) => STMatrix s e -> ST s ()
clearMatrix a = case maybeVectorViewMatrixST a of
    Just x  -> V.clear x
    Nothing -> withColsMatrixST a $ mapM_ V.clear

-- | Add a constant to all entries of a matrix.
shiftToMatrix :: (RMatrix m, VNum e)
              => e -> m e -> STMatrix s e -> ST s ()
shiftToMatrix e = checkMatrixOp2 "shiftToMatrix" $
    matrixVectorOp2 (V.shiftTo e)

-- | Add a vector to the diagonal of a matrix.
shiftDiagToMatrix :: (RVector v, RMatrix m, BLAS1 e)
                  => v e -> m e -> STMatrix s e -> ST s ()
shiftDiagToMatrix s a b
    | V.dim s /= mn || dimMatrix a /= (m,n) = error $
        printf ("shiftDiagToMatrix <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (V.dim s)
                (fst $ dimMatrix a) (snd $ dimMatrix a) m n
    | otherwise = shiftDiagToMatrixWithScale 1 s a b
  where
    (m,n) = dimMatrix b
    mn = min m n

-- | Add a scaled vector to the diagonal of a matrix.
shiftDiagToMatrixWithScale :: (RVector v, RMatrix m, BLAS1 e)
                           => e -> v e -> m e -> STMatrix s e -> ST s ()
shiftDiagToMatrixWithScale e s a b
    | V.dim s /= mn || dimMatrix a /= (m,n) = error $
        printf ("shiftDiagToMatrixWithScale _ <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (V.dim s)
                (fst $ dimMatrix a) (snd $ dimMatrix a) m n
    | otherwise = do
        unsafeCopyToMatrix a b
        unsafeIOToST $
            V.unsafeWith s $ \ps ->
            unsafeWithMatrix b $ \pb ldb ->
                BLAS.axpy mn e ps 1 pb (ldb+1)
  where
    (m,n) = dimMatrix b
    mn = min m n

-- | Add two matrices.
addToMatrix :: (RMatrix m1, RMatrix m2, VNum e)
            => m1 e -> m2 e -> STMatrix s e -> ST s ()
addToMatrix = checkMatrixOp3 "addToMatrix" $ matrixVectorOp3 V.addTo

-- | Add two matrices with the given scales.
addToMatrixWithScales :: (RMatrix m1, RMatrix m2, VNum e)
                      => e -> m1 e -> e -> m2 e -> STMatrix s e -> ST s ()
addToMatrixWithScales alpha a beta b c =
    (checkMatrixOp3 "addToMatrixWithScales" $
        matrixVectorOp3 (\x y z -> V.addToWithScales alpha x beta y z))
        a b c

-- | Subtract two matrices.
subToMatrix :: (RMatrix m1, RMatrix m2, VNum e)
            => m1 e -> m2 e -> STMatrix s e -> ST s ()
subToMatrix = checkMatrixOp3 "subToMatrix" $ matrixVectorOp3 V.subTo

-- | Conjugate the entries of a matrix.
conjToMatrix :: (RMatrix m, VNum e)
             => m e -> STMatrix s e -> ST s ()
conjToMatrix = checkMatrixOp2 "conjToMatrix" $
    matrixVectorOp2 V.conjTo

-- | Negate the entries of a matrix.
negateToMatrix :: (RMatrix m, VNum e)
               => m e -> STMatrix s e -> ST s ()
negateToMatrix = checkMatrixOp2 "negateToMatrix" $
    matrixVectorOp2 V.negateTo

-- | Scale the entries of a matrix by the given value.
scaleToMatrix :: (RMatrix m, VNum e)
              => e -> m e -> STMatrix s e -> ST s ()
scaleToMatrix e = checkMatrixOp2 "scaleToMatrix" $
    matrixVectorOp2 (V.scaleTo e)

-- | Scale the rows of a matrix; @scaleRowsToMatrix s a c@ sets
-- @c := diag(s) * a@.
scaleRowsToMatrix :: (RVector v, RMatrix m, VNum e)
                  => v e -> m e -> STMatrix s e -> ST s ()
scaleRowsToMatrix s a b
    | V.dim s /= m || dimMatrix a /= (m,n) = error $
        printf ("scaleRowsToMatrix <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (V.dim s)
                (fst $ dimMatrix a) (snd $ dimMatrix a) m n
    | otherwise =
        withColsMatrix   a $ \xs ->
        withColsMatrixST b $ \ys ->
            zipWithM_ (V.mulTo s) xs ys
  where
    (m,n) = dimMatrix b

-- | Scale the columns of a matrix; @scaleColsToMatrix s a c@ sets
-- @c := a * diag(s)@.
scaleColsToMatrix :: (RVector v, RMatrix m, VNum e)
                  => v e -> m e -> STMatrix s e -> ST s ()
scaleColsToMatrix s a b 
    | V.dim s /= n || dimMatrix a /= (m,n) = error $
        printf ("scaleColsToMatrix <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (V.dim s)
                (fst $ dimMatrix a) (snd $ dimMatrix a) m n
    | otherwise =
        V.getElems   s >>= \es ->
        withColsMatrix   a $ \xs ->
        withColsMatrixST b $ \ys ->
            sequence_ [ V.scaleTo e x y
                      | (e,x,y) <- zip3 es xs ys
                      ]
  where
    (m,n) = dimMatrix b

-- | @rank1UpdateToMatrix alpha x y a@ sets @a := alpha * x * y^H + a@.
rank1UpdateToMatrix :: (RVector v1, RVector v2, BLAS2 e)
                    => e -> v1 e -> v2 e -> STMatrix s e -> ST s ()
rank1UpdateToMatrix alpha x y a
    | V.dim x /= m || V.dim y /= n = error $
        printf ("rank1UpdateToMatrix _ <vector with dim %d>"
                ++ " <vector with dim %d> <matrix with dim (%d,%d)>:"
                ++ " dimension mismatch")
                (V.dim x) (V.dim y) m n
    | otherwise = do
        unsafeIOToST $
            V.unsafeWith x $ \px ->
            V.unsafeWith y $ \py ->
            unsafeWithMatrix a $ \pa lda ->
                BLAS.gerc m n alpha px 1 py 1 pa lda
  where
    (m,n) = dimMatrix a


-- | @transToMatrix a c@ sets @c := trans(a)@.
transToMatrix :: (RMatrix m, BLAS1 e)
              => m e
              -> STMatrix s e
              -> ST s ()
transToMatrix a a'
    | (ma,na) /= (na',ma') = error $
        printf ("transToMatrix <matrix with dim (%d,%d)>"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
    | otherwise = unsafeIOToST $
        unsafeWithMatrix a $ \pa lda ->
        unsafeWithMatrix a' $ \pa' lda' -> let
            go j px py | j == n = return ()
                       | otherwise = do
                           BLAS.copy m px 1 py lda'
                           go (j+1) (px `advancePtr` lda) (py `advancePtr` 1)
            in go 0 pa pa'
  where
    (ma,na) = dimMatrix a
    (ma',na') = dimMatrix a'
    (m,n) = (ma,na)

-- | @conjTransToMatrix a c@ sets @c := conj(trans(a))@.
conjTransToMatrix :: (RMatrix m, BLAS1 e)
                  => m e
                  -> STMatrix s e
                  -> ST s ()
conjTransToMatrix a a' = do
    transToMatrix a a'
    conjToMatrix a' a'

-- | @mulMatrixToVector transa a x y@
-- sets @y := op(a) * x@, where @op(a)@ is determined by @transa@.                   
mulMatrixToVector :: (RMatrix m, RVector v, BLAS2 e)
                  => Trans -> m e
                  -> v e
                  -> STVector s e
                  -> ST s ()
mulMatrixToVector = mulMatrixToVectorWithScale 1

-- | @mulMatrixToVectorWithScale alpha transa a x y@
-- sets @y := alpha * op(a) * x@, where @op(a)@ is determined by @transa@.                   
mulMatrixToVectorWithScale :: (RMatrix m, RVector v, BLAS2 e)
                           => e
                           -> Trans -> m e
                           -> v e
                           -> STVector s e
                           -> ST s ()
mulMatrixToVectorWithScale alpha t a x y =
    mulMatrixAddToVectorWithScales alpha t a x 0 y

-- | @mulMatrixAddToVectorWithScales alpha transa a x beta y@
-- sets @y := alpha * op(a) * x + beta * y@, where @op(a)@ is
-- determined by @transa@.
mulMatrixAddToVectorWithScales :: (RMatrix m, RVector v, BLAS2 e)
                               => e
                               -> Trans -> m e
                               -> v e
                               -> e
                               -> STVector s e
                               -> ST s ()
mulMatrixAddToVectorWithScales alpha transa a x beta y
    | (not . and) [ case transa of NoTrans -> (ma,na) == (m,n)
                                   _       -> (ma,na) == (n,m)
                  , nx == n
                  , ny == m
                  ] = error $
        printf ("mulMatrixAddToVectorWithScales _"
                ++ " %s <matrix with dim (%d,%d)>" 
                ++ " <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>: dimension mismatch")
               (show transa) ma na
               nx ny
    | otherwise =
        unsafeIOToST $
            unsafeWithMatrix a $ \pa lda ->
            V.unsafeWith x $ \px ->
            V.unsafeWith y $ \py ->
                if n == 0
                    then BLAS.scal m beta py 1
                    else BLAS.gemv transa ma na alpha pa lda px 1 beta py 1
  where
    (ma,na) = dimMatrix a
    nx = V.dim x
    ny = V.dim y
    (m,n) = (ny,nx)

-- | @mulMatrixToMatrix transa a transb b c@
-- sets @c := op(a) * op(b)@, where @op(a)@ and @op(b)@ are determined
-- by @transa@ and @transb@.                   
mulMatrixToMatrix :: (RMatrix m1, RMatrix m2, BLAS3 e)
                  => Trans -> m1 e
                  -> Trans -> m2 e
                  -> STMatrix s e
                  -> ST s ()
mulMatrixToMatrix = mulMatrixToMatrixWithScale 1

-- | @mulMatrixToMatrixWithScale alpha transa a transb b c@
-- sets @c := alpha * op(a) * op(b)@, where @op(a)@ and @op(b)@ are determined
-- by @transa@ and @transb@.                   
mulMatrixToMatrixWithScale :: (RMatrix m1, RMatrix m2, BLAS3 e)
                           => e
                           -> Trans -> m1 e
                           -> Trans -> m2 e
                           -> STMatrix s e
                           -> ST s ()
mulMatrixToMatrixWithScale alpha ta a tb b c =
    mulMatrixAddToMatrixWithScales alpha ta a tb b 0 c

-- | @mulMatrixAddToMatrixWithScales alpha transa a transb b beta c@
-- sets @c := alpha * op(a) * op(b) + beta * c@, where @op(a)@ and
-- @op(b)@ are determined by @transa@ and @transb@.
mulMatrixAddToMatrixWithScales :: (RMatrix m1, RMatrix m2, BLAS3 e)
                               => e
                               -> Trans -> m1 e
                               -> Trans -> m2 e
                               -> e
                               -> STMatrix s e
                               -> ST s ()
mulMatrixAddToMatrixWithScales alpha transa a transb b beta c
    | (not . and) [ case transa of NoTrans -> (ma,na) == (m,k)
                                   _       -> (ma,na) == (k,m)
                  , case transb of NoTrans -> (mb,nb) == (k,n)
                                   _       -> (mb,nb) == (n,k)
                  , (mc, nc) == (m,n)
                  ] = error $
        printf ("mulMatrixAddToMatrixWithScales _"
                ++ " %s <matrix with dim (%d,%d)>" 
                ++ " %s <matrix with dim (%d,%d)>"
                ++ " _"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
               (show transa) ma na
               (show transb) mb nb
               mc nc
    | otherwise =
        unsafeIOToST $
            unsafeWithMatrix a $ \pa lda ->
            unsafeWithMatrix b $ \pb ldb ->
            unsafeWithMatrix c $ \pc ldc ->
                BLAS.gemm transa transb m n k alpha pa lda pb ldb beta pc ldc
  where
    (ma,na) = dimMatrix a
    (mb,nb) = dimMatrix b
    (mc,nc) = dimMatrix c
    (m,n) = dimMatrix c
    k = case transa of NoTrans -> na
                       _       -> ma


checkMatrixOp2 :: (RMatrix x, RMatrix y, Storable e, Storable f)
               => String
               -> (x e -> y f -> a)
               -> x e
               -> y f
               -> a
checkMatrixOp2 str f x y
    | (m1,n1) /= (m2,n2) = error $
        printf ("%s <matrix with dim (%d,%d)> <matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") str m1 n1 m2 n2
    | otherwise =
        f x y
  where
    (m1,n1) = dimMatrix x
    (m2,n2) = dimMatrix y        
{-# INLINE checkMatrixOp2 #-}

checkMatrixOp3 :: (RMatrix x, RMatrix y, RMatrix z, Storable e, Storable f, Storable g)
               => String
               -> (x e -> y f -> z g -> a)
               -> x e
               -> y f
               -> z g
               -> a
checkMatrixOp3 str f x y z
    | (m1,n1) /= (m2,n2) || (m1,n1) /= (m3,n3) = error $
        printf ("%s <matrix with dim (%d,%d)> <matrix with dim (%d,%d)>:"
                ++ " <matrix with dim (%d,%d)> dimension mismatch")
               str m1 n1 m2 n2 m3 n3
    | otherwise =
        f x y z
  where
    (m1,n1) = dimMatrix x
    (m2,n2) = dimMatrix y
    (m3,n3) = dimMatrix z
{-# INLINE checkMatrixOp3 #-}

matrixVectorOp2 :: (RMatrix m, Storable e, Storable f)
                => (forall v . RVector v => v e -> STVector s f -> ST s ())
                -> m e -> STMatrix s f -> ST s ()
matrixVectorOp2 f a c =
    fromMaybe colwise $ maybeWithVectorViewMatrix   a $ \x ->
    fromMaybe colwise $ maybeWithVectorViewMatrixST c $ \z ->
        f x z
  where
    colwise = withColsMatrix   a $ \xs ->
              withColsMatrixST c $ \zs ->
                  sequence_ [ f x z | (x,z) <- zip xs zs ]
{-# INLINE matrixVectorOp2 #-}

matrixVectorOp3 :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
                => (forall v1 v2 . (RVector v1, RVector v2) => 
                    v1 e1 -> v2 e2 -> STVector s f -> ST s ())
                -> m1 e1 -> m2 e2 -> STMatrix s f -> ST s ()
matrixVectorOp3 f a b c =
    fromMaybe colwise $ maybeWithVectorViewMatrix   a $ \x ->
    fromMaybe colwise $ maybeWithVectorViewMatrix   b $ \y ->
    fromMaybe colwise $ maybeWithVectorViewMatrixST c $ \z ->
        f x y z
  where
    colwise = withColsMatrix   a $ \xs ->
              withColsMatrix   b $ \ys ->
              withColsMatrixST c $ \zs ->
                  sequence_ [ f x y z | (x,y,z) <- zip3 xs ys zs ]
{-# INLINE matrixVectorOp3 #-}

newResultMatrix :: (RMatrix m, Storable e, Storable f)
                => (m e -> STMatrix s f -> ST s a)
                -> m e
                -> ST s (STMatrix s f)
newResultMatrix f a = do
    c <- newMatrix_ (dimMatrix a)
    _ <- f a c
    return c
{-# INLINE newResultMatrix #-}

newResultMatrix2 :: (RMatrix m1, RMatrix m2, Storable e, Storable f, Storable g)
                 => (m1 e -> m2 f -> STMatrix s g -> ST s a)
                 -> m1 e
                 -> m2 f
                 -> ST s (STMatrix s g)
newResultMatrix2 f a1 a2 = do
    c <- newMatrix_ (dimMatrix a1)
    _ <- f a1 a2 c
    return c
{-# INLINE newResultMatrix2 #-}
