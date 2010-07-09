{-# LANGUAGE DeriveDataTypeable, TypeFamilies #-}
{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Matrix.STBase
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : asinerimental
--

module BLAS.Matrix.STBase
    where
      
import Control.Monad( forM_ )
import Control.Monad.ST( ST, unsafeIOToST )
import Data.Typeable( Typeable )
import Foreign( ForeignPtr, Ptr, peekElemOff, pokeElemOff )
import Text.Printf( printf )

import BLAS.Elem
import BLAS.Vector.STBase

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
  deriving (Eq, Typeable)

class HasVectorView m where
    type VectorView m :: * -> *

instance HasVectorView (STMatrix s) where
    type VectorView (STMatrix s) = STVector s

-- | Read-only matrices
class (HasVectorView m, RVector (VectorView m)) => RMatrix m where
    -- | The dimensions of the matrix (number of rows and columns).
    dimMatrix :: m e -> (Int,Int)
    
    unsafeSpliceMatrix :: (Storable e)
                       => m e -> (Int,Int) -> (Int,Int) -> m e

    -- | Execute an 'IO' action with a pointer to the first element in the
    -- matrix and the leading dimension (lda).
    unsafeWithMatrix :: m e -> (Ptr e -> Int -> IO a) -> IO a

    unsafeVectorToMatrix :: (Int,Int) -> VectorView m e -> m e
    maybeMatrixToVector :: m e -> Maybe (VectorView m e)

    unsafeColMatrix :: (Storable e)
                    => m e
                    -> Int
                    -> VectorView m e

instance RMatrix (STMatrix s) where
    dimMatrix (STMatrix _ m n _) = (m,n)
    {-# INLINE dimMatrix #-}

    unsafeSpliceMatrix (STMatrix a _ _ lda) (i,j) (m',n') = let
        o = i + j*lda
        a' = unsafeSpliceVector a o (dimVector a - o)
        in STMatrix a' m' n' lda
    {-# INLINE unsafeSpliceMatrix #-}
    
    unsafeWithMatrix (STMatrix a _ _ lda) f =
        unsafeWithVector a $ \p -> f p lda
    {-# INLINE unsafeWithMatrix #-}

    unsafeVectorToMatrix (m,n) x = (STMatrix x m n lda)
        where lda = max 1 m
    {-# INLINE unsafeVectorToMatrix #-}

    maybeMatrixToVector (STMatrix a m n lda) | lda == max 1 m = Just a
                                             | n == 1         = Just a
                                             | otherwise      = Nothing
    {-# INLINE maybeMatrixToVector #-}

    unsafeColMatrix (STMatrix a m _ lda) j = let
        o = j * lda
        x = unsafeSpliceVector a o m
        in x
    {-# INLINE unsafeColMatrix #-}


vectorToMatrix :: (RMatrix m) => (Int,Int) -> VectorView m e -> m e
vectorToMatrix mn@(m,n) x 
    | dimVector x /= m*n = error $
        printf ("vectorToMatrix (%d,%d) <vector with dim %d>:"
                ++ " incorrect number of elements") m n (dimVector x)
    | otherwise =
        unsafeVectorToMatrix mn x
{-# INLINE vectorToMatrix #-}

colMatrix :: (RMatrix m, Storable e)
          => m e
          -> Int
          -> VectorView m e
colMatrix a j
    | j < 0 || j >= n = error $
        printf ("colMatrix <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n j
    | otherwise =
        unsafeColMatrix a j
  where
    (m,n) = dimMatrix a
{-# INLINE colMatrix #-}

-- | @spliceMatrix a (i,j) (m,n)@ creates a submatrix view of @a@ starting at
-- element @(i,j)@ and having dimensions @(m,n)@.
spliceMatrix :: (RMatrix m, Storable e)
             => m e
             -> (Int,Int)
             -> (Int,Int)
             -> m e
spliceMatrix a (i,j) (m',n')
    | (i < 0 || m' < 0 || i + m' > m 
       || j < 0 || n' < 0 || j + n' > n) = error $
        printf ("spliceMatrix <matrix with dim (%d,%d)> (%d,%d) (%d,%d):"
                ++ " index out of range") m n i j m' n'
    | otherwise =
        unsafeSpliceMatrix a (i,j) (m',n')
  where
    (m,n) = dimMatrix a
{-# INLINE spliceMatrix #-}

-- | View an array in memory as a matrix.
matrixViewArray :: (Storable e)
                => ForeignPtr e 
                -> Int          -- ^ offset
                -> (Int,Int)    -- ^ dimension
                -> STMatrix s e
matrixViewArray f o (m,n) = let
    x   = vectorViewArray f o (m*n)
    lda = max 1 m
    in STMatrix x m n lda
{-# INLINE matrixViewArray #-}

-- | Create a new matrix of given shape, but do not initialize the elements.
newMatrix_ :: (Storable e) => (Int,Int) -> ST s (STMatrix s e)
newMatrix_ (m,n) 
    | m < 0 || n < 0 = error $
        printf "newMatrix_ (%d,%d): invalid dimensions" m n
    | otherwise =  do
        x <- newVector_ (m*n)
        return $ unsafeVectorToMatrix (m,n) x

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
unsafeCopyToMatrix a b = 
    case (maybeMatrixToVector a, maybeMatrixToVector b) of
        (Just x, Just y) -> unsafeCopyToVector x y
        _ ->
            forM_ [ 0..n-1 ] $ \j ->
                unsafeCopyToVector (colMatrix a j) (colMatrix b j)
  where
    n = (snd . dimMatrix) b
{-# INLINE unsafeCopyToMatrix #-}

-- | Split a matrix into two blocks and returns views into the blocks.  In
-- @(a1, a2) = splitRowsMatrixAt i a@, we have that
-- @a1 = spliceMatrix a (0,0) (i,n)@ and
-- @a2 = spliceMatrix a (i,0) (m-i,n)@, where @(m,n)@ is the dimension of @a@.
splitRowsMatrixAt :: (RMatrix m, Storable e) => Int -> m e -> (m e, m e)
splitRowsMatrixAt i a
    | i < 0 || i > m = error $
        printf ("splitRowsMatrixAt %d <matrix with dim (%d,%d)>:"
                ++ " invalid index") i m n
    | otherwise = let
        a1 = unsafeSpliceMatrix a (0,0) (i,n)
        a2 = unsafeSpliceMatrix a (i,0) (m-i,n)
    in (a1,a2)
  where
    (m,n) = dimMatrix a
{-# INLINE splitRowsMatrixAt #-}

-- | Split a matrix into two blocks and returns views into the blocks.  In
-- @(a1, a2) = splitColsMatrixAt j a@, we have that
-- @a1 = spliceMatrix a (0,0) (m,j)@ and
-- @a2 = spliceMatrix a (0,j) (m,n-j)@, where @(m,n)@ is the dimension of @a@.
splitColsMatrixAt :: (RMatrix m, Storable e) => Int -> m e -> (m e, m e)
splitColsMatrixAt j a
    | j < 0 || j > n = error $
        printf ("splitColsMatrixAt %d <matrix with dim (%d,%d)>:"
                ++ " invalid index") j m n
    | otherwise = let
        a1 = unsafeSpliceMatrix a (0,0) (m,j)
        a2 = unsafeSpliceMatrix a (0,j) (m,n-j)
    in (a1,a2)
  where
    (m,n) = dimMatrix a
{-# INLINE splitColsMatrixAt #-}

-- | Get the indices of the elements in the matrix, in column-major order.
indicesMatrix :: (RMatrix m) => m e -> [(Int,Int)]
indicesMatrix a = [ (i,j) | j <- [ 0..n-1 ], i <- [ 0..m-1 ] ]
  where (m,n) = dimMatrix a
  
-- | Lazily get the elements of the matrix, in column-major order.  
getElemsMatrix :: (RMatrix m, Storable e) => m e -> ST s [e]
getElemsMatrix a = case maybeMatrixToVector a of
    Just x -> getElemsVector x
    Nothing -> concat `fmap` sequence [ getElemsVector (colMatrix a j)
                                      | j <- [ 0..n-1 ]
                                      ]
  where
    n = (snd . dimMatrix) a

-- | Get the elements of the matrix, in column-major order.
getElemsMatrix' :: (RMatrix m, Storable e) => m e -> ST s [e]
getElemsMatrix' a = case maybeMatrixToVector a of
    Just x -> getElemsVector' x
    Nothing -> concat `fmap` sequence [ getElemsVector' (colMatrix a j)
                                      | j <- [ 0..n-1 ]
                                      ]
  where
    n = (snd . dimMatrix) a

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
setElemsMatrix a es = case maybeMatrixToVector a of
    Just x  -> setElemsVector x es
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
            setElemsVector (colMatrix a j) es1'
            go (j+1) es2'

-- | Set the given values in the matrix.  If an index is repeated twice,
-- the value is implementation-defined.
setAssocsMatrix :: (Storable e) => STMatrix s e -> [((Int,Int),e)] -> ST s ()
setAssocsMatrix a ies =
    sequence_ [ writeMatrix a i e | (i,e) <- ies ]

unsafeSetAssocsMatrix :: (Storable e) => STMatrix s e -> [((Int,Int),e)] -> ST s ()
unsafeSetAssocsMatrix a ies =
    sequence_ [ unsafeWriteMatrix a i e | (i,e) <- ies ]

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

unsafeReadMatrix :: (RMatrix m, Storable e) => m e -> (Int,Int) -> ST s e
unsafeReadMatrix a (i,j) = unsafeIOToST $
    unsafeWithMatrix a $ \p lda ->
        peekElemOff p (i + j * lda)

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

unsafeWriteMatrix :: (Storable e)
                  => STMatrix s e -> (Int,Int) -> e -> ST s ()
unsafeWriteMatrix a (i,j) e = unsafeIOToST $
    unsafeWithMatrix a $ \p lda ->
        pokeElemOff p (i + j * lda) e


{-
add,
shift,
scale,
rank1 update
scaleRows
scaleCols
rowSums
colSums
rowMeans
colMeans
transpose/herm
-}

checkMatrixOp2 :: (RMatrix x, RMatrix y)
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
