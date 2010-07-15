{-# LANGUAGE DeriveDataTypeable, FlexibleContexts, TypeFamilies #-}
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
import Control.Monad.ST( ST, unsafeIOToST )
import Data.Typeable( Typeable )
import Foreign( ForeignPtr, Ptr, peekElemOff, pokeElemOff )
import Text.Printf( printf )

import Numeric.LinearAlgebra.Elem
import qualified Numeric.LinearAlgebra.Elem.BLAS1 as BLAS
import qualified Numeric.LinearAlgebra.Elem.BLAS2 as BLAS
import Numeric.LinearAlgebra.Types( HasVectorView(..) )
import Numeric.LinearAlgebra.Vector.STBase

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

    unsafeMatrixViewVector :: VectorView m e -> (Int,Int) -> m e
    
    -- | Possibly view a matrix as a vector.  This only succeeds if the
    -- matrix is stored contiguously in memory, i.e. if the matrix contains
    -- a single column or the "lda" of the matrix is equal to the number
    -- of rows.
    maybeVectorViewMatrix :: (Storable e) => m e -> Maybe (VectorView m e)

    unsafeColMatrix :: (Storable e)
                    => m e
                    -> Int
                    -> VectorView m e

instance RMatrix (STMatrix s) where
    dimMatrix (STMatrix _ m n _) = (m,n)
    {-# INLINE dimMatrix #-}

    unsafeSpliceMatrix (STMatrix a _ _ lda) (i,j) (m',n') = let
        o = i + j*lda
        l = if m' == 0 || n' == 0
                then 0
                else lda * (n' - 1) + m'
        a' = unsafeSpliceVector a o l
        in STMatrix a' m' n' lda
    {-# INLINE unsafeSpliceMatrix #-}
    
    unsafeWithMatrix (STMatrix a _ _ lda) f =
        unsafeWithVector a $ \p -> f p lda
    {-# INLINE unsafeWithMatrix #-}

    unsafeMatrixViewVector x (m,n) = (STMatrix x m n lda)
        where lda = max 1 m
    {-# INLINE unsafeMatrixViewVector #-}

    maybeVectorViewMatrix (STMatrix a m n lda) | lda == max 1 m = Just a
                                               | n == 1 = Just a
                                               | otherwise = Nothing
    {-# INLINE maybeVectorViewMatrix #-}

    unsafeColMatrix (STMatrix a m _ lda) j = let
        o = j * lda
        x = unsafeSpliceVector a o m
        in x
    {-# INLINE unsafeColMatrix #-}

-- | View a vector as a matrix with the given dimensions.
matrixViewVector :: (RMatrix m, Storable e)
                 => VectorView m e -> (Int,Int) -> m e
matrixViewVector x mn@(m,n)
    | dimVector x /= m*n = error $
        printf ("matrixViewVector <vector with dim %d> (%d,%d):"
                ++ " incorrect number of elements") (dimVector x) m n
    | otherwise =
        unsafeMatrixViewVector x mn
{-# INLINE matrixViewVector #-}

-- | Get a view of a matrix column.
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

-- | Get a list of views of the matrix columns.
colsMatrix :: (RMatrix m, Storable e)
           => m e
           -> [VectorView m e]
colsMatrix a = [ unsafeColMatrix a j | j <- [ 0..n-1 ] ]
  where n = (snd . dimMatrix) a

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
        return $ unsafeMatrixViewVector x (m,n)

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
unsafeCopyToMatrix = matrixVectorOp2 unsafeCopyToVector
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
getElemsMatrix a = case maybeVectorViewMatrix a of
    Just x -> getElemsVector x
    Nothing -> concat `fmap` sequence [ getElemsVector (colMatrix a j)
                                      | j <- [ 0..n-1 ]
                                      ]
  where
    n = (snd . dimMatrix) a

-- | Get the elements of the matrix, in column-major order.
getElemsMatrix' :: (RMatrix m, Storable e) => m e -> ST s [e]
getElemsMatrix' a = case maybeVectorViewMatrix a of
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
setElemsMatrix a es = case maybeVectorViewMatrix a of
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

-- | Set the specified row of the matrix to the given vector.
setRowMatrix :: (RVector v, Storable e)
             => STMatrix s e -> Int -> v e -> ST s ()
setRowMatrix a i x
    | i < 0 || i >= m = error $
        printf ("setRowMatrix <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n i
    | dimVector x /= n = error $
        printf ("setRowMatrix <matrix with dim (%d,%d)> _"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (dimVector x)
    | otherwise =
        unsafeSetRowMatrix a i x
  where
    (m,n) = dimMatrix a
{-# INLINE setRowMatrix #-}

unsafeSetRowMatrix :: (RVector v, Storable e)
                   => STMatrix s e -> Int -> v e -> ST s ()
unsafeSetRowMatrix a i x = do
    jes <- getAssocsVector x
    sequence_ [ unsafeWriteMatrix a (i,j) e | (j,e) <- jes ]
{-# INLINE unsafeSetRowMatrix #-}

-- | Copy the specified row of the matrix to the vector.
getRowMatrix :: (RMatrix m, Storable e)
             => m e -> Int -> STVector s e -> ST s ()
getRowMatrix a i x
    | i < 0 || i >= m = error $
        printf ("getRowMatrix <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n i
    | dimVector x /= n = error $
        printf ("getRowMatrix <matrix with dim (%d,%d)> _"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (dimVector x)
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
        unsafeWriteVector x j e
  where
    (_,n) = dimMatrix a
{-# INLINE unsafeGetRowMatrix #-}

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
    case (maybeVectorViewMatrix a, maybeVectorViewMatrix c) of
        (Just x, Just z) -> unsafeMapToVector f x z
        _ -> sequence_ [ unsafeMapToVector f x z 
                       | (x,z) <- zip (colsMatrix a) (colsMatrix c)
                       ]
{-# INLINE unsafeMapToMatrix #-}

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
    case (maybeVectorViewMatrix a, maybeVectorViewMatrix b,
          maybeVectorViewMatrix c) of
        (Just x, Just y, Just z) -> unsafeZipWithToVector f x y z
        _ -> sequence_ [ unsafeZipWithToVector f x y z
                       | (x,y,z) <- zip3 (colsMatrix a) (colsMatrix b) (colsMatrix c)
                       ]
{-# INLINE unsafeZipWithToMatrix #-}

-- | Set every element in the matrix to a default value.  For
-- standard numeric types (including 'Double', 'Complex Double', and 'Int'),
-- the default value is '0'.
clearMatrix :: (Storable e) => STMatrix s e -> ST s ()
clearMatrix a = case maybeVectorViewMatrix a of
    Just x -> clearVector x
    Nothing -> sequence_ [ clearVector x | x <- colsMatrix a ]

-- | Add a constant to all entries of a matrix.
shiftToMatrix :: (RMatrix m, VNum e)
              => e -> m e -> STMatrix s e -> ST s ()
shiftToMatrix e = checkMatrixOp2 "shiftToMatrix" $
    matrixVectorOp2 (shiftToVector e)

-- | Add a vector to the diagonal of a matrix.
shiftDiagToMatrix :: (RVector v, RMatrix m, BLAS1 e)
                  => v e -> m e -> STMatrix s e -> ST s ()
shiftDiagToMatrix s a b
    | dimVector s /= mn || dimMatrix a /= (m,n) = error $
        printf ("shiftDiagToMatrix <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (dimVector s)
                (fst $ dimMatrix a) (snd $ dimMatrix a) m n
    | otherwise = shiftDiagToMatrixWithScale 1 s a b
  where
    (m,n) = dimMatrix b
    mn = min m n

-- | Add a scaled vector to the diagonal of a matrix.
shiftDiagToMatrixWithScale :: (RVector v, RMatrix m, BLAS1 e)
                           => e -> v e -> m e -> STMatrix s e -> ST s ()
shiftDiagToMatrixWithScale e s a b
    | dimVector s /= mn || dimMatrix a /= (m,n) = error $
        printf ("shiftDiagToMatrixWithScale _ <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (dimVector s)
                (fst $ dimMatrix a) (snd $ dimMatrix a) m n
    | otherwise = do
        unsafeCopyToMatrix a b
        unsafeIOToST $
            unsafeWithVector s $ \ps ->
            unsafeWithMatrix b $ \pb ldb ->
                BLAS.axpy mn e ps 1 pb (ldb+1)
  where
    (m,n) = dimMatrix b
    mn = min m n

-- | Add two matrices.
addToMatrix :: (RMatrix m1, RMatrix m2, VNum e)
            => m1 e -> m2 e -> STMatrix s e -> ST s ()
addToMatrix = checkMatrixOp3 "addToMatrix" $ matrixVectorOp3 addToVector

-- | Add two matrices with the given scales.
addToMatrixWithScale :: (RMatrix m1, RMatrix m2, VNum e)
                     => e -> m1 e -> e -> m2 e -> STMatrix s e -> ST s ()
addToMatrixWithScale alpha a beta b c =
    (checkMatrixOp3 "addToMatrixWithScale" $
        matrixVectorOp3 (\x y z -> addToVectorWithScale alpha x beta y z))
        a b c

-- | Subtract two matrices.
subToMatrix :: (RMatrix m1, RMatrix m2, VNum e)
            => m1 e -> m2 e -> STMatrix s e -> ST s ()
subToMatrix = checkMatrixOp3 "subToMatrix" $ matrixVectorOp3 subToVector

-- | Negate the entries of a matrix.
negateToMatrix :: (RMatrix m, VNum e)
               => m e -> STMatrix s e -> ST s ()
negateToMatrix = checkMatrixOp2 "negateToMatrix" $
    matrixVectorOp2 negateToVector

-- | Scale the entries of a matrix by the given value.
scaleToMatrix :: (RMatrix m, VNum e)
              => e -> m e -> STMatrix s e -> ST s ()
scaleToMatrix e = checkMatrixOp2 "scaleToMatrix" $
    matrixVectorOp2 (scaleToVector e)

-- | Scale the rows of a matrix; @scaleRowsToMatrix s a c@ sets
-- @c := diag(s) * a@.
scaleRowsToMatrix :: (RVector v, RMatrix m, VNum e)
                  => v e -> m e -> STMatrix s e -> ST s ()
scaleRowsToMatrix s a b
    | dimVector s /= m || dimMatrix a /= (m,n) = error $
        printf ("scaleRowsToMatrix <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (dimVector s)
                (fst $ dimMatrix a) (snd $ dimMatrix a) m n
    | otherwise = do
        zipWithM_ (mulToVector s) (colsMatrix a) (colsMatrix b)
  where
    (m,n) = dimMatrix b

-- | Scale the columns of a matrix; @scaleColsToMatrix s a c@ sets
-- @c := a * diag(s)@.
scaleColsToMatrix :: (RVector v, RMatrix m, VNum e)
                  => v e -> m e -> STMatrix s e -> ST s ()
scaleColsToMatrix s a b 
    | dimVector s /= n || dimMatrix a /= (m,n) = error $
        printf ("scaleColsToMatrix <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (dimVector s)
                (fst $ dimMatrix a) (snd $ dimMatrix a) m n
    | otherwise = do
        es <- getElemsVector s
        sequence_ [ scaleToVector e x y
                  | (e,x,y) <- zip3 es (colsMatrix a) (colsMatrix b)
                  ]
  where
    (m,n) = dimMatrix b

-- | @rank1UpdateToMatrix alpha x y a c@ sets @c := alpha * x * y^H + a@.
rank1UpdateToMatrix :: (RVector v1, RVector v2, RMatrix m, BLAS2 e)
                    => e -> v1 e -> v2 e -> m e -> STMatrix s e -> ST s ()
rank1UpdateToMatrix alpha x y a c
    | dimVector x /= m || dimVector y /= n || dimMatrix a /= (m,n) = error $
        printf ("rank1UpdateToMatrix _ <vector with dim %d>"
                ++ " <vector with dim %d> <matrix with dim (%d,%d)>"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
                (dimVector x) (dimVector y) (fst $ dimMatrix a)
                (snd $ dimMatrix a) m n
    | otherwise = do
        unsafeCopyToMatrix a c
        unsafeIOToST $
            unsafeWithVector x $ \px ->
            unsafeWithVector y $ \py ->
            unsafeWithMatrix c $ \pc ldc ->
                BLAS.gerc m n alpha px 1 py 1 pc ldc
  where
    (m,n) = dimMatrix c

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

checkMatrixOp3 :: (RMatrix x, RMatrix y, RMatrix z)
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
                => (VectorView m e -> STVector s f -> ST s ())
                -> m e -> STMatrix s f -> ST s ()
matrixVectorOp2 f a b =
    case (maybeVectorViewMatrix a, maybeVectorViewMatrix b) of
        (Just x, Just y) -> f x y
        _ -> zipWithM_ f (colsMatrix a) (colsMatrix b)
{-# INLINE matrixVectorOp2 #-}

matrixVectorOp3 :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
                => (VectorView m1 e1 -> VectorView m2 e2 -> STVector s f -> ST s ())
                -> m1 e1 -> m2 e2 -> STMatrix s f -> ST s ()
matrixVectorOp3 f a b c =
    case (maybeVectorViewMatrix a, maybeVectorViewMatrix b,
            maybeVectorViewMatrix c) of
        (Just x, Just y, Just z) -> f x y z
        _ -> sequence_
            [ f x y z
            | (x,y,z) <- zip3 (colsMatrix a) (colsMatrix b) (colsMatrix c) ]
{-# INLINE matrixVectorOp3 #-}

newResultMatrix :: (RMatrix m, Storable f)
                => (m e -> STMatrix s f -> ST s a)
                -> m e
                -> ST s (STMatrix s f)
newResultMatrix f a = do
    c <- newMatrix_ (dimMatrix a)
    _ <- f a c
    return c
{-# INLINE newResultMatrix #-}

newResultMatrix2 :: (RMatrix m1, RMatrix m2, Storable g)
                 => (m1 e -> m2 f -> STMatrix s g -> ST s a)
                 -> m1 e
                 -> m2 f
                 -> ST s (STMatrix s g)
newResultMatrix2 f a1 a2 = do
    c <- newMatrix_ (dimMatrix a1)
    _ <- f a1 a2 c
    return c
{-# INLINE newResultMatrix2 #-}
