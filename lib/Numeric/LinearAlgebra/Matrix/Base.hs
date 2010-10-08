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

module Numeric.LinearAlgebra.Matrix.Base (
    Matrix(..),
    dim,
    
    fromList,
    fromRow,
    fromCol,
    zero,
    constant,
    
    at,
    unsafeAt,

    indices,
    elems,
    assocs,

    col,
    unsafeCol,
    cols,

    replace,
    unsafeReplace,
    accum,
    unsafeAccum,

    map,
    zipWith,
    unsafeZipWith,

    
    slice,
    unsafeSlice,
    
    splitRowsAt,
    dropRows,
    takeRows,

    splitColsAt,
    dropCols,
    takeCols,

    fromVector,
    toVector,

    isContig,
    unsafeFromForeignPtr,
    unsafeToForeignPtr,
    unsafeWith,

    ) where

import Prelude hiding ( read, map, zipWith )
import qualified Prelude as P

import Data.AEq( AEq(..) )
import Data.Typeable( Typeable )
import Foreign( ForeignPtr, Ptr, Storable )
import Text.Printf( printf )

import Numeric.LinearAlgebra.Vector( Vector )
import qualified Numeric.LinearAlgebra.Vector as V

-- | Immutable dense matrices.
data Matrix e =
      Matrix {-# UNPACK #-} !(Vector e) -- matrix data
             {-# UNPACK #-} !Int        -- row dimension
             {-# UNPACK #-} !Int        -- column dimension
             {-# UNPACK #-} !Int        -- leading dimension
  deriving (Typeable)

-- | Get the matrix dimensions (number of rows and number of columns).
dim :: (Storable e) => Matrix e -> (Int,Int)
dim (Matrix _ m n _) = (m,n)
{-# INLINE dim #-}

-- | Indicates if the elements of the matrix are stored contigously
isContig :: (Storable e) => Matrix e -> Bool
isContig (Matrix _ m _ lda) = lda == m || m == 0
{-# INLINE isContig #-}

-- | Create a matrix of the given dimension with elements initialized
-- to the values from the list, in column major order.
fromList :: (Storable e) => (Int,Int) -> [e] -> Matrix e
fromList (m,n) es
    | m < 0 || n < 0 = error $
        printf "fromList (%d,%d): negative dimension" m n
    | otherwise = let
        v = V.fromList (m*n) es
        lda = max 1 m
        in Matrix v m n lda
{-# INLINE fromList #-}

-- | Create a matrix of the given dimension with all elements initialized
-- to the given value
constant :: (Storable e) => (Int,Int) -> e -> Matrix e
constant (m,n) e
    | m < 0 || n < 0 = error $
        printf "constant (%d,%d): negative dimension" m n
    | otherwise = let
        v = V.constant (m*n) e
        lda = max 1 m
        in Matrix v m n lda
{-# INLINE constant #-}

-- | Create a zero of the given dimension with all elements initialized
-- to zero.
zero :: (Storable e, Num e) => (Int,Int) -> Matrix e
zero (m,n)
    | m < 0 || n < 0 = error $
        printf "zero (%d,%d): negative dimension" m n
    | otherwise = let
        v = V.zero (m*n)
        lda = max 1 m
        in Matrix v m n lda
{-# INLINE zero #-}

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
unsafeAt (Matrix v _ _ lda) (i,j) = 
    V.unsafeAt v (i + j * lda)
{-# INLINE unsafeAt #-}

-- | Get the indices of the elements in the matrix, in column-major order.
indices :: (Storable e) => Matrix e -> [(Int,Int)]
indices a = [ (i,j) | j <- [ 0..n-1 ], i <- [ 0..m-1 ] ]
  where (m,n) = dim a

-- | Returns a list of the elements of a matrix, in the same order as their
-- indices.
elems :: (Storable e) => Matrix e -> [e]
elems (Matrix v m _ lda)
    | lda == m  = V.elems v
    | otherwise = let
        breakCols [] = []
        breakCols es = let (c,es') = splitAt lda es in c:(breakCols es')

        dropJunk c = take m c

        in concatMap dropJunk $ breakCols (V.elems v)
{-# INLINE elems #-}

-- | Returns the contents of a matrix as a list of associations.
assocs :: (Storable e) => Matrix e -> [((Int,Int),e)]
assocs x = zip (indices x) (elems x)
{-# INLINE assocs #-}

-- | Version of 'replace' that doesn't range-check indices.
unsafeReplace :: (Storable e) => Matrix e -> [((Int,Int),e)] -> Matrix e
unsafeReplace (Matrix v m n lda) ijes = let
    ies = [ (i + j * lda, e) | ((i,j),e) <- ijes ]
    v' = V.unsafeReplace v ies
    lda' = max 1 m
    in Matrix v' m n lda'

-- | Create a new matrix by replacing the values at the specified indices.
replace :: (Storable e) => Matrix e -> [((Int,Int),e)] -> Matrix e
replace (Matrix v m n lda) ijes = let
    ies = [ if i < 0 || i >= m || j < 0 || j >= n
                then error $ printf
                         ("replace"
                         ++ " <matrix with dim (%d,%d)>"
                         ++ " [ ..((%d,%d),_).. ]"
                         ++ ": invalid index")
                         m n i j
                else (i + j * lda, e)
          | ((i,j),e) <- ijes
          ]
    v' = V.unsafeReplace v ies
    lda' = max 1 m
    in Matrix v' m n lda'

-- | Same as 'accum' but does not range-check indices.
unsafeAccum :: (Storable e)
            => (e -> e' -> e)
            -> Matrix e
            -> [((Int,Int), e')]
            -> Matrix e
unsafeAccum f (Matrix v m n lda) ijes = let
    ies = [ (i + j * lda, e) | ((i,j),e) <- ijes ]
    v' = V.unsafeAccum f v ies
    lda' = max 1 m
    in Matrix v' m n lda'

-- | @accum f@ takes a matrix and an association list and accumulates
-- pairs from the list into the matrix with the accumulating function @f@.
accum :: (Storable e)
      => (e -> e' -> e) 
      -> Matrix e
      -> [((Int,Int), e')]
      -> Matrix e
accum f (Matrix v m n lda) ijes = let
    ies = [ if i < 0 || i >= m || j < 0 || j >= n
                then error $ printf
                         ("accum"
                         ++ " <matrix with dim (%d,%d)>"
                         ++ " [ ..((%d,%d),_).. ]"
                         ++ ": invalid index")
                         m n i j
                else (i + j * lda, e)
          | ((i,j),e) <- ijes
          ]
    v' = V.unsafeAccum f v ies
    lda' = max 1 m
    in Matrix v' m n lda'

-- | Construct a new matrix by applying a function to every element of
-- a matrix.
map :: (Storable e, Storable e')
    => (e -> e')
    -> Matrix e
    -> Matrix e'
map f a = fromList (dim a) $ P.map f (elems a)
{-# INLINE map #-}

-- | Construct a new matrix by applying a function to every pair of elements
-- of two matrices.  The two matrices must have identical dimensions.
zipWith :: (Storable e, Storable e', Storable f)
        => (e -> e' -> f)
        -> Matrix e
        -> Matrix e'
        -> Matrix f
zipWith f a a'
    | mn /= mn' = error $
        printf ("zipWith"
                ++ " <matrix with dim %s> "
                ++ " <matrix with dim %s>"
                ++ ": dimension mismatch"
                ) (show mn) (show mn')
    | otherwise =
        unsafeZipWith f a a'
  where
    mn  = dim a
    mn' = dim a'    
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

-- | Version of 'col' that doesn't range-check indices.
unsafeCol :: (Storable e) => Matrix e -> Int -> Vector e
unsafeCol (Matrix v m _ lda) j = 
    V.unsafeSlice (j*lda) m v
{-# INLINE unsafeCol #-}

-- | Get a list of the columns of the matrix.
cols :: (Storable e) => Matrix e -> [Vector e]
cols a = P.map (unsafeCol a) [ 0..n-1 ]
  where
    (_,n) = dim a
{-# INLINE cols #-}

-- | @slice (i,j) (m,n) a@ creates a submatrix view of @a@ starting at
-- element @(i,j)@ and having dimensions @(m,n)@.
slice :: (Storable e)
      => (Int,Int)
      -> (Int,Int)
      -> Matrix e
      -> Matrix e
slice (i,j) (m',n') a
    | (i < 0 || m' < 0 || i + m' > m 
       || j < 0 || n' < 0 || j + n' > n) = error $
        printf ( "slice"
               ++ " (%d,%d)"
               ++ " (%d,%d)"
               ++ " <matrix with dim (%d,%d)>"
               ++ ": index out of range"
               ) i j m' n' m n
    | otherwise =
        unsafeSlice (i,j) (m',n') a
  where
    (m,n) = dim a
{-# INLINE slice #-}

-- | Version of 'slice' that doesn't range-check indices.
unsafeSlice :: (Storable e)
            => (Int,Int) -> (Int,Int) -> Matrix e -> Matrix e
unsafeSlice (i,j) (m',n') (Matrix v _ _ lda) = let
    o = i + j*lda
    l = if m' == 0 || n' == 0
            then 0
            else lda * (n' - 1) + m'
    v' = V.unsafeSlice o l v
    in Matrix v' m' n' lda
{-# INLINE unsafeSlice #-}

-- | Create a view of a matrix by taking the initial rows.
takeRows :: (Storable e) => Int -> Matrix e -> Matrix e
takeRows i a = slice (0,0) (i,n) a
  where
    (_,n) = dim a

-- | Create a view of a matrix by dropping the initial rows.
dropRows :: (Storable e) => Int -> Matrix e -> Matrix e
dropRows i a = slice (i,0) (m-i,n) a
  where
    (m,n) = dim a

-- | Split a matrix into two blocks and returns views into the blocks.  If
-- @(a1, a2) = splitRowsAt i a@, then
-- @a1 = slice (0,0) (i,n) a@ and
-- @a2 = slice (i,0) (m-i,n) a@, where @(m,n)@ is the dimension of @a@.
splitRowsAt :: (Storable e) => Int -> Matrix e -> (Matrix e, Matrix e)
splitRowsAt i a
    | i < 0 || i > m = error $
        printf ("splitRowsAt %d <matrix with dim (%d,%d)>:"
                ++ " invalid index") i m n
    | otherwise = let
        a1 = unsafeSlice (0,0) (i,n)   a
        a2 = unsafeSlice (i,0) (m-i,n) a
    in (a1,a2)
  where
    (m,n) = dim a
{-# INLINE splitRowsAt #-}

-- | Create a view of a matrix by taking the initial columns.
takeCols :: (Storable e) => Int -> Matrix e -> Matrix e
takeCols j a = slice (0,0) (m,j) a
  where
    (m,_) = dim a

-- | Create a view of a matrix by dropping the initial columns.
dropCols :: (Storable e) => Int -> Matrix e -> Matrix e
dropCols j a = slice (0,j) (m,n-j) a
  where
    (m,n) = dim a

-- | Split a matrix into two blocks and returns views into the blocks.  If
-- @(a1, a2) = splitColsAt j a@, then
-- @a1 = slice (0,0) (m,j) a@ and
-- @a2 = slice (0,j) (m,n-j) a@, where @(m,n)@ is the dimension of @a@.
splitColsAt :: (Storable e) => Int -> Matrix e -> (Matrix e, Matrix e)
splitColsAt j a
    | j < 0 || j > n = error $
        printf ("splitColsAt %d <matrix with dim (%d,%d)>:"
                ++ " invalid index") j m n
    | otherwise = let
        a1 = unsafeSlice (0,0) (m,j)   a
        a2 = unsafeSlice (0,j) (m,n-j) a
    in (a1,a2)
  where
    (m,n) = dim a
{-# INLINE splitColsAt #-}



      

-- | Convert a matrix to a vector by stacking its columns.
toVector :: (Storable e)
         => Matrix e
         -> Vector e
toVector a@(Matrix v m n _)
    | isContig a = v
    | otherwise  = V.fromList (m*n) $ elems a

-- | Cast a vector to a matrix of the given shape.
fromVector :: (Storable e)
           => (Int,Int)
           -> Vector e
           -> Matrix e
fromVector (m,n) v
    | nv /= m * n = error $
        printf ("fromVector"
               ++ " (%d,%d)"
               ++ " <vector with dim %d>"
               ++ ": dimension mismatch"
               ) m n nv
    | otherwise =
        Matrix v m n (max 1 m)
  where
    nv = V.dim v
{-# INLINE fromVector #-}

-- | Cast a vector to a matrix with one column.
fromCol :: (Storable e)
        => Vector e
        -> Matrix e
fromCol v = Matrix v m 1 (max 1 m)
  where
    m = V.dim v

-- | Cast a vector to a matrix with one row.
fromRow :: (Storable e)
        => Vector e
        -> Matrix e
fromRow v = Matrix v 1 n 1
  where
    n = V.dim v

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

compareWith :: (Storable e, Storable e')
            => (e -> e' -> Bool)
            -> Matrix e
            -> Matrix e'
            -> Bool
compareWith cmp a a' =
    dim a == dim a'
    && and (P.zipWith cmp (elems a) (elems a'))
{-# INLINE compareWith #-}

-- | Create a matrix from a 'ForeignPtr' with offset, dimensions, and lda. The
-- data may not be modified through the ForeignPtr afterwards.
unsafeFromForeignPtr :: (Storable e)
                     => ForeignPtr e -- ^ pointer
                     -> Int	         -- ^ offset
                     -> (Int,Int)    -- ^ dimensions
                     -> Int          -- ^ leading dimension (lda)
                     -> Matrix e
unsafeFromForeignPtr p o (m,n) lda = let
    v = V.unsafeFromForeignPtr p o (m*lda)
    in Matrix v m n lda
{-# INLINE unsafeFromForeignPtr #-}

-- | Yield the underlying 'ForeignPtr' together with the offset to the data
-- the matrix dimensions, and the lda. The data may not be modified through
-- the 'ForeignPtr'.
unsafeToForeignPtr :: (Storable e)
                   => Matrix e
                   -> (ForeignPtr e, Int, (Int,Int), Int)
unsafeToForeignPtr (Matrix v m n lda) = let
    (f,o,_) = V.unsafeToForeignPtr v
    in (f, o, (m,n), lda)
{-# INLINE unsafeToForeignPtr #-}

-- | Execute an 'IO' action with a pointer to the first element in the
-- matrix and the leading dimension (lda).
unsafeWith :: (Storable e) => Matrix e -> (Ptr e -> Int -> IO a) -> IO a
unsafeWith (Matrix v _ _ lda) f =
    V.unsafeWith v $ \p -> f p lda
{-# INLINE unsafeWith #-}
