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
      
import Control.Monad( forM_ )
import Control.Monad.ST( ST, RealWorld, unsafeIOToST )
import Data.Maybe( fromMaybe )
import Data.Typeable( Typeable )
import Foreign( ForeignPtr, Ptr, advancePtr, peekElemOff, pokeElemOff,
    mallocForeignPtrArray )
import Text.Printf( printf )
import Unsafe.Coerce( unsafeCoerce )

import Numeric.LinearAlgebra.Types
import qualified Foreign.BLAS as BLAS
import Numeric.LinearAlgebra.Vector( STVector, RVector )
import qualified Numeric.LinearAlgebra.Vector as V

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
    dim :: (Storable e) => m e -> (Int,Int)
    
    unsafeSlice :: (Storable e)
                => (Int,Int) -> (Int,Int) -> m e -> m e

    unsafeWithColView :: (Storable e)
                      => m e
                      -> Int 
                      -> (forall v . RVector v => v e -> ST s a)
                      -> ST s a

    -- | Perform an action with a list of views of the matrix columns.
    withColViews :: (Storable e)
                 => m e
                 -> (forall v . RVector v => [v e] -> ST s a)
                 -> ST s a

    -- | Possibly view a matrix as a vector and perform an action on the
    -- view.  This only succeeds if the matrix is stored contiguously in
    -- memory, i.e. if the matrix contains a single column or the \"lda\"
    -- of the matrix is equal to the number of rows.
    maybeWithVectorView :: (Storable e)
                        => m e
                        -> (forall v . RVector v => v e -> ST s a)
                        -> Maybe (ST s a)

    -- | Execute an 'IO' action with a pointer to the first element in the
    -- matrix and the leading dimension (lda).
    unsafeWith :: (Storable e) => m e -> (Ptr e -> Int -> IO a) -> IO a

    -- | Convert a vector to a @ForeignPtr@, and offset, dimensions,
    -- and leading dimension (lda).
    unsafeToForeignPtr :: (Storable e)
                       => m e -> (ForeignPtr e, Int, (Int,Int), Int)
                             
    -- | Cast a @ForeignPtr@ to a matrix.
    unsafeFromForeignPtr :: (Storable e)
                         => ForeignPtr e -- ^ the pointer
                         -> Int          -- ^ the offset
                         -> (Int,Int)    -- ^ the dimensions
                         -> Int          -- ^ the leading dimension
                         -> m e


instance RMatrix (STMatrix s) where
    dim (STMatrix _ m n _) = (m,n)
    {-# INLINE dim #-}

    unsafeWithColView a j f = cast $ unsafeWithSTColView a j (cast . f)
      where
        cast :: ST s' a -> ST s a
        cast = unsafeCoerce
    {-# INLINE unsafeWithColView #-}

    withColViews a f = cast $ withSTColViews a (cast . f)
      where
        cast :: ST s' a -> ST s a
        cast = unsafeCoerce
    {-# INLINE withColViews #-}

    unsafeSlice (i,j) (m',n') (STMatrix a _ _ lda) = let
        o = i + j*lda
        l = if m' == 0 || n' == 0
                then 0
                else lda * (n' - 1) + m'
        a' = V.unsafeSlice o l a
        in STMatrix a' m' n' lda
    {-# INLINE unsafeSlice #-}
    
    maybeWithVectorView a f =
        cast `fmap` maybeWithSTVectorView a (cast . f)
      where
        cast :: ST s' a -> ST s a
        cast = unsafeCoerce
    {-# INLINE maybeWithVectorView #-}

    unsafeWith (STMatrix a _ _ lda) f =
        V.unsafeWith a $ \p -> f p lda
    {-# INLINE unsafeWith #-}

    unsafeToForeignPtr (STMatrix a m n lda) = let
        (f,o,_) = V.unsafeToForeignPtr a
        in (f, o, (m,n), lda)
    {-# INLINE unsafeToForeignPtr #-}

    unsafeFromForeignPtr f o (m,n) lda = let
        d = if m == 0 then 0 else lda * n
        a = V.unsafeFromForeignPtr f o d
        in (STMatrix a m n lda)
    {-# INLINE unsafeFromForeignPtr #-}


unsafeSTColView :: (Storable e)
                => STMatrix s e
                -> Int
                -> STVector s e
unsafeSTColView (STMatrix a m _ lda) j = let
    o = j * lda
    in V.unsafeSlice o m a
{-# INLINE unsafeSTColView #-}

-- | Perform an action with a view of a matrix column (no index checking).
unsafeWithSTColView :: (Storable e)
                    => STMatrix s e
                    -> Int
                    -> (STVector s e -> ST s a)
                    -> ST s a
unsafeWithSTColView a j f = f $ unsafeSTColView a j
{-# INLINE unsafeWithSTColView #-}

-- | Get a list of views of the matrix columns.
stCols :: (Storable e)
       => STMatrix s e
       -> [STVector s e]
stCols (STMatrix a m n lda) = let
    os = take n $ [ 0, lda .. ]
    xs = [ V.unsafeSlice o m a | o <- os ]
    in xs
{-# INLINE stCols #-}

-- | Perform an action with a list of views of the matrix columns.  See
-- also 'withColViews'.
withSTColViews :: (Storable e)
               => STMatrix s e
               -> ([STVector s e] -> ST s a)
               -> ST s a
withSTColViews a f = f $ stCols a
{-# INLINE withSTColViews #-}

-- | Possibly view a matrix as a vector.  This only succeeds if the
-- matrix is stored contiguously in memory, i.e. if the matrix contains
-- a single column or the \"lda\" of the matrix is equal to the number
-- of rows.
maybeSTVectorView :: (Storable e)
                  => STMatrix s e
                  -> Maybe (STVector s e)
maybeSTVectorView (STMatrix a m n lda) | lda == max 1 m = Just a
                                       | n   == 1       = Just a
                                       | otherwise      = Nothing
{-# INLINE maybeSTVectorView #-}

-- | Possibly view a matrix as a vector and perform an action on the
-- view.  This succeeds when the matrix is stored contiguously in memory,
-- i.e. if the matrix contains a single column or the \"lda\" of the matrix
-- is equal to the number of rows.  See also 'maybeWithVectorView'.
maybeWithSTVectorView :: (Storable e)
                      => STMatrix s e
                      -> (STVector s e -> ST s a)
                      -> Maybe (ST s a)
maybeWithSTVectorView a f = f `fmap` maybeSTVectorView a
{-# INLINE maybeWithSTVectorView #-}

-- | Cast a vector to a matrix of the given shape.  See also
-- 'withViewFromVector'.
fromSTVector :: (Storable e)
                 => (Int,Int)
                 -> STVector s e
                 -> STMatrix s e
fromSTVector mn@(m,n) v
    | V.dim v /= m*n = error $
        printf ("fromSTVector (%d,%d) <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim v)
    | otherwise =
        cast v
  where
    cast x = let
        (fptr,o,_) = V.unsafeToForeignPtr x
        lda = max 1 m
        in unsafeFromForeignPtr fptr o mn lda
{-# INLINE fromSTVector #-}

-- | Cast a vector to a matrix with one column.  See also
-- 'withViewFromCol'.
fromSTCol :: (Storable e)
              => STVector s e
              -> STMatrix s e
fromSTCol v = fromSTVector (V.dim v, 1) v
{-# INLINE fromSTCol #-}

-- | Cast a vector to a matrix with one row.  See also
-- 'withViewFromRow'.
fromSTRow :: (Storable e)
              => STVector s e
              -> STMatrix s e
fromSTRow v = fromSTVector (V.dim v, 1) v
{-# INLINE fromSTRow #-}

-- | View a vector as a matrix of the given shape and pass it to
-- the specified function.
withViewFromVector :: (RVector v, Storable e)
                   => (Int,Int)
                   -> v e
                   -> (forall m . RMatrix m => m e -> a)
                   -> a
withViewFromVector mn@(m,n) v f
    | V.dim v /= m*n = error $
        printf ("withViewFromVector (%d,%d) <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim v)
    | otherwise =
        f (cast v)
  where
    cast :: (RVector v, Storable e) => v e -> STMatrix s e
    cast x = let
        (fptr,o,_) = V.unsafeToForeignPtr x
        lda = max 1 m
        in unsafeFromForeignPtr fptr o mn lda
{-# INLINE withViewFromVector #-}

-- | View a mutable vector as a mutable matrix of the given shape and pass it
-- to the specified function.
withViewFromSTVector :: (Storable e)
                     => (Int,Int)
                     -> STVector s e
                     -> (STMatrix s e -> ST s a)
                     -> ST s a
withViewFromSTVector mn@(m,n) v f
    | V.dim v /= m*n = error $
        printf ("withViewFromSTVector (%d,%d) <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim v)
    | otherwise =
        f (cast v)
  where
    cast :: (Storable e) => STVector s e -> STMatrix s e
    cast x = let
        (fptr,o,_) = V.unsafeToForeignPtr x
        lda = max 1 m
        in unsafeFromForeignPtr fptr o mn lda
{-# INLINE withViewFromSTVector #-}

-- | View a vector as a matrix with one column and pass it to
-- the specified function.
withViewFromCol :: (RVector v, Storable e)
                 => v e
                 -> (forall m . RMatrix m => m e -> a)
                 -> a
withViewFromCol v = withViewFromVector (V.dim v, 1) v
{-# INLINE withViewFromCol #-}

-- | View a mutable vector as a mutable matrix with one column and pass it to
-- the specified function.
withViewFromSTCol :: (Storable e)
                 => STVector s e
                 -> (STMatrix s e -> ST s a)
                 -> ST s a
withViewFromSTCol v = withViewFromSTVector (V.dim v, 1) v
{-# INLINE withViewFromSTCol #-}

-- | View a vector as a matrix with one row and pass it to
-- the specified function.
withViewFromRow :: (RVector v, Storable e)
                 => v e
                 -> (forall m . RMatrix m => m e -> a)
                 -> a
withViewFromRow v = withViewFromVector (1, V.dim v) v
{-# INLINE withViewFromRow #-}

-- | View a mutable vector as a mutable matrix with one row and pass it to
-- the specified function.
withViewFromSTRow :: (Storable e)
                  => STVector s e
                  -> (STMatrix s e -> ST s a)
                  -> ST s a
withViewFromSTRow v = withViewFromSTVector (1, V.dim v) v
{-# INLINE withViewFromSTRow #-}

-- | Perform an action with a view of a matrix column.
withColView :: (RMatrix m, Storable e)
            => m e
            -> Int
            -> (forall v . RVector v => v e -> ST s a)
            -> ST s a
withColView a j f
    | j < 0 || j >= n = error $
        printf ("withColView <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n j
    | otherwise =
        unsafeWithColView a j f
  where
    (m,n) = dim a
{-# INLINE withColView #-}

-- | Perform an action with a view of a mutable matrix column.
withSTColView :: (Storable e)
              => STMatrix s e
              -> Int
              -> (STVector s e -> ST s a)
              -> ST s a
withSTColView a j f
    | j < 0 || j >= n = error $
        printf ("withSTColView <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n j
    | otherwise =
        unsafeWithSTColView a j f
  where
    (m,n) = dim a
{-# INLINE withSTColView #-}

-- | @slice (i,j) (m,n) a@ creates a submatrix view of @a@ starting at
-- element @(i,j)@ and having dimensions @(m,n)@.
slice :: (RMatrix m, Storable e)
      => (Int,Int)
      -> (Int,Int)
      -> m e
      -> m e
slice (i,j) (m',n') a
    | (i < 0 || m' < 0 || i + m' > m 
       || j < 0 || n' < 0 || j + n' > n) = error $
        printf ("slice (%d,%d) (%d,%d) <matrix with dim (%d,%d)>:"
                ++ " index out of range") i j m' n' m n
    | otherwise =
        unsafeSlice (i,j) (m',n') a
  where
    (m,n) = dim a
{-# INLINE slice #-}

-- | Create a new matrix of given shape, but do not initialize the elements.
new_ :: (Storable e) => (Int,Int) -> ST s (STMatrix s e)
new_ (m,n) 
    | m < 0 || n < 0 = error $
        printf "new_ (%d,%d): invalid dimensions" m n
    | otherwise = unsafeIOToST $ do
        f <- mallocForeignPtrArray (m*n)
        return $ unsafeFromForeignPtr f 0 (m,n) (max 1 m)

-- | Create a matrix with every element initialized to the same value.
new :: (Storable e) => (Int,Int) -> e -> ST s (STMatrix s e)
new (m,n) e = do
    a <- new_ (m,n)
    setElems a $ replicate (m*n) e
    return a

-- | Creates a new matrix by copying another one.    
newCopy :: (RMatrix m, Storable e) => m e -> ST s (STMatrix s e)
newCopy a = do
    b <- new_ (dim a)
    unsafeCopyTo b a
    return b

-- | @copyTo dst src@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyTo :: (RMatrix m, Storable e) => STMatrix s e -> m e -> ST s ()
copyTo = checkOp2 "copyTo" unsafeCopyTo
{-# INLINE copyTo #-}

unsafeCopyTo :: (RMatrix m, Storable e) => STMatrix s e -> m e -> ST s ()
unsafeCopyTo = vectorOp2 V.unsafeCopyTo
{-# INLINE unsafeCopyTo #-}

-- | Create a view of a matrix by taking the initial rows.
takeRows :: (RMatrix m, Storable e) => Int -> m e -> m e
takeRows i a = slice (0,0) (i,n) a
  where
    (_,n) = dim a

-- | Create a view of a matrix by dropping the initial rows.
dropRows :: (RMatrix m, Storable e) => Int -> m e -> m e
dropRows i a = slice (i,0) (m-i,n) a
  where
    (m,n) = dim a

-- | Split a matrix into two blocks and returns views into the blocks.  If
-- @(a1, a2) = splitRowsAt i a@, then
-- @a1 = slice (0,0) (i,n) a@ and
-- @a2 = slice (i,0) (m-i,n) a@, where @(m,n)@ is the dimension of @a@.
splitRowsAt :: (RMatrix m, Storable e) => Int -> m e -> (m e, m e)
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
takeCols :: (RMatrix m, Storable e) => Int -> m e -> m e
takeCols j a = slice (0,0) (m,j) a
  where
    (m,_) = dim a

-- | Create a view of a matrix by dropping the initial columns.
dropCols :: (RMatrix m, Storable e) => Int -> m e -> m e
dropCols j a = slice (0,j) (m,n-j) a
  where
    (m,n) = dim a

-- | Split a matrix into two blocks and returns views into the blocks.  If
-- @(a1, a2) = splitColsAt j a@, then
-- @a1 = slice (0,0) (m,j) a@ and
-- @a2 = slice (0,j) (m,n-j) a@, where @(m,n)@ is the dimension of @a@.
splitColsAt :: (RMatrix m, Storable e) => Int -> m e -> (m e, m e)
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

-- | Get the indices of the elements in the matrix, in column-major order.
indices :: (RMatrix m, Storable e) => m e -> [(Int,Int)]
indices a = [ (i,j) | j <- [ 0..n-1 ], i <- [ 0..m-1 ] ]
  where (m,n) = dim a
  
-- | Lazily get the elements of the matrix, in column-major order.  
getElems :: (RMatrix m, Storable e) => m e -> ST s [e]
getElems a = case maybeWithVectorView a V.getElems of
    Just es -> es
    Nothing -> withColViews a $ \xs ->
                   concat `fmap` mapM V.getElems xs

-- | Get the elements of the matrix, in column-major order.
getElems' :: (RMatrix m, Storable e) => m e -> ST s [e]
getElems' a = case maybeWithVectorView a V.getElems' of
    Just es -> es
    Nothing -> withColViews a $ \xs ->
                   concat `fmap` mapM V.getElems' xs

-- | Lazily get the association list of the matrix, in column-major order.
getAssocs :: (RMatrix m, Storable e) => m e -> ST s [((Int,Int),e)]
getAssocs a = do
    es <- getElems a
    return $ zip (indices a) es

-- | Get the association list of the matrix, in column-major order.
getAssocs' :: (RMatrix m, Storable e) => m e -> ST s [((Int,Int),e)]
getAssocs' a = do
    es <- getElems' a
    return $ zip (indices a) es

-- | Set all of the values of the matrix from the elements in the list,
-- in column-major order.
setElems :: (Storable e) => STMatrix s e -> [e] -> ST s ()
setElems a es =
    case maybeWithSTVectorView a (`V.setElems` es) of
        Just st  -> st
        Nothing -> go 0 es
  where
    (m,n) = dim a
    go j [] | j == n = return ()
    go j [] | j < n = error $ 
        printf ("setElems <matrix with dim (%d,%d>"
                ++ "<list with length %d>: not enough elements)") m n (j*m)
    go j es' =
        let (es1', es2') = splitAt m es'
        in do
            withSTColView a j (`V.setElems` es1')
            go (j+1) es2'

-- | Set the given values in the matrix.  If an index is repeated twice,
-- the value is implementation-defined.
setAssocs :: (Storable e) => STMatrix s e -> [((Int,Int),e)] -> ST s ()
setAssocs a ies =
    sequence_ [ write a i e | (i,e) <- ies ]

unsafeSetAssocs :: (Storable e) => STMatrix s e -> [((Int,Int),e)] -> ST s ()
unsafeSetAssocs a ies =
    sequence_ [ unsafeWrite a i e | (i,e) <- ies ]

-- | Set the specified row of the matrix to the given vector.
setRow :: (RVector v, Storable e)
       => STMatrix s e -> Int -> v e -> ST s ()
setRow a i x
    | i < 0 || i >= m = error $
        printf ("setRow <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n i
    | V.dim x /= n = error $
        printf ("setRow <matrix with dim (%d,%d)> _"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim x)
    | otherwise =
        unsafeSetRow a i x
  where
    (m,n) = dim a
{-# INLINE setRow #-}

unsafeSetRow :: (RVector v, Storable e)
             => STMatrix s e -> Int -> v e -> ST s ()
unsafeSetRow a i x = do
    jes <- V.getAssocs x
    sequence_ [ unsafeWrite a (i,j) e | (j,e) <- jes ]
{-# INLINE unsafeSetRow #-}

-- | Copy the specified row of the matrix to the vector.
getRow :: (RMatrix m, Storable e)
       => m e -> Int -> STVector s e -> ST s ()
getRow a i x
    | i < 0 || i >= m = error $
        printf ("getRow <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n i
    | V.dim x /= n = error $
        printf ("getRow <matrix with dim (%d,%d)> _"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim x)
    | otherwise =
        unsafeGetRow a i x
  where
    (m,n) = dim a
{-# INLINE getRow #-}

unsafeGetRow :: (RMatrix m, Storable e)
             => m e -> Int -> STVector s e -> ST s ()
unsafeGetRow a i x = 
    forM_ [ 0..n-1 ] $ \j -> do
        e <- unsafeRead a (i,j)
        V.unsafeWrite x j e
  where
    (_,n) = dim a
{-# INLINE unsafeGetRow #-}

-- | Set the diagonal of the matrix to the given vector.
setDiag :: (RVector v, Storable e)
        => STMatrix s e -> v e -> ST s ()
setDiag a x
    | V.dim x /= mn = error $
        printf ("setRow <matrix with dim (%d,%d)>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim x)
    | otherwise =
        unsafeSetDiag a x
  where
    (m,n) = dim a
    mn = min m n
{-# INLINE setDiag #-}

unsafeSetDiag :: (RVector v, Storable e)
              => STMatrix s e -> v e -> ST s ()
unsafeSetDiag a x = do
    ies <- V.getAssocs x
    sequence_ [ unsafeWrite a (i,i) e | (i,e) <- ies ]
{-# INLINE unsafeSetDiag #-}

-- | Copy the diagonal of the matrix to the vector.
getDiag :: (RMatrix m, Storable e)
        => m e -> STVector s e -> ST s ()
getDiag a x
    | V.dim x /= mn = error $
        printf ("getDiag <matrix with dim (%d,%d)>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim x)
    | otherwise =
        unsafeGetDiag a x
  where
    (m,n) = dim a
    mn = min m n
{-# INLINE getDiag #-}

unsafeGetDiag :: (RMatrix m, Storable e)
              => m e -> STVector s e -> ST s ()
unsafeGetDiag a x = 
    forM_ [ 0..mn-1 ] $ \i -> do
        e <- unsafeRead a (i,i)
        V.unsafeWrite x i e
  where
    (m,n) = dim a
    mn = min m n
{-# INLINE unsafeGetDiag #-}

-- | Get the element stored at the given index.
read :: (RMatrix m, Storable e) => m e -> (Int,Int) -> ST s e
read a (i,j)
    | i < 0 || i >= m || j < 0 || j >= n = error $
        printf ("read <matrix with dim (%d,%d)> (%d,%d):"
                ++ " index out of range") m n i j
    | otherwise =
        unsafeRead a (i,j)
  where
    (m,n) = dim a
{-# INLINE read #-}

unsafeRead :: (RMatrix m, Storable e) => m e -> (Int,Int) -> ST s e
unsafeRead a (i,j) = unsafeIOToST $
    unsafeWith a $ \p lda ->
        peekElemOff p (i + j * lda)
{-# INLINE unsafeRead #-}

-- | Set the element stored at the given index.
write :: (Storable e)
      => STMatrix s e -> (Int,Int) -> e -> ST s ()
write a (i,j)
    | i < 0 || i >= m || j < 0 || j >= n = error $
        printf ("write <matrix with dim (%d,%d)> (%d,%d):"
                ++ " index out of range") m n i j
    | otherwise =
        unsafeWrite a (i,j)
  where
    (m,n) = dim a
{-# INLINE write #-}

unsafeWrite :: (Storable e)
            => STMatrix s e -> (Int,Int) -> e -> ST s ()
unsafeWrite a (i,j) e = unsafeIOToST $
    unsafeWith a $ \p lda ->
        pokeElemOff p (i + j * lda) e
{-# INLINE unsafeWrite #-}

-- | Modify the element stored at the given index.
modify :: (Storable e)
       => STMatrix s e -> (Int,Int) -> (e -> e) -> ST s ()
modify a (i,j)
    | i < 0 || i >= m || j < 0 || j >= n = error $
        printf ("modify <matrix with dim (%d,%d)> (%d,%d):"
                ++ " index out of range") m n i j
    | otherwise =
        unsafeModify a (i,j)
  where
    (m,n) = dim a
{-# INLINE modify #-}

unsafeModify :: (Storable e)
             => STMatrix s e -> (Int,Int) -> (e -> e) -> ST s ()
unsafeModify a (i,j) f = unsafeIOToST $
    unsafeWith a $ \p lda -> 
        let o = i + j * lda
        in do
            e <- peekElemOff p o
            pokeElemOff p o $ f e
{-# INLINE unsafeModify #-}

-- | @mapTo dst f src@ replaces @dst@ elementwise with @f(src)@.
mapTo :: (RMatrix m, Storable e, Storable f)
      => STMatrix s f
      -> (e -> f)
      -> m e
      -> ST s ()
mapTo dst f src = (checkOp2 "mapTo _" $ \z x -> unsafeMapTo z f x) dst src
{-# INLINE mapTo #-}
                            
unsafeMapTo :: (RMatrix m, Storable e, Storable f)
            => STMatrix s f
            -> (e -> f)
            -> m e
            -> ST s ()
unsafeMapTo dst f src =
    fromMaybe colwise $ maybeWithSTVectorView dst $ \vdst ->
    fromMaybe colwise $ maybeWithVectorView src $ \vsrc ->
        V.unsafeMapTo vdst f vsrc
  where
    colwise = withSTColViews dst $ \zs ->
              withColViews   src $ \xs ->
                  sequence_ [ V.unsafeMapTo z f x
                            | (z,x) <- zip zs xs
                            ]

-- | @zipWithTo dst f x y@ replaces @dst@ elementwise with @f(x, y)@.
zipWithTo :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
          => STMatrix s f
          -> (e1 -> e2 -> f)
          -> m1 e1
          -> m2 e2
          -> ST s ()
zipWithTo dst f x y = 
    (checkOp3 "zipWithTo _" $ \dst1 x1 y1 -> unsafeZipWithTo dst1 f x1 y1)
        dst x y
{-# INLINE zipWithTo #-}

unsafeZipWithTo :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
                => STMatrix s f
                -> (e1 -> e2 -> f)
                -> m1 e1
                -> m2 e2
                -> ST s ()
unsafeZipWithTo dst f x y =
    fromMaybe colwise $ maybeWithSTVectorView dst $ \vdst ->
    fromMaybe colwise $ maybeWithVectorView x $ \vx ->
    fromMaybe colwise $ maybeWithVectorView y $ \vy ->
        V.unsafeZipWithTo vdst f vx vy
  where
    colwise = withSTColViews dst $ \vdsts ->
              withColViews   x $ \vxs ->
              withColViews   y $ \vys ->
              
                  sequence_ [ V.unsafeZipWithTo vdst f vx vy
                            | (vdst,vx,vy) <- zip3 vdsts vxs vys
                            ]

-- | Set every element in the matrix to a default value.  For
-- standard numeric types (including 'Double', 'Complex Double', and 'Int'),
-- the default value is '0'.
clear :: (Storable e) => STMatrix s e -> ST s ()
clear a = case maybeSTVectorView a of
    Just x  -> V.clear x
    Nothing -> withSTColViews a $ mapM_ V.clear

-- | Add a vector to the diagonal of a matrix.
shiftDiagByM_ :: (RVector v, BLAS1 e)
              => v e -> STMatrix s e -> ST s ()
shiftDiagByM_ s a
    | V.dim s /= mn = error $
        printf ("shiftDiagByM_"
                ++ " <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
                (V.dim s)
                m n
    | otherwise = shiftDiagByWithScaleM_ 1 s a
  where
    (m,n) = dim a
    mn = min m n

-- | Add a scaled vector to the diagonal of a matrix.
shiftDiagByWithScaleM_ :: (RVector v, BLAS1 e)
                       => e -> v e -> STMatrix s e -> ST s ()
shiftDiagByWithScaleM_ e s a
    | V.dim s /= mn = error $
        printf ("shiftDiagByWithScaleM_"
                ++ " _"
                ++ " <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
                (V.dim s)
                m n
    | otherwise = do
        unsafeIOToST $
            V.unsafeWith s $ \ps ->
            unsafeWith a $ \pa lda ->
                BLAS.axpy mn e ps 1 pa (lda+1)
  where
    (m,n) = dim a
    mn = min m n

-- | Add two matrices.
addTo :: (RMatrix m1, RMatrix m2, VNum e)
      =>  STMatrix s e -> m1 e -> m2 e -> ST s ()
addTo = checkOp3 "addTo" $ vectorOp3 V.addTo

-- | Subtract two matrices.
subTo :: (RMatrix m1, RMatrix m2, VNum e)
      => STMatrix s e -> m1 e -> m2 e -> ST s ()
subTo = checkOp3 "subTo" $ vectorOp3 V.subTo

-- | Conjugate the entries of a matrix.
conjugateTo :: (RMatrix m, VNum e)
            => STMatrix s e -> m e -> ST s ()
conjugateTo = checkOp2 "conjugateTo" $
    vectorOp2 V.conjugateTo

-- | Negate the entries of a matrix.
negateTo :: (RMatrix m, VNum e)
         => STMatrix s e -> m e -> ST s ()
negateTo = checkOp2 "negateTo" $
    vectorOp2 V.negateTo

-- | Scale the entries of a matrix by the given value.
scaleByM_ :: (BLAS1 e)
          => e -> STMatrix s e -> ST s ()
scaleByM_ e = vectorOp (V.scaleByM_ e)

-- | @addWithScaleM_ alpha x y@ sets @y := alpha * x + y@.
addWithScaleM_ :: (RMatrix m, BLAS1 e)
               => e -> m e -> STMatrix s e -> ST s ()
addWithScaleM_ e = checkOp2 "addWithScaleM_" $
    unsafeAddWithScaleM_ e

unsafeAddWithScaleM_ :: (RMatrix m, BLAS1 e)
                     => e -> m e -> STMatrix s e -> ST s ()
unsafeAddWithScaleM_ alpha x y =
    fromMaybe colwise $ maybeWithVectorView   x $ \vx ->
    fromMaybe colwise $ maybeWithSTVectorView y $ \vy ->
        V.unsafeAddWithScaleM_ alpha vx vy
  where
    colwise = withColViews   x $ \vxs ->
              withSTColViews y $ \vys ->
                  sequence_ [ V.unsafeAddWithScaleM_ alpha vx vy
                            | (vx,vy) <- zip vxs vys ]                

-- | Scale the rows of a matrix; @scaleRowsTo dst s a@ sets
-- @dst := diag(s) * a@.
scaleRowsTo :: (RVector v, RMatrix m, VNum e)
            => STMatrix s e -> v e -> m e -> ST s ()
scaleRowsTo dst s a
    | V.dim s /= m || dim a /= (m,n) = error $
        printf ("scaleRowsTo"
                ++ " <matrix with dim (%d,%d)>"
                ++ " <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
                m n
                (V.dim s)
                (fst $ dim a) (snd $ dim a)
    | otherwise =
        withSTColViews dst $ \ys ->
        withColViews   a $ \xs ->
            sequence_ [ V.mulTo y s x | (y,x) <- zip ys xs ]
  where
    (m,n) = dim dst

-- | Scale the columns of a matrix; @scaleCols_ a s@ sets
-- @a := a * diag(s)@.
scaleCols_ :: (RVector v, BLAS1 e)
            => STMatrix s e -> v e -> ST s ()
scaleCols_ a s
    | V.dim s /= n  = error $
        printf ("scaleCols_"
                ++ " <matrix with dim (%d,%d)>"        
                ++ " <vector with dim %d>"
                ++ ": dimension mismatch") 
                m n
                (V.dim s)
    | otherwise =
        V.getElems   s >>= \es ->
        withSTColViews a $ \xs ->
            sequence_ [ V.scaleByM_ e x
                      | (e,x) <- zip es xs
                      ]
  where
    (m,n) = dim a

-- | @rank1UpdateTo dst alpha x y a@ sets @dst := alpha * x * y^H + a@.
-- Arguments @dst@ and @a@ can be the same; other forms of aliasing give
-- undefined results.
rank1UpdateTo :: (RVector v1, RVector v2, RMatrix m, BLAS2 e)
              => STMatrix s e -> e -> v1 e -> v2 e -> m e -> ST s ()
rank1UpdateTo dst alpha x y a
    | V.dim x /= m || V.dim y /= n || (mdst,ndst) /= (ma,na) = error $
        printf ("rank1UpdateTo"
                ++ "<matrix with dim (%d,%d)>"
                ++ " _"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>"
                ++ ": dimension mismatch"
                ++ "<matrix with dim (%d,%d)>"                )
                mdst ndst
                (V.dim x)
                (V.dim y)
                ma na 
    | otherwise = do
        unsafeCopyTo dst a
        unsafeIOToST $
            V.unsafeWith x $ \px ->
            V.unsafeWith y $ \py ->
            unsafeWith dst $ \pdst lddst ->
                BLAS.gerc m n alpha px 1 py 1 pdst lddst
  where
    (mdst,ndst) = dim dst
    (ma,na) = dim a
    (m,n) = dim a


-- | @transTo dst a@ sets @dst := trans(a)@.
transTo :: (RMatrix m, BLAS1 e)
        => STMatrix s e
        -> m e
        -> ST s ()
transTo a' a
    | (ma,na) /= (na',ma') = error $
        printf ( "transTo"
               ++ " <matrix with dim (%d,%d)>"
               ++ " <matrix with dim (%d,%d)>"
               ++ ": dimension mismatch"
               )
               ma' na'
               ma na
    | otherwise = unsafeIOToST $
        unsafeWith a' $ \pa' lda' ->
        unsafeWith a $ \pa lda -> let
            go j px py | j == n = return ()
                       | otherwise = do
                           BLAS.copy m px 1 py lda'
                           go (j+1) (px `advancePtr` lda) (py `advancePtr` 1)
            in go 0 pa pa'
  where
    (ma,na) = dim a
    (ma',na') = dim a'
    (m,n) = (ma,na)

-- | @conjTransTo dst a@ sets @dst := conjugate(trans(a))@.
conjTransTo :: (RMatrix m, BLAS1 e)
            => STMatrix s e
            -> m e
            -> ST s ()
conjTransTo a' a = do
    transTo a' a
    conjugateTo a' a'

-- | @mulVectorTo dst transa a x@
-- sets @dst := op(a) * x@, where @op(a)@ is determined by @transa@.                   
mulVectorTo :: (RMatrix m, RVector v, BLAS2 e)
            => STVector s e
            -> Trans -> m e
            -> v e
            -> ST s ()
mulVectorTo dst = mulVectorWithScaleTo dst 1

-- | @mulVectorWithScaleTo dst alpha transa a x@
-- sets @dst := alpha * op(a) * x@, where @op(a)@ is determined by @transa@.                   
mulVectorWithScaleTo :: (RMatrix m, RVector v, BLAS2 e)
                     => STVector s e
                     -> e
                     -> Trans -> m e
                     -> v e
                     -> ST s ()
mulVectorWithScaleTo dst alpha t a x =
    mulAddVectorWithScalesTo dst alpha t a x 0 dst

-- | @mulAddVectorWithScalesTo dst alpha transa a x beta y@
-- sets @dst := alpha * op(a) * x + beta * y@, where @op(a)@ is
-- determined by @transa@.  Arguments @dst@ and @y@ can be the same;
-- other forms of aliasing give undefined results.
mulAddVectorWithScalesTo :: (RMatrix m, RVector v1, RVector v2, BLAS2 e)
                         => STVector s e
                         -> e
                         -> Trans -> m e
                         -> v1 e
                         -> e
                         -> v2 e
                         -> ST s ()
mulAddVectorWithScalesTo dst alpha transa a x beta y
    | (not . and) [ case transa of NoTrans -> (ma,na) == (m,n)
                                   _       -> (ma,na) == (n,m)
                  , nx == n
                  , ny == m
                  , ndst == ny
                  ] = error $
        printf ("mulAddVectorWithScalesTo"
                ++ " <vector with dim %d>"        
                ++ " _"
                ++ " %s"
                ++ " <matrix with dim (%d,%d)>" 
                ++ " <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>"
                ++ ": dimension mismatch")
               ndst
               (show transa)
               ma na
               nx
               ny
    | otherwise = do
        V.unsafeCopyTo dst y
        unsafeIOToST $
            unsafeWith a $ \pa lda ->
            V.unsafeWith x $ \px ->
            V.unsafeWith dst $ \pdst ->
                if n == 0
                    then BLAS.scal m beta pdst 1
                    else BLAS.gemv transa ma na alpha pa lda px 1 beta pdst 1
  where
    ndst = V.dim dst
    (ma,na) = dim a
    nx = V.dim x
    ny = V.dim y
    (m,n) = (ny,nx)

-- | @mulMatrixTo dst transa a transb b@
-- sets @dst := op(a) * op(b)@, where @op(a)@ and @op(b)@ are determined
-- by @transa@ and @transb@.                   
mulMatrixTo :: (RMatrix m1, RMatrix m2, BLAS3 e)
            => STMatrix s e
            -> Trans -> m1 e
            -> Trans -> m2 e
            -> ST s ()
mulMatrixTo dst = mulMatrixWithScaleTo dst 1

-- | @mulMatrixWithScaleTo alpha transa a transb b c@
-- sets @c := alpha * op(a) * op(b)@, where @op(a)@ and @op(b)@ are determined
-- by @transa@ and @transb@.                   
mulMatrixWithScaleTo :: (RMatrix m1, RMatrix m2, BLAS3 e)
                     => STMatrix s e
                     -> e
                     -> Trans -> m1 e
                     -> Trans -> m2 e
                     -> ST s ()
mulMatrixWithScaleTo dst alpha ta a tb b =
    mulAddMatrixWithScalesTo dst alpha ta a tb b 0 dst

-- | @mulAddMatrixWithScalesTo dst alpha transa a transb b beta c@
-- sets @dst := alpha * op(a) * op(b) + beta * c@, where @op(a)@ and
-- @op(b)@ are determined by @transa@ and @transb@.  Arguments @dst@ and
-- @c@ can be the same; other forms of aliasing give undefined results.
mulAddMatrixWithScalesTo :: (RMatrix m1, RMatrix m2, RMatrix m3, BLAS3 e)
                         => STMatrix s e
                         -> e
                         -> Trans -> m1 e
                         -> Trans -> m2 e
                         -> e
                         -> m3 e
                         -> ST s ()
mulAddMatrixWithScalesTo dst alpha transa a transb b beta c
    | (not . and) [ case transa of NoTrans -> (ma,na) == (m,k)
                                   _       -> (ma,na) == (k,m)
                  , case transb of NoTrans -> (mb,nb) == (k,n)
                                   _       -> (mb,nb) == (n,k)
                  , (mc, nc) == (m,n)
                  , (mdst, ndst) == (mc, nc)
                  ] = error $
        printf ("mulAddMatrixWithScalesTo"
                ++ " <matrix with dim (%d,%d)>"        
                ++ " _"
                ++ " %s <matrix with dim (%d,%d)>" 
                ++ " %s <matrix with dim (%d,%d)>"
                ++ " _"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
                mdst ndst
               (show transa) ma na
               (show transb) mb nb
               mc nc
    | otherwise = do
        unsafeCopyTo dst c
        unsafeIOToST $
            unsafeWith a $ \pa lda ->
            unsafeWith b $ \pb ldb ->
            unsafeWith dst $ \pdst lddst ->
                BLAS.gemm transa transb m n k alpha pa lda pb ldb beta pdst lddst
  where
    (ma,na) = dim a
    (mb,nb) = dim b
    (mc,nc) = dim c
    (mdst,ndst) = dim dst
    (m,n) = dim c
    k = case transa of NoTrans -> na
                       _       -> ma


checkOp2 :: (RMatrix x, RMatrix y, Storable e, Storable f)
         => String
         -> (x e -> y f -> a)
         -> x e
         -> y f
         -> a
checkOp2 str f x y
    | (m1,n1) /= (m2,n2) = error $
        printf ("%s <matrix with dim (%d,%d)> <matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") str m1 n1 m2 n2
    | otherwise =
        f x y
  where
    (m1,n1) = dim x
    (m2,n2) = dim y        
{-# INLINE checkOp2 #-}

checkOp3 :: (RMatrix x, RMatrix y, RMatrix z, Storable e, Storable f, Storable g)
         => String
         -> (x e -> y f -> z g -> a)
         -> x e
         -> y f
         -> z g
         -> a
checkOp3 str f x y z
    | (m1,n1) /= (m2,n2) || (m1,n1) /= (m3,n3) = error $
        printf ("%s <matrix with dim (%d,%d)> <matrix with dim (%d,%d)>:"
                ++ " <matrix with dim (%d,%d)> dimension mismatch")
               str m1 n1 m2 n2 m3 n3
    | otherwise =
        f x y z
  where
    (m1,n1) = dim x
    (m2,n2) = dim y
    (m3,n3) = dim z
{-# INLINE checkOp3 #-}

vectorOp :: (Storable e)
         => (STVector s e -> ST s ())
         -> STMatrix s e -> ST s ()
vectorOp f x =
    fromMaybe colwise $ maybeWithSTVectorView x $ \vx -> f vx
  where
    colwise = withSTColViews x $ \vxs ->
                  sequence_ [ f vx | vx <- vxs ]

vectorOp2 :: (RMatrix m, Storable e, Storable f)
          => (forall v . RVector v => STVector s f -> v e -> ST s ())
          -> STMatrix s f -> m e -> ST s ()
vectorOp2 f dst x =
    fromMaybe colwise $ maybeWithSTVectorView dst $ \vdst ->
    fromMaybe colwise $ maybeWithVectorView   x   $ \vx ->
        f vdst vx
  where
    colwise = withSTColViews dst $ \vdsts ->
              withColViews   x   $ \vxs ->
                  sequence_ [ f vdst vx | (vdst,vx) <- zip vdsts vxs ]
{-# INLINE vectorOp2 #-}

vectorOp3 :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
          => (forall v1 v2 . (RVector v1, RVector v2) => 
                  STVector s f -> v1 e1 -> v2 e2 -> ST s ())
          -> STMatrix s f -> m1 e1 -> m2 e2 -> ST s ()
vectorOp3 f dst x y =
    fromMaybe colwise $ maybeWithSTVectorView dst $ \vdst ->
    fromMaybe colwise $ maybeWithVectorView   x $ \vx ->
    fromMaybe colwise $ maybeWithVectorView   y $ \vy ->
        f vdst vx vy
  where
    colwise = withSTColViews dst $ \vdsts ->
              withColViews   x   $ \vxs ->
              withColViews   y   $ \vys ->
                  sequence_ [ f vdst vx vy
                            | (vdst,vx,vy) <- zip3 vdsts vxs vys ]
{-# INLINE vectorOp3 #-}

newResult :: (RMatrix m, Storable e, Storable f)
          => (STMatrix s f -> m e -> ST s a)
          -> m e
          -> ST s (STMatrix s f)
newResult f a = do
    c <- new_ (dim a)
    _ <- f c a
    return c
{-# INLINE newResult #-}

newResult2 :: (RMatrix m1, RMatrix m2, Storable e, Storable f, Storable g)
           => (STMatrix s g -> m1 e -> m2 f -> ST s a)
           -> m1 e
           -> m2 f
           -> ST s (STMatrix s g)
newResult2 f a1 a2 = do
    c <- new_ (dim a1)
    _ <- f c a1 a2
    return c
{-# INLINE newResult2 #-}
