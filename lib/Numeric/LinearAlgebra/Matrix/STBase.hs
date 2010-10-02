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

-- | Perform an action with a list of views of the matrix columns.
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
-- view.  This only succeeds when 'maybeSTVectorView' does.
maybeWithSTVectorView :: (Storable e)
                      => STMatrix s e
                      -> (STVector s e -> ST s a)
                      -> Maybe (ST s a)
maybeWithSTVectorView a f = f `fmap` maybeSTVectorView a
{-# INLINE maybeWithSTVectorView #-}

-- | Cast a vector to a matrix of the given shape.  See also
-- 'withViewFromVector'.
viewFromSTVector :: (Storable e)
                 => (Int,Int)
                 -> STVector s e
                 -> STMatrix s e
viewFromSTVector mn@(m,n) v
    | V.dim v /= m*n = error $
        printf ("viewFromSTVector (%d,%d) <vector with dim %d>:"
                ++ " dimension mismatch") m n (V.dim v)
    | otherwise =
        cast v
  where
    cast x = let
        (fptr,o,_) = V.unsafeToForeignPtr x
        lda = max 1 m
        in unsafeFromForeignPtr fptr o mn lda
{-# INLINE viewFromSTVector #-}

-- | Cast a vector to a matrix with one column.  See also
-- 'withViewFromCol'.
viewFromSTCol :: (Storable e)
              => STVector s e
              -> STMatrix s e
viewFromSTCol v = viewFromSTVector (V.dim v, 1) v
{-# INLINE viewFromSTCol #-}

-- | Cast a vector to a matrix with one row.  See also
-- 'withViewFromRow'.
viewFromSTRow :: (Storable e)
              => STVector s e
              -> STMatrix s e
viewFromSTRow v = viewFromSTVector (V.dim v, 1) v
{-# INLINE viewFromSTRow #-}

-- | Cast a vector to a matrix of the given shape and pass it to
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

-- | Cast a vector to a matrix with one column and pass it to
-- the specified function.
withViewFromCol :: (RVector v, Storable e)
                 => v e
                 -> (forall m . RMatrix m => m e -> a)
                 -> a
withViewFromCol v = withViewFromVector (V.dim v, 1) v
{-# INLINE withViewFromCol #-}

-- | Cast a vector to a matrix with one row and pass it to
-- the specified function.
withViewFromRow :: (RVector v, Storable e)
                 => v e
                 -> (forall m . RMatrix m => m e -> a)
                 -> a
withViewFromRow v = withViewFromVector (1, V.dim v) v
{-# INLINE withViewFromRow #-}

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
    unsafeCopyTo a b
    return b

-- | @copyTo src dst@ replaces the values in @dst@ with those in
-- source.  The operands must be the same shape.
copyTo :: (RMatrix m, Storable e) => m e -> STMatrix s e -> ST s ()
copyTo = checkOp2 "copyTo" unsafeCopyTo
{-# INLINE copyTo #-}

unsafeCopyTo :: (RMatrix m, Storable e) => m e -> STMatrix s e -> ST s ()
unsafeCopyTo = vectorOp2 V.unsafeCopyTo
{-# INLINE unsafeCopyTo #-}

-- | Create a view of a matrix by taking the initial rows.
takeRows :: (RMatrix m, Storable e) => Int -> m e -> m e
takeRows i a = slice (0,0) (i,n) a
  where
    (_,n) = dim a

-- | Create a view of a matrix by droping the initial rows.
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

-- | Create a view of a matrix by droping the initial columns.
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

-- | @mapTo f a c@ replaces @c@ elementwise with @f(a)@.
mapTo :: (RMatrix m, Storable e, Storable f)
      => (e -> f)
      -> m e
      -> STMatrix s f
      -> ST s ()
mapTo f = checkOp2 "mapTo _" $ unsafeMapTo f
{-# INLINE mapTo #-}
                            
unsafeMapTo :: (RMatrix m, Storable e, Storable f)
            => (e -> f)
            -> m e
            -> STMatrix s f
            -> ST s ()
unsafeMapTo f a c =
    fromMaybe colwise $ maybeWithVectorView a $ \x ->
    fromMaybe colwise $ maybeWithSTVectorView c $ \z ->
        V.unsafeMapTo f x z
  where
    colwise = withColViews   a $ \xs ->
              withSTColViews c $ \zs ->
                  sequence_ [ V.unsafeMapTo f x z
                            | (x,z) <- zip xs zs
                            ]

-- | @zipWithTo f a b c@ replaces @c@ elementwise with @f(a, b)@.
zipWithTo :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
          => (e1 -> e2 -> f)
          -> m1 e1
          -> m2 e2
          -> STMatrix s f
          -> ST s ()
zipWithTo f = checkOp3 "zipWithTo _" $
    unsafeZipWithTo f
{-# INLINE zipWithTo #-}

unsafeZipWithTo :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
                => (e1 -> e2 -> f)
                -> m1 e1
                -> m2 e2
                -> STMatrix s f
                -> ST s ()
unsafeZipWithTo f a b c =
    fromMaybe colwise $ maybeWithVectorView a $ \x ->
    fromMaybe colwise $ maybeWithVectorView b $ \y ->
    fromMaybe colwise $ maybeWithSTVectorView c $ \z ->
        V.unsafeZipWithTo f x y z
  where
    colwise = withColViews   a $ \xs ->
              withColViews   b $ \ys ->
              withSTColViews c $ \zs ->
                  sequence_ [ V.unsafeZipWithTo f x y z
                            | (x,y,z) <- zip3 xs ys zs
                            ]

-- | Set every element in the matrix to a default value.  For
-- standard numeric types (including 'Double', 'Complex Double', and 'Int'),
-- the default value is '0'.
clear :: (Storable e) => STMatrix s e -> ST s ()
clear a = case maybeSTVectorView a of
    Just x  -> V.clear x
    Nothing -> withSTColViews a $ mapM_ V.clear

-- | Add a constant to all entries of a matrix.
shiftTo :: (RMatrix m, VNum e)
        => e -> m e -> STMatrix s e -> ST s ()
shiftTo e = checkOp2 "shiftTo" $
    vectorOp2 (V.shiftTo e)

-- | Add a vector to the diagonal of a matrix.
shiftDiagTo :: (RVector v, RMatrix m, BLAS1 e)
            => v e -> m e -> STMatrix s e -> ST s ()
shiftDiagTo s a b
    | V.dim s /= mn || dim a /= (m,n) = error $
        printf ("shiftDiagTo <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (V.dim s)
                (fst $ dim a) (snd $ dim a) m n
    | otherwise = shiftDiagToWithScale 1 s a b
  where
    (m,n) = dim b
    mn = min m n

-- | Add a scaled vector to the diagonal of a matrix.
shiftDiagToWithScale :: (RVector v, RMatrix m, BLAS1 e)
                     => e -> v e -> m e -> STMatrix s e -> ST s ()
shiftDiagToWithScale e s a b
    | V.dim s /= mn || dim a /= (m,n) = error $
        printf ("shiftDiagToWithScale _ <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (V.dim s)
                (fst $ dim a) (snd $ dim a) m n
    | otherwise = do
        unsafeCopyTo a b
        unsafeIOToST $
            V.unsafeWith s $ \ps ->
            unsafeWith b $ \pb ldb ->
                BLAS.axpy mn e ps 1 pb (ldb+1)
  where
    (m,n) = dim b
    mn = min m n

-- | Add two matrices.
addTo :: (RMatrix m1, RMatrix m2, VNum e)
      => m1 e -> m2 e -> STMatrix s e -> ST s ()
addTo = checkOp3 "addTo" $ vectorOp3 V.addTo

-- | Add two matrices with the given scales.
addToWithScales :: (RMatrix m1, RMatrix m2, VNum e)
                => e -> m1 e -> e -> m2 e -> STMatrix s e -> ST s ()
addToWithScales alpha a beta b c =
    (checkOp3 "addToWithScales" $
        vectorOp3 (\x y z -> V.addToWithScales alpha x beta y z))
        a b c

-- | Subtract two matrices.
subTo :: (RMatrix m1, RMatrix m2, VNum e)
      => m1 e -> m2 e -> STMatrix s e -> ST s ()
subTo = checkOp3 "subTo" $ vectorOp3 V.subTo

-- | Conjugate the entries of a matrix.
conjTo :: (RMatrix m, VNum e)
       => m e -> STMatrix s e -> ST s ()
conjTo = checkOp2 "conjTo" $
    vectorOp2 V.conjTo

-- | Negate the entries of a matrix.
negateTo :: (RMatrix m, VNum e)
         => m e -> STMatrix s e -> ST s ()
negateTo = checkOp2 "negateTo" $
    vectorOp2 V.negateTo

-- | Scale the entries of a matrix by the given value.
scaleTo :: (RMatrix m, VNum e)
        => e -> m e -> STMatrix s e -> ST s ()
scaleTo e = checkOp2 "scaleTo" $
    vectorOp2 (V.scaleTo e)

-- | Scale the rows of a matrix; @scaleRowsTo s a c@ sets
-- @c := diag(s) * a@.
scaleRowsTo :: (RVector v, RMatrix m, VNum e)
            => v e -> m e -> STMatrix s e -> ST s ()
scaleRowsTo s a b
    | V.dim s /= m || dim a /= (m,n) = error $
        printf ("scaleRowsTo <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (V.dim s)
                (fst $ dim a) (snd $ dim a) m n
    | otherwise =
        withColViews   a $ \xs ->
        withSTColViews b $ \ys ->
            zipWithM_ (V.mulTo s) xs ys
  where
    (m,n) = dim b

-- | Scale the columns of a matrix; @scaleColsTo s a c@ sets
-- @c := a * diag(s)@.
scaleColsTo :: (RVector v, RMatrix m, VNum e)
            => v e -> m e -> STMatrix s e -> ST s ()
scaleColsTo s a b 
    | V.dim s /= n || dim a /= (m,n) = error $
        printf ("scaleColsTo <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)> matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") (V.dim s)
                (fst $ dim a) (snd $ dim a) m n
    | otherwise =
        V.getElems   s >>= \es ->
        withColViews   a $ \xs ->
        withSTColViews b $ \ys ->
            sequence_ [ V.scaleTo e x y
                      | (e,x,y) <- zip3 es xs ys
                      ]
  where
    (m,n) = dim b

-- | @rank1UpdateTo alpha x y a@ sets @a := alpha * x * y^H + a@.
rank1UpdateTo :: (RVector v1, RVector v2, BLAS2 e)
              => e -> v1 e -> v2 e -> STMatrix s e -> ST s ()
rank1UpdateTo alpha x y a
    | V.dim x /= m || V.dim y /= n = error $
        printf ("rank1UpdateTo _ <vector with dim %d>"
                ++ " <vector with dim %d> <matrix with dim (%d,%d)>:"
                ++ " dimension mismatch")
                (V.dim x) (V.dim y) m n
    | otherwise = do
        unsafeIOToST $
            V.unsafeWith x $ \px ->
            V.unsafeWith y $ \py ->
            unsafeWith a $ \pa lda ->
                BLAS.gerc m n alpha px 1 py 1 pa lda
  where
    (m,n) = dim a


-- | @transTo a c@ sets @c := trans(a)@.
transTo :: (RMatrix m, BLAS1 e)
        => m e
        -> STMatrix s e
        -> ST s ()
transTo a a'
    | (ma,na) /= (na',ma') = error $
        printf ("transTo <matrix with dim (%d,%d)>"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
    | otherwise = unsafeIOToST $
        unsafeWith a $ \pa lda ->
        unsafeWith a' $ \pa' lda' -> let
            go j px py | j == n = return ()
                       | otherwise = do
                           BLAS.copy m px 1 py lda'
                           go (j+1) (px `advancePtr` lda) (py `advancePtr` 1)
            in go 0 pa pa'
  where
    (ma,na) = dim a
    (ma',na') = dim a'
    (m,n) = (ma,na)

-- | @conjTransTo a c@ sets @c := conj(trans(a))@.
conjTransTo :: (RMatrix m, BLAS1 e)
            => m e
            -> STMatrix s e
            -> ST s ()
conjTransTo a a' = do
    transTo a a'
    conjTo a' a'

-- | @mulToVector transa a x y@
-- sets @y := op(a) * x@, where @op(a)@ is determined by @transa@.                   
mulToVector :: (RMatrix m, RVector v, BLAS2 e)
            => Trans -> m e
            -> v e
            -> STVector s e
            -> ST s ()
mulToVector = mulToVectorWithScale 1

-- | @mulToVectorWithScale alpha transa a x y@
-- sets @y := alpha * op(a) * x@, where @op(a)@ is determined by @transa@.                   
mulToVectorWithScale :: (RMatrix m, RVector v, BLAS2 e)
                     => e
                     -> Trans -> m e
                     -> v e
                     -> STVector s e
                     -> ST s ()
mulToVectorWithScale alpha t a x y =
    mulAddToVectorWithScales alpha t a x 0 y

-- | @mulAddToVectorWithScales alpha transa a x beta y@
-- sets @y := alpha * op(a) * x + beta * y@, where @op(a)@ is
-- determined by @transa@.
mulAddToVectorWithScales :: (RMatrix m, RVector v, BLAS2 e)
                         => e
                         -> Trans -> m e
                         -> v e
                         -> e
                         -> STVector s e
                         -> ST s ()
mulAddToVectorWithScales alpha transa a x beta y
    | (not . and) [ case transa of NoTrans -> (ma,na) == (m,n)
                                   _       -> (ma,na) == (n,m)
                  , nx == n
                  , ny == m
                  ] = error $
        printf ("mulAddToVectorWithScales _"
                ++ " %s <matrix with dim (%d,%d)>" 
                ++ " <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>: dimension mismatch")
               (show transa) ma na
               nx ny
    | otherwise =
        unsafeIOToST $
            unsafeWith a $ \pa lda ->
            V.unsafeWith x $ \px ->
            V.unsafeWith y $ \py ->
                if n == 0
                    then BLAS.scal m beta py 1
                    else BLAS.gemv transa ma na alpha pa lda px 1 beta py 1
  where
    (ma,na) = dim a
    nx = V.dim x
    ny = V.dim y
    (m,n) = (ny,nx)

-- | @mulToMatrix transa a transb b c@
-- sets @c := op(a) * op(b)@, where @op(a)@ and @op(b)@ are determined
-- by @transa@ and @transb@.                   
mulToMatrix :: (RMatrix m1, RMatrix m2, BLAS3 e)
            => Trans -> m1 e
            -> Trans -> m2 e
            -> STMatrix s e
            -> ST s ()
mulToMatrix = mulToMatrixWithScale 1

-- | @mulToMatrixWithScale alpha transa a transb b c@
-- sets @c := alpha * op(a) * op(b)@, where @op(a)@ and @op(b)@ are determined
-- by @transa@ and @transb@.                   
mulToMatrixWithScale :: (RMatrix m1, RMatrix m2, BLAS3 e)
                     => e
                     -> Trans -> m1 e
                     -> Trans -> m2 e
                     -> STMatrix s e
                     -> ST s ()
mulToMatrixWithScale alpha ta a tb b c =
    mulAddToMatrixWithScales alpha ta a tb b 0 c

-- | @mulAddToMatrixWithScales alpha transa a transb b beta c@
-- sets @c := alpha * op(a) * op(b) + beta * c@, where @op(a)@ and
-- @op(b)@ are determined by @transa@ and @transb@.
mulAddToMatrixWithScales :: (RMatrix m1, RMatrix m2, BLAS3 e)
                               => e
                               -> Trans -> m1 e
                               -> Trans -> m2 e
                               -> e
                               -> STMatrix s e
                               -> ST s ()
mulAddToMatrixWithScales alpha transa a transb b beta c
    | (not . and) [ case transa of NoTrans -> (ma,na) == (m,k)
                                   _       -> (ma,na) == (k,m)
                  , case transb of NoTrans -> (mb,nb) == (k,n)
                                   _       -> (mb,nb) == (n,k)
                  , (mc, nc) == (m,n)
                  ] = error $
        printf ("mulAddToMatrixWithScales _"
                ++ " %s <matrix with dim (%d,%d)>" 
                ++ " %s <matrix with dim (%d,%d)>"
                ++ " _"
                ++ " <matrix with dim (%d,%d)>: dimension mismatch")
               (show transa) ma na
               (show transb) mb nb
               mc nc
    | otherwise =
        unsafeIOToST $
            unsafeWith a $ \pa lda ->
            unsafeWith b $ \pb ldb ->
            unsafeWith c $ \pc ldc ->
                BLAS.gemm transa transb m n k alpha pa lda pb ldb beta pc ldc
  where
    (ma,na) = dim a
    (mb,nb) = dim b
    (mc,nc) = dim c
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

vectorOp2 :: (RMatrix m, Storable e, Storable f)
          => (forall v . RVector v => v e -> STVector s f -> ST s ())
          -> m e -> STMatrix s f -> ST s ()
vectorOp2 f a c =
    fromMaybe colwise $ maybeWithVectorView   a $ \x ->
    fromMaybe colwise $ maybeWithSTVectorView c $ \z ->
        f x z
  where
    colwise = withColViews   a $ \xs ->
              withSTColViews c $ \zs ->
                  sequence_ [ f x z | (x,z) <- zip xs zs ]
{-# INLINE vectorOp2 #-}

vectorOp3 :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
          => (forall v1 v2 . (RVector v1, RVector v2) => 
                  v1 e1 -> v2 e2 -> STVector s f -> ST s ())
          -> m1 e1 -> m2 e2 -> STMatrix s f -> ST s ()
vectorOp3 f a b c =
    fromMaybe colwise $ maybeWithVectorView   a $ \x ->
    fromMaybe colwise $ maybeWithVectorView   b $ \y ->
    fromMaybe colwise $ maybeWithSTVectorView c $ \z ->
        f x y z
  where
    colwise = withColViews   a $ \xs ->
              withColViews   b $ \ys ->
              withSTColViews c $ \zs ->
                  sequence_ [ f x y z | (x,y,z) <- zip3 xs ys zs ]
{-# INLINE vectorOp3 #-}

newResult :: (RMatrix m, Storable e, Storable f)
          => (m e -> STMatrix s f -> ST s a)
          -> m e
          -> ST s (STMatrix s f)
newResult f a = do
    c <- new_ (dim a)
    _ <- f a c
    return c
{-# INLINE newResult #-}

newResult2 :: (RMatrix m1, RMatrix m2, Storable e, Storable f, Storable g)
           => (m1 e -> m2 f -> STMatrix s g -> ST s a)
           -> m1 e
           -> m2 f
           -> ST s (STMatrix s g)
newResult2 f a1 a2 = do
    c <- new_ (dim a1)
    _ <- f a1 a2 c
    return c
{-# INLINE newResult2 #-}
