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
      
import Control.Monad( forM_, when )
import Control.Monad.ST( ST, RealWorld, runST, unsafeInterleaveST,
    unsafeIOToST )
import Data.Maybe( fromMaybe )
import Data.Typeable( Typeable )
import Foreign( Ptr, advancePtr, peek, peekElemOff, pokeElemOff,
    mallocForeignPtrArray )
import Text.Printf( printf )
import Unsafe.Coerce( unsafeCoerce )

import Numeric.LinearAlgebra.Types
import qualified Foreign.BLAS as BLAS
import Numeric.LinearAlgebra.Matrix.Base hiding ( unsafeWith,
    unsafeToForeignPtr, unsafeFromForeignPtr, )
import qualified Numeric.LinearAlgebra.Matrix.Base as M
import Numeric.LinearAlgebra.Vector( STVector, RVector )
import qualified Numeric.LinearAlgebra.Vector as V

-- | Mutable dense matrices in the 'ST' monad.
newtype STMatrix s e = STMatrix { unSTMatrix :: Matrix e }
    deriving (Typeable)

-- | Mutable dense matrices in the 'IO' monad.
type IOMatrix = STMatrix RealWorld

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


-- | Read-only matrices
class RMatrix m where
    -- | Get the dimensions of the matrix (number of rows and columns).
    getDim :: (Storable e) => m e -> ST s (Int,Int)
    
    unsafeWithColView :: (Storable e)
                      => m e
                      -> Int 
                      -> (forall v. RVector v => v e -> ST s a)
                      -> ST s a

    -- | Perform an action with a list of views of the matrix columns.
    withColViews :: (Storable e)
                 => m e
                 -> (forall v . RVector v => [v e] -> ST s a)
                 -> ST s a

    unsafeWithSliceView :: (Storable e)
                        => (Int,Int)
                        -> (Int,Int)
                        -> m e
                        -> (forall m'. RMatrix m' => m' e -> ST s a)
                        -> ST s a

    -- | Possibly view a matrix as a vector and perform an action on the
    -- view.  This only succeeds if the matrix is stored contiguously in
    -- memory, i.e. if the matrix contains a single column or the \"lda\"
    -- of the matrix is equal to the number of rows.
    maybeWithVectorView :: (Storable e)
                        => m e
                        -> (forall v . RVector v => v e -> ST s a)
                        -> Maybe (ST s a)

    -- | Unsafe cast from a read-only matrix to a mutable matrix.
    unsafeThaw :: (Storable e)
               => m e -> ST s (STMatrix s e)

    -- | Execute an 'IO' action with a pointer to the first element in the
    -- matrix and the leading dimension (lda).
    unsafeWith :: (Storable e) => m e -> (Ptr e -> Int -> IO a) -> IO a

instance RMatrix Matrix where
    getDim = return . dim
    {-# INLINE getDim #-}
    unsafeWithColView a j f = f (unsafeCol a j)
    {-# INLINE unsafeWithColView #-}
    withColViews a f = f (cols a)
    {-# INLINE withColViews #-}
    unsafeWithSliceView ij mn a f = f (unsafeSlice ij mn a)
    {-# INLINE unsafeWithSliceView #-}
    maybeWithVectorView a f | isContig a = Just $ f (toVector a)
                            | otherwise = Nothing
    {-# INLINE maybeWithVectorView #-}
    unsafeWith = M.unsafeWith
    {-# INLINE unsafeWith #-}
    unsafeThaw = return . STMatrix
    {-# INLINE unsafeThaw #-}


instance RMatrix (STMatrix s) where
    getDim = return . dim . unSTMatrix
    {-# INLINE getDim #-}
    unsafeWithColView = unsafeWithColView . unSTMatrix
    {-# INLINE unsafeWithColView #-}
    withColViews = withColViews . unSTMatrix
    {-# INLINE withColViews #-}
    unsafeWithSliceView ij mn = unsafeWithSliceView ij mn . unSTMatrix
    {-# INLINE unsafeWithSliceView #-}
    maybeWithVectorView = maybeWithVectorView . unSTMatrix
    {-# INLINE maybeWithVectorView #-}
    unsafeWith = unsafeWith . unSTMatrix
    {-# INLINE unsafeWith #-}
    unsafeThaw v = return $ cast v
      where
        cast :: STMatrix s e -> STMatrix s' e
        cast = unsafeCoerce
    {-# INLINE unsafeThaw #-}



-- | Perform an action with a view of a matrix column (no index checking).
unsafeWithSTColView :: (Storable e)
                    => STMatrix s e
                    -> Int
                    -> (STVector s e -> ST s a)
                    -> ST s a
unsafeWithSTColView a j f = 
    unsafeWithColView a j $ \c -> do
        mc <- V.unsafeThaw c
        f mc
{-# INLINE unsafeWithSTColView #-}

-- | Perform an action with a list of views of the matrix columns.  See
-- also 'withColViews'.
withSTColViews :: (Storable e)
               => STMatrix s e
               -> ([STVector s e] -> ST s a)
               -> ST s a
withSTColViews a f =
    withColViews a $ \cs -> do
        mcs <- thawVecs cs
        f mcs
  where
    thawVecs [] = return []
    thawVecs (c:cs) = unsafeInterleaveST $ do
        mc <- V.unsafeThaw c
        mcs <- thawVecs cs
        return $ mc:mcs
{-# INLINE withSTColViews #-}


-- | Possibly view a matrix as a vector and perform an action on the
-- view.  This succeeds when the matrix is stored contiguously in memory,
-- i.e. if the matrix contains a single column or the \"lda\" of the matrix
-- is equal to the number of rows.  See also 'maybeWithVectorView'.
maybeWithSTVectorView :: (Storable e)
                      => STMatrix s e
                      -> (STVector s e -> ST s a)
                      -> Maybe (ST s a)
maybeWithSTVectorView a f = 
    maybeWithVectorView a $ \v -> do
        mv <- V.unsafeThaw v
        f mv
{-# INLINE maybeWithSTVectorView #-}


-- | View a vector as a matrix of the given shape and pass it to
-- the specified function.
withViewFromVector :: (RVector v, Storable e)
                   => (Int,Int)
                   -> v e
                   -> (forall m . RMatrix m => m e -> ST s a)
                   -> ST s a
withViewFromVector mn@(m,n) v f = do
    nv <- V.getDim v
    when (nv /= m*n) $ error $
        printf ("withViewFromVector (%d,%d) <vector with dim %d>:"
                ++ " dimension mismatch") m n nv
    iv <- V.unsafeFreeze v
    f $ fromVector mn iv
{-# INLINE withViewFromVector #-}


-- | View a mutable vector as a mutable matrix of the given shape and pass it
-- to the specified function.
withViewFromSTVector :: (Storable e)
                     => (Int,Int)
                     -> STVector s e
                     -> (STMatrix s e -> ST s a)
                     -> ST s a
withViewFromSTVector mn@(m,n) v f = do
    nv <- V.getDim v
    when (nv /= m*n) $ error $
        printf ("withViewFromSTVector (%d,%d) <vector with dim %d>:"
                ++ " dimension mismatch") m n nv
    withViewFromVector mn v $ \a -> do
        ma <- unsafeThaw a
        f ma
{-# INLINE withViewFromSTVector #-}


-- | View a vector as a matrix with one column and pass it to
-- the specified function.
withViewFromCol :: (RVector v, Storable e)
                => v e
                -> (forall m . RMatrix m => m e -> ST s a)
                -> ST s a
withViewFromCol v f = do
    m <- V.getDim v
    withViewFromVector (m,1) v f
{-# INLINE withViewFromCol #-}


-- | View a mutable vector as a mutable matrix with one column and pass it to
-- the specified function.
withViewFromSTCol :: (Storable e)
                  => STVector s e
                  -> (STMatrix s e -> ST s a)
                  -> ST s a
withViewFromSTCol v f = do
    m <- V.getDim v
    withViewFromSTVector (m, 1) v f
{-# INLINE withViewFromSTCol #-}


-- | View a vector as a matrix with one row and pass it to
-- the specified function.
withViewFromRow :: (RVector v, Storable e)
                => v e
                -> (forall m . RMatrix m => m e -> ST s a)
                -> ST s a
withViewFromRow v f = do
    n <- V.getDim v
    withViewFromVector (1,n) v f
{-# INLINE withViewFromRow #-}

-- | View a mutable vector as a mutable matrix with one row and pass it to
-- the specified function.
withViewFromSTRow :: (Storable e)
                  => STVector s e
                  -> (STMatrix s e -> ST s a)
                  -> ST s a
withViewFromSTRow v f = do
    n <- V.getDim v
    withViewFromSTVector (1,n) v f
{-# INLINE withViewFromSTRow #-}

-- | Perform an action with a view of a matrix column.
withColView :: (RMatrix m, Storable e)
            => m e
            -> Int
            -> (forall v . RVector v => v e -> ST s a)
            -> ST s a
withColView a j f = do
    (m,n) <- getDim a
    when (j < 0 || j >= n) $ error $
        printf ("withColView <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n j

    unsafeWithColView a j f
{-# INLINE withColView #-}

-- | Perform an action with a view of a mutable matrix column.
withSTColView :: (Storable e)
              => STMatrix s e
              -> Int
              -> (STVector s e -> ST s a)
              -> ST s a
withSTColView a j f = do
    (m,n) <- getDim a
    when (j < 0 || j >= n) $ error $
        printf ("withSTColView <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n j

    unsafeWithSTColView a j f
{-# INLINE withSTColView #-}

-- | Create a new matrix of given shape, but do not initialize the elements.
new_ :: (Storable e) => (Int,Int) -> ST s (STMatrix s e)
new_ (m,n) 
    | m < 0 || n < 0 = error $
        printf "new_ (%d,%d): invalid dimensions" m n
    | otherwise = unsafeIOToST $ do
        f <- mallocForeignPtrArray (m*n)
        return $ STMatrix $ M.unsafeFromForeignPtr f 0 (m,n) (max 1 m)

-- | Create a matrix with every element initialized to the same value.
new :: (Storable e) => (Int,Int) -> e -> ST s (STMatrix s e)
new (m,n) e = do
    a <- new_ (m,n)
    setElems a $ replicate (m*n) e
    return a

-- | Creates a new matrix by copying another one.    
newCopy :: (RMatrix m, Storable e) => m e -> ST s (STMatrix s e)
newCopy a = do
    mn <- getDim a
    b <- new_ mn
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

-- | Get the indices of the elements in the matrix, in column-major order.
getIndices :: (RMatrix m, Storable e) => m e -> ST s [(Int,Int)]
getIndices a = do
    (m,n) <- getDim a
    return $ [ (i,j) | j <- [ 0..n-1 ], i <- [ 0..m-1 ] ]
  
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
    is <- getIndices a
    es <- getElems a
    return $ zip is es

-- | Get the association list of the matrix, in column-major order.
getAssocs' :: (RMatrix m, Storable e) => m e -> ST s [((Int,Int),e)]
getAssocs' a = do
    is <- getIndices a
    es <- getElems' a
    return $ zip is es

-- | Set all of the values of the matrix from the elements in the list,
-- in column-major order.
setElems :: (Storable e) => STMatrix s e -> [e] -> ST s ()
setElems a es =
    case maybeWithSTVectorView a (`V.setElems` es) of
        Just st  -> st
        Nothing -> do
            (m,n) <- getDim a
            go m n 0 es
  where
    go _ n j [] | j == n = return ()
    go m n j [] | j < n = error $ 
        printf ("setElems <matrix with dim (%d,%d>"
                ++ "<list with length %d>: not enough elements)") m n (j*m)
    go m n j es' =
        let (es1', es2') = splitAt m es'
        in do
            withSTColView a j (`V.setElems` es1')
            go m n (j+1) es2'

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
setRow a i x = do
    (m,n) <- getDim a
    nx <- V.getDim x
    
    when (i < 0 || i >= m) $ error $
        printf ("setRow <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n i
    when (nx /= n) $ error $
        printf ("setRow <matrix with dim (%d,%d)> _"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n nx

    unsafeSetRow a i x
{-# INLINE setRow #-}

unsafeSetRow :: (RVector v, Storable e)
             => STMatrix s e -> Int -> v e -> ST s ()
unsafeSetRow a i x = do
    jes <- V.getAssocs x
    sequence_ [ unsafeWrite a (i,j) e | (j,e) <- jes ]
{-# INLINE unsafeSetRow #-}

-- | Copy the specified row of the matrix to the vector.
rowTo :: (RMatrix m, Storable e)
         => STVector s e -> m e -> Int -> ST s ()
rowTo x a i = do
    (m,n) <- getDim a
    nx <- V.getDim x
    when (i < 0 || i >= m) $ error $
        printf ("rowTo"
               ++ " _"
               ++ " <matrix with dim (%d,%d)>"
               ++ " %d:"
               ++ ": index out of range"
               ) m n i
    when (nx /= n) $ error $
        printf ("rowTo"
                ++ " <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)>"
                ++ " _"
                ++ ": dimension mismatch") nx m n

    unsafeRowTo x a i
{-# INLINE rowTo #-}

unsafeRowTo :: (RMatrix m, Storable e)
             =>  STVector s e -> m e -> Int ->ST s ()
unsafeRowTo x a i = do
    (_,n) <- getDim a
    forM_ [ 0..n-1 ] $ \j -> do
        e <- unsafeRead a (i,j)
        V.unsafeWrite x j e
{-# INLINE unsafeRowTo #-}

-- | Set the diagonal of the matrix to the given vector.
setDiag :: (RVector v, Storable e)
        => STMatrix s e -> v e -> ST s ()
setDiag a x = do
    (m,n) <- getDim a
    nx <- V.getDim x
    let mn = min m n
    
    when (nx /= mn) $ error $
        printf ("setRow <matrix with dim (%d,%d)>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n nx

    unsafeSetDiag a x
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
getDiag a x = do
    (m,n) <- getDim a
    nx <- V.getDim x
    let mn = min m n
    
    when (nx /= mn) $ error $
        printf ("getDiag <matrix with dim (%d,%d)>"
                ++ " <vector with dim %d>:"
                ++ " dimension mismatch") m n nx

    unsafeGetDiag a x
{-# INLINE getDiag #-}

unsafeGetDiag :: (RMatrix m, Storable e)
              => m e -> STVector s e -> ST s ()
unsafeGetDiag a x = do
    (m,n) <- getDim a
    let mn = min m n
    forM_ [ 0..mn-1 ] $ \i -> do
        e <- unsafeRead a (i,i)
        V.unsafeWrite x i e
{-# INLINE unsafeGetDiag #-}

-- | Get the element stored at the given index.
read :: (RMatrix m, Storable e) => m e -> (Int,Int) -> ST s e
read a (i,j) = do
    (m,n) <- getDim a
    when (i < 0 || i >= m || j < 0 || j >= n) $ error $
        printf ("read <matrix with dim (%d,%d)> (%d,%d):"
                ++ " index out of range") m n i j
    unsafeRead a (i,j)
{-# INLINE read #-}

unsafeRead :: (RMatrix m, Storable e) => m e -> (Int,Int) -> ST s e
unsafeRead a (i,j) = unsafeIOToST $
    unsafeWith a $ \p lda ->
        peekElemOff p (i + j * lda)
{-# INLINE unsafeRead #-}

-- | Set the element stored at the given index.
write :: (Storable e)
      => STMatrix s e -> (Int,Int) -> e -> ST s ()
write a (i,j) e = do
    (m,n) <- getDim a
    when (i < 0 || i >= m || j < 0 || j >= n) $ error $
        printf ("write <matrix with dim (%d,%d)> (%d,%d):"
                ++ " index out of range") m n i j
    unsafeWrite a (i,j) e
    
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
modify a (i,j) f = do
    (m,n) <- getDim a
    when (i < 0 || i >= m || j < 0 || j >= n) $ error $
        printf ("modify <matrix with dim (%d,%d)> (%d,%d):"
                ++ " index out of range") m n i j
    unsafeModify a (i,j) f
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
clear a = fromMaybe colwise $ maybeWithSTVectorView a V.clear
  where
    colwise = withSTColViews a $ mapM_ V.clear

-- | Add a vector to the diagonal of a matrix.
shiftDiagByM_ :: (RVector v, BLAS1 e)
              => v e -> STMatrix s e -> ST s ()
shiftDiagByM_ s a = do
    (m,n) <- getDim a
    ns <- V.getDim s
    let mn = min m n
    
    when (ns /= mn) $ error $
        printf ("shiftDiagByM_"
                ++ " <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
                ns
                m n
    
    shiftDiagByWithScaleM_ 1 s a
    

-- | Add a scaled vector to the diagonal of a matrix.
shiftDiagByWithScaleM_ :: (RVector v, BLAS1 e)
                       => e -> v e -> STMatrix s e -> ST s ()
shiftDiagByWithScaleM_ e s a = do
    (m,n) <- getDim a
    ns <- V.getDim s
    let mn = min m n

    when (ns /= mn) $ error $
        printf ("shiftDiagByWithScaleM_"
                ++ " _"
                ++ " <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
                ns
                m n

    unsafeIOToST $
        V.unsafeWith s $ \ps ->
        unsafeWith a $ \pa lda ->
            BLAS.axpy mn e ps 1 pa (lda+1)


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

-- | Scale the rows of a matrix; @scaleRowsByM_ s a@ sets
-- @a := diag(s) * a@.
scaleRowsByM_ :: (RVector v, BLAS1 e)
              => v e -> STMatrix s e -> ST s ()
scaleRowsByM_  s a = do
    (m,n) <- getDim a
    ns <- V.getDim s
    when (ns /= m) $ error $
        printf ("scaleRowsByM_"
                ++ " <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
                ns
                m n

    unsafeIOToST $
        V.unsafeWith s $ \ps ->
        unsafeWith a   $ \pa lda ->
            go m n lda pa ps 0
  where
    go m n lda pa ps i | i == m    = return ()
                       | otherwise = do
                             e <- peek ps
                             BLAS.scal n e pa lda
                             go m n lda (pa `advancePtr` 1)
                                        (ps `advancePtr` 1)
                                        (i+1)

-- | Scale the columns of a matrix; @scaleColBysM_ s a@ sets
-- @a := a * diag(s)@.
scaleColsByM_ :: (RVector v, BLAS1 e)
            => v e -> STMatrix s e -> ST s ()
scaleColsByM_ s a = do
    (m,n) <- getDim a
    ns <- V.getDim s
    when (ns /= n) $ error $
        printf ("scaleColsByM_"
                ++ " <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)>"        
                ++ ": dimension mismatch") 
                ns
                m n

    es <- V.getElems s
    withSTColViews a $ \xs ->
        sequence_ [ V.scaleByM_ e x
                  | (e,x) <- zip es xs
                  ]


-- | @rank1UpdateM_ alpha x y a@ sets @a := alpha * x * y^H + a@.
rank1UpdateM_ :: (RVector v1, RVector v2, BLAS2 e)
              => e -> v1 e -> v2 e -> STMatrix s e -> ST s ()
rank1UpdateM_ alpha x y a = do
    (m,n) <- getDim a    
    nx <- V.getDim x
    ny <- V.getDim y
    
    when (nx /= m || ny /= n) $ error $
        printf ("rank1UpdateTo"
                ++ " _"
                ++ " <vector with dim %d>"
                ++ " <vector with dim %d>"
                ++ ": dimension mismatch"
                ++ "<matrix with dim (%d,%d)>"                )
                nx
                ny
                m n
    
    unsafeIOToST $
        V.unsafeWith x $ \px ->
        V.unsafeWith y $ \py ->
        unsafeWith a $ \pa lda ->
            BLAS.gerc m n alpha px 1 py 1 pa lda


-- | @transTo dst a@ sets @dst := trans(a)@.
transTo :: (RMatrix m, BLAS1 e)
        => STMatrix s e
        -> m e
        -> ST s ()
transTo a' a = do
    (ma,na) <- getDim a
    (ma',na') <- getDim a'
    let (m,n) = (ma,na)

    when ((ma,na) /= (na',ma')) $ error $
        printf ( "transTo"
               ++ " <matrix with dim (%d,%d)>"
               ++ " <matrix with dim (%d,%d)>"
               ++ ": dimension mismatch"
               )
               ma' na'
               ma na
    
    unsafeIOToST $
        unsafeWith a' $ \pa' lda' ->
        unsafeWith a $ \pa lda -> let
            go j px py | j == n = return ()
                       | otherwise = do
                           BLAS.copy m px 1 py lda'
                           go (j+1) (px `advancePtr` lda) (py `advancePtr` 1)
            in go 0 pa pa'


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
    addMulVectorWithScalesM_ alpha t a x 0 dst

-- | @addMulVectorWithScalesM_ alpha transa a x beta y@
-- sets @y := alpha * op(a) * x + beta * y@, where @op(a)@ is
-- determined by @transa@.
addMulVectorWithScalesM_ :: (RMatrix m, RVector v, BLAS2 e)
                         => e
                         -> Trans -> m e
                         -> v e
                         -> e
                         -> STVector s e
                         -> ST s ()
addMulVectorWithScalesM_ alpha transa a x beta y = do
    (ma,na) <- getDim a
    nx <- V.getDim x
    ny <- V.getDim y
    let (m,n) = (ny,nx)

    when ((not . and) [ case transa of NoTrans -> (ma,na) == (m,n)
                                       _       -> (ma,na) == (n,m)
                      , nx == n
                      , ny == m
                      ]) $ error $
        printf ("addMulVectorWithScalesTo"
                ++ " _"
                ++ " %s"
                ++ " <matrix with dim (%d,%d)>" 
                ++ " <vector with dim %d>"
                ++ " _"
                ++ " <vector with dim %d>"
                ++ ": dimension mismatch")
               (show transa)
               ma na
               nx
               ny

    unsafeIOToST $
        unsafeWith a $ \pa lda ->
        V.unsafeWith x $ \px ->
        V.unsafeWith y $ \py ->
            if n == 0
                then BLAS.scal m beta py 1
                else BLAS.gemv transa ma na alpha pa lda px 1 beta py 1

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
    addMulMatrixWithScalesM_ alpha ta a tb b 0 dst

-- | @addMulMatrixWithScalesM_ alpha transa a transb b beta c@
-- sets @c := alpha * op(a) * op(b) + beta * c@, where @op(a)@ and
-- @op(b)@ are determined by @transa@ and @transb@.
addMulMatrixWithScalesM_ :: (RMatrix m1, RMatrix m2, BLAS3 e)
                         => e
                         -> Trans -> m1 e
                         -> Trans -> m2 e
                         -> e
                         -> STMatrix s e
                         -> ST s ()
addMulMatrixWithScalesM_ alpha transa a transb b beta c = do
    (ma,na) <- getDim a
    (mb,nb) <- getDim b
    (mc,nc) <- getDim c
    let (m,n) = (mc,nc)
        k = case transa of NoTrans -> na
                           _       -> ma

    when ((not . and) [ case transa of NoTrans -> (ma,na) == (m,k)
                                       _       -> (ma,na) == (k,m)
                      , case transb of NoTrans -> (mb,nb) == (k,n)
                                       _       -> (mb,nb) == (n,k)
                      , (mc, nc) == (m,n)
                      ]) $ error $
        printf ("addMulMatrixWithScalesM_"
                ++ " _"
                ++ " %s <matrix with dim (%d,%d)>" 
                ++ " %s <matrix with dim (%d,%d)>"
                ++ " _"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
               (show transa) ma na
               (show transb) mb nb
               mc nc

    unsafeIOToST $
        unsafeWith a $ \pa lda ->
        unsafeWith b $ \pb ldb ->
        unsafeWith c $ \pc ldc ->
            BLAS.gemm transa transb m n k alpha pa lda pb ldb beta pc ldc

checkOp2 :: (RMatrix x, RMatrix y, Storable e, Storable f)
         => String
         -> (x e -> y f -> ST s a)
         -> x e
         -> y f
         -> ST s a
checkOp2 str f x y = do
    (m1,n1) <- getDim x
    (m2,n2) <- getDim y
    when ((m1,n1) /= (m2,n2)) $ error $
        printf ("%s <matrix with dim (%d,%d)> <matrix with dim (%d,%d)>:"
                ++ " dimension mismatch") str m1 n1 m2 n2
    f x y
{-# INLINE checkOp2 #-}

checkOp3 :: (RMatrix x, RMatrix y, RMatrix z, Storable e, Storable f, Storable g)
         => String
         -> (x e -> y f -> z g -> ST s a)
         -> x e
         -> y f
         -> z g
         -> ST s a
checkOp3 str f x y z = do
    (m1,n1) <- getDim x
    (m2,n2) <- getDim y
    (m3,n3) <- getDim z
    when((m1,n1) /= (m2,n2) || (m1,n1) /= (m3,n3)) $ error $
        printf ("%s <matrix with dim (%d,%d)> <matrix with dim (%d,%d)>:"
                ++ " <matrix with dim (%d,%d)> dimension mismatch")
               str m1 n1 m2 n2 m3 n3
    f x y z
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
