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
import Control.Monad.ST( ST, RealWorld, runST )
import Control.Monad.ST.Unsafe( unsafeInterleaveST, unsafeIOToST )
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
freeze :: (RMatrix m, Storable e) => m e -> ST s (Matrix e)
freeze a = do
    a' <- newCopy a
    unsafeFreeze a'
{-# INLINE freeze #-}


-- | Read-only matrices
class RMatrix m where
    -- | Get the dimensions of the matrix (number of rows and columns).
    getDim :: (Storable e) => m e -> ST s (Int,Int)
    
    -- | Same as 'withCol' but does not range-check index.
    unsafeWithCol :: (Storable e)
                  => m e
                  -> Int 
                  -> (forall v. RVector v => v e -> ST s a)
                  -> ST s a

    -- | Perform an action with a list of views of the matrix columns.
    withCols :: (Storable e)
                 => m e
                 -> (forall v . RVector v => [v e] -> ST s a)
                 -> ST s a

    -- | Same as 'withSlice' but does not range-check index.
    unsafeWithSlice :: (Storable e)
                    => (Int,Int)
                    -> (Int,Int)
                    -> m e
                    -> (forall m'. RMatrix m' => m' e -> ST s a)
                    -> ST s a

    -- | Possibly view a matrix as a vector and perform an action on the
    -- view.  This only succeeds if the matrix is stored contiguously in
    -- memory, i.e. if the matrix contains a single column or the \"lda\"
    -- of the matrix is equal to the number of rows.
    maybeWithVector :: (Storable e)
                    => m e
                    -> (forall v . RVector v => v e -> ST s a)
                    -> Maybe (ST s a)

    -- | Converts a read-only matrix into an immutable matrix. This simply
    -- casts the matrix from one type to the other without copying.
    -- Note that because the matrix is possibly not copied, any subsequent
    -- modifications made to the read-only version of the matrix may be shared
    -- with the immutable version. It is safe to use, therefore, if the
    -- read-only version is never modified after the freeze operation.
    unsafeFreeze :: (Storable e) => m e -> ST s (Matrix e)

    -- | Unsafe cast from a read-only matrix to a mutable matrix.
    unsafeThaw :: (Storable e)
               => m e -> ST s (STMatrix s e)

    -- | Execute an 'IO' action with a pointer to the first element in the
    -- matrix and the leading dimension (lda).
    unsafeWith :: (Storable e) => m e -> (Ptr e -> Int -> IO a) -> IO a

instance RMatrix Matrix where
    getDim = return . dim
    {-# INLINE getDim #-}
    unsafeWithCol a j f = f (unsafeCol a j)
    {-# INLINE unsafeWithCol #-}
    withCols a f = f (cols a)
    {-# INLINE withCols #-}
    unsafeWithSlice ij mn a f = f (unsafeSlice ij mn a)
    {-# INLINE unsafeWithSlice #-}
    maybeWithVector a f | isContig a = Just $ f (toVector a)
                        | otherwise = Nothing
    {-# INLINE maybeWithVector #-}
    unsafeWith = M.unsafeWith
    {-# INLINE unsafeWith #-}
    unsafeFreeze = return
    {-# INLINE unsafeFreeze #-}
    unsafeThaw = return . STMatrix
    {-# INLINE unsafeThaw #-}


instance RMatrix (STMatrix s) where
    getDim = return . dim . unSTMatrix
    {-# INLINE getDim #-}
    unsafeWithCol st = unsafeWithCol $ unSTMatrix st
    {-# INLINE unsafeWithCol #-}
    withCols st = withCols $ unSTMatrix st
    {-# INLINE withCols #-}
    unsafeWithSlice ij mn st = unsafeWithSlice ij mn $ unSTMatrix st
    {-# INLINE unsafeWithSlice #-}
    maybeWithVector st = maybeWithVector $ unSTMatrix st
    {-# INLINE maybeWithVector #-}
    unsafeWith = unsafeWith . unSTMatrix
    {-# INLINE unsafeWith #-}
    unsafeFreeze = return . unSTMatrix
    {-# INLINE unsafeFreeze #-}
    unsafeThaw v = return $ cast v
      where
        cast :: STMatrix s e -> STMatrix s' e
        cast = unsafeCoerce
    {-# INLINE unsafeThaw #-}



-- | Perform an action with a view of a mutable matrix column
-- (no index checking).
unsafeWithColM :: (Storable e)
               => STMatrix s e
               -> Int
               -> (STVector s e -> ST s a)
               -> ST s a
unsafeWithColM a j f = 
    unsafeWithCol a j $ \c -> do
        mc <- V.unsafeThaw c
        f mc
{-# INLINE unsafeWithColM #-}

-- | Perform an action with a list of views of the mutable matrix columns. See
-- also 'withCols'.
withColsM :: (Storable e)
               => STMatrix s e
               -> ([STVector s e] -> ST s a)
               -> ST s a
withColsM a f =
    withCols a $ \cs -> do
        mcs <- thawVecs cs
        f mcs
  where
    thawVecs [] = return []
    thawVecs (c:cs) = unsafeInterleaveST $ do
        mc <- V.unsafeThaw c
        mcs <- thawVecs cs
        return $ mc:mcs
{-# INLINE withColsM #-}


-- | Possibly view a matrix as a vector and perform an action on the
-- view.  This succeeds when the matrix is stored contiguously in memory,
-- i.e. if the matrix contains a single column or the \"lda\" of the matrix
-- is equal to the number of rows.  See also 'maybeWithVector'.
maybeWithVectorM :: (Storable e)
                 => STMatrix s e
                 -> (STVector s e -> ST s a)
                 -> Maybe (ST s a)
maybeWithVectorM a f = 
    maybeWithVector a $ \v -> do
        mv <- V.unsafeThaw v
        f mv
{-# INLINE maybeWithVectorM #-}


-- | View a vector as a matrix of the given shape and pass it to
-- the specified function.
withFromVector :: (RVector v, Storable e)
                   => (Int,Int)
                   -> v e
                   -> (forall m . RMatrix m => m e -> ST s a)
                   -> ST s a
withFromVector mn@(m,n) v f = do
    nv <- V.getDim v
    when (nv /= m*n) $ error $
        printf ("withFromVector (%d,%d) <vector with dim %d>:"
                ++ " dimension mismatch") m n nv
    iv <- V.unsafeFreeze v
    f $ fromVector mn iv
{-# INLINE withFromVector #-}


-- | View a mutable vector as a mutable matrix of the given shape and pass it
-- to the specified function.
withFromVectorM :: (Storable e)
                     => (Int,Int)
                     -> STVector s e
                     -> (STMatrix s e -> ST s a)
                     -> ST s a
withFromVectorM mn@(m,n) v f = do
    nv <- V.getDim v
    when (nv /= m*n) $ error $
        printf ("withFromVectorM (%d,%d) <vector with dim %d>:"
                ++ " dimension mismatch") m n nv
    withFromVector mn v $ \a -> do
        ma <- unsafeThaw a
        f ma
{-# INLINE withFromVectorM #-}


-- | View a vector as a matrix with one column and pass it to
-- the specified function.
withFromCol :: (RVector v, Storable e)
                => v e
                -> (forall m . RMatrix m => m e -> ST s a)
                -> ST s a
withFromCol v f = do
    m <- V.getDim v
    withFromVector (m,1) v f
{-# INLINE withFromCol #-}


-- | View a mutable vector as a mutable matrix with one column and pass it to
-- the specified function.
withFromColM :: (Storable e)
                  => STVector s e
                  -> (STMatrix s e -> ST s a)
                  -> ST s a
withFromColM v f = do
    m <- V.getDim v
    withFromVectorM (m, 1) v f
{-# INLINE withFromColM #-}


-- | View a vector as a matrix with one row and pass it to
-- the specified function.
withFromRow :: (RVector v, Storable e)
                => v e
                -> (forall m . RMatrix m => m e -> ST s a)
                -> ST s a
withFromRow v f = do
    n <- V.getDim v
    withFromVector (1,n) v f
{-# INLINE withFromRow #-}

-- | View a mutable vector as a mutable matrix with one row and pass it to
-- the specified function.
withFromRowM :: (Storable e)
                  => STVector s e
                  -> (STMatrix s e -> ST s a)
                  -> ST s a
withFromRowM v f = do
    n <- V.getDim v
    withFromVectorM (1,n) v f
{-# INLINE withFromRowM #-}

-- | Perform an action with a view of a matrix column.
withCol :: (RMatrix m, Storable e)
            => m e
            -> Int
            -> (forall v . RVector v => v e -> ST s a)
            -> ST s a
withCol a j f = do
    (m,n) <- getDim a
    when (j < 0 || j >= n) $ error $
        printf ("withCol <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n j

    unsafeWithCol a j f
{-# INLINE withCol #-}

-- | Like 'withCol', but perform the action with a mutable view.
withColM :: (Storable e)
         => STMatrix s e
         -> Int
         -> (STVector s e -> ST s a)
         -> ST s a
withColM a j f = do
    (m,n) <- getDim a
    when (j < 0 || j >= n) $ error $
        printf ("withColM <matrix with dim (%d,%d)> %d:"
                ++ " index out of range") m n j

    unsafeWithColM a j f
{-# INLINE withColM #-}

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

-- | Same as 'copyTo' but does not range-check indices.
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
getElems a = case maybeWithVector a V.getElems of
    Just es -> es
    Nothing -> withCols a $ \xs ->
                   concat `fmap` mapM V.getElems xs

-- | Get the elements of the matrix, in column-major order.
getElems' :: (RMatrix m, Storable e) => m e -> ST s [e]
getElems' a = case maybeWithVector a V.getElems' of
    Just es -> es
    Nothing -> withCols a $ \xs ->
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
    case maybeWithVectorM a (`V.setElems` es) of
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
            withColM a j (`V.setElems` es1')
            go m n (j+1) es2'

-- | Set the given values in the matrix.  If an index is repeated twice,
-- the value is implementation-defined.
setAssocs :: (Storable e) => STMatrix s e -> [((Int,Int),e)] -> ST s ()
setAssocs a ies =
    sequence_ [ write a i e | (i,e) <- ies ]

-- | Same as 'setAssocs' but does not range-check indices.
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

-- | Same as 'setRow' but does not range-check index or check
-- vector dimension.
unsafeSetRow :: (RVector v, Storable e)
             => STMatrix s e -> Int -> v e -> ST s ()
unsafeSetRow a i x = do
    jes <- V.getAssocs x
    sequence_ [ unsafeWrite a (i,j) e | (j,e) <- jes ]
{-# INLINE unsafeSetRow #-}

-- | Exchange corresponding elements in the given rows.
swapRows :: (BLAS1 e)
         => STMatrix s e -> Int -> Int -> ST s ()
swapRows a i1 i2 = do
    (m,n) <- getDim a
    when (i1 < 0 || i1 >= m || i2 < 0 || i2 >= m) $ error $
        printf ("swapRows <matrix with dim (%d,%d)> %d %d"
                ++ ": index out of range") m n i1 i2
    unsafeSwapRows a i1 i2

-- | Same as 'swapRows' but does not range-check indices.
unsafeSwapRows :: (BLAS1 e)
               => STMatrix s e -> Int -> Int -> ST s ()
unsafeSwapRows a i1 i2 = when (i1 /= i2) $ do
    (_,n) <- getDim a
    unsafeIOToST $
        unsafeWith a $ \pa lda ->
            let px = pa `advancePtr` i1
                py = pa `advancePtr` i2
                incx = lda
                incy = lda
            in
                BLAS.swap n px incx py incy

-- | Exchange corresponding elements in the given columns.
swapCols :: (BLAS1 e)
         => STMatrix s e -> Int -> Int -> ST s ()
swapCols a j1 j2 = do
    (m,n) <- getDim a
    when (j1 < 0 || j1 >= n || j2 < 0 || j2 >= n) $ error $
        printf ("swapCols <matrix with dim (%d,%d)> %d %d"
                ++ ": index out of range") m n j1 j2
    unsafeSwapCols a j1 j2

-- | Same as 'swapCols' but does not range-check indices.
unsafeSwapCols :: (BLAS1 e)
               => STMatrix s e -> Int -> Int -> ST s ()
unsafeSwapCols a j1 j2 = when (j1 /= j2) $ do
    (m,_) <- getDim a
    unsafeIOToST $
        unsafeWith a $ \pa lda ->
            let px = pa `advancePtr` (j1*lda)
                py = pa `advancePtr` (j2*lda)
                incx = 1
                incy = 1
            in
                BLAS.swap m px incx py incy

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

-- | Same as 'rowTo' but does not range-check index or check dimension.
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

-- | Same as 'setDiag' but does not range-check index or check dimension.
unsafeSetDiag :: (RVector v, Storable e)
              => STMatrix s e -> v e -> ST s ()
unsafeSetDiag a x = do
    ies <- V.getAssocs x
    sequence_ [ unsafeWrite a (i,i) e | (i,e) <- ies ]
{-# INLINE unsafeSetDiag #-}

-- | Copy the diagonal of the matrix to the vector.
diagTo :: (RMatrix m, Storable e)
       => STVector s e -> m e -> ST s ()
diagTo x a = do
    nx <- V.getDim x
    (m,n) <- getDim a
    let mn = min m n
    
    when (nx /= mn) $ error $
        printf ("diagTo"
                ++ " <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch") nx m n

    unsafeDiagTo x a
{-# INLINE diagTo #-}

-- | Same as 'diagTo' but does not range-check index or check dimensions.
unsafeDiagTo :: (RMatrix m, Storable e)
             => STVector s e -> m e -> ST s ()
unsafeDiagTo x a = do
    (m,n) <- getDim a
    let mn = min m n
    forM_ [ 0..mn-1 ] $ \i -> do
        e <- unsafeRead a (i,i)
        V.unsafeWrite x i e
{-# INLINE unsafeDiagTo #-}

-- | Get the element stored at the given index.
read :: (RMatrix m, Storable e) => m e -> (Int,Int) -> ST s e
read a (i,j) = do
    (m,n) <- getDim a
    when (i < 0 || i >= m || j < 0 || j >= n) $ error $
        printf ("read <matrix with dim (%d,%d)> (%d,%d):"
                ++ " index out of range") m n i j
    unsafeRead a (i,j)
{-# INLINE read #-}

-- | Same as 'read' but does not range-check index.
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

-- | Same as 'write' but does not range-check index.
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

-- | Same as 'modify' but does not range-check index.
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
             
-- | Same as 'mapTo' but does not check dimensions.
unsafeMapTo :: (RMatrix m, Storable e, Storable f)
            => STMatrix s f
            -> (e -> f)
            -> m e
            -> ST s ()
unsafeMapTo dst f src =
    fromMaybe colwise $ maybeWithVectorM dst $ \vdst ->
    fromMaybe colwise $ maybeWithVector  src $ \vsrc ->
        V.unsafeMapTo vdst f vsrc
  where
    colwise = withColsM dst $ \zs ->
              withCols   src $ \xs ->
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

-- | Same as 'zipWithTo' but does not check dimensions.
unsafeZipWithTo :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
                => STMatrix s f
                -> (e1 -> e2 -> f)
                -> m1 e1
                -> m2 e2
                -> ST s ()
unsafeZipWithTo dst f x y =
    fromMaybe colwise $ maybeWithVectorM dst $ \vdst ->
    fromMaybe colwise $ maybeWithVector    x $ \vx ->
    fromMaybe colwise $ maybeWithVector    y $ \vy ->
        V.unsafeZipWithTo vdst f vx vy
  where
    colwise = withColsM dst $ \vdsts ->
              withCols   x $ \vxs ->
              withCols   y $ \vys ->
              
                  sequence_ [ V.unsafeZipWithTo vdst f vx vy
                            | (vdst,vx,vy) <- zip3 vdsts vxs vys
                            ]

-- | Set every element in the matrix to a default value.  For
-- standard numeric types (including 'Double', 'Complex Double', and 'Int'),
-- the default value is '0'.
clear :: (Storable e) => STMatrix s e -> ST s ()
clear a = fromMaybe colwise $ maybeWithVectorM a V.clear
  where
    colwise = withColsM a $ mapM_ V.clear

-- | @withSlice (i,j) (m,n) a@ performs an action with a view of the
-- submatrix of @a@ starting at index @(i,j)@ and having dimension @(m,n)@.
withSlice :: (RMatrix m, Storable e)
          => (Int,Int)
          -> (Int,Int)
          -> m e
          -> (forall m'. RMatrix m' => m' e -> ST s a)
          -> ST s a
withSlice ij mn a f = do
    ia <- unsafeFreeze a
    f $ slice ij mn ia

-- | Like 'withSlice', but perform the action with a mutable view.
withSliceM :: (Storable e)
           => (Int,Int)
           -> (Int,Int)
           -> STMatrix s e
           -> (STMatrix s e -> ST s a)
           -> ST s a
withSliceM ij mn a f =
    withSlice ij mn a $ \a' -> do
        ma <- unsafeThaw a'
        f ma

-- | Perform an action with a view gotten from taking the given number of
-- rows from the start of the matrix.
withTakeRows :: (RMatrix m, Storable e)
             => Int
             -> m e
             -> (forall m'. RMatrix m' => m' e -> ST s a)
             -> ST s a
withTakeRows i a f = do
    ia <- unsafeFreeze a
    f $ takeRows i ia

-- | Like 'withTakeRows', but perform the action with a mutable view.
withTakeRowsM :: (Storable e)
              => Int
              -> STMatrix s e
              -> (STMatrix s e -> ST s a)
              -> ST s a
withTakeRowsM i a f =
    withTakeRows i a $ \a' -> do
        ma <- unsafeThaw a'
        f ma

-- | Perform an action with a view gotten from dropping the given number of
-- rows from the start of the matrix.
withDropRows :: (RMatrix m, Storable e)
             => Int
             -> m e
             -> (forall m'. RMatrix m' => m' e -> ST s a)
             -> ST s a
withDropRows n a f = do
    ia <- unsafeFreeze a
    f $ dropRows n ia

-- | Like 'withDropRows', but perform the action with a mutable view.
withDropRowsM :: (Storable e)
              => Int
              -> STMatrix s e
              -> (STMatrix s e -> ST s a)
              -> ST s a
withDropRowsM i a f =
    withDropRows i a $ \a' -> do
        ma <- unsafeThaw a'
        f ma

-- | Perform an action with views from splitting the matrix rows at the given
-- index.
withSplitRowsAt :: (RMatrix m, Storable e)
                => Int
                -> m e
                -> (forall m1 m2. (RMatrix m1, RMatrix m2) => m1 e -> m2 e -> ST s a)
                -> ST s a
withSplitRowsAt i a f = do
    ia <- unsafeFreeze a
    uncurry f $ splitRowsAt i ia

-- | Like 'withSplitRowsAt', but perform the action with a mutable view.
withSplitRowsAtM :: (Storable e)
                 => Int
                 -> STMatrix s e
                 -> (STMatrix s e -> STMatrix s e -> ST s a)
                 -> ST s a
withSplitRowsAtM i a f =
    withSplitRowsAt i a $ \a1' a2' -> do
        ma1 <- unsafeThaw a1'
        ma2 <- unsafeThaw a2'        
        f ma1 ma2
    
-- | Perform an action with a view gotten from taking the given number of
-- columns from the start of the matrix.
withTakeCols :: (RMatrix m, Storable e)
             => Int
             -> m e
             -> (forall m'. RMatrix m' => m' e -> ST s a)
             -> ST s a
withTakeCols i a f = do
    ia <- unsafeFreeze a
    f $ takeCols i ia

-- | Like 'withTakeCols', but perform the action with a mutable view.
withTakeColsM :: (Storable e)
              => Int
              -> STMatrix s e
              -> (STMatrix s e -> ST s a)
              -> ST s a
withTakeColsM i a f =
    withTakeCols i a $ \a' -> do
        ma <- unsafeThaw a'
        f ma

-- | Perform an action with a view gotten from dropping the given number of
-- columns from the start of the matrix.
withDropCols :: (RMatrix m, Storable e)
             => Int
             -> m e
             -> (forall m'. RMatrix m' => m' e -> ST s a)
             -> ST s a
withDropCols n a f = do
    ia <- unsafeFreeze a
    f $ dropCols n ia

-- | Like 'withDropCols', but perform the action with a mutable view.
withDropColsM :: (Storable e)
              => Int
              -> STMatrix s e
              -> (STMatrix s e -> ST s a)
              -> ST s a
withDropColsM i a f =
    withDropCols i a $ \a' -> do
        ma <- unsafeThaw a'
        f ma

-- | Perform an action with views from splitting the matrix columns at the given
-- index.
withSplitColsAt :: (RMatrix m, Storable e)
                => Int
                -> m e
                -> (forall m1 m2. (RMatrix m1, RMatrix m2) => m1 e -> m2 e -> ST s a)
                -> ST s a
withSplitColsAt i a f = do
    ia <- unsafeFreeze a
    uncurry f $ splitColsAt i ia

-- | Like 'withSplitColsAt', but perform the action with mutable views.    
withSplitColsAtM :: (Storable e)
                 => Int
                 -> STMatrix s e
                 -> (STMatrix s e -> STMatrix s e -> ST s a)
                 -> ST s a
withSplitColsAtM i a f =
    withSplitColsAt i a $ \a1' a2' -> do
        ma1 <- unsafeThaw a1'
        ma2 <- unsafeThaw a2'        
        f ma1 ma2


-- | Add a vector to the diagonal of a matrix.
shiftDiagM_ :: (RVector v, BLAS1 e)
              => v e -> STMatrix s e -> ST s ()
shiftDiagM_ s a = do
    (m,n) <- getDim a
    ns <- V.getDim s
    let mn = min m n
    
    when (ns /= mn) $ error $
        printf ("shiftDiagM_"
                ++ " <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)>"
                ++ ": dimension mismatch")
                ns
                m n
    
    shiftDiagWithScaleM_ 1 s a
    

-- | Add a scaled vector to the diagonal of a matrix.
shiftDiagWithScaleM_ :: (RVector v, BLAS1 e)
                       => e -> v e -> STMatrix s e -> ST s ()
shiftDiagWithScaleM_ e s a = do
    (m,n) <- getDim a
    ns <- V.getDim s
    let mn = min m n

    when (ns /= mn) $ error $
        printf ("shiftDiagWithScaleM_"
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
scaleM_ :: (BLAS1 e)
          => e -> STMatrix s e -> ST s ()
scaleM_ e = vectorOp (V.scaleM_ e)

-- | @addWithScaleM_ alpha x y@ sets @y := alpha * x + y@.
addWithScaleM_ :: (RMatrix m, BLAS1 e)
               => e -> m e -> STMatrix s e -> ST s ()
addWithScaleM_ e = checkOp2 "addWithScaleM_" $
    unsafeAddWithScaleM_ e

unsafeAddWithScaleM_ :: (RMatrix m, BLAS1 e)
                     => e -> m e -> STMatrix s e -> ST s ()
unsafeAddWithScaleM_ alpha x y =
    fromMaybe colwise $ maybeWithVector  x $ \vx ->
    fromMaybe colwise $ maybeWithVectorM y $ \vy ->
        V.unsafeAddWithScaleM_ alpha vx vy
  where
    colwise = withCols   x $ \vxs ->
              withColsM y $ \vys ->
                  sequence_ [ V.unsafeAddWithScaleM_ alpha vx vy
                            | (vx,vy) <- zip vxs vys ]                

-- | Scale the rows of a matrix; @scaleRowsM_ s a@ sets
-- @a := diag(s) * a@.
scaleRowsM_ :: (RVector v, BLAS1 e)
              => v e -> STMatrix s e -> ST s ()
scaleRowsM_  s a = do
    (m,n) <- getDim a
    ns <- V.getDim s
    when (ns /= m) $ error $
        printf ("scaleRowsM_"
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
scaleColsM_ :: (RVector v, BLAS1 e)
            => v e -> STMatrix s e -> ST s ()
scaleColsM_ s a = do
    (m,n) <- getDim a
    ns <- V.getDim s
    when (ns /= n) $ error $
        printf ("scaleColsM_"
                ++ " <vector with dim %d>"
                ++ " <matrix with dim (%d,%d)>"        
                ++ ": dimension mismatch") 
                ns
                m n

    es <- V.getElems s
    withColsM a $ \xs ->
        sequence_ [ V.scaleM_ e x
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
    fromMaybe colwise $ maybeWithVectorM x $ \vx -> f vx
  where
    colwise = withColsM x $ \vxs ->
                  sequence_ [ f vx | vx <- vxs ]

vectorOp2 :: (RMatrix m, Storable e, Storable f)
          => (forall v . RVector v => STVector s f -> v e -> ST s ())
          -> STMatrix s f -> m e -> ST s ()
vectorOp2 f dst x =
    fromMaybe colwise $ maybeWithVectorM dst $ \vdst ->
    fromMaybe colwise $ maybeWithVector    x $ \vx ->
        f vdst vx
  where
    colwise = withColsM dst $ \vdsts ->
              withCols   x   $ \vxs ->
                  sequence_ [ f vdst vx | (vdst,vx) <- zip vdsts vxs ]
{-# INLINE vectorOp2 #-}

vectorOp3 :: (RMatrix m1, RMatrix m2, Storable e1, Storable e2, Storable f)
          => (forall v1 v2 . (RVector v1, RVector v2) => 
                  STVector s f -> v1 e1 -> v2 e2 -> ST s ())
          -> STMatrix s f -> m1 e1 -> m2 e2 -> ST s ()
vectorOp3 f dst x y =
    fromMaybe colwise $ maybeWithVectorM dst $ \vdst ->
    fromMaybe colwise $ maybeWithVector    x $ \vx ->
    fromMaybe colwise $ maybeWithVector    y $ \vy ->
        f vdst vx vy
  where
    colwise = withColsM dst $ \vdsts ->
              withCols   x   $ \vxs ->
              withCols   y   $ \vys ->
                  sequence_ [ f vdst vx vy
                            | (vdst,vx,vy) <- zip3 vdsts vxs vys ]
{-# INLINE vectorOp3 #-}
