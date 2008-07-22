{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Dense.Internal
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Dense.Internal (
    -- * Dense matrix data types
    DMatrix(..),
    IOMatrix,
    Matrix,

    module BLAS.Matrix.Base,
    module BLAS.Tensor,

    -- * Converting to and from foreign pointers
    toForeignPtr,
    fromForeignPtr,
    ldaOf,
    isHerm,

    -- * Creating new matrices
    newMatrix,
    newMatrix_,
    newListMatrix,
    newColsMatrix,
    newRowsMatrix,
    
    listMatrix,

    -- * Special matrices
    newIdentity,
    setIdentity,
    
    -- * Row and column views
    row,
    col,
    rows,
    cols,
    
    -- * Diagonal views
    diag,
    
    -- * Matrix views
    submatrix,

    -- * Converting to/from vectors
    maybeFromRow,
    maybeFromCol,
    maybeToVector,
    
    -- * Lifting scalar and vector operations
    liftV,
    liftV2,
    
    -- * Casting matrices
    coerceMatrix,
    
    -- * Unsafe operations
    unsafeThaw,
    unsafeFreeze,
    unsafeNewMatrix,
    unsafeWithElemPtr,
    unsafeRow,
    unsafeCol,
    unsafeDiag,
    unsafeSubmatrix,
    ) where

import Control.Monad ( forM_, zipWithM_ )
import Data.Ix ( inRange, range )
import Foreign
import System.IO.Unsafe   
import Unsafe.Coerce

import Data.AEq

import Data.Vector.Dense.Internal hiding ( toForeignPtr, fromForeignPtr,
    unsafeFreeze, unsafeThaw, fptr, offset, unsafeWithElemPtr )
import qualified Data.Vector.Dense.Internal as V
import qualified Data.Vector.Dense.Operations as V

import BLAS.Access
import BLAS.Internal ( inlinePerformIO, checkedRow, checkedCol, checkedDiag,
    checkedSubmatrix, diagStart, diagLen )
import BLAS.Elem ( Elem, BLAS1 )
import qualified BLAS.Elem as E
import BLAS.Matrix.Base hiding ( Matrix )
import qualified BLAS.Matrix.Base as C
import BLAS.Tensor
import BLAS.Types


-- | The mutable dense matrix data type.  It can either store elements in 
-- column-major order, or provide a view into another matrix.  The view 
-- transposes and conjugates the underlying matrix.
data DMatrix t mn e =
      DM { fptr   :: !(ForeignPtr e) -- ^ a pointer to the storage region
         , offset :: !Int            -- ^ an offset (in elements, not bytes) to the first element in the matrix. 
         , size1  :: !Int            -- ^ the number of rows in the matrix
         , size2  :: !Int            -- ^ the number of columns in the matrix
         , lda    :: !Int            -- ^ the leading dimension size of the matrix
         }
    | H !(DMatrix t mn e)           -- ^ a transposed and conjugated matrix

type Matrix = DMatrix Imm
type IOMatrix = DMatrix Mut

unsafeFreeze :: DMatrix t mn e -> Matrix mn e
unsafeFreeze = unsafeCoerce

unsafeThaw :: DMatrix t mn e -> IOMatrix mn e
unsafeThaw = unsafeCoerce

-- | Coerce the phantom shape type from one type to another.
coerceMatrix :: DMatrix t mn e -> DMatrix t kl e
coerceMatrix = unsafeCoerce

-- | @fromForeignPtr f o mn l@ creates a matrix view of the data pointed to
-- by @f@ starting at offset @o@ and having shape @mn@ and lda @l@.
fromForeignPtr :: ForeignPtr e -> Int -> (Int,Int) -> Int -> DMatrix t (m,n) e
fromForeignPtr f o (m,n) l = DM f o m n l

-- | Convert a dense matrix to a pointer, offset, size, and lda.  Note that this
-- does not give the conjugacy/transpose information.  For that, use 'isHerm'.
toForeignPtr :: DMatrix t (m,n) e -> (ForeignPtr e, Int, (Int,Int), Int)
toForeignPtr   (H a)            = toForeignPtr a
toForeignPtr a@(DM _ _ _ _ _)   = (fptr a, offset a, (size1 a, size2 a), lda a)

-- | Get the lda of a matrix, defined as the number of elements in the underlying
-- array that separate two consecutive elements in the same row of the matrix.
ldaOf :: DMatrix t (m,n) e -> Int
ldaOf   (H a)          = ldaOf a
ldaOf a@(DM _ _ _ _ _) = lda a


indexOf :: DMatrix t (m,n) e -> (Int,Int) -> Int
indexOf   (H a)          (i,j) = indexOf a (j,i)
indexOf a@(DM _ _ _ _ _) (i,j) = 
    let o = offset a
        l = lda a
    in o + i + j*l
    
-- | Get the storage order of the matrix.  If 'isTrans' is true, this 
-- will be 'RowMajor'.  Otherwise, it will be 'ColMajor'.
orderOf :: DMatrix t (m,n) e -> Order
orderOf (H a)          = flipOrder (orderOf a)
orderOf (DM _ _ _ _ _) = ColMajor

-- | Get whether or not the matrix is transposed and conjugated.
isHerm :: DMatrix t (m,n) e -> Bool
isHerm (H a)          = not (isHerm a)
isHerm (DM _ _ _ _ _) = False

-- | Create a new matrix of the given size and initialize the given elements to
-- the given values.  All other elements get initialized to zero.
newMatrix :: (BLAS1 e) => (Int,Int) -> [((Int,Int), e)] -> IO (DMatrix t (m,n) e)
newMatrix = newMatrixHelp writeElem

-- | Same as 'newMatrix' but do not do any bounds-checking.
unsafeNewMatrix :: (BLAS1 e) => (Int,Int) -> [((Int,Int), e)] -> IO (DMatrix t (m,n) e)
unsafeNewMatrix = newMatrixHelp unsafeWriteElem

newMatrixHelp :: (BLAS1 e) => 
       (IOMatrix (m,n) e -> (Int,Int) -> e -> IO ()) 
    -> (Int,Int) -> [((Int,Int),e)] -> IO (DMatrix t (m,n) e)
newMatrixHelp set mn ijes = do
    x <- newZero mn
    io <- unsafeInterleaveIO $ mapM_ (uncurry $ set $ unsafeThaw x) ijes
    return $ io `seq` x

-- | Create a new matrix of given shape, but do not initialize the elements.
newMatrix_ :: (Elem e) => (Int,Int) -> IO (DMatrix t (m,n) e)
newMatrix_ (m,n) 
    | m < 0 || n < 0 =
        ioError $ userError $ 
            "Tried to create a matrix with shape `" ++ show (m,n) ++ "'"
    | otherwise = do
        f <- mallocForeignPtrArray (m*n)
        return $ fromForeignPtr f 0 (m,n) (max 1 m)

-- | Create a new matrix with the given elements in column-major order.
newListMatrix :: (Elem e) => (Int,Int) -> [e] -> IO (DMatrix t (m,n) e)
newListMatrix (m,n) es = do
    a <- newMatrix_ (m,n)
    unsafeWithElemPtr a (0,0) $ flip pokeArray (take (m*n) es)
    return a

-- | Create a new matrix with the given elements in row-major order.
listMatrix :: (Elem e) => (Int,Int) -> [e] -> Matrix (m,n) e
listMatrix mn es = unsafePerformIO $ newListMatrix mn es
{-# NOINLINE listMatrix #-}

-- | Create a new matrix of the given shape with ones along the diagonal, 
-- and zeros everywhere else.
newIdentity :: (BLAS1 e) => (Int,Int) -> IO (DMatrix t (m,n) e)
newIdentity mn = do
    a <- newMatrix_ mn
    setIdentity (unsafeThaw a)
    return a

-- | Set the diagonal to ones, and set everywhere else to zero.
setIdentity :: (BLAS1 e) => IOMatrix (m,n) e -> IO ()
setIdentity a = do
    s <- getSize a
    case s of
        0 -> return ()
        _ -> setZero a >>
             setConstant 1 (diag a 0)

-- | Form a matrix from a list of column vectors.
newColsMatrix :: (BLAS1 e) => (Int,Int) -> [DVector t m e] -> IO (DMatrix r (m,n) e)
newColsMatrix (m,n) cs = do
    a <- newZero (m,n)
    forM_ (zip [0..(n-1)] cs) $ \(j,c) ->
        V.copyVector (unsafeCol (unsafeThaw a) j) c
    return a

-- | Form a matrix from a list of row vectors.
newRowsMatrix :: (BLAS1 e) => (Int,Int) -> [DVector t n e] -> IO (DMatrix r (m,n) e)
newRowsMatrix (m,n) rs = do
    a <- newZero (m,n)
    forM_ (zip [0..(m-1)] rs) $ \(i,r) ->
        V.copyVector (unsafeRow (unsafeThaw a) i) r
    return a

-- | Evaluate a function with a pointer to the raw storage for the element
-- at the given index.  It may be necessary to conjugate or scale values before
-- reading or writing to or from the location.
unsafeWithElemPtr :: (Elem e) => DMatrix t (m,n) e -> (Int,Int) -> (Ptr e -> IO a) -> IO a
unsafeWithElemPtr   (H a) (i,j) f = unsafeWithElemPtr a (j,i) f
unsafeWithElemPtr a@(DM _ _ _ _ _) ij f =
    withForeignPtr (fptr a) $ \ptr ->
        let ptr' = ptr `advancePtr` (indexOf a ij)
        in f ptr'
        
-- | Get a vector view of the given row in a matrix.
row :: (Elem e) => DMatrix t (m,n) e -> Int -> DVector t n e
row a = checkedRow (shape a) (unsafeRow a)

-- | Get a list of vector views of the rows of the matrix.
rows :: (Elem e) => DMatrix t (m,n) e -> [DVector t n e]
rows a = [ unsafeRow a i | i <- [0..numRows a - 1] ]

-- | Get a list of vector views of the columns of the matrix.
cols :: (Elem e) => DMatrix t (m,n) e -> [DVector t m e]
cols a = [ unsafeCol a j | j <- [0..numCols a - 1] ]

-- | Get a vector view of the given column in a matrix.
col :: (Elem e) => DMatrix t (m,n) e -> Int -> DVector t m e
col a = checkedCol (shape a) (unsafeCol a)

-- | Same as 'row', but does not do any bounds checking.
unsafeRow ::  (Elem e) => DMatrix t (m,n) e -> Int -> DVector t n e
unsafeRow a@(H _)          i = conj $ unsafeCol (herm a) i
unsafeRow a@(DM _ _ _ _ _) i =
    let f = fptr a
        o = indexOf a (i,0)
        n = size2 a
        s = lda a
    in V.fromForeignPtr f o n s

-- | Same as 'col', but does not do any bounds checking.
unsafeCol :: (Elem e) => DMatrix t (m,n) e -> Int -> DVector t m e
unsafeCol a@(H _)          j = conj $ unsafeRow (herm a) j
unsafeCol a@(DM _ _ _ _ _) j =
    let f = fptr a
        o = indexOf a (0,j)
        m = size1 a
        s = 1
    in V.fromForeignPtr f o m s

-- | @diag a 0@ gets a vector view of the main diagonal of @a@.  @diag a k@ for 
-- @k@ positive gets a view of the @k@th superdiagonal.  For @k@ negative, it
-- gets a view of the @(-k)@th subdiagonal.
diag :: (Elem e) => DMatrix t (m,n) e -> Int -> DVector t k e
diag a = checkedDiag (shape a) (unsafeDiag a)

-- | Same as 'diag', but does not do any bounds checking.
unsafeDiag :: (Elem e) => DMatrix t (m,n) e -> Int -> DVector t k e
unsafeDiag   (H a)          i = conj $ unsafeDiag a (negate i)
unsafeDiag a@(DM _ _ _ _ _) i =
    let f = fptr a
        o = indexOf a (diagStart i)
        n = diagLen (shape a) i
        s = lda a + 1
    in V.fromForeignPtr f o n s

-- | @submatrix a ij mn@ returns a view of the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrix :: (Elem e) => DMatrix t (m,n) e -> (Int,Int) -> (Int,Int) -> DMatrix t (k,l) e
submatrix a = checkedSubmatrix (shape a) (unsafeSubmatrix a)

-- | Same as 'submatrix' but does not do any bounds checking.
unsafeSubmatrix :: (Elem e) => DMatrix t (m,n) e -> (Int,Int) -> (Int,Int) -> DMatrix t (k,l) e
unsafeSubmatrix a@(H _)          (i,j) (m',n') = herm $ unsafeSubmatrix (herm a) (j,i) (n',m')
unsafeSubmatrix a@(DM _ _ _ _ _) (i,j) mn' =
    let f = fptr a
        o = indexOf a (i,j)
        l = lda a
    in fromForeignPtr f o mn' l
    

-- | Create a matrix view of a row vector.  This will fail if the
-- stride is not @1@ and the vector is conjugated.
maybeFromRow :: (Elem e) => DVector t m e -> Maybe (DMatrix t (one,m) e)
maybeFromRow (V.C   (V.C x))   = maybeFromRow x
maybeFromRow (V.C x@(V.DV _ _ _ _))
    | V.stride x == 1 =
        let f = V.fptr x
            o = V.offset x
            n = V.dim x
            l = max 1 n
        in Just $ herm $ fromForeignPtr f o (n,1) l
    | otherwise =
        Nothing
maybeFromRow x@(V.DV _ _ _ _) =
    let f = V.fptr x
        o = V.offset x
        n = V.dim x
        s = V.stride x
        l = max 1 s
    in Just $ fromForeignPtr f o (1,n) l


-- | Possibly create a matrix view of a column vector.  This will fail
-- if the stride of the vector is not @1@ and the vector is not conjugated.
maybeFromCol :: (Elem e) => DVector t n e -> Maybe (DMatrix t (n,one) e)
maybeFromCol   (V.C x)   = maybeFromRow x >>= return . herm
maybeFromCol x@(V.DV _ _ _ _)
    | V.stride x == 1 =
        let f = V.fptr x
            o = V.offset x
            m = dim x
            l = max 1 m
        in Just $ fromForeignPtr f o (m,1) l
    | otherwise =
        Nothing

maybeToVector :: (Elem e) => DMatrix t (m,n) e -> Maybe (Order, DVector t k e)
maybeToVector (H a)   = maybeToVector a >>= (\(o,x) -> return (flipOrder o, conj x))
maybeToVector (DM f o m n ld)
    | ld == m =
        Just $ (ColMajor, V.fromForeignPtr f o (m*n) 1)
    | m == 1 =
        Just $ (ColMajor, V.fromForeignPtr f o n    ld)
    | otherwise =
        Nothing

-- | Modify each element in-place by applying a function to it.
-- modifyWith :: (Elem e) => (e -> e) -> IOMatrix (m,n) e -> IO ()


-- | Take a unary elementwise vector operation and apply it to the elements of a matrix.
liftV :: (Elem e) => (DVector t k e -> IO ()) -> DMatrix t (m,n) e -> IO ()
liftV f a =
    case maybeToVector a of
        Just (_,x) -> f x
        _ -> 
            let xs  = case orderOf a of
                          RowMajor -> rows (coerceMatrix a)
                          ColMajor -> cols (coerceMatrix a)
            in mapM_ f xs

-- | Take a binary elementwise vector operation and apply it to the elements of a pair
-- of matrices.
liftV2 :: (Elem e) => (DVector s k e -> DVector t k e -> IO ()) 
       -> DMatrix s (m,n) e -> DMatrix t (m,n) e -> IO ()
liftV2 f a b =
    case (maybeToVector a, maybeToVector b) of
        (Just (RowMajor,x), Just (RowMajor,y)) -> f x y
        (Just (ColMajor,x), Just (ColMajor,y)) -> f x y
        _ -> 
            let (xs,ys) = case orderOf a of
                               RowMajor -> (rows (coerceMatrix a), rows (coerceMatrix b))
                               ColMajor -> (cols (coerceMatrix a), cols (coerceMatrix b))
            in zipWithM_ f xs ys


instance C.Matrix (DMatrix t) where
    numRows = fst . shape
    numCols = snd . shape

    herm a = case a of
        (H a')   -> coerceMatrix a'
        _        -> H (coerceMatrix a)
    
instance Tensor (DMatrix t (m,n)) (Int,Int) e where
    shape a = case a of
        (H a')   -> case shape a' of (m,n) -> (n,m)
        _        -> (size1 a, size2 a)
    
    bounds a = let (m,n) = shape a in ((0,0), (m-1,n-1))

instance (BLAS1 e) => ITensor (DMatrix Imm (m,n)) (Int,Int) e where
    size a = (numRows a * numCols a)
    
    unsafeAt a = inlinePerformIO . unsafeReadElem a
    {-# INLINE unsafeAt #-}
    
    indices a = [ (i,j) | j <- range (0,n-1), i <- range (0,m-1) ]
      where (m,n) = shape a
    {-# INLINE indices #-}
    
    elems = inlinePerformIO . getElems
    {-# INLINE elems #-}
    
    assocs = inlinePerformIO . getAssocs
    {-# INLINE assocs #-}

    (//) = replaceHelp writeElem
    unsafeReplace = replaceHelp unsafeWriteElem

    amap f a = listMatrix (shape a) (map f $ elems a)


replaceHelp :: (BLAS1 e) => 
       (IOMatrix (m,n) e -> (Int,Int) -> e -> IO ())
    -> Matrix (m,n) e -> [((Int,Int), e)] -> Matrix (m,n) e
replaceHelp set x ies = 
    unsafeFreeze $ unsafePerformIO $ do
        y  <- newCopy (unsafeThaw x)
        mapM_ (uncurry $ set y) ies
        return y
{-# NOINLINE replaceHelp #-}


instance (BLAS1 e) => IDTensor (DMatrix Imm (m,n)) (Int,Int) e where
    zero = unsafePerformIO . newZero
    {-# NOINLINE zero #-}
     
    constant mn = unsafePerformIO . newConstant mn
    {-# NOINLINE constant #-}
     
    azipWith f a b
        | shape b /= mn =
            error ("azipWith: matrix shapes differ; first has shape `"
                   ++ show mn ++ "' and second has shape `" 
                   ++ show (shape b) ++ "'")
        | otherwise =
            listMatrix mn (zipWith f (elems a) (elems b))
        where
            mn = shape a
     
     
instance (BLAS1 e) => RTensor (DMatrix t (m,n)) (Int,Int) e IO where
    getSize a = return (numRows a * numCols a)
    
    newCopy a = case a of
        (H a')   -> newCopy a' >>= return . H
        _        -> do
            a' <- newMatrix_ (shape a)
            liftV2 V.copyVector (unsafeThaw a') a
            return a'
    
    unsafeReadElem a (i,j) = case a of
        (H a')   -> unsafeReadElem a' (j,i) >>= return . E.conj
        _        -> withForeignPtr (fptr a) $ \ptr ->
                        peekElemOff ptr (indexOf a (i,j))
    {-# INLINE unsafeReadElem #-}
    
    getIndices = return . indices . unsafeFreeze
    {-# INLINE getIndices #-}

    getElems a = return $ go (cols a)
        where go cs | cs `seq` False = undefined
              go []     = []
              go (c:cs) =
                  let e  = inlinePerformIO $ getElems c
                      es = go cs
                  in e ++ es
    {-# NOINLINE getElems #-}

    getAssocs a = return $ go (cols a) 0
        where go cs j | cs `seq` j `seq` False = undefined
              go []     _ = []
              go (c:cs) j =
                  let ie   = inlinePerformIO $ getAssocs c
                      ije  = map (\(i,e) -> ((i,j),e)) ie
                      ijes = go cs (j+1)
                  in ije ++ ijes
    {-# NOINLINE getAssocs #-}


instance (BLAS1 e) => RDTensor (DMatrix t (m,n)) (Int,Int) e IO where
    newZero mn = do
        a <- newMatrix_ mn
        setZero (unsafeThaw a)
        return a
    
    newConstant mn e = do
        a <- newMatrix_ mn
        setConstant e (unsafeThaw a)
        return a
    
    
instance (BLAS1 e) => MTensor (DMatrix Mut (m,n)) (Int,Int) e IO where 
    setZero = liftV setZero
    setConstant e = liftV (setConstant e)
    
    canModifyElem a ij =
        return $ inRange (bounds a) ij
    {-# INLINE canModifyElem #-}

    unsafeWriteElem a (i,j) e = case a of
        (H a')   -> unsafeWriteElem a' (j,i) $ E.conj e
        _        -> withForeignPtr (fptr a) $ \ptr ->
                        pokeElemOff ptr (indexOf a (i,j)) e
    
    modifyWith f = liftV (modifyWith f)       


instance (BLAS1 e, Show e) => Show (DMatrix Imm (m,n) e) where
    show a = case a of
        (H a')   -> "herm (" ++ show a' ++ ")"
        _        -> "listMatrix " ++ show (shape a) ++ " " ++ show (elems a)
        
compareHelp :: (BLAS1 e) => 
    (e -> e -> Bool) -> Matrix (m,n) e -> Matrix (m,n) e -> Bool
compareHelp cmp x y
    | isHerm x && isHerm y =
        compareHelp cmp (herm x) (herm y)
compareHelp cmp x y =
    (shape x == shape y) && (and $ zipWith cmp (elems x) (elems y))

instance (BLAS1 e, Eq e) => Eq (DMatrix Imm (m,n) e) where
    (==) = compareHelp (==)

instance (BLAS1 e, AEq e) => AEq (DMatrix Imm (m,n) e) where
    (===) = compareHelp (===)
    (~==) = compareHelp (~==)
    