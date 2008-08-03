{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}
{-# OPTIONS_GHC -fglasgow-exts #-}
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
    unsafeFreeze, unsafeThaw, storageOf, offsetOf, unsafeWithElemPtr )
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
      DM { storageOf :: {-# UNPACK #-} !(ForeignPtr e) -- ^ a pointer to the storage region
         , offsetOf  :: {-# UNPACK #-} !Int            -- ^ an offset (in elements, not bytes) to the first element in the matrix. 
         , size1     :: {-# UNPACK #-} !Int            -- ^ the number of rows in the matrix
         , size2     :: {-# UNPACK #-} !Int            -- ^ the number of columns in the matrix
         , ldaOf     :: {-# UNPACK #-} !Int            -- ^ the leading dimension size of the matrix
         , isHerm    :: {-# UNPACK #-} !Bool           -- ^ indicates whether or not the matrix is transposed and conjugated
         }

type Matrix = DMatrix Imm
type IOMatrix = DMatrix Mut

unsafeFreeze :: DMatrix t mn e -> Matrix mn e
unsafeFreeze = unsafeCoerce

unsafeThaw :: DMatrix t mn e -> IOMatrix mn e
unsafeThaw = unsafeCoerce

-- | Coerce the phantom shape type from one type to another.
coerceMatrix :: DMatrix t mn e -> DMatrix t kl e
coerceMatrix = unsafeCoerce

-- | @fromForeignPtr f o mn l h@ creates a matrix view of the data pointed to
-- by @f@ starting at offset @o@ and having shape @mn@ and lda @l@.  If @h@
-- is @True@ the matrix is interpreted as transposed and conjugated.
fromForeignPtr :: ForeignPtr e -> Int -> (Int,Int) -> Int -> Bool -> DMatrix t (m,n) e
fromForeignPtr f o (m,n) l h = DM f o m n l h

-- | Convert a dense matrix to a pointer, offset, size, lda, herm.
toForeignPtr :: DMatrix t (m,n) e -> (ForeignPtr e, Int, (Int,Int), Int, Bool)
toForeignPtr a = (storageOf a, offsetOf a, (size1 a, size2 a), ldaOf a, isHerm a)

indexOf :: DMatrix t (m,n) e -> (Int,Int) -> Int
indexOf a (i,j) = 
    let (i',j') = case isHerm a of
                        True  -> (j,i)
                        False -> (i,j)
        o = offsetOf a
        l = ldaOf a
    in o + i' + j'*l
{-# INLINE indexOf #-}
    
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
        return $ fromForeignPtr f 0 (m,n) (max 1 m) False

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
unsafeWithElemPtr a ij f =
    withForeignPtr (storageOf a) $ \ptr ->
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
unsafeRow a i
    | isHerm a =
        conj $ unsafeCol (herm a) i
    | otherwise =
        let f = storageOf a
            o = indexOf a (i,0)
            n = numCols a
            s = ldaOf a
            c = False
        in V.fromForeignPtr f o n s c

-- | Same as 'col', but does not do any bounds checking.
unsafeCol :: (Elem e) => DMatrix t (m,n) e -> Int -> DVector t m e
unsafeCol a j 
    | isHerm a =
        conj $ unsafeRow (herm a) j
    | otherwise =
        let f = storageOf a
            o = indexOf a (0,j)
            m = numRows a
            s = 1
            c = False
        in V.fromForeignPtr f o m s c

-- | @diag a 0@ gets a vector view of the main diagonal of @a@.  @diag a k@ for 
-- @k@ positive gets a view of the @k@th superdiagonal.  For @k@ negative, it
-- gets a view of the @(-k)@th subdiagonal.
diag :: (Elem e) => DMatrix t (m,n) e -> Int -> DVector t k e
diag a = checkedDiag (shape a) (unsafeDiag a)

-- | Same as 'diag', but does not do any bounds checking.
unsafeDiag :: (Elem e) => DMatrix t (m,n) e -> Int -> DVector t k e
unsafeDiag a i 
    | isHerm a = 
        conj $ unsafeDiag (herm a) (negate i)
    | otherwise =            
        let f = storageOf a
            o = indexOf a (diagStart i)
            n = diagLen (shape a) i
            s = ldaOf a + 1
            c = False
        in V.fromForeignPtr f o n s c

-- | @submatrix a ij mn@ returns a view of the submatrix of @a@ with element @(0,0)@
-- being element @ij@ in @a@, and having shape @mn@.
submatrix :: (Elem e) => DMatrix t (m,n) e -> (Int,Int) -> (Int,Int) -> DMatrix t (k,l) e
submatrix a = checkedSubmatrix (shape a) (unsafeSubmatrix a)

-- | Same as 'submatrix' but does not do any bounds checking.
unsafeSubmatrix :: (Elem e) => DMatrix t (m,n) e -> (Int,Int) -> (Int,Int) -> DMatrix t (k,l) e
unsafeSubmatrix a (i,j) (m,n)
    | isHerm a  = 
        herm $ unsafeSubmatrix (herm a) (j,i) (n,m)
    | otherwise =
        let f = storageOf a
            o = indexOf a (i,j)
            l = ldaOf a
        in fromForeignPtr f o (m,n) l False
    

-- | Create a matrix view of a row vector.  This will fail if the
-- vector is conjugated and the stride is not @1@.
maybeFromRow :: (Elem e) => DVector t m e -> Maybe (DMatrix t (one,m) e)
maybeFromRow x
    | isConj x && strideOf x == 1 =
        let f = V.storageOf x
            o = V.offsetOf x
            n = dim x
            l = max 1 n
            h = True
        in Just $ fromForeignPtr f o (n,1) l h
    | (not . isConj) x =
        let f = V.storageOf x
            o = V.offsetOf x
            n = dim x
            s = strideOf x
            l = max 1 s
            h = False
        in Just $ fromForeignPtr f o (1,n) l h
    | otherwise =
        Nothing


-- | Possibly create a matrix view of a column vector.  This will fail
-- if the stride of the vector is not @1@ and the vector is not conjugated.
maybeFromCol :: (Elem e) => DVector t n e -> Maybe (DMatrix t (n,one) e)
maybeFromCol x
    | isConj x = maybeFromRow (conj x) >>= return . herm
    | strideOf x == 1 =
        let f = V.storageOf x
            o = V.offsetOf x
            m = dim x
            l = max 1 m
        in Just $ fromForeignPtr f o (m,1) l False
    | otherwise =
        Nothing

maybeToVector :: (Elem e) => DMatrix t (m,n) e -> Maybe (Order, DVector t k e)
maybeToVector a@(DM f o m n ld h) 
    | h = 
        maybeToVector (herm a) >>= (\(ord,x) -> return (flipOrder ord, conj x))
    | ld == m =
        Just $ (ColMajor, V.fromForeignPtr f o (m*n) 1  False)
    | m == 1 =
        Just $ (ColMajor, V.fromForeignPtr f o n     ld False)
    | otherwise =
        Nothing


-- | Take a unary elementwise vector operation and apply it to the elements of a matrix.
liftV :: (Elem e) => (DVector t k e -> IO ()) -> DMatrix t (m,n) e -> IO ()
liftV f a =
    case maybeToVector a of
        Just (_,x) -> f x
        _ -> 
            let xs  = case isHerm a of
                          True  -> rows (coerceMatrix a)
                          False -> cols (coerceMatrix a)
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
            let (xs,ys) = case isHerm a of
                               True  -> (rows (coerceMatrix a), rows (coerceMatrix b))
                               False -> (cols (coerceMatrix a), cols (coerceMatrix b))
            in zipWithM_ f xs ys


instance C.Matrix (DMatrix t) where
    numRows a | isHerm a  = size2 a
              | otherwise = size1 a

    numCols a | isHerm a  = size1 a
              | otherwise = size2 a

    herm a = a{ isHerm=(not . isHerm) a }


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
    
    newCopy a | isHerm a  = 
                  newCopy (herm a) >>= return . herm
              | otherwise = do
                  a' <- newMatrix_ (shape a)
                  liftV2 V.copyVector (unsafeThaw a') a
                  return a'
           
    unsafeReadElem a (i,j)
        | isHerm a  = unsafeReadElem (herm a) (j,i) >>= return . E.conj
        | otherwise = withForeignPtr (storageOf a) $ \ptr ->
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

    unsafeWriteElem a (i,j) e
        | isHerm a  = unsafeWriteElem (herm a) (j,i) $ E.conj e
        | otherwise = withForeignPtr (storageOf a) $ \ptr ->
                          pokeElemOff ptr (indexOf a (i,j)) e
    
    modifyWith f = liftV (modifyWith f)       


instance (BLAS1 e, Show e) => Show (DMatrix Imm (m,n) e) where
    show a | isHerm a = 
                "herm (" ++ show (herm a) ++ ")"
           | otherwise =
                "listMatrix " ++ show (shape a) ++ " " ++ show (elems a)
        
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
    