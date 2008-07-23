{-# LANGUAGE BangPatterns, FlexibleInstances, MultiParamTypeClasses #-}
{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Vector.Dense.Internal
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Vector.Dense.Internal (
    -- * Vector data types
    Vector,
    IOVector,
    DVector(..),

    module BLAS.Vector,
    module BLAS.Tensor,

    -- * Conversion to and from @ForeignPtr@s.
    fromForeignPtr,
    toForeignPtr,
    
    -- * Creating new vectors
    newVector, 
    newVector_,
    newListVector,
    
    -- * Special vectors
    newBasis,
    setBasis,

    -- * Vector views
    subvector,
    subvectorWithStride,
    
    -- * Casting vectors
    coerceVector,
    
    -- * Unsafe operations
    unsafeNewVector,
    unsafeWithElemPtr,
    unsafeSubvector,
    unsafeSubvectorWithStride,
    
    unsafeFreeze,
    unsafeThaw,
    
    ) where
  
  
import Control.Monad
import Data.Ix
import Foreign 
import System.IO.Unsafe       
import Unsafe.Coerce

import Data.AEq

import BLAS.Access
import BLAS.Elem.Base ( Elem )
import qualified BLAS.Elem.Base as E
import BLAS.Vector hiding ( Vector )
import qualified BLAS.Vector as C
import BLAS.Tensor

import BLAS.Internal  ( clearArray, inlinePerformIO, checkedSubvector,
    checkedSubvectorWithStride )
import BLAS.C.Level1  ( BLAS1, copy )


-- | A dense vector.  @t@ is a type that will usually be @Imm@ or @Mut@.  
-- @n@ is a phantom type for the dimension of the vector, and @e@ is the 
-- element type.  A @DVector@ @x@ stores @dim x@ elements.  Indices into
-- the vector are @0@-based.
data DVector t n e = 
      DV { storageOf :: {-# UNPACK #-} !(ForeignPtr e) -- ^ a pointer to the storage region
         , offsetOf  :: {-# UNPACK #-} !Int            -- ^ an offset (in elements, not bytes) to the first element in the vector. 
         , lengthOf  :: {-# UNPACK #-} !Int            -- ^ the length of the vector
         , strideOf  :: {-# UNPACK #-} !Int            -- ^ the stride (in elements, not bytes) between elements.
         , isConj    :: {-# UNPACK #-} !Bool           -- ^ indicates whether or not the vector is conjugated
         }

type Vector n e = DVector Imm n e
type IOVector n e = DVector Mut n e

-- | Cast the phantom length type.
coerceVector :: DVector t n e -> DVector t m e
coerceVector = unsafeCoerce

-- | @fromForeignPtr fptr offset n inc c@ creates a vector view of a
-- region in memory starting at the given offset and having dimension @n@,
-- with a stride of @inc@, and with @isConj@ set to @c@.
fromForeignPtr :: ForeignPtr e -> Int -> Int -> Int -> Bool -> DVector t n e
fromForeignPtr = DV
{-# INLINE fromForeignPtr #-}

-- | Gets the tuple @(fptr,offset,n,inc,c)@, where @n@ is the dimension and 
-- @inc@ is the stride of the vector, and @c@ indicates whether or not the
-- vector is conjugated.
toForeignPtr :: DVector t n e -> (ForeignPtr e, Int, Int, Int, Bool)
toForeignPtr (DV f o n s c) = (f, o, n, s, c)
{-# INLINE toForeignPtr #-}

-- | @subvector x o n@ creates a subvector view of @x@ starting at index @o@ 
-- and having length @n@.
subvector :: DVector t n e -> Int -> Int -> DVector t m e
subvector x = checkedSubvector (dim x) (unsafeSubvector x)

-- | Same as 'subvector' but arguments are not range-checked.
unsafeSubvector :: DVector t n e -> Int -> Int -> DVector t m e
unsafeSubvector = unsafeSubvectorWithStride 1

-- | @subvectorWithStride s x o n@ creates a subvector view of @x@ starting 
-- at index @o@, having length @n@ and stride @s@.
subvectorWithStride :: Int -> DVector t n e -> Int -> Int -> DVector t m e
subvectorWithStride s x = 
    checkedSubvectorWithStride s (dim x) (unsafeSubvectorWithStride s x)
    
-- | Same as 'subvectorWithStride' but arguments are not range-checked.
unsafeSubvectorWithStride :: Int -> DVector t n e -> Int -> Int -> DVector t m e
unsafeSubvectorWithStride s x o n =
    let f  = storageOf x
        o' = indexOf x o
        n' = n
        s' = s * (strideOf x)
        c  = isConj x
    in 
        fromForeignPtr f o' n' s' c

-- | Creates a new vector of the given length.  The elements will be 
-- uninitialized.
newVector_ :: (Elem e) => Int -> IO (DVector t n e)
newVector_ n
    | n < 0 = 
        ioError $ userError $ 
            "Tried to create a vector with `" ++ show n ++ "' elements."
    | otherwise = do
        arr <- mallocForeignPtrArray n
        return $ fromForeignPtr arr 0 n 1 False

-- | Creates a new vector of the given dimension with the given elements.
-- If the list has length less than the passed-in dimenson, the tail of
-- the vector will be uninitialized.
newListVector :: (Elem e) => Int -> [e] -> IO (DVector t n e)
newListVector n es = do
    x <- newVector_ n
    withForeignPtr (storageOf x) $ flip pokeArray $ take n es
    return x

-- | @listVector n es@ is equivalent to @vector n (zip [0..(n-1)] es)@, except
-- that the result is undefined if @length es@ is less than @n@.
listVector :: (Elem e) => Int -> [e] -> Vector n e
listVector n es = unsafeFreeze $ unsafePerformIO $ newListVector n es
{-# NOINLINE listVector #-}


-- | Creates a new vector with the given association list.  Unspecified
-- indices will get initialized to zero.
newVector :: (BLAS1 e) => Int -> [(Int,e)] -> IO (DVector t n e)
newVector =
    newVectorHelp writeElem

-- | Same as 'newVector' but indices are not range-checked.
unsafeNewVector :: (BLAS1 e) => Int -> [(Int,e)] -> IO (DVector t n e)
unsafeNewVector =
    newVectorHelp unsafeWriteElem

newVectorHelp :: (BLAS1 e) => 
       (IOVector n e -> Int -> e -> IO ()) 
    -> Int -> [(Int,e)] -> IO (DVector t n e)
newVectorHelp set n ies = do
    x <- newZero n
    mapM_ (uncurry $ set x) ies
    return (unsafeCoerce x)

-- | @newBasis n i@ creates a vector of length @n@ that is all zero except for
-- at position @i@, where it equal to one.
newBasis :: (BLAS1 e) => Int -> Int -> IO (IOVector n e)
newBasis n i = do
    x <- newVector_ n
    setBasis i x
    return x

-- | @setBasis x i@ sets the @i@th coordinate of @x@ to @1@, and all other
-- coordinates to @0@.  If the vector has been scaled, it is possible that
-- @readVector x i@ will not return exactly @1@.  See 'setElem'.
setBasis :: (BLAS1 e) => Int -> IOVector n e -> IO ()
setBasis i x
    | i < 0 || i >= dim x =
        ioError $ userError $ 
            "tried to set a vector of dimension `" ++ show (dim x) ++ "'"
            ++ " to basis vector `" ++ show i ++ "'"
    | otherwise = do
        setZero x
        unsafeWriteElem x i 1 

indexOf :: DVector t n e -> Int -> Int
indexOf x i = offsetOf x + i * strideOf x
{-# INLINE indexOf #-}

-- | Evaluate a function with a pointer to the value stored at the given
-- index.  Note that the value may need to conjugated before using it.  See
-- 'isConj'.
unsafeWithElemPtr :: (Elem e) => DVector t n e -> Int -> (Ptr e -> IO a) -> IO a
unsafeWithElemPtr x i f =
    withForeignPtr (storageOf x) $ \ptr ->
        let elemPtr = ptr `advancePtr` (indexOf x i)
        in f elemPtr
{-# INLINE unsafeWithElemPtr #-}

-- | Cast the access type to @Imm@.
unsafeFreeze :: DVector t n e -> Vector n e
unsafeFreeze = unsafeCoerce

-- | Cast the access type to @Mut@.
unsafeThaw :: DVector t n e -> IOVector n e
unsafeThaw = unsafeCoerce

instance C.Vector (DVector t) where
    dim = lengthOf
    {-# INLINE dim #-}

    conj x = let c' = (not . isConj) x 
             in x { isConj=c' }
    {-# INLINE conj #-}

{-# RULES 
    "conj/Float"  conj = conjFloat 
    "conj/Double" conj = conjDouble
  #-}
conjFloat :: DVector t n Float -> DVector t n Float
conjFloat = id

conjDouble :: DVector t n Double -> DVector t n Double
conjDouble = id



instance Tensor (DVector t n) Int e where
    shape = dim
    {-# INLINE shape #-}

    bounds x = (0, dim x - 1)
    {-# INLINE bounds #-}

instance (BLAS1 e) => ITensor (DVector Imm n) Int e where
    size = dim
    
    indices = range . bounds
    {-# INLINE indices #-}

    elems  = inlinePerformIO . getElems . unsafeThaw
    assocs = inlinePerformIO . getAssocs . unsafeThaw

    unsafeAt x = inlinePerformIO . unsafeReadElem (unsafeThaw x)
    {-# INLINE unsafeAt #-}
    
    amap f x = listVector (dim x) (map f $ elems x)
    
    (//) = replaceHelp writeElem
    unsafeReplace = replaceHelp unsafeWriteElem

replaceHelp :: (BLAS1 e) => 
       (IOVector n e -> Int -> e -> IO ())
    -> Vector n e -> [(Int, e)] -> Vector n e
replaceHelp set x ies =
    unsafeFreeze $ unsafePerformIO $ do
        y  <- newCopy (unsafeThaw x)
        mapM_ (uncurry $ set y) ies
        return y
{-# NOINLINE replaceHelp #-}


instance (BLAS1 e) => IDTensor (DVector Imm n) Int e where
    zero n = unsafeFreeze $ unsafePerformIO $ newZero n
    {-# NOINLINE zero #-}
    
    constant n e = unsafeFreeze $ unsafePerformIO $ newConstant n e
    {-# NOINLINE constant #-}

    azipWith f x y
        | dim y /= n =
            error ("amap2: vector lengths differ; first has length `" ++
                    show n ++ "' and second has length `" ++
                    show (dim y) ++ "'")
        | otherwise =
            listVector n (zipWith f (elems x) (elems y))
      where
        n = dim x
    

instance (BLAS1 e) => RTensor (DVector t n) Int e IO where
    getSize = return . dim
    
    newCopy x
        | isConj x = 
            newCopy (conj x) >>= return . conj
        | otherwise = do
            y <- newVector_ (dim x)
            unsafeWithElemPtr x 0 $ \pX ->
                unsafeWithElemPtr y 0 $ \pY ->
                    let n    = dim x
                        incX = strideOf x
                        incY = strideOf y
                    in copy n pX incX pY incY >> 
                       return y
                       
    getIndices = return . indices . unsafeFreeze
    {-# INLINE getIndices #-}
    
    unsafeReadElem x i
        | isConj x = 
            unsafeReadElem (conj x) i >>= return . E.conj
        | otherwise =
            withForeignPtr (storageOf x) $ \ptr ->
                peekElemOff ptr (indexOf x i) 

    getAssocs x
        | isConj x =
            getAssocs (conj x) >>= return . map (\(i,e) -> (i,E.conj e))
        | otherwise =
            let (f,o,n,incX,_) = toForeignPtr x
                ptr = (unsafeForeignPtrToPtr f) `advancePtr` o
            in return $ go n f incX ptr 0
      where
            go !n !f !incX !ptr !i 
                | i >= n = 
                     -- This is very important since we are doing unsafe IO.
                     -- Otherwise, the DVector might get discared and the
                     -- memory freed before all of the elements are read
                     inlinePerformIO $ do
                         touchForeignPtr f
                         return []
                | otherwise =
                    let e    = inlinePerformIO $ peek ptr
                        ptr' = ptr `advancePtr` incX
                        i'   = i + 1
                        ies  = go n f incX ptr' i'
                    in e `seq` ((i,e):ies)
    {-# NOINLINE getAssocs #-}


instance (BLAS1 e) => RDTensor (DVector t n) Int e IO where
    newZero n = newVector_ n >>= (\x -> setZero (unsafeThaw x) >> return x)
    
    newConstant n e = newVector_ n >>= (\x -> setConstant e (unsafeThaw x) >> return x)
    
    
instance (BLAS1 e) => MTensor (DVector Mut n) Int e IO where
    setZero x 
        | strideOf x == 1 = unsafeWithElemPtr x 0 $ flip clearArray (dim x)
        | otherwise       = setConstant 0 x

    setConstant e x 
        | isConj x  = setConstant (E.conj e) (conj x)
        | otherwise = unsafeWithElemPtr x 0 $ go (dim x)
      where
        go !n !ptr | n <= 0 = return ()
                   | otherwise = let ptr' = ptr `advancePtr` (strideOf x)
                                     n'   = n - 1
                                 in poke ptr e >> 
                                    go n' ptr'
    
    unsafeWriteElem x i e =
        let e' = if isConj x then E.conj e else e
        in withForeignPtr (storageOf x) $ \ptr -> 
               pokeElemOff ptr (indexOf x i) e'
                        
    canModifyElem x i = return $ inRange (bounds x) i
    {-# INLINE canModifyElem #-}
    
    modifyWith f x
        | isConj x  = modifyWith (E.conj . f . E.conj) (conj x)
        | otherwise = withForeignPtr (storageOf x) $ go (dim x)
      where
        go !n !ptr | n <= 0 = return ()
                   | otherwise = do
                       peek ptr >>= poke ptr . f
                       go (n-1) (ptr `advancePtr` incX)

        incX = strideOf x
    
compareHelp :: (BLAS1 e) => 
    (e -> e -> Bool) -> Vector n e -> Vector n e -> Bool
compareHelp cmp x y
    | isConj x && isConj y =
        compareHelp cmp (conj x) (conj y)
    | otherwise =
        (dim x == dim y) && (and $ zipWith cmp (elems x) (elems y))

instance (BLAS1 e, Eq e) => Eq (DVector Imm n e) where
    (==) = compareHelp (==)

instance (BLAS1 e, AEq e) => AEq (DVector Imm n e) where
    (===) = compareHelp (===)
    (~==) = compareHelp (~==)

instance (BLAS1 e, Show e) => Show (DVector Imm n e) where
    show x
        | isConj x  = "conj (" ++ show (conj x) ++ ")"
        | otherwise = "listVector " ++ show (dim x) ++ " " ++ show (elems x)
    