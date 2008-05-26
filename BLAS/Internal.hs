{-# LANGUAGE CPP, ForeignFunctionInterface #-}
{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Internal
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--


module BLAS.Internal (
    clearArray,
    bzero,
    inlinePerformIO,
    checkedSubvector,
    checkedSubvectorWithStride,
    checkVecVecOp,
    ) where

import Data.Ix     ( inRange )
import Foreign                  ( Ptr, Storable, castPtr, sizeOf )
import Foreign.C.Types          ( CSize )
import Text.Printf ( printf )

#if defined(__GLASGOW_HASKELL__)
import GHC.Base                 ( realWorld# )
import GHC.IOBase               ( IO(IO) )
#else
import System.IO.Unsafe         ( unsafePerformIO )
#endif


clearArray :: Storable e => Ptr e -> Int -> IO ()
clearArray = clearArray' undefined
    where
    clearArray' :: Storable e => e -> Ptr e -> Int -> IO ()
    clearArray' e ptr n =
        let nbytes = (fromInteger . toInteger) (n * sizeOf e)
        in do
            bzero ptr nbytes
{-# INLINE clearArray #-}


bzero :: Ptr a -> Int -> IO ()
bzero ptr n =
    let ptr' = castPtr ptr
        n'   = (fromInteger . toInteger) n
    in bzero_ ptr' n'
        
foreign import ccall "strings.h bzero"
    bzero_ :: Ptr () -> CSize -> IO ()
    

inlinePerformIO :: IO a -> a
#if defined(__GLASGOW_HASKELL__)
inlinePerformIO (IO m) = case m realWorld# of (# _, r #) -> r
#else
inlinePerformIO = unsafePerformIO
#endif
{-# INLINE inlinePerformIO #-}

checkedSubvector :: Int -> (Int -> Int -> v) -> Int -> Int -> v
checkedSubvector n sub o n'
    | (o < 0) && (n' /= 0) = 
        error $ printf 
            "tried to create a subvector starting at a negative offset: `%d'" o
    | n' < 0 = 
        error $ printf 
            "tried to create a subvector with a negative length `%d'" n'
    | n' + o > n = 
        error $ printf
            ("tried to create a subvector of length `%d' and offset `%d' "
             ++ " from a vector of length `%d'") n' o n
    | otherwise =
        sub o n'
        

checkedSubvectorWithStride :: Int -> Int -> (Int -> Int -> v)
    -> Int -> Int -> v
checkedSubvectorWithStride s n sub o n'
    | (o < 0) && (n' /= 0) =
        error $ printf
            "Tried to create a subvector starting at a negative offset: `%d'" o
    | n' < 0 =
        error $ printf
            "Tried to create a subvector with a negative length `%d'" n'
    | s <= 0 =
        error $ printf
            "Tried to create a subvector with non-positive stride `%d'" s
    | not $ inRange (-1,n) (o + s * n') =
        error $ printf
            ("tried to create a subvector of length `%d',  offset `%d',"
             ++ " and stride '%d' from a vector of length `%d'") n' o s n
    | otherwise =
        sub o n'

checkVecVecOp :: String -> Int -> Int -> IO ()
checkVecVecOp name n1 n2
    | n1 /= n2 =
        ioError $ userError $ printf
            ("%s: x and y have different dimensions.  x has dimension `%d',"
             ++ " and y has dimension `%d'") name n1 n2
    | otherwise =
        return ()
