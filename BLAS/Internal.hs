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
    ) where

import Foreign                  ( Ptr, Storable, castPtr, sizeOf )
import Foreign.C.Types          ( CSize )

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

