{-# LANGUAGE ForeignFunctionInterface #-}
{-# CFILES cbits/vmath-zomplex.c #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.VMath.Zomplex
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Foreign.VMath.Zomplex
    where
        
import Data.Complex ( Complex )
import Foreign.Ptr  ( Ptr )

foreign import ccall unsafe "vzConj"
    vzConj :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzShift"
    vzShift :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzAdd"
    vzAdd :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzSub"
    vzSub :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzAxpby"
    vzAxpby :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzScale"
    vzScale :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzMul"
    vzMul :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzNeg"
    vzNeg :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzAbs"
    vzAbs :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzSgn"
    vzSgn :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzInv"
    vzInv :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzDiv"
    vzDiv :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vzExp"
    vzExp :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vzSqrt"    
    vzSqrt :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vzLog"
    vzLog :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vzPow"
    vzPow :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()            
    
foreign import ccall unsafe "vzSin"
    vzSin :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                
    
foreign import ccall unsafe "vzCos"
    vzCos :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vzTan"    
    vzTan :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                        

foreign import ccall unsafe "vzASin"
    vzASin :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                
    
foreign import ccall unsafe "vzACos"
    vzACos :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vzATan"    
    vzATan :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                        

foreign import ccall unsafe "vzSinh"
    vzSinh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                
    
foreign import ccall unsafe "vzCosh"
    vzCosh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vzTanh"    
    vzTanh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                        

foreign import ccall unsafe "vzASinh"
    vzASinh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                
    
foreign import ccall unsafe "vzACosh"
    vzACosh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vzATanh"    
    vzATanh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                        

