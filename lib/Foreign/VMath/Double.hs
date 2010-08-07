{-# LANGUAGE ForeignFunctionInterface #-}
{-# CFILES cbits/vmath-double.c #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.VMath.Double
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Foreign.VMath.Double
    where
        
import Foreign.Ptr( Ptr )

---------------------------- Vector Routines --------------------------------

foreign import ccall unsafe "vdShift"
    vdShift :: Int -> Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vdAdd"
    vdAdd :: Int -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vdSub"
    vdSub :: Int -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vdAxpby"
    vdAxpby :: Int -> Double -> Ptr Double -> Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vdScale"
    vdScale :: Int -> Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vdMul"
    vdMul :: Int -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vdNeg"
    vdNeg :: Int -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vdAbs"
    vdAbs :: Int -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vdSgn"
    vdSgn :: Int -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vdInv"
    vdInv :: Int -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vdDiv"
    vdDiv :: Int -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vdExp"
    vdExp :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vdSqrt"    
    vdSqrt :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vdLog"
    vdLog :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vdPow"
    vdPow :: Int -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()            
    
foreign import ccall unsafe "vdSin"
    vdSin :: Int -> Ptr Double -> Ptr Double -> IO ()                
    
foreign import ccall unsafe "vdCos"
    vdCos :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vdTan"    
    vdTan :: Int -> Ptr Double -> Ptr Double -> IO ()                        

foreign import ccall unsafe "vdASin"
    vdASin :: Int -> Ptr Double -> Ptr Double -> IO ()                
    
foreign import ccall unsafe "vdACos"
    vdACos :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vdATan"    
    vdATan :: Int -> Ptr Double -> Ptr Double -> IO ()                        

foreign import ccall unsafe "vdSinh"
    vdSinh :: Int -> Ptr Double -> Ptr Double -> IO ()                
    
foreign import ccall unsafe "vdCosh"
    vdCosh :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vdTanh"    
    vdTanh :: Int -> Ptr Double -> Ptr Double -> IO ()                        

foreign import ccall unsafe "vdASinh"
    vdASinh :: Int -> Ptr Double -> Ptr Double -> IO ()                
    
foreign import ccall unsafe "vdACosh"
    vdACosh :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vdATanh"    
    vdATanh :: Int -> Ptr Double -> Ptr Double -> IO ()                        

                      
