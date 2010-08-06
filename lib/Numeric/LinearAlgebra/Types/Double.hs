{-# LANGUAGE ForeignFunctionInterface #-}
{-# CFILES cbits/BLAS-double.c cbits/vectorOps-double.c #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Types.Double
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module Numeric.LinearAlgebra.Types.Double 
    where
        
import Foreign.Ptr ( Ptr )
import Numeric.LinearAlgebra.Types.CEnums

---------------------------- Vector Routines --------------------------------

foreign import ccall unsafe "vectorOps.h dVectorShift"
    vdShift :: Int -> Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vectorOps.h dVectorAdd"
    vdAdd :: Int -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vectorOps.h dVectorSub"
    vdSub :: Int -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vectorOps.h dVectorAxpby"
    vdAxpby :: Int -> Double -> Ptr Double -> Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vectorOps.h dVectorScale"
    vdScale :: Int -> Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vectorOps.h dVectorMul"
    vdMul :: Int -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vectorOps.h dVectorNeg"
    vdNeg :: Int -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vectorOps.h dVectorAbs"
    vdAbs :: Int -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vectorOps.h dVectorSgn"
    vdSgn :: Int -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vectorOps.h dVectorInv"
    vdInv :: Int -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vectorOps.h dVectorDiv"
    vdDiv :: Int -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "vectorOps.h dVectorExp"
    vdExp :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vectorOps.h dVectorSqrt"    
    vdSqrt :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vectorOps.h dVectorLog"
    vdLog :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vectorOps.h dVectorPow"
    vdPow :: Int -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()            
    
foreign import ccall unsafe "vectorOps.h dVectorSin"
    vdSin :: Int -> Ptr Double -> Ptr Double -> IO ()                
    
foreign import ccall unsafe "vectorOps.h dVectorCos"
    vdCos :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vectorOps.h dVectorTan"    
    vdTan :: Int -> Ptr Double -> Ptr Double -> IO ()                        

foreign import ccall unsafe "vectorOps.h dVectorASin"
    vdASin :: Int -> Ptr Double -> Ptr Double -> IO ()                
    
foreign import ccall unsafe "vectorOps.h dVectorACos"
    vdACos :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vectorOps.h dVectorATan"    
    vdATan :: Int -> Ptr Double -> Ptr Double -> IO ()                        

foreign import ccall unsafe "vectorOps.h dVectorSinh"
    vdSinh :: Int -> Ptr Double -> Ptr Double -> IO ()                
    
foreign import ccall unsafe "vectorOps.h dVectorCosh"
    vdCosh :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vectorOps.h dVectorTanh"    
    vdTanh :: Int -> Ptr Double -> Ptr Double -> IO ()                        

foreign import ccall unsafe "vectorOps.h dVectorASinh"
    vdASinh :: Int -> Ptr Double -> Ptr Double -> IO ()                
    
foreign import ccall unsafe "vectorOps.h dVectorACosh"
    vdACosh :: Int -> Ptr Double -> Ptr Double -> IO ()
    
foreign import ccall unsafe "vectorOps.h dVectorATanh"    
    vdATanh :: Int -> Ptr Double -> Ptr Double -> IO ()                        

                      
---------------------------- Level 1 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_ddot"
    ddot :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO Double

foreign import ccall unsafe "BLAS.h blas_dnrm2"
    dnrm2  :: Int -> Ptr Double -> Int -> IO Double

foreign import ccall unsafe "BLAS.h blas_dasum"
    dasum  :: Int -> Ptr Double -> Int -> IO Double

foreign import ccall unsafe "BLAS.h blas_idamax"
    idamax :: Int -> Ptr Double -> Int -> IO Int

foreign import ccall unsafe "BLAS.h blas_dscal"
    dscal  :: Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dswap"
    dswap  :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dcopy"
    dcopy  :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_daxpy"
    daxpy  :: Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_drotg"
    drotg  :: Ptr Double -> Ptr Double -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall unsafe "BLAS.h blas_drot"
    drot :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Double -> IO ()

foreign import ccall unsafe "BLAS.h blas_drotmg"
    drotmg :: Ptr Double -> Ptr Double -> Ptr Double -> Double -> Ptr Double -> IO ()

foreign import ccall unsafe "BLAS.h blas_drotm"
    drotm :: Int -> Ptr Double -> Int -> Ptr Double -> Int -> Ptr Double -> IO ()


---------------------------- Level 2 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_dgemv"
    dgemv :: CBLASTrans -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dgbmv"
    dgbmv ::  CBLASTrans -> Int -> Int -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dtrmv"
    dtrmv ::  CBLASUplo -> CBLASTrans -> CBLASDiag -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dtbmv"
    dtbmv ::  CBLASUplo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()
                 
foreign import ccall unsafe "BLAS.h blas_dtrsv"
    dtrsv ::  CBLASUplo -> CBLASTrans -> CBLASDiag -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dtbsv"
    dtbsv ::  CBLASUplo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()
    
foreign import ccall unsafe "BLAS.h blas_dsymv"
    dsymv ::  CBLASUplo -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
    
foreign import ccall unsafe "BLAS.h blas_dsbmv"
    dsbmv ::  CBLASUplo -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
    
foreign import ccall unsafe "BLAS.h blas_dger"
    dger  ::  Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()
        
foreign import ccall unsafe "BLAS.h blas_dsyr"
    dsyr  ::  CBLASUplo -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dsyr2"
    dsyr2 ::  CBLASUplo -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()


---------------------------- Level 3 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_dgemm"
    dgemm  ::  CBLASTrans -> CBLASTrans -> Int -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dsymm"
    dsymm  ::  CBLASSide -> CBLASUplo -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dtrmm"
    dtrmm  ::  CBLASSide -> CBLASUplo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dtrsm"
    dtrsm  ::  CBLASSide -> CBLASUplo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_dsyrk"
    dsyrk  ::  CBLASUplo -> CBLASTrans -> Int -> Int -> Double -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
           
foreign import ccall unsafe "BLAS.h blas_dsyr2k"           
    dsyr2k ::  CBLASUplo -> CBLASTrans -> Int -> Int -> Double -> Ptr Double -> Int -> Ptr Double -> Int -> Double -> Ptr Double -> Int -> IO ()
    