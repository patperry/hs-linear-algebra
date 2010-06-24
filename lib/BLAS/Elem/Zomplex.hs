{-# LANGUAGE ForeignFunctionInterface #-}
{-# CFILES cbits/zomplex.c #-}
-----------------------------------------------------------------------------
-- |
-- Module     : BLAS.Elem.Zomplex
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--

module BLAS.Elem.Zomplex
    where
        
import Data.Complex ( Complex )
import Foreign.Ptr  ( Ptr )
import BLAS.CTypes

foreign import ccall unsafe "vectorOps.h zVectorConj"
    vzConj :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorShift"
    vzShift :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorAdd"
    vzAdd :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorSub"
    vzSub :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorAxpby"
    vzAxpby :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorScale"
    vzScale :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorMul"
    vzMul :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorNeg"
    vzNeg :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorAbs"
    vzAbs :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorSgn"
    vzSgn :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorInv"
    vzInv :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorDiv"
    vzDiv :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "vectorOps.h zVectorExp"
    vzExp :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vectorOps.h zVectorSqrt"    
    vzSqrt :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vectorOps.h zVectorLog"
    vzLog :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vectorOps.h zVectorPow"
    vzPow :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()            
    
foreign import ccall unsafe "vectorOps.h zVectorSin"
    vzSin :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                
    
foreign import ccall unsafe "vectorOps.h zVectorCos"
    vzCos :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vectorOps.h zVectorTan"    
    vzTan :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                        

foreign import ccall unsafe "vectorOps.h zVectorASin"
    vzASin :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                
    
foreign import ccall unsafe "vectorOps.h zVectorACos"
    vzACos :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vectorOps.h zVectorATan"    
    vzATan :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                        

foreign import ccall unsafe "vectorOps.h zVectorSinh"
    vzSinh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                
    
foreign import ccall unsafe "vectorOps.h zVectorCosh"
    vzCosh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vectorOps.h zVectorTanh"    
    vzTanh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                        

foreign import ccall unsafe "vectorOps.h zVectorASinh"
    vzASinh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                
    
foreign import ccall unsafe "vectorOps.h zVectorACosh"
    vzACosh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
    
foreign import ccall unsafe "vectorOps.h zVectorATanh"    
    vzATanh :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()                        


---------------------------- Level 1 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_zdotu_sub"
    zdotu_sub   :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "BLAS.h blas_zdotc_sub"
    zdotc_sub   :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> IO ()


foreign import ccall unsafe "BLAS.h blas_dznrm2"
    znrm2  :: Int -> Ptr (Complex Double) -> Int -> IO Double

foreign import ccall unsafe "BLAS.h blas_dzasum"
    zasum  :: Int -> Ptr (Complex Double) -> Int -> IO Double

foreign import ccall unsafe "BLAS.h blas_izamax"
    izamax :: Int -> Ptr (Complex Double) -> Int -> IO Int

foreign import ccall unsafe "BLAS.h blas_zscal"
    zscal  :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zswap"
    zswap  :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zcopy"
    zcopy  :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zaxpy"
    zaxpy  :: Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zrotg"
    zrotg  :: Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import ccall unsafe "BLAS.h blas_zdrot"
    zdrot :: Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Double -> Double -> IO ()


---------------------------- Level 2 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_zgemv"
    zgemv :: CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zgbmv"
    zgbmv ::  CBLASTrans -> Int -> Int -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_ztrmv"
    ztrmv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_ztbmv"
    ztbmv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()
                 
foreign import ccall unsafe "BLAS.h blas_ztrsv"
    ztrsv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_ztbsv"
    ztbsv ::  CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()
    
foreign import ccall unsafe "BLAS.h blas_zhemv"
    zhemv ::  CBLASUpLo -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zhbmv"
    zhbmv ::  CBLASUpLo -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
    
foreign import ccall unsafe "BLAS.h blas_zgeru"
    zgeru  ::  Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zgerc"
    zgerc  ::  Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()
        
foreign import ccall unsafe "BLAS.h blas_zher"
    zher  ::  CBLASUpLo -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zher2"
    zher2 ::  CBLASUpLo -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()


---------------------------- Level 3 Routines -------------------------------

foreign import ccall unsafe "BLAS.h blas_zgemm"
    zgemm  ::  CBLASTrans -> CBLASTrans -> Int -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zsymm"
    zsymm  ::  CBLASSide -> CBLASUpLo -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zhemm"
    zhemm  ::  CBLASSide -> CBLASUpLo -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_ztrmm"
    ztrmm  ::  CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_ztrsm"
    ztrsm  ::  CBLASSide -> CBLASUpLo -> CBLASTrans -> CBLASDiag -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zsyrk"
    zsyrk  ::  CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
           
foreign import ccall unsafe "BLAS.h blas_zsyr2k"           
    zsyr2k ::  CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()

foreign import ccall unsafe "BLAS.h blas_zherk"
    zherk  ::  CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
           
foreign import ccall unsafe "BLAS.h blas_zher2k"           
    zher2k ::  CBLASUpLo -> CBLASTrans -> Int -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Int -> Ptr (Complex Double) -> Ptr (Complex Double) -> Int -> IO ()
    