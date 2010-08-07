{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.BLAS.Level3
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Matrix-Matrix operations.
--

module Foreign.BLAS.Level3 (
    BLAS3(..),
    ) where
     
import Data.Complex 
import Foreign( Ptr, Storable, with )

import Foreign.BLAS.Types
import Foreign.BLAS.Level2
import Foreign.BLAS.Double
import Foreign.BLAS.Zomplex

-- | Types with matrix-matrix operations.        
class (BLAS2 a) => BLAS3 a where
    gemm  :: Trans -> Trans -> Int -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    symm  :: Side -> Uplo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    hemm  :: Side -> Uplo -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    trmm  :: Side -> Uplo -> Trans -> Diag -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    trsm  :: Side -> Uplo -> Trans -> Diag -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> IO ()
    syrk  :: Uplo -> Trans -> Int -> Int -> a -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    syr2k :: Uplo -> Trans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    herk  :: Uplo -> Trans -> Int -> Int -> a -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
    her2k :: Uplo -> Trans -> Int -> Int -> a -> Ptr a -> Int -> Ptr a -> Int -> a -> Ptr a -> Int -> IO ()
  

withEnum :: (Enum a, Storable a) => Int -> (Ptr a -> IO b) -> IO b
withEnum = with . toEnum
{-# INLINE withEnum #-}
  
    
instance BLAS3 Double where
    gemm transa transb m n k alpha pa lda pb ldb beta pc ldc =
        withTrans transa $ \ptransa ->
        withTrans transb $ \ptransb ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
        with beta $ \pbeta ->
        withEnum ldc $ \pldc ->
            dgemm ptransa ptransb pm pn pk palpha pa plda pb pldb pbeta pc pldc
    {-# INLINE gemm #-}

    symm side uplo m n alpha pa lda pb ldb beta pc ldc =
        withSide side $ \pside ->
        withUplo uplo $ \puplo ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
        with beta $ \pbeta ->
        withEnum ldc $ \pldc ->
            dsymm pside puplo pm pn palpha pa plda pb pldb pbeta pc pldc
    {-# INLINE symm #-}

    hemm = symm
    {-# INLINE hemm #-}

    trmm side uplo transa diag m n alpha pa lda pb ldb =
        withSide side $ \pside ->
        withUplo uplo $ \puplo ->
        withTrans transa $ \ptransa ->
        withDiag diag $ \pdiag ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
            dtrmm pside puplo ptransa pdiag pm pn palpha pa plda pb pldb
    {-# INLINE trmm #-}

    trsm side uplo transa diag m n alpha pa lda pb ldb =
        withSide side $ \pside ->
        withUplo uplo $ \puplo ->
        withTrans transa $ \ptransa ->
        withDiag diag $ \pdiag ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
            dtrsm pside puplo ptransa pdiag pm pn palpha pa plda pb pldb
    {-# INLINE trsm #-}

    syrk uplo transa n k alpha pa lda beta pc ldc =
        withUplo uplo $ \puplo ->
        withTrans transa $ \ptransa ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        with beta $ \pbeta ->
        withEnum ldc $ \pldc ->
            dsyrk puplo ptransa pn pk palpha pa plda pbeta pc pldc
    {-# INLINE syrk #-}

    syr2k uplo transa n k alpha pa lda pb ldb beta pc ldc =
        withUplo uplo $ \puplo ->
        withTrans transa $ \ptransa ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
        with beta $ \pbeta ->
        withEnum ldc $ \pldc ->
            dsyr2k puplo ptransa pn pk palpha pa plda pb pldb pbeta pc pldc
    {-# INLINE syr2k #-}

    herk = syrk
    {-# INLINE herk #-}

    her2k = syr2k
    {-# INLINE her2k #-}

    
instance BLAS3 (Complex Double) where
    gemm transa transb m n k alpha pa lda pb ldb beta pc ldc =
        withTrans transa $ \ptransa ->
        withTrans transb $ \ptransb ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
        with beta $ \pbeta ->
        withEnum ldc $ \pldc ->
            zgemm ptransa ptransb pm pn pk palpha pa plda pb pldb pbeta pc pldc
    {-# INLINE gemm #-}

    symm side uplo m n alpha pa lda pb ldb beta pc ldc =
        withSide side $ \pside ->
        withUplo uplo $ \puplo ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
        with beta $ \pbeta ->
        withEnum ldc $ \pldc ->
            zsymm pside puplo pm pn palpha pa plda pb pldb pbeta pc pldc
    {-# INLINE symm #-}

    hemm side uplo m n alpha pa lda pb ldb beta pc ldc =
        withSide side $ \pside ->
        withUplo uplo $ \puplo ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
        with beta $ \pbeta ->
        withEnum ldc $ \pldc ->
            zhemm pside puplo pm pn palpha pa plda pb pldb pbeta pc pldc
    {-# INLINE hemm #-}

    trmm side uplo transa diag m n alpha pa lda pb ldb =
        withSide side $ \pside ->
        withUplo uplo $ \puplo ->
        withTrans transa $ \ptransa ->
        withDiag diag $ \pdiag ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
            ztrmm pside puplo ptransa pdiag pm pn palpha pa plda pb pldb
    {-# INLINE trmm #-}

    trsm side uplo transa diag m n alpha pa lda pb ldb =
        withSide side $ \pside ->
        withUplo uplo $ \puplo ->
        withTrans transa $ \ptransa ->
        withDiag diag $ \pdiag ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
            ztrsm pside puplo ptransa pdiag pm pn palpha pa plda pb pldb
    {-# INLINE trsm #-}

    syrk uplo transa n k alpha pa lda beta pc ldc =
        withUplo uplo $ \puplo ->
        withTrans transa $ \ptransa ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        with beta $ \pbeta ->
        withEnum ldc $ \pldc ->
            zsyrk puplo ptransa pn pk palpha pa plda pbeta pc pldc
    {-# INLINE syrk #-}

    syr2k uplo transa n k alpha pa lda pb ldb beta pc ldc =
        withUplo uplo $ \puplo ->
        withTrans transa $ \ptransa ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
        with beta $ \pbeta ->
        withEnum ldc $ \pldc ->
            zsyr2k puplo ptransa pn pk palpha pa plda pb pldb pbeta pc pldc
    {-# INLINE syr2k #-}

    herk uplo transa n k alpha pa lda beta pc ldc =
        withUplo uplo $ \puplo ->
        withTrans transa $ \ptransa ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        with beta $ \pbeta ->
        withEnum ldc $ \pldc ->
            zherk puplo ptransa pn pk palpha pa plda pbeta pc pldc
    {-# INLINE herk #-}

    her2k uplo transa n k alpha pa lda pb ldb beta pc ldc =
        withUplo uplo $ \puplo ->
        withTrans transa $ \ptransa ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        with alpha $ \palpha ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
        with beta $ \pbeta ->
        withEnum ldc $ \pldc ->
            zher2k puplo ptransa pn pk palpha pa plda pb pldb pbeta pc pldc
    {-# INLINE her2k #-}
