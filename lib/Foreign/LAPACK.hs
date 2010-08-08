{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Foreign.LAPACK
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Types with a Linear Algebra PACKage.
--
module Foreign.LAPACK (
    LAPACK(..),
    ) where

import Control.Exception( assert )
import Data.Complex( Complex )
import Foreign
import Foreign.BLAS
import Foreign.LAPACK.Double
import Foreign.LAPACK.Zomplex


-- | Types with LAPACK operations.
class (BLAS3 e) => LAPACK e where
    geqrf :: Int -> Int -> Ptr e -> Int -> Ptr e -> IO ()
    gelqf :: Int -> Int -> Ptr e -> Int -> Ptr e -> IO ()
    larfg :: Int -> Ptr e -> Ptr e -> Int -> IO e
    potrf :: Uplo -> Int -> Ptr e -> Int -> IO Int
    potrs :: Uplo -> Int -> Int -> Ptr e -> Int -> Ptr e -> Int -> IO ()
    
    unmqr :: Side -> Trans -> Int -> Int -> Int -> Ptr e -> Int -> Ptr e -> Ptr e -> Int -> IO ()
    unmlq :: Side -> Trans -> Int -> Int -> Int -> Ptr e -> Int -> Ptr e -> Ptr e -> Int -> IO ()


callWithWork :: (Storable e) => (Ptr e -> Ptr LAInt -> Ptr LAInt -> IO ()) -> IO LAInt
callWithWork call =
    alloca $ \pquery ->
    with (-1) $ \pnone ->
    alloca $ \pinfo -> do
        call pquery pnone pinfo
        lwork <- (ceiling . max 1) `fmap` (peek (castPtr pquery) :: IO Double)
        allocaArray lwork $ \pwork ->
            withEnum lwork $ \plwork -> do
                call pwork plwork pinfo
                peek pinfo

checkInfo :: (Num info) => info -> IO ()
checkInfo info = assert (info == 0) $ return ()

withEnum :: (Enum a, Storable a) => Int -> (Ptr a -> IO b) -> IO b
withEnum = with . toEnum
{-# INLINE withEnum #-}
        

instance LAPACK Double where
    geqrf m n pa lda ptau =
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        withEnum lda $ \plda ->
            checkInfo =<< callWithWork (dgeqrf pm pn pa plda ptau)

    gelqf m n pa lda ptau =
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        withEnum lda $ \plda ->
            checkInfo =<< callWithWork (dgelqf pm pn pa plda ptau)

    unmqr side trans m n k pa lda ptau pc ldc =
        withSide side $ \pside ->
        withTrans trans $ \ptrans ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        withEnum lda $ \plda ->
        withEnum ldc $ \pldc ->
            checkInfo =<< callWithWork (dormqr pside ptrans pm pn pk pa plda ptau pc pldc)

    unmlq side trans m n k pa lda ptau pc ldc =
        withSide side $ \pside ->
        withTrans trans $ \ptrans ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        withEnum lda $ \plda ->
        withEnum ldc $ \pldc ->
            checkInfo =<< callWithWork (dormlq pside ptrans pm pn pk pa plda ptau pc pldc)

    larfg n palpha px incx = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
        alloca $ \ptau -> do
            dlarfg pn palpha px pincx ptau
            peek ptau

    potrf uplo n pa lda =
        withUplo uplo $ \puplo ->
        withEnum n $ \pn ->
        withEnum lda $ \plda ->
        alloca $ \pinfo -> do
            dpotrf puplo pn pa plda pinfo
            info <- fromEnum `fmap` peek pinfo
            assert (info >= 0) $ return info            

    potrs uplo n nrhs pa lda pb ldb =
        withUplo uplo $ \puplo ->
        withEnum n $ \pn ->
        withEnum nrhs $ \pnrhs ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
        alloca $ \pinfo -> do
            dpotrs puplo pn pnrhs pa plda pb pldb pinfo
            checkInfo =<< peek pinfo


instance LAPACK (Complex Double) where
    geqrf m n pa lda ptau =
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        withEnum lda $ \plda ->
            checkInfo =<< callWithWork (zgeqrf pm pn pa plda ptau)

    gelqf m n pa lda ptau =
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        withEnum lda $ \plda ->
            checkInfo =<< callWithWork (zgelqf pm pn pa plda ptau)

    unmqr side trans m n k pa lda ptau pc ldc =
        withSide side $ \pside ->
        withTrans trans $ \ptrans ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        withEnum lda $ \plda ->
        withEnum ldc $ \pldc ->
            checkInfo =<< callWithWork (zunmqr pside ptrans pm pn pk pa plda ptau pc pldc)

    unmlq side trans m n k pa lda ptau pc ldc =
        withSide side $ \pside ->
        withTrans trans $ \ptrans ->
        withEnum m $ \pm ->
        withEnum n $ \pn ->
        withEnum k $ \pk ->
        withEnum lda $ \plda ->
        withEnum ldc $ \pldc ->
            checkInfo =<< callWithWork (zunmlq pside ptrans pm pn pk pa plda ptau pc pldc)

    larfg n palpha px incx = 
        withEnum n $ \pn ->
        withEnum incx $ \pincx ->
        alloca $ \ptau -> do
            zlarfg pn palpha px pincx ptau
            peek ptau

    potrf uplo n pa lda =
        withUplo uplo $ \puplo ->
        withEnum n $ \pn ->
        withEnum lda $ \plda ->
        alloca $ \pinfo -> do
            zpotrf puplo pn pa plda pinfo
            info <- fromEnum `fmap` peek pinfo
            assert (info >= 0) $ return info

    potrs uplo n nrhs pa lda pb ldb =
        withUplo uplo $ \puplo ->
        withEnum n $ \pn ->
        withEnum nrhs $ \pnrhs ->
        withEnum lda $ \plda ->
        withEnum ldb $ \pldb ->
        alloca $ \pinfo -> do
            zpotrs puplo pn pnrhs pa plda pb pldb pinfo
            checkInfo =<< peek pinfo
