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
    checkedRow,
    checkedCol,
    checkedDiag,
    checkedSubmatrix,
    checkMatMatOp,
    checkMatVecMult,
    checkMatMatMult,
    diagStart,
    diagLen,    
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

checkedRow ::  (Int,Int) -> (Int -> v) -> Int -> v
checkedRow (m,n) row i 
    | i < 0 || i >= m =
        error $ printf
            "Error in row index.  Tried to get row `%d' in a matrix with shape `(%d,%d)'" i m n
    | otherwise =
        row i

checkedCol :: (Int,Int) -> (Int -> v) -> Int -> v
checkedCol (m,n) col j 
    | j < 0 || j >= n =
        error $ printf
            "Error in column index.  Tried to get column `%d' in a matrix with shape `(%d,%d)'" j m n
    | otherwise =
        col j

checkedDiag :: (Int,Int) -> (Int -> v) ->  Int -> v
checkedDiag (m,n) diag i
    | i < 0 && negate i >= m =
        error $ printf
            "Tried to get sub-diagonal `%d' of a matrix with shape `(%d,%d)'" (negate i) m n
    | i > 0 && i >= n =
        error $ printf
            "Tried to get super-diagonal `%d' of a matrix with shape `(%d,%d)'" i m n        
    | otherwise = 
        diag i

diagStart :: Int -> (Int,Int)
diagStart i
    | i <= 0 =
        (negate i, 0)
    | otherwise =
        (0, i)
        
diagLen :: (Int,Int) -> Int -> Int
diagLen (m,n) i
    | m <= n =
        if i <= 0 
            then max (m + i) 0
            else min (n - i) m
    | otherwise =
        if i > 0
            then max (n - i) 0
            else min (m + i) n

checkedSubmatrix :: (Int,Int) -> ((Int,Int) -> (Int,Int) -> a) -> (Int,Int) -> (Int,Int) -> a
checkedSubmatrix (m,n) sub (i,j) (m',n')
    | or [ i < 0, m' < 0, i + m' > m, 
           j < 0, n' < 0, j + n' > n ] =
        error $ printf ("tried to create submatrix of a `(%d,%d)' matrix " ++
                        " using offset `(%d,%d)' and shape (%d,%d)") m n i j m' n'
    | otherwise =
        sub (i,j) (m',n')


checkMatMatOp :: String -> (Int,Int) -> (Int,Int) -> IO ()
checkMatMatOp name mn1 mn2
    | mn1 /= mn2 =
        ioError $ userError $ printf
            ("%s: x and y have different shapes.  x has shape `%s',"
             ++ " and y has shape `%s'") name (show mn1) (show mn2)
    | otherwise =
        return ()
        
checkMatVecMult :: (Int,Int) -> Int -> IO ()
checkMatVecMult mn n
    | snd mn /= n =
        ioError $ userError $ printf
            ("Tried to multiply a matrix with shape `%s' by a vector of dimension `%d'")
            (show mn) n
    | otherwise =
        return ()
        
checkMatMatMult :: (Int,Int) -> (Int,Int) -> IO ()
checkMatMatMult mk kn
    | snd mk /= fst kn =
        ioError $ userError $ printf
            ("Tried to multiply a matrix with shape `%s' by a matrix with shape `%s'")
            (show mk) (show kn)
    | otherwise =
        return ()

        
