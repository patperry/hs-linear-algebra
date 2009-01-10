{-# LANGUAGE ScopedTypeVariables #-}
module LU ( luFactorize ) where

import Data.Elem.BLAS( BLAS3 )

import Control.Monad
import Control.Monad.ST

import Data.Matrix.Dense
import Data.Matrix.Dense.ST
import Data.Matrix.Tri
import Data.Vector.Dense.ST


lu :: (BLAS3 e) => Matrix (n,p) e -> Either Int (Matrix (n,p) e, [Int])
lu (a :: Matrix (n,p) e) = runST $ do
    ma <- thawMatrix a :: ST s (STMatrix s (n,p) e)
    luFactorize ma >>=
        either (return . Left) (\pivots -> do
            a' <- unsafeFreezeMatrix ma
            return $ Right (a',pivots)
            )


{-
 - Recursive LU factorization with row pivoting.  Takes a matrix
 - A and factors it as P A = L U, where P is a permutation matrix, 
 - L is a lower triangular matrix with ones along the diagonal, and 
 - U is an upper triangular matrix.  On successful return, the values of
 - L and U are stored in A, and a list of the row swaps are returned.
 - On failure, the index of the failing column is returned.
 -}      
luFactorize :: (WriteMatrix a e m) => a (n,p) e -> m (Either Int [Int])
luFactorize a
    | mn > 1 =
        let nleft = mn `div` 2
            (a_1, a_2) = splitColsAt nleft a
            (a11, a21) = splitRowsAt nleft a_1
            (a12, a22) = splitRowsAt nleft a_2
        in luFactorize a_1 >>=
               either (return . Left) (\pivots -> do
                   zipWithM_ (swapRows a_2) [0..] pivots
                   doSolveMat_ (lowerU a11) a12
                   doSApplyAddMat (-1) a21 a12 1 a22
                   luFactorize a22 >>=
                       either (return . Left . (nleft+)) (\pivots' -> do
                           zipWithM_ (swapRows a21) [0..] pivots'
                           return $ Right (pivots ++ map (nleft+) pivots')
                       )
               )
    | mn == 1 = 
        let x = colView a 0
        in getWhichMaxAbs x >>= \(i,e) ->
            if (e /= 0) 
                then do
                    scaleBy (1/e) x
                    readElem x 0 >>= writeElem x i
                    writeElem x 0 e
                    return $ Right [i]
                else
                    return $ Left 0
    | otherwise =
        return $ Right []
  where
    (m,n) = shape a
    mn    = min m n
