module LU ( luFactor, luFactorM ) where

import Control.Monad( zipWithM_ )
import Control.Monad.ST( ST, runST )

import Numeric.LinearAlgebra
import qualified Numeric.LinearAlgebra.Matrix as M
import qualified Numeric.LinearAlgebra.Vector as V


luFactor :: (BLAS3 e) => Matrix e -> Either Int (Matrix e, [Int])
luFactor a = runST $ do
    ma <- M.newCopy a
    luFactorM ma >>=
        either (return . Left) (\pivots -> do
            a' <- M.unsafeFreeze ma
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
luFactorM :: (BLAS3 e) => STMatrix s e -> ST s (Either Int [Int])
luFactorM a = do
    (m,n) <- M.getDim a
    let mn = min m n
        nleft = mn `div` 2
    
    case undefined of
        _ | mn > 1 ->     
            M.withSplitColsAtM nleft a $ \a_1 a_2 ->
            M.withSplitRowsAtM nleft a_1 $ \a11 a21 ->
            M.withSplitRowsAtM nleft a_2 $ \a12 a22 ->
            luFactorM a_1 >>=
                either (return . Left) (\pivots -> do
                    zipWithM_ (M.swapRows a_2) [ 0.. ] pivots
                    M.triSolvMatrixM_ LeftSide NoTrans (Tri Lower Unit a11) a12
                    M.addMulMatrixWithScalesM_ (-1) NoTrans a21 NoTrans a12 1 a22
                    luFactorM a22 >>=
                        either (return . Left . (nleft+)) (\pivots' -> do
                            zipWithM_ (M.swapRows a21) [ 0.. ] pivots'
                            return $ Right (pivots ++ map (nleft+) pivots')
                        )
                )
        
        _ | mn == 1 ->
            M.withColM a 0 $ \x ->
                V.getWhichMaxAbs x >>= \(i,e) ->
                    if (e /= 0) 
                        then do
                            V.scaleByM_ (1/e) x
                            V.read x 0 >>= V.write x i
                            V.write x 0 e
                            return $ Right [i]
                        else
                            return $ Left 0

        _ | otherwise ->
            return $ Right []
