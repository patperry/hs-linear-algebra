module Statistics (
    tests_Statistics
    ) where

import Control.Monad( replicateM )
import Data.AEq
import Data.List( foldl' )
import Debug.Trace
import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Numeric.LinearAlgebra
import qualified Numeric.LinearAlgebra.Packed as P
import qualified Numeric.LinearAlgebra.Matrix as M
import qualified Numeric.LinearAlgebra.Vector as V

import Test.QuickCheck.LinearAlgebra( VectorList(..), WeightedVectorList(..) )
import qualified Test.QuickCheck.LinearAlgebra as Test

import Typed

testAEq actual expected =
    if actual ~== expected
        then True
        else trace ("expected: " ++ show expected ++ "\nactual: " ++ show actual)
                   False



tests_Statistics = testGroup "Statistics"
    [ testPropertyD "sum" prop_sum
    , testPropertyD "weightedSum" prop_weightedSum
    , testPropertyD "mean" prop_mean
    , testPropertyD "weightedMean" prop_weightedMean
    , testPropertyD "weightedMean (equal weights)" prop_weightedMean_eqw
    , testPropertyD "Matrix cov" prop_Matrix_cov
    , testPropertyD "Matrix covWithMean" prop_Matrix_covWithMean
    , testPropertyD "Matrix weightedCov (equal weights)" prop_Matrix_weightedCov_eqw
    , testPropertyD "Matrix weightedCov" prop_Matrix_weightedCov
    , testPropertyD "Matrix weightedCovWithMean" prop_Matrix_weightedCovWithMean
    , testPropertyD "Packed cov" prop_Packed_cov
    , testPropertyD "Packed covWithMean" prop_Packed_covWithMean
    , testPropertyD "Packed weightedCov (equal weights)" prop_Packed_weightedCov_eqw
    , testPropertyD "Packed weightedCov" prop_Packed_weightedCov
    , testPropertyD "Packed weightedCovWithMean" prop_Packed_weightedCovWithMean
    ]


prop_sum t (VectorList p xs) =
    V.sum p xs ~== foldl' V.add (V.constant p 0) xs
  where
    _ = typed t (head xs)

prop_weightedSum t (WeightedVectorList p wxs) =
    V.weightedSum p wxs
            ~== V.sum p [ V.scale x w | (w,x) <- wxs ]
  where
    n = length wxs
    _ = typed t (snd $ head wxs)

prop_mean t (VectorList p xs) =
    V.mean p xs ~== V.scale (V.sum p xs) (1/n)
  where
    n = fromIntegral $ max (length xs) 1
    _ = typed t (head xs)

prop_weightedMean_eqw t (VectorList p xs) = let
    wxs = zip (repeat 1) xs
    in V.weightedMean p wxs ~== V.mean p xs
  where
    _ = typed t (head xs)
    
prop_weightedMean t (WeightedVectorList p wxs) = let
        w_sum = (sum . fst . unzip) wxs
        in V.weightedMean p wxs
            ~== if w_sum == 0 then V.constant p 0
                              else V.scale (V.weightedSum p wxs) (1/w_sum)
  where
    n = length wxs
    _ = typed t (snd $ head wxs)

prop_Matrix_cov t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = V.mean p xs
        ys = [ V.sub x xbar | x <- xs ]
        scale = case method of { UnbiasedCov -> 1/(n-1) ; MLCov -> 1/n }
        cov' = foldl' (flip $ \y -> M.rank1Update scale y y)
                      (M.constant (p,p) 0)
                      ys
        cov = M.cov p method xs

        in M.hermMulVector cov z ~== M.mulVector NoTrans cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_Matrix_covWithMean t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = V.mean p xs
        cov' = M.covWithMean xbar method xs
        cov = M.cov p method xs
        in M.hermMulVector cov z ~== M.hermMulVector cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_Matrix_weightedCov_eqw t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        wxs = zip (repeat 1) xs
        cov = M.weightedCov p method wxs
        cov' = M.cov p method xs
        in M.hermMulVector cov z ~== M.hermMulVector cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_Matrix_weightedCov t (WeightedVectorList p wxs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        (ws,xs) = unzip wxs
        
        w_sum = sum ws
        ws' = [ w / w_sum | w <- ws ]
        w2_sum = sum [ w*w | w <- ws' ]
        scale = case method of { UnbiasedCov -> 1/(1-w2_sum) ; MLCov -> 1 }

        xbar = V.weightedMean p wxs
        wys = [ (w, V.sub x xbar) | (w,x) <- zip ws' xs ]
        cov' = if w_sum == 0
                    then M.constant (p,p) 0
                    else foldl' (flip $ \(w,y) -> M.rank1Update (scale*w) y y)
                                (M.constant (p,p) 0)
                                wys

        cov = M.weightedCov p method wxs

        in M.hermMulVector cov z ~== M.mulVector NoTrans cov' z
  where
    n = fromIntegral $ length wxs
    _ = typed t $ snd $ head wxs

prop_Matrix_weightedCovWithMean t (WeightedVectorList p wxs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = V.weightedMean p wxs
        cov' = M.weightedCov p method wxs
        cov = M.weightedCovWithMean xbar method wxs
        in M.hermMulVector cov z ~== M.hermMulVector cov' z
  where
    n = fromIntegral $ length wxs
    _ = typed t $ snd $ head wxs

prop_Packed_cov t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        cov' = M.cov p method xs
        cov = P.cov p method xs
        in P.hermMulVector cov z ~== M.hermMulVector cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_Packed_covWithMean t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = V.mean p xs
        cov = P.covWithMean xbar method xs
        cov' = P.cov p method xs
        in P.hermMulVector cov z ~== P.hermMulVector cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_Packed_weightedCov_eqw t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        wxs = zip (repeat 1) xs
        cov = P.weightedCov p method wxs
        cov' = P.cov p method xs
        in P.hermMulVector cov z ~== P.hermMulVector cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_Packed_weightedCov t (WeightedVectorList p wxs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        cov' = M.weightedCov p method wxs
        cov = P.weightedCov p method wxs
        in P.hermMulVector cov z ~== M.hermMulVector cov' z
  where
    n = fromIntegral $ length wxs
    _ = typed t $ snd $ head wxs

prop_Packed_weightedCovWithMean t (WeightedVectorList p wxs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = V.weightedMean p wxs
        cov' = P.weightedCov p method wxs
        cov = P.weightedCovWithMean xbar method wxs
        in P.hermMulVector cov z ~== P.hermMulVector cov' z
  where
    n = fromIntegral $ length wxs
    _ = typed t $ snd $ head wxs
