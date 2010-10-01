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
import qualified Numeric.LinearAlgebra.Matrix.Packed as P
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
    [ testPropertyD "sumVector" prop_sumVector
    , testPropertyD "weightedSumVector" prop_weightedSumVector
    , testPropertyD "meanVector" prop_meanVector
    , testPropertyD "weightedMeanVector" prop_weightedMeanVector
    , testPropertyD "weightedMeanVector (equal weights)" prop_weightedMeanVector_eqw
    , testPropertyD "covMatrix" prop_covMatrix
    , testPropertyD "covMatrixWithMean" prop_covMatrixWithMean
    , testPropertyD "weightedCovMatrix (equal weights)" prop_weightedCovMatrix_eqw
    , testPropertyD "weightedCovMatrix" prop_weightedCovMatrix
    , testPropertyD "weightedCovMatrixWithMean" prop_weightedCovMatrixWithMean
    , testPropertyD "covPacked" prop_covPacked
    , testPropertyD "covPackedWithMean" prop_covPackedWithMean
    , testPropertyD "weightedCovPacked (equal weights)" prop_weightedCovPacked_eqw
    , testPropertyD "weightedCovPacked" prop_weightedCovPacked
    , testPropertyD "weightedCovPackedWithMean" prop_weightedCovPackedWithMean
    ]


prop_sumVector t (VectorList p xs) =
    sumVector p xs ~== foldl' V.add (V.constant p 0) xs
  where
    _ = typed t (head xs)

prop_weightedSumVector t (WeightedVectorList p wxs) =
    weightedSumVector p wxs
            ~== sumVector p (map (uncurry V.scale) wxs)
  where
    n = length wxs
    _ = typed t (snd $ head wxs)

prop_meanVector t (VectorList p xs) =
    meanVector p xs ~== V.scale (1/n) (sumVector p xs)
  where
    n = fromIntegral $ max (length xs) 1
    _ = typed t (head xs)

prop_weightedMeanVector_eqw t (VectorList p xs) = let
    wxs = zip (repeat 1) xs
    in weightedMeanVector p wxs ~== meanVector p xs
  where
    _ = typed t (head xs)
    
prop_weightedMeanVector t (WeightedVectorList p wxs) = let
        w_sum = (sum . fst . unzip) wxs
        in weightedMeanVector p wxs
            ~== if w_sum == 0 then V.constant p 0
                              else V.scale (1/w_sum) (weightedSumVector p wxs)
  where
    n = length wxs
    _ = typed t (snd $ head wxs)

prop_covMatrix t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = meanVector p xs
        ys = [ V.sub x xbar | x <- xs ]
        scale = case method of { UnbiasedCov -> 1/(n-1) ; MLCov -> 1/n }
        cov' = foldl' (flip $ \y -> M.rank1Update scale y y)
                      (M.constant (p,p) 0)
                      ys
        cov = covMatrix p method xs

        in mulHermMatrixVector cov z ~== M.mulVector NoTrans cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_covMatrixWithMean t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = meanVector p xs
        cov' = covMatrixWithMean xbar method xs
        cov = covMatrix p method xs
        in mulHermMatrixVector cov z ~== mulHermMatrixVector cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_weightedCovMatrix_eqw t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        wxs = zip (repeat 1) xs
        cov = weightedCovMatrix p method wxs
        cov' = covMatrix p method xs
        in mulHermMatrixVector cov z ~== mulHermMatrixVector cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_weightedCovMatrix t (WeightedVectorList p wxs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        (ws,xs) = unzip wxs
        
        w_sum = sum ws
        ws' = [ w / w_sum | w <- ws ]
        w2_sum = sum [ w*w | w <- ws' ]
        scale = case method of { UnbiasedCov -> 1/(1-w2_sum) ; MLCov -> 1 }

        xbar = weightedMeanVector p wxs
        wys = [ (w, V.sub x xbar) | (w,x) <- zip ws' xs ]
        cov' = if w_sum == 0
                    then M.constant (p,p) 0
                    else foldl' (flip $ \(w,y) -> M.rank1Update (scale*w) y y)
                                (M.constant (p,p) 0)
                                wys

        cov = weightedCovMatrix p method wxs

        in mulHermMatrixVector cov z ~== M.mulVector NoTrans cov' z
  where
    n = fromIntegral $ length wxs
    _ = typed t $ snd $ head wxs

prop_weightedCovMatrixWithMean t (WeightedVectorList p wxs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = weightedMeanVector p wxs
        cov' = weightedCovMatrix p method wxs
        cov = weightedCovMatrixWithMean xbar method wxs
        in mulHermMatrixVector cov z ~== mulHermMatrixVector cov' z
  where
    n = fromIntegral $ length wxs
    _ = typed t $ snd $ head wxs

prop_covPacked t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        cov' = covMatrix p method xs
        cov = covPacked p method xs
        in P.hermMulVector cov z ~== mulHermMatrixVector cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_covPackedWithMean t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = meanVector p xs
        cov = covPackedWithMean xbar method xs
        cov' = covPacked p method xs
        in P.hermMulVector cov z ~== P.hermMulVector cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_weightedCovPacked_eqw t (VectorList p xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        wxs = zip (repeat 1) xs
        cov = weightedCovPacked p method wxs
        cov' = covPacked p method xs
        in P.hermMulVector cov z ~== P.hermMulVector cov' z
  where
    n = fromIntegral $ length xs
    _ = typed t $ head xs

prop_weightedCovPacked t (WeightedVectorList p wxs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        cov' = weightedCovMatrix p method wxs
        cov = weightedCovPacked p method wxs
        in P.hermMulVector cov z ~== mulHermMatrixVector cov' z
  where
    n = fromIntegral $ length wxs
    _ = typed t $ snd $ head wxs

prop_weightedCovPackedWithMean t (WeightedVectorList p wxs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = weightedMeanVector p wxs
        cov' = weightedCovPacked p method wxs
        cov = weightedCovPackedWithMean xbar method wxs
        in P.hermMulVector cov z ~== P.hermMulVector cov' z
  where
    n = fromIntegral $ length wxs
    _ = typed t $ snd $ head wxs
