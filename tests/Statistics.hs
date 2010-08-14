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

import Test.QuickCheck.LinearAlgebra( NonEmptyVectorList(..) )
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
    ]


prop_sumVector t (NonEmptyVectorList xs) =
    sumVector xs ~== foldl' addVector (constantVector p 0) xs
  where
    p = dimVector (head xs)
    _ = typed t (head xs)

prop_weightedSumVector t (NonEmptyVectorList xs) =
    forAll (replicateM n $ fmap abs arbitrary) $ \ws ->
        weightedSumVector (zip ws xs)
            ~== sumVector (zipWith scaleVector ws xs)
  where
    n = length xs
    _ = typed t (head xs)

prop_meanVector t (NonEmptyVectorList xs) =
    meanVector xs ~== scaleVector (1/n) (sumVector xs)
  where
    n = fromIntegral $ length xs
    _ = typed t (head xs)

prop_weightedMeanVector_eqw t (NonEmptyVectorList xs) = let
    wxs = zip (repeat 1) xs
    in weightedMeanVector wxs ~== meanVector xs
  where
    _ = typed t (head xs)
    
prop_weightedMeanVector t (NonEmptyVectorList xs) =
    forAll (replicateM n $ fmap abs arbitrary) $ \ws -> let
        wxs = zip ws xs
        w_sum = sum ws
        in weightedMeanVector wxs
            ~== if w_sum == 0 then constantVector p 0
                              else scaleVector (1/w_sum) (weightedSumVector wxs)
  where
    n = length xs
    p = dimVector (head xs)
    _ = typed t (head xs)

prop_covMatrix t (NonEmptyVectorList xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = meanVector xs
        ys = [ subVector x xbar | x <- xs ]
        scale = case method of { UnbiasedCov -> 1/(n-1) ; MLCov -> 1/n }
        cov' = foldl' (flip $ \y -> rank1UpdateMatrix scale y y)
                      (constantMatrix (p,p) 0)
                      ys
        cov = covMatrix method xs

        in mulHermMatrixVector cov z ~== mulMatrixVector NoTrans cov' z
  where
    n = fromIntegral $ length xs
    p = dimVector $ head xs
    _ = typed t $ head xs

prop_covMatrixWithMean t (NonEmptyVectorList xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        xbar = meanVector xs
        cov' = covMatrixWithMean xbar method xs
        cov = covMatrix method xs
        in mulHermMatrixVector cov z ~== mulHermMatrixVector cov' z
  where
    n = fromIntegral $ length xs
    p = dimVector $ head xs
    _ = typed t $ head xs

prop_weightedCovMatrix_eqw t (NonEmptyVectorList xs) =
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        wxs = zip (repeat 1) xs
        cov = weightedCovMatrix method wxs
        cov' = covMatrix method xs
        in mulHermMatrixVector cov z ~== mulHermMatrixVector cov' z
  where
    n = fromIntegral $ length xs
    p = dimVector $ head xs
    _ = typed t $ head xs

prop_weightedCovMatrix t (NonEmptyVectorList xs) =
    forAll (replicateM n $ fmap abs arbitrary) $ \ws ->
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        wxs = zip ws xs
        
        w_sum = sum ws
        ws' = [ w / w_sum | w <- ws ]
        w2_sum = sum [ w*w | w <- ws' ]
        scale = case method of { UnbiasedCov -> 1/(1-w2_sum) ; MLCov -> 1 }

        xbar = weightedMeanVector wxs
        wys = [ (w, subVector x xbar) | (w,x) <- zip ws' xs ]
        cov' = foldl' (flip $ \(w,y) -> rank1UpdateMatrix (scale*w) y y)
                      (constantMatrix (p,p) 0)
                      wys

        cov = weightedCovMatrix method wxs

        in mulHermMatrixVector cov z ~== mulMatrixVector NoTrans cov' z
  where
    n = fromIntegral $ length xs
    p = dimVector $ head xs
    _ = typed t $ head xs

prop_weightedCovMatrixWithMean t (NonEmptyVectorList xs) =
    forAll (replicateM n $ fmap abs arbitrary) $ \ws ->
    forAll (Test.vector p) $ \z ->
    forAll (elements [ UnbiasedCov, MLCov ]) $ \method -> let
        wxs = zip ws xs
        xbar = weightedMeanVector wxs
        cov' = weightedCovMatrix method wxs
        cov = weightedCovMatrixWithMean xbar method wxs
        in mulHermMatrixVector cov z ~== mulHermMatrixVector cov' z
  where
    n = fromIntegral $ length xs
    p = dimVector $ head xs
    _ = typed t $ head xs
