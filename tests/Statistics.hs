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
    n = realToFrac $ length xs
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



