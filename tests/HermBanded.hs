
module HermBanded
    where

import Driver
import Monadic
import qualified Test.QuickCheck.BLAS as Test

import Data.Matrix.Herm
import Data.Matrix.Banded
import Data.Vector.Dense
import Data.Vector.Dense.ST
import Data.Matrix.Dense.ST

import Test.Matrix.Herm.Banded

type V = Vector E
type B = Banded E
type HB = Herm Banded E


prop_herm_apply (HermBandedMV (h :: HB) a x) =
    h <*> x ~== a <*> x

prop_herm_sapply k (HermBandedMV (h :: HB) a x) =
    sapplyVector k h x ~== sapplyVector k a x

prop_herm_herm_apply (HermBandedMV (h :: HB) a x) =
    herm h <*> x ~== h <*> x

prop_doSapplyAddVector alpha beta (HermBandedMV (a :: HB) _ x) = monadicST $ do
    forAllM (Test.vector (numRows a)) $ \y -> do
        y'  <- run $ (unsafeThawVector y :: ST s (STVector s E))
        y'' <- run $ freezeVector y'
        run $ doSApplyAddVector alpha a x beta y'
        assert $ y ~== a <*> (alpha *> x) + (beta *> y'')

prop_herm_applyMatrix (HermBandedMM (h :: HB) a b) =
    h <**> b ~== a <**> b

prop_herm_sapplyMatrix k (HermBandedMM (h :: HB) a b) =
    sapplyMatrix k h b ~== sapplyMatrix k a b

prop_herm_herm_applyMatrix (HermBandedMM (h :: HB) _ b) =
    herm h <**> b ~== h <**> b

prop_doSapplyAddMatrix alpha beta (HermBandedMM (a :: HB) _ b) = monadicST $ do
    forAllM (Test.matrix (numRows a, numCols b)) $ \c -> do
        c'  <- run $ unsafeThawMatrix c
        c'' <- run $ freezeMatrix c'
        run $ doSApplyAddMatrix alpha a b beta c'
        assert $ c ~== a <**> (alpha *> b) + (beta *> c'')


tests_HermBanded =
    [ testProperty "herm apply" prop_herm_apply
    , testProperty "herm sapply" prop_herm_sapply
    , testProperty "herm herm apply" prop_herm_herm_apply
    , testProperty "doSApplyAddVector" prop_doSapplyAddVector    

    , testProperty "herm applyMatrix" prop_herm_applyMatrix
    , testProperty "herm sapplyMatrix" prop_herm_sapplyMatrix
    , testProperty "herm herm applyMatrix" prop_herm_herm_applyMatrix
    , testProperty "doSApplyAddMatrix" prop_doSapplyAddMatrix            
    ]
