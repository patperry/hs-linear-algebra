
module HermMatrix
    where

import Driver

import Data.Matrix.Herm
import Data.Matrix.Dense
import Data.Vector.Dense

import Test.Matrix.Herm.Dense
import qualified Test.QuickCheck.BLAS as Test

type V  = Vector Int E
type M  = Matrix (Int,Int) E
type HM = Herm Matrix (Int,Int) E

prop_herm_col (Index n j) =
    forAll (Test.hermMatrix n) $ \(a :: HM) ->
        col a j ~== a <*> basisVector n j

prop_herm_row (Index m i) =
    forAll (Test.hermMatrix m) $ \(a :: HM) ->
        row a i ~== conj (herm a <*> basisVector m i)

prop_herm_apply (HermMatrixMV (h :: HM) a x) =
    h <*> x ~== a <*> x

prop_herm_sapply k (HermMatrixMV (h :: HM) a x) =
    sapplyVector k h x ~== sapplyVector k a x

prop_herm_herm_apply (HermMatrixMV (h :: HM) a x) =
    herm h <*> x ~== h <*> x

prop_herm_applyMatrix (HermMatrixMM (h :: HM) a b) =
    h <**> b ~== a <**> b

prop_herm_sapplyMatrix k (HermMatrixMM (h :: HM) a b) =
    sapplyMatrix k h b ~== sapplyMatrix k a b

prop_herm_herm_applyMatrix (HermMatrixMM (h :: HM) _ b) =
    herm h <**> b ~== h <**> b

tests_HermMatrix =
    [ ("herm col"              , mytest prop_herm_col)
    , ("herm row"              , mytest prop_herm_row)
    , ("herm apply"            , mytest prop_herm_apply)
    , ("herm sapply"           , mytest prop_herm_sapply)
    , ("herm herm apply"       , mytest prop_herm_herm_apply)

    , ("herm applyMatrix"         , mytest prop_herm_applyMatrix)
    , ("herm sapplyMatrix"        , mytest prop_herm_sapplyMatrix)
    , ("herm herm applyMatrix"    , mytest prop_herm_herm_applyMatrix)
    ]
