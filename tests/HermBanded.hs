
module HermBanded( tests_HermBanded ) where

import Driver

import Data.Matrix.Herm
import Data.Matrix.Banded
import Data.Vector.Dense

import Generators.Matrix.Herm.Banded

type V = Vector Int E
type B = Banded (Int,Int) E
type HB = Herm Banded (Int,Int) E


prop_herm_apply (HermBandedMV (h :: HB) a x) =
    h <*> x ~== a <*> x

prop_herm_sapply k (HermBandedMV (h :: HB) a x) =
    sapply k h x ~== sapply k a x

prop_herm_herm_apply (HermBandedMV (h :: HB) a x) =
    herm h <*> x ~== h <*> x

prop_herm_applyMat (HermBandedMM (h :: HB) a b) =
    h <**> b ~== a <**> b

prop_herm_sapplyMat k (HermBandedMM (h :: HB) a b) =
    sapplyMat k h b ~== sapplyMat k a b

prop_herm_herm_applyMat (HermBandedMM (h :: HB) _ b) =
    herm h <**> b ~== h <**> b


tests_HermBanded =
    [ ("herm apply"            , mytest prop_herm_apply)
    , ("herm sapply"           , mytest prop_herm_sapply)
    , ("herm herm apply"       , mytest prop_herm_herm_apply)

    , ("herm applyMat"         , mytest prop_herm_applyMat)
    , ("herm sapplyMat"        , mytest prop_herm_sapplyMat)
    , ("herm herm applyMat"    , mytest prop_herm_herm_applyMat)
    ]
