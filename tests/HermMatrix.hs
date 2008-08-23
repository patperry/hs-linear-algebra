
module HermMatrix( tests_HermMatrix ) where

import Driver

import Data.Matrix.Herm
import Data.Matrix.Dense
import Data.Vector.Dense

import Generators.Matrix.Herm.Dense

type V  = Vector Int E
type M  = Matrix (Int,Int) E
type HM = Herm Matrix (Int,Int) E

prop_herm_apply (HermMatrixMV (h :: HM) a x) =
    h <*> x ~== a <*> x

prop_herm_sapply k (HermMatrixMV (h :: HM) a x) =
    sapply k h x ~== sapply k a x

prop_herm_herm_apply (HermMatrixMV (h :: HM) a x) =
    herm h <*> x ~== h <*> x

prop_herm_applyMat (HermMatrixMM (h :: HM) a b) =
    h <**> b ~== a <**> b

prop_herm_sapplyMat k (HermMatrixMM (h :: HM) a b) =
    sapplyMat k h b ~== sapplyMat k a b

prop_herm_herm_applyMat (HermMatrixMM (h :: HM) _ b) =
    herm h <**> b ~== h <**> b

tests_HermMatrix =
    [ ("herm apply"            , mytest prop_herm_apply)
    , ("herm sapply"           , mytest prop_herm_sapply)
    , ("herm herm apply"       , mytest prop_herm_herm_apply)

    , ("herm applyMat"         , mytest prop_herm_applyMat)
    , ("herm sapplyMat"        , mytest prop_herm_sapplyMat)
    , ("herm herm applyMat"    , mytest prop_herm_herm_applyMat)
    ]
