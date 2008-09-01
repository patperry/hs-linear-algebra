module Diag( tests_Diag ) where
    
import Driver

import Generators.Matrix.Diag

import qualified BLAS.Elem as E
import Data.Vector.Dense
import Data.Matrix.Dense
import Data.Matrix.Diag

type V = Vector Int E
type M = Matrix (Int,Int) E
type D = Diag Vector (Int,Int) E


prop_diag_apply (DiagMV (d :: D) a x) =
    d <*> x ~== a <*> x

prop_diag_sapply k (DiagMV (d :: D) a x) =
    sapply k d x ~== sapply k a x

prop_diag_applyMat (DiagMM (d :: D) a b) =
    d <**> b ~== a <**> b

prop_diag_sapplyMat k (DiagMM (d :: D) a b) =
    sapplyMat k d b ~== sapplyMat k a b


prop_diag_solve (DiagSV (d :: D) y) =
    let x = d <\> y
    in d <*> x ~== y

prop_diag_ssolve k (DiagSV (d :: D) y) =
    ssolve k d y ~== d <\> (k *> y)

prop_diag_solveMat (DiagSM (d :: D) b) =
    let a = d <\\> b
    in d <**> a ~== b

prop_diag_ssolveMat k (DiagSM (d :: D) b) =
    ssolveMat k d b ~== d <\\> (k *> b)


tests_Diag =
    [ ("diag apply"             , mytest prop_diag_apply)
    , ("diag sapply"            , mytest prop_diag_sapply)
    , ("diag applyMat"          , mytest prop_diag_applyMat)
    , ("diag sapplyMat"         , mytest prop_diag_sapplyMat)

    , ("diag solve"             , mytest prop_diag_solve)
    , ("diag ssolve"            , mytest prop_diag_ssolve)
    , ("diag solveMat"          , mytest prop_diag_solveMat)
    , ("diag ssolveMat"         , mytest prop_diag_ssolveMat)
    
    ]
