-----------------------------------------------------------------------------
-- |
-- Module     : TriMatrix
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module TriMatrix
    where

import Driver
import Test.Matrix.Tri.Dense
import qualified Test.QuickCheck.BLAS as Test

import Data.Vector.Dense
import Data.Matrix.Dense
import Data.Matrix.Tri

type V = Vector Int E
type M = Matrix (Int,Int) E
type TM = Tri Matrix (Int,Int) E

prop_tri_col (Index2 (m,n) (_,j)) =
    forAll (Test.triMatrix (m,n)) $ \(a :: TM) ->
        col a j ~== a <*> basisVector n j

prop_tri_row (Index2 (m,n) (i,_)) =
    forAll (Test.triMatrix (m,n)) $ \(a :: TM) ->
        row a i ~== conj (herm a <*> basisVector m i)

prop_tri_apply (TriMatrixMV (t :: TM) a x) =
    t <*> x ~== a <*> x

prop_tri_sapply k (TriMatrixMV (t :: TM) a x) =
    sapply k t x ~== sapply k a x

prop_tri_applyMat (TriMatrixMM (t :: TM) a b) =
    t <**> b ~== a <**> b

prop_tri_sapplyMat k (TriMatrixMM (t :: TM) a b) =
    sapplyMat k t b ~== sapplyMat k a b

prop_tri_solve (TriMatrixSV (t :: TM) y) =
    let x = t <\> y
    in t <*> x ~== y

prop_tri_ssolve k (TriMatrixSV (t :: TM) y) =
    ssolve k t y ~== t <\> (k *> y)

prop_tri_solveMat (TriMatrixSM (t :: TM) b) =
    let a = t <\\> b
    in t <**> a ~== b

prop_tri_ssolveMat k (TriMatrixSM (t :: TM) b) =
    ssolveMat k t b ~== t <\\> (k *> b)


tests_TriMatrix =
    [ ("tri col"               , mytest prop_tri_col)
    , ("tri row"               , mytest prop_tri_row)
    , ("tri apply"             , mytest prop_tri_apply)
    , ("tri sapply"            , mytest prop_tri_sapply)
    , ("tri applyMat"          , mytest prop_tri_applyMat)
    , ("tri sapplyMat"         , mytest prop_tri_sapplyMat)

    , ("tri solve"             , mytest prop_tri_solve)
    , ("tri ssolve"            , mytest prop_tri_ssolve)
    , ("tri solveMat"          , mytest prop_tri_solveMat)
    , ("tri ssolveMat"         , mytest prop_tri_ssolveMat)
    ]
