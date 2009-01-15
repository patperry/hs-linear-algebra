
module TriBanded
    where

import Driver

import Data.Vector.Dense
import Data.Matrix.Dense
import Data.Matrix.Banded
import Data.Matrix.Tri  


import Test.Matrix.Tri.Banded


type V = Vector Int E
type M = Matrix (Int,Int) E
type B = Banded (Int,Int) E
type TB = Tri Banded (Int,Int) E


prop_tri_apply (TriBandedMV (t :: TB) a x) =
    t <*> x ~== a <*> x
        
prop_herm_tri_apply (TriBandedMV (t :: TB) a x) =
    herm t <*> x ~== herm a <*> x

prop_scale_tri_apply k (TriBandedMV (t :: TB) a x) =
    sapplyVector k t x ~== sapplyVector k a x

prop_scale_herm_tri_apply k (TriBandedMV (t :: TB) a x) =
    sapplyVector k (herm t) x ~== sapplyVector k (herm a) x


prop_tri_compose (TriBandedMM (t :: TB) a b) =
    t <**> b ~== a <**> b

prop_herm_tri_compose (TriBandedMM (t :: TB) a b) =
    herm t <**> b ~== herm a <**> b

prop_scale_tri_compose k (TriBandedMM (t :: TB) a b) =
    sapplyMatrix k t b ~== sapplyMatrix k a b

prop_scale_herm_tri_compose k (TriBandedMM (t :: TB) a b) =
    sapplyMatrix k (herm t) b ~== sapplyMatrix k (herm a) b


prop_tri_solve (TriBandedSV (t :: TB) y) =
    let x = t <\> y
    in t <*> x ~== y

prop_tri_ssolve k (TriBandedSV (t :: TB) y) =
    ssolveVector k t y ~== t <\> (k*>y)

prop_tri_solveMatrix (TriBandedSM (t :: TB) c) =
    let b = t <\\> c
    in t <**> b ~== c

prop_tri_ssolveMatrix k (TriBandedSM (t :: TB) c) =
    ssolveMatrix k t c ~== t <\\> (k*>c)

tests_TriBanded =
    [ ("tri apply"             , mytest prop_tri_apply)
    , ("scale tri apply"       , mytest prop_scale_tri_apply)
    , ("herm tri apply"        , mytest prop_herm_tri_apply)
    , ("scale herm tri apply"  , mytest prop_scale_herm_tri_apply)
    
    , ("tri compose"           , mytest prop_tri_compose)
    , ("scale tri compose"     , mytest prop_scale_tri_compose)
    , ("herm tri compose"      , mytest prop_herm_tri_compose)
    , ("scale herm tri compose", mytest prop_scale_herm_tri_compose)
    
    , ("tri solve"             , mytest prop_tri_solve)
    , ("tri solveMatrix"          , mytest prop_tri_solveMatrix)
    , ("tri ssolve"            , mytest prop_tri_ssolve)
    , ("tri ssolveMatrix"         , mytest prop_tri_ssolveMatrix)

    ]
