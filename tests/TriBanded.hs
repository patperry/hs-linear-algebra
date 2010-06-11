
module TriBanded
    where

import Driver
import Monadic
import qualified Test.QuickCheck.BLAS as Test

import Data.Vector.Dense
import Data.Vector.Dense.ST
import Data.Matrix.Dense
import Data.Matrix.Dense.ST
import Data.Matrix.Banded
import Data.Matrix.Tri  


import Test.Matrix.Tri.Banded


type V = Vector E
type M = Matrix E
type B = Banded E
type TB = Tri Banded E


prop_tri_apply (TriBandedMV (t :: TB) a x) =
    t <*> x ~== a <*> x
        
prop_herm_tri_apply (TriBandedMV (t :: TB) a x) =
    herm t <*> x ~== herm a <*> x

prop_scale_tri_apply k (TriBandedMV (t :: TB) a x) =
    sapplyVector k t x ~== sapplyVector k a x

prop_scale_herm_tri_apply k (TriBandedMV (t :: TB) a x) =
    sapplyVector k (herm t) x ~== sapplyVector k (herm a) x

prop_doSapplyAddVector alpha beta (TriBandedMV (a :: TB) _ x) = monadicST $ do
    forAllM (Test.vector (numRows a)) $ \y -> do
        y'  <- run $ (unsafeThawVector y :: ST s (STVector s E))
        y'' <- run $ freezeVector y'
        run $ doSApplyAddVector alpha a x beta y'
        assert $ y ~== a <*> (alpha *> x) + (beta *> y'')

prop_tri_compose (TriBandedMM (t :: TB) a b) =
    t <**> b ~== a <**> b

prop_herm_tri_compose (TriBandedMM (t :: TB) a b) =
    herm t <**> b ~== herm a <**> b

prop_scale_tri_compose k (TriBandedMM (t :: TB) a b) =
    sapplyMatrix k t b ~== sapplyMatrix k a b

prop_scale_herm_tri_compose k (TriBandedMM (t :: TB) a b) =
    sapplyMatrix k (herm t) b ~== sapplyMatrix k (herm a) b

prop_doSapplyAddMatrix alpha beta (TriBandedMM (a :: TB) _ b) = monadicST $ do
    forAllM (Test.matrix (numRows a, numCols b)) $ \c -> do
        c'  <- run $ unsafeThawMatrix c
        c'' <- run $ freezeMatrix c'
        run $ doSApplyAddMatrix alpha a b beta c'
        assert $ c ~== a <**> (alpha *> b) + (beta *> c'')

prop_tri_solve (TriBandedSV (t :: TB) y) =
    let x = t <\> y
    in t <*> x ~== y

prop_tri_ssolve k (TriBandedSV (t :: TB) y) =
    ssolveVector k t y ~== t <\> (k*>y)

prop_doSSolveVector alpha (TriBandedSV (a :: TB) x) = monadicST $ do
    forAllM (Test.vector (numCols a)) $ \y -> do
        y'  <- run $ (unsafeThawVector y :: ST s (STVector s E))
        run $ doSSolveVector alpha a x y'
        assert $ y ~== a <\> (alpha *> x)

prop_tri_solveMatrix (TriBandedSM (t :: TB) c) =
    let b = t <\\> c
    in t <**> b ~== c

prop_tri_ssolveMatrix k (TriBandedSM (t :: TB) c) =
    ssolveMatrix k t c ~== t <\\> (k*>c)

prop_doSSolveMatrix alpha (TriBandedSM (a :: TB) b) = monadicST $ do
    forAllM (Test.matrix (numCols a, numCols b)) $ \c -> do
        c'  <- run $ unsafeThawMatrix c
        run $ doSSolveMatrix alpha a b c'
        assert $ c ~== a <\\> (alpha *> b)


tests_TriBanded =
    [ testProperty "tri apply" prop_tri_apply
    , testProperty "scale tri apply" prop_scale_tri_apply
    , testProperty "herm tri apply" prop_herm_tri_apply
    , testProperty "scale herm tri apply" prop_scale_herm_tri_apply
    , testProperty "tri doSApplyAddVector" prop_doSapplyAddVector 
        
    , testProperty "tri compose" prop_tri_compose
    , testProperty "scale tri compose" prop_scale_tri_compose
    , testProperty "herm tri compose" prop_herm_tri_compose
    , testProperty "scale herm tri compose" prop_scale_herm_tri_compose
    , testProperty "tri doSApplyAddMatrix" prop_doSapplyAddMatrix            
    
    , testProperty "tri solve" prop_tri_solve
    , testProperty "tri ssolve" prop_tri_ssolve
    , testProperty "tri doSSolveVector" prop_doSSolveVector
    
    , testProperty "tri solveMatrix" prop_tri_solveMatrix
    , testProperty "tri ssolveMatrix" prop_tri_ssolveMatrix
    , testProperty "tri doSSolveMatrix" prop_doSSolveMatrix   
    ]
