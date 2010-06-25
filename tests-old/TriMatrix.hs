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
import Monadic
import Test.Matrix.Tri.Dense
import qualified Test.QuickCheck.BLAS as Test

import Data.Vector.Dense
import Data.Vector.Dense.ST
import Data.Matrix.Dense
import Data.Matrix.Dense.ST
import Data.Matrix.Tri

type V = Vector E
type M = Matrix E
type TM = Tri Matrix E

prop_tri_col (Index2 (m,n) (_,j)) =
    forAll (Test.triMatrix (m,n)) $ \(a :: TM) ->
        col a j ~== a <*> basisVector n j

prop_tri_row (Index2 (m,n) (i,_)) =
    forAll (Test.triMatrix (m,n)) $ \(a :: TM) ->
        row a i ~== conj (herm a <*> basisVector m i)

prop_tri_apply (TriMatrixMV (t :: TM) a x) =
    t <*> x ~== a <*> x

prop_tri_sapply k (TriMatrixMV (t :: TM) a x) =
    sapplyVector k t x ~== sapplyVector k a x

prop_doSapplyAddVector alpha beta (TriMatrixMV (a :: TM) _ x) = monadicST $ do
    forAllM (Test.vector (numRows a)) $ \y -> do
        y'  <- run $ (unsafeThawVector y :: ST s (STVector s E))
        y'' <- run $ freezeVector y'
        run $ doSApplyAddVector alpha a x beta y'
        assert $ y ~== a <*> (alpha *> x) + (beta *> y'')

prop_doSSolveVector alpha (TriMatrixSV (a :: TM) x) = monadicST $ do
    forAllM (Test.vector (numCols a)) $ \y -> do
        y'  <- run $ (unsafeThawVector y :: ST s (STVector s E))
        run $ doSSolveVector alpha a x y'
        assert $ y ~== a <\> (alpha *> x)

prop_tri_applyMatrix (TriMatrixMM (t :: TM) a b) =
    t <**> b ~== a <**> b

prop_tri_sapplyMatrix k (TriMatrixMM (t :: TM) a b) =
    sapplyMatrix k t b ~== sapplyMatrix k a b

prop_doSapplyAddMatrix alpha beta (TriMatrixMM (a :: TM) _ b) = monadicST $ do
    forAllM (Test.matrix (numRows a, numCols b)) $ \c -> do
        c'  <- run $ unsafeThawMatrix c
        c'' <- run $ freezeMatrix c'
        run $ doSApplyAddMatrix alpha a b beta c'
        assert $ c ~== a <**> (alpha *> b) + (beta *> c'')

prop_tri_solve (TriMatrixSV (t :: TM) y) =
    let x = t <\> y
    in t <*> x ~== y

prop_tri_ssolve k (TriMatrixSV (t :: TM) y) =
    ssolveVector k t y ~== t <\> (k *> y)

prop_tri_solveMatrix (TriMatrixSM (t :: TM) b) =
    let a = t <\\> b
    in t <**> a ~== b

prop_tri_ssolveMatrix k (TriMatrixSM (t :: TM) b) =
    ssolveMatrix k t b ~== t <\\> (k *> b)

prop_doSSolveMatrix alpha (TriMatrixSM (a :: TM) b) = monadicST $ do
    forAllM (Test.matrix (numCols a, numCols b)) $ \c -> do
        c'  <- run $ unsafeThawMatrix c
        run $ doSSolveMatrix alpha a b c'
        assert $ c ~== a <\\> (alpha *> b)

tests_TriMatrix =
    [ testProperty "tri col" prop_tri_col
    , testProperty "tri row" prop_tri_row
    , testProperty "tri apply" prop_tri_apply
    , testProperty "tri sapply" prop_tri_sapply
    , testProperty "tri doSApplyAddVector" prop_doSapplyAddVector    
    , testProperty "tri applyMatrix" prop_tri_applyMatrix
    , testProperty "tri sapplyMatrix" prop_tri_sapplyMatrix
    , testProperty "tri doSApplyAddMatrix" prop_doSapplyAddMatrix        

    , testProperty "tri solve" prop_tri_solve
    , testProperty "tri ssolve" prop_tri_ssolve
    , testProperty "tri doSSolveVector" prop_doSSolveVector        
    , testProperty "tri solveMatrix" prop_tri_solveMatrix
    , testProperty "tri ssolveMatrix" prop_tri_ssolveMatrix
    , testProperty "tri doSSolveMatrix" prop_doSSolveMatrix            
    ]
