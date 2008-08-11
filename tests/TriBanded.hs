{-# OPTIONS -fglasgow-exts -fno-excess-precision -cpp #-}
-----------------------------------------------------------------------------
-- |
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

import System.Environment ( getArgs )
import Test.QuickCheck.Parallel hiding ( vector )
import qualified Test.QuickCheck as QC

import Data.Complex ( Complex(..) )

import qualified BLAS.Elem as E
import Data.Vector.Dense
import Data.Matrix.Dense
import Data.Matrix.Banded
import Data.Matrix.Tri.Banded  

import Data.AEq
import Numeric.IEEE

import Test.QuickCheck.Complex
import Test.QuickCheck.Matrix.Tri.Banded

isUndefR x = isNaN x || isInfinite x
isUndefC (x :+ y) = isUndefR x || isUndefR y
        
#ifdef COMPLEX
field = "Complex Double"
type E = Complex Double
isUndef = isUndefC
#else
field = "Double"
type E = Double
isUndef = isUndefR
#endif        

type V = Vector Int E
type M = Matrix (Int,Int) E
type B = Banded (Int,Int) E
type TB = Tri Banded (Int,Int) E

instance (Arbitrary e, RealFloat e) => Arbitrary (Complex e) where
    arbitrary   = arbitrary >>= \(TestComplex x) -> return x
    coarbitrary = coarbitrary . TestComplex

prop_tri_apply (TriBandedMV (t :: TB) a x) =
    t <*> x ~== a <*> x
        
prop_herm_tri_apply (TriBandedMV (t :: TB) a x) =
    herm t <*> x ~== herm a <*> x

prop_scale_tri_apply k (TriBandedMV (t :: TB) a x) =
    sapply k t x ~== sapply k a x

prop_scale_herm_tri_apply k (TriBandedMV (t :: TB) a x) =
    sapply k (herm t) x ~== sapply k (herm a) x


prop_tri_compose (TriBandedMM (t :: TB) a b) =
    t <**> b ~== a <**> b

prop_herm_tri_compose (TriBandedMM (t :: TB) a b) =
    herm t <**> b ~== herm a <**> b

prop_scale_tri_compose k (TriBandedMM (t :: TB) a b) =
    sapplyMat k t b ~== sapplyMat k a b

prop_scale_herm_tri_compose k (TriBandedMM (t :: TB) a b) =
    sapplyMat k (herm t) b ~== sapplyMat k (herm a) b


prop_tri_solve (TriBandedSV (t :: TB) y) =
    let x = t <\> y
    in t <*> x ~== y || (any isUndef $ elems x)

prop_tri_ssolve k (TriBandedSV (t :: TB) y) =
    ssolve k t y ~== t <\> (k*>y)

prop_tri_solveMat (TriBandedSM (t :: TB) c) =
    let b = t <\\> c
    in t <**> b ~== c || (any isUndef $ elems b)

prop_tri_ssolveMat k (TriBandedSM (t :: TB) c) =
    ssolveMat k t c ~== t <\\> (k*>c)

properties =
    [ ("tri apply"             , pDet prop_tri_apply)
    , ("scale tri apply"       , pDet prop_scale_tri_apply)
    , ("herm tri apply"        , pDet prop_herm_tri_apply)
    , ("scale herm tri apply"  , pDet prop_scale_herm_tri_apply)
    
    , ("tri compose"           , pDet prop_tri_compose)
    , ("scale tri compose"     , pDet prop_scale_tri_compose)
    , ("herm tri compose"      , pDet prop_herm_tri_compose)
    , ("scale herm tri compose", pDet prop_scale_herm_tri_compose)
    
    , ("tri solve"             , pDet prop_tri_solve)
    , ("tri solveMat"          , pDet prop_tri_solveMat)
    , ("tri ssolve"            , pDet prop_tri_ssolve)
    , ("tri ssolveMat"         , pDet prop_tri_ssolveMat)

    ]


main = do
    args <- getArgs
    n <- case args of
             (a:_) -> readIO a
             _     -> return 1
    main' n

main' n = do
    putStrLn $ "Running tests for " ++ field
    pRun n 400 properties
