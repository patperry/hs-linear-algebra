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
import Data.Matrix.Tri.Dense  

import Data.AEq
import Numeric.IEEE

import Test.QuickCheck.Complex
import Test.QuickCheck.Matrix.Tri.Dense

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
type TM = Tri (Matrix) (Int,Int) E

instance (Arbitrary e, RealFloat e) => Arbitrary (Complex e) where
    arbitrary   = arbitrary >>= \(TestComplex x) -> return x
    coarbitrary = coarbitrary . TestComplex

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
    in t <*> x ~== y || (any isUndef $ elems x)

prop_tri_ssolve k (TriMatrixSV (t :: TM) y) =
    ssolve k t y ~== t <\> (k *> y)

prop_tri_solveMat (TriMatrixSM (t :: TM) b) =
    let a = t <\\> b
    in t <**> a ~== b || (any isUndef $ elems a)

prop_tri_ssolveMat k (TriMatrixSM (t :: TM) b) =
    ssolveMat k t b ~== t <\\> (k *> b)


properties =
    [ ("tri apply"             , pDet prop_tri_apply)
    , ("tri sapply"            , pDet prop_tri_sapply)
    , ("tri applyMat"          , pDet prop_tri_applyMat)
    , ("tri sapplyMat"         , pDet prop_tri_sapplyMat)

    , ("tri solve"             , pDet prop_tri_solve)
    , ("tri ssolve"            , pDet prop_tri_ssolve)
    , ("tri solveMat"          , pDet prop_tri_solveMat)
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
    pRun n 200 properties
