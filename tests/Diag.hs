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
import Data.Matrix.Diag

import Data.AEq
import Numeric.IEEE

import Debug.Trace

import Test.QuickCheck.Complex
import Test.QuickCheck.Matrix.Diag

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
type D = Diag (Int,Int) E

instance (Arbitrary e, RealFloat e) => Arbitrary (Complex e) where
    arbitrary   = arbitrary >>= \(TestComplex x) -> return x
    coarbitrary = coarbitrary . TestComplex

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
    in d <*> x ~== y || (any isUndef $ elems x)

prop_diag_ssolve k (DiagSV (d :: D) y) =
    ssolve k d y ~== d <\> (k *> y)

prop_diag_solveMat (DiagSM (d :: D) b) =
    let a = d <\\> b
    in d <**> a ~== b || (any isUndef $ elems a)

prop_diag_ssolveMat k (DiagSM (d :: D) b) =
    ssolveMat k d b ~== d <\\> (k *> b)


properties =
    [ ("diag apply"             , pDet prop_diag_apply)
    , ("diag sapply"            , pDet prop_diag_sapply)
    , ("diag applyMat"          , pDet prop_diag_applyMat)
    , ("diag sapplyMat"         , pDet prop_diag_sapplyMat)

    , ("diag solve"             , pDet prop_diag_solve)
    , ("diag ssolve"            , pDet prop_diag_ssolve)
    , ("diag solveMat"          , pDet prop_diag_solveMat)
    , ("diag ssolveMat"         , pDet prop_diag_ssolveMat)
    
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
