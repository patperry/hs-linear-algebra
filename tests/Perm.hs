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
import Data.Matrix.Perm

import Data.Permutation ( Permutation, permutation )
import qualified Data.Permutation as P

import Data.AEq
import Numeric.IEEE

import Test.QuickCheck.Complex
import Test.QuickCheck.Matrix.Perm

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
type P = Perm (Int,Int) E

instance (Arbitrary e, RealFloat e) => Arbitrary (Complex e) where
    arbitrary   = arbitrary >>= \(TestComplex x) -> return x
    coarbitrary = coarbitrary . TestComplex

prop_perm_herm (TestPerm (p :: P)) =
    toPermutation (herm p) == P.inverse (toPermutation p)

prop_perm_apply_basis (PermMBasis (p :: P) i) =
    n > 0 ==> p <*> (basis n i) === basis n (P.apply (toPermutation p) i)
  where
    n = numCols p

prop_perm_herm_apply (PermMV (p :: P) x) =
    p <*> herm p <*> x === x

prop_herm_perm_apply (PermMV (p :: P) x) =
    herm p <*> p <*> x === x

prop_perm_solve (PermMV (p :: P) x) =
    p <\> x === herm p <*> x

prop_perm_applyMat_cols (PermMM (p :: P) a) =
    cols (p <**> a) === map (p <*>) (cols a)

prop_perm_herm_applyMat (PermMM (p :: P) a) =
    p <**> herm p <**> a === a

prop_herm_perm_applyMat (PermMM (p :: P) a) =
    herm p <**> p <**> a === a

prop_perm_solveMat_cols (PermMM (p :: P) a) =
    cols (p <\\> a) === map (p <\>) (cols a)

prop_perm_solveMat (PermMM (p :: P) a) =
    p <\\> a === herm p <**> a
    

properties =
    [ ("perm herm"             , pDet prop_perm_herm)
    , ("perm apply basis"      , pDet prop_perm_apply_basis)
    , ("perm herm apply"       , pDet prop_perm_herm_apply)
    , ("herm perm apply"       , pDet prop_herm_perm_apply)
    , ("perm solve"            , pDet prop_perm_solve)
    , ("perm applyMat cols"    , pDet prop_perm_applyMat_cols)
    , ("perm herm applyMat"    , pDet prop_perm_herm_applyMat)
    , ("herm perm applyMat"    , pDet prop_herm_perm_applyMat)
    , ("perm solveMat cols"    , pDet prop_perm_solveMat_cols)
    , ("perm solveMat"         , pDet prop_perm_solveMat)
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
