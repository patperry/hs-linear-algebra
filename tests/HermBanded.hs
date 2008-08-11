{-# OPTIONS -fglasgow-exts -fno-excess-precision -cpp #-}
-----------------------------------------------------------------------------
-- |
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

import Data.AEq
import Numeric.IEEE

import System.Environment ( getArgs )
import Test.QuickCheck.Parallel hiding ( vector )
import qualified Test.QuickCheck as QC


import qualified BLAS.Elem as E
import Data.Complex ( Complex(..) )

import Data.Matrix.Herm.Banded
import Data.Matrix.Banded
import Data.Vector.Dense

import Test.QuickCheck.Complex
import Test.QuickCheck.Matrix.Herm.Banded

#ifdef COMPLEX
field = "Complex Double"
type E = Complex Double
#else
field = "Double"
type E = Double
#endif

type V = Vector Int E
type B = Banded (Int,Int) E
type HB = Herm Banded (Int,Int) E

instance (Arbitrary e, RealFloat e) => Arbitrary (Complex e) where
    arbitrary   = arbitrary >>= \(TestComplex x) -> return x
    coarbitrary = coarbitrary . TestComplex

prop_herm_apply (HermBandedMV (h :: HB) a x) =
    h <*> x ~== a <*> x

prop_herm_sapply k (HermBandedMV (h :: HB) a x) =
    sapply k h x ~== sapply k a x

prop_herm_herm_apply (HermBandedMV (h :: HB) a x) =
    herm h <*> x ~== h <*> x

prop_herm_applyMat (HermBandedMM (h :: HB) a b) =
    h <**> b ~== a <**> b

prop_herm_sapplyMat k (HermBandedMM (h :: HB) a b) =
    sapplyMat k h b ~== sapplyMat k a b

prop_herm_herm_applyMat (HermBandedMM (h :: HB) _ b) =
    herm h <**> b ~== h <**> b


properties =
    [ ("herm apply"            , pDet prop_herm_apply)
    , ("herm sapply"           , pDet prop_herm_sapply)
    , ("herm herm apply"       , pDet prop_herm_herm_apply)

    , ("herm applyMat"         , pDet prop_herm_applyMat)
    , ("herm sapplyMat"        , pDet prop_herm_sapplyMat)
    , ("herm herm applyMat"    , pDet prop_herm_herm_applyMat)
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
