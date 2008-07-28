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

instance (Arbitrary e, RealFloat e) => Arbitrary (Complex e) where
    arbitrary   = arbitrary >>= \(TestComplex x) -> return x
    coarbitrary = coarbitrary . TestComplex

prop_herm_apply (HermBandedMV u h (a :: B) x) =
    case u of
        Lower -> hermL h <*> x ~== a <*> x
        Upper -> hermU h <*> x ~== a <*> x

prop_scale_herm_apply k (HermBandedMV u h (a :: B) x) =
    case u of
        Lower -> (k *> hermL h) <*> x ~== (k *> (a <*> x))
        Upper -> (k *> hermU h) <*> x ~== (k *> (a <*> x))

prop_herm_herm_apply (HermBandedMV u h (a :: B) x) =
    case u of
        Lower -> hermU (herm h) <*> x ~== a <*> x
        Upper -> hermL (herm h) <*> x ~== a <*> x

prop_herm_compose (HermBandedMM u h (a :: B) b) =
    case u of
        Lower -> hermL h <**> b ~== a <**> b
        Upper -> hermU h <**> b ~== a <**> b

prop_scale_herm_compose k (HermBandedMM u h (a :: B) b) =
    case u of
        Lower -> (k *> hermL h) <**> b ~== (k *> (a <**> b))
        Upper -> (k *> hermU h) <**> b ~== (k *> (a <**> b))

prop_herm_herm_compose (HermBandedMM u h (a :: B) b) =
    case u of
        Lower -> hermU (herm h) <**> b ~== a <**> b
        Upper -> hermL (herm h) <**> b ~== a <**> b


properties =
    [ 
      ("herm apply"            , pDet prop_herm_apply)
    , ("scale herm apply"      , pDet prop_scale_herm_apply)
    , ("herm herm apply"       , pDet prop_herm_herm_apply)

    , ("herm compose"          , pDet prop_herm_compose)
    , ("scale herm compose"    , pDet prop_scale_herm_compose)
    , ("herm herm compose"     , pDet prop_herm_herm_compose)
    
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
