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

prop_tri_apply (TriMatrixMV u d (a :: M) x) =
    case (u,d) of
        (Lower, NonUnit) -> (lower  a) <*> x ~== a <*> x
        (Lower, Unit   ) -> (lowerU a) <*> x ~== a <*> x + x
        (Upper, NonUnit) -> (upper  a) <*> x ~== a <*> x
        (Upper, Unit   ) -> (upperU a) <*> x ~== a <*> x + x

prop_herm_tri_apply (TriMatrixMV u d (a :: M) x) =
    case (u,d) of
        (Lower, NonUnit) -> (herm $ lower  a) <*> x ~== herm a <*> x
        (Lower, Unit   ) -> (herm $ lowerU a) <*> x ~== herm a <*> x + x
        (Upper, NonUnit) -> (herm $ upper  a) <*> x ~== herm a <*> x
        (Upper, Unit   ) -> (herm $ upperU a) <*> x ~== herm a <*> x + x

prop_scale_tri_apply k (TriMatrixMV u d (a :: M) x) =
    case (u,d) of
        (Lower, NonUnit) -> (k *> lower  a) <*> x ~== (k *> a) <*> x
        (Lower, Unit   ) -> (k *> lowerU a) <*> x ~== (k *> a) <*> x + (k *> x)
        (Upper, NonUnit) -> (k *> upper  a) <*> x ~== (k *> a) <*> x
        (Upper, Unit   ) -> (k *> upperU a) <*> x ~== (k *> a) <*> x + (k *> x)

prop_scale_herm_tri_apply k (TriMatrixMV u d (a :: M) x) =
    case (u,d) of
        (Lower, NonUnit) -> (k *> (herm $ lower  a)) <*> x ~== (k *> herm a) <*> x
        (Lower, Unit   ) -> (k *> (herm $ lowerU a)) <*> x ~== (k *> herm a) <*> x + (k *> x)
        (Upper, NonUnit) -> (k *> (herm $ upper  a)) <*> x ~== (k *> herm a) <*> x
        (Upper, Unit   ) -> (k *> (herm $ upperU a)) <*> x ~== (k *> herm a) <*> x + (k *> x)

prop_herm_scale_tri_apply k (TriMatrixMV u d (a :: M) x) =
    case (u,d) of
        (Lower, NonUnit) -> (herm $ k *> lower  a) <*> x ~== (herm $ k *> a) <*> x
        (Lower, Unit   ) -> (herm $ k *> lowerU a) <*> x ~== (herm $ k *> a) <*> x + ((E.conj k) *> x)
        (Upper, NonUnit) -> (herm $ k *> upper  a) <*> x ~== (herm $ k *> a) <*> x
        (Upper, Unit   ) -> (herm $ k *> upperU a) <*> x ~== (herm $ k *> a) <*> x + ((E.conj k) *> x)


prop_tri_compose (TriMatrixMM u d (a :: M) b) =
    case (u,d) of
        (Lower, NonUnit) -> (lower  a) <**> b ~== a <**> b
        (Lower, Unit   ) -> (lowerU a) <**> b ~== a <**> b + b
        (Upper, NonUnit) -> (upper  a) <**> b ~== a <**> b
        (Upper, Unit   ) -> (upperU a) <**> b ~== a <**> b + b

prop_herm_tri_compose (TriMatrixMM u d (a :: M) b) =
    case (u,d) of
        (Lower, NonUnit) -> (herm $ lower  a) <**> b ~== herm a <**> b
        (Lower, Unit   ) -> (herm $ lowerU a) <**> b ~== herm a <**> b + b
        (Upper, NonUnit) -> (herm $ upper  a) <**> b ~== herm a <**> b
        (Upper, Unit   ) -> (herm $ upperU a) <**> b ~== herm a <**> b + b

prop_scale_tri_compose k (TriMatrixMM u d (a :: M) b) =
    case (u,d) of
        (Lower, NonUnit) -> (k *> lower  a) <**> b ~== (k *> a) <**> b
        (Lower, Unit   ) -> (k *> lowerU a) <**> b ~== (k *> a) <**> b + (k *> b)
        (Upper, NonUnit) -> (k *> upper  a) <**> b ~== (k *> a) <**> b
        (Upper, Unit   ) -> (k *> upperU a) <**> b ~== (k *> a) <**> b + (k *> b)

prop_scale_herm_tri_compose k (TriMatrixMM u d (a :: M) b) =
    case (u,d) of
        (Lower, NonUnit) -> (k *> (herm $ lower  a)) <**> b ~== (k *> herm a) <**> b
        (Lower, Unit   ) -> (k *> (herm $ lowerU a)) <**> b ~== (k *> herm a) <**> b + (k *> b)
        (Upper, NonUnit) -> (k *> (herm $ upper  a)) <**> b ~== (k *> herm a) <**> b
        (Upper, Unit   ) -> (k *> (herm $ upperU a)) <**> b ~== (k *> herm a) <**> b + (k *> b)

prop_herm_scale_tri_compose k (TriMatrixMM u d (a :: M) b) =
    case (u,d) of
        (Lower, NonUnit) -> (herm $ k *> lower  a) <**> b ~== (herm $ k *> a) <**> b
        (Lower, Unit   ) -> (herm $ k *> lowerU a) <**> b ~== (herm $ k *> a) <**> b + ((E.conj k) *> b)
        (Upper, NonUnit) -> (herm $ k *> upper  a) <**> b ~== (herm $ k *> a) <**> b
        (Upper, Unit   ) -> (herm $ k *> upperU a) <**> b ~== (herm $ k *> a) <**> b + ((E.conj k) *> b)


prop_tri_solve (TriMatrixSV (t :: TM) y) =
    let x = t <\> y
    in t <*> x ~== y || (any isUndef $ elems x)

prop_tri_invCompose (TriMatrixSM (t :: TM) b) =
    let a = t <\\> b
    in t <**> a ~== b || (any isUndef $ elems a)


properties =
    [ ("tri apply"             , pDet prop_tri_apply)
    , ("scale tri apply"       , pDet prop_scale_tri_apply)
    , ("herm tri apply"        , pDet prop_herm_tri_apply)
    , ("scale herm tri apply"  , pDet prop_scale_herm_tri_apply)
    , ("herm scale tri apply"  , pDet prop_herm_scale_tri_apply)
    
    , ("tri compose"           , pDet prop_tri_compose)
    , ("scale tri compose"     , pDet prop_scale_tri_compose)
    , ("herm tri compose"      , pDet prop_herm_tri_compose)
    , ("scale herm tri compose", pDet prop_scale_herm_tri_compose)
    , ("herm scale tri compose", pDet prop_herm_scale_tri_compose)
    
    , ("tri solve"             , pDet prop_tri_solve)
    , ("tri invCompose"        , pDet prop_tri_invCompose)
    
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
