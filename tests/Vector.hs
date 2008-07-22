{-# OPTIONS -fglasgow-exts -cpp #-}
-----------------------------------------------------------------------------
-- |
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

import qualified Data.Array as Array
import Data.Ix   ( inRange )
import Data.List ( nub, sortBy )
import Data.Ord  ( comparing )
import System.Environment ( getArgs )
import Test.QuickCheck.Parallel hiding ( vector )
import qualified Test.QuickCheck as QC

import Data.Vector.Dense
import BLAS.Elem ( BLAS1 )
import qualified BLAS.Elem as E
import Data.Complex ( Complex(..) )

import Data.AEq
import Numeric.IEEE

import Test.QuickCheck.Complex
import Test.QuickCheck.Vector
import Test.QuickCheck.Vector.Dense


#ifdef COMPLEX
field = "Complex Double"
type E = Complex Double
#else
field = "Double"
type E = Double
#endif

type V = Vector Index E

instance (Arbitrary e, RealFloat e) => Arbitrary (Complex e) where
    arbitrary   = arbitrary >>= \(TestComplex x) -> return x
    coarbitrary = coarbitrary . TestComplex

instance (Arbitrary e, BLAS1 e) => Arbitrary (Vector n e) where
    arbitrary   = arbitrary >>= \(TestVector x) -> return x
    coarbitrary = coarbitrary . TestVector

assocsEq x ies =
    let is = fst $ unzip ies
    in filter (\(i,e) -> i `elem` is) (assocs x) === sortBy (comparing fst) ies
       && (all (==0) $ map snd $ filter (\(i,e) -> not $ i `elem` is) $ assocs x)

prop_vector_dim (Assocs n ies) =
    dim (vector n ies :: V) == n
prop_vector_assocs (Assocs n ies) =
    (vector n ies :: V) `assocsEq` ies

prop_listVector_dim es =
    let n = length es
    in dim (listVector n es :: V) == n
prop_listVector_assocs es =
    let n = length es
    in (listVector n es :: V) `assocsEq` (zip [0..] es)

prop_zero_dim (Index n) =
    dim (zero n :: V) == n
prop_zero_elems (Index n) =
    elems (zero n :: V) == replicate n 0
    
prop_constant_dim (Index n) e =
    dim (constant n e :: V) == n
prop_constant_elems (Index n) (e :: E) =
    elems (constant n e :: V) === replicate n e

prop_basis_dim (Basis n i) =
    dim (basis n i :: V) == n
prop_basis_elems (Basis n i) =
    elems (basis n i :: V) == (replicate i 0) ++ [1] ++ (replicate (n-i-1) 0)

prop_replace_elems (x :: V) (Assocs _ ies) =
    let ies' = filter (\(i,e) -> i < dim x) ies
        x'   = x // ies'
        n    = dim x
    in and $ zipWith (\(i1,e1) (i2,e2) -> (i1 == i2) && (e1 === e2))
                     (assocs x')
                     (Array.assocs $ (Array.//) (Array.array (0,n-1) $ assocs x) ies')

prop_subvector_dim (SubVector _ (x :: V) o n) =
    dim (subvector x o n) == n
prop_subvector_elems (SubVector _ (x :: V) o n) =
    elems (subvector x o n) === (take n $ drop o $ elems x)
    
prop_subvectorWithStride_dim (SubVector s (x :: V) o n) =
    dim (subvectorWithStride s x o n) == n
    
prop_subvectorWithStride_elems (SubVector s (x :: V) o n) =
    let expected = (map snd $ filter (\(i,_) -> (i - o >= 0) 
                              && ((i - o) `mod` s == 0) 
                              && ((i - o) `div` s < n))
                      (assocs x))
        actual = elems (subvectorWithStride s x o n) 
    in expected === actual
        

prop_dim (x :: V) = 
    dim x == length (elems x)
prop_bounds (x :: V) =
    bounds x == (0, dim x - 1)
prop_at (x :: V) (Index i) = 
    i < dim x ==> (x!i) === ((elems x) !! i)


prop_indices (x :: V) =
    indices x == [0..(dim x - 1)]
prop_assocs (x :: V) = 
    assocs x === zip (indices x) (elems x)

prop_scale_elems k (x :: V) =
    (elems $ k *> x) ~== (map (k*) $ elems x)
prop_conj_elems (x :: V) =
    and $ zipWith (===) (elems (conj x)) (map (E.conj) (elems x))
prop_conj_scale k (x :: V) =
    conj (k *> x) ===  (E.conj k *> (conj x))

prop_to_from_list es =
    toList (fromList es :: V) === es

prop_sumAbs (x :: V) =
    sumAbs x ~== (sum $ map E.norm1 $ elems x)
prop_norm2 (x :: V) =
    norm2 x ~== (sqrt $ sum $ map (^2) $ map E.norm $ elems x)
prop_whichMaxAbs1 (x :: V) =
    (dim x > 0) && all (not . isNaN) (map E.norm1 $ elems x) ==>
        let (i,e) = whichMaxAbs x
        in x ! i === e
prop_whichMaxAbs2 (x :: V) =
    (dim x > 0) && all (not . isNaN) (map E.norm1 $ elems x) ==>
        let a = E.norm1 $ snd $ whichMaxAbs x
        in all (<= a) (map E.norm1 $ elems x)
        
prop_dot_self (x :: V) =
    (sqrt $ x <.> x) ~== (E.fromReal $ norm2 x)
prop_dot_conj (Pair (x :: V) y) =
    (x <.> y) ~== (E.conj $ y <.> x)
prop_dot_scale1 k (Pair (x :: V) y) =
    (x <.> (k *> y)) ~== k * (x <.> y)
prop_dot_scale2 k (Pair (x :: V) y) =
    ((k *> x) <.> y) ~== (E.conj k) * (x <.> y)
prop_dot_linear1 (Triple (x :: V) y z) =
    (x <.> (y + z)) ~== (x <.> y + x <.> z)
prop_dot_linear2 (Triple (x :: V) y z) =
    ((x + y) <.> z) ~== (x <.> z + y <.> z)

prop_shift k (x :: V) =
    shift k x ~== x + constant (dim x) k
prop_scale k (x :: V) =
    scale k x ~== x * constant (dim x) k
prop_invScale k (x :: V) =
    invScale k x ~== x / constant (dim x) k

prop_plus (Pair (x :: V) y) =
    elems (x + y) ~== zipWith (+) (elems x) (elems y)
prop_minus (Pair (x :: V) y) =
    elems (x - y) ~== zipWith (-) (elems x) (elems y)
prop_times (Pair (x :: V) y) =
    elems (x * y) ~== zipWith (*) (elems x) (elems y)
prop_divide (Pair (x :: V) y) =
    elems (x / y) ~== zipWith (/) (elems x) (elems y)

prop_negate (x :: V) =
    negate x ~== (-1) *> x
prop_abs (x :: V) =
    elems (abs x) ~== map abs (elems x)
prop_signum (x :: V) =
    elems (signum x) === map signum (elems x)
prop_recip (x :: V) =
    elems (recip x) ~== (map recip $ elems x)

properties =
    [ ("dim of vector"        , pDet prop_vector_dim)
    , ("assocs of vector"     , pDet prop_vector_assocs)
    
    , ("dim of listVector"    , pDet prop_listVector_dim)
    , ("assocs of listVector" , pDet prop_listVector_assocs)
    
    , ("dim of zero"          , pDet prop_zero_dim)
    , ("elems of zero"        , pDet prop_zero_elems)
    
    , ("dim of constant"      , pDet prop_constant_dim)
    , ("elems of constant"    , pDet prop_constant_elems)

    , ("dim of basis"         , pDet prop_basis_dim)
    , ("elems of basis"       , pDet prop_basis_elems)

    , ("dim of subvector"     , pDet prop_subvector_dim)
    , ("elems of subvector"   , pDet prop_subvector_elems)
    
    , ("dim of subvectorWithStride"  
                              , pDet prop_subvectorWithStride_dim)
    , ("elems of subvectorWithStride"   
                              , pDet prop_subvectorWithStride_elems)
    
    , ("elems of replace"     , pDet prop_replace_elems)
    , ("elems of scale"       , pDet prop_scale_elems)
    , ("elems of conj"        , pDet prop_conj_elems)
    , ("conj/scale"           , pDet prop_conj_scale)
    
    , ("dim"                  , pDet prop_dim)
    , ("bounds"               , pDet prop_bounds)
    , ("at"                   , pDet prop_at)

    , ("indices"              , pDet prop_indices)
    , ("assocs"               , pDet prop_assocs)
    
    , ("to/from list"         , pDet prop_to_from_list)
    
    , ("sumAbs"               , pDet prop_sumAbs)
    , ("norm2"                , pDet prop_norm2)
    , ("whichMaxAbs1"         , pDet prop_whichMaxAbs1)
    , ("whichMaxAbs2"         , pDet prop_whichMaxAbs2)
    
    , ("dot self"             , pDet prop_dot_self)
    , ("dot conj"             , pDet prop_dot_conj)
    , ("dot scale1"           , pDet prop_dot_scale1)
    , ("dot scale2"           , pDet prop_dot_scale2)
    , ("dot linear1"          , pDet prop_dot_linear1)
    , ("dot linear2"          , pDet prop_dot_linear2)
    
    , ("shift"                , pDet prop_shift)
    , ("scale"                , pDet prop_scale)
    , ("invScale"             , pDet prop_invScale)
    
    , ("plus"                 , pDet prop_plus)
    , ("minus"                , pDet prop_minus)
    , ("times"                , pDet prop_times)
    , ("divide"               , pDet prop_divide)
    
    , ("negate"               , pDet prop_negate)
    , ("abs"                  , pDet prop_abs)
    , ("signum"               , pDet prop_signum)
    , ("recip"                , pDet prop_recip)
    ]


main = do
    args <- getArgs
    n <- case args of
             (a:_) -> readIO a
             _     -> return 1
    main' n

main' n = do
    putStrLn $ "Running tests for " ++ field
    pRun n 500 properties
