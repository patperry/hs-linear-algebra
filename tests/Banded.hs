{-# OPTIONS -fglasgow-exts -fno-excess-precision -cpp #-}
-----------------------------------------------------------------------------
-- |
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--


import qualified Data.Array as Array
import Data.Ix   ( inRange, range )
import Data.List ( nub, sortBy, foldl' )
import Data.Ord  ( comparing )
import System.Environment ( getArgs )
import Test.QuickCheck.Parallel hiding ( vector )
import qualified Test.QuickCheck as QC

import Data.Matrix.Banded

import BLAS.Elem ( BLAS1 )
import qualified BLAS.Elem as E

import Data.Complex ( Complex(..) )
import Data.Vector.Dense

import Data.AEq
import Numeric.IEEE

import Test.QuickCheck.Complex
import Test.QuickCheck.Vector hiding ( Assocs )
import Test.QuickCheck.Vector.Dense hiding ( Pair )
import Test.QuickCheck.Matrix
import Test.QuickCheck.Matrix.Banded
        
import Debug.Trace

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

instance (Arbitrary e, BLAS1 e) => Arbitrary (Vector n e) where
    arbitrary   = arbitrary >>= \(TestVector x) -> return x
    coarbitrary = coarbitrary . TestVector

instance (Arbitrary e, BLAS1 e) => Arbitrary (Banded (m,n) e) where
    arbitrary   = arbitrary >>= \(TestBanded x) -> return x
    coarbitrary = coarbitrary . TestBanded


assocsEq :: (BLAS1 e, AEq e) => Banded (m,n) e -> [((Int,Int), e)] -> Bool
assocsEq x ijes =
    let ijs = fst $ unzip ijes
    in filter (\(ij,e) -> ij `elem` ijs) (sortBy (comparing fst) $ assocs x) === sortBy (comparing fst) ijes
       && (all (==0) $ map snd $ filter (\(ij,e) -> not $ ij `elem` ijs) $ assocs x)

fromAssocs (Assocs mn ijes) =
    let kl = foldl' max 0 $ map (\((i,j),_) -> (i-j)) ijes
        ku = foldl' max 0 $ map (\((i,j),_) -> (j-i)) ijes
    in banded mn (kl,ku) ijes 
    
prop_banded_shape a@(Assocs mn _) =
    shape (fromAssocs a :: B) == mn
    
prop_banded_assocs a@(Assocs mn ijes) =
    (fromAssocs a :: B) `assocsEq` ijes

prop_listsBanded_shape (ListsBanded mn lu ds) =
    shape (listsBanded mn lu ds :: B) == mn
    
prop_listsBanded_toLists (ListsBanded mn lu ds) =
    toLists (listsBanded mn lu ds :: B) === (mn,lu,ds)



prop_replace_elems (a :: B) (Assocs _ ijes) =
    let ijes' = filter (\((i,j),_) -> i < numRows a 
                                      && j < numCols a
                                      && (i-j) <= numLower a
                                      && (j-i) <= numUpper a) ijes
        old = filter (\(ij,_) -> not $ ij `elem` (fst . unzip) ijes') $ assocs a
        expected = sortBy (comparing fst) $ old ++ ijes'
        actual   = sortBy (comparing fst) $ assocs (a // ijes')
    in expected === actual


prop_shape (a :: B) = 
    shape a == (numRows a, numCols a)

prop_bandwidth (a :: B) =
    bandwidth a == ((negate . numLower) a, numUpper a)

prop_size (a :: B) =
    size a == (sum $ map diagLen (range $ bandwidth a))
  where
    (m,n) = shape a
    diagLen i | i <= 0    = min (m+i) n
              | otherwise = min m     (n-i) 

prop_bounds (a :: B) =
    bounds a == ((0,0), (numRows a - 1, numCols a - 1))
    
prop_at (BandedAt (a :: B) ij@(i,j)) =
    (a!ij) === expected
    
  where
    (_,(kl,ku),ds) = toLists a
    d = ds !! (kl + (j-i))
    expected = case undefined of
                   _ | i - j > kl -> 0
                   _ | j - i > ku -> 0
                   _ | otherwise  -> d!!i
    
prop_row_dim (BandedAt (a :: B) (i,_)) =
    dim (row a i) == numCols a
prop_col_dim (BandedAt (a :: B) (_,j)) =
    dim (col a j) == numRows a
prop_rows_len (a :: B) =
    length (rows a) == numRows a
prop_cols_len (a :: B) =
    length (cols a) == numCols a
prop_rows_dims (a :: B) =
    map dim (rows a) == replicate (numRows a) (numCols a)
prop_cols_dims (a :: B) =
    map dim (cols a) == replicate (numCols a) (numRows a)

prop_indices_length (a :: B) =
    length (indices a) == size a
prop_indices_lower (a :: B) =
    all (\(i,j) -> i - j <= numLower a) $ indices a
prop_indices_upper (a :: B) =
    all (\(i,j) -> j - i <= numUpper a) $ indices a

prop_elems_length (a :: B) =
    length (elems a) == size a

prop_assocs (a :: B) = 
    assocs a === zip (indices a) (elems a)
prop_assocs_at (a :: B) =
    all (\(ij,e) -> a!ij === e) $ assocs a

{-
prop_scale_elems (a :: M) k =
    and $ zipWith (~==) (elems (k *> a)) (map (k*) (elems a))
prop_herm_elem (BandedAt (a :: M) (i,j)) =
    (herm a) ! (j,i) == E.conj (a!(i,j))
prop_herm_scale (a :: M) k =
    herm (k *> a) === (E.conj k) *> (herm a)

prop_herm_shape (a :: M) =
    shape (herm a) == (numCols a, numRows a)
prop_herm_rows (a :: M) =
    rows (herm a) === map conj (cols a)
prop_herm_cols (a :: M) = 
    cols (herm a) === map conj (rows a)

prop_herm_herm (a :: M) =
    herm (herm a) === a

prop_diag_herm1 (BandedAt (a :: M) (k,_)) =
    diag a (-k) === conj (diag (herm a) k)
prop_diag_herm2 (BandedAt (a :: M) (_,k)) =
    diag a k === conj (diag (herm a) (-k))

prop_fromRow_shape (x :: V) =
    shape (fromRow x :: M) == (1,dim x)
prop_fromRow_elems (x :: V) =
    elems (fromRow x :: M) === elems x

prop_fromCol_shape (x :: V) =
    shape (fromCol x :: M) == (dim x,1)
prop_fromCol_elems (x :: V) =
    elems (fromCol x :: M) === elems x


prop_apply_basis (BandedAt (a :: M) (_,j)) =
    a <*> (basis (numCols a) j :: V) ~== col a j
prop_apply_herm_basis (BandedAt (a :: M) (i,_)) =
    (herm a) <*> (basis (numRows a) i :: V) ~== conj (row a i)
prop_apply_scale k (MultMV (a :: M) x) =
    a <*> (k *> x) ~== k *> (a <*> x)
prop_apply_linear (MultMVPair (a :: M) x y) =
    a <*> (x + y) ~== a <*> x + a <*> y

prop_compose_id_right (a :: M) =
    let n = numCols a
    in a <**> (identity (n,n) :: M) ~== a
prop_compose_id_left (a :: M) =
    let m = numRows a
    in (identity (m,m) :: M) <**> a ~== a
prop_compose_scale_left (MultMM (a:: M) b) k =
    a <**> (k *> b) ~== k *> (a <**> b)    
prop_compose_scale_right (MultMM (a:: M) b) k =
    (k *> a) <**> b ~== k *> (a <**> b)
prop_compose_linear (MultMMPair (a :: M) b c) =
    a <**> (b + c) ~== a <**> b + a <**> c
prop_compose_herm (MultMM (a :: M) b) =
    herm b <**> herm a ~== herm (a <**> b)
prop_compose_cols (MultMM (a :: M) b) =
    cols (a <**> b) ~== map (a <*> ) (cols b)

prop_scale k (a :: M) =
    k *> a ~== a * constant (shape a) k
prop_invScale k (a :: M) =
    invScale k a ~== a / constant (shape a) k

prop_plus (Pair (a :: M) b) =
    elems (a + b) ~== zipWith (+) (elems a) (elems b)
prop_minus (Pair (a :: M) b) =
    elems (a - b) ~== zipWith (-) (elems a) (elems b)
prop_times (Pair (a :: M) b) =
    elems (a * b) ~== zipWith (*) (elems a) (elems b)
prop_divide (Pair (a :: M) b) =
    elems (a / b) ~== zipWith (/) (elems a) (elems b)

prop_negate (a :: M) =
    negate a ~== (-1) *> a
    
prop_abs (a :: M) =
    elems (abs a) ~== map abs (elems a)
prop_signum (a :: M) =
    elems (signum a) === map signum (elems a)
prop_recip (a :: M) =
    elems (recip a) ~== (map recip $ elems a)

-}

properties =
    [ ("shape of banded"       , pDet prop_banded_shape)
    , ("assocs of banded"      , pDet prop_banded_assocs)
    , ("shape of listsBanded"  , pDet prop_listsBanded_shape)
    , ("listsBanded/toLists"   , pDet prop_listsBanded_toLists)

    , ("elems of replace"      , pDet prop_replace_elems)
    
    , ("numRows/numCols"       , pDet prop_shape)
    , ("numLower/numUpper"     , pDet prop_bandwidth)
    , ("size"                  , pDet prop_size)
    , ("bounds"                , pDet prop_bounds)

    , ("at"                    , pDet prop_at)
    , ("row dim"               , pDet prop_row_dim)
    , ("col dim"               , pDet prop_col_dim)
    , ("rows length"           , pDet prop_rows_len)
    , ("cols length"           , pDet prop_cols_len)
    , ("rows dims"             , pDet prop_rows_dims)
    , ("cols dims"             , pDet prop_cols_dims)

    , ("indices length"        , pDet prop_indices_length)
    , ("indices low bw"        , pDet prop_indices_lower)
    , ("indices up bw"         , pDet prop_indices_upper)

    , ("elems length"          , pDet prop_elems_length)

    , ("assocs"                , pDet prop_assocs)
    , ("assocs/at"             , pDet prop_assocs_at)

{-    
    , ("elems of scale"        , pDet prop_scale_elems)
    , ("elem of herm"          , pDet prop_herm_elem)
    , ("herm/scale"            , pDet prop_herm_scale)
                               
    , ("shape . herm"          , pDet prop_herm_shape)
    , ("rows . herm"           , pDet prop_herm_rows)
    , ("cols . herm"           , pDet prop_herm_cols)
                               
    , ("herm . herm == id"     , pDet prop_herm_herm)
                               
    , ("subdiag . herm"        , pDet prop_diag_herm1)
    , ("superdiag . herm"      , pDet prop_diag_herm2)
                               
    , ("shape . fromRow"       , pDet prop_fromRow_shape)
    , ("elems . fromRow"       , pDet prop_fromRow_elems)
    , ("shape . fromCol"       , pDet prop_fromCol_shape)
    , ("elems . fromCol"       , pDet prop_fromCol_elems)
                               
    , ("apply basis"           , pDet prop_apply_basis)
    , ("apply herm basis"      , pDet prop_apply_herm_basis)
    , ("apply scale"           , pDet prop_apply_scale)
    , ("apply linear"          , pDet prop_apply_linear)
    
    , ("compose id left"       , pDet prop_compose_id_left)
    , ("compose id right"      , pDet prop_compose_id_right)
    , ("compose scale left"    , pDet prop_compose_scale_left)
    , ("compose scale right"   , pDet prop_compose_scale_right)
    , ("compose linear"        , pDet prop_compose_linear)
    , ("compose herm"          , pDet prop_compose_herm)
    , ("compose cols"          , pDet prop_compose_cols)
    
    , ("scale"                 , pDet prop_scale)
    , ("invScale"              , pDet prop_invScale)
    
    , ("plus"                  , pDet prop_plus)
    , ("minus"                 , pDet prop_minus)
    , ("times"                 , pDet prop_times)
    , ("divide"                , pDet prop_divide)
    
    , ("negate"                , pDet prop_negate)
    , ("abs"                   , pDet prop_abs)
    , ("signum"                , pDet prop_signum)
    , ("recip"                 , pDet prop_recip)
-}
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
