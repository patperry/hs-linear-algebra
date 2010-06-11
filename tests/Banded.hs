{-# OPTIONS -fglasgow-exts -fno-excess-precision -cpp #-}
-----------------------------------------------------------------------------
-- |
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Banded( tests_Banded ) where

import Driver
import Monadic
import qualified Test.QuickCheck.BLAS as Test

import Data.Elem.BLAS
import Data.Vector.Dense 
import Data.Vector.Dense.ST
import Data.Matrix.Dense ( Matrix, identityMatrix )
import Data.Matrix.Dense.ST
import Data.Matrix.Banded

import Test.Matrix.Banded hiding ( banded )
        

listsFromBanded :: (BLAS1 e) => Banded e -> ((Int,Int), (Int,Int),[[e]])
listsFromBanded a = ( (m,n)
            , (kl,ku)
            , map paddedDiag [(-kl)..ku]
            )
  where
    (m,n)   = shape a
    (kl,ku) = bandwidths a
    
    padBegin i   = replicate (max (-i) 0)    0
    padEnd   i   = replicate (max (m-n+i) 0) 0
    paddedDiag i = (  padBegin i
                   ++ elems (diagBanded a i) 
                   ++ padEnd i 
                   )
                   
type V = Vector E
type M = Matrix E
type B = Banded E

fromAssocs (Assocs2 mn ijes) =
    let kl = foldl' max 0 $ map (\((i,j),_) -> (i-j)) ijes
        ku = foldl' max 0 $ map (\((i,j),_) -> (j-i)) ijes
    in banded mn (kl,ku) ijes 
    
prop_banded_shape a@(Assocs2 mn _) =
    shape (fromAssocs a :: B) == mn
    
prop_banded_assocs a@(Assocs2 (m,n) ijes) =
    let as    = (assocs (fromAssocs a :: B))
        ijes' = (reverse . (nubBy ((==) `on` fst)) . reverse) ijes
    in (all (==0) $ map snd (as \\ ijes'))
       && (length (as `intersect` ijes') === length ijes')

prop_listsBanded_shape (ListsBanded mn lu ds) =
    shape (listsBanded mn lu ds :: B) == mn
    
prop_listsBanded_listsFromBanded (ListsBanded mn lu ds) =
    listsFromBanded (listsBanded mn lu ds :: B) === (mn,lu,ds)

prop_replace_elems (a :: B) (Assocs2 _ ijes) =
    let ijes' = filter (\((i,j),_) -> i < numRows a 
                                   && j < numCols a
                                   && (i-j) <= numLower a
                                   && (j-i) <= numUpper a) ijes
        old = filter (\(ij,_) -> not $ ij `elem` (fst . unzip) ijes') $ assocs a
        new = (reverse . nubBy ((==) `on` fst) . reverse) ijes'
        expected = sortBy (comparing fst) $ old ++ new
        actual   = sortBy (comparing fst) $ assocs (a // ijes')
    in expected === actual


prop_shape (a :: B) = 
    shape a == (numRows a, numCols a)

prop_bandwidths (a :: B) =
    bandwidths a == (numLower a, numUpper a)

prop_size (a :: B) =
    size a == (sum $ map diagLen (range (-kl,ku)))
  where
    (m,n) = shape a
    (kl,ku) = bandwidths a
    diagLen i | i <= 0    = min (m+i) n
              | otherwise = min m     (n-i) 

prop_bounds (a :: B) =
    bounds a == ((0,0), (numRows a - 1, numCols a - 1))
  
prop_at (BandedAt (a :: B) ij@(i,j)) =
    (a!ij) === (diagBanded a (j-i))!(min i j)
    
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
prop_row_at (BandedAt (a :: B) (i,j)) =
    (row a i)!j === a!(i,j)
prop_col_at (BandedAt (a :: B) (i,j)) =
    (col a j)!i === a!(i,j)

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

prop_scale_elems (a :: B) k =
    (assocs (k *> a)) `assocsEq` (map (second (*k)) (assocs a))
prop_herm_elem (BandedAt (a :: B) (i,j)) =
    (herm a) ! (j,i) === conjugate (a!(i,j))
prop_herm_scale (a :: B) k =
    herm (k *> a) === (conjugate k) *> (herm a)

prop_herm_shape (a :: B) =
    shape (herm a) == (numCols a, numRows a)
prop_herm_rows (a :: B) =
    rows (herm a) === map conj (cols a)
prop_herm_cols (a :: B) = 
    cols (herm a) === map conj (rows a)

prop_herm_herm (a :: B) =
    herm (herm a) === a

prop_diag_herm1 (BandedAt (a :: B) (k,_)) =
    diagBanded a (-k) === conj (diagBanded (herm a) k)
prop_diag_herm2 (BandedAt (a :: B) (_,k)) =
    diagBanded a k === conj (diagBanded (herm a) (-k))


prop_apply_basis (BandedAt (a :: B) (_,j)) =
    a <*> (basisVector (numCols a) j :: V) ~== col a j
prop_apply_herm_basis (BandedAt (a :: B) (i,_)) =
    (herm a) <*> (basisVector (numRows a) i :: V) ~== conj (row a i)
prop_apply_scale k (BandedMV (a :: B) x) =
    a <*> (k *> x) ~== k *> (a <*> x)
prop_apply_linear (BandedMVPair (a :: B) x y) =
    a <*> (x + y) ~== a <*> x + a <*> y
prop_doSapplyAddVector alpha beta (BandedMV (a :: B) x) = monadicST $ do
    forAllM (Test.vector (numRows a)) $ \y -> do
        y'  <- run $ (unsafeThawVector y :: ST s (STVector s E))
        y'' <- run $ freezeVector y'
        run $ doSApplyAddVector alpha a x beta y'
        assert $ y ~== a <*> (alpha *> x) + (beta *> y'')

prop_applyMatrix_scale_left (BandedMM (a:: B) b) k =
    a <**> (k *> b) ~== k *> (a <**> b)    
prop_applyMatrix_scale_right (BandedMM (a:: B) b) k =
    (k *> a) <**> b ~== k *> (a <**> b)
prop_applyMatrix_linear (BandedMMPair (a :: B) b c) =
    (a <**> (b + c) ~== a <**> b + a <**> c)
prop_applyMatrix_cols (BandedMM (a :: B) b) =
    cols (a <**> b) ~== map (a <*> ) (cols b)
prop_doSapplyAddMatrix alpha beta (BandedMM (a :: B) b) = monadicST $ do
    forAllM (Test.matrix (numRows a, numCols b)) $ \c -> do
        c'  <- run $ unsafeThawMatrix c
        c'' <- run $ freezeMatrix c'
        run $ doSApplyAddMatrix alpha a b beta c'
        assert $ c ~== a <**> (alpha *> b) + (beta *> c'')

prop_scale k (a :: B) =
    k *> a ~== tmap (\e -> e * k) a

tests_Banded =
    [ testProperty "shape of banded" prop_banded_shape
    , testProperty "assocs of banded" prop_banded_assocs
    , testProperty "shape of listsBanded" prop_listsBanded_shape
    , testProperty "listsFromBanded/listsBanded" prop_listsBanded_listsFromBanded

    , testProperty "elems of replace" prop_replace_elems
    
    , testProperty "numRows/numCols" prop_shape
    , testProperty "numLower/numUpper" prop_bandwidths
    , testProperty "size" prop_size
    , testProperty "bounds" prop_bounds

    , testProperty "at" prop_at
    , testProperty "row dim" prop_row_dim
    , testProperty "col dim" prop_col_dim
    , testProperty "rows length" prop_rows_len
    , testProperty "cols length" prop_cols_len
    , testProperty "rows dims" prop_rows_dims
    , testProperty "cols dims" prop_cols_dims
    , testProperty "row at" prop_row_at
    , testProperty "col at" prop_col_at

    , testProperty "indices length" prop_indices_length
    , testProperty "indices low bw" prop_indices_lower
    , testProperty "indices up bw" prop_indices_upper

    , testProperty "elems length" prop_elems_length

    , testProperty "assocs" prop_assocs
    , testProperty "assocs/at" prop_assocs_at

    , testProperty "elems of scale" prop_scale_elems
    , testProperty "elem of herm" prop_herm_elem
    , testProperty "herm/scale" prop_herm_scale
                               
    , testProperty "shape . herm" prop_herm_shape
    , testProperty "rows . herm" prop_herm_rows
    , testProperty "cols . herm" prop_herm_cols
                               
    , testProperty "herm . herm == id" prop_herm_herm
                               
    , testProperty "subdiag . herm" prop_diag_herm1
    , testProperty "superdiag . herm" prop_diag_herm2
                               
    , testProperty "apply basis" prop_apply_basis
    , testProperty "apply herm basis" prop_apply_herm_basis
    , testProperty "apply scale" prop_apply_scale
    , testProperty "apply linear" prop_apply_linear
    , testProperty "doSApplyAddVector" prop_doSapplyAddVector
    
    , testProperty "applyMatrix scale left" prop_applyMatrix_scale_left
    , testProperty "applyMatrix scale right" prop_applyMatrix_scale_right
    , testProperty "applyMatrix linear" prop_applyMatrix_linear
    , testProperty "applyMatrix cols" prop_applyMatrix_cols
    , testProperty "doSApplyAddMatrix" prop_doSapplyAddMatrix    
    
    , testProperty "scale" prop_scale
    
    ]

assocsEq :: [((Int,Int), E)] -> [((Int,Int), E)] -> Bool
assocsEq ies ies' = ordered ies ~== ordered ies'
  where
    ordered = sortAssocs . nubAssocs
    nubAssocs = reverse . nubBy ((==) `on` fst) . reverse      
    sortAssocs = sortBy (comparing fst)
