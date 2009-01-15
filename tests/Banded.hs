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

import Data.Elem.BLAS
import Data.Vector.Dense 
import Data.Matrix.Dense ( Matrix, identityMatrix )
import Data.Matrix.Banded

import Test.Matrix.Banded hiding ( banded )
        

listsFromBanded :: (BLAS1 e) => Banded np e -> ((Int,Int), (Int,Int),[[e]])
listsFromBanded a = ( (m,n)
            , (kl,ku)
            , map paddedDiag [(-kl)..ku]
            )
  where
    (m,n)   = shape a
    (kl,ku) = bandwidths (coerceBanded a)
    
    padBegin i   = replicate (max (-i) 0)    0
    padEnd   i   = replicate (max (m-n+i) 0) 0
    paddedDiag i = (  padBegin i
                   ++ elems (diagBanded (coerceBanded a) i) 
                   ++ padEnd i 
                   )
                   
type V = Vector Int E
type M = Matrix (Int,Int) E
type B = Banded (Int,Int) E

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

prop_applyMatrix_scale_left (BandedMM (a:: B) b) k =
    a <**> (k *> b) ~== k *> (a <**> b)    
prop_applyMatrix_scale_right (BandedMM (a:: B) b) k =
    (k *> a) <**> b ~== k *> (a <**> b)
prop_applyMatrix_linear (BandedMMPair (a :: B) b c) =
    (a <**> (b + c) ~== a <**> b + a <**> c)
prop_applyMatrix_cols (BandedMM (a :: B) b) =
    cols (a <**> b) ~== map (a <*> ) (cols b)

prop_scale k (a :: B) =
    k *> a ~== tmap (\e -> e * k) a

tests_Banded =
    [ ("shape of banded"       , mytest prop_banded_shape)
    , ("assocs of banded"      , mytest prop_banded_assocs)
    , ("shape of listsBanded"  , mytest prop_listsBanded_shape)
    , ("listsFromBanded/listsBanded"
                               , mytest prop_listsBanded_listsFromBanded)

    , ("elems of replace"      , mytest prop_replace_elems)
    
    , ("numRows/numCols"       , mytest prop_shape)
    , ("numLower/numUpper"     , mytest prop_bandwidths)
    , ("size"                  , mytest prop_size)
    , ("bounds"                , mytest prop_bounds)

    , ("at"                    , mytest prop_at)
    , ("row dim"               , mytest prop_row_dim)
    , ("col dim"               , mytest prop_col_dim)
    , ("rows length"           , mytest prop_rows_len)
    , ("cols length"           , mytest prop_cols_len)
    , ("rows dims"             , mytest prop_rows_dims)
    , ("cols dims"             , mytest prop_cols_dims)
    , ("row at"                , mytest prop_row_at)
    , ("col at"                , mytest prop_col_at)

    , ("indices length"        , mytest prop_indices_length)
    , ("indices low bw"        , mytest prop_indices_lower)
    , ("indices up bw"         , mytest prop_indices_upper)

    , ("elems length"          , mytest prop_elems_length)

    , ("assocs"                , mytest prop_assocs)
    , ("assocs/at"             , mytest prop_assocs_at)

    , ("elems of scale"        , mytest prop_scale_elems)
    , ("elem of herm"          , mytest prop_herm_elem)
    , ("herm/scale"            , mytest prop_herm_scale)
                               
    , ("shape . herm"          , mytest prop_herm_shape)
    , ("rows . herm"           , mytest prop_herm_rows)
    , ("cols . herm"           , mytest prop_herm_cols)
                               
    , ("herm . herm == id"     , mytest prop_herm_herm)
                               
    , ("subdiag . herm"        , mytest prop_diag_herm1)
    , ("superdiag . herm"      , mytest prop_diag_herm2)
                               
    , ("apply basis"           , mytest prop_apply_basis)
    , ("apply herm basis"      , mytest prop_apply_herm_basis)
    , ("apply scale"           , mytest prop_apply_scale)
    , ("apply linear"          , mytest prop_apply_linear)
    
    , ("applyMatrix scale left"   , mytest prop_applyMatrix_scale_left)
    , ("applyMatrix scale right"  , mytest prop_applyMatrix_scale_right)
    , ("applyMatrix linear"       , mytest prop_applyMatrix_linear)
    , ("applyMatrix cols"         , mytest prop_applyMatrix_cols)
    
    , ("scale"                 , mytest prop_scale)
    
    ]

assocsEq :: [((Int,Int), E)] -> [((Int,Int), E)] -> Bool
assocsEq ies ies' = ordered ies ~== ordered ies'
  where
    ordered = sortAssocs . nubAssocs
    nubAssocs = reverse . nubBy ((==) `on` fst) . reverse      
    sortAssocs = sortBy (comparing fst)
