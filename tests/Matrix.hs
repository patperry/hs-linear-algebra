{-# OPTIONS -fglasgow-exts -fno-excess-precision -cpp #-}
-----------------------------------------------------------------------------
-- |
-- Copyright  : Copyright (c) 2008, Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

import Debug.Trace ( trace )

import qualified Data.Array as Array
import Data.Ix   ( inRange, range )
import Data.List ( nub, sortBy )
import Data.Ord  ( comparing )
import System.Environment ( getArgs )
import Test.QuickCheck.Parallel hiding ( vector )
import qualified Test.QuickCheck as QC

import BLAS.Access
import Data.Matrix.Dense
import Data.Vector.Dense.Internal ( DVector )
import Data.Matrix.Dense.Internal ( DMatrix )
import Data.Vector.Dense hiding ( shift, scale, invScale )
import qualified Data.Vector.Dense as V
import BLAS.Elem ( Elem, BLAS1 )
import qualified BLAS.Elem as E
import Data.Complex ( Complex(..) )

import Data.AEq
import Numeric.IEEE

import Test.QuickCheck.Complex
import Test.QuickCheck.Vector hiding ( Assocs )
import Test.QuickCheck.Vector.Dense hiding ( Pair )
import Test.QuickCheck.Matrix
import Test.QuickCheck.Matrix.Dense

import Debug.Trace
        

#ifdef COMPLEX
field = "Complex Double"
type E = Complex Double
#else
field = "Double"
type E = Double
#endif

type V = Vector Int E
type M = Matrix (Int,Int) E

instance (Arbitrary e, RealFloat e) => Arbitrary (Complex e) where
    arbitrary   = arbitrary >>= \(TestComplex x) -> return x
    coarbitrary = coarbitrary . TestComplex

instance (Arbitrary e, BLAS1 e) => Arbitrary (DVector Imm n e) where
    arbitrary   = arbitrary >>= \(TestVector x) -> return x
    coarbitrary = coarbitrary . TestVector

instance (Arbitrary e, BLAS1 e) => Arbitrary (DMatrix Imm (m,n) e) where
    arbitrary   = arbitrary >>= \(TestMatrix x) -> return x
    coarbitrary = coarbitrary . TestMatrix


assocsEq :: (BLAS1 e, AEq e) => Matrix (m,n) e -> [((Int,Int), e)] -> Bool
assocsEq x ijes =
    let ijs = fst $ unzip ijes
    in filter (\(ij,e) -> ij `elem` ijs) (sortBy (comparing fst) $ assocs x) === sortBy (comparing fst) ijes
       && (all (==0) $ map snd $ filter (\(ij,e) -> not $ ij `elem` ijs) $ assocs x)

prop_matrix_shape (Assocs mn ijes) =
    shape (matrix mn ijes :: M) == mn
prop_matrix_assocs (Assocs mn ijes) =
    (matrix mn ijes :: M) `assocsEq` ijes

prop_listMatrix_shape (IndexPair mn) es =
    shape (listMatrix mn es :: M) == mn
prop_listMatrix_assocs (IndexPair (m,n)) es =
    let es' = repeat es
    in assocs (listMatrix (m,n) es' :: M) === zip [(i,j) | j <- range (0,n-1), i <- range (0,m-1)] es'

prop_zero_shape (IndexPair mn) =
    shape (zero mn :: M) == mn
prop_zero_elems (IndexPair (m,n)) =
    elems (zero (m,n) :: M) == replicate (m*n) 0

prop_constant_shape (IndexPair mn) (e :: E) =
    shape (constant mn e :: M) == mn
prop_constant_elems (IndexPair (m,n)) (e :: E) =
    elems (constant (m,n) e :: M) == replicate (m*n) e

prop_identity_shape (IndexPair mn) =
    shape (identity mn :: M) == mn
prop_identity_diag (IndexPair (m,n)) =
    diag (identity (m,n) :: M) 0 === (constant (min m n) 1)
prop_identity_row (Basis m i) (Basis n _) =
    if i < min m n 
        then row (identity (m,n) :: M) i === V.basis n i
        else row (identity (m,n) :: M) i === V.zero n
prop_identity_col (Basis m _) (Basis n j) =
    if j < min m n
        then col (identity (m,n) :: M) j === V.basis m j
        else col (identity (m,n) :: M) j === V.zero m

prop_replace_elems (a :: M) (Assocs _ ijes) =
    let ijes' = filter (\((i,j),_) -> i < numRows a && j < numCols a) ijes
        a'   = a // ijes'
        mn   = (numRows a - 1, numCols a - 1)
    in and $ zipWith (\(ij1,e1) (ij2,e2) -> (ij1 == ij2) && (e1 === e2))
                     (sortBy (comparing fst) $ assocs a')
                     (Array.assocs $ (Array.//) (Array.array ((0,0),mn) $ assocs a) ijes')


prop_submatrix_shape (SubMatrix a ij mn) =
    shape (submatrix a ij mn :: M) == mn
prop_submatrix_rows (SubMatrix a (i,j) (m,n)) =
    rows (submatrix a (i,j) (m,n) :: M) === map (\k -> V.subvector (row a (i+k)) j n) [0..(m-1)]
prop_submatrix_cols (SubMatrix a (i,j) (m,n)) (Index l) =
    cols (submatrix a (i,j) (m,n) :: M) === map (\l -> V.subvector (col a (j+l)) i m) [0..(n-1)]

prop_shape (a :: M) = 
    shape a == (numRows a, numCols a)
prop_size (a :: M) =
    size a == numRows a * numCols a
prop_bounds (a :: M) =
    bounds a == ((0,0), (numRows a - 1, numCols a - 1))
    
prop_at (MatAt (a :: M) (i,j)) =
    let ij = (i,j)
        k  = i + j * numRows a
    in (a!ij) === ((elems a) !! k)
    
prop_row_dim (MatAt (a :: M) (i,_)) =
    V.dim (row a i) == numCols a
prop_col_dim (MatAt (a :: M) (_,j)) =
    V.dim (col a j) == numRows a
prop_rows_len (a :: M) =
    length (rows a) == numRows a
prop_cols_len (a :: M) =
    length (cols a) == numCols a
prop_rows_dims (a :: M) =
    map (V.dim) (rows a) == replicate (numRows a) (numCols a)
prop_cols_dims (a :: M) =
    map (V.dim) (cols a) == replicate (numCols a) (numRows a)

prop_indices (a :: M) =
    let (m,n) = shape a
    in indices a == [(i,j) | j <- range (0,n-1), i <- range(0,m-1)]
prop_elems (a :: M) =
    and $ zipWith (===) (elems a) $ concatMap V.elems (cols a)
prop_assocs (a :: M) = 
    assocs a === zip (indices a) (elems a)

prop_scale_elems (a :: M) k =
    and $ zipWith (~==) (elems (k *> a)) (map (k*) (elems a))
prop_herm_elem (MatAt (a :: M) (i,j)) =
    (herm a) ! (j,i) == E.conj (a!(i,j))
prop_herm_scale (a :: M) k =
    herm (scale k a) === scale (E.conj k) (herm a)

prop_herm_shape (a :: M) =
    shape (herm a) == (numCols a, numRows a)
prop_herm_rows (a :: M) =
    rows (herm a) === map (V.conj) (cols a)
prop_herm_cols (a :: M) = 
    cols (herm a) === map (V.conj) (rows a)

prop_herm_herm (a :: M) =
    herm (herm a) === a

prop_diag_herm1 (MatAt (a :: M) (k,_)) =
    diag a (-k) === V.conj (diag (herm a) k)
prop_diag_herm2 (MatAt (a :: M) (_,k)) =
    diag a k === V.conj (diag (herm a) (-k))

prop_fromRow_shape (x :: V) =
    shape (fromRow x :: M) == (1,V.dim x)
prop_fromRow_elems (x :: V) =
    elems (fromRow x :: M) === V.elems x

prop_fromCol_shape (x :: V) =
    shape (fromCol x :: M) == (V.dim x,1)
prop_fromCol_elems (x :: V) =
    elems (fromCol x :: M) === V.elems x


prop_apply_basis (MatAt (a :: M) (_,j)) =
    a <*> (V.basis (numCols a) j :: V) ~== col a j
prop_apply_herm_basis (MatAt (a :: M) (i,_)) =
    (herm a) <*> (V.basis (numRows a) i :: V) ~== V.conj (row a i)
prop_apply_scale k (MultMV (a :: M) x) =
    a <*> (V.scale k x) ~== V.scale k (a <*> x)
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

prop_shift k (a :: M) =
    shift k a ~== a + constant (shape a) k
prop_scale k (a :: M) =
    scale k a ~== a * constant (shape a) k
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
    negate a ~== scale (-1) a
    
prop_abs (a :: M) =
    elems (abs a) ~== map abs (elems a)
prop_signum (a :: M) =
    elems (signum a) === map signum (elems a)
prop_recip (a :: M) =
    elems (recip a) ~== (map recip $ elems a)

properties =
    [ ("shape of matrix"       , pDet prop_matrix_shape)
    , ("assocs of matrix"      , pDet prop_matrix_assocs)
    , ("shape of listMatrix"   , pDet prop_listMatrix_shape)
    , ("assocs of listMatrix"  , pDet prop_listMatrix_assocs)
    , ("shape of zero"         , pDet prop_zero_shape)
    , ("elems of zero"         , pDet prop_zero_elems)
    
    , ("shape of constant"     , pDet prop_constant_shape)
    , ("elems of constant"     , pDet prop_constant_elems)
    
    , ("shape of identity"     , pDet prop_identity_shape)
    , ("diag of identity"      , pDet prop_identity_diag)
    , ("row of identity"       , pDet prop_identity_row)
    , ("col of identity"       , pDet prop_identity_col)
    
    , ("elems of replace"      , pDet prop_replace_elems)
    
    , ("numRows/numCols"       , pDet prop_shape)
    , ("size"                  , pDet prop_size)
    , ("bounds"                , pDet prop_bounds)
    , ("at"                    , pDet prop_at)
    , ("row dim"               , pDet prop_row_dim)
    , ("col dim"               , pDet prop_col_dim)
    , ("rows length"           , pDet prop_rows_len)
    , ("cols length"           , pDet prop_cols_len)
    , ("rows dims"             , pDet prop_rows_dims)
    , ("cols dims"             , pDet prop_cols_dims)

    , ("indices"               , pDet prop_indices)
    , ("elems"                 , pDet prop_elems)
    , ("assocs"                , pDet prop_assocs)
    
    , ("shape of submatrix"    , pDet prop_submatrix_shape)
    , ("rows of submatrix"     , pDet prop_submatrix_rows)
    , ("col of submatrix"      , pDet prop_submatrix_cols)
    
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
    
    , ("shift"                 , pDet prop_shift)
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
