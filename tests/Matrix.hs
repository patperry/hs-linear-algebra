{-# LANGUAGE ScopedTypeVariables #-}

module Matrix( tests_Matrix ) where

import Driver
import qualified Data.Array as Array

import Data.Elem.BLAS
import Data.Matrix.Dense
import Data.Vector.Dense

import Test.Matrix.Dense hiding ( matrix )
        

type V = Vector Int E
type M = Matrix (Int,Int) E


prop_matrix_shape (Assocs2 mn ijes) =
    shape (matrix mn ijes :: M) == mn
prop_matrix_assocs (Assocs2 (m,n) ijes) =
    assocs (matrix (m,n) ijes :: M) 
    `assocsEq` 
    (zip (range ((0,0),(m-1,n-1))) (repeat 0) ++ ijes)

prop_listMatrix_shape (Nat2 mn) es =
    shape (listMatrix mn es :: M) == mn
prop_listMatrix_assocs (Nat2 (m,n)) es =
    let es' = repeat es
    in assocs (listMatrix (m,n) es' :: M) 
       === 
       zip [(i,j) | j <- range (0,n-1), i <- range (0,m-1)] es'

prop_zero_shape (Nat2 mn) =
    shape (zeroMatrix mn :: M) == mn
prop_zero_elems (Nat2 (m,n)) =
    elems (zeroMatrix (m,n) :: M) == replicate (m*n) 0

prop_constant_shape (Nat2 mn) (e :: E) =
    shape (constantMatrix mn e :: M) == mn
prop_constant_elems (Nat2 (m,n)) (e :: E) =
    elems (constantMatrix (m,n) e :: M) == replicate (m*n) e

prop_identityMatrix_shape (Nat2 mn) =
    shape (identityMatrix mn :: M) == mn
prop_identityMatrix_diag (Nat2 (m,n)) =
    diag (identityMatrix (m,n) :: M) 0 === (constantVector (min m n) 1)
prop_identityMatrix_row (Index i m) (Index _ n) =
    if i < min m n 
        then row (identityMatrix (m,n) :: M) i === basisVector n i
        else row (identityMatrix (m,n) :: M) i === zeroVector n
prop_identityMatrix_col (Index _ m) (Index j n) =
    if j < min m n
        then col (identityMatrix (m,n) :: M) j === basisVector m j
        else col (identityMatrix (m,n) :: M) j === zeroVector m

prop_replace_elems (a :: M) (Assocs2 _ ijes) =
    let ijes' = filter (\((i,j),_) -> i < numRows a && j < numCols a) ijes
        a'   = a // ijes'
        mn   = (numRows a - 1, numCols a - 1)
    in and $ zipWith (\(ij1,e1) (ij2,e2) -> (ij1 == ij2) && (e1 === e2))
                     (sortBy (comparing fst) $ assocs a')
                     (Array.assocs $ (Array.//) (Array.array ((0,0),mn) $ assocs a) ijes')


prop_submatrix_shape (SubMatrix a ij mn) =
    shape (submatrix a ij mn :: M) == mn
prop_submatrix_rows (SubMatrix a (i,j) (m,n)) =
    rows (submatrix a (i,j) (m,n) :: M) === map (\k -> subvector (row a (i+k)) j n) [0..(m-1)]
prop_submatrix_cols (SubMatrix a (i,j) (m,n)) =
    cols (submatrix a (i,j) (m,n) :: M) === map (\l -> subvector (col a (j+l)) i m) [0..(n-1)]

prop_shape (a :: M) = 
    shape a == (numRows a, numCols a)
prop_size (a :: M) =
    size a == numRows a * numCols a
prop_bounds (a :: M) =
    bounds a == ((0,0), (numRows a - 1, numCols a - 1))
    
prop_at (MatrixAt (a :: M) (i,j)) =
    let ij = (i,j)
        k  = if isHermMatrix a then j + i * numCols a
                               else i + j * numRows a
    in (a!ij) === ((elems a) !! k)
    
prop_row_dim (MatrixAt (a :: M) (i,_)) =
    dim (row a i) == numCols a
prop_col_dim (MatrixAt (a :: M) (_,j)) =
    dim (col a j) == numRows a
    
prop_rows_len (a :: M) =
    length (rows a) == numRows a
prop_cols_len (a :: M) =
    length (cols a) == numCols a
prop_rows_dims (a :: M) =
    map dim (rows a) == replicate (numRows a) (numCols a)
prop_cols_dims (a :: M) =
    map dim (cols a) == replicate (numCols a) (numRows a)

prop_indices (a :: M)
    | isHermMatrix a =
        indices a == [(i,j) | i <- range (0,m-1), j <- range(0,n-1)]
    | otherwise =
        indices a == [(i,j) | j <- range (0,n-1), i <- range(0,m-1)]
  where (m,n) = shape a
  
prop_elems (a :: M) 
    | isHermMatrix a =
        elems a === concatMap elems (rows a)        
    | otherwise =
        elems a === concatMap elems (cols a)
        
prop_assocs (a :: M) = 
    assocs a === zip (indices a) (elems a)

prop_scale_elems (a :: M) k =
    and $ zipWith (~==) (elems (k *> a)) (map (k*) (elems a))
prop_herm_elem (MatrixAt (a :: M) (i,j)) =
    (herm a) ! (j,i) == conjugate (a!(i,j))
prop_herm_scale (a :: M) k =
    herm (k *> a) === (conjugate k) *> (herm a)

prop_herm_shape (a :: M) =
    shape (herm a) == (numCols a, numRows a)
prop_herm_rows (a :: M) =
    rows (herm a) === map conj (cols a)
prop_herm_cols (a :: M) = 
    cols (herm a) === map conj (rows a)

prop_herm_herm (a :: M) =
    herm (herm a) === a

prop_diag_herm1 (MatrixAt (a :: M) (k,_)) =
    diag a (-k) === conj (diag (herm a) k)
prop_diag_herm2 (MatrixAt (a :: M) (_,k)) =
    diag a k === conj (diag (herm a) (-k))

prop_matrixFromRow_shape (x :: V) =
    shape (matrixFromRow x :: M) == (1,dim x)
prop_matrixFromRow_elems (x :: V) =
    elems (matrixFromRow x :: M) === elems x

prop_matrixFromCol_shape (x :: V) =
    shape (matrixFromCol x :: M) == (dim x,1)
prop_matrixFromCol_elems (x :: V) =
    elems (matrixFromCol x :: M) === elems x

prop_apply_basis (MatrixAt (a :: M) (_,j)) =
    a <*> (basisVector (numCols a) j :: V) ~== col a j
prop_apply_herm_basis (MatrixAt (a :: M) (i,_)) =
    (herm a) <*> (basisVector (numRows a) i :: V) ~== conj (row a i)
prop_apply_scale k (MatrixMV (a :: M) x) =
    sapply k a x ~== k *> (a <*> x)
prop_apply_linear (MatrixMVPair (a :: M) x y) =
    a <*> (x + y) ~== a <*> x + a <*> y

prop_applyMat_id_right (a :: M) =
    let n = numCols a
    in a <**> (identityMatrix (n,n) :: M) ~== a
prop_applyMat_id_left (a :: M) =
    let m = numRows a
    in (identityMatrix (m,m) :: M) <**> a ~== a
prop_applyMat_scale_left (MatrixMM (a:: M) b) k =
    a <**> (k *> b) ~== k *> (a <**> b)    
prop_applyMat_scale_right (MatrixMM (a:: M) b) k =
    (k *> a) <**> b ~== k *> (a <**> b)
prop_applyMat_linear (MatrixMMPair (a :: M) b c) =
    a <**> (b + c) ~== a <**> b + a <**> c
prop_applyMat_herm (MatrixMM (a :: M) b) =
    herm b <**> herm a ~== herm (a <**> b)
prop_applyMat_cols (MatrixMM (a :: M) b) =
    cols (a <**> b) ~== map (a <*> ) (cols b)

prop_shift k (a :: M) =
    shift k a ~== a + constantMatrix (shape a) k
prop_scale k (a :: M) =
    k *> a ~== a * constantMatrix (shape a) k

prop_plus (MatrixPair (a :: M) b) =
    colElems (a + b) ~== zipWith (+) (colElems a) (colElems b)
prop_minus (MatrixPair (a :: M) b) =
    colElems (a - b) ~== zipWith (-) (colElems a) (colElems b)
prop_times (MatrixPair (a :: M) b) =
    colElems (a * b) ~== zipWith (*) (colElems a) (colElems b)
prop_divide (MatrixPair (a :: M) b) =
    colElems (a / b) ~== zipWith (/) (colElems a) (colElems b)

colElems a = concatMap elems (cols a)

prop_negate (a :: M) =
    negate a ~== (-1) *> a
    
prop_abs (a :: M) =
    elems (abs a) ~== map abs (elems a)
prop_signum (a :: M) =
    elems (signum a) === map signum (elems a)
prop_recip (a :: M) =
    elems (recip a) ~== (map recip $ elems a)

tests_Matrix =
    [ ("shape of matrix"       , mytest prop_matrix_shape)
    , ("assocs of matrix"      , mytest prop_matrix_assocs)
    , ("shape of listMatrix"   , mytest prop_listMatrix_shape)
    , ("assocs of listMatrix"  , mytest prop_listMatrix_assocs)
    , ("shape of zero"         , mytest prop_zero_shape)
    , ("elems of zero"         , mytest prop_zero_elems)
    
    , ("shape of constant"     , mytest prop_constant_shape)
    , ("elems of constant"     , mytest prop_constant_elems)
    
    , ("shape of identityMatrix"     , mytest prop_identityMatrix_shape)
    , ("diag of identityMatrix"      , mytest prop_identityMatrix_diag)
    , ("row of identityMatrix"       , mytest prop_identityMatrix_row)
    , ("col of identityMatrix"       , mytest prop_identityMatrix_col)
    
    , ("elems of replace"      , mytest prop_replace_elems)
    
    , ("numRows/numCols"       , mytest prop_shape)
    , ("size"                  , mytest prop_size)
    , ("bounds"                , mytest prop_bounds)
    , ("at"                    , mytest prop_at)
    , ("row dim"               , mytest prop_row_dim)
    , ("col dim"               , mytest prop_col_dim)
    , ("rows length"           , mytest prop_rows_len)
    , ("cols length"           , mytest prop_cols_len)
    , ("rows dims"             , mytest prop_rows_dims)
    , ("cols dims"             , mytest prop_cols_dims)

    , ("indices"               , mytest prop_indices)
    , ("elems"                 , mytest prop_elems)
    , ("assocs"                , mytest prop_assocs)
    
    , ("shape of submatrix"    , mytest prop_submatrix_shape)
    , ("rows of submatrix"     , mytest prop_submatrix_rows)
    , ("col of submatrix"      , mytest prop_submatrix_cols)
    
    , ("elems of scale"        , mytest prop_scale_elems)
    , ("elem of herm"          , mytest prop_herm_elem)
    , ("herm/scale"            , mytest prop_herm_scale)
                               
    , ("shape . herm"          , mytest prop_herm_shape)
    , ("rows . herm"           , mytest prop_herm_rows)
    , ("cols . herm"           , mytest prop_herm_cols)
                               
    , ("herm . herm == id"     , mytest prop_herm_herm)
                               
    , ("subdiag . herm"        , mytest prop_diag_herm1)
    , ("superdiag . herm"      , mytest prop_diag_herm2)
                               
    , ("shape . matrixFromRow"       , mytest prop_matrixFromRow_shape)
    , ("elems . matrixFromRow"       , mytest prop_matrixFromRow_elems)
    , ("shape . matrixFromCol"       , mytest prop_matrixFromCol_shape)
    , ("elems . matrixFromCol"       , mytest prop_matrixFromCol_elems)

    , ("apply basis"           , mytest prop_apply_basis)
    , ("apply herm basis"      , mytest prop_apply_herm_basis)
    , ("apply scale"           , mytest prop_apply_scale)
    , ("apply linear"          , mytest prop_apply_linear)
    
    , ("applyMat id left"       , mytest prop_applyMat_id_left)
    , ("applyMat id right"      , mytest prop_applyMat_id_right)
    , ("applyMat scale left"    , mytest prop_applyMat_scale_left)
    , ("applyMat scale right"   , mytest prop_applyMat_scale_right)
    , ("applyMat linear"        , mytest prop_applyMat_linear)
    , ("applyMat herm"          , mytest prop_applyMat_herm)
    , ("applyMat cols"          , mytest prop_applyMat_cols)
    , ("shift"                 , mytest prop_shift)
    , ("scale"                 , mytest prop_scale)
    
    , ("plus"                  , mytest prop_plus)
    , ("minus"                 , mytest prop_minus)
    , ("times"                 , mytest prop_times)
    , ("divide"                , mytest prop_divide)
    
    , ("negate"                , mytest prop_negate)
    , ("abs"                   , mytest prop_abs)
    , ("signum"                , mytest prop_signum)
    , ("recip"                 , mytest prop_recip)
    ]


assocsEq :: [((Int,Int), E)] -> [((Int,Int), E)] -> Bool
assocsEq ies ies' = ordered ies === ordered ies'
  where
    ordered = sortAssocs . nubAssocs
    nubAssocs = reverse . nubBy ((==) `on` fst) . reverse      
    sortAssocs = sortBy (comparing fst)
