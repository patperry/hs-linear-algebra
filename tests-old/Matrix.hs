{-# LANGUAGE ScopedTypeVariables #-}

module Matrix
    where

import Driver
import Monadic
import qualified Data.Array as Array

import Data.Elem.BLAS
import Data.Matrix.Dense
import Data.Matrix.Dense.ST
import Data.Vector.Dense
import Data.Vector.Dense.ST

import Test.Matrix.Dense hiding ( matrix )
import qualified Test.QuickCheck.BLAS as Test

type V = Vector E
type M = Matrix E


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
prop_identityMatrix_row (Index m i) (Index _ n) =
    if i < min m n 
        then row (identityMatrix (m,n) :: M) i === basisVector n i
        else row (identityMatrix (m,n) :: M) i === zeroVector n
prop_identityMatrix_col (Pos m) (Index n j) =
    if j < min m n
        then col (identityMatrix (m,n) :: M) j === basisVector m j
        else col (identityMatrix (m,n) :: M) j === zeroVector m

prop_replace_elems (a :: M) =
    forAll (Test.assocs2 (shape a)) $ \ies ->
        let a' = a // ies
            mn = (numRows a - 1, numCols a - 1)
        in sortBy (comparing fst) (assocs a')
           ===
           (Array.assocs $ (Array.//) (Array.array ((0,0),mn) $ assocs a) ies)


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
    sapplyVector k a x ~== k *> (a <*> x)
prop_apply_linear (MatrixMVPair (a :: M) x y) =
    a <*> (x + y) ~== a <*> x + a <*> y
prop_doSapplyAddVector alpha beta (MatrixMV (a ::M) x) = monadicST $ do
    forAllM (Test.vector (numRows a)) $ \y -> do
        y' <- run $ (unsafeThawVector y :: ST s (STVector s E))
        y'' <- run $ freezeVector y'
        run $ doSApplyAddVector alpha a x beta y'
        assert $ y ~== a <*> (alpha *> x) + (beta *> y'')

prop_applyMatrix_id_right (a :: M) =
    let n = numCols a
    in a <**> (identityMatrix (n,n) :: M) ~== a
prop_applyMatrix_id_left (a :: M) =
    let m = numRows a
    in (identityMatrix (m,m) :: M) <**> a ~== a
prop_applyMatrix_scale_left (MatrixMM (a:: M) b) k =
    a <**> (k *> b) ~== k *> (a <**> b)    
prop_applyMatrix_scale_right (MatrixMM (a:: M) b) k =
    (k *> a) <**> b ~== k *> (a <**> b)
prop_applyMatrix_linear (MatrixMMPair (a :: M) b c) =
    a <**> (b + c) ~== a <**> b + a <**> c
prop_applyMatrix_herm (MatrixMM (a :: M) b) =
    herm b <**> herm a ~== herm (a <**> b)
prop_applyMatrix_cols (MatrixMM (a :: M) b) =
    cols (a <**> b) ~== map (a <*> ) (cols b)
prop_doSapplyAddMatrix alpha beta (MatrixMM (a ::M) b) = monadicST $ do
    forAllM (Test.matrix (numRows a, numCols b)) $ \c -> do
        c'  <- run $ unsafeThawMatrix c
        c'' <- run $ freezeMatrix c'
        run $ doSApplyAddMatrix alpha a b beta c'
        assert $ c ~== a <**> (alpha *> b) + (beta *> c'')

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
    [ testProperty "shape of matrix" prop_matrix_shape
    , testProperty "assocs of matrix" prop_matrix_assocs
    , testProperty "shape of listMatrix" prop_listMatrix_shape
    , testProperty "assocs of listMatrix" prop_listMatrix_assocs
    , testProperty "shape of zero" prop_zero_shape
    , testProperty "elems of zero" prop_zero_elems
    
    , testProperty "shape of constant" prop_constant_shape
    , testProperty "elems of constant" prop_constant_elems
    
    , testProperty "shape of identityMatrix" prop_identityMatrix_shape
    , testProperty "diag of identityMatrix" prop_identityMatrix_diag
    , testProperty "row of identityMatrix" prop_identityMatrix_row
    , testProperty "col of identityMatrix" prop_identityMatrix_col
    
    , testProperty "elems of replace" prop_replace_elems
    
    , testProperty "numRows/numCols" prop_shape
    , testProperty "size" prop_size
    , testProperty "bounds" prop_bounds
    , testProperty "at" prop_at
    , testProperty "row dim" prop_row_dim
    , testProperty "col dim" prop_col_dim
    , testProperty "rows length" prop_rows_len
    , testProperty "cols length" prop_cols_len
    , testProperty "rows dims" prop_rows_dims
    , testProperty "cols dims" prop_cols_dims

    , testProperty "indices" prop_indices
    , testProperty "elems" prop_elems
    , testProperty "assocs" prop_assocs
    
    , testProperty "shape of submatrix" prop_submatrix_shape
    , testProperty "rows of submatrix" prop_submatrix_rows
    , testProperty "col of submatrix" prop_submatrix_cols
    
    , testProperty "elems of scale" prop_scale_elems
    , testProperty "elem of herm" prop_herm_elem
    , testProperty "herm/scale" prop_herm_scale
                               
    , testProperty "shape . herm" prop_herm_shape
    , testProperty "rows . herm" prop_herm_rows
    , testProperty "cols . herm" prop_herm_cols
                               
    , testProperty "herm . herm == id" prop_herm_herm
                               
    , testProperty "subdiag . herm" prop_diag_herm1
    , testProperty "superdiag . herm" prop_diag_herm2
                               
    , testProperty "shape . matrixFromRow" prop_matrixFromRow_shape
    , testProperty "elems . matrixFromRow" prop_matrixFromRow_elems
    , testProperty "shape . matrixFromCol" prop_matrixFromCol_shape
    , testProperty "elems . matrixFromCol" prop_matrixFromCol_elems

    , testProperty "apply basis" prop_apply_basis
    , testProperty "apply herm basis" prop_apply_herm_basis
    , testProperty "apply scale" prop_apply_scale
    , testProperty "apply linear" prop_apply_linear
    , testProperty "doSApplyAddVector" prop_doSapplyAddVector
    
    , testProperty "applyMatrix id left" prop_applyMatrix_id_left
    , testProperty "applyMatrix id right" prop_applyMatrix_id_right
    , testProperty "applyMatrix scale left" prop_applyMatrix_scale_left
    , testProperty "applyMatrix scale right" prop_applyMatrix_scale_right
    , testProperty "applyMatrix linear" prop_applyMatrix_linear
    , testProperty "applyMatrix herm" prop_applyMatrix_herm
    , testProperty "applyMatrix cols" prop_applyMatrix_cols
    , testProperty "doSApplyAddMatrix" prop_doSapplyAddMatrix    
    , testProperty "shift" prop_shift
    , testProperty "scale" prop_scale
    
    , testProperty "plus" prop_plus
    , testProperty "minus" prop_minus
    , testProperty "times" prop_times
    , testProperty "divide" prop_divide
    
    , testProperty "negate" prop_negate
    , testProperty "abs" prop_abs
    , testProperty "signum" prop_signum
    , testProperty "recip" prop_recip
    ]


assocsEq :: [((Int,Int), E)] -> [((Int,Int), E)] -> Bool
assocsEq ies ies' = ordered ies === ordered ies'
  where
    ordered = sortAssocs . nubAssocs
    nubAssocs = reverse . nubBy ((==) `on` fst) . reverse      
    sortAssocs = sortBy (comparing fst)
