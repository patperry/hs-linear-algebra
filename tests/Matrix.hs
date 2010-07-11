module Matrix (
    tests_Matrix
    ) where

import Control.Monad( zipWithM_ )
import Data.AEq
import Debug.Trace
import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Numeric.LinearAlgebra.Elem
import Numeric.LinearAlgebra.Vector
import Numeric.LinearAlgebra.Matrix

import Test.QuickCheck.LinearAlgebra( TestElem(..), Dim2(..), Index2(..),
    Assocs2(..), MatrixPair(..) )
import qualified Test.QuickCheck.LinearAlgebra as Test

import Typed


tests_Matrix = testGroup "Matrix"
    [ testPropertyI "dim/matrix" prop_dim_matrix
    , testPropertyI "at/matrix" prop_at_matrix
    , testPropertyI "listMatrix" prop_listMatrix
    , testPropertyI "constantMatrix" prop_constantMatrix
    , testPropertyI "indices" prop_indices
    , testPropertyI "elems" prop_elems
    , testPropertyI "assocs" prop_assocs
    , testPropertyI "replace" prop_replace
    , testPropertyI "accum" prop_accum
    , testPropertyI "map" prop_map
    , testPropertyI "zipWith" prop_zipWith
    , testPropertyI "col" prop_col
    , testPropertyI "cols" prop_cols
    , testPropertyI "splice" prop_splice
    , testPropertyI "splitRowsAt" prop_splitRowsAt
    , testPropertyI "splitColsAt" prop_splitColsAt
    , testPropertyDZ "shift" prop_shift prop_shift
    , testPropertyDZ "shiftDiag" prop_shiftDiag prop_shiftDiag
    , testPropertyDZ "shiftDiagWithScale" prop_shiftDiagWithScale prop_shiftDiagWithScale
    , testPropertyDZ "add" prop_add prop_add
    , testPropertyDZ "addWithScale" prop_addWithScale prop_addWithScale
    , testPropertyDZ "sub" prop_sub prop_sub
    , testPropertyDZ "scale" prop_scale prop_scale
    , testPropertyDZ "scaleRows" prop_scaleRows prop_scaleRows
    , testPropertyDZ "scaleCols" prop_scaleCols prop_scaleCols
    , testPropertyDZ "negate" prop_negate prop_negate
    ]



------------------------- Matrix Construction ------------------------------

prop_dim_matrix t (Assocs2 mn ies) =
    dimMatrix (matrix mn ies) === mn
  where
    _ = typed t $ matrix mn ies

prop_at_matrix t (Assocs2 mn ies) = let
    x = matrix mn ies
    is = (fst . unzip) ies
    in and [ atMatrix x i `elem` [ e | (i',e) <- ies, i' == i ]
           | i <- is]
  where
    _ = typed t $ matrix mn ies

prop_listMatrix t (Dim2 (m,n)) =
    forAll (QC.vector $ m*n) $ \es ->
        listMatrix (m,n) es === (typed t $ matrix (m,n) $ 
            zip [ (i,j) | j <- [ 0..n-1], i <- [ 0..m-1 ] ] es)

prop_constantMatrix t (Dim2 (m,n)) e =
    constantMatrix (m,n) e === listMatrix (m,n) (replicate (m*n) e)
  where
    _ = typed t [e]


-------------------------- Accessing Matrices ------------------------------

prop_indices t x =
    indicesMatrix x === [ (i,j) | j <- [ 0..n-1 ], i <- [ 0..m-1 ] ]
  where
    (m,n) = dimMatrix x
    _ = immutableMatrix x
    _ = typed t x

prop_elems t x =
    elemsMatrix x === [ atMatrix x i | i <- indicesMatrix x ]
  where
    _ = typed t x
    
prop_assocs t x =
    assocsMatrix x === zip (indicesMatrix x) (elemsMatrix x)
  where
    _ = typed t x


------------------------- Incremental Updates ------------------------------
    
prop_replace t (Assocs2 mn ies) =
    forAll (typed t `fmap` Test.matrix mn) $ \x -> let
        x' = replaceMatrix x ies
        is = indicesMatrix x
        is1 = (fst . unzip) ies
        is0 = [ i | i <- is, i `notElem` is1 ]
        in and $
            [ atMatrix x' i `elem` [ e | (i',e) <- ies, i' == i ]
            | i <- is1
            ] ++
            [ atMatrix x' i === atMatrix x i
            | i <- is0
            ]

prop_accum t (Blind f) (Assocs2 mn ies) =
    forAll (typed t `fmap` Test.matrix mn) $ \x -> let
        x' = accumMatrix f x ies
        in x' === listMatrix mn [ foldl f e [ e' | (i',e') <- ies, i' == i]
                                | (i,e) <- assocsMatrix x ]
  where
      _ = typed t $ (snd . unzip) ies

-------------------------- Derived Matrices ------------------------------
     
prop_map t (Blind f) x =
    mapMatrix f x === listMatrix (dimMatrix x) (map f $ elemsMatrix x)
  where
    _ = typed t x
    _ = typed t $ mapMatrix f x

prop_zipWith t (Blind f) (MatrixPair x y) =
    zipWithMatrix f x y === (listMatrix (dimMatrix x) $
                                zipWith f (elemsMatrix x) (elemsMatrix y))
  where
    _ = typed t x
    _ = typed t y    
    _ = typed t $ zipWithMatrix f x y
     

------------------------------ Matrix Views --------------------------------

prop_col t (Index2 (m,n) (_,j)) =
    forAll (typed t `fmap` Test.matrix (m,n)) $ \a ->
        colMatrix a j === listVector m [ atMatrix a (i,j) | i <- [ 0..m-1 ] ]

prop_cols t a =
    colsMatrix a === [ colMatrix a j | j <- [ 0..n-1 ] ]
  where
    (_,n) = dimMatrix a
    _ = typed t $ immutableMatrix a

prop_splice t a =
    forAll (choose (0,m)) $ \m' ->
    forAll (choose (0,n)) $ \n' ->
    forAll (choose (0,m-m')) $ \i ->
    forAll (choose (0,n-n')) $ \j ->
        spliceMatrix a (i,j) (m',n')
            === colListMatrix (m',n') [ spliceVector (colMatrix a j') i m'
                                      | j' <- [ j..j+n'-1 ] ]
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_splitRowsAt t a =
    forAll (choose (0,m)) $ \i ->
        splitRowsMatrixAt i a
            === ( spliceMatrix a (0,0) (i,n)
                , spliceMatrix a (i,0) (m-i,n)
                )
  where
    (m,n) = dimMatrix a
    _  = typed t $ immutableMatrix a

prop_splitColsAt t a =
    forAll (choose (0,n)) $ \j ->
        splitColsMatrixAt j a
            === ( spliceMatrix a (0,0) (m,j)
                , spliceMatrix a (0,j) (m,n-j)
                )
  where
    (m,n) = dimMatrix a
    _  = typed t $ immutableMatrix a
    

-------------------------- Num Matrix Operations --------------------------

prop_shift t k a =
    k `shiftMatrix` a === mapMatrix (k+) a
  where
    _ = typed t a

prop_shiftDiag t a =
    forAll (Test.vector (min m n)) $ \d ->
        d `shiftDiagMatrix` a
            === accumMatrix (+) a [ ((i,i),e) | (i,e) <- assocsVector d ]
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_shiftDiagWithScale t k a =
    forAll (Test.vector (min m n)) $ \d ->
        shiftDiagMatrixWithScale k d a
            ~== accumMatrix (+) a [ ((i,i),k * e) | (i,e) <- assocsVector d ]
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_add t (MatrixPair x y) =
    x `addMatrix` y === zipWithMatrix (+) x y
  where
    _ = typed t x

prop_addWithScale t a b (MatrixPair x y) =
    addMatrixWithScale a x b y ~== 
        zipWithMatrix (+) (mapMatrix (a*) x) (mapMatrix (b*) y)
  where
    _ = typed t x

prop_sub t (MatrixPair x y) =
    x `subMatrix` y === zipWithMatrix (-) x y
  where
    _ = typed t x

prop_scale t k x =
    k `scaleMatrix` x === mapMatrix (k*) x
  where
    _ = typed t x

prop_scaleRows t a =
    forAll (Test.vector m) $ \s ->
        scaleRowsMatrix s a
            === colListMatrix (m,n) [ mulVector s x | x <- colsMatrix a ]
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_scaleCols t a =
    forAll (Test.vector n) $ \s ->
        scaleColsMatrix s a
            === colListMatrix (m,n)
                    [ scaleVector e x 
                    | (e,x) <- zip (elemsVector s) (colsMatrix a) ]
  where
    (m,n) = dimMatrix a
    _ = typed t a
    
prop_negate t x =
    negateMatrix x === mapMatrix negate x
  where
    _ = typed t x
