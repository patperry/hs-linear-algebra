module Matrix (
    tests_Matrix
    ) where

import Data.AEq
import Debug.Trace
import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import BLAS.Elem
import BLAS.Vector
import BLAS.Matrix

import Test.QuickCheck.BLAS( TestElem(..), Dim2(..), Assocs2(..),
    MatrixPair(..) )
import qualified Test.QuickCheck.BLAS as Test

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
    {-, testPropertyI "splice" prop_splice
    , testPropertyI "splitRowsAt" prop_splitRowsAt
    , testPropertyI "splitColsAt" prop_splitColsAt
    , testPropertyDZ "shift" prop_shift prop_shift
    , testPropertyDZ "shiftDiag" prop_shiftDiag prop_shiftDiag
    , testPropertyDZ "shiftDiagWithScale" prop_shiftDiagWithScale prop_shiftDiagWithScale
    , testPropertyDZ "add" prop_add prop_add
    , testPropertyDZ "addWithScale" prop_addWithScale prop_addWithScale
    , testPropertyDZ "sub" prop_sub prop_sub
    , testPropertyDZ "scale" prop_scale prop_scale
    , testPropertyDZ "negate" prop_negate prop_negate -}
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
     

------------------------------ Matrix Views-- --------------------------------

{-
prop_splice t x = 
    forAll (choose (0,n)) $ \n' ->
    forAll (choose (0,n-n')) $ \o ->
        spliceVector x o n' === listVector n' (take n' $ drop o $ es)
  where
    n  = dimVector x
    es = elemsVector x
    _  = typed t x

prop_splitAt t x =
    forAll (choose (0,n)) $ \k ->
        splitVectorAt k x === (listVector k $ take k es,
                               listVector (n-k) $ drop k es)
  where
    n  = dimVector x
    es = elemsVector x
    _  = typed t x
    


-------------------------- Num Vector Operations --------------------------

prop_shift t k x =
    k `shiftVector` x === mapVector (k+) x
  where
    _ = typed t x

prop_add t (VectorPair x y) =
    x `addVector` y === zipWithVector (+) x y
  where
    _ = typed t x

prop_addWithScale t a b (VectorPair x y) =
    addVectorWithScale a x b y ~== 
        zipWithVector (+) (mapVector (a*) x) (mapVector (b*) y)
  where
    _ = typed t x

prop_sub t (VectorPair x y) =
    x `subVector` y === zipWithVector (-) x y
  where
    _ = typed t x

prop_scale t k x =
    k `scaleVector` x === mapVector (k*) x
  where
    _ = typed t x

prop_negate t x =
    negateVector x === mapVector negate x
  where
    _ = typed t x

-}
