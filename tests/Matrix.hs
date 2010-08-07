module Matrix (
    tests_Matrix
    ) where

import Control.Monad( replicateM, zipWithM_ )
import Data.AEq
import Data.List( transpose )
import Debug.Trace
import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Numeric.LinearAlgebra

import Test.QuickCheck.LinearAlgebra( TestElem(..), Dim2(..), Index2(..),
    Assocs2(..), MatrixPair(..) )
import qualified Test.QuickCheck.LinearAlgebra as Test

import Typed


tests_Matrix = testGroup "Matrix"
    [ testPropertyI "dim/matrix" prop_dim_matrix
    , testPropertyI "at/matrix" prop_at_matrix
    , testPropertyI "listMatrix" prop_listMatrix
    , testPropertyI "colListMatrix" prop_colListMatrix    
    , testPropertyI "rowListMatrix" prop_rowListMatrix        
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
    , testPropertyI "row" prop_row
    , testPropertyI "rows" prop_rows
    , testPropertyI "slice" prop_slice
    , testPropertyI "splitRowsAt" prop_splitRowsAt
    , testPropertyI "splitColsAt" prop_splitColsAt
    , testPropertyDZ "shift" prop_shift prop_shift
    , testPropertyDZ "shiftDiag" prop_shiftDiag prop_shiftDiag
    , testPropertyDZ "shiftDiagWithScale" prop_shiftDiagWithScale prop_shiftDiagWithScale
    , testPropertyDZ "add" prop_add prop_add
    , testPropertyDZ "addWithScales" prop_addWithScales prop_addWithScales
    , testPropertyDZ "sub" prop_sub prop_sub
    , testPropertyDZ "scale" prop_scale prop_scale
    , testPropertyDZ "scaleRows" prop_scaleRows prop_scaleRows
    , testPropertyDZ "scaleCols" prop_scaleCols prop_scaleCols
    , testPropertyDZ "negate" prop_negate prop_negate
    , testPropertyDZ "conj" prop_conj prop_conj
    , testPropertyDZ "trans" prop_trans prop_trans
    , testPropertyDZ "conjTrans" prop_conjTrans prop_conjTrans
    , testPropertyDZ "rank1Update" prop_rank1Update prop_rank1Update
    , testPropertyDZ "mulMatrixVector" prop_mulMatrixVector prop_mulMatrixVector
    , testPropertyDZ "mulMatrixVectorWithScale" prop_mulMatrixVectorWithScale prop_mulMatrixVectorWithScale
    , testPropertyDZ "mulMatrixAddVectorWithScales" prop_mulMatrixAddVectorWithScales prop_mulMatrixAddVectorWithScales    
    , testPropertyDZ "mulMatrixMatrix" prop_mulMatrixMatrix prop_mulMatrixMatrix
    , testPropertyDZ "mulMatrixMatrixWithScale" prop_mulMatrixMatrixWithScale prop_mulMatrixMatrixWithScale
    , testPropertyDZ "mulMatrixAddMatrixWithScales" prop_mulMatrixAddMatrixWithScales prop_mulMatrixAddMatrixWithScales    

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

prop_colListMatrix t (Dim2 (m,n)) =
    forAll (replicateM n $ Test.vector m) $ \cs ->
        colListMatrix (m,n) cs === (typed t $ listMatrix (m,n) $
            concatMap elemsVector cs)

prop_rowListMatrix t (Dim2 (m,n)) =
    forAll (replicateM m $ Test.vector n) $ \rs ->
        rowListMatrix (m,n) rs === (typed t $ listMatrix (m,n) $
            concat $ transpose $ map elemsVector rs)

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

prop_row t (Index2 (m,n) (i,_)) =
    forAll (typed t `fmap` Test.matrix (m,n)) $ \a ->
        rowMatrix a i === listVector n [ atMatrix a (i,j) | j <- [ 0..n-1 ] ]

prop_rows t a =
    rowsMatrix a === [ rowMatrix a i | i <- [ 0..m-1 ] ]
  where
    (m,_) = dimMatrix a
    _ = typed t $ immutableMatrix a

prop_slice t a =
    forAll (choose (0,m)) $ \m' ->
    forAll (choose (0,n)) $ \n' ->
    forAll (choose (0,m-m')) $ \i ->
    forAll (choose (0,n-n')) $ \j ->
        sliceMatrix (i,j) (m',n') a
            === colListMatrix (m',n') [ sliceVector i m' (colMatrix a j')
                                      | j' <- [ j..j+n'-1 ] ]
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_splitRowsAt t a =
    forAll (choose (0,m)) $ \i ->
        splitRowsMatrixAt i a
            === ( sliceMatrix (0,0) (i,n) a
                , sliceMatrix (i,0) (m-i,n) a
                )
  where
    (m,n) = dimMatrix a
    _  = typed t $ immutableMatrix a

prop_splitColsAt t a =
    forAll (choose (0,n)) $ \j ->
        splitColsMatrixAt j a
            === ( sliceMatrix (0,0) (m,j) a
                , sliceMatrix (0,j) (m,n-j) a
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

prop_addWithScales t a b (MatrixPair x y) =
    addMatrixWithScales a x b y ~== 
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

prop_conj t x =
    conjMatrix x === mapMatrix conj x
  where
    _ = typed t x


-------------------------- Linear Algebra --------------------------

prop_trans t a =
    transMatrix a
        ===
        matrix (swap $ dimMatrix a) [ (swap ij, e) | (ij,e) <- assocsMatrix a ]
  where
    swap (i,j) = (j,i)
    _ = typed t a
    
prop_conjTrans t a =
    conjTransMatrix a === conjMatrix (transMatrix a)
  where
    _ = typed t a

prop_rank1Update t alpha a =
    forAll (Test.vector m) $ \x ->
    forAll (Test.vector n) $ \y -> let y' = conjVector y in
        rank1UpdateMatrix alpha x y a
            ~==
            matrix (m,n) [ ((i,j), alpha * atVector x i * atVector y' j + e)
                         | ((i,j),e) <- assocsMatrix a
                         ]
  where
    (m,n)= dimMatrix a
    _ = typed t a

data MulMatrixAddVector e =
    MulMatrixAddVector Trans (Matrix e) (Vector e) (Vector e) deriving (Show)
    
instance (Storable e, Arbitrary e) => Arbitrary (MulMatrixAddVector e) where
    arbitrary = do
        transa <- arbitrary
        a <- arbitrary
        let (ma,na) = dimMatrix a
            (m,n) = case transa of NoTrans -> (ma,na)
                                   _       -> (na,ma)
        x <- Test.vector n
        y <- Test.vector m
        return $ MulMatrixAddVector transa a x y
        
data MulMatrixVector e =
    MulMatrixVector Trans (Matrix e) (Vector e) deriving (Show)
instance (Storable e, Arbitrary e) => Arbitrary (MulMatrixVector e) where
    arbitrary = do
        (MulMatrixAddVector transa a x _) <- arbitrary
        return $ MulMatrixVector transa a x
        
prop_mulMatrixVector t (MulMatrixVector transa a x) =
    mulMatrixVector transa a x
        ~==
        case transa of
            NoTrans   -> listVector (fst $ dimMatrix a)
                                    [ dotVector x (conjVector r)
                                    | r <- rowsMatrix a ]

            Trans     -> listVector (snd $ dimMatrix a)
                                    [ dotVector x (conjVector c)
                                    | c <- colsMatrix a ]
                                    
            ConjTrans -> listVector (snd $ dimMatrix a)
                                    [ dotVector x c
                                    | c <- colsMatrix a ]
  where
    _ = typed t a

prop_mulMatrixVectorWithScale t alpha (MulMatrixVector transa a x) =
    mulMatrixVectorWithScale alpha transa a x
        ~==
        mulMatrixVector transa a (scaleVector alpha x)
  where
    _ = typed t a

prop_mulMatrixAddVectorWithScales t alpha beta (MulMatrixAddVector transa a x y) =
    mulMatrixAddVectorWithScales alpha transa a x beta y
        ~==
        addVector (mulMatrixVectorWithScale alpha transa a x)
                  (scaleVector beta y)
  where
    _ = typed t a

data MulMatrixAddMatrix e =
    MulMatrixAddMatrix Trans (Matrix e) Trans (Matrix e) (Matrix e) deriving (Show)
    
instance (Storable e, Arbitrary e) => Arbitrary (MulMatrixAddMatrix e) where
    arbitrary = do
        transa <- arbitrary
        transb <- arbitrary
        c <- arbitrary
        k <- fst `fmap` Test.dim2
        
        let (m,n) = dimMatrix c
            (ma,na) = case transa of NoTrans -> (m,k)
                                     _       -> (k,m)
            (mb,nb) = case transb of NoTrans -> (k,n)
                                     _       -> (n,k)
        a <- Test.matrix (ma,na)
        b <- Test.matrix (mb,nb)
        
        return $ MulMatrixAddMatrix transa a transb b c

data MulMatrixMatrix e =
    MulMatrixMatrix Trans (Matrix e) Trans (Matrix e) deriving (Show)
instance (Storable e, Arbitrary e) => Arbitrary (MulMatrixMatrix e) where
    arbitrary = do
        (MulMatrixAddMatrix transa a transb b _) <- arbitrary
        return $ MulMatrixMatrix transa a transb b

prop_mulMatrixMatrix t (MulMatrixMatrix transa a transb b) =
    mulMatrixMatrix transa a transb b
        ~==
        colListMatrix (m,n) [ mulMatrixVector transa a x | x <- colsMatrix b' ]
  where
    m = case transa of NoTrans -> (fst $ dimMatrix a)
                       _       -> (snd $ dimMatrix a)
    n = case transb of NoTrans -> (snd $ dimMatrix b)
                       _       -> (fst $ dimMatrix b)
    b' = case transb of NoTrans   -> b
                        Trans     -> transMatrix b
                        ConjTrans -> conjTransMatrix b
    _ = typed t a

prop_mulMatrixMatrixWithScale t alpha (MulMatrixMatrix transa a transb b) =
    mulMatrixMatrixWithScale alpha transa a transb b
        ~==
        scaleMatrix alpha (mulMatrixMatrix transa a transb b)
  where
    _ = typed t a

prop_mulMatrixAddMatrixWithScales t alpha beta (MulMatrixAddMatrix transa a transb b c) =
    mulMatrixAddMatrixWithScales alpha transa a transb b beta c
        ~==
        addMatrixWithScales alpha (mulMatrixMatrix transa a transb b)
                            beta c
  where
    _ = typed t a



testAEq a b =
    if a ~== b then True
               else trace ("expected: " ++ show b ++ "\nactual: " ++ show a) False
