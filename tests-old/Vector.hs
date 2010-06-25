{-# LANGUAGE ScopedTypeVariables #-}
module Vector
    where

import Data.Elem.BLAS
import Data.Vector.Dense
import Data.Tensor.Class


import Driver
import qualified Data.Array as Array
import Test.Vector.Dense hiding ( vector )
import qualified Test.Vector.Dense as Test
import System.IO.Unsafe

type V = Vector E


---------------------------- Creating Vectors --------------------------------

prop_vector_dim (Assocs n ies) =
    dim (vector n ies :: V) == n
prop_vector_assocs (Assocs n ies) =
    (zip [0..(n-1)] (repeat 0) ++ ies) 
    `assocsEq` 
    (assocs (vector n ies :: V))

prop_listVector_dim es =
    let n = length es
    in dim (listVector n es :: V) == n
prop_listVector_assocs es =
    let n = length es
    in assocs (listVector n es :: V) `assocsEq` (zip [0..] es)


---------------------- Reading and Writing Elements --------------------------

prop_dim (x :: V) = 
    dim x == length (elems x)
prop_bounds (x :: V) =
    bounds x == (0, dim x - 1)
prop_at (Index n i) =
    forAll (Test.vector n) $ \(x :: V) -> 
        (x!i) === ((elems x) !! i)

prop_indices (x :: V) =
    indices x == [0..(dim x - 1)]
prop_assocs (x :: V) = 
    assocs x === zip (indices x) (elems x)

prop_replace_elems (Assocs n ies) =
    forAll (Test.vector n) $ \x ->
        let x'   = x // ies
            ies' = Array.assocs $ (Array.//) (Array.array (0,n-1) $ assocs x) ies
        in ies' `assocsEq` assocs x'

------------------------------ Vector Views-- --------------------------------

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
        

----------------------------- Special Vectors --------------------------------

prop_zeroVector_dim (Nat n) =
    dim (zeroVector n :: V) == n
    
prop_zeroVector_elems (Nat n) =
    elems (zeroVector n :: V) == replicate n 0
    
prop_constantVector_dim (Nat n) e =
    dim (constantVector n e :: V) == n
    
prop_constantVector_elems (Nat n) (e :: E) =
    elems (constantVector n e :: V) === replicate n e

prop_basisVector_dim (Index n i) =
    dim (basisVector n i :: V) == n
    
prop_basisVector_elems (Index n i) =
    elems (basisVector n i :: V) == (replicate i 0) ++ [1] ++ (replicate (n-i-1) 0)


-------------------------- Unsary Vector Operations --------------------------

prop_shift k (x :: V) =
    elems (shift k x) ~== map (k+) (elems x)

prop_scale k (x :: V) =
    elems (k *> x) ~== map (k*) (elems x)

prop_conj_elems (x :: V) =
    and $ zipWith (===) (elems $ conj x) (map conjugate $ elems x)

prop_conj_scale k (x :: V) =
    conj (k *> x) ===  (conjugate k *> (conj x))

prop_negate (x :: V) =
    negate x ~== (-1) *> x

prop_abs (x :: V) =
    elems (abs x) ~== map abs (elems x)

prop_signum (x :: V) =
    elems (signum x) ~== map signum (elems x)

prop_recip (x :: V) =
    elems (recip x) ~== (map recip $ elems x)


------------------------- Binary Vector Operations ---------------------------

prop_plus (VectorPair (x :: V) y) =
    elems (x + y) ~== zipWith (+) (elems x) (elems y)
    
prop_minus (VectorPair (x :: V) y) =
    elems (x - y) ~== zipWith (-) (elems x) (elems y)
    
prop_times (VectorPair (x :: V) y) =
    elems (x * y) ~== zipWith (*) (elems x) (elems y)
    
prop_divide (VectorPair (x :: V) y) =
    elems (x / y) ~== zipWith (/) (elems x) (elems y)


-------------------------- Vector Properties ---------------------------------

prop_sumAbs (x :: V) =
    sumAbs x ~== (sum $ map norm1 $ elems x)
    
prop_norm2 (x :: V) =
    norm2 x ~== (sqrt $ sum $ map (^2) $ map norm $ elems x)
    
prop_whichMaxAbs1 (x :: V) =
    (dim x > 0) && all (not . isNaN) (map norm1 $ elems x) ==>
        let (i,e) = whichMaxAbs x
        in x ! i === e
        
prop_whichMaxAbs2 (x :: V) =
    (dim x > 0) && all (not . isNaN) (map norm1 $ elems x) ==>
        let a = norm1 $ snd $ whichMaxAbs x
        in all (<= a) (map norm1 $ elems x)
        
prop_dot_self (x :: V) =
    (sqrt $ x <.> x) ~== (fromReal $ norm2 x)
    
prop_dot_conj (VectorPair (x :: V) y) =
    (x <.> y) ~== (conjugate $ y <.> x)
    
prop_dot_scale1 k (VectorPair (x :: V) y) =
    (x <.> (k *> y)) ~== k * (x <.> y)
    
prop_dot_scale2 k (VectorPair (x :: V) y) =
    ((k *> x) <.> y) ~== (conjugate k) * (x <.> y)
    
prop_dot_linear1 (VectorTriple (x :: V) y z) =
    (x <.> (y + z)) ~== (x <.> y + x <.> z)
    
prop_dot_linear2 (VectorTriple (x :: V) y z) =
    ((x + y) <.> z) ~== (x <.> z + y <.> z)

------------------------------------------------------------------------------
tests_Vector =
    [ testProperty "vector/dim" prop_vector_dim
    , testProperty "vector/assocs" prop_vector_assocs
    , testProperty "listVector/dim" prop_listVector_dim
    , testProperty "listVector/assocs" prop_listVector_assocs
    
    , testProperty "dim" prop_dim
    , testProperty "bounds" prop_bounds
    , testProperty "(!)" prop_at
    , testProperty "indices" prop_indices
    , testProperty "assocs" prop_assocs
    , testProperty "(//)" prop_replace_elems
    
    , testProperty "subvector/dim" prop_subvector_dim
    , testProperty "subvector/elems" prop_subvector_elems
    , testProperty "subvectorWithStride/dim" prop_subvectorWithStride_dim
    , testProperty "subvectorWithStride/elems" prop_subvectorWithStride_elems
    
    , testProperty "zeroVector/dim" prop_zeroVector_dim
    , testProperty "zeroVector/elems" prop_zeroVector_elems
    , testProperty "constantVector/dim" prop_constantVector_dim
    , testProperty "constantVector/elems" prop_constantVector_elems
    , testProperty "basisVector/dim" prop_basisVector_dim
    , testProperty "basisVector/elems" prop_basisVector_elems

    , testProperty "conj" prop_conj_elems
    , testProperty "(*>)" prop_scale
    , testProperty "shift" prop_shift
    , testProperty "conj . (*>)" prop_conj_scale
    , testProperty "negate" prop_negate
    , testProperty "abs" prop_abs
    , testProperty "signum" prop_signum
    , testProperty "recip" prop_recip
    
    , testProperty "(+)" prop_plus
    , testProperty "(-)" prop_minus
    , testProperty "(*)" prop_times
    , testProperty "(/)" prop_divide
    
    , testProperty "sumAbs" prop_sumAbs
    , testProperty "norm2" prop_norm2
    , testProperty "whichMaxAbs1" prop_whichMaxAbs1
    , testProperty "whichMaxAbs2" prop_whichMaxAbs2
    , testProperty "dot self" prop_dot_self
    , testProperty "dot conj" prop_dot_conj
    , testProperty "dot scale1" prop_dot_scale1
    , testProperty "dot scale2" prop_dot_scale2
    , testProperty "dot linear1" prop_dot_linear1
    , testProperty "dot linear2" prop_dot_linear2
    
    ]



assocsEq :: [(Int,E)] -> [(Int,E)] -> Bool
assocsEq ies ies' = 
    if (ordered ies === ordered ies')
        then True
        else unsafePerformIO $ do
            zipWithM_ (\(i1,e1) (i2,e2) -> do
                putStr $ show i1 ++ ": " ++ show e1 ++ " " ++ show e2
                unless (e1 ~== e2) $ putStr " **** "
                putStrLn ""
                )
                (ordered ies)
                (ordered ies')
            return False
  where
    ordered = sortAssocs . nubAssocs
    nubAssocs = reverse . nubBy ((==) `on` fst) . reverse      
    sortAssocs = sortBy (comparing fst)

