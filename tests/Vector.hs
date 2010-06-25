module Vector (
    tests_Vector
    ) where

import Data.AEq
import BLAS.Elem
import BLAS.Vector

import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck hiding ( vector )
import Test.QuickCheck.BLAS( Dim(..), Assocs(..), VectorPair(..) )
import qualified Test.QuickCheck.BLAS as Test


tests_Vector = testGroup "Vector"
    [ testPropertyI "dim/vector" prop_dim_vector
    , testPropertyI "at/vector" prop_at_vector
    , testPropertyI "listVector" prop_listVector
    , testPropertyI "constantVector" prop_constantVector
    , testPropertyI "indices" prop_indices
    , testPropertyI "elems" prop_elems
    , testPropertyI "assocs" prop_assocs
    , testPropertyI "replace" prop_replace
    , testPropertyI "accum" prop_accum
    , testPropertyI "map" prop_map
    , testPropertyI "zipWith" prop_zipWith
    , testPropertyI "splice" prop_splice
    , testPropertyI "splitAt" prop_splitAt
    , testPropertyDZ "(+)" prop_plus prop_plus
    , testPropertyDZ "(*)" prop_times prop_times
    , testPropertyDZ "(-)" prop_minus prop_minus
    , testPropertyDZ "negate" prop_negate prop_negate
    , testPropertyDZ "abs" prop_abs prop_abs
    , testPropertyDZ "signum" prop_signum prop_signum
    , testPropertyDZ "fromInteger" prop_fromInteger prop_fromInteger
    , testPropertyDZ "(/)" prop_div prop_div
    , testPropertyDZ "(recip)" prop_recip prop_recip
    , testPropertyDZ "fromRational" prop_fromRational prop_fromRational
    , testPropertyDZ "pi" prop_pi prop_pi
    , testPropertyDZ "exp" prop_exp prop_exp
    , testPropertyDZ "sqrt" prop_sqrt prop_sqrt
    , testPropertyDZ "log" prop_log prop_log
    , testPropertyDZ "(**)" prop_pow prop_pow
    , testPropertyDZ "logBase" prop_logBase prop_logBase
    , testPropertyDZ "sin" prop_sin prop_sin
    , testPropertyDZ "tan" prop_tan prop_tan
    , testPropertyDZ "cos" prop_cos prop_cos
    , testPropertyDZ "asin" prop_asin prop_asin    
    , testPropertyDZ "atan" prop_atan prop_atan    
    , testPropertyDZ "acos" prop_acos prop_acos    
    , testPropertyDZ "sinh" prop_sinh prop_sinh
    , testPropertyDZ "tanh" prop_tanh prop_tanh
    , testPropertyDZ "cosh" prop_cosh prop_cosh
    , testPropertyDZ "asinh" prop_asinh prop_asinh
    , testPropertyDZ "atanh" prop_atanh prop_atanh
    , testPropertyDZ "acosh" prop_acosh prop_acosh
    ]

typed :: e -> a e -> a e
typed _ = id

immutableVector :: Vector e -> Vector e
immutableVector = id


testPropertyI :: (Testable a)
              => TestName
              -> (Int -> a)
              -> Test
testPropertyI str prop =
    testProperty str $ prop undefined

testPropertyDZ :: (Testable a, Testable b)
               => TestName
               -> (Double -> a)
               -> (Complex Double -> b)
               -> Test
testPropertyDZ str propd propz =
    testGroup str
        [ testProperty "Double" $ propd undefined
        , testProperty "Complex Double" $ propz undefined
        ]
    

------------------------- Vector Construction ------------------------------

prop_dim_vector t (Assocs n ies) =
    dimVector (vector n ies) == n
  where
    _ = typed t $ vector n ies

prop_at_vector t (Assocs n ies) = let
    x = vector n ies
    is = (fst . unzip) ies
    in and [ atVector x i `elem` [ e | (i',e) <- ies, i' == i ]
           | i <- is]
  where
    _ = typed t $ vector n ies

prop_listVector t (Dim n) es =
    listVector n es == vector n (zip [ 0..n-1 ] es)
  where
    _ = typed t es

prop_constantVector t (Dim n) e =
    constantVector n e == listVector n (replicate n e)
  where
    _ = typed t [e]


-------------------------- Accessing Vectors ------------------------------

prop_indices t x =
    indicesVector x == [ 0..((dimVector x) - 1) ]
  where
    _ = immutableVector x
    _ = typed t x

prop_elems t x =
    elemsVector x == [ atVector x i | i <- indicesVector x ]
  where
    _ = typed t x
    
prop_assocs t x =
    assocsVector x == zip (indicesVector x) (elemsVector x)
  where
    _ = typed t x


------------------------- Incremental Updates ------------------------------
    
prop_replace t (Assocs n ies) =
    forAll (typed t `fmap` Test.vector n) $ \x -> let
        x' = replaceVector x ies
        is = indicesVector x
        is1 = (fst . unzip) ies
        is0 = [ i | i <- is, i `notElem` is1 ]
        in and $
            [ atVector x' i `elem` [ e | (i',e) <- ies, i' == i ]
            | i <- is1
            ] ++
            [ atVector x' i == atVector x i
            | i <- is0
            ]

prop_accum t (Blind f) (Assocs n ies) =
    forAll (typed t `fmap` Test.vector n) $ \x -> let
        x' = accumVector f x ies
        in x' == listVector n [ foldl f e [ e' | (i',e') <- ies, i' == i]
                              | (i,e) <- assocsVector x ]
  where
      _ = typed t $ (snd . unzip) ies
     
     
-------------------------- Derived Vectors ------------------------------
     
prop_map t (Blind f) x =
    mapVector f x == listVector (dimVector x) (map f $ elemsVector x)
  where
    _ = typed t x
    _ = typed t $ mapVector f x

prop_zipWith t (Blind f) (VectorPair x y) =
    zipWithVector f x y == (listVector (dimVector x) $
                                zipWith f (elemsVector x) (elemsVector y))
  where
    _ = typed t x
    _ = typed t y    
    _ = typed t $ zipWithVector f x y


------------------------------ Vector Views-- --------------------------------

prop_splice t x = 
    forAll (choose (0,n)) $ \n' ->
    forAll (choose (0,n-n')) $ \o ->
        spliceVector x o n' == listVector n' (take n' $ drop o $ es)
  where
    n  = dimVector x
    es = elemsVector x
    _  = typed t x

prop_splitAt t x =
    forAll (choose (0,n)) $ \k ->
        splitVectorAt k x == (listVector k $ take k es,
                              listVector (n-k) $ drop k es)
  where
    n  = dimVector x
    es = elemsVector x
    _  = typed t x
    


-------------------------- Num Vector Operations --------------------------

prop_plus t (VectorPair x y) =
    x + y === zipWithVector (+) x y
  where
    _ = typed t x

prop_times t (VectorPair x y) =
    x * y === zipWithVector (*) x y
  where
    _ = typed t x

prop_minus t (VectorPair x y) =
    x - y === zipWithVector (-) x y
  where
    _ = typed t x

prop_negate t x =
    negate x === mapVector negate x
  where
    _ = typed t x

prop_abs t x =
    abs x === mapVector abs x
  where
    _ = typed t x

prop_signum t x =
    signum x == mapVector signum x
  where
    _ = typed t x

prop_fromInteger t n =
    fromInteger n == typed t (listVector 1 [fromInteger n])


---------------------- Fractional Vector Operations ------------------------

prop_div t (VectorPair x y) =
    x / y ~== zipWithVector (/) x y
  where
    _ = typed t x
    
prop_recip t x =
    recip x === constantVector (dimVector x) 1 / x
  where
    _ = typed t x
    
prop_fromRational t q =
    fromRational q == typed t (listVector 1 [fromRational q])


---------------------- Floating Vector Operations ------------------------

prop_pi t =
    pi == typed t (listVector 1 [pi])

prop_exp t x =
    exp x ~== mapVector exp x
  where
    _ = typed t x

prop_sqrt t x =
    sqrt x ~== mapVector sqrt x
  where
    _ = typed t x

prop_log t x =
    log x ~== mapVector log x
  where
    _ = typed t x

prop_pow t (VectorPair x y) =
    x ** y ~== zipWithVector (**) x y
  where
    _ = typed t x

prop_logBase t (VectorPair x y) =
    logBase x y ~== zipWithVector logBase x y
  where
    _ = typed t x

prop_sin t x =
    sin x ~== mapVector sin x
  where
    _ = typed t x

prop_tan t x =
    tan x ~== mapVector tan x
  where
    _ = typed t x

prop_cos t x =
    cos x ~== mapVector cos x
  where
    _ = typed t x

prop_asin t x =
    asin x ~== mapVector asin x
  where
    _ = typed t x

prop_atan t x =
    atan x ~== mapVector atan x
  where
    _ = typed t x

prop_acos t x =
    acos x ~== mapVector acos x
  where
    _ = typed t x

prop_sinh t x =
    sin x ~== mapVector sinh x
  where
    _ = typed t x

prop_tanh t x =
    tan x ~== mapVector tanh x
  where
    _ = typed t x

prop_cosh t x =
    cos x ~== mapVector cosh x
  where
    _ = typed t x

prop_asinh t x =
    asin x ~== mapVector asinh x
  where
    _ = typed t x

prop_atanh t x =
    atan x ~== mapVector atanh x
  where
    _ = typed t x

prop_acosh t x =
    acos x ~== mapVector acosh x
  where
    _ = typed t x





{-


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


-}
