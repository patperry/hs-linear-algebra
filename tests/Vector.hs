module Vector (
    tests_Vector
    ) where

import Data.AEq
import Debug.Trace
import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Numeric.LinearAlgebra.Elem
import Numeric.LinearAlgebra.Vector

import Test.QuickCheck.LinearAlgebra( TestElem(..), Dim(..), Assocs(..),
    VectorPair(..) )
import qualified Test.QuickCheck.LinearAlgebra as Test

import Typed


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
    , testPropertyDZ "sum" prop_sum prop_sum
    , testPropertyDZ "sumAbs" prop_sumAbs prop_sumAbs    
    , testPropertyDZ "norm2" prop_norm2 prop_norm2
    , testPropertyDZ "whichMaxAbs1" prop_whichMaxAbs1 prop_whichMaxAbs1
    , testPropertyDZ "whichMaxAbs2" prop_whichMaxAbs1 prop_whichMaxAbs2
    , testPropertyDZ "dot" prop_dot prop_dot
    , testPropertyDZ "kronecker" prop_kronecker prop_kronecker
    , testPropertyDZ "shift" prop_shift prop_shift
    , testPropertyDZ "add" prop_add prop_add
    , testPropertyDZ "addWithScale" prop_addWithScale prop_addWithScale
    , testPropertyDZ "sub" prop_sub prop_sub
    , testPropertyDZ "scale" prop_scale prop_scale
    , testPropertyDZ "mul" prop_mul prop_mul
    , testPropertyDZ "negate" prop_negate prop_negate
    , testPropertyDZ "conj" prop_conj prop_conj
    , testPropertyDZ "abs" prop_abs prop_abs
    , testPropertyDZ "signum" prop_signum prop_signum
    , testPropertyDZ "div" prop_div prop_div
    , testPropertyDZ "recip" prop_recip prop_recip
    , testPropertyDZ "sqrt" prop_sqrt prop_sqrt
    , testPropertyDZ "exp" prop_exp prop_exp
    , testPropertyDZ "log" prop_log prop_log
    , testPropertyDZ "pow" prop_pow prop_pow
    , testPropertyDZ "sin" prop_sin prop_sin
    , testPropertyDZ "cos" prop_cos prop_cos
    , testPropertyDZ "tan" prop_tan prop_tan
    , testPropertyDZ "asin" prop_asin prop_asin    
    , testPropertyDZ "acos" prop_acos prop_acos    
    , testPropertyDZ "atan" prop_atan prop_atan    
    , testPropertyDZ "sinh" prop_sinh prop_sinh
    , testPropertyDZ "cosh" prop_cosh prop_cosh
    , testPropertyDZ "tanh" prop_tanh prop_tanh
    , testPropertyDZ "asinh" prop_asinh prop_asinh
    , testPropertyDZ "acosh" prop_acosh prop_acosh
    , testPropertyDZ "atanh" prop_atanh prop_atanh
    ]



------------------------- Vector Construction ------------------------------

prop_dim_vector t (Assocs n ies) =
    dimVector (vector n ies) === n
  where
    _ = typed t $ vector n ies

prop_at_vector t (Assocs n ies) = let
    x = vector n ies
    is = (fst . unzip) ies
    in and [ atVector x i `elem` [ e | (i',e) <- ies, i' == i ]
           | i <- is]
  where
    _ = typed t $ vector n ies

prop_listVector t (Dim n) =
    forAll (QC.vector n) $ \es ->
        listVector n es === (typed t $ vector n $ zip [ 0..n-1 ] es)

prop_constantVector t (Dim n) e =
    constantVector n e === listVector n (replicate n e)
  where
    _ = typed t [e]


-------------------------- Accessing Vectors ------------------------------

prop_indices t x =
    indicesVector x === [ 0..((dimVector x) - 1) ]
  where
    _ = immutableVector x
    _ = typed t x

prop_elems t x =
    elemsVector x === [ atVector x i | i <- indicesVector x ]
  where
    _ = typed t x
    
prop_assocs t x =
    assocsVector x === zip (indicesVector x) (elemsVector x)
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
            [ atVector x' i === atVector x i
            | i <- is0
            ]

prop_accum t (Blind f) (Assocs n ies) =
    forAll (typed t `fmap` Test.vector n) $ \x -> let
        x' = accumVector f x ies
        in x' === listVector n [ foldl f e [ e' | (i',e') <- ies, i' == i]
                               | (i,e) <- assocsVector x ]
  where
      _ = typed t $ (snd . unzip) ies
     
     
-------------------------- Derived Vectors ------------------------------
     
prop_map t (Blind f) x =
    mapVector f x === listVector (dimVector x) (map f $ elemsVector x)
  where
    _ = typed t x
    _ = typed t $ mapVector f x

prop_zipWith t (Blind f) (VectorPair x y) =
    zipWithVector f x y === (listVector (dimVector x) $
                                zipWith f (elemsVector x) (elemsVector y))
  where
    _ = typed t x
    _ = typed t y    
    _ = typed t $ zipWithVector f x y


------------------------------ Vector Views-- --------------------------------

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

prop_mul t (VectorPair x y) =
    x `mulVector` y === zipWithVector (*) x y
  where
    _ = typed t x

prop_negate t x =
    negateVector x === mapVector negate x
  where
    _ = typed t x

prop_conj t x =
    conjVector x === mapVector conj x
  where
    _ = typed t x

prop_abs t x =
    absVector x === mapVector abs x
  where
    _ = typed t x

prop_signum t x =
    signumVector x === mapVector signum x
  where
    _ = typed t x


---------------------- Fractional Vector Operations ------------------------

prop_div t (VectorPair x y) =
    x `divVector` y ~== zipWithVector (/) x y
  where
    _ = typed t x
    
prop_recip t x =
    recipVector x ~== mapVector (1/) x
  where
    _ = typed t x
    

---------------------- Floating Vector Operations ------------------------

prop_exp t x =
    expVector x ~== mapVector exp x
  where
    _ = typed t x

prop_sqrt t x =
    sqrtVector x ~== mapVector sqrt x
  where
    _ = typed t x

prop_log t x =
    logVector x ~== mapVector log x
  where
    _ = typed t x

prop_pow t (VectorPair x y) =
    x `powVector` y ~== zipWithVector (**) x y
  where
    _ = typed t x

prop_sin t x =
    sinVector x ~== mapVector sin x
  where
    _ = typed t x

prop_cos t x =
    cosVector x ~== mapVector cos x
  where
    _ = typed t x

prop_tan t x =
    tanVector x ~== mapVector tan x
  where
    _ = typed t x

prop_asin t x =
    -- trace (show (asinVector x) ++ "\n" ++ (show $ mapVector asin x)) $    
    asinVector x ~== mapVector asin x
  where
    _ = typed t x

prop_acos t x =
    acosVector x ~== mapVector acos x
  where
    _ = typed t x

prop_atan t x =
    atanVector x ~== mapVector atan x
  where
    _ = typed t x

prop_sinh t x =
    sinhVector x ~== mapVector sinh x
  where
    _ = typed t x

prop_cosh t x =
    coshVector x ~== mapVector cosh x
  where
    _ = typed t x

prop_tanh t x =
    tanhVector x ~== mapVector tanh x
  where
    _ = typed t x

prop_asinh t x =
    -- trace (show (asinhVector x) ++ "\n" ++ (show $ mapVector asinh x)) $
    asinhVector x ~== mapVector asinh x
  where
    _ = typed t x

prop_acosh t x =
    acoshVector x ~== mapVector acosh x
  where
    _ = typed t x

prop_atanh t x =
    atanhVector x ~== mapVector atanh x
  where
    _ = typed t x


-------------------------- Vector Properties ---------------------------------

prop_sum t x =
    sumVector x ~== (sum $ elemsVector x)
  where
    _ = typed t x

prop_sumAbs t x =
    sumAbsVector x ~== (sum $ map norm1 $ elemsVector x)
  where
    _ = typed t x

prop_norm2 t x =
    norm2Vector x ~== (sqrt $ sum $ map (^^2) $ map norm $ elemsVector x)
  where
    _ = typed t x

prop_whichMaxAbs1 t x =
    (dimVector x > 0) && all (not . isNaN) (map norm1 $ elemsVector x) ==>
        atVector x i === e
  where
    (i,e) = whichMaxAbsVector x      
    _     = typed t x

prop_whichMaxAbs2 t x =
    (dimVector x > 0) && all (not . isNaN) (map norm1 $ elemsVector x) ==>
        all (<= norm1 e) $ map norm1 (elemsVector x)
  where
    (_,e) = whichMaxAbsVector x
    _     = typed t x

prop_dot t (VectorPair x y) =
    dotVector x y ~== sum (x * conj y)
  where
    sum  = sumVector
    conj = conjVector
    (*)  = mulVector
    _    = typed t x

prop_kronecker t x y =
    x `kroneckerVector` y ===
        listVector (dimVector x * dimVector y)
                   [ e*f | e <- elemsVector x, f <- elemsVector y ]
  where
    _ = typed t x
