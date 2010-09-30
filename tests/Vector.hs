module Vector (
    tests_Vector
    ) where

import Data.AEq
import Debug.Trace
import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Numeric.LinearAlgebra
import qualified Numeric.LinearAlgebra.Vector as V

import Test.QuickCheck.LinearAlgebra( TestElem(..), Dim(..), Assocs(..),
    VectorPair(..) )
import qualified Test.QuickCheck.LinearAlgebra as Test

import Typed


tests_Vector = testGroup "Vector"
    [ testPropertyI "dim/fromAssocs" prop_dim_fromAssocs
    , testPropertyI "at/fromAssocs" prop_at_fromAssocs
    , testPropertyI "fromList" prop_fromList
    , testPropertyI "constant" prop_constant
    , testPropertyI "indices" prop_indices
    , testPropertyI "elems" prop_elems
    , testPropertyI "assocs" prop_assocs
    , testPropertyI "replace" prop_replace
    , testPropertyI "accum" prop_accum
    , testPropertyI "map" prop_map
    , testPropertyI "zipWith" prop_zipWith
    , testPropertyI "concat" prop_concat
    , testPropertyI "slice" prop_slice
    , testPropertyI "splitAt" prop_splitAt
    , testPropertyDZ "sumAbs" prop_sumAbs prop_sumAbs    
    , testPropertyDZ "norm2" prop_norm2 prop_norm2
    , testPropertyDZ "whichMaxAbs1" prop_whichMaxAbs1 prop_whichMaxAbs1
    , testPropertyDZ "whichMaxAbs2" prop_whichMaxAbs1 prop_whichMaxAbs2
    , testPropertyDZ "dot" prop_dot prop_dot
    , testPropertyDZ "kronecker" prop_kronecker prop_kronecker
    , testPropertyDZ "shift" prop_shift prop_shift
    , testPropertyDZ "add" prop_add prop_add
    , testPropertyDZ "addWithScales" prop_addWithScales prop_addWithScales
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

prop_dim_fromAssocs t (Assocs n ies) =
    V.dim (V.fromAssocs n ies) === n
  where
    _ = typed t $ V.fromAssocs n ies

prop_at_fromAssocs t (Assocs n ies) = let
    x = V.fromAssocs n ies
    is = (fst . unzip) ies
    in and [ V.at x i `elem` [ e | (i',e) <- ies, i' == i ]
           | i <- is]
  where
    _ = typed t $ V.fromAssocs n ies

prop_fromList t (Dim n) =
    forAll (QC.vector n) $ \es ->
        V.fromList n es === (typed t $ V.fromAssocs n $ zip [ 0..n-1 ] es)

prop_constant t (Dim n) e =
    V.constant n e === V.fromList n (replicate n e)
  where
    _ = typed t [e]


-------------------------- Accessing Vectors ------------------------------

prop_indices t x =
    V.indices x === [ 0..((V.dim x) - 1) ]
  where
    _ = immutableVector x
    _ = typed t x

prop_elems t x =
    V.elems x === [ V.at x i | i <- V.indices x ]
  where
    _ = typed t x
    
prop_assocs t x =
    V.assocs x === zip (V.indices x) (V.elems x)
  where
    _ = typed t x


------------------------- Incremental Updates ------------------------------
    
prop_replace t (Assocs n ies) =
    forAll (typed t `fmap` Test.vector n) $ \x -> let
        x' = V.replace x ies
        is = V.indices x
        is1 = (fst . unzip) ies
        is0 = [ i | i <- is, i `notElem` is1 ]
        in and $
            [ V.at x' i `elem` [ e | (i',e) <- ies, i' == i ]
            | i <- is1
            ] ++
            [ V.at x' i === V.at x i
            | i <- is0
            ]

prop_accum t (Blind f) (Assocs n ies) =
    forAll (typed t `fmap` Test.vector n) $ \x -> let
        x' = V.accum f x ies
        in x' === V.fromList n [ foldl f e [ e' | (i',e') <- ies, i' == i]
                               | (i,e) <- V.assocs x ]
  where
      _ = typed t $ (snd . unzip) ies
     
     
-------------------------- Derived Vectors ------------------------------
     
prop_map t (Blind f) x =
    V.map f x === V.fromList (V.dim x) (map f $ V.elems x)
  where
    _ = typed t x
    _ = typed t $ V.map f x

prop_zipWith t (Blind f) (VectorPair x y) =
    V.zipWith f x y === (V.fromList (V.dim x) $
                                zipWith f (V.elems x) (V.elems y))
  where
    _ = typed t x
    _ = typed t y    
    _ = typed t $ V.zipWith f x y

prop_concat t xs =
    V.elems (V.concat xs) === concatMap V.elems xs
  where
    _ = typed t $ head xs

------------------------------ Vector Views-- --------------------------------

prop_slice t x = 
    forAll (choose (0,n)) $ \n' ->
    forAll (choose (0,n-n')) $ \o ->
        V.slice o n' x === V.fromList n' (take n' $ drop o $ es)
  where
    n  = V.dim x
    es = V.elems x
    _  = typed t x

prop_splitAt t x =
    forAll (choose (0,n)) $ \k ->
        V.splitAt k x === (V.fromList k $ take k es,
                           V.fromList (n-k) $ drop k es)
  where
    n  = V.dim x
    es = V.elems x
    _  = typed t x
    


-------------------------- Num Vector Operations --------------------------

prop_shift t k x =
    k `V.shift` x === V.map (k+) x
  where
    _ = typed t x

prop_add t (VectorPair x y) =
    x `V.add` y === V.zipWith (+) x y
  where
    _ = typed t x

prop_addWithScales t a b (VectorPair x y) =
    V.addWithScales a x b y ~== 
        V.zipWith (+) (V.map (a*) x) (V.map (b*) y)
  where
    _ = typed t x

prop_sub t (VectorPair x y) =
    x `V.sub` y === V.zipWith (-) x y
  where
    _ = typed t x

prop_scale t k x =
    k `V.scale` x === V.map (k*) x
  where
    _ = typed t x

prop_mul t (VectorPair x y) =
    x `V.mul` y === V.zipWith (*) x y
  where
    _ = typed t x

prop_negate t x =
    V.negate x === V.map negate x
  where
    _ = typed t x

prop_conj t x =
    V.conj x === V.map conj x
  where
    _ = typed t x

prop_abs t x =
    V.abs x === V.map abs x
  where
    _ = typed t x

prop_signum t x =
    V.signum x === V.map signum x
  where
    _ = typed t x


---------------------- Fractional Vector Operations ------------------------

prop_div t (VectorPair x y) =
    x `V.div` y ~== V.zipWith (/) x y
  where
    _ = typed t x
    
prop_recip t x =
    V.recip x ~== V.map (1/) x
  where
    _ = typed t x
    

---------------------- Floating Vector Operations ------------------------

prop_exp t x =
    V.exp x ~== V.map exp x
  where
    _ = typed t x

prop_sqrt t x =
    V.sqrt x ~== V.map sqrt x
  where
    _ = typed t x

prop_log t x =
    V.log x ~== V.map log x
  where
    _ = typed t x

prop_pow t (VectorPair x y) =
    x `V.pow` y ~== V.zipWith (**) x y
  where
    _ = typed t x

prop_sin t x =
    V.sin x ~== V.map sin x
  where
    _ = typed t x

prop_cos t x =
    V.cos x ~== V.map cos x
  where
    _ = typed t x

prop_tan t x =
    V.tan x ~== V.map tan x
  where
    _ = typed t x

prop_asin t x =
    -- trace (show (V.asin x) ++ "\n" ++ (show $ V.map asin x)) $    
    V.asin x ~== V.map asin x
  where
    _ = typed t x

prop_acos t x =
    V.acos x ~== V.map acos x
  where
    _ = typed t x

prop_atan t x =
    V.atan x ~== V.map atan x
  where
    _ = typed t x

prop_sinh t x =
    V.sinh x ~== V.map sinh x
  where
    _ = typed t x

prop_cosh t x =
    V.cosh x ~== V.map cosh x
  where
    _ = typed t x

prop_tanh t x =
    V.tanh x ~== V.map tanh x
  where
    _ = typed t x

prop_asinh t x =
    -- trace (show (V.asinh x) ++ "\n" ++ (show $ V.map asinh x)) $
    V.asinh x ~== V.map asinh x
  where
    _ = typed t x

prop_acosh t x =
    V.acosh x ~== V.map acosh x
  where
    _ = typed t x

prop_atanh t x =
    V.atanh x ~== V.map atanh x
  where
    _ = typed t x


-------------------------- Vector Properties ---------------------------------

prop_sumAbs t x =
    V.sumAbs x ~== (sum $ map norm1 $ V.elems x)
  where
    _ = typed t x

prop_norm2 t x =
    V.norm2 x ~== (sqrt $ sum $ map (^^2) $ map norm $ V.elems x)
  where
    _ = typed t x

prop_whichMaxAbs1 t x =
    (V.dim x > 0) && all (not . isNaN) (map norm1 $ V.elems x) ==>
        V.at x i === e
  where
    (i,e) = V.whichMaxAbs x      
    _     = typed t x

prop_whichMaxAbs2 t x =
    (V.dim x > 0) && all (not . isNaN) (map norm1 $ V.elems x) ==>
        all (<= norm1 e) $ map norm1 (V.elems x)
  where
    (_,e) = V.whichMaxAbs x
    _     = typed t x

prop_dot t (VectorPair x y) =
    V.dot x y ~== sum (V.elems (x * conj y))
  where
    conj = V.conj
    (*)  = V.mul
    _    = typed t x

prop_kronecker t x y =
    x `V.kronecker` y ===
        V.fromList (V.dim x * V.dim y)
                   [ e*f | e <- V.elems x, f <- V.elems y ]
  where
    _ = typed t x
