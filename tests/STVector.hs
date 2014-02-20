{-# LANGUAGE NoMonomorphismRestriction, Rank2Types #-}
module STVector (
    tests_STVector,
    mutatesToVector,
    readOnlyVector,
    ) where

import Control.Monad
import Control.Monad.ST
import Data.AEq
import Data.Complex( magnitude )
import Debug.Trace
import Foreign( Storable )
import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Numeric.LinearAlgebra.Vector( Vector, STVector )
import qualified Numeric.LinearAlgebra.Vector as V

import Test.QuickCheck.LinearAlgebra( TestElem(..), Dim(..), Index(..),
    Assocs(..), VectorPair(..), VectorTriple(..) )
import qualified Test.QuickCheck.LinearAlgebra as Test
import Test.QuickCheck.Property(rejected)

import Typed


--- this could be a bad idea ..
instance Testable () where
  property _ = property rejected


tests_STVector = testGroup "STVector"
    [ testPropertyI "new_" prop_new_
    , testPropertyI "new" prop_new
    , testPropertyI "newCopy" prop_newCopy
    , testPropertyI "copyTo" prop_copyTo
    , testPropertyDZ "swap" prop_swap prop_swap
    , testPropertyI "read" prop_read
    , testPropertyI "write" prop_write
    , testPropertyI "modify" prop_modify
    , testPropertyI "getElems" prop_getElems
    , testPropertyI "getElems'" prop_getElems'
    , testPropertyI "getAssocs" prop_getAssocs
    , testPropertyI "getAssocs'" prop_getAssocs'
    , testPropertyI "setElems" prop_setElems
    , testPropertyI "setAssocs" prop_setAssocs
    , testPropertyI "mapTo" prop_mapTo
    , testPropertyI "zipWithTo" prop_zipWithTo
    , testPropertyDZ "getSumAbs" prop_getSumAbs prop_getSumAbs
    , testPropertyDZ "getNorm2" prop_getNorm2 prop_getNorm2
    , testPropertyDZ "getDot" prop_getDot prop_getDot
    , testPropertyDZ "kroneckerTo" prop_kroneckerTo prop_kroneckerTo
    , testPropertyDZ "addTo" prop_addTo prop_addTo
    , testPropertyDZ "subTo" prop_subTo prop_subTo
    , testPropertyDZ "mulTo" prop_mulTo prop_mulTo
    , testPropertyDZ "negateTo" prop_negateTo prop_negateTo
    , testPropertyDZ "conjugateTo" prop_conjugateTo prop_conjugateTo
    , testPropertyDZ "absTo" prop_absTo prop_absTo
    , testPropertyDZ "signumTo" prop_signumTo prop_signumTo
    , testPropertyDZ "divTo" prop_divTo prop_divTo
    , testPropertyDZ "recipTo" prop_recipTo prop_recipTo
    , testPropertyDZ "sqrtTo" prop_sqrtTo prop_sqrtTo
    , testPropertyDZ "expTo" prop_expTo prop_expTo
    , testPropertyDZ "logTo" prop_logTo prop_logTo
    , testPropertyDZ "powTo" prop_powTo prop_powTo
    , testPropertyDZ "sinTo" prop_sinTo prop_sinTo
    , testPropertyDZ "cosTo" prop_cosTo prop_cosTo
    , testPropertyDZ "tanTo" prop_tanTo prop_tanTo
    , testPropertyDZ "asinTo" prop_asinTo prop_asinTo
    , testPropertyDZ "acosTo" prop_acosTo prop_acosTo
    , testPropertyDZ "atanTo" prop_atanTo prop_atanTo
    , testPropertyDZ "sinhTo" prop_sinhTo prop_sinhTo
    , testPropertyDZ "coshTo" prop_coshTo prop_coshTo
    , testPropertyDZ "tanhTo" prop_tanhTo prop_tanhTo
    , testPropertyDZ "asinhTo" prop_asinhTo prop_asinhTo
    , testPropertyDZ "acoshTo" prop_acoshTo prop_acoshTo
    , testPropertyDZ "atanhTo" prop_atanhTo prop_atanhTo
    ]
    

{-
  getSumAbs:

    Double: [Failed]
*** Failed! Falsifiable (after 34 tests and 1 shrink): 
fromList [48.00809219086489,-28.969358622663314,63.04503159639685,5.714612973491514,-5.312024428882535,60.56484370204321,20.339282130713237,-25.609691263383986,11.776052643355653,-72.22738635520177]
LHS
(used seed -5655945031033394469)

    Complex Double: [Failed]
*** Failed! Falsifiable (after 46 tests): 
fromList [(-126.63684010249744) :+ (-25.732604105949363),48.30271669448071 :+ (-39.257907001851365),(-742.581416038894) :+ (-40.94772884905545),(-21.033310146590043) :+ 11.040319311574454,(-32.47272728995846) :+ (-75.38889146104881),34.61215019691849 :+ 63.02796474654257,25.495435934718248 :+ 14.86496960897311,285.1569736733574 :+ 20.592531068434027,(-75.03887610996705) :+ 62.568157533798285,13.719037128011735 :+ 49.186886152517744,(-23.33031638959207) :+ (-41.802953657765805),46.325514491999364 :+ (-40.3054052139413),(-97.8445174356864) :+ 62.27686352261467,(-210.08994152230045) :+ 23.4837419831381,(-47.12015960420664) :+ 9.502863461401958]
LHS

    Double: [Failed]
*** Failed! Falsifiable (after 61 tests): 
VectorPair (fromList [24.074823693680273,166.2322847046224,-2.8103216076602617,6.544349669641704,157.92295537758915,-52.792112116283036,-56.334929255514766,12.079209274206134,-51.61912866052746,163.39401573569245,-49.34434955256682,-17.020978906318945,116.41247221011058]) (fromList [91.62365332716433,-46.23902424675873,151.71325376854423,34.663441501788725,-2.5921683575991397,-48.720845208466315,51.34009695413845,41.16706632522662,-52.1781413047511,18.0562293933656,-198.82058422020668,-4.213635327473868,40.59717128981193])
LHS
LHS
(used seed -6941761855370166739)

    Complex Double: [OK, passed 100 tests]

-}    
    
prop_new_ t (Dim n) = 
    (V.dim $ typed t $ V.create $ V.new_ n) === n

prop_new t (Dim n) e = 
    (V.create $ V.new n e) === (typed t $ V.constant n e)
    
prop_newCopy t x = 
    (V.create $ V.newCopy x) === x
  where
    _ = typed t x
        
prop_copyTo t (VectorPair x y) = runST $
    x `readOnlyVector` \mx ->
    y `mutatesToVector` x $ \my ->
        V.copyTo my mx
  where
    _ = typed t x

prop_swap t (VectorPair x y) = runST $
    x `mutatesToVector` y $ \mx ->
    y `mutatesToVector` x $ \my ->
        V.swap mx my
  where
    _ = typed t x
    
prop_read t (Index n i) =
    forAll (typed t `fmap` Test.vector n) $ \x -> runST $
        x `readOnlyVector` \mx -> do
            e <- V.read mx i
            return $ e === V.at x i

prop_write t (Index n i) e =
    forAll (Test.vector n) $ \x -> runST $
        x `mutatesToVector` (x `V.update` [(i,e)]) $ \mx -> do
            V.write mx i e
  where
    _ = e == t

prop_modify t (Index n i) (Blind f) =
    forAll (Test.vector n) $ \x -> runST $
        x `mutatesToVector`
            (typed t $ x `V.update` [(i, f $ V.at x i)]) $ \mx ->
                V.modify mx i f

prop_getElems t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (V.map f x) $ \mx -> do
            es <- V.getElems mx
            V.mapTo mx f mx
            return $ es === V.elems (V.map f x)
  where
    _ = typed t x

prop_getElems' t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (V.map f x) $ \mx -> do
            es <- V.getElems' mx
            V.mapTo mx f mx
            return $ es === V.elems x
  where
    _ = typed t x

prop_getAssocs t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (V.map f x) $ \mx -> do
            ies <- V.getAssocs mx
            V.mapTo mx f mx
            return $ ies === V.assocs (V.map f x)
  where
    _ = typed t x

prop_getAssocs' t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (V.map f x) $ \mx -> do
            ies <- V.getAssocs' mx
            V.mapTo mx f mx
            return $ ies === V.assocs x
  where
    _ = typed t x

prop_setElems t x =
    forAll (QC.vector n) $ \es -> runST $
        x `mutatesToVector` (V.fromList n es) $ \mx ->
            V.setElems mx es
  where
    n = V.dim x
    _ = typed t x

prop_setAssocs t (Assocs n ies) =
    forAll (Test.vector n) $ \x -> runST $
        x `mutatesToVector` (typed t $ V.update x ies) $ \mx ->
            V.setAssocs mx ies

prop_mapTo t (Blind f) = binaryProp t
    (\x -> V.map f x)
    (\dst mx -> V.mapTo dst f mx)

prop_zipWithTo t (Blind f) = ternaryProp t
    (\x y -> V.zipWith f x y)
    (\dst mx my -> V.zipWithTo dst f mx my)

prop_getSumAbs t x = runST $
    x `readOnlyVector` \mx -> do
        s <- V.getSumAbs mx
        return $ s === V.sumAbs x
  where
    _ = typed t x

prop_getNorm2 t x = runST $
    x `readOnlyVector` \mx -> do
        n <- V.getNorm2 mx
        return $ n === V.norm2 x
  where
    _ = typed t x

prop_getDot t (VectorPair x y) = runST $
    x `readOnlyVector` \mx ->
    y `readOnlyVector` \my -> do
        d <- V.getDot mx my
        return $ d === V.dot x y
  where
    _ = typed t x

prop_kroneckerTo t = sized $ \s -> resize (s `div` 2) $
    forAll arbitrary $ \x ->
    forAll arbitrary $ \y ->
    forAll (Test.vector $ V.dim x * V.dim y) $ \z -> runST $
        x `readOnlyVector` \mx ->
        y `readOnlyVector` \my ->
        z `mutatesToVector` (typed t $ V.kronecker x y) $ \mz ->
            V.kroneckerTo mz mx my

prop_addTo t = ternaryProp t V.add V.addTo
prop_subTo t = ternaryProp t V.sub V.subTo
prop_mulTo t = ternaryProp t V.mul V.mulTo
prop_negateTo t = binaryProp t V.negate V.negateTo
prop_conjugateTo t = binaryProp t V.conjugate V.conjugateTo
prop_absTo t = binaryProp t V.abs V.absTo
prop_signumTo t = binaryProp t V.signum V.signumTo

prop_divTo t = ternaryProp t V.div V.divTo
prop_recipTo t = binaryProp t V.recip V.recipTo

prop_sqrtTo t = binaryProp t V.sqrt V.sqrtTo
prop_expTo t = binaryProp t V.exp V.expTo
prop_logTo t = binaryProp t V.log V.logTo
prop_powTo t = ternaryProp t V.pow V.powTo
prop_sinTo t = binaryProp t V.sin V.sinTo
prop_cosTo t = binaryProp t V.cos V.cosTo
prop_tanTo t = binaryProp t V.tan V.tanTo
prop_asinTo t = binaryProp t V.asin V.asinTo
prop_acosTo t = binaryProp t V.acos V.acosTo
prop_atanTo t = binaryProp t V.atan V.atanTo
prop_sinhTo t = binaryProp t V.sinh V.sinhTo
prop_coshTo t = binaryProp t V.cosh V.coshTo
prop_tanhTo t = binaryProp t V.tanh V.tanhTo
prop_asinhTo t = binaryProp t V.asinh V.asinhTo
prop_acoshTo t = binaryProp t V.acosh V.acoshTo
prop_atanhTo t = binaryProp t V.atanh V.atanhTo


binaryProp :: (AEq e, Arbitrary e, Show e, Storable e, Testable a)
            => e
            -> (Vector e -> Vector e)
            -> (forall s . STVector s e -> STVector s e -> ST s a)
            -> Property
binaryProp t imm_f f = let
    prop2 t (VectorPair x y) = runST $ let _ = typed t x in
        x `readOnlyVector` \mx ->
        y `mutatesToVector` (imm_f x) $ \my -> f my mx
      where _ = typed t x
    in label "" prop2
    
ternaryProp :: (AEq e, Arbitrary e, Show e, Storable e, Testable a)
            => e
            -> (Vector e -> Vector e -> Vector e)
            -> (forall s . STVector s e -> STVector s e -> STVector s e -> ST s a)
            -> Property
ternaryProp t imm_f f = let
    prop3 t (VectorTriple x y z) = runST $
        x `readOnlyVector` \mx ->
        y `readOnlyVector` \my ->
        z `mutatesToVector` (imm_f x y) $ \mz -> f mz mx my
      where _ = typed t x
    in label "" prop3
       


readOnlyVector :: (Storable e, AEq e, Testable prop, Show e)
               => Vector e
               -> (STVector s e -> ST s prop)
               -> ST s Property
readOnlyVector x = mutatesToVector x x

mutatesToVector :: (Storable e, AEq e, Testable prop, Show e)
          => Vector e
          -> Vector e
          -> (STVector s e -> ST s prop)
          -> ST s Property
mutatesToVector x x_new f = do
    mx <- V.newCopy x
    prop <- f mx
    x' <- V.freeze mx
    return $ prop .&. (x' === x_new)
