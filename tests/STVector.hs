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

import Numeric.LinearAlgebra.Vector( Vector )
import qualified Numeric.LinearAlgebra.Vector as V
import Numeric.LinearAlgebra.Vector.ST( STVector )
import qualified Numeric.LinearAlgebra.Vector.ST as V

import Test.QuickCheck.LinearAlgebra( TestElem(..), Dim(..), Index(..),
    Assocs(..), VectorPair(..), VectorTriple(..) )
import qualified Test.QuickCheck.LinearAlgebra as Test

import Typed


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
    , testPropertyDZ "shiftTo" prop_shiftTo prop_shiftTo
    , testPropertyDZ "addTo" prop_addTo prop_addTo
    , testPropertyDZ "addToWithScales" prop_addToWithScales prop_addToWithScales
    , testPropertyDZ "subTo" prop_subTo prop_subTo
    , testPropertyDZ "scaleTo" prop_scaleTo prop_scaleTo
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
        V.copyTo mx my
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
        x `mutatesToVector` (x `V.replace` [(i,e)]) $ \mx -> do
            V.write mx i e
  where
    _ = e == t

prop_modify t (Index n i) (Blind f) =
    forAll (Test.vector n) $ \x -> runST $
        x `mutatesToVector`
            (typed t $ x `V.replace` [(i, f $ V.at x i)]) $ \mx ->
                V.modify mx i f

prop_getElems t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (V.map f x) $ \mx -> do
            es <- V.getElems mx
            V.mapTo f mx mx
            return $ es === V.elems (V.map f x)
  where
    _ = typed t x

prop_getElems' t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (V.map f x) $ \mx -> do
            es <- V.getElems' mx
            V.mapTo f mx mx
            return $ es === V.elems x
  where
    _ = typed t x

prop_getAssocs t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (V.map f x) $ \mx -> do
            ies <- V.getAssocs mx
            V.mapTo f mx mx
            return $ ies === V.assocs (V.map f x)
  where
    _ = typed t x

prop_getAssocs' t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (V.map f x) $ \mx -> do
            ies <- V.getAssocs' mx
            V.mapTo f mx mx
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
        x `mutatesToVector` (typed t $ V.replace x ies) $ \mx ->
            V.setAssocs mx ies

prop_mapTo t (Blind f) = binaryProp t
    (\x -> V.map f x)
    (\mx my -> V.mapTo f mx my)

prop_zipWithTo t (Blind f) = ternaryProp t
    (\x y -> V.zipWith f x y)
    (\mx my mz -> V.zipWithTo f mx my mz)

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
            V.kroneckerTo mx my mz

prop_shiftTo t e = binaryProp t
    (\x -> V.shift e x)
    (\mx my -> V.shiftTo e mx my)
    
prop_addTo t = ternaryProp t V.add V.addTo

prop_addToWithScales t e f = ternaryProp t
    (\x y -> V.addWithScales e x f y)
    (\mx my mz -> V.addToWithScales e mx f my mz)
    
prop_subTo t = ternaryProp t V.sub V.subTo

prop_scaleTo t e = binaryProp t
    (\x -> V.scale e x)
    (\mx my -> V.scaleTo e mx my)

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
    prop1 t x = runST $
        x `mutatesToVector` (imm_f x) $ \mx -> f mx mx
      where _ = typed t x
    prop2 t (VectorPair x y) = runST $ let _ = typed t x in
        x `readOnlyVector` \mx ->
        y `mutatesToVector` (imm_f x) $ \my -> f mx my
      where _ = typed t x
    in (    label "1" prop1
        .&. label "2" prop2
       )
    
ternaryProp :: (AEq e, Arbitrary e, Show e, Storable e, Testable a)
            => e
            -> (Vector e -> Vector e -> Vector e)
            -> (forall s . STVector s e -> STVector s e -> STVector s e -> ST s a)
            -> Property
ternaryProp t imm_f f = let
    prop1 t x = runST $
        x `mutatesToVector` (imm_f x x) $ \mx -> f mx mx mx
      where _ = typed t x
    prop2a t (VectorPair x y) = runST $
        x `readOnlyVector` \mx ->
        y `mutatesToVector` (imm_f x y) $ \my -> f mx my my
      where _ = typed t x
    prop2b t (VectorPair x y) = runST $
        x `mutatesToVector` (imm_f x y) $ \mx ->
        y `readOnlyVector` \my -> f mx my mx
      where _ = typed t x        
    prop3 t (VectorTriple x y z) = runST $
        x `readOnlyVector` \mx ->
        y `readOnlyVector` \my ->
        z `mutatesToVector` (imm_f x y) $ \mz -> f mx my mz
      where _ = typed t x
    in (    label "1" prop1
        .&. label "2a" prop2a
        .&. label "2b" prop2b
        .&. label "3" prop3
       )


readOnlyVector :: (Storable e, AEq e, Testable prop, Show e)
               => Vector e
               -> (STVector s e ->  ST s prop)
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
