{-# LANGUAGE NoMonomorphismRestriction #-}
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
import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Numeric.LinearAlgebra

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
    , testPropertyI "update" prop_update
    , testPropertyI "getElems" prop_getElems
    , testPropertyI "getElems'" prop_getElems'
    , testPropertyI "getAssocs" prop_getAssocs
    , testPropertyI "getAssocs'" prop_getAssocs'
    , testPropertyI "setElems" prop_setElems
    , testPropertyI "setAssocs" prop_setAssocs
    , testPropertyI "mapTo" prop_mapTo
    , testPropertyI "zipWithTo" prop_zipWithTo
    , testPropertyDZ "getSum" prop_getSum prop_getSum
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
    , testPropertyDZ "conjTo" prop_conjTo prop_conjTo
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
    (dimVector $ typed t $ runVector $ newVector_ n) === n

prop_new t (Dim n) e = 
    (runVector $ newVector n e) === (typed t $ constantVector n e)
    
prop_newCopy t x = 
    (runVector $ newCopyVector x) === x
  where
    _ = typed t x
        
prop_copyTo t (VectorPair x y) = runST $
    x `readOnlyVector` \mx ->
    y `mutatesToVector` x $ \my ->
        copyToVector mx my
  where
    _ = typed t x

prop_swap t (VectorPair x y) = runST $
    x `mutatesToVector` y $ \mx ->
    y `mutatesToVector` x $ \my ->
        swapVector mx my
  where
    _ = typed t x
    
prop_read t (Index n i) =
    forAll (typed t `fmap` Test.vector n) $ \x -> runST $
        x `readOnlyVector` \mx -> do
            e <- readVector mx i
            return $ e === atVector x i

prop_write t (Index n i) e =
    forAll (Test.vector n) $ \x -> runST $
        x `mutatesToVector` (x `replaceVector` [(i,e)]) $ \mx -> do
            writeVector mx i e
  where
    _ = e == t

prop_update t (Index n i) (Blind f) =
    forAll (Test.vector n) $ \x -> runST $
        x `mutatesToVector`
            (typed t $ x `replaceVector` [(i, f $ atVector x i)]) $ \mx ->
                updateVector mx i f

prop_getElems t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (mapVector f x) $ \mx -> do
            es <- getElemsVector mx
            mapToVector f mx mx
            return $ es === elemsVector (mapVector f x)
  where
    _ = typed t x

prop_getElems' t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (mapVector f x) $ \mx -> do
            es <- getElemsVector' mx
            mapToVector f mx mx
            return $ es === elemsVector x
  where
    _ = typed t x

prop_getAssocs t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (mapVector f x) $ \mx -> do
            ies <- getAssocsVector mx
            mapToVector f mx mx
            return $ ies === assocsVector (mapVector f x)
  where
    _ = typed t x

prop_getAssocs' t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToVector` (mapVector f x) $ \mx -> do
            ies <- getAssocsVector' mx
            mapToVector f mx mx
            return $ ies === assocsVector x
  where
    _ = typed t x

prop_setElems t x =
    forAll (QC.vector n) $ \es -> runST $
        x `mutatesToVector` (listVector n es) $ \mx ->
            setElemsVector mx es
  where
    n = dimVector x
    _ = typed t x

prop_setAssocs t (Assocs n ies) =
    forAll (Test.vector n) $ \x -> runST $
        x `mutatesToVector` (typed t $ replaceVector x ies) $ \mx ->
            setAssocsVector mx ies

prop_mapTo t (Blind f) = binaryProp t
    (\x -> mapVector f x)
    (\mx my -> mapToVector f mx my)

prop_zipWithTo t (Blind f) = ternaryProp t
    (\x y -> zipWithVector f x y)
    (\mx my mz -> zipWithToVector f mx my mz)

prop_getSum t x = runST $
    x `readOnlyVector` \mx -> do
        s <- getSumVector mx
        return $ s === sumVector x
  where
    _ = typed t x

prop_getSumAbs t x = runST $
    x `readOnlyVector` \mx -> do
        s <- getSumAbsVector mx
        return $ s === sumAbsVector x
  where
    _ = typed t x

prop_getNorm2 t x = runST $
    x `readOnlyVector` \mx -> do
        n <- getNorm2Vector mx
        return $ n === norm2Vector x
  where
    _ = typed t x

prop_getDot t (VectorPair x y) = runST $
    x `readOnlyVector` \mx ->
    y `readOnlyVector` \my -> do
        d <- getDotVector mx my
        return $ d === dotVector x y
  where
    _ = typed t x

prop_kroneckerTo t = sized $ \s -> resize (s `div` 2) $
    forAll arbitrary $ \x ->
    forAll arbitrary $ \y ->
    forAll (Test.vector $ dimVector x * dimVector y) $ \z -> runST $
        x `readOnlyVector` \mx ->
        y `readOnlyVector` \my ->
        z `mutatesToVector` (typed t $ kroneckerVector x y) $ \mz ->
            kroneckerToVector mx my mz

prop_shiftTo t e = binaryProp t
    (\x -> shiftVector e x)
    (\mx my -> shiftToVector e mx my)
    
prop_addTo t = ternaryProp t addVector addToVector

prop_addToWithScales t e f = ternaryProp t
    (\x y -> addVectorWithScales e x f y)
    (\mx my mz -> addToVectorWithScales e mx f my mz)
    
prop_subTo t = ternaryProp t subVector subToVector

prop_scaleTo t e = binaryProp t
    (\x -> scaleVector e x)
    (\mx my -> scaleToVector e mx my)

prop_mulTo t = ternaryProp t mulVector mulToVector
prop_negateTo t = binaryProp t negateVector negateToVector
prop_conjTo t = binaryProp t conjVector conjToVector
prop_absTo t = binaryProp t absVector absToVector
prop_signumTo t = binaryProp t signumVector signumToVector

prop_divTo t = ternaryProp t divVector divToVector
prop_recipTo t = binaryProp t recipVector recipToVector

prop_sqrtTo t = binaryProp t sqrtVector sqrtToVector
prop_expTo t = binaryProp t expVector expToVector
prop_logTo t = binaryProp t logVector logToVector
prop_powTo t = ternaryProp t powVector powToVector
prop_sinTo t = binaryProp t sinVector sinToVector
prop_cosTo t = binaryProp t cosVector cosToVector
prop_tanTo t = binaryProp t tanVector tanToVector
prop_asinTo t = binaryProp t asinVector asinToVector
prop_acosTo t = binaryProp t acosVector acosToVector
prop_atanTo t = binaryProp t atanVector atanToVector
prop_sinhTo t = binaryProp t sinhVector sinhToVector
prop_coshTo t = binaryProp t coshVector coshToVector
prop_tanhTo t = binaryProp t tanhVector tanhToVector
prop_asinhTo t = binaryProp t asinhVector asinhToVector
prop_acoshTo t = binaryProp t acoshVector acoshToVector
prop_atanhTo t = binaryProp t atanhVector atanhToVector


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
    mx <- newCopyVector x
    prop <- f mx
    x' <- freezeVector mx
    return $ prop .&. (x' === x_new)
