{-# LANGUAGE NoMonomorphismRestriction #-}
module STMatrix (
    tests_STMatrix,
    mutatesToMatrix,
    readOnlyMatrix,
    
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

import BLAS.Elem
import BLAS.Matrix
import BLAS.Matrix.ST
import BLAS.Vector
import BLAS.Vector.ST

import Test.QuickCheck.BLAS( TestElem(..), Dim2(..), Index2(..), Assocs2(..),
    MatrixPair(..), MatrixTriple(..) )
import qualified Test.QuickCheck.BLAS as Test

import STVector( readOnlyVector )
import Typed


tests_STMatrix = testGroup "STMatrix"
    [ testPropertyI "new_" prop_new_
    , testPropertyI "new" prop_new
    , testPropertyI "newCopy" prop_newCopy
    , testPropertyI "copyTo" prop_copyTo
    , testPropertyI "read" prop_read
    , testPropertyI "write" prop_write
    , testPropertyI "getElems" prop_getElems
    , testPropertyI "getElems'" prop_getElems'
    , testPropertyI "getAssocs" prop_getAssocs
    , testPropertyI "getAssocs'" prop_getAssocs'
    , testPropertyI "setElems" prop_setElems
    , testPropertyI "setAssocs" prop_setAssocs
    , testPropertyI "mapTo" prop_mapTo
    , testPropertyI "zipWithTo" prop_zipWithTo
    , testPropertyDZ "shiftTo" prop_shiftTo prop_shiftTo
    , testPropertyDZ "shiftDiagTo (1)" prop_shiftDiagTo1 prop_shiftDiagTo1
    , testPropertyDZ "shiftDiagTo (2)" prop_shiftDiagTo2 prop_shiftDiagTo2
    , testPropertyDZ "shiftDiagToWithScale (1)"
        prop_shiftDiagToWithScale1 prop_shiftDiagToWithScale1
    , testPropertyDZ "shiftDiagToWithScale (2)"
        prop_shiftDiagToWithScale2 prop_shiftDiagToWithScale2
    , testPropertyDZ "addTo" prop_addTo prop_addTo
    , testPropertyDZ "addToWithScale" prop_addToWithScale prop_addToWithScale
    , testPropertyDZ "subTo" prop_subTo prop_subTo
    , testPropertyDZ "scaleTo" prop_scaleTo prop_scaleTo
    , testPropertyDZ "scaleRowsTo (1)" prop_scaleRowsTo1 prop_scaleRowsTo1
    , testPropertyDZ "scaleRowsTo (2)" prop_scaleRowsTo2 prop_scaleRowsTo2
    , testPropertyDZ "scaleColsTo (1)" prop_scaleColsTo1 prop_scaleColsTo1
    , testPropertyDZ "scaleColsTo (2)" prop_scaleColsTo2 prop_scaleColsTo2
    , testPropertyDZ "negateTo" prop_negateTo prop_negateTo
    ]
    
    
prop_new_ t (Dim2 n) = 
    (dimMatrix $ typed t $ runMatrix $ newMatrix_ n) === n

prop_new t (Dim2 n) e = 
    (runMatrix $ newMatrix n e) === (typed t $ constantMatrix n e)
    
prop_newCopy t x = 
    (runMatrix $ newCopyMatrix x) === x
  where
    _ = typed t x
        
prop_copyTo t (MatrixPair x y) = runST $
    x `readOnlyMatrix` \mx ->
    y `mutatesToMatrix` x $ \my ->
        copyToMatrix mx my
  where
    _ = typed t x

prop_read t (Index2 n i) =
    forAll (typed t `fmap` Test.matrix n) $ \x -> runST $
        x `readOnlyMatrix` \mx -> do
            e <- readMatrix mx i
            return $ e === atMatrix x i

prop_write t (Index2 n i) e =
    forAll (Test.matrix n) $ \x -> runST $
        x `mutatesToMatrix` (x `replaceMatrix` [(i,e)]) $ \mx -> do
            writeMatrix mx i e
  where
    _ = e == t

prop_getElems t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToMatrix` (mapMatrix f x) $ \mx -> do
            es <- getElemsMatrix mx
            mapToMatrix f mx mx
            return $ es === elemsMatrix (mapMatrix f x)
  where
    _ = typed t x

prop_getElems' t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToMatrix` (mapMatrix f x) $ \mx -> do
            es <- getElemsMatrix' mx
            mapToMatrix f mx mx
            return $ es === elemsMatrix x
  where
    _ = typed t x

prop_getAssocs t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToMatrix` (mapMatrix f x) $ \mx -> do
            ies <- getAssocsMatrix mx
            mapToMatrix f mx mx
            return $ ies === assocsMatrix (mapMatrix f x)
  where
    _ = typed t x

prop_getAssocs' t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToMatrix` (mapMatrix f x) $ \mx -> do
            ies <- getAssocsMatrix' mx
            mapToMatrix f mx mx
            return $ ies === assocsMatrix x
  where
    _ = typed t x

prop_setElems t x =
    forAll (QC.vector $ m*n) $ \es -> runST $
        x `mutatesToMatrix` (listMatrix (m,n) es) $ \mx ->
            setElemsMatrix mx es
  where
    (m,n) = dimMatrix x
    _ = typed t x

prop_setAssocs t (Assocs2 mn ies) =
    forAll (Test.matrix mn) $ \x -> runST $
        x `mutatesToMatrix` (typed t $ replaceMatrix x ies) $ \mx ->
            setAssocsMatrix mx ies

prop_mapTo t (Blind f) = binaryProp t
    (\x -> mapMatrix f x)
    (\mx my -> mapToMatrix f mx my)

prop_zipWithTo t (Blind f) = ternaryProp t
    (\x y -> zipWithMatrix f x y)
    (\mx my mz -> zipWithToMatrix f mx my mz)

prop_shiftTo t e = binaryProp t
    (\x -> shiftMatrix e x)
    (\mx my -> shiftToMatrix e mx my)

prop_shiftDiagTo1 t a =
    forAll (Test.vector (min m n)) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `mutatesToMatrix` (shiftDiagMatrix s a) $ \ma ->
            shiftDiagToMatrix ms ma ma
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_shiftDiagTo2 t (MatrixPair a b) =
    forAll (Test.vector (min m n)) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `readOnlyMatrix` \ma ->
        b `mutatesToMatrix` (shiftDiagMatrix s a) $ \mb ->
            shiftDiagToMatrix ms ma mb
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_shiftDiagToWithScale1 t e a =
    forAll (Test.vector (min m n)) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `mutatesToMatrix` (shiftDiagMatrixWithScale e s a) $ \ma ->
            shiftDiagToMatrixWithScale e ms ma ma
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_shiftDiagToWithScale2 t e (MatrixPair a b) =
    forAll (Test.vector (min m n)) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `readOnlyMatrix` \ma ->
        b `mutatesToMatrix` (shiftDiagMatrixWithScale e s a) $ \mb ->
            shiftDiagToMatrixWithScale e ms ma mb
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_addTo t = ternaryProp t addMatrix addToMatrix

prop_addToWithScale t e f = ternaryProp t
    (\x y -> addMatrixWithScale e x f y)
    (\mx my mz -> addToMatrixWithScale e mx f my mz)
    
prop_subTo t = ternaryProp t subMatrix subToMatrix

prop_scaleTo t e = binaryProp t
    (\x -> scaleMatrix e x)
    (\mx my -> scaleToMatrix e mx my)

prop_scaleRowsTo1 t a =
    forAll (Test.vector m) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `mutatesToMatrix` (scaleRowsMatrix s a) $ \ma ->
            scaleRowsToMatrix ms ma ma
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_scaleRowsTo2 t (MatrixPair a b) =
    forAll (Test.vector m) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `readOnlyMatrix` \ma ->
        b `mutatesToMatrix` (scaleRowsMatrix s a) $ \mb ->
            scaleRowsToMatrix ms ma mb
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_scaleColsTo1 t a =
    forAll (Test.vector n) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `mutatesToMatrix` (scaleColsMatrix s a) $ \ma ->
            scaleColsToMatrix ms ma ma
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_scaleColsTo2 t (MatrixPair a b) =
    forAll (Test.vector n) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `readOnlyMatrix` \ma ->
        b `mutatesToMatrix` (scaleColsMatrix s a) $ \mb ->
            scaleColsToMatrix ms ma mb
  where
    (m,n) = dimMatrix a
    _ = typed t a

prop_negateTo t = binaryProp t negateMatrix negateToMatrix


binaryProp :: (AEq e, Arbitrary e, Show e, Storable e, Testable a)
            => e
            -> (Matrix e -> Matrix e)
            -> (forall s . STMatrix s e -> STMatrix s e -> ST s a)
            -> Property
binaryProp t imm_f f = let
    prop1 t x = runST $
        x `mutatesToMatrix` (imm_f x) $ \mx -> f mx mx
      where _ = typed t x
    prop2 t (MatrixPair x y) = runST $ let _ = typed t x in
        x `readOnlyMatrix` \mx ->
        y `mutatesToMatrix` (imm_f x) $ \my -> f mx my
      where _ = typed t x
    in (    label "1" prop1
        .&. label "2" prop2
       )
    
ternaryProp :: (AEq e, Arbitrary e, Show e, Storable e, Testable a)
            => e
            -> (Matrix e -> Matrix e -> Matrix e)
            -> (forall s . STMatrix s e -> STMatrix s e -> STMatrix s e -> ST s a)
            -> Property
ternaryProp t imm_f f = let
    prop1 t x = runST $
        x `mutatesToMatrix` (imm_f x x) $ \mx -> f mx mx mx
      where _ = typed t x
    prop2a t (MatrixPair x y) = runST $
        x `readOnlyMatrix` \mx ->
        y `mutatesToMatrix` (imm_f x y) $ \my -> f mx my my
      where _ = typed t x
    prop2b t (MatrixPair x y) = runST $
        x `mutatesToMatrix` (imm_f x y) $ \mx ->
        y `readOnlyMatrix` \my -> f mx my mx
      where _ = typed t x        
    prop3 t (MatrixTriple x y z) = runST $
        x `readOnlyMatrix` \mx ->
        y `readOnlyMatrix` \my ->
        z `mutatesToMatrix` (imm_f x y) $ \mz -> f mx my mz
      where _ = typed t x
    in (    label "1" prop1
        .&. label "2a" prop2a
        .&. label "2b" prop2b
        .&. label "3" prop3
       )


readOnlyMatrix :: (Storable e, AEq e, Testable prop, Show e)
               => Matrix e
               -> (STMatrix s e ->  ST s prop)
               -> ST s Property
readOnlyMatrix x = mutatesToMatrix x x

mutatesToMatrix :: (Storable e, AEq e, Testable prop, Show e)
                => Matrix e
                -> Matrix e
                -> (STMatrix s e -> ST s prop)
                -> ST s Property
mutatesToMatrix x x_new f = do
    mx <- thawMatrix x
    prop <- f mx
    x' <- freezeMatrix mx
    return $ prop .&. (x' === x_new)
