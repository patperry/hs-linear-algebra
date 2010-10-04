{-# LANGUAGE NoMonomorphismRestriction, Rank2Types #-}
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

import Numeric.LinearAlgebra
import qualified Numeric.LinearAlgebra.Matrix as M
import qualified Numeric.LinearAlgebra.Matrix.ST as M

import Test.QuickCheck.LinearAlgebra( TestElem(..), Dim2(..), Index2(..),
    Assocs2(..), MatrixPair(..), MatrixTriple(..) )
import qualified Test.QuickCheck.LinearAlgebra as Test

import STVector( readOnlyVector )
import Typed


tests_STMatrix = testGroup "STMatrix"
    [ testPropertyI "new_" prop_new_
    , testPropertyI "new" prop_new
    , testPropertyI "newCopy" prop_newCopy
    , testPropertyI "copyTo" prop_copyTo
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
    , testPropertyDZ "shiftDiagByM_" prop_shiftDiagByM_ prop_shiftDiagByM_
    , testPropertyDZ "shiftDiagByWithScaleToM_"
        prop_shiftDiagByWithScaleM_ prop_shiftDiagByWithScaleM_
    , testPropertyDZ "addTo" prop_addTo prop_addTo
    , testPropertyDZ "subTo" prop_subTo prop_subTo
    , testPropertyDZ "scaleByM_" prop_scaleByM_ prop_scaleByM_
    , testPropertyDZ "addWithScaleM_" prop_addWithScaleM_ prop_addWithScaleM_
    , testPropertyDZ "scaleRowsTo (1)" prop_scaleRowsTo1 prop_scaleRowsTo1
    , testPropertyDZ "scaleRowsTo (2)" prop_scaleRowsTo2 prop_scaleRowsTo2
    , testPropertyDZ "scaleCols_" prop_scaleCols_ prop_scaleCols_
    , testPropertyDZ "negateTo" prop_negateTo prop_negateTo
    , testPropertyDZ "conjugateTo" prop_conjugateTo prop_conjugateTo
    ]
    
    
prop_new_ t (Dim2 n) = 
    (M.dim $ typed t $ M.create $ M.new_ n) === n

prop_new t (Dim2 n) e = 
    (M.create $ M.new n e) === (typed t $ M.constant n e)
    
prop_newCopy t x = 
    (M.create $ M.newCopy x) === x
  where
    _ = typed t x
        
prop_copyTo t (MatrixPair x y) = runST $
    x `readOnlyMatrix` \mx ->
    y `mutatesToMatrix` x $ \my ->
        M.copyTo my mx
  where
    _ = typed t x

prop_read t (Index2 n i) =
    forAll (typed t `fmap` Test.matrix n) $ \x -> runST $
        x `readOnlyMatrix` \mx -> do
            e <- M.read mx i
            return $ e === M.at x i

prop_write t (Index2 n i) e =
    forAll (Test.matrix n) $ \x -> runST $
        x `mutatesToMatrix` (x `M.replace` [(i,e)]) $ \mx -> do
            M.write mx i e
  where
    _ = e == t

prop_modify t (Index2 n i) (Blind f) =
    forAll (Test.matrix n) $ \x -> runST $
        x `mutatesToMatrix`
            (typed t $ x `M.replace` [(i, f $ M.at x i)]) $ \mx ->
                M.modify mx i f

prop_getElems t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToMatrix` (M.map f x) $ \mx -> do
            es <- M.getElems mx
            M.mapTo mx f mx
            return $ es === M.elems (M.map f x)
  where
    _ = typed t x

prop_getElems' t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToMatrix` (M.map f x) $ \mx -> do
            es <- M.getElems' mx
            M.mapTo mx f mx
            return $ es === M.elems x
  where
    _ = typed t x

prop_getAssocs t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToMatrix` (M.map f x) $ \mx -> do
            ies <- M.getAssocs mx
            M.mapTo mx f mx
            return $ ies === M.assocs (M.map f x)
  where
    _ = typed t x

prop_getAssocs' t x =
    forAll arbitrary $ \(Blind f) -> runST $
        x `mutatesToMatrix` (M.map f x) $ \mx -> do
            ies <- M.getAssocs' mx
            M.mapTo mx f mx
            return $ ies === M.assocs x
  where
    _ = typed t x

prop_setElems t x =
    forAll (QC.vector $ m*n) $ \es -> runST $
        x `mutatesToMatrix` (M.fromList (m,n) es) $ \mx ->
            M.setElems mx es
  where
    (m,n) = M.dim x
    _ = typed t x

prop_setAssocs t (Assocs2 mn ies) =
    forAll (Test.matrix mn) $ \x -> runST $
        x `mutatesToMatrix` (typed t $ M.replace x ies) $ \mx ->
            M.setAssocs mx ies

prop_mapTo t (Blind f) = binaryProp t
    (\x -> M.map f x)
    (\dst mx -> M.mapTo dst f mx)

prop_zipWithTo t (Blind f) = ternaryProp t
    (\x y -> M.zipWith f x y)
    (\dst mx my -> M.zipWithTo dst f mx my)

prop_shiftDiagByM_ t a =
    forAll (Test.vector (min m n)) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `mutatesToMatrix` (M.shiftDiagBy s a) $ \ma ->
            M.shiftDiagByM_ ms ma
  where
    (m,n) = M.dim a
    _ = typed t a

prop_shiftDiagByWithScaleM_ t e a =
    forAll (Test.vector (min m n)) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `mutatesToMatrix` (M.shiftDiagByWithScale e s a) $ \ma ->
            M.shiftDiagByWithScaleM_ e ms ma
  where
    (m,n) = M.dim a
    _ = typed t a

prop_addTo t = ternaryProp t M.add M.addTo

prop_subTo t = ternaryProp t M.sub M.subTo

prop_scaleByM_ t e a = runST $
    a `mutatesToMatrix` (M.scaleBy e a) $ \ma ->
        M.scaleByM_ e ma
  where
    _ = typed t a

prop_addWithScaleM_ t alpha (MatrixPair x y) = runST $
    x `readOnlyMatrix` \mx ->
    y `mutatesToMatrix` (M.addWithScale alpha x y) $ \my ->
        M.addWithScaleM_ alpha mx my
  where
    _ = typed t x

prop_scaleRowsTo1 t a =
    forAll (Test.vector m) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `mutatesToMatrix` (M.scaleRows s a) $ \ma ->
            M.scaleRowsTo ma ms ma
  where
    (m,n) = M.dim a
    _ = typed t a

prop_scaleRowsTo2 t (MatrixPair a b) =
    forAll (Test.vector m) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `readOnlyMatrix` \ma ->
        b `mutatesToMatrix` (M.scaleRows s a) $ \mb ->
            M.scaleRowsTo mb ms ma
  where
    (m,n) = M.dim a
    _ = typed t a

prop_scaleCols_ t a =
    forAll (Test.vector n) $ \s -> runST $
        s `readOnlyVector` \ms ->
        a `mutatesToMatrix` (M.scaleCols s a) $ \ma ->
            M.scaleCols_ ma ms
  where
    (m,n) = M.dim a
    _ = typed t a

prop_negateTo t = binaryProp t M.negate M.negateTo

prop_conjugateTo t = binaryProp t M.conjugate M.conjugateTo


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
        y `mutatesToMatrix` (imm_f x) $ \my -> f my mx
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
        y `mutatesToMatrix` (imm_f x y) $ \my -> f my mx my
      where _ = typed t x
    prop2b t (MatrixPair x y) = runST $
        x `mutatesToMatrix` (imm_f x y) $ \mx ->
        y `readOnlyMatrix` \my -> f mx mx my
      where _ = typed t x        
    prop3 t (MatrixTriple x y z) = runST $
        x `readOnlyMatrix` \mx ->
        y `readOnlyMatrix` \my ->
        z `mutatesToMatrix` (imm_f x y) $ \mz -> f mz mx my
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
    mx <- M.newCopy x
    prop <- f mx
    x' <- M.freeze mx
    return $ prop .&. (x' === x_new)
