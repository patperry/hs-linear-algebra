module STVector (
    tests_STVector
    ) where

import Control.Monad
import Control.Monad.ST
import Data.AEq
import Debug.Trace
import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import BLAS.Elem
import BLAS.Vector
import BLAS.Vector.ST

import Test.QuickCheck.BLAS( TestElem(..), Dim(..), Index(..), Assocs(..),
    VectorPair(..), VectorTriple(..) )
import qualified Test.QuickCheck.BLAS as Test

import Typed


tests_STVector = testGroup "STVector"
    [ testPropertyI "new_" prop_new_
    , testPropertyI "new" prop_new
    , testPropertyI "newCopy" prop_newCopy
    , testPropertyI "copyTo" prop_copyTo
    , testPropertyDZ "swap" prop_swap prop_swap
    , testPropertyI "read" prop_read
    , testPropertyI "write" prop_write
    , testPropertyI "getElems" prop_getElems
    , testPropertyI "getElems'" prop_getElems'
    , testPropertyI "getAssocs" prop_getAssocs
    , testPropertyI "getAssocs'" prop_getAssocs'
    , testPropertyI "setElems" prop_setElems
    , testPropertyI "setAssocs" prop_setAssocs
    , testPropertyI "mapTo 1" prop_mapTo1
    , testPropertyI "mapTo 2" prop_mapTo2
    , testPropertyI "zipWithTo 1" prop_zipWithTo1
    , testPropertyI "zipWithTo 2a" prop_zipWithTo2a
    , testPropertyI "zipWithTo 2b" prop_zipWithTo2b
    , testPropertyI "zipWithTo 3" prop_zipWithTo3
    , testPropertyDZ "getSum" prop_getSum prop_getSum
    , testPropertyDZ "getSumAbs" prop_getSumAbs prop_getSumAbs
    , testPropertyDZ "getNorm2" prop_getNorm2 prop_getNorm2
    , testPropertyDZ "getDot" prop_getDot prop_getDot
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

prop_mapTo1 t (Blind f) x = runST $
    x `mutatesToVector` (mapVector f x) $ \mx ->
        mapToVector f mx mx
  where
    _ = typed t x

prop_mapTo2 t (Blind f) (VectorPair x y) = runST $
    x `readOnlyVector` \mx ->
    y `mutatesToVector` (mapVector f x) $ \my ->
        mapToVector f mx my
  where
    _ = typed t x
    _ = typed t y


prop_zipWithTo1 t (Blind f) x = runST $
    x `mutatesToVector` (zipWithVector f x x) $ \mx ->
        zipWithToVector f mx mx mx
  where
    _ = typed t x

prop_zipWithTo2a t (Blind f) (VectorPair x y) = runST $
    x `mutatesToVector` (zipWithVector f x y) $ \mx ->
    y `readOnlyVector` \my ->
        zipWithToVector f mx my mx
  where
    _ = typed t x
    _ = typed t y    

prop_zipWithTo2b t (Blind f) (VectorPair x y) = runST $
    x `readOnlyVector` \mx ->
    y `mutatesToVector` (zipWithVector f x y) $ \my ->
        zipWithToVector f mx my my
  where
    _ = typed t x
    _ = typed t y        

prop_zipWithTo3 t (Blind f) (VectorTriple x y z) = runST $
    x `readOnlyVector` \mx ->
    y `readOnlyVector` \my ->
    z `mutatesToVector` (zipWithVector f x y) $ \mz ->
        zipWithToVector f mx my mz
  where
    _ = typed t x
    _ = typed t y
    _ = typed t z        

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

readOnlyVector :: (Storable e, AEq e, Testable prop)
               => Vector e
               -> (STVector s e ->  ST s prop)
               -> ST s Property
readOnlyVector x = mutatesToVector x x

mutatesToVector :: (Storable e, AEq e, Testable prop)
                => Vector e
                -> Vector e
                -> (STVector s e -> ST s prop)
                -> ST s Property
mutatesToVector x x_new f = do
    mx <- thawVector x
    prop <- f mx
    x' <- freezeVector mx
    return $ prop .&. (x' === x_new)
