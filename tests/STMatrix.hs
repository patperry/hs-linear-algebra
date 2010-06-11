{-# LANGUAGE Rank2Types, ScopedTypeVariables #-}
module STMatrix
    where

import Data.Elem.BLAS
import Data.Vector.Dense
import Data.Matrix.Dense
import Data.Matrix.Dense.ST

import qualified Test.Matrix.Dense as Test

import Driver



-------------------------- Creating Matrices --------------------------------

newMatrix_S = matrix

prop_NewMatrix (Assocs2 mn ijes) =
    newMatrix mn ijes `equivalent` newMatrix_S mn ijes

newListMatrix_S = listMatrix

prop_NewListMatrix (Nat2 mn) es =
    newListMatrix mn es `equivalent` newListMatrix_S mn es

---------------------- Reading and Writing Elements --------------------------

getSize_S a = ( size a, a )
prop_GetSize = getSize `implements` getSize_S

readElem_S a ij = ( a!ij, a )
prop_ReadElem (Index2 mn ij) =
    implementsFor mn (`readElem` ij) (`readElem_S` ij)

canModifyElem_S a ij = ( True, a )
prop_CanModifyElem ij = (`canModifyElem` ij) `implements` (`canModifyElem_S` ij)

writeElem_S a ij e = ( (), a // [(ij,e)] )
prop_WriteElem (Index2 mn ij) e =
    implementsFor mn (\a -> writeElem a ij e) (\a -> writeElem_S a ij e)

modifyElem_S a ij f = writeElem_S a ij $ f (a!ij)
prop_ModifyElem (Index2 mn ij) f =
    implementsFor mn (\a -> modifyElem a ij f) (\a -> modifyElem_S a ij f)

getIndices_S a = ( indices a, a )

prop_GetIndicesLazy   = getIndices `implements`  getIndices_S
prop_GetIndicesStrict = getIndices' `implements` getIndices_S

getElems_S a = ( elems a, a )
prop_GetElemsLazy   = getElems  `implements` getElems_S
prop_GetElemsStrict = getElems' `implements` getElems_S
    
getElemsLazyModifyWith_S f a = ( elems a', a' ) where a' = tmap f a
prop_GetElemsLazyModifyWith f =
    (\a -> do { es <- getElems a ; modifyWith f a ; return es })
    `implements `
    (getElemsLazyModifyWith_S f)

getElemsStrictModifyWith_S f a = ( elems a, a' ) where a' = tmap f a
prop_GetElemsStrictModifyWith f =
    (\a -> do { es <- getElems' a ; modifyWith f a ; return es })
    `implements `
    (getElemsStrictModifyWith_S f)

getAssocsLazyModifyWith_S f a = ( assocs a', a' ) where a' = tmap f a
prop_GetAssocsLazyModifyWith f =
    (\a -> do { ijes <- getAssocs a ; modifyWith f a ; return ijes })
    `implements` 
    getAssocsLazyModifyWith_S f

getAssocsStrictModifyWith_S f a = ( assocs a, a' ) where a' = tmap f a
prop_GetAssocsStrictModifyWith f =
    (\a -> do { ijes <- getAssocs' a ; modifyWith f a ; return ijes })
    `implements` 
    getAssocsStrictModifyWith_S f


----------------------------- Special Matrices --------------------------------

newZeroMatrix_S = zeroMatrix
prop_NewZeroMatrix (Nat2 mn) = 
    newZeroMatrix mn `equivalent` newZeroMatrix_S mn

setZeroMatrix_S a = ( (), newZeroMatrix_S (shape a) )
prop_SetZeroMatrix = setZeroMatrix `implements` setZeroMatrix_S

newConstantMatrix_S mn e = constantMatrix mn e
prop_NewConstantMatrix (Nat2 mn) e = 
    newConstantMatrix mn e `equivalent` newConstantMatrix_S mn e

setConstantMatrix_S e a = ( (), newConstantMatrix_S (shape a) e )
prop_SetConstantMatrix e = setConstantMatrix e `implements` setConstantMatrix_S e

newIdentityMatrix_S = identityMatrix
prop_NewIdentityMatrix (Nat2 mn) = 
    newIdentityMatrix mn `equivalent` newIdentityMatrix_S mn

setIdentityMatrix_S a = ( (), newIdentityMatrix_S (shape a) )
prop_SetIdentityMatrix =
    setIdentityMatrix `implements` setIdentityMatrix_S


---------------------------- Copying Matrices --------------------------------

newCopyMatrix_S a = ( a, a )
prop_NewCopyMatrix = 
    (\a -> newCopyMatrix a >>= abstract) `implements` newCopyMatrix_S

copyMatrix_S a b = ( (), b, b )
prop_CopyMatrix = copyMatrix `implements2` copyMatrix_S

swapMatrix_S a b = ( (), b, a )
prop_SwapMatrix = swapMatrix `implements2` swapMatrix_S


-------------------------- Unsary Matrix Operations --------------------------

doConj_S x = ( (), tmap conjugate x )
prop_DoConj = doConj `implements` doConj_S

scaleBy_S k x = ( (), tmap (k*) x )
prop_ScaleBy k = scaleBy k `implements` scaleBy_S k

shiftBy_S k x = ( (), tmap (k+) x )
prop_ShiftBy k = shiftBy k `implements` shiftBy_S k

modifyWith_S f x = ( (), tmap f x )
prop_ModifyWith f = modifyWith f `implements` modifyWith_S f

rank1UpdateMatrix_S a alpha x y = 
    let a' = a + matrixFromCol x <**> matrixFromRow (alpha *> (conj y))
    in ( (), a')
prop_Rank1UpdateMatrix alpha x y =
    implementsFor (dim x, dim y)
        (\a -> rank1UpdateMatrix   a alpha x y) 
        (\a -> rank1UpdateMatrix_S a alpha x y)     

getConjMatrix_S x = ( tmap conjugate x, x )
prop_GetConjMatrix = 
    (\x -> getConjMatrix x >>= abstract) `implements` getConjMatrix_S

getScaledMatrix_S k x = ( tmap (k*) x, x )
prop_GetScaledMatrix k = 
    (\x -> getScaledMatrix k x >>= abstract) `implements` (getScaledMatrix_S k)

getShiftedMatrix_S k x = ( tmap (k+) x, x )
prop_GetShiftedMatrix k = 
    (\x -> getShiftedMatrix k x >>= abstract) `implements` (getShiftedMatrix_S k)


------------------------- Binary Matrix Operations ---------------------------

addMatrix_S x y = ( (), x + y, y )
prop_AddMatrix = addMatrix `implements2` addMatrix_S

subMatrix_S x y = ( (), x - y, y )
prop_SubMatrix = subMatrix `implements2` subMatrix_S

axpyMatrix_S alpha x y = ( (), x, alpha *> x + y )
prop_AxpyMatrix alpha = axpyMatrix alpha `implements2` axpyMatrix_S alpha

mulMatrix_S x y = ( (), x * y, y )
prop_MulMatrix = mulMatrix `implements2` mulMatrix_S

divMatrix_S x y = ( (), x / y, y )
prop_DivMatrix = divMatrix `implements2` divMatrix_S

getAddMatrix_S x y = ( x + y, x, y )
prop_GetAddMatrix =
    (\x y -> getAddMatrix x y >>= abstract) `implements2` getAddMatrix_S

getSubMatrix_S x y = ( x - y, x, y )
prop_GetSubMatrix =
    (\x y -> getSubMatrix x y >>= abstract) `implements2` getSubMatrix_S

getMulMatrix_S x y = ( x * y, x, y )
prop_GetMulMatrix =
    (\x y -> getMulMatrix x y >>= abstract) `implements2` getMulMatrix_S

getDivMatrix_S x y = ( x / y, x, y )
prop_GetDivMatrix =
    (\x y -> getDivMatrix x y >>= abstract) `implements2` getDivMatrix_S

------------------------------------------------------------------------
-- 
-- The specification language
--
    
abstract :: (BLAS3 e) => STMatrix s e -> ST s (Matrix e)
abstract = freezeMatrix

commutes :: (AEq a, Show a, AEq e, BLAS3 e) =>
    STMatrix s e -> (STMatrix s e -> ST s a) ->
        (Matrix e -> (a,Matrix e)) -> ST s Bool
commutes x a f = do
    old <- abstract x
    r <- a x
    new <- abstract x
    let s      = f old
        s'     = (r,new)
        passed = s ~== s'
        
    when (not passed) $
        trace (printf ("expected `%s' but got `%s'") (show s) (show s'))
              return ()
              
    return passed

commutes2 :: (AEq a, Show a, AEq e, BLAS3 e) =>
    STMatrix s e -> STMatrix s e -> 
    (STMatrix s e ->  STMatrix s e -> ST s a) ->
        (Matrix e -> Matrix e -> (a,Matrix e,Matrix e)) -> ST s Bool
commutes2 x y a f = do
    oldX <- abstract x
    oldY <- abstract y
    r <- a x y
    newX <- abstract x
    newY <- abstract y
    let s      = f oldX oldY
        s'     = (r,newX,newY)
        passed = s ~== s'
        
    when (not passed) $
        trace (printf ("expected `%s' but got `%s'") (show s) (show s'))
            return ()
            
    return passed

equivalent :: (forall s . ST s (STMatrix s E)) -> Matrix E -> Bool
equivalent x s = runST $ do
    x' <- (x >>= abstract)
    when (not $ x' === s) $
        trace (printf ("expected `%s' but got `%s'") (show s) (show x'))
            return ()
    return (x' === s)
    
implements :: (AEq a, Show a) =>
    (forall s . STMatrix s E -> ST s a) ->
    (Matrix E -> (a,Matrix E)) -> 
        Property
a `implements` f =
    forAll arbitrary $ \(Nat2 mn) ->
        implementsFor mn a f

implements2 :: (AEq a, Show a) =>
    (forall s . STMatrix s E -> STMatrix s E -> ST s a) ->
    (Matrix E -> Matrix E -> (a,Matrix E,Matrix E)) -> 
        Property
a `implements2` f =
    forAll arbitrary $ \(Nat2 mn) ->
        implementsFor2 mn a f

implementsFor :: (AEq a, Show a) =>
    (Int,Int) ->
    (forall s . STMatrix s E -> ST s a) ->
    (Matrix E -> (a,Matrix E)) -> 
        Property
implementsFor mn a f =
    forAll (Test.matrix mn) $ \x ->
        runST $ do
            x' <- unsafeThawMatrix x
            commutes x' a f

implementsFor2 :: (AEq a, Show a) =>
    (Int,Int) ->
    (forall s . STMatrix s E -> STMatrix s E -> ST s a) ->
    (Matrix E -> Matrix E -> (a,Matrix E,Matrix E)) -> 
        Property
implementsFor2 mn a f =
    forAll (Test.matrix mn) $ \x ->
    forAll (Test.matrix mn) $ \y ->
        runST $ do
            x' <- unsafeThawMatrix x
            y' <- unsafeThawMatrix y
            commutes2 x' y' a f

implementsIf :: (AEq a, Show a) =>
    (forall s . STMatrix s E -> ST s Bool) ->
    (forall s . STMatrix s E -> ST s a) ->
    (Matrix E -> (a,Matrix E)) -> 
        Property
implementsIf pre a f =
    forAll arbitrary $ \(Nat2 mn) ->
    forAll (Test.matrix mn) $ \x ->
        runST ( do
            x' <- thawMatrix x
            pre x') ==>
        runST ( do
            x' <- unsafeThawMatrix x
            commutes x' a f )

implementsIf2 :: (AEq a, Show a) =>
    (forall s . STMatrix s E -> STMatrix s E -> ST s Bool) ->
    (forall s . STMatrix s E -> STMatrix s E -> ST s a) ->
    (Matrix E -> Matrix E -> (a,Matrix E,Matrix E)) -> 
        Property
implementsIf2 pre a f =
    forAll arbitrary $ \(Nat2 mn) ->
    forAll (Test.matrix mn) $ \x ->
    forAll (Test.matrix mn) $ \y ->
        runST ( do
            x' <- thawMatrix x
            y' <- thawMatrix y
            pre x' y') ==>
        runST ( do
            x' <- unsafeThawMatrix x
            y' <- unsafeThawMatrix y
            commutes2 x' y' a f )
            
------------------------------------------------------------------------

tests_STMatrix =
    [ testProperty "newMatrix" prop_NewMatrix
    , testProperty "newListMatrix" prop_NewListMatrix
    
    , testProperty "getSize" prop_GetSize
    , testProperty "readElem" prop_ReadElem
    , testProperty "canModifyElem" prop_CanModifyElem
    , testProperty "writeElem" prop_WriteElem
    , testProperty "modifyElem" prop_ModifyElem

    , testProperty "getIndices" prop_GetIndicesLazy
    , testProperty "getIndices'" prop_GetIndicesStrict
    , testProperty "getElems" prop_GetElemsLazy
    , testProperty "getElems'" prop_GetElemsStrict

    , testProperty "getElems . modifyWith" prop_GetElemsLazyModifyWith
    , testProperty "getElems' . modifyWith" prop_GetElemsStrictModifyWith
    , testProperty "getAssocs . modifyWith" prop_GetAssocsLazyModifyWith
    , testProperty "getAssocs' . modifyWith" prop_GetAssocsStrictModifyWith

    , testProperty "newZeroMatrix" prop_NewZeroMatrix
    , testProperty "setZeroMatrix" prop_SetZeroMatrix
    , testProperty "newConstantMatrix" prop_NewConstantMatrix
    , testProperty "setConstantMatrix" prop_SetConstantMatrix
    , testProperty "newIdentityMatrix" prop_NewIdentityMatrix
    , testProperty "setIdentityMatrix" prop_SetIdentityMatrix
    
    , testProperty "newCopyMatrix" prop_NewCopyMatrix
    , testProperty "copyMatrix" prop_CopyMatrix
    , testProperty "swapMatrix" prop_SwapMatrix

    , testProperty "doConj" prop_DoConj
    , testProperty "scaleBy" prop_ScaleBy
    , testProperty "shiftBy" prop_ShiftBy
    , testProperty "modifyWith" prop_ModifyWith
    , testProperty "rank1UpdateMatrix" prop_Rank1UpdateMatrix
    
    , testProperty "getConjMatrix" prop_GetConjMatrix
    , testProperty "getScaledMatrix" prop_GetScaledMatrix
    , testProperty "getShiftedMatrix" prop_GetShiftedMatrix

    , testProperty "axpyMatrix" prop_AxpyMatrix
    , testProperty "addMatrix" prop_AddMatrix
    , testProperty "subMatrix" prop_SubMatrix
    , testProperty "mulMatrix" prop_MulMatrix
    , testProperty "divMatrix" prop_DivMatrix
    
    , testProperty "getAddMatrix" prop_GetAddMatrix
    , testProperty "getSubMatrix" prop_GetSubMatrix
    , testProperty "getMulMatrix" prop_GetMulMatrix
    , testProperty "getDivMatrix" prop_GetDivMatrix

    ]
