{-# LANGUAGE Rank2Types, PatternSignatures #-}
module STMatrix ( tests_STMatrix ) where


import BLAS.Elem
import Data.Matrix.Dense
import Data.Matrix.Dense.ST

import qualified Generators.Matrix.Dense as Test

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
prop_ReadElem (Index2 ij mn) =
    implementsFor mn (`readElem` ij) (`readElem_S` ij)

canModifyElem_S a ij = ( True, a )
prop_CanModifyElem ij = (`canModifyElem` ij) `implements` (`canModifyElem_S` ij)

writeElem_S a ij e = ( (), a // [(ij,e)] )
prop_WriteElem (Index2 ij mn) e =
    implementsFor mn (\a -> writeElem a ij e) (\a -> writeElem_S a ij e)

modifyElem_S a ij f = writeElem_S a ij $ f (a!ij)
prop_ModifyElem (Index2 ij mn) f =
    implementsFor mn (\a -> modifyElem a ij f) (\a -> modifyElem_S a ij f)

getIndices_S a = ( indices a, a )

prop_GetIndicesLazy   = getIndices `implements`  getIndices_S
prop_GetIndicesStrict = getIndices' `implements` getIndices_S
    
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


----------------------------- Copying Vectors --------------------------------

{-
newCopyVector_S x = ( x, x )
prop_NewCopyVector = 
    (\x -> newCopyVector x >>= abstract) `implements` newCopyVector_S

copyVector_S x y = ( (), y, y )
prop_CopyVector = copyVector `implements2` copyVector_S

swap_S x y = ( (), y, x )
prop_Swap = swap `implements2` swap_S


-------------------------- Unsary Vector Operations --------------------------

doConj_S x = ( (), map conj x )
prop_DoConj = doConj `implements` doConj_S

scaleBy_S k x = ( (), map (k*) x )
prop_ScaleBy k = scaleBy k `implements` scaleBy_S k

shiftBy_S k x = ( (), map (k+) x )
prop_ShiftBy k = shiftBy k `implements` shiftBy_S k

modifyWith_S f x = ( (), map f x )
prop_ModifyWith f = modifyWith f `implements` modifyWith_S f

getConj_S x = ( map conj x, x )
prop_GetConj = 
    (\x -> getConj x >>= abstract) `implements` getConj_S

getScaled_S k x = ( map (k*) x, x )
prop_GetScaled k = 
    (\x -> getScaled k x >>= abstract) `implements` (getScaled_S k)

getShifted_S k x = ( map (k+) x, x )
prop_GetShifted k = 
    (\x -> getShifted k x >>= abstract) `implements` (getShifted_S k)


------------------------- Binary Vector Operations ---------------------------

addEquals_S x y = ( (), zipWith (+) x y, y )
prop_AddEquals = (+=) `implements2` addEquals_S

subEquals_S x y = ( (), zipWith (-) x y, y )
prop_SubEquals = (-=) `implements2` subEquals_S

mulEquals_S x y = ( (), zipWith (*) x y, y )
prop_MulEquals = (*=) `implements2` mulEquals_S

divEquals_S x y = ( (), zipWith (/) x y, y )
prop_DivEquals = (//=) `implements2` divEquals_S

axpy_S alpha x y = ( (), x, zipWith (\xi yi -> alpha * xi + yi) x y )
prop_Axpy alpha = axpy alpha `implements2` axpy_S alpha

getAdd_S x y = ( zipWith (+) x y, x, y )
prop_GetAdd =
    (\x y -> getAdd x y >>= abstract) `implements2` getAdd_S

getSub_S x y = ( zipWith (-) x y, x, y )
prop_GetSub =
    (\x y -> getSub x y >>= abstract) `implements2` getSub_S

getMul_S x y = ( zipWith (*) x y, x, y )
prop_GetMul =
    (\x y -> getMul x y >>= abstract) `implements2` getMul_S

getDiv_S x y = ( zipWith (/) x y, x, y )
prop_GetDiv =
    (\x y -> getDiv x y >>= abstract) `implements2` getDiv_S


-------------------------- Vector Properties ---------------------------------

getSumAbs_S x = ( foldl' (+) 0 $ map norm1 x, x )
prop_GetSumAbs = getSumAbs `implements` getSumAbs_S

getNorm2_S x = ( sqrt $ foldl' (+) 0 $ map (**2) $ map norm x, x )
prop_GetNorm2 = getNorm2 `implements` getNorm2_S

getWhichMaxAbs_S x = ( (i,x!!i), x) 
    where i = fst $ maximumBy (comparing snd) $ reverse $ zip [0..] $ map norm1 x
prop_GetWhichMaxAbs = 
    implementsIf (return . (>0) . dim) getWhichMaxAbs getWhichMaxAbs_S

getDot_S x y = ( foldl' (+) 0 $ zipWith (*) (map conj x) y, x, y )
prop_GetDot = getDot `implements2` getDot_S

-}
------------------------------------------------------------------------
-- 
-- The specification language
--
    
abstract :: (BLAS1 e) => STMatrix s (m,n) e -> ST s (Matrix (m,n) e)
abstract = freezeMatrix

commutes :: (AEq a, Show a, AEq e, BLAS1 e) =>
    STMatrix s (m,n) e -> (STMatrix s (m,n) e -> ST s a) ->
        (Matrix (m,n) e -> (a,Matrix (m,n) e)) -> ST s Bool
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

{-
commutes2 :: (AEq a, Show a, AEq e, Elem e) =>
    STVector s n e -> STVector s n e -> 
    (STVector s n e ->  STVector s n e -> ST s a) ->
        ([e] -> [e] -> (a,[e],[e])) -> ST s Bool
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
-}

equivalent :: (forall s . ST s (STMatrix s (m,n) E)) -> Matrix (m,n) E -> Bool
equivalent x s = runST $ do
    x' <- (x >>= abstract)
    when (not $ x' === s) $
        trace (printf ("expected `%s' but got `%s'") (show s) (show x'))
            return ()
    return (x' === s)
    
implements :: (AEq a, Show a) =>
    (forall s . STMatrix s (m,n) E -> ST s a) ->
    (Matrix (m,n) E -> (a,Matrix (m,n) E)) -> 
        Property
a `implements` f =
    forAll arbitrary $ \(Nat2 mn) ->
        implementsFor mn a f

{-
implements2 :: (AEq a, Show a) =>
    (forall s . STVector s n E -> STVector s n E -> ST s a) ->
    ([E] -> [E] -> (a,[E],[E])) -> 
        Property
a `implements2` f =
    forAll arbitrary $ \(Nat n) ->
        implementsFor2 n a f
-}

implementsFor :: (AEq a, Show a) =>
    (Int,Int) ->
    (forall s . STMatrix s (m,n) E -> ST s a) ->
    (Matrix (m,n) E -> (a,Matrix (m,n) E)) -> 
        Property
implementsFor mn a f =
    forAll (Test.matrix mn) $ \x ->
        runST $ do
            commutes (unsafeThawMatrix x) a f

{-
implementsFor2 :: (AEq a, Show a) =>
    Int ->
    (forall s . STVector s n E -> STVector s n E -> ST s a) ->
    ([E] -> [E] -> (a,[E],[E])) -> 
        Property
implementsFor2 n a f =
    forAll (Test.vector n) $ \x ->
    forAll (Test.vector n) $ \y ->
        runST $ do
            commutes2 (unsafeThawVector x) (unsafeThawVector y) a f
-}

implementsIf :: (AEq a, Show a) =>
    (forall s . STMatrix s (m,n) E -> ST s Bool) ->
    (forall s . STMatrix s (m,n) E -> ST s a) ->
    (Matrix (m,n) E -> (a,Matrix (m,n) E)) -> 
        Property
implementsIf pre a f =
    forAll arbitrary $ \(Nat2 mn) ->
    forAll (Test.matrix mn) $ \x ->
        runST ( do
            x' <- thawMatrix x
            pre x') ==>
        runST ( do
            commutes (unsafeThawMatrix x) a f )

{-
implementsIf2 :: (AEq a, Show a) =>
    (forall s . STVector s n E -> STVector s n E -> ST s Bool) ->
    (forall s . STVector s n E -> STVector s n E -> ST s a) ->
    ([E] -> [E] -> (a,[E],[E])) -> 
        Property
implementsIf2 pre a f =
    forAll arbitrary $ \(Nat n) ->
    forAll (Test.vector n) $ \x ->
    forAll (Test.vector n) $ \y ->
        runST ( do
            x' <- thawVector x
            y' <- thawVector y
            pre x' y') ==>
        runST ( do
            commutes2 (unsafeThawVector x) (unsafeThawVector y) a f )
-}
            
------------------------------------------------------------------------

tests_STMatrix =
    [ ("newMatrix", mytest prop_NewMatrix)
    , ("newListMatrix", mytest prop_NewListMatrix)
    
    , ("getSize", mytest prop_GetSize)
    , ("readElem", mytest prop_ReadElem)
    , ("canModifyElem", mytest prop_CanModifyElem)
    , ("writeElem", mytest prop_WriteElem)
    , ("modifyElem", mytest prop_ModifyElem)

    , ("getIndices", mytest prop_GetIndicesLazy)
    , ("getIndices'", mytest prop_GetIndicesStrict)
    , ("getElems", mytest prop_GetElemsLazyModifyWith)
    , ("getElems'", mytest prop_GetElemsStrictModifyWith)
    , ("getAssocs", mytest prop_GetAssocsLazyModifyWith)
    , ("getAssocs'", mytest prop_GetAssocsStrictModifyWith)

    , ("newZeroMatrix", mytest prop_NewZeroMatrix)
    , ("setZeroMatrix", mytest prop_SetZeroMatrix)
    , ("newConstantMatrix", mytest prop_NewConstantMatrix)
    , ("setConstantMatrix", mytest prop_SetConstantMatrix)
    , ("newIdentityMatrix", mytest prop_NewIdentityMatrix)
    , ("setIdentityMatrix", mytest prop_SetIdentityMatrix)
    
 {-    
    , ("newCopyVector", mytest prop_NewCopyVector)
    , ("copyVector", mytest prop_CopyVector)
    , ("swap", mytest prop_Swap)

    , ("doConj", mytest prop_DoConj)
    , ("scaleBy", mytest prop_ScaleBy)
    , ("shiftBy", mytest prop_ShiftBy)
    , ("modifyWith", mytest prop_ModifyWith)
    
    , ("getConj", mytest prop_GetConj)
    , ("getScaled", mytest prop_GetScaled)
    , ("getShifted", mytest prop_GetShifted)

    , ("axpy", mytest prop_Axpy)
    , ("(+=)", mytest prop_AddEquals)
    , ("(-=)", mytest prop_SubEquals)
    , ("(*=)", mytest prop_MulEquals)
    , ("(//=)", mytest prop_DivEquals)
    
    , ("getAdd", mytest prop_GetAdd)
    , ("getSub", mytest prop_GetSub)
    , ("getMul", mytest prop_GetMul)
    , ("getDiv", mytest prop_GetDiv)

    , ("getSumAbs", mytest prop_GetSumAbs)
    , ("getNorm2", mytest prop_GetNorm2)
    , ("getWhichMaxAbs", mytest prop_GetWhichMaxAbs)
    , ("getDot", mytest prop_GetDot)
-}
    ]
