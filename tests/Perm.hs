
module Perm( tests_Perm ) where
    
import Driver

import Data.Vector.Dense
import Data.Matrix.Dense
import Data.Matrix.Perm

import Data.Permutation ( Permutation, permutation )
import qualified Data.Permutation as P

import Generators.Matrix.Perm

type V = Vector Int E
type M = Matrix (Int,Int) E
type P = Perm (Int,Int) E

prop_perm_herm (TestPerm (p :: P)) =
    toPermutation (herm p) == P.inverse (toPermutation p)

prop_perm_apply_basis (PermMBasis (p :: P) i) =
    n > 0 ==> p <*> (basisVector n i) === basisVector n (P.apply (toPermutation p) i)
  where
    n = numCols p

prop_perm_herm_apply (PermMV (p :: P) x) =
    p <*> herm p <*> x === x

prop_herm_perm_apply (PermMV (p :: P) x) =
    herm p <*> p <*> x === x

prop_perm_solve (PermMV (p :: P) x) =
    p <\> x === herm p <*> x

prop_perm_applyMat_cols (PermMM (p :: P) a) =
    cols (p <**> a) === map (p <*>) (cols a)

prop_perm_herm_applyMat (PermMM (p :: P) a) =
    p <**> herm p <**> a === a

prop_herm_perm_applyMat (PermMM (p :: P) a) =
    herm p <**> p <**> a === a

prop_perm_solveMat_cols (PermMM (p :: P) a) =
    cols (p <\\> a) === map (p <\>) (cols a)

prop_perm_solveMat (PermMM (p :: P) a) =
    p <\\> a === herm p <**> a
    

tests_Perm =
    [ ("perm herm"             , mytest prop_perm_herm)
    , ("perm apply basis"      , mytest prop_perm_apply_basis)
    , ("perm herm apply"       , mytest prop_perm_herm_apply)
    , ("herm perm apply"       , mytest prop_herm_perm_apply)
    , ("perm solve"            , mytest prop_perm_solve)
    , ("perm applyMat cols"    , mytest prop_perm_applyMat_cols)
    , ("perm herm applyMat"    , mytest prop_perm_herm_applyMat)
    , ("herm perm applyMat"    , mytest prop_herm_perm_applyMat)
    , ("perm solveMat cols"    , mytest prop_perm_solveMat_cols)
    , ("perm solveMat"         , mytest prop_perm_solveMat)
    ]
