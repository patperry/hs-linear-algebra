{-# LANGUAGE FlexibleInstances #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.BLAS
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
--
-- Test generators for BLAS types.
--

module Test.QuickCheck.BLAS (
    -- * Testable element types
    TestElem(..),
    
    -- * Generating random objects
    -- ** Shapes
    dim,
    shape,

    -- ** Indices
    index,
    index2,
    
    -- ** Elements
    elem,
    realElem,
    elems,
    realElems,
    
    -- ** Association lists
    assocs,
    assocs2,
    bandedAssocs,
    
    -- ** Vectors
    vector,
    
    -- ** Dense matrices
    matrix,
    hermMatrix,
    triMatrix,

    -- ** Banded matrices
    bandwidths,
    banded,
    bandedWith,
    hermBanded,
    triBanded,
    
    -- * Convenience types
    Nat(..),
    Nat2(..),
    Pos(..),
    Pos2(..),
    Index(..),
    Index2(..),
    Assocs(..),
    Assocs2(..),
    BandedAssocs(..),
    
    ) where

import Prelude hiding ( elem )

import BLAS.Types( UpLoEnum(..), DiagEnum(..) )
import Control.Monad
import Data.Maybe( fromJust, mapMaybe )

import Test.QuickCheck hiding ( vector, elements )
import qualified Test.QuickCheck as QC

import Data.Vector.Dense( Vector, listVector, subvectorWithStride,
    conj )
import Data.Matrix.Dense( Matrix, listMatrix, herm, submatrix  )
import Data.Matrix.Dense.ST( runSTMatrix, setElems, diagView,
    unsafeThawMatrix )
import Data.Matrix.Banded( Banded, maybeBandedFromMatrixStorage )
import Data.Matrix.Banded.ST( runSTBanded, unsafeThawBanded, 
    diagViewBanded )
import Data.Matrix.Herm
import Data.Matrix.Tri
import Data.Elem.BLAS

-- | Element types that can be tested with QuickCheck properties.
class (Elem e, Arbitrary e, CoArbitrary e) => TestElem e where
    -- | Inicates whether or not the value should be used in tests.  For
    -- 'Double's, @isTestElemElem e@ is defined as 
    -- @not (isNaN e || isInfinite e || isDenormalized e)@.
    isTestElemElem :: e -> Bool
    
instance TestElem Double where
    isTestElemElem e = not (isNaN e || isInfinite e || isDenormalized e)
    {-# INLINE isTestElemElem #-}

instance Arbitrary (Complex Double) where
    arbitrary = liftM2 (:+) arbitrary arbitrary
    {-# INLINE arbitrary #-}

instance CoArbitrary (Complex Double) where
    coarbitrary (x:+y) = coarbitrary (x,y)
    {-# INLINE coarbitrary #-}

instance TestElem (Complex Double) where
    isTestElemElem (x :+ y) = isTestElemElem x && isTestElemElem y
    {-# INLINE isTestElemElem #-}

-- | Generate a random element.
elem :: (TestElem e) => Gen e
elem = do
    e <- arbitrary
    if isTestElemElem e then return e
                        else elem

-- | Generate a random element that has no imaginary part.
realElem :: (TestElem e) => Gen e
realElem = liftM fromReal elem

-- | Generate a list of elements suitable for testing with.
elems :: (TestElem e) => Int -> Gen [e]
elems n = replicateM n elem

-- | Generate a list of elements for testing that have no imaginary part.
realElems :: (TestElem e) => Int -> Gen [e]
realElems n = replicateM n realElem

-- | Get an appropriate dimension for a random vector
dim :: Gen Int
dim = sized $ \n -> 
      resize (n `div` 4) $ liftM abs arbitrary

-- | Given a dimension generate a valid index.  The dimension must be positive.
index :: Int -> Gen Int
index n | n <= 0 = 
            error $  "index " ++ (show n) ++ ":"
                  ++ " dimension must be positive (QuickCheck error)"
        | otherwise =
            choose (0,n-1)

-- | Generate a random vector of the given size.
vector :: (TestElem e) => Int -> Gen (Vector e)
vector n =
    frequency [ (3, rawVector n)  
              , (2, conjVector n)
              , (1, subVector n    >>= \(SubVector s x o _) -> 
                    return $ subvectorWithStride s x o n)
              ]    

data SubVector e = 
    SubVector Int 
              (Vector e) 
              Int 
              Int 
    deriving (Show)

instance (TestElem e) => Arbitrary (SubVector e) where
    arbitrary = sized $ \m -> 
        choose (0,m) >>= subVector
        
rawVector :: (TestElem e) => Int -> Gen (Vector e)
rawVector n = do
    es <- elems n
    return $ listVector n es

conjVector :: (TestElem e) => Int -> Gen (Vector e)
conjVector n = do
    x <- vector n
    return $ (conj x)

subVector :: (TestElem e) => Int -> Gen (SubVector e)
subVector n = do
    o <- choose (0,5)
    s <- choose (1,5)
    e <- choose (0,5)
    x <- vector (o + s*n + e)
    return (SubVector s x o n)

data SubMatrix e = 
    SubMatrix (Matrix e) 
              (Int,Int) 
              (Int,Int) 
    deriving (Show)

-- | Generate an appropriate shape for a random matrix.
shape :: Gen (Int,Int)
shape = sized $ \s ->
    let s' = (ceiling . sqrt) (fromIntegral s :: Double)
    in liftM2 (,) (choose (0,s')) (choose (0,1+s'))

-- | Given a shape generate a valid index.  The shape must be have nonzero size.
index2 :: (Int,Int) -> Gen (Int,Int)
index2 (m,n) = liftM2 (,) (index m) (index n)

-- | Generate a random matrix of the given shape.
matrix :: (TestElem e) => (Int,Int) -> Gen (Matrix e)
matrix mn = frequency [ (3, rawMatrix mn)  
                      , (2, hermedMatrix mn)
                      , (1, subMatrix mn >>= \(SubMatrix a ij _) -> 
                                 return $ submatrix a ij mn)
                      ]

-- | Generate a triangular dense matrix.
triMatrix :: (TestElem e) => (Int,Int) -> Gen (Tri Matrix e)
triMatrix (m,n) = do
    a <- matrix (m,n)
    u <- QC.elements [ Lower, Upper  ]
    d <- QC.elements [ Unit, NonUnit ]
    return $ Tri u d a
    
-- | Generate a Hermitian dense matrix.
hermMatrix :: (TestElem e) => Int -> Gen (Herm Matrix e)
hermMatrix n = do
    a <- matrix (n,n)
    d <- realElems n
    let a' = runSTMatrix $ do
                 ma <- unsafeThawMatrix a
                 setElems (diagView ma 0) d
                 return ma
    u <- QC.elements [ Lower, Upper ]
    return $ Herm u a'

rawMatrix :: (TestElem e) => (Int,Int) -> Gen (Matrix e)
rawMatrix (m,n) = do
    es <- elems (m*n)
    return $ listMatrix (m,n) es

hermedMatrix :: (TestElem e) => (Int,Int) -> Gen (Matrix e)
hermedMatrix (m,n) = do
    x <- matrix (n,m)
    return $ (herm x)

subMatrix :: (TestElem e) => (Int,Int) -> Gen (SubMatrix e)
subMatrix (m,n) = 
    oneof [ rawSubMatrix (m,n)
          , rawSubMatrix (n,m) >>= \(SubMatrix a (i,j) (m',n')) ->
                return $ SubMatrix (herm a) (j,i) (n',m')
          ]

rawSubMatrix :: (TestElem e) => (Int,Int) -> Gen (SubMatrix e)
rawSubMatrix (m,n) = do
    i <- choose (0,5)
    j <- choose (0,5)
    e <- choose (0,5)
    f <- choose (0,5)
    x <- rawMatrix (i+m+e, j+n+f)
    return $ SubMatrix x (i,j) (m,n)

instance (TestElem e) => Arbitrary (SubMatrix e) where
    arbitrary = do
        (m,n) <- shape
        (SubMatrix a ij mn) <- subMatrix (m,n)
        return $ SubMatrix a ij mn

-- | Generate valid bandwidth for a given matrix dimension size
bandwidth :: Int -> Gen Int
bandwidth n = if n == 0 then return 0 else choose (0,n-1)
    
-- | Generate valid bandwidths for the given matrix shape.
bandwidths :: (Int,Int) -> Gen (Int,Int)
bandwidths (m,n) = liftM2 (,) (bandwidth m) (bandwidth n)

-- | Generate a random banded matrix of the given shape.
banded :: (TestElem e) => (Int,Int) -> Gen (Banded e)
banded mn = do
    lu <- bandwidths mn
    bandedWith lu mn

-- | Generate a random banded matrix with the given bandwidths.
bandedWith :: (TestElem e) 
           => (Int,Int) -> (Int,Int) -> Gen (Banded e)
bandedWith lu mn = frequency [ (3, rawBanded mn lu)  
                             , (2, hermedBanded mn lu)
                             ]

-- | Generate a triangular banded matrix.
triBanded :: (TestElem e) => Int -> Gen (Tri Banded e)
triBanded n = do
    a <- banded (n,n)
    u <- QC.elements [ Lower, Upper  ]
    d <- QC.elements [ Unit, NonUnit ]
    return $ Tri u d a
    
-- | Generate a Hermitian banded matrix.
hermBanded :: (TestElem e) => Int -> Gen (Herm Banded e)
hermBanded n = do
    a <- banded (n,n)
    d <- realElems n
    let a' = runSTBanded $ do
                 ma <- unsafeThawBanded a
                 setElems (diagViewBanded ma 0) d
                 return ma
    u <- QC.elements [ Lower, Upper ]
    return $ Herm u a'    

rawBanded :: (TestElem e) => 
    (Int,Int) -> (Int,Int) -> Gen (Banded e)
rawBanded (m,n) (kl,ku) = 
    let bw = kl+ku+1
    in do
        a <- frequency [ (2, rawMatrix (bw,n))
                       , (1, rawSubMatrix (bw,n) >>= \(SubMatrix b ij _) ->
                                 return $ submatrix b ij (bw,n))     
                       ]
        return $ fromJust (maybeBandedFromMatrixStorage (m,n) (kl,ku) a)

hermedBanded :: (TestElem e) => 
    (Int,Int) -> (Int,Int) -> Gen (Banded e)
hermedBanded (m,n) (kl,ku) = do
    x <- rawBanded (n,m) (ku,kl)
    return $ herm x

-- | Generate an associations list for a vector of the given dimension.
assocs :: (TestElem e) => Int -> Gen [(Int,e)]
assocs n | n == 0    = return []
         | otherwise = do
    (Nat l) <- arbitrary
    is <- replicateM l $ index n
    es <- elems l
    return $ zip is es
    
-- | Generate an associations list for a matrix of the given shape.
assocs2 :: (TestElem e) => (Int,Int) -> Gen [((Int,Int),e)]
assocs2 (m,n) | m*n == 0  = return []
              | otherwise = do
    (Nat l) <- arbitrary
    is <- replicateM l $ index2 (m,n)
    es <- elems l
    return $ zip is es

-- | Generate an associations list for a banded matrix of the given shape
-- and bandwidths.
bandedAssocs :: (TestElem e) => (Int,Int) -> (Int,Int) -> Gen [((Int,Int),e)]
bandedAssocs (m,n) (kl,ku) | m*n == 0  = return []
                           | otherwise = do
    (Nat l) <- arbitrary
    ijs     <- replicateM l $ index2 (kl+1+ku,n)
    let ijs' = mapMaybe (\(i,j) -> let i' = i - j - ku in
                                   if 0 <= i' && i' < m then Just (i',j)
                                                        else Nothing    ) ijs
    es      <- replicateM l elem
    return $ zip ijs' es

-- | A non-negative integer.
newtype Nat = Nat Int deriving (Eq,Show)
instance Arbitrary Nat where
    arbitrary           = liftM Nat dim

instance CoArbitrary Nat where    
    coarbitrary (Nat n) = coarbitrary n


-- | A pair of non-negative integers
newtype Nat2 = Nat2 (Int,Int) deriving (Eq,Show)
instance Arbitrary Nat2 where
    arbitrary            = liftM Nat2 shape

instance CoArbitrary Nat2 where
    coarbitrary (Nat2 n) = coarbitrary n


-- | A positive integer.
newtype Pos = Pos Int deriving (Eq,Show)
instance Arbitrary Pos where
    arbitrary           = liftM (Pos . (1+)) dim

instance CoArbitrary Pos where
    coarbitrary (Pos n) = coarbitrary n


-- | A pair of positive integers
newtype Pos2 = Pos2 (Int,Int) deriving (Eq,Show)
instance Arbitrary Pos2 where
    arbitrary = do
        (m,n) <- shape
        return $ Pos2 (m+1,n+1)
        
instance CoArbitrary Pos2 where        
    coarbitrary (Pos2 mn) =
        coarbitrary mn


-- | A dimension and a valid index for it.
data Index = Index Int Int deriving (Eq,Show)
instance Arbitrary Index where
    arbitrary = do
        n <- dim
        let n' = n+1
        i <- index n'
        return $ Index n' i

instance CoArbitrary Index where        
    coarbitrary (Index n i) = 
        coarbitrary (n,i)


-- | A shape and a valid index for it.
data Index2 = Index2 (Int,Int) (Int,Int) deriving (Eq,Show)
instance Arbitrary Index2 where
    arbitrary = do
        (m,n) <- shape
        let mn' = (m+1,n+1)
        ij <- index2 mn'
        return $ Index2 mn' ij

instance CoArbitrary Index2 where        
    coarbitrary (Index2 mn ij) = 
        coarbitrary (mn,ij)


-- | A dimension and an associations list.
data Assocs e = Assocs Int [(Int,e)] deriving (Eq,Show)
instance (TestElem e) => Arbitrary (Assocs e) where
    arbitrary = do
        n   <- dim
        ies <- assocs n
        return $ Assocs n ies
        
instance (TestElem e) => CoArbitrary (Assocs e) where        
    coarbitrary (Assocs n ies) =
        coarbitrary (n,ies)

        
-- | A shape and an associations list.
data Assocs2 e = Assocs2 (Int,Int) [((Int,Int),e)] deriving (Eq,Show)
instance (TestElem e) => Arbitrary (Assocs2 e) where
    arbitrary = do
        n    <- shape
        ies  <- assocs2 n
        return $ Assocs2 n ies

instance (TestElem e) => CoArbitrary (Assocs2 e) where
    coarbitrary (Assocs2 n ies) =
        coarbitrary (n,ies)


-- | A shape, bandwidths, and an associations list.
data BandedAssocs e = BandedAssocs (Int,Int) (Int,Int) [((Int,Int),e)] deriving (Eq,Show)
instance (TestElem e) => Arbitrary (BandedAssocs e) where
    arbitrary = do
        mn  <- shape
        bw  <- bandwidths mn
        ies <- bandedAssocs mn bw
        return $ BandedAssocs mn bw ies

instance (TestElem e) => CoArbitrary (BandedAssocs e) where
    coarbitrary (BandedAssocs mn bw ies) =
        coarbitrary (mn,bw,ies)
