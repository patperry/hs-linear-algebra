{-# LANGUAGE FlexibleInstances, GeneralizedNewtypeDeriving #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.LinearAlgebra
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
--
-- Test generators for linear algebra types.
--

module Test.QuickCheck.LinearAlgebra (
    -- Testable element types
    TestElem(..),
    
    -- * Generating random objects
    -- ** Dimensions
    dim,
    dim2,
    Dim(..),
    Dim2(..),
    
    -- ** Indices
    index,
    index2,    
    Index(..),
    Index2(..),     
    
    -- ** Elements
    elem,
    elems,
    -- realElem,
    -- realElems,
    
    -- ** Association lists
    assocs,
    assocs2,
    Assocs(..),    
    Assocs2(..),
    -- bandedAssocs,
    
    -- ** Vectors
    vector,
    VectorPair(..),
    VectorTriple(..),
    VectorList(..),
    WeightedVectorList(..),
    NonEmptyVectorList(..),
    NonEmptyWeightedVectorList(..),
    
    --  ** Matrices
    matrix,
    MatrixPair(..),
    MatrixTriple(..),
    -- hermMatrix,
    -- triMatrix,

    --  Banded matrices
    -- bandwidths,
    -- banded,
    -- bandedWith,
    -- hermBanded,
    -- triBanded,
    
    -- BandedAssocs(..),

    ) where

import Prelude hiding ( elem )

-- import BLAS.Types( UpLoEnum(..), DiagEnum(..) )
import Control.Monad
import Data.Complex( magnitude )
import Data.Maybe( fromJust )

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Vector( Vector )
import qualified Numeric.LinearAlgebra.Vector as V
import Numeric.LinearAlgebra.Matrix( Matrix )
import qualified Numeric.LinearAlgebra.Matrix as M
-- import Data.Matrix.Banded( Banded, maybeBandedFromMatrixStorage )
-- import Data.Matrix.Banded.ST( runSTBanded, unsafeThawBanded, 
--     diagViewBanded )


class TestElem e where
    maybeToReal :: e -> Maybe Double
    toReal :: e -> Double
    norm1 :: e -> Double
    norm :: e -> Double
    conjugate :: e -> e
    
    toReal = fromJust . maybeToReal
    
instance (TestElem Double) where
    maybeToReal = Just
    toReal = id
    norm1 = abs
    norm = abs
    conjugate = id
    
instance (TestElem (Complex Double)) where
    maybeToReal (x :+ y) | y == 0    = Just x
                         | otherwise = Nothing
    norm1 (x :+ y) = abs x + abs y
    norm z = magnitude z
    conjugate (x :+ y) = x :+ (-y)

{-
-- | Element types that can be tested with QuickCheck properties.
class (Storable e, Arbitrary e, CoArbitrary e) => TestElem e where
    -- | Inicates whether or not the value should be used in tests.  For
    -- 'Double's, @isTestElem e@ is defined as 
    -- @not (isNaN e || isInfinite e || isDenormalized e)@.
    isTestElem :: e -> Bool

instance TestElem Double where
    isTestElem e = not (isNaN e || isInfinite e || isDenormalized e)
    {-# INLINE isTestElem #-}

instance TestElem (Complex Double) where
    isTestElem (x :+ y) = isTestElem x && isTestElem y
    {-# INLINE isTestElem #-}
-}

-- | Generate a random element.
elem :: (Arbitrary e) => Gen e
elem = arbitrary

-- | Generate a list of elements suitable for testing with.
elems :: (Arbitrary e) => Int -> Gen [e]
elems n = replicateM n elem

-- Generate a random element that has no imaginary part.
-- realElem :: (TestElem e) => Gen e
-- realElem = liftM fromReal elem

-- Generate a list of elements for testing that have no imaginary part.
-- realElems :: (TestElem e) => Int -> Gen [e]
-- realElems n = replicateM n realElem

-- | Get an appropriate dimension for a random vector
dim :: Gen Int
dim = sized $ \s -> do
    (NonNegative n) <- resize (s `div` 4) $ arbitrary
    return n
 
-- | Get an appropriate dimension for a random matrix
dim2 :: Gen (Int,Int)
dim2 = sized $ \s -> do
    m <- resize (s `div` 2) $ dim
    n <- resize (s `div` 2) $ dim
    return (m,n)
 
-- | A vector dimension.   
newtype Dim = Dim Int
  deriving (Eq, Ord, Num, Integral, Real, Enum, Show, Read)

instance Arbitrary Dim where
    arbitrary = sized $ \s -> do
        (NonNegative n) <- resize (s `div` 4) $ arbitrary
        return $ Dim n
    
    shrink (Dim a) = 
        [ Dim a' | (NonNegative a') <- shrink (NonNegative a) ]

-- | A matrix dimension.   
newtype Dim2 = Dim2 (Int,Int)
  deriving (Eq, Show, Read)

instance Arbitrary Dim2 where
    arbitrary = do
        mn <- dim2
        return $ Dim2 mn


-- | Given a dimension generate a valid index.  The dimension must be positive.
index :: Int -> Gen Int
index n | n <= 0 = 
            error $ "index " ++ (show n) ++ ":"
                  ++ " dimension must be positive (QuickCheck error)"
        | otherwise =
            choose (0,n-1)

-- | Given a matrix dimension generate a valid index.
index2 :: (Int,Int) -> Gen (Int,Int)
index2 (m,n) = do
    i <- index m
    j <- index n
    return (i,j)

-- | A dimension and a valid index for it.
data Index = Index Int Int deriving (Eq,Show)
instance Arbitrary Index where
    arbitrary = do
        n <- (1+) `fmap` dim
        i <- index n
        return $ Index n i

-- | A matrix dimension and a valid index for it.
data Index2 = Index2 (Int,Int) (Int,Int) deriving (Eq,Show)
instance Arbitrary Index2 where
    arbitrary = do
        (m',n') <- dim2
        let mn = (m' + 1, n' + 1)
        ij <- index2 mn
        return $ Index2 mn ij
    
    
-- | Generate an associations list for a vector of the given dimension.
assocs :: (Arbitrary e) => Int -> Gen [(Int,e)]
assocs n | n == 0    = return []
         | otherwise = do
    l <- choose(0, 2*n)
    is <- replicateM l $ index n
    es <- elems l
    return $ zip is es

-- | Generate an associations list for a matrix of the given shape.
assocs2 :: (Arbitrary e) => (Int,Int) -> Gen [((Int,Int),e)]
assocs2 (m,n) | m*n == 0  = return []
              | otherwise = do
    l <- choose(0, 2*m*n)
    is <- replicateM l $ index2 (m,n)
    es <- elems l
    return $ zip is es


-- | A dimension and an associations list.
data Assocs e = Assocs Int [(Int,e)] deriving (Eq,Show)

instance (Arbitrary e) => Arbitrary (Assocs e) where
    arbitrary = do
        n   <- dim
        ies <- assocs n
        return $ Assocs n ies
    
    shrink (Assocs n ies) =
        [ Assocs n' $ filter ((< n') . fst) ies
        | n' <- shrink n
        ] ++
        [ Assocs n ies'
        | ies' <- shrink ies
        ]

-- | A shape and an associations list.
data Assocs2 e = Assocs2 (Int,Int) [((Int,Int),e)] deriving (Eq,Show)
instance (Arbitrary e) => Arbitrary (Assocs2 e) where
    arbitrary = do
        mn   <- dim2
        ies  <- assocs2 mn
        return $ Assocs2 mn ies


-- | Generate a random vector of the given size.
vector :: (Arbitrary e, Storable e) => Int -> Gen (Vector e)
vector n = do
    es <- elems n
    return $ V.fromList n es

instance (Arbitrary e, Storable e) => Arbitrary (Vector e) where
    arbitrary = dim >>= vector
    
    shrink x =
        [ V.slice 0 n x
        | (NonNegative n) <- shrink (NonNegative $ V.dim x)
        ]

-- | Two vectors with the same dimension.
data VectorPair e f = 
    VectorPair (Vector e) (Vector f) deriving (Eq, Show)
instance (Arbitrary e, Storable e, Arbitrary f, Storable f) =>
    Arbitrary (VectorPair e f) where
        arbitrary = do
            x <- arbitrary
            y <- vector (V.dim x)
            return $ VectorPair x y
            
        shrink (VectorPair x y) =
            [ VectorPair (V.slice 0 n' x) (V.slice 0 n' y)
            | n' <- shrink (V.dim x)
            ]

-- | Three vectors with the same dimension.
data VectorTriple e f g =
    VectorTriple (Vector e) (Vector f) (Vector g) deriving (Eq, Show)
instance (Arbitrary e, Storable e, Arbitrary f, Storable f,
          Arbitrary g, Storable g) =>
    Arbitrary (VectorTriple e f g) where
        arbitrary = do
            x <- arbitrary
            y <- vector (V.dim x)
            z <- vector (V.dim x)
            return $ VectorTriple x y z

-- | A nonempty list of vectors with the same dimension.
data NonEmptyVectorList e = NonEmptyVectorList Int [Vector e] deriving (Eq, Show)
instance (Arbitrary e, Storable e) => Arbitrary (NonEmptyVectorList e) where
    arbitrary = do
        x <- arbitrary
        n <- choose (0,20)
        let p = V.dim x
        xs <- replicateM n $ vector p
        return $ NonEmptyVectorList p $ x:xs

-- | A nonempty list of (weight, vector) pairs, with the weights all non-negative
-- and the vectors all having the same dimension.
data NonEmptyWeightedVectorList e = NonEmptyWeightedVectorList Int [(e, Vector e)]
    deriving (Eq, Show)
instance (Arbitrary e, Storable e, Num e) => Arbitrary (NonEmptyWeightedVectorList e) where
    arbitrary = do
        (NonEmptyVectorList p xs) <- arbitrary
        ws <- replicateM (length xs) $ fmap abs arbitrary
        return $ NonEmptyWeightedVectorList p $ zip ws xs
        
-- | A list of vectors with the same dimension.
data VectorList e = VectorList Int [Vector e] deriving (Eq, Show)
instance (Arbitrary e, Storable e) => Arbitrary (VectorList e) where
    arbitrary = do
        (NonEmptyVectorList p (_:xs)) <- arbitrary
        return $ VectorList p xs

-- | A list of (weight, vector) pairs, with the weights all non-negative
-- and the vectors all having the same dimension.
data WeightedVectorList e = WeightedVectorList Int [(e, Vector e)]
    deriving (Eq, Show)
instance (Arbitrary e, Storable e, Num e) => Arbitrary (WeightedVectorList e) where
    arbitrary = do
        (NonEmptyWeightedVectorList p (_:wxs)) <- arbitrary
        return $ WeightedVectorList p wxs

-- | Generate a random matrix of the given size.
matrix :: (Arbitrary e, Storable e) => (Int,Int) -> Gen (Matrix e)
matrix (m,n) = 
    oneof [ raw, sub ]
  where
    raw = do
        es <- elems $ m * n
        return $ M.fromList (m,n) es
    sub = do
        m' <- choose (m, 2*m)
        es <- elems $ m' * n
        return $ M.slice (0,0) (m,n) (M.fromList (m',n) es)

instance (Arbitrary e, Storable e) => Arbitrary (Matrix e) where
    arbitrary = dim2 >>= matrix
    
-- | Two matrices with the same dimension.
data MatrixPair e f = 
    MatrixPair (Matrix e) (Matrix f) deriving (Eq, Show)
instance (Arbitrary e, Storable e, Arbitrary f, Storable f) =>
    Arbitrary (MatrixPair e f) where
        arbitrary = do
            x <- arbitrary
            y <- matrix (M.dim x)
            return $ MatrixPair x y

-- | Three matrices with the same dimension.
data MatrixTriple e f g =
    MatrixTriple (Matrix e) (Matrix f) (Matrix g) deriving (Eq, Show)
instance (Arbitrary e, Storable e, Arbitrary f, Storable f,
          Arbitrary g, Storable g) =>
    Arbitrary (MatrixTriple e f g) where
        arbitrary = do
            x <- arbitrary
            y <- matrix (M.dim x)
            z <- matrix (M.dim x)
            return $ MatrixTriple x y z

instance Arbitrary Trans where
    arbitrary = elements [ NoTrans, Trans, ConjTrans ]
{-

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
    return $ M.fromList (m,n) es


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


-- | A shape, bandwidths, and an associations list.
data BandedAssocs e = BandedAssocs (Int,Int) (Int,Int) [((Int,Int),e)] deriving (Eq,Show)
instance (TestElem e) => Arbitrary (BandedAssocs e) where
    arbitrary = do
        mn  <- shape
        bw  <- bandwidths mn
        ies <- bandedAssocs mn bw
        return $ BandedAssocs mn bw ies
-}
