-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Statistics
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Basic multivariate statistics.
--

module Numeric.LinearAlgebra.Statistics (
    CovMethod(..),
    defaultCovUplo,

    -- * Immutable interface
    
    -- ** Sums and means
    sumVector,
    meanVector,
    weightedSumVector,
    weightedMeanVector,

    -- ** Covariance matrices
    covMatrix,
    covMatrixWithMean,
    weightedCovMatrix,
    weightedCovMatrixWithMean,

    -- ** Covariance matrices in packed form
    covPacked,
    covPackedWithMean,
    weightedCovPacked,
    weightedCovPackedWithMean,

    -- * Mutable interface
    
    -- ** Sums and means
    sumToVector,
    meanToVector,
    weightedSumToVector,
    weightedMeanToVector,
        
    -- ** Covariance matrices
    covToMatrix,
    covToMatrixWithMean,
    weightedCovToMatrix,
    weightedCovToMatrixWithMean,
    
    -- ** Covariance matrices in packed form
    covToPacked,
    covToPackedWithMean,
    weightedCovToPacked,
    weightedCovToPackedWithMean,

    ) where

import Control.Monad( forM_ )
import Control.Monad.ST( ST )
import Data.List( foldl' )
import Text.Printf( printf )

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Matrix.Herm
import Numeric.LinearAlgebra.Matrix.Packed( Packed, STPacked )
import qualified Numeric.LinearAlgebra.Matrix.Packed as P

import Numeric.LinearAlgebra.Vector( Vector )
import qualified Numeric.LinearAlgebra.Vector as V

import Numeric.LinearAlgebra.Vector.ST( STVector, RVector )
import qualified Numeric.LinearAlgebra.Vector.ST as V

import Numeric.LinearAlgebra.Matrix( Matrix )
import qualified Numeric.LinearAlgebra.Matrix as M

import Numeric.LinearAlgebra.Matrix.ST( STMatrix )
import qualified Numeric.LinearAlgebra.Matrix.ST as M


-- | The method of scaling the sample covariance matrix.
data CovMethod =
      UnbiasedCov -- ^ This is the default behavior. Corresponds to a
                  -- scaling of @n/(n-1)@ in the unweighed case, and
                  -- @1/(1 - \\sum w_i^2)@ in the weighted case, where @w_i@
                  -- is the normalized weight. Note the unweighted and
                  -- weighted cases agree when @w_i = 1/n@.
                  
    | MLCov       -- ^ Returns the centered second moment matrix without
                  -- scaling the result.
    deriving (Eq, Show)

-- | Returns the default storage scheme for covariance matrices.
defaultCovUplo :: Uplo
defaultCovUplo = Lower

-- | Returns the sum of the vectors.  The first argument gives the dimension
-- of the vectors.
sumVector :: (VNum e) => Int -> [Vector e] -> Vector e
sumVector p xs = V.create $ do
    s <- V.new_ p
    sumToVector xs s
    return s

-- | Returns the mean of the vectors.  The first argument gives the dimension
-- of the vectors.
meanVector :: (VNum e, Fractional e) => Int -> [Vector e] -> Vector e
meanVector p xs = V.create $ do
      m <- V.new_ p
      meanToVector xs m
      return m

-- | Returns the weighted sum of the vectors.  The first argument gives the
-- dimension of the vectors.
weightedSumVector :: (VNum e) => Int -> [(e, Vector e)] -> Vector e
weightedSumVector p wxs = V.create $ do
    s <- V.new_ p
    weightedSumToVector wxs s
    return s

-- | Returns the weighted mean of the vectors.  The first argument gives the
-- dimension of the vectors.
weightedMeanVector :: (VNum e, Fractional e) => Int -> [(Double, Vector e)] -> Vector e
weightedMeanVector p wxs = V.create $ do
       s <- V.new_ p
       weightedMeanToVector wxs s
       return s

-- | Sets the target vector to the sum of the vectors.
sumToVector :: (RVector v, VNum e) => [v e] -> STVector s e -> ST s ()
sumToVector = weightedSumToVector . zip (repeat 1)

-- | Sets the target vector to the mean of the vectors.
meanToVector :: (RVector v, VNum e, Fractional e)
             => [v e] -> STVector s e -> ST s()
meanToVector = weightedMeanToVector . zip (repeat 1)

-- | Sets the target vector to the weigthed sum of the vectors.
weightedSumToVector :: (RVector v, VNum e) => [(e, v e)] -> STVector s e -> ST s ()
weightedSumToVector wxs s = do
    err <- V.new n 0
    old_s <- V.new_ n
    diff <- V.new_ n
    val <- V.new_ n
    
    V.setElems s (replicate n 0)
    forM_ wxs $ \(w,x) -> do
        V.unsafeCopyTo s old_s -- old_s := s
        V.scaleTo w x val      -- val := w * x
        V.addTo err val err    -- err := err + val
        V.addTo s err s        -- s := s + err
        
        V.subTo old_s s diff   -- diff := old_s - s
        V.addTo diff val err   -- err := diff + val
  where
    n = V.dim s

-- | Sets the target vector to the weighted mean of the vectors.
weightedMeanToVector :: (RVector v, VNum e, Fractional e)
                     => [(Double, v e)] -> STVector s e -> ST s ()
weightedMeanToVector wxs m = let
    go _ _ [] = return ()
    go diff w_sum ((w,x):wxs') | w == 0    = go diff w_sum wxs'
                               | otherwise = let w_sum' = w_sum + w
                                             in do
                                    V.subTo x m diff
                                    V.addToWithScales
                                        (realToFrac $ w/w_sum') diff 1 m m
                                    go diff w_sum' wxs'
    in do
        diff <- V.new_ n
        V.setElems m (replicate n 0)
        go diff 0 wxs
  where
    n = V.dim m

-- | Returns the sample covariance matrix as a hermitian matrix with storage
-- scheme equal to 'defaultCovUplo'.  The first argument gives the dimension
-- of the vectors.
covMatrix :: (BLAS3 e)
          => Int -> CovMethod -> [Vector e] -> Herm Matrix e
covMatrix p t xs = runHermMatrix $ do
    cov <- Herm uplo `fmap` M.new_ (p,p)
    covToMatrix t xs cov
    return cov
  where
    uplo = defaultCovUplo

-- | Returns the sample covariance matrix hermitian matrix (in packed form)
-- with storage scheme equal to 'defaultCovUplo'.  The first argument gives
-- the dimension of the vectors.
covPacked :: (BLAS2 e)
          => Int -> CovMethod -> [Vector e] -> Herm Packed e
covPacked p t xs = P.hermCreate $ do
    cov <- (Herm uplo . P.unsafeViewFromSTVector p) `fmap` V.new_ (p*(p+1) `div` 2)
    covToPacked t xs cov
    return cov
  where
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the sample covariance matrix
-- with storage scheme equal to 'defaultCovUplo'.
covMatrixWithMean :: (BLAS3 e)
                  => Vector e -> CovMethod -> [Vector e] -> Herm Matrix e
covMatrixWithMean mu t xs = runHermMatrix $ do
    cov <- Herm uplo `fmap` M.new_ (p,p)
    covToMatrixWithMean mu t xs cov
    return cov
  where
    p = V.dim mu
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the sample covariance matrix
-- (in packed form) with storage scheme equal to 'defaultCovUplo'.
covPackedWithMean :: (BLAS2 e)
                  => Vector e -> CovMethod -> [Vector e] -> Herm Packed e
covPackedWithMean mu t xs = P.hermCreate $ do
    cov <- (Herm uplo . P.unsafeViewFromSTVector p) `fmap` V.new_ (p*(p+1) `div` 2)
    covToPackedWithMean mu t xs cov
    return cov
  where
    p = V.dim mu
    uplo = defaultCovUplo

-- | Returns the weighed sample covariance matrix with storage scheme equal
-- to 'defaultCovUplo'. The first argument gives the dimension of the vectors.
weightedCovMatrix :: (BLAS3 e)
                  => Int -> CovMethod -> [(Double, Vector e)] -> Herm Matrix e
weightedCovMatrix p t wxs = runHermMatrix $ do
    cov <- Herm uplo `fmap` M.new_ (p,p)
    weightedCovToMatrix t wxs cov
    return cov
  where
    uplo = defaultCovUplo

-- | Returns the weighed sample covariance matrix (in packed form) with
-- storage scheme equal to 'defaultCovUplo'. The first argument gives the
-- dimension of the vectors.
weightedCovPacked :: (BLAS2 e)
                  => Int -> CovMethod -> [(Double, Vector e)] -> Herm Packed e
weightedCovPacked p t wxs = P.hermCreate $ do
    cov <- (Herm uplo . P.unsafeViewFromSTVector p) `fmap` V.new_ (p*(p+1) `div` 2)
    weightedCovToPacked t wxs cov
    return cov
  where
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the weighed sample covariance matrix
-- with storage scheme equal to 'defaultCovUplo'.
weightedCovMatrixWithMean :: (BLAS3 e)
                          => Vector e -> CovMethod -> [(Double, Vector e)]
                          -> Herm Matrix e
weightedCovMatrixWithMean mu t wxs = runHermMatrix $ do
    cov <- Herm uplo `fmap` M.new_ (p,p)
    weightedCovToMatrixWithMean mu t wxs cov
    return cov
  where
    p = V.dim mu
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the weighed sample covariance matrix
-- (in packed form) with storage scheme equal to 'defaultCovUplo'.
weightedCovPackedWithMean :: (BLAS2 e)
                          => Vector e -> CovMethod -> [(Double, Vector e)]
                          -> Herm Packed e
weightedCovPackedWithMean mu t wxs = P.hermCreate $ do
    cov <- (Herm uplo . P.unsafeViewFromSTVector p) `fmap` V.new_ (p*(p+1) `div` 2)
    weightedCovToPackedWithMean mu t wxs cov
    return cov
  where
    p = V.dim mu
    uplo = defaultCovUplo

-- | Computes and copies the sample covariance matrix to the given
-- destination.
covToMatrix :: (RVector v, BLAS3 e)
            => CovMethod -> [v e] -> Herm (STMatrix s) e -> ST s ()
covToMatrix t xs cov@(Herm _ a) = do
    mu <- V.new p 1
    meanToVector xs mu
    covToMatrixWithMean mu t xs cov
  where
    (p,_) = M.dim a

-- | Computes and copies the sample covariance matrix (in packed form)
-- to the given destination.
covToPacked :: (RVector v, BLAS2 e)
            => CovMethod -> [v e] -> Herm (STPacked s) e -> ST s ()
covToPacked t xs cov@(Herm _ a) = do
    mu <- V.new p 1
    meanToVector xs mu
    covToPackedWithMean mu t xs cov
  where
    p = P.dim a

-- | Given the pre-computed mean, computes and copies the sample covariance
-- matrix to the given destination.
covToMatrixWithMean :: (RVector v1, RVector v2, BLAS3 e)
                    => v1 e -> CovMethod -> [v2 e] -> Herm (STMatrix s) e
                    -> ST s ()
covToMatrixWithMean mu t xs cov@(Herm _ a)
    | M.dim a /= (p,p) = error $
        printf ("covToMatrixWithMean <vector with dim %d> _ _"
                ++ " <matrix with dim %s>: dimension mismatch")
               n (show $ M.dim a)
    | otherwise = do
        xt <- M.new_ (p,n)
        M.withColsST xt $ \xs' ->
            sequence_ [ V.subTo x mu x'
                      | (x,x') <- zip xs xs'
                      ]
        rankKUpdateToHermMatrix (1/df) NoTrans xt 0 cov
  where
    p = V.dim mu
    n = length xs
    df = fromIntegral $ case t of { MLCov -> n ; UnbiasedCov -> n - 1 }

-- | Given the pre-computed mean, computes and copies the sample covariance
-- matrix (in packed form) to the given destination.
covToPackedWithMean :: (RVector v1, RVector v2, BLAS2 e)
                    => v1 e -> CovMethod -> [v2 e] -> Herm (STPacked s) e
                    -> ST s ()
covToPackedWithMean mu t xs cov@(Herm _ a)
    | P.dim a /= p = error $
        printf ("covToPackedWithMean <vector with dim %d> _ _"
                ++ " (Herm _ <packed matrix with dim %d>):"
                ++ " dimension mismatch")
               n (P.dim a)
    | otherwise = do
        xt <- M.new_ (p,n)
        M.withColsST xt $ \xs' ->
            sequence_ [ V.subTo x mu x'
                      | (x,x') <- zip xs xs'
                      ]
        P.withSTVectorView a V.clear
        M.withColsST xt $ \xs' ->
            sequence_ [ P.hermRank1UpdateTo scale x' cov | x' <- xs' ]
  where
    p = V.dim mu
    n = length xs
    df = fromIntegral $ case t of { MLCov -> n ; UnbiasedCov -> n - 1 }
    scale = 1/df

-- | Computes and copies the weighed sample covariance matrix to the
-- given destination.
weightedCovToMatrix :: (RVector v, BLAS3 e)
                    => CovMethod -> [(Double, v e)] -> Herm (STMatrix s) e
                    -> ST s ()
weightedCovToMatrix t wxs cov@(Herm _ a) = do
    mu <- V.new p 1
    weightedMeanToVector wxs mu
    weightedCovToMatrixWithMean mu t wxs cov
  where
    (p,_) = M.dim a

-- | Computes and copies the weighed sample covariance matrix (in packed
-- form) to the given destination.
weightedCovToPacked :: (RVector v, BLAS2 e)
                    => CovMethod -> [(Double, v e)] -> Herm (STPacked s) e
                    -> ST s ()
weightedCovToPacked t wxs cov@(Herm _ a) = do
    mu <- V.new p 1
    weightedMeanToVector wxs mu
    weightedCovToPackedWithMean mu t wxs cov
  where
    p = P.dim a

-- | Given the pre-computed mean, computes and copies the weighed sample
-- covariance matrix to the given destination.
weightedCovToMatrixWithMean :: (RVector v1, RVector v2, BLAS3 e)
                            => v1 e -> CovMethod -> [(Double, v2 e)]
                            -> Herm (STMatrix s) e -> ST s ()
weightedCovToMatrixWithMean mu t wxs cov@(Herm _ a)
    | M.dim a /= (p,p) = error $
        printf ("weightedCovToMatrixWithMean <vector with dim %d> _ _"
                ++ " (Herm _ <matrix with dim %s>):"
                ++ " dimension mismatch")
               n (show $ M.dim a)
    | otherwise = do
        xt <- M.new_ (p,n)
        M.withColsST xt $ \xs' ->
            sequence_ [  V.subTo x mu x'
                      >> V.scaleTo (realToFrac $ sqrt (w / invscale)) x' x'
                      |  (w,x,x') <- zip3 ws xs xs'
                      ]
        rankKUpdateToHermMatrix 1 NoTrans xt 0 cov
  where
    (ws0,xs) = unzip wxs
    w_sum = foldl' (+) 0 ws0
    ws = if w_sum == 0 then ws0 else map (/w_sum) ws0
    w2s_sum = foldl' (+) 0 $ map (^^(2::Int)) ws
    invscale = case t of 
                   MLCov -> 1
                   UnbiasedCov -> (1 - w2s_sum)
    n = length ws0
    p = V.dim mu

-- | Given the pre-computed mean, computes and copies the weighed sample
-- covariance matrix (in packed form) to the given destination.
weightedCovToPackedWithMean :: (RVector v1, RVector v2, BLAS2 e)
                            => v1 e -> CovMethod -> [(Double, v2 e)]
                            -> Herm (STPacked s) e -> ST s ()
weightedCovToPackedWithMean mu t wxs cov@(Herm _ a)
    | P.dim a /= p = error $
        printf ("weightedCovToPackedWithMean <vector with dim %d> _ _"
                ++ " (Herm _ <packed matrix with dim %d>):"
                ++ " dimension mismatch")
               n (P.dim a)
    | otherwise = do
        xt <- M.new_ (p,n)
        M.withColsST xt $ \xs' ->
            sequence_ [  V.subTo x mu x'
                      >> V.scaleTo (realToFrac $ sqrt (w / invscale)) x' x'
                      |  (w,x,x') <- zip3 ws xs xs'
                      ]
        P.withSTVectorView a V.clear                      
        M.withCols xt $ \xs' ->
            sequence_ [ P.hermRank1UpdateTo 1 x' cov | x' <- xs' ]
  where
    (ws0,xs) = unzip wxs
    w_sum = foldl' (+) 0 ws0
    ws = if w_sum == 0 then ws0 else map (/w_sum) ws0
    w2s_sum = foldl' (+) 0 $ map (^^(2::Int)) ws
    invscale = case t of 
                   MLCov -> 1
                   UnbiasedCov -> (1 - w2s_sum)
    n = length ws0
    p = V.dim mu
