-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Packed.Statistics
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Basic multivariate statistics.
--

module Numeric.LinearAlgebra.Packed.Statistics (
    defaultCovUplo,

    -- * Immutable interface
    cov,
    covWithMean,
    weightedCov,
    weightedCovWithMean,

    -- * Mutable interface
    covTo,
    covWithMeanTo,
    weightedCovTo,
    weightedCovWithMeanTo,

    ) where

import Control.Monad.ST( ST )
import Data.List( foldl' )
import Text.Printf( printf )

import Numeric.LinearAlgebra.Types
import Numeric.LinearAlgebra.Packed.Base( Packed, STPacked )
import qualified Numeric.LinearAlgebra.Packed.Base as P

import Numeric.LinearAlgebra.Vector( Vector, RVector )
import qualified Numeric.LinearAlgebra.Vector as V

import qualified Numeric.LinearAlgebra.Matrix as M


-- | Returns the default storage scheme for covariance matrices.
defaultCovUplo :: Uplo
defaultCovUplo = Lower

-- | Returns the sample covariance matrix hermitian matrix (in packed form)
-- with storage scheme equal to 'defaultCovUplo'.  The first argument gives
-- the dimension of the vectors.
cov :: (BLAS2 e)
    => Int -> CovMethod -> [Vector e] -> Herm Packed e
cov p t xs = P.hermCreate $ do
    c <- Herm uplo `fmap` P.new_ p
    covTo t xs c
    return c
  where
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the sample covariance matrix
-- (in packed form) with storage scheme equal to 'defaultCovUplo'.
covWithMean :: (BLAS2 e)
            => Vector e -> CovMethod -> [Vector e] -> Herm Packed e
covWithMean mu t xs = P.hermCreate $ do
    c <- Herm uplo `fmap` P.new_ p
    covWithMeanTo mu t xs c
    return c
  where
    p = V.dim mu
    uplo = defaultCovUplo

-- | Returns the weighed sample covariance matrix (in packed form) with
-- storage scheme equal to 'defaultCovUplo'. The first argument gives the
-- dimension of the vectors.
weightedCov :: (BLAS2 e)
            => Int -> CovMethod -> [(Double, Vector e)] -> Herm Packed e
weightedCov p t wxs = P.hermCreate $ do
    c <- Herm uplo `fmap` P.new_ p
    weightedCovTo t wxs c
    return c
  where
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the weighed sample covariance matrix
-- (in packed form) with storage scheme equal to 'defaultCovUplo'.
weightedCovWithMean :: (BLAS2 e)
                    => Vector e -> CovMethod -> [(Double, Vector e)]
                    -> Herm Packed e
weightedCovWithMean mu t wxs = P.hermCreate $ do
    c <- Herm uplo `fmap` P.new_ p
    weightedCovWithMeanTo mu t wxs c
    return c
  where
    p = V.dim mu
    uplo = defaultCovUplo

-- | Computes and copies the sample covariance matrix (in packed form)
-- to the given destination.
covTo :: (RVector v, BLAS2 e)
      => CovMethod -> [v e] -> Herm (STPacked s) e -> ST s ()
covTo t xs c@(Herm _ a) = do
    mu <- V.new p 1
    V.meanTo xs mu
    covWithMeanTo mu t xs c
  where
    p = P.dim a

-- | Given the pre-computed mean, computes and copies the sample covariance
-- matrix (in packed form) to the given destination.
covWithMeanTo :: (RVector v1, RVector v2, BLAS2 e)
              => v1 e -> CovMethod -> [v2 e] -> Herm (STPacked s) e
              -> ST s ()
covWithMeanTo mu t xs c@(Herm _ a)
    | P.dim a /= p = error $
        printf ("covWithMeanTo <vector with dim %d> _ _"
                ++ " (Herm _ <packed matrix with dim %d>):"
                ++ " dimension mismatch")
               n (P.dim a)
    | otherwise = do
        xt <- M.new_ (p,n)
        M.withSTColViews xt $ \xs' ->
            sequence_ [ V.subTo x' mu x
                      | (x,x') <- zip xs xs'
                      ]
        P.withSTVectorView a V.clear
        M.withSTColViews xt $ \xs' ->
            sequence_ [ P.hermRank1UpdateM_ scale x' c | x' <- xs' ]
  where
    p = V.dim mu
    n = length xs
    df = fromIntegral $ case t of { MLCov -> n ; UnbiasedCov -> n - 1 }
    scale = 1/df

-- | Computes and copies the weighed sample covariance matrix (in packed
-- form) to the given destination.
weightedCovTo :: (RVector v, BLAS2 e)
              => CovMethod -> [(Double, v e)] -> Herm (STPacked s) e
              -> ST s ()
weightedCovTo t wxs c@(Herm _ a) = do
    mu <- V.new p 1
    V.weightedMeanTo wxs mu
    weightedCovWithMeanTo mu t wxs c
  where
    p = P.dim a

-- | Given the pre-computed mean, computes and copies the weighed sample
-- covariance matrix (in packed form) to the given destination.
weightedCovWithMeanTo :: (RVector v1, RVector v2, BLAS2 e)
                      => v1 e -> CovMethod -> [(Double, v2 e)]
                      -> Herm (STPacked s) e -> ST s ()
weightedCovWithMeanTo mu t wxs c@(Herm _ a)
    | P.dim a /= p = error $
        printf ("weightedCovWithMeanTo <vector with dim %d> _ _"
                ++ " (Herm _ <packed matrix with dim %d>):"
                ++ " dimension mismatch")
               n (P.dim a)
    | otherwise = do
        xt <- M.new_ (p,n)
        M.withSTColViews xt $ \xs' ->
            sequence_ [  V.subTo x' mu x
                      >> V.scaleByM_ (realToFrac $ sqrt (w / invscale)) x'
                      |  (w,x,x') <- zip3 ws xs xs'
                      ]
        P.withSTVectorView a V.clear                      
        M.withColViews xt $ \xs' ->
            sequence_ [ P.hermRank1UpdateM_ 1 x' c | x' <- xs' ]
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
