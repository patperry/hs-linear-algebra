-----------------------------------------------------------------------------
-- |
-- Module     : Numeric.LinearAlgebra.Matrix.Statistics
-- Copyright  : Copyright (c) 2010, Patrick Perry <patperry@gmail.com>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@gmail.com>
-- Stability  : experimental
--
-- Basic multivariate statistics.
--

module Numeric.LinearAlgebra.Matrix.Statistics (
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

import Numeric.LinearAlgebra.Vector( Vector, RVector )
import qualified Numeric.LinearAlgebra.Vector as V

import Numeric.LinearAlgebra.Matrix.Base( Matrix )
import Numeric.LinearAlgebra.Matrix.STBase( STMatrix )
import qualified Numeric.LinearAlgebra.Matrix.STBase as M
import qualified Numeric.LinearAlgebra.Matrix.Herm as M


-- | Returns the default storage scheme for covariance matrices.
defaultCovUplo :: Uplo
defaultCovUplo = Lower


-- | Returns the sample covariance matrix as a hermitian matrix with storage
-- scheme equal to 'defaultCovUplo'.  The first argument gives the dimension
-- of the vectors.
cov :: (BLAS3 e)
          => Int -> CovMethod -> [Vector e] -> Herm Matrix e
cov p t xs = M.hermCreate $ do
    c <- Herm uplo `fmap` M.new_ (p,p)
    covTo t xs c
    return c
  where
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the sample covariance matrix
-- with storage scheme equal to 'defaultCovUplo'.
covWithMean :: (BLAS3 e)
                  => Vector e -> CovMethod -> [Vector e] -> Herm Matrix e
covWithMean mu t xs = M.hermCreate $ do
    c <- Herm uplo `fmap` M.new_ (p,p)
    covWithMeanTo mu t xs c
    return c
  where
    p = V.dim mu
    uplo = defaultCovUplo

-- | Returns the weighed sample covariance matrix with storage scheme equal
-- to 'defaultCovUplo'. The first argument gives the dimension of the vectors.
weightedCov :: (BLAS3 e)
                  => Int -> CovMethod -> [(Double, Vector e)] -> Herm Matrix e
weightedCov p t wxs = M.hermCreate $ do
    c <- Herm uplo `fmap` M.new_ (p,p)
    weightedCovTo t wxs c
    return c
  where
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the weighed sample covariance matrix
-- with storage scheme equal to 'defaultCovUplo'.
weightedCovWithMean :: (BLAS3 e)
                          => Vector e -> CovMethod -> [(Double, Vector e)]
                          -> Herm Matrix e
weightedCovWithMean mu t wxs = M.hermCreate $ do
    c <- Herm uplo `fmap` M.new_ (p,p)
    weightedCovWithMeanTo mu t wxs c
    return c
  where
    p = V.dim mu
    uplo = defaultCovUplo

-- | Computes and copies the sample covariance matrix to the given
-- destination.
covTo :: (RVector v, BLAS3 e)
            => CovMethod -> [v e] -> Herm (STMatrix s) e -> ST s ()
covTo t xs c@(Herm _ a) = do
    mu <- V.new p 1
    V.meanTo xs mu
    covWithMeanTo mu t xs c
  where
    (p,_) = M.dim a

-- | Given the pre-computed mean, computes and copies the sample covariance
-- matrix to the given destination.
covWithMeanTo :: (RVector v1, RVector v2, BLAS3 e)
                    => v1 e -> CovMethod -> [v2 e] -> Herm (STMatrix s) e
                    -> ST s ()
covWithMeanTo mu t xs c@(Herm _ a)
    | M.dim a /= (p,p) = error $
        printf ("covWithMeanTo <vector with dim %d> _ _"
                ++ " <matrix with dim %s>: dimension mismatch")
               n (show $ M.dim a)
    | otherwise = do
        xt <- M.new_ (p,n)
        M.withSTColViews xt $ \xs' ->
            sequence_ [ V.subTo x' mu x
                      | (x,x') <- zip xs xs'
                      ]
        M.hermRankKUpdateTo (1/df) NoTrans xt 0 c
  where
    p = V.dim mu
    n = length xs
    df = fromIntegral $ case t of { MLCov -> n ; UnbiasedCov -> n - 1 }

-- | Computes and copies the weighed sample covariance matrix to the
-- given destination.
weightedCovTo :: (RVector v, BLAS3 e)
                    => CovMethod -> [(Double, v e)] -> Herm (STMatrix s) e
                    -> ST s ()
weightedCovTo t wxs c@(Herm _ a) = do
    mu <- V.new p 1
    V.weightedMeanTo wxs mu
    weightedCovWithMeanTo mu t wxs c
  where
    (p,_) = M.dim a

-- | Given the pre-computed mean, computes and copies the weighed sample
-- covariance matrix to the given destination.
weightedCovWithMeanTo :: (RVector v1, RVector v2, BLAS3 e)
                            => v1 e -> CovMethod -> [(Double, v2 e)]
                            -> Herm (STMatrix s) e -> ST s ()
weightedCovWithMeanTo mu t wxs c@(Herm _ a)
    | M.dim a /= (p,p) = error $
        printf ("weightedCovWithMeanTo <vector with dim %d> _ _"
                ++ " (Herm _ <matrix with dim %s>):"
                ++ " dimension mismatch")
               n (show $ M.dim a)
    | otherwise = do
        xt <- M.new_ (p,n)
        M.withSTColViews xt $ \xs' ->
            sequence_ [  V.subTo x' mu x
                      >> V.scaleTo x' (realToFrac $ sqrt (w / invscale)) x'
                      |  (w,x,x') <- zip3 ws xs xs'
                      ]
        M.hermRankKUpdateTo 1 NoTrans xt 0 c
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
