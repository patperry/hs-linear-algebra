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

import Control.Monad( when )
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
    covTo c t xs
    return c
  where
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the sample covariance matrix
-- with storage scheme equal to 'defaultCovUplo'.
covWithMean :: (BLAS3 e)
                  => Vector e -> CovMethod -> [Vector e] -> Herm Matrix e
covWithMean mu t xs = M.hermCreate $ do
    p <- V.getDim mu    
    c <- Herm uplo `fmap` M.new_ (p,p)
    covWithMeanTo c mu t xs
    return c
  where
    uplo = defaultCovUplo

-- | Returns the weighed sample covariance matrix with storage scheme equal
-- to 'defaultCovUplo'. The first argument gives the dimension of the vectors.
weightedCov :: (BLAS3 e)
                  => Int -> CovMethod -> [(Double, Vector e)] -> Herm Matrix e
weightedCov p t wxs = M.hermCreate $ do
    c <- Herm uplo `fmap` M.new_ (p,p)
    weightedCovTo c t wxs
    return c
  where
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the weighed sample covariance matrix
-- with storage scheme equal to 'defaultCovUplo'.
weightedCovWithMean :: (BLAS3 e)
                    => Vector e -> CovMethod -> [(Double, Vector e)]
                    -> Herm Matrix e
weightedCovWithMean mu t wxs = M.hermCreate $ do
    p <- V.getDim mu    
    c <- Herm uplo `fmap` M.new_ (p,p)
    weightedCovWithMeanTo c mu t wxs
    return c
  where
    uplo = defaultCovUplo

-- | Computes and copies the sample covariance matrix to the given
-- destination.
covTo :: (RVector v, BLAS3 e)
      => Herm (STMatrix s) e -> CovMethod -> [v e] -> ST s ()
covTo c@(Herm _ a) t xs = do
    (p,_) <- M.getDim a
    mu <- V.new p 1
    V.meanTo mu xs
    covWithMeanTo c mu t xs
    

-- | Given the pre-computed mean, computes and copies the sample covariance
-- matrix to the given destination.
covWithMeanTo :: (RVector v1, RVector v2, BLAS3 e)
              => Herm (STMatrix s) e -> v1 e -> CovMethod -> [v2 e] -> ST s ()
covWithMeanTo c@(Herm _ a) mu t xs = do
    (ma,na) <- M.getDim a
    p <- V.getDim mu
    
    when ((ma,na) /= (p,p)) $ error $
        printf ("covWithMeanTo"
                ++ " (Herm _ <matrix with dim (%d,%d)>)"
                ++ " <vector with dim %d>"
                ++ " _ _"
                ++ ": dimension mismatch")
               ma na p

    xt <- M.new_ (p,n)
    M.withSTColViews xt $ \xs' ->
        sequence_ [ V.subTo x' mu x
                  | (x,x') <- zip xs xs'
                  ]
    M.hermRankKUpdateM_ (1/df) NoTrans xt 0 c
  where

    n = length xs
    df = fromIntegral $ case t of { MLCov -> n ; UnbiasedCov -> n - 1 }

-- | Computes and copies the weighed sample covariance matrix to the
-- given destination.
weightedCovTo :: (RVector v, BLAS3 e)
              => Herm (STMatrix s) e -> CovMethod -> [(Double, v e)] -> ST s ()
weightedCovTo c@(Herm _ a) t wxs = do
    (p,_) <- M.getDim a
    mu <- V.new p 1
    V.weightedMeanTo mu wxs
    weightedCovWithMeanTo c mu t wxs
    

-- | Given the pre-computed mean, computes and copies the weighed sample
-- covariance matrix to the given destination.
weightedCovWithMeanTo :: (RVector v1, RVector v2, BLAS3 e)
                      => Herm (STMatrix s) e
                      -> v1 e -> CovMethod -> [(Double, v2 e)]
                      -> ST s ()
weightedCovWithMeanTo c@(Herm _ a) mu t wxs = do
    (ma,na) <- M.getDim a
    p <- V.getDim mu

    when ((ma,na) /= (p,p)) $ error $
        printf ("weightedCovWithMeanTo"
                ++ " (Herm _ <matrix with dim (%d,%d)>):"
                ++ " <vector with dim %d>"
                ++ " _ _"
                ++ " dimension mismatch")
               ma na p

    xt <- M.new_ (p,n)
    M.withSTColViews xt $ \xs' ->
        sequence_ [  V.subTo x' mu x
                  >> V.scaleByM_ (realToFrac $ sqrt (w / invscale)) x'
                  |  (w,x,x') <- zip3 ws xs xs'
                  ]
    M.hermRankKUpdateM_ 1 NoTrans xt 0 c
  where
    (ws0,xs) = unzip wxs
    w_sum = foldl' (+) 0 ws0
    ws = if w_sum == 0 then ws0 else map (/w_sum) ws0
    w2s_sum = foldl' (+) 0 $ map (^^(2::Int)) ws
    invscale = case t of 
                   MLCov -> 1
                   UnbiasedCov -> (1 - w2s_sum)
    n = length ws0
