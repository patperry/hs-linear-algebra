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

import Control.Monad( when )
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
    covTo c t xs
    return c
  where
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the sample covariance matrix
-- (in packed form) with storage scheme equal to 'defaultCovUplo'.
covWithMean :: (BLAS2 e)
            => Vector e -> CovMethod -> [Vector e] -> Herm Packed e
covWithMean mu t xs = P.hermCreate $ do
    c <- Herm uplo `fmap` P.new_ p
    covWithMeanTo c mu t xs
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
    weightedCovTo c t wxs
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
    weightedCovWithMeanTo c mu t wxs
    return c
  where
    p = V.dim mu
    uplo = defaultCovUplo

-- | Computes and copies the sample covariance matrix (in packed form)
-- to the given destination.
covTo :: (RVector v, BLAS2 e)
      => Herm (STPacked s) e -> CovMethod -> [v e] -> ST s ()
covTo c@(Herm _ a) t xs = do
    p <- P.getDim a
    mu <- V.new p 1
    V.meanTo mu xs
    covWithMeanTo c mu t xs


-- | Given the pre-computed mean, computes and copies the sample covariance
-- matrix (in packed form) to the given destination.
covWithMeanTo :: (RVector v1, RVector v2, BLAS2 e)
              => Herm (STPacked s) e
              -> v1 e -> CovMethod -> [v2 e]
              -> ST s ()
covWithMeanTo c@(Herm _ a) mu t xs = do
    pa <- P.getDim a
    p <- V.getDim mu   
     
    when (pa /= p) $ error $
        printf ("covWithMeanTo"
                ++ " (Herm _ <packed matrix with dim %d>)"
                ++ " <vector with dim %d>"
                ++ " _ _"
                ++ ": dimension mismatch")
               pa p

    xt <- M.new_ (p,n)
    M.withColsM xt $ \xs' ->
        sequence_ [ V.subTo x' mu x
                  | (x,x') <- zip xs xs'
                  ]
    P.withVectorM a V.clear
    M.withColsM xt $ \xs' ->
        sequence_ [ P.hermRank1UpdateM_ scale x' c | x' <- xs' ]
  where
    n = length xs
    df = fromIntegral $ case t of { MLCov -> n ; UnbiasedCov -> n - 1 }
    scale = 1/df


-- | Computes and copies the weighed sample covariance matrix (in packed
-- form) to the given destination.
weightedCovTo :: (RVector v, BLAS2 e)
              => Herm (STPacked s) e
              -> CovMethod -> [(Double, v e)] 
              -> ST s ()
weightedCovTo c@(Herm _ a) t wxs = do
    p <- P.getDim a
    mu <- V.new p 1
    V.weightedMeanTo mu wxs
    weightedCovWithMeanTo c mu t wxs


-- | Given the pre-computed mean, computes and copies the weighed sample
-- covariance matrix (in packed form) to the given destination.
weightedCovWithMeanTo :: (RVector v1, RVector v2, BLAS2 e)
                      => Herm (STPacked s) e
                      -> v1 e -> CovMethod -> [(Double, v2 e)]
                      -> ST s ()
weightedCovWithMeanTo c@(Herm _ a) mu t wxs = do
    pa <- P.getDim a
    p <- V.getDim mu
    
    when (pa /= p) $ error $
        printf ("weightedCovWithMeanTo"
                ++ " (Herm _ <packed matrix with dim %d>)"
                ++ " <vector with dim %d>"
                ++ " _ _"
                ++ ": dimension mismatch")
               pa p

    xt <- M.new_ (p,n)
    M.withColsM xt $ \xs' ->
        sequence_ [  V.subTo x' mu x
                  >> V.scaleM_ (realToFrac $ sqrt (w / invscale)) x'
                  |  (w,x,x') <- zip3 ws xs xs'
                  ]
    P.withVectorM a V.clear                      
    M.withCols xt $ \xs' ->
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
