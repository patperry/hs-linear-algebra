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
import Numeric.LinearAlgebra.Vector
import Numeric.LinearAlgebra.Vector.ST
import Numeric.LinearAlgebra.Matrix
import Numeric.LinearAlgebra.Matrix.Herm
import Numeric.LinearAlgebra.Matrix.Packed
import Numeric.LinearAlgebra.Matrix.ST


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
defaultCovUplo = Upper

-- | Returns the sum of the vectors.  The first argument gives the dimension
-- of the vectors.
sumVector :: (VNum e) => Int -> [Vector e] -> Vector e
sumVector p xs = runVector $ do
    s <- newVector_ p
    sumToVector xs s
    return s

-- | Returns the mean of the vectors.  The first argument gives the dimension
-- of the vectors.
meanVector :: (VNum e, Fractional e) => Int -> [Vector e] -> Vector e
meanVector p xs = runVector $ do
      m <- newVector_ p
      meanToVector xs m
      return m

-- | Returns the weighted sum of the vectors.  The first argument gives the
-- dimension of the vectors.
weightedSumVector :: (VNum e) => Int -> [(e, Vector e)] -> Vector e
weightedSumVector p wxs = runVector $ do
    s <- newVector_ p
    weightedSumToVector wxs s
    return s

-- | Returns the weighted mean of the vectors.  The first argument gives the
-- dimension of the vectors.
weightedMeanVector :: (VNum e, Fractional e) => Int -> [(Double, Vector e)] -> Vector e
weightedMeanVector p wxs = runVector $ do
       s <- newVector_ p
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
    err <- newVector n 0
    old_s <- newVector_ n
    diff <- newVector_ n
    val <- newVector_ n
    
    setElemsVector s (replicate n 0)
    forM_ wxs $ \(w,x) -> do
        unsafeCopyToVector s old_s -- old_s := s
        scaleToVector w x val      -- val := w * x
        addToVector err val err    -- err := err + val
        addToVector s err s        -- s := s + err
        
        subToVector old_s s diff   -- diff := old_s - s
        addToVector diff val err   -- err := diff + val
  where
    n = dimVector s

-- | Sets the target vector to the weighted mean of the vectors.
weightedMeanToVector :: (RVector v, VNum e, Fractional e)
                     => [(Double, v e)] -> STVector s e -> ST s ()
weightedMeanToVector wxs m = let
    go _ _ [] = return ()
    go diff w_sum ((w,x):wxs') | w == 0    = go diff w_sum wxs'
                               | otherwise = let w_sum' = w_sum + w
                                             in do
                                    subToVector x m diff
                                    addToVectorWithScales
                                        (realToFrac $ w/w_sum') diff 1 m m
                                    go diff w_sum' wxs'
    in do
        diff <- newVector_ n
        setElemsVector m (replicate n 0)
        go diff 0 wxs
  where
    n = dimVector m

-- | Returns the sample covariance matrix as a hermitian matrix with storage
-- scheme equal to 'defaultCovUplo'.  The first argument gives the dimension
-- of the vectors.
covMatrix :: (VNum e, BLAS3 e)
          => Int -> CovMethod -> [Vector e] -> Herm Matrix e
covMatrix p t xs = runHermMatrix $ do
    cov <- Herm uplo `fmap` newMatrix_ (p,p)
    covToMatrix t xs cov
    return cov
  where
    uplo = defaultCovUplo

-- | Returns the sample covariance matrix hermitian matrix (in packed form)
-- with storage scheme equal to 'defaultCovUplo'.  The first argument gives
-- the dimension of the vectors.
covPacked :: (VNum e, BLAS2 e)
          => Int -> CovMethod -> [Vector e] -> Herm Packed e
covPacked p t xs = runHermPacked $ do
    cov <- (Herm uplo . unsafeToSTPacked p) `fmap` newVector_ (p*(p+1) `div` 2)
    covToPacked t xs cov
    return cov
  where
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the sample covariance matrix
-- with storage scheme equal to 'defaultCovUplo'.
covMatrixWithMean :: (VNum e, BLAS3 e)
                  => Vector e -> CovMethod -> [Vector e] -> Herm Matrix e
covMatrixWithMean mu t xs = runHermMatrix $ do
    cov <- Herm uplo `fmap` newMatrix_ (p,p)
    covToMatrixWithMean mu t xs cov
    return cov
  where
    p = dimVector mu
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the sample covariance matrix
-- (in packed form) with storage scheme equal to 'defaultCovUplo'.
covPackedWithMean :: (VNum e, BLAS2 e)
                  => Vector e -> CovMethod -> [Vector e] -> Herm Packed e
covPackedWithMean mu t xs = runHermPacked $ do
    cov <- (Herm uplo . unsafeToSTPacked p) `fmap` newVector_ (p*(p+1) `div` 2)
    covToPackedWithMean mu t xs cov
    return cov
  where
    p = dimVector mu
    uplo = defaultCovUplo

-- | Returns the weighed sample covariance matrix with storage scheme equal
-- to 'defaultCovUplo'. The first argument gives the dimension of the vectors.
weightedCovMatrix :: (VFloating e, BLAS3 e)
                  => Int -> CovMethod -> [(Double, Vector e)] -> Herm Matrix e
weightedCovMatrix p t wxs = runHermMatrix $ do
    cov <- Herm uplo `fmap` newMatrix_ (p,p)
    weightedCovToMatrix t wxs cov
    return cov
  where
    uplo = defaultCovUplo

-- | Returns the weighed sample covariance matrix (in packed form) with
-- storage scheme equal to 'defaultCovUplo'. The first argument gives the
-- dimension of the vectors.
weightedCovPacked :: (VFloating e, BLAS2 e)
                  => Int -> CovMethod -> [(Double, Vector e)] -> Herm Packed e
weightedCovPacked p t wxs = runHermPacked $ do
    cov <- (Herm uplo . unsafeToSTPacked p) `fmap` newVector_ (p*(p+1) `div` 2)
    weightedCovToPacked t wxs cov
    return cov
  where
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the weighed sample covariance matrix
-- with storage scheme equal to 'defaultCovUplo'.
weightedCovMatrixWithMean :: (VFloating e, BLAS3 e)
                          => Vector e -> CovMethod -> [(Double, Vector e)]
                          -> Herm Matrix e
weightedCovMatrixWithMean mu t wxs = runHermMatrix $ do
    cov <- Herm uplo `fmap` newMatrix_ (p,p)
    weightedCovToMatrixWithMean mu t wxs cov
    return cov
  where
    p = dimVector mu
    uplo = defaultCovUplo

-- | Given the pre-computed mean, returns the weighed sample covariance matrix
-- (in packed form) with storage scheme equal to 'defaultCovUplo'.
weightedCovPackedWithMean :: (VFloating e, BLAS2 e)
                          => Vector e -> CovMethod -> [(Double, Vector e)]
                          -> Herm Packed e
weightedCovPackedWithMean mu t wxs = runHermPacked $ do
    cov <- (Herm uplo . unsafeToSTPacked p) `fmap` newVector_ (p*(p+1) `div` 2)
    weightedCovToPackedWithMean mu t wxs cov
    return cov
  where
    p = dimVector mu
    uplo = defaultCovUplo

-- | Computes and copies the sample covariance matrix to the given
-- destination.
covToMatrix :: (RVector v, VNum e, BLAS3 e)
            => CovMethod -> [v e] -> Herm (STMatrix s) e -> ST s ()
covToMatrix t xs cov@(Herm _ a) = do
    mu <- newVector p 1
    meanToVector xs mu
    covToMatrixWithMean mu t xs cov
  where
    (p,_) = dimMatrix a

-- | Computes and copies the sample covariance matrix (in packed form)
-- to the given destination.
covToPacked :: (RVector v, VNum e, BLAS2 e)
            => CovMethod -> [v e] -> Herm (STPacked s) e -> ST s ()
covToPacked t xs cov@(Herm _ a) = do
    mu <- newVector p 1
    meanToVector xs mu
    covToPackedWithMean mu t xs cov
  where
    p = dimPacked a

-- | Given the pre-computed mean, computes and copies the sample covariance
-- matrix to the given destination.
covToMatrixWithMean :: (RVector v1, RVector v2, VNum e, BLAS3 e)
                    => v1 e -> CovMethod -> [v2 e] -> Herm (STMatrix s) e
                    -> ST s ()
covToMatrixWithMean mu t xs cov@(Herm _ a)
    | dimMatrix a /= (p,p) = error $
        printf ("covToMatrixWithMean <vector with dim %d> _ _"
                ++ " <matrix with dim %s>: dimension mismatch")
               n (show $ dimMatrix a)
    | otherwise = do
        xt <- newMatrix_ (p,n)
        withColsMatrixST xt $ \xs' ->
            sequence_ [ subToVector x mu x'
                      | (x,x') <- zip xs xs'
                      ]
        rankKUpdateToHermMatrix (1/df) NoTrans xt 0 cov
  where
    p = dimVector mu
    n = length xs
    df = fromIntegral $ case t of { MLCov -> n ; UnbiasedCov -> n - 1 }

-- | Given the pre-computed mean, computes and copies the sample covariance
-- matrix (in packed form) to the given destination.
covToPackedWithMean :: (RVector v1, RVector v2, VNum e, BLAS2 e)
                    => v1 e -> CovMethod -> [v2 e] -> Herm (STPacked s) e
                    -> ST s ()
covToPackedWithMean mu t xs cov@(Herm _ a)
    | dimPacked a /= p = error $
        printf ("covToPackedWithMean <vector with dim %d> _ _"
                ++ " (Herm _ <packed matrix with dim %d>):"
                ++ " dimension mismatch")
               n (dimPacked a)
    | otherwise = do
        xt <- newMatrix_ (p,n)
        withColsMatrixST xt $ \xs' ->
            sequence_ [ subToVector x mu x'
                      | (x,x') <- zip xs xs'
                      ]
        withVectorPackedST a clearVector
        withColsMatrixST xt $ \xs' ->
            sequence_ [ rank1UpdateToHermPacked scale x' cov | x' <- xs' ]
  where
    p = dimVector mu
    n = length xs
    df = fromIntegral $ case t of { MLCov -> n ; UnbiasedCov -> n - 1 }
    scale = 1/df

-- | Computes and copies the weighed sample covariance matrix to the
-- given destination.
weightedCovToMatrix :: (RVector v, VFloating e, BLAS3 e)
                    => CovMethod -> [(Double, v e)] -> Herm (STMatrix s) e
                    -> ST s ()
weightedCovToMatrix t wxs cov@(Herm _ a) = do
    mu <- newVector p 1
    weightedMeanToVector wxs mu
    weightedCovToMatrixWithMean mu t wxs cov
  where
    (p,_) = dimMatrix a

-- | Computes and copies the weighed sample covariance matrix (in packed
-- form) to the given destination.
weightedCovToPacked :: (RVector v, VFloating e, BLAS2 e)
                    => CovMethod -> [(Double, v e)] -> Herm (STPacked s) e
                    -> ST s ()
weightedCovToPacked t wxs cov@(Herm _ a) = do
    mu <- newVector p 1
    weightedMeanToVector wxs mu
    weightedCovToPackedWithMean mu t wxs cov
  where
    p = dimPacked a

-- | Given the pre-computed mean, computes and copies the weighed sample
-- covariance matrix to the given destination.
weightedCovToMatrixWithMean :: (RVector v1, RVector v2, VFloating e, BLAS3 e)
                            => v1 e -> CovMethod -> [(Double, v2 e)]
                            -> Herm (STMatrix s) e -> ST s ()
weightedCovToMatrixWithMean mu t wxs cov@(Herm _ a)
    | dimMatrix a /= (p,p) = error $
        printf ("weightedCovToMatrixWithMean <vector with dim %d> _ _"
                ++ " (Herm _ <matrix with dim %s>):"
                ++ " dimension mismatch")
               n (show $ dimMatrix a)
    | otherwise = do
        xt <- newMatrix_ (p,n)
        withColsMatrixST xt $ \xs' ->
            sequence_ [  subToVector x mu x'
                      >> scaleToVector (realToFrac $ sqrt (w / invscale)) x' x'
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
    p = dimVector mu

-- | Given the pre-computed mean, computes and copies the weighed sample
-- covariance matrix (in packed form) to the given destination.
weightedCovToPackedWithMean :: (RVector v1, RVector v2, VFloating e, BLAS2 e)
                            => v1 e -> CovMethod -> [(Double, v2 e)]
                            -> Herm (STPacked s) e -> ST s ()
weightedCovToPackedWithMean mu t wxs cov@(Herm _ a)
    | dimPacked a /= p = error $
        printf ("weightedCovToPackedWithMean <vector with dim %d> _ _"
                ++ " (Herm _ <packed matrix with dim %d>):"
                ++ " dimension mismatch")
               n (dimPacked a)
    | otherwise = do
        xt <- newMatrix_ (p,n)
        withColsMatrixST xt $ \xs' ->
            sequence_ [  subToVector x mu x'
                      >> scaleToVector (realToFrac $ sqrt (w / invscale)) x' x'
                      |  (w,x,x') <- zip3 ws xs xs'
                      ]
        withVectorPackedST a clearVector                      
        withColsMatrix xt $ \xs' ->
            sequence_ [ rank1UpdateToHermPacked 1 x' cov | x' <- xs' ]
  where
    (ws0,xs) = unzip wxs
    w_sum = foldl' (+) 0 ws0
    ws = if w_sum == 0 then ws0 else map (/w_sum) ws0
    w2s_sum = foldl' (+) 0 $ map (^^(2::Int)) ws
    invscale = case t of 
                   MLCov -> 1
                   UnbiasedCov -> (1 - w2s_sum)
    n = length ws0
    p = dimVector mu
