{-# OPTIONS_HADDOCK hide #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Test.Matrix.Herm.Dense
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.Matrix.Herm.Dense( 
    HermMatrix(..),
    HermMatrixMV(..),
    HermMatrixMM(..),
    ) where

import Control.Monad( liftM )

import Test.QuickCheck hiding ( vector )
import Test.QuickCheck.BLAS ( TestElem )
import qualified Test.QuickCheck.BLAS as Test

import Data.Elem.BLAS
import Data.Vector.Dense hiding ( vector )
import Data.Matrix.Dense hiding ( matrix )
import Data.Matrix.Dense.ST
import Data.Matrix.Herm


hermRawMatrix :: (TestElem e) => Int -> Gen (Matrix e)
hermRawMatrix n = do
    a <- Test.matrix (n,n)
    d <- Test.realElems n
    return $ runSTMatrix $ do
        h <- unsafeThawMatrix a
        sequence_ [ writeElem h (i,j) (conjugate $ a!(j,i)) 
                  | i <- [ 0..n-1 ], j <- [ 0..i-1 ] ]
        setElems (diagView h 0) d
        return h

data HermMatrix e = 
    HermMatrix (Herm Matrix e)
               (Matrix e)
    deriving Show

instance (TestElem e) => Arbitrary (HermMatrix e) where
    arbitrary = do
        n <- liftM fst Test.shape
        a <- hermRawMatrix n
        
        junk <- Test.elems (n*n)
        let (u ,b ) = (Upper, a // zip (filter (uncurry (>)) $ indices a) junk)
            (u',b') = (Lower, a // zip (filter (uncurry (<)) $ indices a) junk)

        h <- elements [ hermFromBase u             b
                      , hermFromBase (flipUpLo u)  (herm b)
                      , hermFromBase u'            b'
                      , hermFromBase (flipUpLo u') (herm b')
                      ]
        return $ HermMatrix h a


data HermMatrixMV e = 
    HermMatrixMV (Herm Matrix e) 
                 (Matrix e) 
                 (Vector e) 
    deriving Show

instance (TestElem e) => Arbitrary (HermMatrixMV e) where
    arbitrary = do
        (HermMatrix h a) <- arbitrary
        x <- Test.vector (numCols a)
        return $ HermMatrixMV h a x

    
data HermMatrixMM e = 
    HermMatrixMM (Herm Matrix e) 
                 (Matrix e) 
                 (Matrix e) 
    deriving Show
    
instance (TestElem e) => Arbitrary (HermMatrixMM e) where
    arbitrary = do
        (HermMatrix h a) <- arbitrary
        n <- liftM fst Test.shape
        b <- Test.matrix (numCols a,n)
        return $ HermMatrixMM h a b
        