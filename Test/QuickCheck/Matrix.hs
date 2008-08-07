-----------------------------------------------------------------------------
-- |
-- Module     : Test.QuickCheck.Matrix
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Test.QuickCheck.Matrix
    where
                
import Data.List ( nub )
import Test.QuickCheck.Vector ( Index(..) )

import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC

newtype IndexPair = IndexPair (Int,Int) deriving (Eq, Show)
instance Arbitrary IndexPair where
    arbitrary = sized $ \k ->
        let k' = floor (sqrt $ fromInteger $ toInteger k :: Double)
        in resize k' $
            two arbitrary >>= (\((Index i), (Index j)) -> return $ IndexPair (i,j))
    coarbitrary (IndexPair ij) = coarbitrary ij
    
data Assocs e = Assocs (Int,Int) [((Int,Int),e)] deriving (Eq, Show)
instance Arbitrary e => Arbitrary (Assocs e) where
    arbitrary = sized $ \k -> 
        let k' = floor (sqrt $ fromInteger $ toInteger k :: Double)
        in do
            m <- choose (0,k')
            n <- choose (0,k')
            ijs <- QC.vector (m*n)
            let ijs' = filter (\(i,j) -> i < m && j < n) 
                           $ nub 
                           $ map (\(IndexPair ij) -> ij) ijs
            es <- QC.vector (length ijs')
            return $ Assocs (m,n) (zip ijs' es)
        
    coarbitrary (Assocs mn ijes) = coarbitrary (mn,ijes)

matrixSized :: (Int -> Gen a) -> Gen a
matrixSized f = 
    sized $ \s ->
        let s' = ceiling (sqrt $ fromInteger $ toInteger s :: Double)
        in f s'
        