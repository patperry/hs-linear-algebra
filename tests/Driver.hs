{-# LANGUAGE CPP #-}
module Driver (
    E,
    
    Nat(..),
    Nat2(..),
    Pos(..),
    Pos2(..),
    Index(..),
    Index2(..),
    Assocs(..),
    Assocs2(..),
    BandedAssocs(..),
    
    module Control.Arrow,
    module Control.Monad,
    module Control.Monad.ST,
    
    AEq,
    (===),
    (~==),
    module Data.Function,
    module Data.Ix,
    module Data.List,
    module Data.Ord,
    
    module Debug.Trace,
    
    module Test.QuickCheck,
    
    module Text.Printf,

    field,

    module Test.Framework,
    module Test.Framework.Providers.QuickCheck2,
    ) where

import Control.Arrow
import Control.Monad
import Control.Monad.ST

import Data.AEq( AEq )
import qualified Data.AEq as AEq
import Data.Complex
import Data.Ix
import Data.Function
import Data.List
import Data.Ord

import Debug.Trace

import System.IO
import System.Random

import Test.Vector.Dense()
import Test.Matrix.Dense()
import Test.Matrix.Banded()

import Test.QuickCheck hiding ( vector )
import Test.QuickCheck.BLAS( Nat(..), Nat2(..), Pos(..), Pos2(..), 
    Index(..), Index2(..), Assocs(..), Assocs2(..), BandedAssocs(..) )
import qualified Test.QuickCheck.BLAS as Test

import Test.Framework
import Test.Framework.Providers.QuickCheck2

import Text.Printf
import Text.Show.Functions

#ifdef COMPLEX
field = "Complex Double"
type E = Complex Double 
#else
field = "Double"
type E = Double
#endif

infix 4 ===, ~==

x === y | (AEq.===) x y = True
        | otherwise = trace (printf "expected `%s', but got `%s'" (show y) (show x)) False

x ~== y | (AEq.~==) x y = True
        | otherwise = trace (printf "expected `%s', but got `%s'" (show y) (show x)) False
