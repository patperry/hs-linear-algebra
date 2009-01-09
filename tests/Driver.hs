{-# LANGUAGE CPP #-}
module Driver (
    E,
    
    Natural(..),
    Index(..),
    Assocs(..),
    
    Natural2(..),
    Index2(..),
    Assocs2(..),
    
    module Control.Arrow,
    module Control.Monad,
    module Control.Monad.ST,
    
    module Data.AEq,
    module Data.Function,
    module Data.Ix,
    module Data.List,
    module Data.Ord,
    
    module Debug.Trace,
    
    module Test.QuickCheck,
    
    module Text.Printf,

    field,

    mytest,
    mycheck,
    mytests,
    done,

    ) where

import Control.Arrow
import Control.Monad
import Control.Monad.ST

import Data.AEq
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
import qualified Test.QuickCheck.BLAS as Test

import Text.Printf
import Text.Show.Functions

#ifdef COMPLEX
field = "Complex Double"
type E = Complex Double 
#else
field = "Double"
type E = Double
#endif

newtype Natural = Nat Int deriving (Eq,Show)
instance Arbitrary Natural where
    arbitrary = do
        n <- arbitrary
        return $ Nat (abs n)
        
    coarbitrary = undefined

newtype Natural2 = Nat2 (Int,Int) deriving (Eq,Show)
instance Arbitrary Natural2 where
    arbitrary = liftM Nat2 Test.shape
    coarbitrary = undefined

data Index = Index Int Int deriving (Eq,Show)
instance Arbitrary Index where
    arbitrary = do
        n <- arbitrary
        i <- choose (0, abs n)
        return $ Index i (abs n + 1)
        
    coarbitrary = undefined

data Index2 = Index2 (Int,Int) (Int,Int) deriving (Eq,Show)
instance Arbitrary Index2 where
    arbitrary = do
        (m',n') <- Test.shape
        i <- choose (0,m')
        j <- choose (0,n')
        return $ Index2 (i,j) (m'+1,n'+1)
    
    coarbitrary = undefined

data Assocs = Assocs Int [(Int,E)] deriving (Eq,Show)
instance Arbitrary Assocs where
    arbitrary = do
        (Nat n) <- arbitrary
        (Nat s) <- if n == 0 then return (Nat 0) else arbitrary
        ies <- replicateM s $ liftM2 (,) (choose (0,n-1)) arbitrary
        return $ Assocs n ies
        
    coarbitrary = undefined

data Assocs2 = Assocs2 (Int,Int) [((Int,Int),E)] deriving (Eq,Show)
instance Arbitrary Assocs2 where
    arbitrary = do
        (Nat2 (m,n)) <- arbitrary
        (Nat s) <- if m == 0 || n == 0 then return (Nat 0) else arbitrary
        ijes <- replicateM s $ liftM2 (,) (liftM2 (,) (choose (0,m-1)) 
                                                      (choose (0,n-1)))
                                          arbitrary
        return $ Assocs2 (m,n) ijes
        
    coarbitrary = undefined


------------------------------------------------------------------------
--
-- QC driver ( taken from xmonad-0.6 )
--

debug = False

mytest :: Testable a => a -> Int -> IO (Bool, Int)
mytest a n = mycheck defaultConfig
    { configMaxTest=n
    , configEvery   = \n args -> let s = show n in s ++ [ '\b' | _ <- s ] } a
 -- , configEvery= \n args -> if debug then show n ++ ":\n" ++ unlines args else [] } a

mycheck :: Testable a => Config -> a -> IO (Bool, Int)
mycheck config a = do
    rnd <- newStdGen
    mytests config (evaluate a) rnd 0 0 []

mytests :: Config -> Gen Result -> StdGen -> Int -> Int -> [[String]] -> IO (Bool, Int)
mytests config gen rnd0 ntest nfail stamps
    | ntest == configMaxTest config = done "OK," ntest stamps >> return (True, ntest)
    | nfail == configMaxFail config = done "Arguments exhausted after" ntest stamps >> return (True, ntest)
    | otherwise               =
      do putStr (configEvery config ntest (arguments result)) >> hFlush stdout
         case ok result of
           Nothing    ->
             mytests config gen rnd1 ntest (nfail+1) stamps
           Just True  ->
             mytests config gen rnd1 (ntest+1) nfail (stamp result:stamps)
           Just False ->
             putStr ( "Falsifiable after "
                   ++ show ntest
                   ++ " tests:\n"
                   ++ unlines (arguments result)
                    ) >> hFlush stdout >> return (False, ntest)
     where
      result      = generate (configSize config ntest) rnd2 gen
      (rnd1,rnd2) = split rnd0

done :: String -> Int -> [[String]] -> IO ()
done mesg ntest stamps = putStr ( mesg ++ " " ++ show ntest ++ " tests" ++ table )
  where
    table = display
            . map entry
            . reverse
            . sort
            . map pairLength
            . group
            . sort
            . filter (not . null)
            $ stamps

    display []  = ".\n"
    display [x] = " (" ++ x ++ ").\n"
    display xs  = ".\n" ++ unlines (map (++ ".") xs)

    pairLength xss@(xs:_) = (length xss, xs)
    entry (n, xs)         = percentage n ntest
                       ++ " "
                       ++ concat (intersperse ", " xs)

    percentage n m        = show ((100 * n) `div` m) ++ "%"

------------------------------------------------------------------------

