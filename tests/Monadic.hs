{-# LANGUAGE RankNTypes #-}
-- | Allows testing of monadic values.  Taken from QuickCheck2
module Monadic( 
    PropertyM,
    assert,
    pre,
    run,
    pick,
    wp,
    forAllM,
    monitor,
    monadicIO,
    monadicST,
    ) where

import Control.Monad( liftM )
import Control.Monad.ST
import System.IO.Unsafe( unsafePerformIO )
import Test.QuickCheck

instance (Testable prop, Show prop) => Testable (Gen prop) where
    property mp = forAll mp property

instance Show Property where
    show = const "<property>"

--------------------------------------------------------------------------
-- type PropertyM

newtype PropertyM m a =
  MkPropertyM { unPropertyM :: (a -> Gen (m Property)) -> Gen (m Property) }

instance Functor (PropertyM m) where
  fmap f (MkPropertyM m) = MkPropertyM (\k -> m (k . f))

instance Monad m => Monad (PropertyM m) where
  return x            = MkPropertyM (\k -> k x)
  MkPropertyM m >>= f = MkPropertyM (\k -> m (\a -> unPropertyM (f a) k))
  fail s              = MkPropertyM (\k -> return (return (property result)))
   where
    result = Result (Just False) [s] []

-- should think about strictness/exceptions here
--assert :: Testable prop => prop -> PropertyM m ()
assert :: Monad m => Bool -> PropertyM m ()
assert b = MkPropertyM $ \k ->
  if b
    then k ()
    else return (return (property False))

-- should think about strictness/exceptions here
pre :: Monad m => Bool -> PropertyM m ()
pre b = MkPropertyM $ \k ->
  if b
    then k ()
    else return (return (property ()))

-- should be called lift?
run :: Monad m => m a -> PropertyM m a
run m = MkPropertyM (liftM (m >>=) . promote)

pick :: (Monad m, Show a) => Gen a -> PropertyM m a
pick gen = MkPropertyM $ \k ->
  do a <- gen
     mp <- k a
     return (do p <- mp
                return (forAll (return a) (const p)))

wp :: Monad m => m a -> (a -> PropertyM m b) -> PropertyM m b
wp m k = run m >>= k

forAllM :: (Monad m, Show a) => Gen a -> (a -> PropertyM m b) -> PropertyM m b
forAllM gen k = pick gen >>= k

monitor :: Monad m => (Property -> Property) -> PropertyM m ()
monitor f = MkPropertyM (\k -> (f `liftM`) `fmap` (k ()))

-- run functions

-- Can't make this work in any other way... :-(
monadicIO :: PropertyM IO a -> Property
monadicIO (MkPropertyM m) =
  property $
    unsafePerformIO `fmap`
      m (const (return (return (property True))))

unsafePropertySTToPropertyIO :: (forall s. PropertyM (ST s) a) -> PropertyM IO a
unsafePropertySTToPropertyIO pm = MkPropertyM $ \f ->
   let m  = unPropertyM pm
       f' = \a -> liftM unsafeIOToST $ f a
   in liftM unsafeSTToIO $ m f'

monadicST :: (forall s. PropertyM (ST s) a) -> Property
monadicST = monadicIO . unsafePropertySTToPropertyIO
