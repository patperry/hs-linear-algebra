{-# LANGUAGE RankNTypes #-}
-- | Allows testing of monadic values.
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

import Test.QuickCheck.Monadic
import Test.QuickCheck.Property

unsafePropertySTToPropertyIO :: (forall s. PropertyM (ST s) a) -> PropertyM IO a
unsafePropertySTToPropertyIO pm = MkPropertyM $ \f ->
   let m  = unPropertyM pm
       f' = \a -> liftM unsafeIOToST $ f a
   in liftM unsafeSTToIO $ m f'

monadicST :: (forall s. PropertyM (ST s) a) -> Property
monadicST = monadicIO . unsafePropertySTToPropertyIO
