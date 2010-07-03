module STVector (
    tests_STVector
    ) where

import Control.Monad
import Control.Monad.ST
import Data.AEq
import Debug.Trace
import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck hiding ( vector )

import BLAS.Elem
import BLAS.Vector
import BLAS.Vector.ST

import Test.QuickCheck.BLAS( TestElem(..), Dim(..), Index(..), Assocs(..),
    VectorPair(..) )
import qualified Test.QuickCheck.BLAS as Test

import Typed


tests_STVector = testGroup "STVector"
    [ testPropertyI "new_" prop_new_
    , testPropertyI "new" prop_new
    , testPropertyI "newCopy" prop_newCopy
    , testPropertyI "copyTo" prop_copyTo
    , testPropertyDZ "swap" prop_swap prop_swap
    , testPropertyI "read" prop_read
    , testPropertyI "write" prop_write
    ]
    
    
prop_new_ t (Dim n) = 
    (dimVector $ typed t $ runVector $ newVector_ n) === n

prop_new t (Dim n) e = 
    (runVector $ newVector n e) === (typed t $ constantVector n e)
    
prop_newCopy t x = 
    (runVector $ newCopyVector x) === x
  where
    _ = typed t x
        
prop_copyTo t (VectorPair x y) = runST $
    x `readOnlyVector` \mx ->
    y `mutatesToVector` x $ \my ->
        copyToVector mx my
  where
    _ = typed t x

prop_swap t (VectorPair x y) = runST $
    x `mutatesToVector` y $ \mx ->
    y `mutatesToVector` x $ \my ->
        swapVector mx my
  where
    _ = typed t x
    
prop_read t (Index n i) =
    forAll (typed t `fmap` Test.vector n) $ \x -> runST $
        x `readOnlyVector` \mx -> do
            e <- readVector mx i
            return $ e === atVector x i

prop_write t (Index n i) e =
    forAll (Test.vector n) $ \x -> runST $
        x `mutatesToVector` (x `replaceVector` [(i,e)]) $ \mx -> do
            writeVector mx i e
  where
    _ = e == t


readOnlyVector :: (Storable e, AEq e, Testable prop)
               => Vector e
               -> (STVector s e ->  ST s prop)
               -> ST s Property
readOnlyVector x = mutatesToVector x x

mutatesToVector :: (Storable e, AEq e, Testable prop)
                => Vector e
                -> Vector e
                -> (STVector s e -> ST s prop)
                -> ST s Property
mutatesToVector x x_new f = do
    mx <- thawVector x
    prop <- f mx
    x' <- freezeVector mx
    return $ prop .&. (x' === x_new)
