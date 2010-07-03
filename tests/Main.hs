
import Test.Framework
import Test.Framework.Providers.QuickCheck2

import Vector
import STVector

main :: IO ()
main = defaultMain tests
  where
    tests = [ tests_Vector
            , tests_STVector
            ]
