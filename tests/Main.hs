
import Test.Framework
import Test.Framework.Providers.QuickCheck2

import Vector

main :: IO ()
main = defaultMain tests
  where
    tests = [ tests_Vector
            ]
