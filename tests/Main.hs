
import Test.Framework
import Test.Framework.Providers.QuickCheck2

import Vector
import STVector
-- import Matrix
-- import STMatrix
-- import Statistics

main :: IO ()
main = defaultMain tests
  where
    tests = [ tests_Vector
            , tests_STVector
            -- , tests_Matrix
            -- , tests_STMatrix
            -- , tests_Statistics
            ]
