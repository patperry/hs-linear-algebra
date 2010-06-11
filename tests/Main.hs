
import System.Environment

import Driver

import Vector     ( tests_Vector     )
import STVector   ( tests_STVector   )
import Matrix     ( tests_Matrix     )
import STMatrix   ( tests_STMatrix   )
import HermMatrix ( tests_HermMatrix )
import TriMatrix  ( tests_TriMatrix  )
import Banded     ( tests_Banded     )
import HermBanded ( tests_HermBanded )
import TriBanded  ( tests_TriBanded  )


main :: IO ()
main = defaultMain tests
  where
    tests = [ testGroup "Vector" tests_Vector
            , testGroup "STVector" tests_STVector
            , testGroup "Matrix" tests_Matrix
            , testGroup "STMatrix" tests_STMatrix
            , testGroup "Herm Matrix" tests_HermMatrix
            , testGroup "Tri Matrix" tests_TriMatrix
            , testGroup "Banded" tests_Banded
            , testGroup "Herm Banded" tests_HermBanded
            , testGroup "Tri Banded" tests_TriBanded
            ]
