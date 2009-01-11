
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
main = do
    args <- getArgs
    let n = if null args then 100 else read (head args)
    
    printf "\nRunnings tests for field `%s'\n" field
    
    (results, passed) <- liftM unzip $
        foldM ( \prev (name,subtests) -> do
                     printf "\n%s\n" name
                     printf "%s\n" $ replicate (length name) '-'
                     cur <- mapM (\(s,a) -> printf "%-30s: " s >> a n) subtests
                     return (prev ++ cur)
              )
              []
              tests
               
    printf "\nPassed %d tests!\n\n" (sum passed)
    when (not . and $ results) $ fail "\nNot all tests passed!"
 where

    tests = [ ("Vector"      , tests_Vector)
            , ("STVector"    , tests_STVector)
            , ("Matrix"      , tests_Matrix)
            , ("STMatrix"    , tests_STMatrix)
            , ("Herm Matrix" , tests_HermMatrix)
            , ("Tri Matrix"  , tests_TriMatrix)
            , ("Banded"      , tests_Banded)
            , ("Herm Banded" , tests_HermBanded)
            , ("Tri Banded"  , tests_TriBanded)
            ]
