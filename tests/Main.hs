
import System.Environment

import Driver

import Vector
import STVector
import Matrix
import STMatrix
import Banded
import HermMatrix
import HermBanded
import TriMatrix
import TriBanded

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
            , ("Banded"      , tests_Banded)
            , ("Herm Matrix" , tests_HermMatrix)
            , ("Herm Banded" , tests_HermBanded)
            , ("Tri Matrix"  , tests_TriMatrix)
            , ("Tri Banded"  , tests_TriBanded)
            ]
