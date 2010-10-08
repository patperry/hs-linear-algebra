module Typed(
    typed,
    
    testPropertyI,
    testPropertyD,
    testPropertyZ,
    testPropertyDZ,
    
    immutableVector,
    immutableMatrix,
    ) where


import Data.Complex( Complex )
import Test.Framework
import Test.Framework.Providers.QuickCheck2
import Test.QuickCheck

import Numeric.LinearAlgebra.Vector( Vector )
import Numeric.LinearAlgebra.Matrix( Matrix )
        

typed :: e -> a e -> a e
typed _ = id

immutableVector :: Vector e -> Vector e
immutableVector = id

immutableMatrix :: Matrix e -> Matrix e
immutableMatrix = id

testPropertyI :: (Testable a)
              => TestName
              -> (Int -> a)
              -> Test
testPropertyI str prop =
    testProperty str $ prop undefined

testPropertyDZ :: (Testable a, Testable b)
               => TestName
               -> (Double -> a)
               -> (Complex Double -> b)
               -> Test
testPropertyDZ str propd propz =
    testGroup str
        [ testProperty "Double" $ propd undefined
        , testProperty "Complex Double" $ propz undefined
        ]

testPropertyD :: (Testable a)
              => TestName
              -> (Double -> a)
              -> Test
testPropertyD str prop =
    testProperty str $ prop undefined

testPropertyZ :: (Testable a)
              => TestName
              -> (Complex Double -> a)
              -> Test
testPropertyZ str prop =
    testProperty str $ prop undefined
        