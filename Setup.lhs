#!/usr/bin/env runhaskell
> import Distribution.Simple
> import System.Cmd
>
> testing _ _ _ _ = do
>     system "runhaskell -lblas -DREAL    tests/Vector.hs"
>     system "runhaskell -lblas -DCOMPLEX tests/Vector.hs"
>     system "runhaskell -lblas -DREAL    tests/Matrix.hs"
>     system "runhaskell -lblas -DCOMPLEX tests/Matrix.hs"
>     system "runhaskell -lblas -DREAL    tests/TriMatrix.hs"
>     system "runhaskell -lblas -DCOMPLEX tests/TriMatrix.hs"
>     return ()
>
> main = defaultMainWithHooks defaultUserHooks{ runTests=testing }
>
