#!/usr/bin/env runhaskell
> import Distribution.Simple
> import System.Cmd
>
> testing _ _ _ _ = do
>     system "runhaskell -lblas -DREAL    tests/Vector.hs"
>     system "runhaskell -lblas -DCOMPLEX tests/Vector.hs"
>     return ()
>
> main = defaultMainWithHooks defaultUserHooks{ runTests=testing }
>
