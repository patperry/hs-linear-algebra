#!/usr/bin/env runhaskell
> import Distribution.Simple
> import System.Cmd
>
> testing _ _ _ _ = do
>     system "runhaskell -lcblas -DREAL    tests/Vector.hs"
>     system "runhaskell -lcblas -DCOMPLEX tests/Vector.hs"
>     system "runhaskell -lcblas -DREAL    tests/Matrix.hs"
>     system "runhaskell -lcblas -DCOMPLEX tests/Matrix.hs"
>     system "runhaskell -lcblas -DREAL    tests/HermMatrix.hs"
>     system "runhaskell -lcblas -DCOMPLEX tests/HermMatrix.hs"
>     system "runhaskell -lcblas -DREAL    tests/TriMatrix.hs"
>     system "runhaskell -lcblas -DCOMPLEX tests/TriMatrix.hs"
>     system "runhaskell -lcblas -DREAL    tests/Banded.hs"
>     system "runhaskell -lcblas -DCOMPLEX tests/Banded.hs"
>     system "runhaskell -lcblas -DREAL    tests/HermBanded.hs"
>     system "runhaskell -lcblas -DCOMPLEX tests/HermBanded.hs"
>     system "runhaskell -lcblas -DREAL    tests/TriBanded.hs"
>     system "runhaskell -lcblas -DCOMPLEX tests/TriBanded.hs"
>     system "runhaskell -lcblas -DREAL    tests/Perm.hs"
>     system "runhaskell -lcblas -DCOMPLEX tests/Perm.hs"
>     system "runhaskell -lcblas -DREAL    tests/Diag.hs"
>     system "runhaskell -lcblas -DCOMPLEX tests/Diag.hs"
>     return ()
>
> main = defaultMainWithHooks defaultUserHooks{ runTests=testing }
>
