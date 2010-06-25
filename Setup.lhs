#!/usr/bin/env runhaskell
> import Distribution.Simple
> import Distribution.Simple.LocalBuildInfo
> import Distribution.PackageDescription
> import System.Cmd
> import System.FilePath
>
> main = defaultMainWithHooks autoconfUserHooks { runTests = runTests' }
>
> runTests' :: Args -> Bool -> PackageDescription -> LocalBuildInfo -> IO ()
> runTests' _ _ _ lbi = system testprog >> return ()
>   where testprog = (buildDir lbi) </> "test" </> "test"
