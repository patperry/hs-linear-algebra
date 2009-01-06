-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Banded.Class.Copying
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.Banded.Class.Copying (
    -- * Copying Banded matrices
    newCopyBanded,
    copyBanded,
    unsafeCopyBanded,
    ) where

import Data.Elem.BLAS.Level1( BLAS1 )
import qualified Data.Elem.BLAS.Level1 as BLAS
import BLAS.UnsafeIOToM

import Control.Monad( zipWithM_ )
import Data.Ix( range )
import Foreign( advancePtr )

import Data.Matrix.Class
import Data.Matrix.Banded.Class.Internal
import Data.Matrix.Banded.Class.Views
import Data.Vector.Dense.Base( unsafeCopyVector )

newCopyBanded :: (ReadBanded a e m, WriteBanded b e m) => 
    a mn e -> m (b mn e)
newCopyBanded a 
    | isHermBanded a =
        newCopyBanded ((herm . coerceBanded) a) >>= 
            return . coerceBanded . herm
    | otherwise = do
        a' <- newBanded_ (shapeBanded a) (numLower a, numUpper a)
        unsafeCopyBanded a' a
        return a'

copyBanded :: (WriteBanded b e m, ReadBanded a e m) =>
    b mn e -> a mn e -> m ()
copyBanded dst src
    | shapeBanded dst /= shapeBanded src =
        error "Shape mismatch in copyBanded."
    | bandwidth dst /= bandwidth src =
        error "Bandwidth mismatch in copyBanded."
    | otherwise =
        unsafeCopyBanded dst src

unsafeCopyBanded :: (WriteBanded b e m, ReadBanded a e m) =>
    b mn e -> a mn e -> m ()
unsafeCopyBanded dst src
    | isHermBanded dst = 
        unsafeCopyBanded ((herm . coerceBanded) dst) 
                         ((herm . coerceBanded) src)
                         
    | (not . isHermBanded) src =
        unsafeIOToM $
        withBandedPtr dst $ \pDst ->
        withBandedPtr src $ \pSrc ->
            if ldDst == m && ldSrc == m
                then copyBlock pDst pSrc
                else copyCols  pDst pSrc n
                
    | otherwise =
        zipWithM_ unsafeCopyVector (diagViews dst) (diagViews src)
        
  where
    m     = numLower dst + numUpper dst + 1 -- we can be sure dst is not herm
    n     = numCols dst
    ldDst = ldaOfBanded dst
    ldSrc = ldaOfBanded src

    copyBlock pDst pSrc =
        BLAS.copy (m*n) pSrc 1 pDst 1

    copyCols pDst pSrc nleft
        | nleft == 0 = return ()
        | otherwise = do
            BLAS.copy m pSrc 1 pDst 1
            copyCols (pDst `advancePtr` ldDst) (pSrc `advancePtr` ldSrc) 
                     (nleft-1)

    diagViews a = map (unsafeDiagViewBanded a) $ (range . bandwidth) a
    