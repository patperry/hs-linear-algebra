# To build and run the tests, make sure you have blas installed and that you
# also have the "pqc" package installed.


# Change NTHREADS to the number of cores on your machine to run the tests
# in parallel.
NTHREADS = 1
THREADS  = $(NTHREADS) +RTS -N$(NTHREADS) -RTS

all:
	ghc --make -fforce-recomp -O -DREAL -threaded Vector.hs -o vector-real
	./vector-real $(THREADS)
	
	ghc --make -fforce-recomp -O -DCOMPLEX -threaded Vector.hs -o vector-cplx
	./vector-cplx $(THREADS)
	
	ghc --make -fforce-recomp -O -DREAL -threaded Matrix.hs -o matrix-real
	./matrix-real $(THREADS)
	
	ghc --make -fforce-recomp -O -DCOMPLEX -threaded Matrix.hs -o matrix-cplx
	./matrix-cplx $(THREADS)
	
	ghc --make -fforce-recomp -O -DREAL -threaded HermMatrix.hs -o herm-matrix-real
	./herm-matrix-real $(THREADS)
	
	ghc --make -fforce-recomp -O -DCOMPLEX -threaded HermMatrix.hs -o herm-matrix-cplx
	./herm-matrix-cplx $(THREADS)
	
	ghc --make -fforce-recomp -O -DREAL -threaded TriMatrix.hs -o tri-matrix-real
	./tri-matrix-real $(THREADS)
	
	ghc --make -fforce-recomp -O -DCOMPLEX -threaded TriMatrix.hs -o tri-matrix-cplx
	./tri-matrix-cplx $(THREADS)

	ghc --make -fforce-recomp -O -DREAL -threaded Banded.hs -o banded-real
	./banded-real $(THREADS)
	
	ghc --make -fforce-recomp -O -DCOMPLEX -threaded Banded.hs -o banded-cplx
	./banded-cplx $(THREADS)

	ghc --make -fforce-recomp -O -DREAL -threaded HermBanded.hs -o herm-banded-real
	./herm-banded-real $(THREADS)
	
	ghc --make -fforce-recomp -O -DCOMPLEX -threaded HermBanded.hs -o herm-banded-cplx
	./herm-banded-cplx $(THREADS)

	ghc --make -fforce-recomp -O -DREAL -threaded TriBanded.hs -o tri-banded-real
	./tri-banded-real $(THREADS)
	
	ghc --make -fforce-recomp -O -DCOMPLEX -threaded TriBanded.hs -o tri-banded-cplx
	./tri-banded-cplx $(THREADS)

	ghc --make -fforce-recomp -O -DREAL -threaded Perm.hs -o perm-real
	./perm-real $(THREADS)
	
	ghc --make -fforce-recomp -O -DCOMPLEX -threaded Perm.hs -o perm-cplx
	./perm-cplx $(THREADS)


clean:
	rm -f *~ *.hi *.o vector-real vector-cplx matrix-real matrix-cplx \
	herm-matrix-real herm-matrix-cplx tri-matrix-real tri-matrix-cplx \
	banded-real banded-cplx herm-banded-real herm-banded-cplx \
	tri-banded-real tri-banded-cplx perm-real perm-cplx
