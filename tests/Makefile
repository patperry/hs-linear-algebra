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

clean:
	rm -f *~ *.hi *.o vector-real vector-cplx matrix-real matrix-cplx \
	herm-matrix-real herm-matrix-cplx tri-matrix-real tri-matrix-cplx
