

all:
	ghc --make -fforce-recomp -XScopedTypeVariables -DREAL Main.hs -o test-real
	./test-real
	
	ghc --make -fforce-recomp -XScopedTypeVariables -DCOMPLEX Main.hs -o test-complex
	./test-complex
	
hpc:
	ghc --make -fforce-recomp -i. -i.. -XScopedTypeVariables -lcblas -fhpc -DCOMPLEX Main.hs -o test-complex
	./test-complex
	hpc markup test-complex
	
clean:
	find ../BLAS ../Data . -name '*.hi' | xargs rm -f
	find ../BLAS ../Data . -name '*.o'  | xargs rm -f
	find . -name '*.html' | xargs rm -f
	rm -f test-real test-complex test-complex.tix
