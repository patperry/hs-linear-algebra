

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
	find ../lib . -name '*.hi' | xargs rm -f
	find ../lib . -name '*.o'  | xargs rm -f
	find . -name '*.html' | xargs rm -f
	rm -f test-real test-complex test-complex.tix
