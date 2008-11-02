

all:
	ghc --make -fforce-recomp -XPatternSignatures -DREAL Main.hs -o test-real
	./test-real
	
	ghc --make -fforce-recomp -XPatternSignatures -DCOMPLEX Main.hs -o test-complex
	./test-complex
	
hpc:
	ghc --make -fforce-recomp -i. -i.. -XPatternSignatures -lcblas -fhpc -DCOMPLEX Main.hs -o test-complex
	./test-complex
	hpc markup test-complex
	
clean:
	find .. -name '*.hi' | xargs rm -f
	find .. -name '*.o'  | xargs rm -f
	find . -name '*.html' | xargs rm -f
	rm -f test-real test-complex test-complex.tix
