

all:
	ghc --make -fforce-recomp -DREAL Main.hs -o test-real
	./test-real
	
	ghc --make -fforce-recomp -DCOMPLEX Main.hs -o test-complex
	./test-complex
	
clean:
	find . -name '*.hi' | xargs rm -f
	find . -name '*.o'  | xargs rm -f
	rm -f test-real test-complex
