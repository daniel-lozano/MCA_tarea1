.PHONY: graficar
graficar:
	python Shock_analityc.py
	c99 shocktube.c -lm 
	./a.out
	python lee.py datos.txt UpWindGoudonovfinal.dat 


clean:
	rm *.txt *.dat a.out 
