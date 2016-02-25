.PHONY: graficar
graficar:
	python Shock_analityc.py
	c99 shocktube.c -lm 
	./a.out
	python lee.py datos.txt UpwindGodunov_step_0.dat 


clean:
	rm *.txt *.dat a.out 
