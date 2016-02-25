.PHONY: graficar
graficar:
	cc 3body_proble.c
	./a.out 0.35
	python lee.py evolucion.dat 2 4 6 8 q3 p3 si

	./a.out 0.4325 
	python lee.py evolucion.dat 2 4 6 8 q3 p3 si 

	./a.out 0.425
	python lee.py evolucion.dat 2 4 6 8 q3 p3 no

	./a.out 0.46	
	python lee.py evolucion.dat 2 4 6 8 q3 p3 no

clean:
	rm *.dat *~ *# *.png a.out 
