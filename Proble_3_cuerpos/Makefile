.PHONY: graficar
graficar:
	cc 3body_proble.c
	./a.out 0.4325 
	python lee.py evolucion.dat 2 4 6 8 q3 p3 

	cc 3body_proble.c		
	./a.out 0.35
	python lee.py evolucion.dat 2 4 6 8 q3 p3 

	cc 3body_proble.c
	./a.out 0.45
	python lee.py evolucion.dat 2 4 6 8 q3 p3 

clean:
	rm *.dat *~ *# *.png a.out 
