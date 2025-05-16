

#makefile pour le projet


main: main.o moindre_carrees.o #dependances
	gfortran -o main main.o moindre_carrees.o #commande

main.o: main.f90 moindre_carrees.o
	gfortran  -c main.f90  

moindre_carrees.o: moindre_carrees.f90 
	gfortran  -c moindre_carrees.f90

clean:
	rm *.o *.mod 

clean_exec: 
	rm main

clean_data:
	rm moindre_carrees.txt
