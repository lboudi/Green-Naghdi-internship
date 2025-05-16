

#makefile pour le projet


regularisation: regularisation.o module_regularisation.o #dependances
	gfortran -o regularisation regularisation.o module_regularisation.o #commande

regularisation.o: regularisation.f90 module_regularisation.o
	gfortran  -c regularisation.f90  

module_regularisation.o: module_regularisation.f90 
	gfortran  -c module_regularisation.f90

clean:
	rm *.o *.mod 

clean_exec: 
	rm regularisation

clean_data:
	rm *.txt
