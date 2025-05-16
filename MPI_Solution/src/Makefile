REPDIR = ../obj
REP=src
MOD = ../mod
SOURCEf ="f90"
RANLIB = ranlib
AR = ar cr


F95 = mpif90

F95FLAGS = -O1 -J$(MOD)  -ftree-vectorize -ftree-vectorizer-verbose=0 -lm # -O3 -J$(MOD) -L$(LIBGFORT64) -ftree-vectorize -m64 -fdefault-integer-8 -ftree-vectorizer-verbose=0

# delete and rename a file
RM = rm -f
MV = mv -f
CD = cd
EXE = ../surf.exe
LIB	=	lib.a
ARCH	=	ar -ruv


OBJS  = mod_TYPES.o\
	mod_CST.o\
	mod_COMMON.o\
	mod_TOOLS.o\
	mod_ANALYTIC.o\
	mod_INTEGRATION.o\
	mod_BOUNDARY.o\
	mod_SOLVER.o\
	mod_MUSCL.o\
	mod_WENO.o\
	mod_FLUX.o\
	mod_MESH.o\
	mod_POST.o\
	mod_OUTPUTS.o\
	mod_INITIALISE.o\
	mod_TIME.o\
	surf.o  

OBJS := $(addprefix $(REPDIR)/, $(OBJS))

all: $(EXE)
	@echo "...fin de la compilation de l'executable $(EXE)"	

$(EXE): $(OBJS)
	@echo "fin compilation dans $(REP)"
	@echo "chargement des librairies dans $(REP)"

# MAINLY FOR PGFORTRAN
#	$(F95) $(F95FLAGS) -o $@ $^ -lm $(LIBGOTO32) $(UMFLIB)

	$(F95) $(F95FLAGS) -o $@ $^ 
	@echo "fin creation de l'executable $(EXE) "




#$(REPDIR)/%.o: %.c
#	$(CC) $(CFLAGS) -o $@ -c $<
$(REPDIR)/%.o: %.f90
	$(F95) $(F95FLAGS) -o $@ -c $<
$(REPDIR)/%.o: %.f
	$(F95) $(F95FLAGS) -o $@ -c $<
clean:
	@echo "nettoyage du repertoire $(REP)"
	rm -f core *~ *.o *.a `grep -l built *.f90`
