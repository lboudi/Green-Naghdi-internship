REPS = src 
REPSCLEAN = src obj
.PHONY:reps $(REPS)
reps: $(REPS)
$(REPS):
	@cd $@; make SOURCEf="f90" SOURCEoldf="f"

clean:
	@for dir in $(REPSCLEAN); do \
	cd $$dir; make clean; cd ..; \
	done
	rm -f surf.exe core *~ MODULES/*
reset:
	@cd res; make reset; cd ..;

