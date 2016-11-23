LIBNAME=event

OBJS=	 ksurvb.o ksurvg.o survkit.o

$(LIBNAME)$(SHLIB_EXT): $(OBJS)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) -o $@ $(OBJS) $(FLIBS)

clean:
	@rm -f *.o *.$(SHLIB_EXT)

realclean: clean
