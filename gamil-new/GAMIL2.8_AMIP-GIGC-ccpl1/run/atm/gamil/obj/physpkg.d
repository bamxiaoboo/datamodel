physpkg.o physpkg.d : physpkg.F90
physpkg.o : misc.h
physpkg.o : params.h
physpkg.o : shr_kind_mod.o
physpkg.o : pmgrid.o
physpkg.o : ppgrid.o
physpkg.o : buffer.o
physpkg.o : comsrf.o
physpkg.o : mpishorthand.o
physpkg.o : phys_grid.o
physpkg.o : physics_types.o
physpkg.o : diagnostics.o
physpkg.o : coupling_chemistry_model_mod.o
physpkg.o : phys_buffer.o
physpkg.o : comctl.h
physpkg.o : comsol.h
