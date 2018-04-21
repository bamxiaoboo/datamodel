wrap_mpi.o wrap_mpi.d : wrap_mpi.F90
wrap_mpi.o : misc.h
wrap_mpi.o : shr_kind_mod.o
wrap_mpi.o : mpishorthand.o
wrap_mpi.o : pmgrid.o
wrap_mpi.o : spmd_dyn.o
