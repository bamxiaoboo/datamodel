mpi_gamil.o mpi_gamil.d : mpi_gamil.F90
mpi_gamil.o : mpishorthand.o
mpi_gamil.o : spmd_dyn.o
mpi_gamil.o : pmgrid.o
