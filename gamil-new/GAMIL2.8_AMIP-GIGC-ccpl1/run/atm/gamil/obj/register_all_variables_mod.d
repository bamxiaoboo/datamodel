register_all_variables_mod.o register_all_variables_mod.d : register_all_variables_mod.F90
register_all_variables_mod.o : register_private_variables_mod.o
register_all_variables_mod.o : surface_subroutines_with_coupler_mod.o
