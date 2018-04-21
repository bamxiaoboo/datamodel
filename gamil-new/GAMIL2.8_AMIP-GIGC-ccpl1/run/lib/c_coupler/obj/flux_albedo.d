flux_albedo.o flux_albedo.d : flux_albedo.F90
flux_albedo.o : shr_orb_mod.o
flux_albedo.o : c_coupler_interface_mod.o
