pwd
cpp -traditional-cpp -traditional -DHAVE_CONFIG_H at_experiments.fpp at_experiments.f90
mpif90  -I. -I../../src  -I../lib -O3 -ffast-math -march=native -ffree-line-length-none -fopenmp  -c at_experiments.f90
cpp -traditional-cpp -traditional -DHAVE_CONFIG_H at_md_vverlet.fpp at_md_vverlet.f90
mpif90  -I. -I../../src  -I../lib -O3 -ffast-math -march=native -ffree-line-length-none -fopenmp  -c at_md_vverlet.f90
cpp -traditional-cpp -traditional -DHAVE_CONFIG_H at_dynamics.fpp at_dynamics.f90
mpif90  -I. -I../../src  -I../lib -O3 -ffast-math -march=native -ffree-line-length-none -fopenmp  -c at_dynamics.f90
cpp -traditional-cpp -traditional -DHAVE_CONFIG_H atdyn.fpp atdyn.f90
mpif90  -I. -I../../src  -I../lib -O3 -ffast-math -march=native -ffree-line-length-none -fopenmp  -c atdyn.f90
mpif90 -o atdyn at_energy_str.o at_enefunc_str.o at_pairlist_str.o at_boundary_str.o at_constraints_str.o at_experiments_str.o at_restraints_str.o at_ensemble_str.o at_dynvars_str.o at_dynamics_str.o at_minimize_str.o at_vibration_str.o at_output_str.o at_remd_str.o at_rpath_str.o at_qmmm.o at_experiments.o at_pairlist.o at_boundary.o at_energy_bonds.o at_energy_angles.o at_energy_dihedrals.o at_energy_table_cubic.o at_energy_table_linear.o at_energy_table_linear_bondcorr.o at_energy_pme.o at_energy_eef1.o at_energy_gbsa.o at_energy_nonbonds.o at_energy_restraints.o at_energy_go.o at_energy_gamd.o at_energy.o at_enefunc_restraints.o at_enefunc_gamd.o at_enefunc_gbsa.o at_enefunc_table.o at_enefunc_pme.o at_enefunc_charmm.o at_enefunc_amber.o at_enefunc_gromacs.o at_enefunc_go.o at_enefunc.o at_constraints.o at_restraints.o at_ensemble.o at_dynvars.o at_input.o at_output.o at_md_leapfrog.o at_md_vverlet.o at_minimize.o at_vibration.o at_remd.o at_rpath_mep.o at_rpath_fep.o at_rpath.o at_gamd.o at_dynamics.o at_control.o at_setup_atdyn.o at_setup_mpi.o atdyn.o ../lib/lib.a -fopenmp  -llapack -lblas 
mv atdyn ../../bin
exit

