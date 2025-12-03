# A fluid dynamics simulation

*Authors* Pedde, Veraldi, Zhaokun

## MILESTONES

- [x] add documentation doxygen style
- [x] CMAKE
- [x] add openMP TEST IT
- [x] put on leonardo and do benchmarks for serial
- [x] MPI implementation (WHEN PUT A VARIABLE PUT IT ON PARALLEL IN MODULE PARAMETERS FILES)
	- [x] model.f90 init MPI (1D GRID, IMPLEMENT REST)
		- rank
		- size
  		- ierr 	
		- x_global (not needed)
		- x_local  (not needed)
		- z_global
		- z_local
		- nz_loc (MEMORY ALONG X IS CONTIGOUS)
		- comm
	- [x] module_physics.f90 MPI runge-kutta (all inside it)
	- [x] module_output.F90 MPI (read doc netCDF for the communicator, put right cursor position, block)
	- [x] module_types.F90 MPI (exchange halos, trend)

- [x] TEST MPI
- [x] Timer
- [x] Benchmarks MPI with and without threads
- [x] add namelist for fortran parameters
- [x] openACC
- [x] CMAKE should have option for using threads (compile with openMP) or openACC (nvfortran)
- [ ] Benchmarks MPI with openacc
- [ ] README update of code how works
- [ ] Doxygen generation of docs
- [ ] presentation (see google)

