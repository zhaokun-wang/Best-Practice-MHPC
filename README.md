# A fluid dynamics simulation

*Authors* Pedde, Veraldi, Zhaokun

## MILESTONES

- [x] add documentation doxygen style
- [x] CMAKE **Z**
- [x] add openMP **E** TEST IT
- [x] put on leonardo and do benchmarks for serial (1-2-4-8-16-32 nodes) **G**
- [ ] MPI implementation (WHEN PUT A VARIABLE PUT IT ON PARALLEL IN MODULE PARAMETERS FILES)
	- [ ] model.f90 init MPI (1D GRID, IMPLEMENT REST)
		- rank
		- size
		- x_global (not needed)
		- x_local  (not needed)
		- z_global
		- z_local
		- nz_loc (MEMORY ALONG X IS CONTIGOUS)
		- comm
	- [ ] module_physics.f90 MPI runge-kutta (all inside it)
	- [ ] module_output.F90 MPI (read doc netCDF for the communicator, put right cursor position, block)
	- [ ] module_types.F90 MPI (exchange halos, trend)

- [ ] TEST MPI
- [ ] Benchmarks MPI with and without threads
- [ ] openACC
- [ ] CMAKE should have option for using threads (compile with openMP) or openACC
- [ ] Doxygen generation of docs

