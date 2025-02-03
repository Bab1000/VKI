This tool interpolates a solution of the VKI stagnation-line code over a new mesh.

It takes as input the following files:
     * solution_old.dat: the restart.dat file of the simulation on the old mesh writen by the VKI stagnation-line code
     * mesh_old.dat:     the old mesh file
     * mesh_new.dat:     the new mesh file

It generates the following output:
     * solution_new.dat: the file to be used as restart file for the simulation with the new mesh
     * solution_new.tec: file identical to the solution_new.dat but including the new mesh in the first column. To be used to check the interpolation results.
     * solution_new.tec: file identical to the solution_old.dat but including the old mesh in the first column. To be used to compare with the interpolated result (solution_new.tec).
