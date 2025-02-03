This tool generates a 1-D mesh in the right format for the VKI stagnation-line code.

It takes as input the followings:

    * stagnation point radius (distance from body center): real 
    * free-stream boundary radius (distance from body center): real
    * number of volumes for the mesh: integer
    * stretching function (parabola, hyperbolic tangent and hyperbolic sine are implemented): string
    * stretching intensity (~0 for equispaced mesh, >0 to a refinement next to the wall): real
         note: the stretching intensity has to be greater than 0 (use 0.000000001 to have an equispaced mesh) 

It generates the following output:

   * mesh.dat: is the file containing the mesh to be used with the VKI stagnation-line code
   * mesh.tec: is a file containing two columns: the cell index and the stencil size. Use this file to check the mesh spacing.
