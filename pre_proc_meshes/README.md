Workflow for downloaded 3D mesh from [NASA Turbulence Modeling Resource](
https://turbmodels.larc.nasa.gov/naca0012numerics_grids.html)

1) Convert mesh from .p3d format to .msh format using [`p3d2gmsh.py`](https://github.com/mrklein/p3d2gmsh)
The easiest way to do this is by renaming both the neutral map file (`.nmf`) and mesh to the same name, such that only the file ending differs. 
Then, call `python3 p3d2gmsh.py SHAREDFILENAME.p3d` which assumes that there is a neutral map file with the same name and the file ending `.nmf` in the same directory present.

This results in a ([Gmsh](https://gmsh.info/)-ready) mesh file `SHAREDFILENAME.msh`.

2) Since the NASA TMR currently provides only neutral map files for the 3D meshes, we need to perform the truncation to 2D ourselves.
For this, the scripts `AssignGeoIDs_SwapCoords.jl` and `Abaqus3Dto2D.jl` are provided.
The first script basically removes the extrusion operation in z-direction and then flips $y$ and $z$ coordinates to get the traditional view on an airfoil mesh.
Additionally, it assigns the correct geometric ID (in addition to the physical ID) of the elements to be able to specify boundary conditions later on.

Thus, to do the first step, call `AssignGeoIDs_SwapCoords.jl` which results in `SHAREDFILENAME_2D.msh`.

3) Now, we need to convert the 2D mesh to `.inp` format for further processing. 
For this we use the gmsh graphical interface to `Export` the mesh with options `Save all elements` and `Save groups of nodes` selected.

This results in a file `SHAREDFILENAME_2D.inp`.

4) Finally, it remains to convert the `.inp` file to truly 2D. 
For this we have the script `Abaqus3Dto2D.jl` which removes the previously 2D elements and changes previous 3D elements to 2D elements.

This results in `SHAREDFILENAME_2D_unique.inp` which is ready for usage in Trixi.jl.