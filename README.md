# pipeMeshNek

f90 program for the creation of a nek5000 .rea file for the creation of a good
mesh for a straight pipe

The software generates a mesh based on the parameters in the file
*INPUTgeometry*. One example file is provided.

Both a 2D and a 3D mesh are generated in *base2d.rea* and *base.rea*
respectively.

## Testing:
The created mesh can be visualised by running the python script *plotMesh.py*.
The script directly reads the .rea file *base2d.rea* and can be run simply
by executing

`python plotMesh.py`

The [pymech](https://github.com/jcanton/pymech) suite (automatically fetched)
is used for reading the .rea file.
