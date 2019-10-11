# Table of Contents

  * [petibmpy](#petibmpy)
  * [petibmpy.probes](#petibmpy.probes)
    * [ProbeVolume](#petibmpy.probes.ProbeVolume)
    * [ProbePoint](#petibmpy.probes.ProbePoint)
    * [Probe](#petibmpy.probes.Probe)
    * [probes\_yaml\_dump](#petibmpy.probes.probes_yaml_dump)
  * [petibmpy.misc](#petibmpy.misc)
  * [petibmpy.forces](#petibmpy.forces)
    * [read\_forces](#petibmpy.forces.read_forces)
    * [get\_force\_coefficients](#petibmpy.forces.get_force_coefficients)
    * [get\_time\_averaged\_values](#petibmpy.forces.get_time_averaged_values)
  * [petibmpy.version](#petibmpy.version)
  * [petibmpy.field](#petibmpy.field)
    * [read\_field\_hdf5](#petibmpy.field.read_field_hdf5)
    * [write\_field\_hdf5](#petibmpy.field.write_field_hdf5)
    * [linear\_interpolation](#petibmpy.field.linear_interpolation)
    * [interpolate3d](#petibmpy.field.interpolate3d)
    * [interpolate2d](#petibmpy.field.interpolate2d)
  * [petibmpy.extrude](#petibmpy.extrude)
    * [extrude2d](#petibmpy.extrude.extrude2d)
  * [petibmpy.rotate](#petibmpy.rotate)
    * [rotate2d](#petibmpy.rotate.rotate2d)
    * [rotate3d](#petibmpy.rotate.rotate3d)
    * [rotate3d\_vec](#petibmpy.rotate.rotate3d_vec)
  * [petibmpy.qcriterion](#petibmpy.qcriterion)
    * [qcriterion](#petibmpy.qcriterion.qcriterion)
  * [petibmpy.regularize](#petibmpy.regularize)
    * [get\_perimeter](#petibmpy.regularize.get_perimeter)
    * [regularize2d](#petibmpy.regularize.regularize2d)
  * [petibmpy.grid](#petibmpy.grid)
    * [CartesianGrid](#petibmpy.grid.CartesianGrid)
    * [GridLine](#petibmpy.grid.GridLine)
    * [Segment](#petibmpy.grid.Segment)
    * [read\_grid\_hdf5](#petibmpy.grid.read_grid_hdf5)
    * [write\_grid\_hdf5](#petibmpy.grid.write_grid_hdf5)
  * [petibmpy.bod](#petibmpy.bod)
    * [write\_body](#petibmpy.bod.write_body)
    * [read\_body](#petibmpy.bod.read_body)
  * [petibmpy.createxdmf](#petibmpy.createxdmf)
    * [write\_xdmf](#petibmpy.createxdmf.write_xdmf)
    * [write\_xdmf\_multi](#petibmpy.createxdmf.write_xdmf_multi)

<h1 id="petibmpy"><code>petibmpy</code></h1>


<h1 id="petibmpy.probes"><code>petibmpy.probes</code></h1>

Module for PetIBM probes.

<h2 id="petibmpy.probes.ProbeVolume"><code>ProbeVolume</code> Objects</h2>

Class for a volume probe (monitoring solution in sub-volume).

<h3 id="petibmpy.probes.ProbeVolume.adjust_box"><code>ProbeVolume.adjust_box(self, grid, box=None)</code></h3>

Adjust the box so that limits lie between two grid points.

Parameters
----------
grid : list or numpy.ndarray
    The grid of the field to minotor
box : list or numpy.ndarray, optional
    Estimated limits of the box, by default None

<h3 id="petibmpy.probes.ProbeVolume.read_hdf5"><code>ProbeVolume.read_hdf5(self, filepath, time, ndigits=6)</code></h3>

Read the probe from a HDF5 file at a given time.

Parameters
----------
filepath : pathlib.Path or str
    Path of file with the solution of the probe
time : float
    Time value
ndigits : int, optional
    Number of digits to round the time value, by default 6

Returns
-------
tuple
    The mesh grid of the probe
numpy.ndarray
    The probe values

<h2 id="petibmpy.probes.ProbePoint"><code>ProbePoint</code> Objects</h2>

Class to monitor a field at a single point.

<h3 id="petibmpy.probes.ProbePoint.set_loc"><code>ProbePoint.set_loc(self, loc)</code></h3>

Set the coordinates of the point to monitor.

Parameters
----------
loc : list or numpy.ndarray
    Coordinates of the point

<h2 id="petibmpy.probes.Probe"><code>Probe(ptype, args, *,, ,, kwargs)</code></h2>

Create a probe.

Parameters
----------
ptype : str
    Type of the probe, choices are 'VOLUME' or 'POINT'

Returns
-------
ProbeVolume or ProbePoint
    The probe

Raises
------
ValueError
    Type is neither 'VOLUME' nor 'POINT'

<h2 id="petibmpy.probes.probes_yaml_dump"><code>probes_yaml_dump(probes, filepath, mode='w')</code></h2>

Save the probes configuration in a YAML file.

Parameters
----------
probes : list
    The list of probes
filepath : pathlib.Path or str
    Path of the YAML file
mode : str, optional
    Mode to open file, choices are 'w' or 'a', by default 'w'

<h1 id="petibmpy.misc"><code>petibmpy.misc</code></h1>

Collection of miscellaneous functions and classes.

<h1 id="petibmpy.forces"><code>petibmpy.forces</code></h1>

Module with functions to process forces..

<h2 id="petibmpy.forces.read_forces"><code>read_forces(filepaths)</code></h2>

Read PetIBM forces from given file(s).

If multiple files are provided, the histories are concatenated.

Parameters
----------
filepaths : tuple of pathlib.Path objects or strings
    Path of the files to load the history from.

Returns
-------
data : numpy.ndarray
    Time followed by the forces in the x, y, and z directions.

<h2 id="petibmpy.forces.get_force_coefficients"><code>get_force_coefficients(forces, *,, ,, =)</code></h2>

Convert forces to force coefficients.

Parameters
----------
forces : tuple of numpy.ndarray objects
    The forces.
coeff : float (optional)
    The scaling coefficient; default: 1.0.

Returns
-------
force_coeffs : tuple of numpy.ndarray objects
    The force coefficients.

<h2 id="petibmpy.forces.get_time_averaged_values"><code>get_time_averaged_values(t, forces, *,, ,, =)</code></h2>

Compute the time-averaged values.

Parameters
----------
t : numpy.ndarray object
    The time values.
forces : tuple of numpy.ndarray objects
    The forces (or force coefficients).
limits : tuple of 2 floats (optional)
    Time limits used to compute the mean; default: (-inf, +inf).

Returns
-------
means : tuple of floats
    The time-averaged values.

<h1 id="petibmpy.version"><code>petibmpy.version</code></h1>

Set up the version.

<h1 id="petibmpy.field"><code>petibmpy.field</code></h1>

Module to read/write a PetIBM field variable.

<h2 id="petibmpy.field.read_field_hdf5"><code>read_field_hdf5(filepath, name)</code></h2>

Read a field from HDF5 file.

Parameters
----------
filepath : string or pathlib.Path object
    Path of the HDF5 file.
name : string
    Name of the field variable.

Returns
-------
field : numpy.ndarray
    The PetIBM field variable as a NumPy array of floats.

<h2 id="petibmpy.field.write_field_hdf5"><code>write_field_hdf5(filepath, name, field)</code></h2>

Write a field to a HDF5 file.

Parameters
----------
filepath : string or pathlib.Path object
    Path of the HDF5 file.
name : string
    Name of the field variable.
field : numpy.ndarray
    The PetIBM field variable as a NumPy array of floats.

<h2 id="petibmpy.field.linear_interpolation"><code>linear_interpolation(u, x, xi)</code></h2>

Perform a linear interpolation along the first axis.

Parameters
----------
u : numpy.ndarray
    Array to interpolate.
x : numpy.ndarray
    Gridline locations.
xi : float
    Target location.

Returns
-------
ui : numpy.ndarray or float
    Interpolated values.

<h2 id="petibmpy.field.interpolate3d"><code>interpolate3d(field, grid1, grid2, kwargs)</code></h2>

Interpolate a 3D field from one grid to another.

Parameters
----------
field : numpy.ndarray
    The 3D field to interpolate.
grid1 : tuple of numpy.ndarray objects
    The grid on which the field is defined.
    The grid should be provided as (x, y, z).
grid2 : tuple of numpy.ndarray objects
    The grid on which to interpolate the field.
    The grid should be provided as (x, y, z).
**kwargs : Arbitrary keyword arguments
    To be passed to scipy.interpolate.interpn.

Returns
-------
field2 : numpy.ndarray
    The interpolated 3D field.

<h2 id="petibmpy.field.interpolate2d"><code>interpolate2d(field, grid1, grid2, kwargs)</code></h2>

Interpolate a 2D field from one grid to another.

Parameters
----------
field : numpy.ndarray
    The 2D field to interpolate.
grid1 : tuple of numpy.ndarray objects
    The grid on which the field is defined.
    The grid should be provided as (x, y).
grid2 : tuple of numpy.ndarray objects
    The grid on which to interpolate the field.
    The grid should be provided as (x, y).
**kwargs : Arbitrary keyword arguments
    To be passed to scipy.interpolate.interpn.

Returns
-------
field2 : numpy.ndarray
    The interpolated 2D field.

<h1 id="petibmpy.extrude"><code>petibmpy.extrude</code></h1>

Module with function to extrude a 2D geometry in the third direction.

<h2 id="petibmpy.extrude.extrude2d"><code>extrude2d(x, y, limits=[-0.5, 0.5], n=None, ds=None, force=False)</code></h2>

Extrude the two-dimensional section along the third direction (z).

Parameters
----------
x : numpy.ndarray
    x-coordinates of the section.
y : numpy.ndarray
    y-coordinates of the section.
limits : 2-list of floats, optional
    Limits of the extrusion;
    default: [-0.5, 0.5].
n : integer, optional
    Number of divisions in the z-direction;
    default: None.
ds : float, optional
    Desired segment-length;
    default: None.
force : boolean, optional
    Forces the extrusion to the limits prescribed;
    default: False.

Returns
-------
x : numpy.ndarray
    x-coordinates of the geometry.
y : numpy.ndarray
    y-coordinates of the geometry.
z : numpy.ndarray
    z-coordinates of the geometry.

<h1 id="petibmpy.rotate"><code>petibmpy.rotate</code></h1>

Module with function to rotate a geometry.

<h2 id="petibmpy.rotate.rotate2d"><code>rotate2d(x, y, center=(0.0, 0.0), angle=0.0, mode='deg')</code></h2>

Rotate (x, y) coordinates around a center.

Parameters
----------
x : numpy.ndarray of floats
    The x-coordinates to rotate.
y : numpy.ndarray of floats
    The y-coordinates to rotate.
center : 2-tuple of floats, optional
    The center of rotation; default: (0.0, 0.0).
angle : float, optional
    The angle of rotation; default: 0.0.
mode : string, optional
    Whether angle is provided in degrees or in radians;
    choices: 'deg', 'rad'; default: 'deg'.

Returns
-------
x_new : numpy.ndarray of floats
    The rotated x-coordinates.
y_new : numpy.ndarray of floats
    The rotated y-coordinates.

<h2 id="petibmpy.rotate.rotate3d"><code>rotate3d(x, y, z, roll=0.0, yaw=0.0, pitch=0.0, center=(0.0, 0.0, 0.0))</code></h2>

Rotate 3D point.

Parameters
----------
x : float
    x-coordinate of point.
y : float
    y-coordinate of point.
z : float
    z-coordinate of point.
roll : float (optional)
    Roll angle (in radians); default: 0.0.
yaw : float (optional)
    Yaw angle (in radians); default: 0.0.
pitch : float (optional)
    Pitch angle (in radians); default: 0.0.
center : tuple of floats
    Coordinates of the center of rotation;
    default: [0.0, 0.0, 0.0].

Returns
-------
xr : float
    x-coordinate of rotated point.
yr : float
    y-coordinate of rotated point.
zr : float
    z-coordinate of rotated point.

<h2 id="petibmpy.rotate.rotate3d_vec"><code>rotate3d_vec</code></h2>


<h1 id="petibmpy.qcriterion"><code>petibmpy.qcriterion</code></h1>

Module with functions to compute the Q-criterion

<h2 id="petibmpy.qcriterion.qcriterion"><code>qcriterion(velocity, grid)</code></h2>

Compute the Q-criterion on a 3D grid.

Parameters
----------
velocity : tuple of numpy.ndarray objects
    The velocity vector field given as (u, v, w).
grid : tuple of numpy.ndarray objects
    The structured Cartesian grid given as (x, y, z).

Returns
-------
qcrit : numpy.ndarray
    Value of the Q-criterion on the 3D grid.

<h1 id="petibmpy.regularize"><code>petibmpy.regularize</code></h1>

Module with function to regularize a 2D curve (with uniform resolution).

<h2 id="petibmpy.regularize.get_perimeter"><code>get_perimeter(x, y)</code></h2>

Return the perimeter of the geometry.

Parameters
----------
x : numpy.ndarray
    x-coordinate of the points along the curve.
y : numpy.ndarray
    y-coordinate of the points along the curve.

Returns
-------
perimeter : float
    The perimeter.

<h2 id="petibmpy.regularize.regularize2d"><code>regularize2d(xo, yo, N=None, ds=None, atol=1.0E-06)</code></h2>

Regularize the geometry.

Parameters
----------
xo: numpy.ndarray of floats
    The x-coordinates of the boundary to regularize.
yo: numpy.ndarray of floats
    The y-coordinates of the boundary to regularize.
N: integer, optional
    Number of divisions;
    default: None.
ds: float, optional
    Desired segment-length;
    default: None.
atol: float, optional
    Desired tolerance for discretization;
    default: 1.0E-06.

Returns
-------
x: numpy.ndarray of floats
    The x-coordinates of the regularized boundary.
y: numpy.ndarray of floats
    The y-coordinates of the regularized boundary.

<h1 id="petibmpy.grid"><code>petibmpy.grid</code></h1>

Module to create/read/write a PetIBM grid.

<h2 id="petibmpy.grid.CartesianGrid"><code>CartesianGrid</code> Objects</h2>

Contain information about a structured Cartesian grid.

<h3 id="petibmpy.grid.CartesianGrid.create"><code>CartesianGrid.create(self, config)</code></h3>

Create the grid.

Parameters
----------
config : dictionary
    Configuration of the grid.

<h3 id="petibmpy.grid.CartesianGrid.get_number_cells"><code>CartesianGrid.get_number_cells(self)</code></h3>

Return the number of cells in the grid.

<h3 id="petibmpy.grid.CartesianGrid.get_gridlines"><code>CartesianGrid.get_gridlines(self)</code></h3>

Return the gridlines as a list of 1D NumPy arrays of floats.

<h3 id="petibmpy.grid.CartesianGrid.write_hdf5"><code>CartesianGrid.write_hdf5(self, filepath)</code></h3>

Save the grid into HDF5 file.

Parameters
----------
filepath : pathlib.Path or string
    Path of the HDF5 file to write into.

<h3 id="petibmpy.grid.CartesianGrid.write_yaml"><code>CartesianGrid.write_yaml(self, filepath, ndigits=6)</code></h3>

Write the YAML configuration node for PetIBM.

Parameters
----------
filepath : pathlib.Path or string
    Path of the YAML file to write into.
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

<h2 id="petibmpy.grid.GridLine"><code>GridLine</code> Objects</h2>

Contain information about a gridline of a structured Cartesian grid.

<h3 id="petibmpy.grid.GridLine.create"><code>GridLine.create(self, config)</code></h3>

Create the gridline.

Parameters
----------
config : dictionary
    Configuration of the gridline.

<h3 id="petibmpy.grid.GridLine.get_size"><code>GridLine.get_size(self)</code></h3>

Return the number of vertices in the gridline.

<h3 id="petibmpy.grid.GridLine.asarray"><code>GridLine.asarray(self)</code></h3>

Return the gridline as a 1D NumPy array of floats.

<h3 id="petibmpy.grid.GridLine.yaml_node"><code>GridLine.yaml_node(self, ndigits=6)</code></h3>

Return the YAML configuration node for PetIBM.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

Returns
-------
node : dictionary
    Configuration node for the gridline.

<h2 id="petibmpy.grid.Segment"><code>Segment</code> Objects</h2>

Contain information about a segment of a gridline.

<h3 id="petibmpy.grid.Segment.create"><code>Segment.create(self, config)</code></h3>

Create the segment.

Parameters
----------
config : dictionary
    Configuration of the segment.

<h3 id="petibmpy.grid.Segment.asarray"><code>Segment.asarray(self)</code></h3>

Return the segment as a 1D NumPy array of floats.

<h3 id="petibmpy.grid.Segment.yaml_node"><code>Segment.yaml_node(self, ndigits=6)</code></h3>

Return the YAML configuration node for PetIBM.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

Returns
-------
node : dictionary
    Configuration node for the segment.

<h2 id="petibmpy.grid.read_grid_hdf5"><code>read_grid_hdf5(filepath, name)</code></h2>

Read a grid from HDF5 file.

Parameters
----------
filepath : string or pathlib.Path object
    Path of the HDF5 file.
name : string
    Name of the grid.

Returns
-------
x : numpy.ndarray
    The x-coordinates along a gridline in the x-direction.
y : numpy.ndarray
    The y-coordinates along a gridline in the y-direction.
z : numpy.ndarray
    The z-coordinates along a gridline in the z-direction.

<h2 id="petibmpy.grid.write_grid_hdf5"><code>write_grid_hdf5(filepath, name, grid)</code></h2>

Write a grid to a HDF5 file.

Parameters
----------
filepath : string or pathlib.Path object
    Path of the HDF5 file.
name : string
    Name of the grid.
grid : tuple of numpy.ndarray objects
    The gridline coordinates as 1D arrays of floats.

<h1 id="petibmpy.bod"><code>petibmpy.bod</code></h1>

Module with I/O functions for immersed body.

<h2 id="petibmpy.bod.write_body"><code>write_body(filepath, coords)</code></h2>

Save the boundary coordinates to a file.

Parameters
----------
filepath : pathlib.Path object or string
    Path of the file to write.
coords : tuple of lists or numpy.ndarray objects
    The x, y, and z coordinates of the boundary.

<h2 id="petibmpy.bod.read_body"><code>read_body(filepath, kwargs)</code></h2>

Read the boundary coordinates from a file.

Parameters
----------
filepath : pathlib.Path object or string
    Path of the file to read.
kwargs : dictionary
    Keyword arguments to pass to numpy.loadtxt.

Returns
-------
coords : numpy.ndarray
    The boundary coordinates.

<h1 id="petibmpy.createxdmf"><code>petibmpy.createxdmf</code></h1>

Module to create a XDMF file for a PetIBM field variable.

<h2 id="petibmpy.createxdmf.write_xdmf"><code>write_xdmf(outpath, datadir, gridpath, name, nstart=None, nt=None, nsave=None, states=None, times=None)</code></h2>

Write a XDMF file to read the solution of a PetIBM variable.

Parameters
----------
outpath : pathlib.Path object
    Path of the XDMF file to create.
datadir : pathlib.Path object
    Data directory.
gridpath : pathlib.Path object
    Path of the file containing the gridline coordinates.
name : string
    Name of the field variable.
nstart : integer (optional)
    Starting time step; default: None.
nt : integer (optional)
    Number of time steps; default: None.
nsave : integer (optional)
    Frequency of saving in number of time steps; default: None.
states : list of integers (optional)
    The list of time-step indices to consider in the XDMF file;
    default: None.
times : list of floats (optional)
    The list of time values; default: None.

<h2 id="petibmpy.createxdmf.write_xdmf_multi"><code>write_xdmf_multi(outpath, config, nstart=None, nt=None, nsave=None, states=None, times=None)</code></h2>

Write a XDMF file to read the solution of multiple PetIBM variables.

Parameters
----------
outpath : pathlib.Path object
    Path of the XDMF file to create.
config : dictionary
    Should contains two keys: 'grid' and 'data'.
    The value mapped to 'grid' is the path of the HDF5 grid file.
    The value mapped to 'data' is a dictionary.
    Each item of the 'data' dictionary is labeled with the name
    of the variable to add to the XDMF file that is mapped to
    the path of the directory that contains the numerical solution
    for that variable.
nstart : integer (optional)
    Starting time step; default: None.
nt : integer (optional)
    Number of time steps; default: None.
nsave : integer (optional)
    Frequency of saving in number of time steps; default: None.
states : list of integers (optional)
    The list of time-step indices to consider in the XDMF file;
    default: None.
times : list of floats (optional)
    The list of time values; default: None.

