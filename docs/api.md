<a name=".petibmpy"></a>
## petibmpy

<a name=".petibmpy.probes"></a>
## petibmpy.probes

Module for PetIBM probes.

<a name=".petibmpy.probes._ProbeBase"></a>
### \_ProbeBase

```python
class _ProbeBase(object)
```

Base class for a probe.

<a name=".petibmpy.probes._ProbeBase._type"></a>
#### \_type

<a name=".petibmpy.probes._ProbeBase.__init__"></a>
#### \_\_init\_\_(name, field, viewer='hdf5', path=None, n\_monitor=None, n\_sum=None, t\_start=None, t\_end=None)

Initialize a base probe.

Parameters
----------
name : str
    Name of the probe
field : str
    Name of the field variable to monitor
viewer : str, optional
    Type of viewer, by default 'hdf5'
path : pathlib.Path or str, optional
    Path of the output file, by default None
n_monitor : int, optional
    Monitoring frequency, by default None
n_sum : int, optional
    Number of time steps to average, by default None
t_start : float, optional
    Starting time of monitoring, by default None
t_end : float, optional
    Ending time of monitoring, by default None

<a name=".petibmpy.probes._ProbeBase.__repr__"></a>
#### \_\_repr\_\_()

Return the string representation.

Returns
-------
str
    The string representation.

<a name=".petibmpy.probes._ProbeBase.set_viewer"></a>
#### set\_viewer(viewer='hdf5', path=None)

Set the output viewer type and path.

The path is relative to the PetIBM output directory.

Parameters
----------
viewer : str, optional
    Type of viewer, choices are 'hdf5' or 'ascii', by default 'hdf5'
path : pathlib.Path or str, optional
    Path of the output file, by default None

<a name=".petibmpy.probes._ProbeBase._get_yaml_node"></a>
#### \_get\_yaml\_node()

<a name=".petibmpy.probes.ProbeVolume"></a>
### ProbeVolume

```python
class ProbeVolume(_ProbeBase)
```

Class for a volume probe (monitoring solution in sub-volume).

<a name=".petibmpy.probes.ProbeVolume._type"></a>
#### \_type

<a name=".petibmpy.probes.ProbeVolume.__init__"></a>
#### \_\_init\_\_(name, field, box=None, adjust\_box=False, grid=None, \*\*kwargs)

Initialize a volume probe.

Parameters
----------
name : str
    Name of the probe
field : str
    Name of the field variable to monitor
box : list or numpy.ndarray, optional
    Limits of the box, by default None
adjust_box : bool, optional
    Adjust the box given a grid, by default False
grid : list or numpy.ndarray, optional
    The grid of the field, by default None
**kwargs: dict, optional
    Optional arguments passed to the base class constructor

<a name=".petibmpy.probes.ProbeVolume.__repr__"></a>
#### \_\_repr\_\_()

Return the string representation.

Returns
-------
str
    The string representation

<a name=".petibmpy.probes.ProbeVolume._check_type"></a>
#### \_check\_type(cls, ptype)

Check if probe type matches class type.

Parameters
----------
ptype : str
    Type of the probe

Returns
-------
bool
    True if type is 'VOLUME'

<a name=".petibmpy.probes.ProbeVolume._get_yaml_node"></a>
#### \_get\_yaml\_node(ndigits=6, \*\*kwargs)

<a name=".petibmpy.probes.ProbeVolume.adjust_box"></a>
#### adjust\_box(grid, box=None)

Adjust the box so that limits lie between two grid points.

Parameters
----------
grid : list or numpy.ndarray
    The grid of the field to minotor
box : list or numpy.ndarray, optional
    Estimated limits of the box, by default None

<a name=".petibmpy.probes.ProbeVolume.read_hdf5"></a>
#### read\_hdf5(filepath, time, ndigits=6)

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

<a name=".petibmpy.probes.ProbeVolume.read_hdf5_deprecated"></a>
#### read\_hdf5\_deprecated(filepath, time, ndigits=6)

Read the probe from a HDF5 file at a given time.

Method is deprecated and will be removed in next release.

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

<a name=".petibmpy.probes.ProbePoint"></a>
### ProbePoint

```python
class ProbePoint(_ProbeBase)
```

Class to monitor a field at a single point.

<a name=".petibmpy.probes.ProbePoint._type"></a>
#### \_type

<a name=".petibmpy.probes.ProbePoint.__init__"></a>
#### \_\_init\_\_(name, field, loc=None, \*\*kwargs)

Initialize a point probe.

Parameters
----------
name : str
    Name of the probe
field : str
    Name of the field to monitor
loc : list or numpy.ndarray, optional
    Coordinates of the point to monitor, by default None
**kwargs: dict, optional
    Optional arguments passed to the base class constructor

<a name=".petibmpy.probes.ProbePoint.__repr__"></a>
#### \_\_repr\_\_()

Return the string representation.

Returns
-------
str
    The string representation

<a name=".petibmpy.probes.ProbePoint._check_type"></a>
#### \_check\_type(cls, ptype)

Check if probe type matches class type.

Parameters
----------
ptype : str
    Type of the probe

Returns
-------
bool
    True if type is 'POINT'

<a name=".petibmpy.probes.ProbePoint.set_loc"></a>
#### set\_loc(loc)

Set the coordinates of the point to monitor.

Parameters
----------
loc : list or numpy.ndarray
    Coordinates of the point

<a name=".petibmpy.probes.ProbePoint._get_yaml_node"></a>
#### \_get\_yaml\_node(\*\*kwargs)

<a name=".petibmpy.probes.Probe"></a>
#### Probe(ptype, \*args, \*\*kwargs)

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

<a name=".petibmpy.probes.probes_yaml_dump"></a>
#### probes\_yaml\_dump(probes, filepath, mode='w')

Save the probes configuration in a YAML file.

Parameters
----------
probes : list
    The list of probes
filepath : pathlib.Path or str
    Path of the YAML file
mode : str, optional
    Mode to open file, choices are 'w' or 'a', by default 'w'

<a name=".petibmpy.misc"></a>
## petibmpy.misc

Collection of miscellaneous functions and classes.

<a name=".petibmpy.misc.check_not_primary_variables"></a>
#### check\_not\_primary\_variables(f)

Check if variable names are not primary variables.

<a name=".petibmpy.misc.delete_datasets_hdf5"></a>
#### delete\_datasets\_hdf5(filepath, names)

Delete datasets from HDF5 file.

If a name if not a dataset, the function moves to the next name.

Parameters
----------
filepath : pathlib.Path or str
    Path of the HDF5 file.
names : list or tuple
    Names of the datasets to delete.

<a name=".petibmpy.misc._Sequence"></a>
### \_Sequence

```python
class _Sequence(list)
```

Dummy class to store list/tuple in YAML file in pretty format.

<a name=".petibmpy.misc._represent_dictionary_order"></a>
#### \_represent\_dictionary\_order(dict\_data)

Pretty output of dictionary to YAML file.

<a name=".petibmpy.misc._represent_limits"></a>
#### \_represent\_limits(data)

Pretty output of list/tuple to YAML file.

<a name=".petibmpy.misc._setup_yaml"></a>
#### \_setup\_yaml()

Configure output format to YAML file.

<a name=".petibmpy.forces"></a>
## petibmpy.forces

Module with functions to process forces.

<a name=".petibmpy.forces.read_forces"></a>
#### read\_forces(\*filepaths)

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

<a name=".petibmpy.forces.get_force_coefficients"></a>
#### get\_force\_coefficients(\*forces, coeff=1.0)

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

<a name=".petibmpy.forces.get_time_averaged_values"></a>
#### get\_time\_averaged\_values(t, \*forces, limits=(-numpy.infty, numpy.infty))

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

<a name=".petibmpy.forces.get_rms_values"></a>
#### get\_rms\_values(t, \*forces, limits=(-numpy.infty, numpy.infty))

Compute the root-mean-square of the signals.

Parameters
----------
t : numpy.ndarray object
    The time values.
forces : tuple of numpy.ndarray objects
    The forces (or force coefficients).
limits : tuple of 2 floats (optional)
    Time limits used to compute the RMS; default: (-inf, +inf).

Returns
-------
rms : tuple of floats
    The RMS values.

<a name=".petibmpy.version"></a>
## petibmpy.version

Set up the version.

<a name=".petibmpy.version._version_major"></a>
#### \_version\_major

<a name=".petibmpy.version._version_minor"></a>
#### \_version\_minor

<a name=".petibmpy.version._version_micro"></a>
#### \_version\_micro

<a name=".petibmpy.version._version_extra"></a>
#### \_version\_extra

<a name=".petibmpy.version._ver"></a>
#### \_ver

<a name=".petibmpy.version.__version__"></a>
#### \_\_version\_\_

<a name=".petibmpy.version.CLASSIFIERS"></a>
#### CLASSIFIERS

<a name=".petibmpy.version.NAME"></a>
#### NAME

<a name=".petibmpy.version.MAINTAINER"></a>
#### MAINTAINER

<a name=".petibmpy.version.MAINTAINER_EMAIL"></a>
#### MAINTAINER\_EMAIL

<a name=".petibmpy.version.DESCRIPTION"></a>
#### DESCRIPTION

<a name=".petibmpy.version.LONG_DESCRIPTION"></a>
#### LONG\_DESCRIPTION

<a name=".petibmpy.version.URL"></a>
#### URL

<a name=".petibmpy.version.DOWNLOAD_URL"></a>
#### DOWNLOAD\_URL

<a name=".petibmpy.version.LICENSE"></a>
#### LICENSE

<a name=".petibmpy.version.AUTHOR"></a>
#### AUTHOR

<a name=".petibmpy.version.AUTHOR_EMAIL"></a>
#### AUTHOR\_EMAIL

<a name=".petibmpy.version.PLATFORMS"></a>
#### PLATFORMS

<a name=".petibmpy.version.MAJOR"></a>
#### MAJOR

<a name=".petibmpy.version.MINOR"></a>
#### MINOR

<a name=".petibmpy.version.MICRO"></a>
#### MICRO

<a name=".petibmpy.version.VERSION"></a>
#### VERSION

<a name=".petibmpy.version.PACKAGES"></a>
#### PACKAGES

<a name=".petibmpy.version.PACKAGE_DATA"></a>
#### PACKAGE\_DATA

<a name=".petibmpy.version.REQUIRES"></a>
#### REQUIRES

<a name=".petibmpy.vorticit"></a>
## petibmpy.vorticit

Module with functions to compute the vorticity.

<a name=".petibmpy.vorticit.gradient"></a>
#### gradient(u, grid, axis=0)

Compute the gradient of u along a given axis.

Parameters
----------
u : numpy.ndarray
    An N-dimensional array.
grid : tuple of N 1-D arrays
    Grid on which u is defined.
axis : int, optional
    Axis along which to compute the gradient, by default 0.

Returns
-------
numpy.ndarray
    The N-dimensional gradient of u along a given axis.
tuple of N 1-D arrays
    The grid on which is defined the gradient.

<a name=".petibmpy.vorticit._gradient"></a>
#### \_gradient(u, x)

<a name=".petibmpy.vorticit.compute_wx"></a>
#### compute\_wx(v, w, grid\_v, grid\_w)

Compute the x-component of the vorticity field.

Parameters
----------
v : numpy.ndarray
    y-component of the velocity field as a N-dimensional array.
w : numpy.ndarray
    z-component of the velocity field as a N-dimensional array.
grid_v : tuple of N 1-D arrays
    Grid on which the y-velocity is defined.
grid_w : tuple of N 1-D arrays
    Grid on which the z-velocity is defined.

Returns
-------
numpy.ndarray
    x-component of the vorticity field.
tuple of N 1-D arrays
    Grid on which the x-vorticity is defined.

<a name=".petibmpy.vorticit.compute_wy"></a>
#### compute\_wy(u, w, grid\_u, grid\_w)

Compute the y-component of the vorticity field.

Parameters
----------
u : numpy.ndarray
    x-component of the velocity field as a N-dimensional array.
w : numpy.ndarray
    z-component of the velocity field as a N-dimensional array.
grid_x : tuple of N 1-D arrays
    Grid on which the x-velocity is defined.
grid_w : tuple of N 1-D arrays
    Grid on which the z-velocity is defined.

Returns
-------
numpy.ndarray
    y-component of the vorticity field.
tuple of N 1-D arrays
    Grid on which the y-vorticity is defined.

<a name=".petibmpy.vorticit.compute_wz"></a>
#### compute\_wz(u, v, grid\_u, grid\_v)

Compute the z-component of the vorticity field.

Parameters
----------
u : numpy.ndarray
    x-component of the velocity field as a N-dimensional array.
v : numpy.ndarray
    y-component of the velocity field as a N-dimensional array.
grid_u : tuple of N 1-D arrays
    Grid on which the x-velocity is defined.
grid_v : tuple of N 1-D arrays
    Grid on which the y-velocity is defined.

Returns
-------
numpy.ndarray
    z-component of the vorticity field.
tuple of N 1-D arrays
    Grid on which the z-vorticity is defined.

<a name=".petibmpy.field"></a>
## petibmpy.field

Module to read/write a PetIBM field variable.

<a name=".petibmpy.field.read_field_hdf5"></a>
#### read\_field\_hdf5(filepath, name)

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

<a name=".petibmpy.field.write_field_hdf5"></a>
#### write\_field\_hdf5(filepath, name, field)

Write a field to a HDF5 file.

Parameters
----------
filepath : string or pathlib.Path object
    Path of the HDF5 file.
name : string
    Name of the field variable.
field : numpy.ndarray
    The PetIBM field variable as a NumPy array of floats.

<a name=".petibmpy.field.linear_interpolation"></a>
#### linear\_interpolation(u, x, xi)

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

<a name=".petibmpy.field.interpolate3d"></a>
#### interpolate3d(field, grid1, grid2, \*\*kwargs)

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

<a name=".petibmpy.field.interpolate2d"></a>
#### interpolate2d(field, grid1, grid2, \*\*kwargs)

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

<a name=".petibmpy.extrude"></a>
## petibmpy.extrude

Module with function to extrude a 2D geometry in the third direction.

<a name=".petibmpy.extrude.extrude2d"></a>
#### extrude2d(x, y, limits, n=None, ds=None, force=False)

Extrude the two-dimensional section along the third direction (z).

Parameters
----------
x : numpy.ndarray
    x-coordinates of the section.
y : numpy.ndarray
    y-coordinates of the section.
limits : 2-list of floats
    Limits of the extrusion.
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

<a name=".petibmpy.rotate"></a>
## petibmpy.rotate

Module with function to rotate a geometry.

<a name=".petibmpy.rotate.rotate2d"></a>
#### rotate2d(x, y, center=(0.0, 0.0), angle=0.0, mode='deg')

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

<a name=".petibmpy.rotate.rotate3d"></a>
#### rotate3d(x, y, z, roll=0.0, yaw=0.0, pitch=0.0, center=(0.0, 0.0, 0.0))

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

<a name=".petibmpy.rotate.rotate3d_vec"></a>
#### rotate3d\_vec

<a name=".petibmpy.qcriterion"></a>
## petibmpy.qcriterion

Module with functions to compute the Q-criterion

<a name=".petibmpy.qcriterion.qcriterion"></a>
#### qcriterion(velocity, grid)

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

<a name=".petibmpy.regularize"></a>
## petibmpy.regularize

Module with function to regularize a 2D curve (with uniform resolution).

<a name=".petibmpy.regularize._get_perimeter"></a>
#### \_get\_perimeter(x, y)

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

<a name=".petibmpy.regularize.regularize2d"></a>
#### regularize2d(xo, yo, N=None, ds=None, atol=1.0E-06)

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

<a name=".petibmpy.grid"></a>
## petibmpy.grid

Module to create/read/write a PetIBM grid.

<a name=".petibmpy.grid.CartesianGrid"></a>
### CartesianGrid

```python
class CartesianGrid()
```

Contain information about a structured Cartesian grid.

<a name=".petibmpy.grid.CartesianGrid.__init__"></a>
#### \_\_init\_\_(config=None)

Initialize the grid.

Parameters
----------
config : dictionary (optional)
    Configuration of the grid to create; default: None.

<a name=".petibmpy.grid.CartesianGrid.__repr__"></a>
#### \_\_repr\_\_(ndigits=6)

Representation of the grid.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

<a name=".petibmpy.grid.CartesianGrid.create"></a>
#### create(config)

Create the grid.

Parameters
----------
config : dictionary
    Configuration of the grid.

<a name=".petibmpy.grid.CartesianGrid.get_number_cells"></a>
#### get\_number\_cells()

Return the number of cells in the grid.

<a name=".petibmpy.grid.CartesianGrid.get_gridlines"></a>
#### get\_gridlines()

Return the gridlines as a list of 1D NumPy arrays of floats.

<a name=".petibmpy.grid.CartesianGrid.write_hdf5"></a>
#### write\_hdf5(filepath)

Save the grid into HDF5 file.

Parameters
----------
filepath : pathlib.Path or string
    Path of the HDF5 file to write into.

<a name=".petibmpy.grid.CartesianGrid.write_yaml"></a>
#### write\_yaml(filepath, ndigits=6)

Write the YAML configuration node for PetIBM.

Parameters
----------
filepath : pathlib.Path or string
    Path of the YAML file to write into.
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

<a name=".petibmpy.grid.CartesianGrid.plot_gridlines"></a>
#### plot\_gridlines(\*\*kwargs)

<a name=".petibmpy.grid.CartesianGrid.plot_gridlines_2d"></a>
#### plot\_gridlines\_2d(figsize=(6.0, 6.0), color='black', xlabel='x', ylabel='y', xrange=(0, None, 1), yrange=(0, None, 1), xlim=(-numpy.infty, numpy.infty), ylim=(-numpy.infty, numpy.infty))

Create a Matplotlib figure with gridlines.

Parameters
----------
figsize : (float, float), optional
    Width and height of the figure in inches; default is (6, 6).
color : str, optional
    Color of the gridlines; default is black.
xlabel : str, optional
    Label along the x axis; default is 'x'.
ylabel : str, optional
    Label along the y axis; default is 'y'.
xrange : (int, int, int), optional
    Index range (min, max, stride) to consider for x gridlines;
    default is to consider all stations (0, None, 1).
yrange : (int, int, int), optional
    Index range (min, max, stride) to consider for y gridlines;
    default is to consider all stations (0, None, 1).
xlim : (float, float), optional
    Limits of the domain in the x direction to plot;
    default is to plot the entire domain.
ylim : (float, float), optional
    Limits of the domain in the y direction to plot;
    default is to plot the entire domain.

Returns
-------
matplotlib.figure.Figure
    Matplotlib Figure.
matplotlib.axes.Axes
    Matplotlib Axes object.

<a name=".petibmpy.grid.CartesianGrid.plot_gridlines_3d"></a>
#### plot\_gridlines\_3d(figsize=(12.0, 6.0), color='black', xlabel='x', ylabel='y', zlabel='z', xrange=(0, None, 1), yrange=(0, None, 1), zrange=(0, None, 1), xlim=(-numpy.infty, numpy.infty), ylim=(-numpy.infty, numpy.infty), zlim=(-numpy.infty, numpy.infty))

Create a Matplotlib figure with gridlines.

Parameters
----------
figsize : (float, float), optional
    Width and height of the figure in inches; default is (12, 6).
color : str, optional
    Color of the gridlines; default is black.
xlabel : str, optional
    Label along the x axis; default is 'x'.
ylabel : str, optional
    Label along the y axis; default is 'y'.
zlabel : str, optional
    Label along the z axis; default is 'z'.
xrange : (int, int, int), optional
    Index range (min, max, stride) to consider for x gridlines;
    default is to consider all stations (0, None, 1).
yrange : (int, int, int), optional
    Index range (min, max, stride) to consider for y gridlines;
    default is to consider all stations (0, None, 1).
zrange : (int, int, int), optional
    Index range (min, max, stride) to consider for z gridlines;
    default is to consider all stations (0, None, 1).
xlim : (float, float), optional
    Limits of the domain in the x direction to plot;
    default is to plot the entire domain.
ylim : (float, float), optional
    Limits of the domain in the y direction to plot;
    default is to plot the entire domain.
zlim : (float, float), optional
    Limits of the domain in the z direction to plot;
    default is to plot the entire domain.

Returns
-------
matplotlib.figure.Figure
    Matplotlib Figure.
array of matplotlib.axes.Axes
    Array of Matplotlib Axes objects.

<a name=".petibmpy.grid.CartesianGrid._plot_gridlines_2d"></a>
#### \_plot\_gridlines\_2d(ax, x, y, color='black', xrange=(0, None, 1), yrange=(0, None, 1), xlim=(-numpy.infty, numpy.infty), ylim=(-numpy.infty, numpy.infty))

<a name=".petibmpy.grid.CartesianGrid.print_info"></a>
#### print\_info()

Print some information about the cell widths.

The method prints the minimum and maximum cell widths
along each direction, as well as max/min ratio across
directions.

<a name=".petibmpy.grid.GridLine"></a>
### GridLine

```python
class GridLine()
```

Contain information about a gridline of a structured Cartesian grid.

<a name=".petibmpy.grid.GridLine.__init__"></a>
#### \_\_init\_\_(config=None)

Initialize the gridline.

Parameters
----------
config : dictionary (optional)
    Configuration of the gridline to create; default: None.

<a name=".petibmpy.grid.GridLine.__repr__"></a>
#### \_\_repr\_\_(ndigits=6)

Representation of the gridline.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

<a name=".petibmpy.grid.GridLine.create"></a>
#### create(config)

Create the gridline.

Parameters
----------
config : dictionary
    Configuration of the gridline.

<a name=".petibmpy.grid.GridLine.get_size"></a>
#### get\_size()

Return the number of vertices in the gridline.

<a name=".petibmpy.grid.GridLine.asarray"></a>
#### asarray(tol=1e-12)

Return the gridline as a 1D NumPy array of floats.

<a name=".petibmpy.grid.GridLine.yaml_node"></a>
#### yaml\_node(ndigits=6)

Return the YAML configuration node for PetIBM.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

Returns
-------
node : dictionary
    Configuration node for the gridline.

<a name=".petibmpy.grid.GridLine._split_needed"></a>
#### \_split\_needed(config)

Check if need to split a configuration into uniform and stretched.

We only to split the configuration is the last width is bigger than
the target maximum width.

Parameters
----------
config : dict
    Configuration of the segment to split.

Returns
-------
bool
    True is splitting is needed.

<a name=".petibmpy.grid.GridLine._split_uniform_and_stretch"></a>
#### \_split\_uniform\_and\_stretch(config)

Split configuration of a stretched portion.

The configuration is split into a stretch portion and a uniform portion
with a cell width equal to the maximum cell width provided.

Parameters
----------
config : dict
    Configuration of the segment to split.

Returns
-------
dict, dict
    Configurations for the stretched and uniform sub-segments.

<a name=".petibmpy.grid.Segment"></a>
### Segment

```python
class Segment()
```

Contain information about a segment of a gridline.

<a name=".petibmpy.grid.Segment.__init__"></a>
#### \_\_init\_\_(config=None)

Initialize the segment.

Parameters
----------
config : dictionary (optional)
    Configuration of the segment to create; default: None.

<a name=".petibmpy.grid.Segment.__repr__"></a>
#### \_\_repr\_\_(ndigits=6)

Representation of the segment.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

<a name=".petibmpy.grid.Segment.create"></a>
#### create(config)

Create the segment.

Parameters
----------
config : dictionary
    Configuration of the segment.

<a name=".petibmpy.grid.Segment.asarray"></a>
#### asarray()

Return the segment as a 1D NumPy array of floats.

<a name=".petibmpy.grid.Segment.yaml_node"></a>
#### yaml\_node(ndigits=6)

Return the YAML configuration node for PetIBM.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

Returns
-------
node : dictionary
    Configuration node for the segment.

<a name=".petibmpy.grid.read_grid_hdf5"></a>
#### read\_grid\_hdf5(filepath, name)

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

<a name=".petibmpy.grid.write_grid_hdf5"></a>
#### write\_grid\_hdf5(filepath, name, \*grid)

Write a grid to a HDF5 file.

Parameters
----------
filepath : string or pathlib.Path object
    Path of the HDF5 file.
name : string
    Name of the grid.
grid : tuple of numpy.ndarray objects
    The gridline coordinates as 1D arrays of floats.

<a name=".petibmpy.bod"></a>
## petibmpy.bod

Module with I/O functions for immersed body.

<a name=".petibmpy.bod.write_body"></a>
#### write\_body(filepath, \*coords)

Save the boundary coordinates to a file.

Parameters
----------
filepath : pathlib.Path object or string
    Path of the file to write.
coords : tuple of lists or numpy.ndarray objects
    The x, y, and z coordinates of the boundary.

<a name=".petibmpy.bod.read_body"></a>
#### read\_body(filepath, \*\*kwargs)

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

<a name=".petibmpy.logview"></a>
## petibmpy.logview

Module to parse a PETSc log view file.

<a name=".petibmpy.logview.PETScLogView"></a>
### PETScLogView

```python
class PETScLogView(object)
```

Parse a PETSc log view file.

<a name=".petibmpy.logview.PETScLogView.__init__"></a>
#### \_\_init\_\_(filepath=None)

Initialize the parser.

<a name=".petibmpy.logview.PETScLogView.parse_log_view"></a>
#### parse\_log\_view(filepath)

Parse a PETSc log view file.

<a name=".petibmpy.logview.PETScLogView._read_walltime"></a>
#### \_read\_walltime()

Parse and return the wall-time clock in seconds.

<a name=".petibmpy.logview.PETScLogView._read_resident_set_size"></a>
#### \_read\_resident\_set\_size(unit='GB')

Parse and return the resident set size.

<a name=".petibmpy.logview.PETScLogView._read_events"></a>
#### \_read\_events()

Parse information about PETSc events.

<a name=".petibmpy.logview.PETScLogView._parse_event"></a>
#### \_parse\_event(line)

Parse information about an event.

<a name=".petibmpy.logview.plot_events_breakdown"></a>
#### plot\_events\_breakdown(ax, runs, ylabel='wall-time (s)', event\_names=None, bar\_width=0.5)

Add a bar chart of the breakdown of events to an axis.

<a name=".petibmpy.createxdmf"></a>
## petibmpy.createxdmf

Module to create a XDMF file for a PetIBM field variable.

<a name=".petibmpy.createxdmf.write_xdmf"></a>
#### write\_xdmf(outpath, datadir, gridpath, name, nstart=None, nt=None, nsave=None, states=None, times=None)

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

<a name=".petibmpy.createxdmf.write_xdmf_multi"></a>
#### write\_xdmf\_multi(outpath, config, nstart=None, nt=None, nsave=None, states=None, times=None)

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

