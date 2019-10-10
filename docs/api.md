# Table of Contents

  * [petibmpy](#petibmpy)
  * [petibmpy.probes](#petibmpy.probes)
    * [read\_probe\_volume\_hdf5](#petibmpy.probes.read_probe_volume_hdf5)
    * [get\_probe\_volume\_limits](#petibmpy.probes.get_probe_volume_limits)
  * [petibmpy.misc](#petibmpy.misc)
    * [Sequence](#petibmpy.misc.Sequence)
    * [\_represent\_dictionary\_order](#petibmpy.misc._represent_dictionary_order)
    * [\_represent\_limits](#petibmpy.misc._represent_limits)
    * [setup\_yaml](#petibmpy.misc.setup_yaml)
  * [petibmpy.forces](#petibmpy.forces)
    * [read\_forces](#petibmpy.forces.read_forces)
    * [get\_force\_coefficients](#petibmpy.forces.get_force_coefficients)
    * [get\_time\_averaged\_values](#petibmpy.forces.get_time_averaged_values)
  * [petibmpy.version](#petibmpy.version)
    * [\_version\_major](#petibmpy.version._version_major)
    * [\_version\_minor](#petibmpy.version._version_minor)
    * [\_version\_micro](#petibmpy.version._version_micro)
    * [\_version\_extra](#petibmpy.version._version_extra)
    * [\_ver](#petibmpy.version._ver)
    * [\_\_version\_\_](#petibmpy.version.__version__)
    * [CLASSIFIERS](#petibmpy.version.CLASSIFIERS)
    * [description](#petibmpy.version.description)
    * [long\_description](#petibmpy.version.long_description)
    * [NAME](#petibmpy.version.NAME)
    * [MAINTAINER](#petibmpy.version.MAINTAINER)
    * [MAINTAINER\_EMAIL](#petibmpy.version.MAINTAINER_EMAIL)
    * [DESCRIPTION](#petibmpy.version.DESCRIPTION)
    * [LONG\_DESCRIPTION](#petibmpy.version.LONG_DESCRIPTION)
    * [URL](#petibmpy.version.URL)
    * [DOWNLOAD\_URL](#petibmpy.version.DOWNLOAD_URL)
    * [LICENSE](#petibmpy.version.LICENSE)
    * [AUTHOR](#petibmpy.version.AUTHOR)
    * [AUTHOR\_EMAIL](#petibmpy.version.AUTHOR_EMAIL)
    * [PLATFORMS](#petibmpy.version.PLATFORMS)
    * [MAJOR](#petibmpy.version.MAJOR)
    * [MINOR](#petibmpy.version.MINOR)
    * [MICRO](#petibmpy.version.MICRO)
    * [VERSION](#petibmpy.version.VERSION)
    * [PACKAGES](#petibmpy.version.PACKAGES)
    * [PACKAGE\_DATA](#petibmpy.version.PACKAGE_DATA)
    * [REQUIRES](#petibmpy.version.REQUIRES)
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

# `petibmpy`


# `petibmpy.probes`

Module for the PetIBM probes.

## `read_probe_volume_hdf5()`

```python
def read_probe_volume_hdf5(filepath, name, time)
```

Read a volume field from file at a given recorded time.

Parameters
----------
filepath : pathlib.Path or string
    Path of the file to read from.
name : string
    Name of the field to read.
time : float
    Time value at which the volume field was recorded.

Returns
-------
mesh : tuple of numpy.ndarray objects
    The gridline coordinates of the volume.
values : numpy.ndarray
    The volume field.

## `get_probe_volume_limits()`

```python
def get_probe_volume_limits(x, loc, neighbors=1, btype='both', decimals=6)
```

Return the interval limits (left, right, or both).

Parameters
----------
x : 1D numpy.ndarray
    The gridline coordinates.
loc : float
    Location for which to request the gridline boundaries.
neighbors : integer (optional)
    Number of neighbors to include as a buffer; default: 1.
btype : string (optional)
    Type of boundary to return;
    choices: 'left', 'right, or 'both';
    default: 'both'.
decimals : integer (optional)
    Number of digits used to round numbers; default: 6.

Returns
-------
'left', 'right', or 'both' gridline boundaries.

# `petibmpy.misc`

Collection of miscellaneous functions and classes.

## `Sequence` Objects

Dummy class to store list/tuple in YAML file in pretty format.

## `_represent_dictionary_order()`

```python
def _represent_dictionary_order(self, dict_data)
```

Pretty output of dictionary to YAML file.

## `_represent_limits()`

```python
def _represent_limits(self, data)
```

Pretty output of list/tuple to YAML file.

## `setup_yaml()`

```python
def setup_yaml()
```

Configure output format to YAML file.

# `petibmpy.forces`

Module with functions to process forces..

## `read_forces()`

```python
def read_forces(filepaths)
```

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

## `get_force_coefficients()`

```python
def get_force_coefficients(forces, *,, ,, =)
```

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

## `get_time_averaged_values()`

```python
def get_time_averaged_values(t, forces, *,, ,, =)
```

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

# `petibmpy.version`

Set up the version.

## `_version_major`

```python
_version_major = 0
```


## `_version_minor`

```python
_version_minor = 1
```


## `_version_micro`

```python
_version_micro = ''
```


## `_version_extra`

```python
_version_extra = ''
```


## `_ver`

```python
_ver = [_version_major, _version_minor]
```


## `__version__`

```python
__version__ = '.'.join(map(str, _ver))
```


## `CLASSIFIERS`

```python
CLASSIFIERS = ['Development Status :: 1 - Alpha',
               'Environment :: Console',
               'License ...
```


## `description`

```python
description = 'PetibmPy: Python processing tools for PetIBM'
```


## `long_description`

```python
long_description = """
PetibmPy
========
Python processing tools for PetIBM.

License
=======
PetibmPy is licensed unde ...
```


## `NAME`

```python
NAME = 'PetibmPy'
```


## `MAINTAINER`

```python
MAINTAINER = 'Olivier Mesnard'
```


## `MAINTAINER_EMAIL`

```python
MAINTAINER_EMAIL = 'mesnardo@gwu.edu'
```


## `DESCRIPTION`

```python
DESCRIPTION = description
```


## `LONG_DESCRIPTION`

```python
LONG_DESCRIPTION = long_description
```


## `URL`

```python
URL = 'https://github.com/mesnardo/petibmpy'
```


## `DOWNLOAD_URL`

```python
DOWNLOAD_URL = ''
```


## `LICENSE`

```python
LICENSE = 'BSD 3-Clause'
```


## `AUTHOR`

```python
AUTHOR = ''
```


## `AUTHOR_EMAIL`

```python
AUTHOR_EMAIL = ''
```


## `PLATFORMS`

```python
PLATFORMS = 'Unix'
```


## `MAJOR`

```python
MAJOR = _version_major
```


## `MINOR`

```python
MINOR = _version_minor
```


## `MICRO`

```python
MICRO = _version_micro
```


## `VERSION`

```python
VERSION = __version__
```


## `PACKAGES`

```python
PACKAGES = ['petibmpy']
```


## `PACKAGE_DATA`

```python
PACKAGE_DATA = {'petibmpy': [os.path.join('styles', '*')]}
```


## `REQUIRES`

```python
REQUIRES = ['h5py', 'lxml', 'matplotlib', 'numpy', 'pyyaml', 'scipy']
```


# `petibmpy.field`

Module to read/write a PetIBM field variable.

## `read_field_hdf5()`

```python
def read_field_hdf5(filepath, name)
```

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

## `write_field_hdf5()`

```python
def write_field_hdf5(filepath, name, field)
```

Write a field to a HDF5 file.

Parameters
----------
filepath : string or pathlib.Path object
    Path of the HDF5 file.
name : string
    Name of the field variable.
field : numpy.ndarray
    The PetIBM field variable as a NumPy array of floats.

## `linear_interpolation()`

```python
def linear_interpolation(u, x, xi)
```

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

## `interpolate3d()`

```python
def interpolate3d(field, grid1, grid2, kwargs)
```

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

## `interpolate2d()`

```python
def interpolate2d(field, grid1, grid2, kwargs)
```

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

# `petibmpy.extrude`

Module with function to extrude a 2D geometry in the third direction.

## `extrude2d()`

```python
def extrude2d(x, y, limits=[-0.5, 0.5], n=None, ds=None, force=False)
```

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

# `petibmpy.rotate`

Module with function to rotate a geometry.

## `rotate2d()`

```python
def rotate2d(x, y, center=(0.0, 0.0), angle=0.0, mode='deg')
```

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

## `rotate3d()`

```python
def rotate3d(x, y, z, roll=0.0, yaw=0.0, pitch=0.0, center=(0.0, 0.0, 0.0))
```

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

## `rotate3d_vec`

```python
rotate3d_vec = numpy.vectorize(rotate3d,
                               excluded=['roll', 'yaw', 'pitch', 'center'] ...
```


# `petibmpy.qcriterion`

Module with functions to compute the Q-criterion

## `qcriterion()`

```python
def qcriterion(velocity, grid)
```

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

# `petibmpy.regularize`

Module with function to regularize a 2D curve (with uniform resolution).

## `get_perimeter()`

```python
def get_perimeter(x, y)
```

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

## `regularize2d()`

```python
def regularize2d(xo, yo, N=None, ds=None, atol=1.0E-06)
```

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

# `petibmpy.grid`

Module to create/read/write a PetIBM grid.

## `CartesianGrid` Objects

```python
def __init__(self, config=None)
```

Contain information about a structured Cartesian grid.

### `CartesianGrid.__init__()`

```python
def __init__(self, config=None)
```

Initialize the grid.

Parameters
----------
config : dictionary (optional)
    Configuration of the grid to create; default: None.

### `CartesianGrid.__repr__()`

```python
def __repr__(self, ndigits=6)
```

Representation of the grid.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

### `CartesianGrid.create()`

```python
def create(self, config)
```

Create the grid.

Parameters
----------
config : dictionary
    Configuration of the grid.

### `CartesianGrid.get_number_cells()`

```python
def get_number_cells(self)
```

Return the number of cells in the grid.

### `CartesianGrid.get_gridlines()`

```python
def get_gridlines(self)
```

Return the gridlines as a list of 1D NumPy arrays of floats.

### `CartesianGrid.write_hdf5()`

```python
def write_hdf5(self, filepath)
```

Save the grid into HDF5 file.

Parameters
----------
filepath : pathlib.Path or string
    Path of the HDF5 file to write into.

### `CartesianGrid.write_yaml()`

```python
def write_yaml(self, filepath, ndigits=6)
```

Write the YAML configuration node for PetIBM.

Parameters
----------
filepath : pathlib.Path or string
    Path of the YAML file to write into.
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

## `GridLine` Objects

```python
def __init__(self, config=None)
```

Contain information about a gridline of a structured Cartesian grid.

### `GridLine.__init__()`

```python
def __init__(self, config=None)
```

Initialize the gridline.

Parameters
----------
config : dictionary (optional)
    Configuration of the gridline to create; default: None.

### `GridLine.__repr__()`

```python
def __repr__(self, ndigits=6)
```

Representation of the gridline.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

### `GridLine.create()`

```python
def create(self, config)
```

Create the gridline.

Parameters
----------
config : dictionary
    Configuration of the gridline.

### `GridLine.get_size()`

```python
def get_size(self)
```

Return the number of vertices in the gridline.

### `GridLine.asarray()`

```python
def asarray(self)
```

Return the gridline as a 1D NumPy array of floats.

### `GridLine.yaml_node()`

```python
def yaml_node(self, ndigits=6)
```

Return the YAML configuration node for PetIBM.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

Returns
-------
node : dictionary
    Configuration node for the gridline.

## `Segment` Objects

```python
def __init__(self, config=None)
```

Contain information about a segment of a gridline.

### `Segment.__init__()`

```python
def __init__(self, config=None)
```

Initialize the segment.

Parameters
----------
config : dictionary (optional)
    Configuration of the segment to create; default: None.

### `Segment.__repr__()`

```python
def __repr__(self, ndigits=6)
```

Representation of the segment.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

### `Segment.create()`

```python
def create(self, config)
```

Create the segment.

Parameters
----------
config : dictionary
    Configuration of the segment.

### `Segment.asarray()`

```python
def asarray(self)
```

Return the segment as a 1D NumPy array of floats.

### `Segment.yaml_node()`

```python
def yaml_node(self, ndigits=6)
```

Return the YAML configuration node for PetIBM.

Parameters
----------
ndigits : integer (optional)
    Number of digits to represent floats; default: 6.

Returns
-------
node : dictionary
    Configuration node for the segment.

## `read_grid_hdf5()`

```python
def read_grid_hdf5(filepath, name)
```

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

## `write_grid_hdf5()`

```python
def write_grid_hdf5(filepath, name, grid)
```

Write a grid to a HDF5 file.

Parameters
----------
filepath : string or pathlib.Path object
    Path of the HDF5 file.
name : string
    Name of the grid.
grid : tuple of numpy.ndarray objects
    The gridline coordinates as 1D arrays of floats.

# `petibmpy.bod`

Module with I/O functions for immersed body.

## `write_body()`

```python
def write_body(filepath, coords)
```

Save the boundary coordinates to a file.

Parameters
----------
filepath : pathlib.Path object or string
    Path of the file to write.
coords : tuple of lists or numpy.ndarray objects
    The x, y, and z coordinates of the boundary.

## `read_body()`

```python
def read_body(filepath, kwargs)
```

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

# `petibmpy.createxdmf`

Module to create a XDMF file for a PetIBM field variable.

## `write_xdmf()`

```python
def write_xdmf(outpath, datadir, gridpath, name, nstart=None, nt=None, nsave=None, states=None, times=None)
```

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

## `write_xdmf_multi()`

```python
def write_xdmf_multi(outpath, config, nstart=None, nt=None, nsave=None, states=None, times=None)
```

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

