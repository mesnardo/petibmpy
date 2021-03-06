"""Module to create a XDMF file for a PetIBM field variable."""

import sys
import pathlib
from lxml import etree

from .grid import read_grid_hdf5


def write_xdmf(outpath, datadir, gridpath, name,
               nstart=None, nt=None, nsave=None,
               states=None, times=None):
    """Write a XDMF file to read the solution of a PetIBM variable.

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

    """
    # Initialize XDMF file.
    xdmf = etree.Element('Xdmf', Version='2.2')
    info = etree.SubElement(xdmf, 'Information',
                            Name='MetaData',
                            Value='ID-23454')
    domain = etree.SubElement(xdmf, 'Domain')
    grid_time_series = etree.SubElement(domain, 'Grid',
                                        Name='TimeSeries',
                                        GridType='Collection',
                                        CollectionType='Temporal')
    # Read grid to get dimension and number of points.
    grid = read_grid_hdf5(gridpath, name)
    dim = len(grid)
    topology_type = '{}DRectMesh'.format(dim)
    geometry_type = 'VXVY' + (dim == 3) * 'VZ'
    components = ('x', 'y', 'z')[:dim]
    gridsize = [len(line) for line in grid]
    number_of_elements = ' '.join(str(n) for n in gridsize[::-1])
    precision = '8'
    # Get time-step indices and time values.
    if states is None:
        states = list(range(nstart, nstart + nt + 1, nsave))
    # Generate the time series.
    for i, state in enumerate(states):
        grid = etree.SubElement(grid_time_series, 'Grid',
                                Name='Grid',
                                GridType='Uniform')
        if times is not None:
            time_value = '{:.6f}'.format(times[i])
        else:
            time_value = '{:0>7}'.format(state)
        time = etree.SubElement(grid, 'Time',
                                Value=time_value)
        topology = etree.SubElement(grid, 'Topology',
                                    TopologyType=topology_type,
                                    NumberOfElements=number_of_elements)
        geometry = etree.SubElement(grid, 'Geometry',
                                    GeometryType=geometry_type)
        # Create XDMF block for the grid. (Use of loop for code-reuse.)
        for component, n in zip(components, gridsize):
            dataitem = etree.SubElement(geometry, 'DataItem',
                                        Dimensions=str(n),
                                        NumberType='Float',
                                        Precision=precision,
                                        Format='HDF')
            dataitem.text = ':/'.join([str(gridpath), name + '/' + component])
        # Create XDMF block for the scalar field variable.
        attribute = etree.SubElement(grid, 'Attribute',
                                     Name=name,
                                     AttributeType='Scalar',
                                     Center='Node')
        dataitem = etree.SubElement(attribute, 'DataItem',
                                    Dimensions=number_of_elements,
                                    NumberType='Float',
                                    Precision=precision,
                                    Format='HDF')
        filepath = datadir / '{:0>7}.h5'.format(state)
        dataitem.text = ':/'.join([str(filepath), name])
    # Write XDMF file.
    tree = etree.ElementTree(xdmf)
    tree.write(str(outpath), pretty_print=True, xml_declaration=True)
    return


def write_xdmf_multi(outpath, config,
                     nstart=None, nt=None, nsave=None,
                     states=None, times=None):
    """Write a XDMF file to read the solution of multiple PetIBM variables.

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

    """
    # Initialize XDMF file.
    xdmf = etree.Element('Xdmf', Version='2.2')
    info = etree.SubElement(xdmf, 'Information',
                            Name='MetaData',
                            Value='ID-23454')
    domain = etree.SubElement(xdmf, 'Domain')
    grid_time_series = etree.SubElement(domain, 'Grid',
                                        Name='TimeSeries',
                                        GridType='Collection',
                                        CollectionType='Temporal')
    # Read grid to get dimension and number of points.
    master_name = list(config['data'].keys())[0]
    gridpath = config['grid']
    grid = read_grid_hdf5(gridpath, master_name)
    dim = len(grid)
    topology_type = '{}DRectMesh'.format(dim)
    geometry_type = 'VXVY' + (dim == 3) * 'VZ'
    components = ('x', 'y', 'z')[:dim]
    gridsize = [len(line) for line in grid]
    number_of_elements = ' '.join(str(n) for n in gridsize[::-1])
    precision = '8'
    # Get time-step indices and time values.
    if states is None:
        states = list(range(nstart, nstart + nt + 1, nsave))
    # Generate the time series.
    for i, state in enumerate(states):
        grid = etree.SubElement(grid_time_series, 'Grid',
                                Name='Grid',
                                GridType='Uniform')
        if times is not None:
            time_value = '{:.6f}'.format(times[i])
        else:
            time_value = '{:0>7}'.format(state)
        time = etree.SubElement(grid, 'Time',
                                Value=time_value)
        topology = etree.SubElement(grid, 'Topology',
                                    TopologyType=topology_type,
                                    NumberOfElements=number_of_elements)
        geometry = etree.SubElement(grid, 'Geometry',
                                    GeometryType=geometry_type)
        # Create XDMF block for the grid. (Use of loop for code-reuse.)
        for component, n in zip(components, gridsize):
            dataitem = etree.SubElement(geometry, 'DataItem',
                                        Dimensions=str(n),
                                        NumberType='Float',
                                        Precision=precision,
                                        Format='HDF')
            dataitem.text = ':/'.join([str(gridpath),
                                       master_name + '/' + component])
        # Create XDMF block for each scalar field variable.
        for name, datadir in config['data'].items():
            attribute = etree.SubElement(grid, 'Attribute',
                                         Name=name,
                                         AttributeType='Scalar',
                                         Center='Node')
            dataitem = etree.SubElement(attribute, 'DataItem',
                                        Dimensions=number_of_elements,
                                        NumberType='Float',
                                        Precision=precision,
                                        Format='HDF')
            filepath = datadir / '{:0>7}.h5'.format(state)
            dataitem.text = ':/'.join([str(filepath), name])
    # Write XDMF file.
    tree = etree.ElementTree(xdmf)
    tree.write(str(outpath), pretty_print=True, xml_declaration=True)
    return
