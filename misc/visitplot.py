"""Module to generate plots of the flow field with VisIt.

Requirements: Python 2, environment variables VISIT_DIR and VISIT_ARCH.
"""

import os
import sys
import yaml

VISIT_DIR = os.environ.get('VISIT_DIR', '')
VISIT_ARCH = os.environ.get('VISIT_ARCH', '')
pkgs_dir = os.path.join(VISIT_DIR, VISIT_ARCH, 'lib', 'site-packages')
if not os.path.exists(pkgs_dir):
    raise ValueError('Set env variables VISIT_DIR and VISIT_ARCH')
sys.path.append(pkgs_dir)
import visit


def visit_initialize():
    visit.LaunchNowin()
    visit_check_version(visit.Version())


def visit_finalize():
    os.remove('visitlog.py')
    visit.CloseComputeEngine()
    visit.Close()


def visit_check_version(version):
    """Check the version of VisIt.

    Parameters
    ----------
    version : string
        VisIt's version.

    """
    # Check version of VisIt.
    script_version = '2.12.1'
    if version != script_version:
        print('[warning] You are using VisIt-{}.'.format(version))
        print('[warning] This script was created with VisIt-{}.'
              .format(script_version))
        print('[warning] It may not work as expected')
    return


def visit_get_view(filepath, dim):
    """Create a view in VisIt.

    Parameters
    ----------
    filepath : string
        Path of the YAML file with the view configuration.
    dim : string
        Either '2D' or '3D'.

    Returns
    -------
    ViewAtts : visit.View2DAttributes or visit.View3DAttributes
        The view attributes.

    """
    dim = dim.upper()
    node = 'View{}Atts'.format(dim)
    with open(str(filepath), 'r') as infile:
        config = yaml.load(infile, Loader=yaml.FullLoader)[node]
    # Set attributes of the view.
    ViewAtts = (visit.View3DAttributes() if dim == '3D'
                else visit.View2DAttributes())
    for key, value in config.items():
        if type(value) is list:
            value = tuple(value)
        setattr(ViewAtts, key, value)
    return ViewAtts


def visit_get_states(state=None, states=None, states_range=[0, None, 1]):
    if state is not None:
        states = [state]
    elif states is None:
        if states_range[1] is None:
            states_range[1] = visit.TimeSliderGetNStates()
        else:
            states_range[1] += 1
        states = range(*states_range)
    return states


def visit_render_save_states(states, out_dir=os.getcwd(), out_prefix=None,
                             figsize=(1024, 1024)):
    visit.Source(os.path.join(VISIT_DIR, VISIT_ARCH, 'bin', 'makemovie.py'))
    visit.ToggleCameraViewMode()

    # Create output directory if necessary.
    if not os.path.isdir(str(out_dir)):
        os.makedirs(str(out_dir))

    # Loop over the states to render and save the plots.
    for i, state in enumerate(states):
        print('[state {}] Rendering and saving figure ...'.format(state))
        visit.SetTimeSliderState(state)

        if i == 0:
            visit.DrawPlots()

        RenderingAtts = visit.RenderingAttributes()
        visit.SetRenderingAttributes(RenderingAtts)

        SaveWindowAtts = visit.SaveWindowAttributes()
        SaveWindowAtts.outputToCurrentDirectory = 0
        SaveWindowAtts.outputDirectory = str(out_dir)
        timestep = int(visit.Query('Time').split(' ')[-1][:-1])
        SaveWindowAtts.fileName = '{}{:0>7}'.format(out_prefix, timestep)
        SaveWindowAtts.family = 0
        SaveWindowAtts.format = SaveWindowAtts.PNG
        SaveWindowAtts.width = figsize[0]
        SaveWindowAtts.height = figsize[1]
        SaveWindowAtts.quality = 100
        SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint
        visit.SetSaveWindowAttributes(SaveWindowAtts)

        visit.SaveWindow()


def visit_plot_pseudocolor_2d(xdmf_path, name,
                              value_range=(-5.0, 5.0),
                              curve2d_paths=None,
                              config_view=None,
                              out_dir=os.getcwd(), out_prefix=None,
                              figsize=(1024, 1024),
                              state=None, states=None,
                              states_range=[0, None, 1]):
    visit_initialize()

    # Create database correlation with optional Curve2D files.
    num_bodies = 0
    databases = [str(xdmf_path)]
    if curve2d_paths is not None:
        num_bodies = len(curve2d_paths)
        databases = [str(path) for path in curve2d_paths]
        databases.append(str(xdmf_path))
    visit.CreateDatabaseCorrelation('common', databases[num_bodies:], 0)

    # Open the file with the coordinates of the immersed boundary.
    if num_bodies > 0:
        for i in range(num_bodies):
            visit.OpenDatabase(databases[i], 0)
            # Add plot the mesh points.
            visit.AddPlot('Curve', 'curve', 1, 1)
            # Set attributes of the curve.
            CurveAtts = visit.CurveAttributes()
            CurveAtts.lineWidth = 1
            CurveAtts.curveColorSource = CurveAtts.Custom
            CurveAtts.curveColor = (0, 0, 0, 255)
            CurveAtts.showLegend = 0
            CurveAtts.showLabels = 0
            visit.SetPlotOptions(CurveAtts)

    # Open the XMF file for the spanwise-averaged z-component of the vorticity.
    visit.OpenDatabase(databases[-1], 0)
    # Add a pseudocolor plot of the scalar field.
    visit.AddPlot('Pseudocolor', name, 1, 1)
    # Set attributes of the pseudocolor.
    PseudocolorAtts = visit.PseudocolorAttributes()
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = value_range[0]
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = value_range[1]
    PseudocolorAtts.colorTableName = 'viridis'
    visit.SetPlotOptions(PseudocolorAtts)

    # Parse the 2D view configuration file.
    if config_view is not None:
        View2DAtts = visit_get_view(config_view, '2D')
    visit.SetView2D(View2DAtts)

    # Remove time and user info.
    AnnotationAtts = visit.AnnotationAttributes()
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.timeInfoFlag = 1
    visit.SetAnnotationAttributes(AnnotationAtts)

    visit.SetActiveWindow(1)

    states = visit_get_states(state=state, states=states,
                              states_range=states_range)
    if out_prefix is None:
        out_prefix = name + '_'
    visit_render_save_states(states, out_dir=out_dir, out_prefix=out_prefix,
                             figsize=figsize)

    visit_finalize()
    return


def visit_plot_contour_3d(xdmf_path, name,
                          value_range=(-5.0, 5.0),
                          p3d_paths=None,
                          config_view=None,
                          out_dir=os.getcwd(), out_prefix=None,
                          figsize=(1024, 1024),
                          state=None, states=None,
                          states_range=[0, None, 1]):
    visit_initialize()

    # Create database correlation with optional Point3D files.
    num_bodies = 0
    databases = [str(xdmf_path)]
    if p3d_paths is not None:
        num_bodies = len(p3d_paths)
        databases = [str(path) for path in p3d_paths]
        databases.append(str(xdmf_path))
    visit.CreateDatabaseCorrelation('common', databases[num_bodies:], 0)

    # Open the file with the coordinates of the immersed boundary.
    if num_bodies > 0:
        for i in range(num_bodies):
            visit.OpenDatabase(databases[i], 0, 'Point3D_1.0')
            # Add plot the mesh points.
            visit.AddPlot('Mesh', 'points', 1, 1)
            # Set attributes of the mesh plot.
            MeshAtts = visit.MeshAttributes()
            MeshAtts.legendFlag = 0
            MeshAtts.meshColor = (255, 204, 0, 1.0 * 255)
            MeshAtts.meshColorSource = MeshAtts.MeshCustom
            MeshAtts.pointSize = 0.05
            MeshAtts.pointType = MeshAtts.Point
            MeshAtts.pointSizePixels = 2
            MeshAtts.opacity = 1
            visit.SetPlotOptions(MeshAtts)

    # Open the XMF file for the z-component of the vorticity.
    visit.OpenDatabase(databases[-1], 0)
    # Add the plot of the contour of the z-component of the vorticity.
    visit.AddPlot('Contour', name, 1, 1)
    # Set attributes of the contour.
    ContourAtts = visit.ContourAttributes()
    ContourAtts.contourNLevels = 2
    ContourAtts.SetMultiColor(0, (0, 51, 102, 0.6 * 255))
    ContourAtts.SetMultiColor(1, (255, 0, 0, 0.6 * 255))
    ContourAtts.legendFlag = 1
    ContourAtts.minFlag = 1
    ContourAtts.maxFlag = 1
    ContourAtts.min = value_range[0]
    ContourAtts.max = value_range[1]
    visit.SetPlotOptions(ContourAtts)

    # Parse the 3D view configuration file.
    if config_view is not None:
        View3DAtts = visit_get_view(config_view, '3D')
    visit.SetView3D(View3DAtts)

    # Remove time and user info.
    AnnotationAtts = visit.AnnotationAttributes()
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.timeInfoFlag = 0
    AnnotationAtts.axes3D.visible = 0
    AnnotationAtts.axes3D.triadFlag = 1
    AnnotationAtts.axes3D.bboxFlag = 0
    visit.SetAnnotationAttributes(AnnotationAtts)

    visit.SetActiveWindow(1)

    states = visit_get_states(state=state, states=states,
                              states_range=states_range)
    if out_prefix is None:
        out_prefix = name + '_'
    visit_render_save_states(states, out_dir=out_dir, out_prefix=out_prefix,
                             figsize=figsize)

    visit_finalize()
    return


def visit_plot_qcrit_wx_3d(xdmf_dir,
                           wx_range=(-5.0, 5.0),
                           q_value=0.1,
                           config_view=None,
                           out_dir=os.getcwd(), out_prefix='qcrit_wx_',
                           figsize=(1024, 1024),
                           state=None, states=None,
                           states_range=[0, None, 1]):
    visit_initialize()

    # Define some variables to get the q_crit and wx_cc.
    p_xdmf_path = os.path.join(str(xdmf_dir), 'p.xmf')
    visit.OpenDatabase(p_xdmf_path, 0)
    # Define cell-centered velocity vector field.
    ux_xdmf_path = os.path.join(str(xdmf_dir), 'u.xmf')
    uy_xdmf_path = os.path.join(str(xdmf_dir), 'v.xmf')
    uz_xdmf_path = os.path.join(str(xdmf_dir), 'w.xmf')
    vel_expr = ('{' +
                'pos_cmfe(<{}[0]id:u>, <p Grid>, 1.0),'.format(ux_xdmf_path) +
                'pos_cmfe(<{}[0]id:v>, <p Grid>, 0.0),'.format(uy_xdmf_path) +
                'pos_cmfe(<{}[0]id:w>, <p Grid>, 0.0)'.format(uz_xdmf_path) +
                '}')
    visit.DefineVectorExpression('velocity', vel_expr)
    # Define Q-criterion.
    qcrit_expr = ('q_criterion(' +
                  'gradient(velocity[0]),' +
                  'gradient(velocity[1]),' +
                  'gradient(velocity[2])' +
                  ')')
    visit.DefineScalarExpression('q_crit', qcrit_expr)
    # Define cell-centered streamwise vorticity.
    wx_xdmf_path = os.path.join(str(xdmf_dir), 'wx.xmf')
    wx_exp = 'pos_cmfe(<{}[0]id:wx>, <p Grid>, 0.0)'.format(wx_xdmf_path)
    visit.DefineScalarExpression('wx_cc', wx_exp)

    # Add a pseudocolor of the cell-centered streamwise vorticity.
    visit.AddPlot('Pseudocolor', 'wx_cc', 1, 1)
    PseudocolorAtts = visit.PseudocolorAttributes()
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = wx_range[0]
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = wx_range[1]
    PseudocolorAtts.colorTableName = 'RdBu'
    PseudocolorAtts.invertColorTable = 1
    PseudocolorAtts.opacityType = PseudocolorAtts.Constant
    PseudocolorAtts.opacity = 0.8
    PseudocolorAtts.legendFlag = 0
    visit.SetPlotOptions(PseudocolorAtts)

    # Add an isosurface of the Q-criterion.
    visit.AddOperator('Isosurface', 1)
    IsosurfaceAtts = visit.IsosurfaceAttributes()
    IsosurfaceAtts.variable = 'q_crit'
    IsosurfaceAtts.contourMethod = IsosurfaceAtts.Value
    IsosurfaceAtts.contourValue = (q_value)
    IsosurfaceAtts.scaling = IsosurfaceAtts.Linear
    visit.SetOperatorOptions(IsosurfaceAtts, 1)

    # Remove info about user, time, database, and legend.
    AnnotationAtts = visit.AnnotationAttributes()
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.timeInfoFlag = 0
    AnnotationAtts.legendInfoFlag = 0
    AnnotationAtts.axes3D.visible = 0
    AnnotationAtts.axes3D.triadFlag = 1
    AnnotationAtts.axes3D.bboxFlag = 0
    visit.SetAnnotationAttributes(AnnotationAtts)

    # Parse the 3D view configuration file.
    if config_view is not None:
        View3DAtts = visit_get_view(config_view, '3D')
    visit.SetView3D(View3DAtts)

    visit.SetActiveWindow(1)

    states = visit_get_states(state=state, states=states,
                              states_range=states_range)
    visit_render_save_states(states, out_dir=out_dir, out_prefix=out_prefix,
                             figsize=figsize)

    visit_finalize()
    return


def visit_plot_qcrit_wx_3d_direct(xdmf_path,
                                  wx_range=(-5.0, 5.0),
                                  q_value=0.1,
                                  config_view=None,
                                  out_dir=os.getcwd(), out_prefix='qcrit_wx_',
                                  figsize=(1024, 1024),
                                  state=None, states=None,
                                  states_range=[0, None, 1]):
    visit_initialize()

    visit.OpenDatabase(str(xdmf_path), 0)

    # Add a pseudocolor of the cell-centered streamwise vorticity.
    visit.AddPlot('Pseudocolor', 'wx_cc', 1, 1)
    PseudocolorAtts = visit.PseudocolorAttributes()
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = wx_range[0]
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = wx_range[1]
    PseudocolorAtts.colorTableName = 'RdBu'
    PseudocolorAtts.invertColorTable = 1
    PseudocolorAtts.opacityType = PseudocolorAtts.Constant
    PseudocolorAtts.opacity = 0.8
    PseudocolorAtts.legendFlag = 0
    visit.SetPlotOptions(PseudocolorAtts)

    # Add an isosurface of the Q-criterion.
    visit.AddOperator('Isosurface', 1)
    IsosurfaceAtts = visit.IsosurfaceAttributes()
    IsosurfaceAtts.variable = 'qcrit'
    IsosurfaceAtts.contourMethod = IsosurfaceAtts.Value
    IsosurfaceAtts.contourValue = (q_value)
    IsosurfaceAtts.scaling = IsosurfaceAtts.Linear
    visit.SetOperatorOptions(IsosurfaceAtts, 1)

    # Remove info about user, time, database, and legend.
    AnnotationAtts = visit.AnnotationAttributes()
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.timeInfoFlag = 0
    AnnotationAtts.legendInfoFlag = 0
    AnnotationAtts.axes3D.visible = 0
    AnnotationAtts.axes3D.triadFlag = 1
    AnnotationAtts.axes3D.bboxFlag = 0
    visit.SetAnnotationAttributes(AnnotationAtts)

    # Parse the 3D view configuration file.
    if config_view is not None:
        View3DAtts = visit_get_view(config_view, '3D')
    visit.SetView3D(View3DAtts)

    visit.SetActiveWindow(1)

    states = visit_get_states(state=state, states=states,
                              states_range=states_range)
    visit_render_save_states(states, out_dir=out_dir, out_prefix=out_prefix,
                             figsize=figsize)

    visit_finalize()
    return
