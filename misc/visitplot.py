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


def visit_plot_qcrit_wx_3d(xdmf_dir,
                           wx_range=(-5.0, 5.0),
                           q_value=0.1,
                           config_view=None,
                           out_dir=os.getcwd(), out_prefix='qcrit_wx_',
                           figsize=(1024, 1024),
                           state=None, states=None,
                           states_range=[0, None, 1]):
    visit.LaunchNowin()
    visit_check_version(visit.Version())

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
        with open(str(config_view), 'r') as infile:
            config_view = yaml.load(infile, Loader=yaml.FullLoader)
            config_view = config_view['View3DAtts']
    # Set attributes of the view.
    View3DAtts = visit.View3DAttributes()
    for key, value in config_view.items():
        if type(value) is list:
            value = tuple(value)
        setattr(View3DAtts, key, value)
    visit.SetView3D(View3DAtts)

    visit.SetActiveWindow(1)

    visit.Source(os.path.join(VISIT_DIR, VISIT_ARCH, 'bin', 'makemovie.py'))
    visit.ToggleCameraViewMode()

    # Create output directory if necessary.
    if not os.path.isdir(str(out_dir)):
        os.makedirs(str(out_dir))

    # Loop over the states to render and save the plots.
    if state is not None:
        states = [state]
    elif states is None:
        if states_range[1] is None:
            states_range[1] = visit.TimeSliderGetNStates()
        else:
            states_range[1] += 1
        states = range(*states_range)

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

    os.remove('visitlog.py')
    visit.CloseComputeEngine()
    visit.Close()
    return
