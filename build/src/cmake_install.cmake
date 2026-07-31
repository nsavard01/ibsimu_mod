# Install script for directory: /home/nsavard/ibsimu_install/ibsimu_mod/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/opt/ibsimu/mod_ibsimu")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/home/nsavard/ibsimu_install/ibsimu_mod/build/src/libibsimu.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ibsimu" TYPE FILE FILES
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/axisymmetricvectorfield.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/callback.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/ccolmatrix.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/cfifo.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/colormap.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/compmath.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/comptime.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/config.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/constants.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/convergence.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/coordmapper.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/coordmatrix.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/crowmatrix.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/csgobject_solid.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/diag_precond.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/dxf_solid.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/empty_precond.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/epot_bicgstabsolver.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/epot_efield.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/epot_field.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/epot_gssolver.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/epot_matrixsolver.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/epot_mgsolver.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/epot_mgsubsolver.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/epot_solver.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/epot_umfpacksolver.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/eqpotgraph.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/error.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/field.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/fielddiagplot.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/fielddiagplotter.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/fieldgraph.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/file.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/fonts.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/frame.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/func_solid.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/geom3dplot.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/geometry.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/geomplot.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/geomplotter.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/glrenderer.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gmg_precond.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/graph.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/graph3d.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkfielddiagdialog.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkfielddiagexportdialog.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkfielddiagwindow.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkframewindow.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkgeom3dwindow.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkgeomwindow.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkhardcopy.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkparticlediagdialog.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkparticlediagexportdialog.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkparticlediagwindow.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkplotter.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkpreferences.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/gtkwindow.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/hbio.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/histogram.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/ibsimu.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/icons.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/id.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/ilu0_precond.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/ilu1_precond.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/interpolation.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/label.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/legend.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/lineclip.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mat3d.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/matrix.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mesh.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/meshcolormap.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/meshgraph.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/meshscalarfield.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/meshvectorfield.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/multimeshvectorfield.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mvector.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxfarc.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxfblocks.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxfcircle.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxfentities.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxffile.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxffont.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxfheader.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxfinsert.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxfline.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxflwpolyline.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxfmtext.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxfspline.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/mydxftables.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/palette.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/parallel_util.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/particledatabase.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/particledatabaseimp.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/particlediagplot.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/particlediagplotter.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/particlegraph.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/particleiterator.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/particles.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/particlestatistics.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/particlestepper.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/plotter.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/polysolver.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/precond.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/random.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/rbsor_precond.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/readascii.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/renderer.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/ruler.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/scalarfield.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/scharge.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/scheduler.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/softwarerenderer.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/solid.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/solidgraph.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/solver.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/sort.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/statusprint.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/stl_solid.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/stlfile.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/symbols.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/timer.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/trajectory.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/trajectorydiagnostics.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/transformation.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/types.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/vec3d.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/vec4d.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/vectorfield.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/vtriangle.hpp"
    "/home/nsavard/ibsimu_install/ibsimu_mod/src/xygraph.hpp"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ibsimu" TYPE FILE FILES "/home/nsavard/ibsimu_install/ibsimu_mod/build/src/config.h")
endif()

