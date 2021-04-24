import fmm3dbie as bie
import fmm3dpy
import numpy as np
import mwaspbie as mw

# Define geometry and number of components
string1 = '../geometries/lens_r00.go3?'
n_components = 1

# compute number of patches and points
npatches,npts = mw.em_solver_wrap_mem(string1,n_components)

# Set translation and scaling of each component
dP = np.zeros((4,n_components),order="F")
dP[3,0] = 1.0

eps = 1e-6

# Get geometry info
[npatches_vect,npts_vect,norders,ixyzs,iptype,srcvals,srccoefs,wts,
sorted_vector,exposed_surfaces] = mw.em_solver_open_geom(string1,dP,npatches,
  npts,eps)

# plot surface info (plots scalar field as z coordinate of the surface)
bie.surf_vtk_plot(norders,ixyzs,iptype,srccoefs,srcvals,'surf.vtk','a')

# plot surface info with prescribed scalar field. We will use
# x coordinate as a demo which is stored in srcvals[0,:]
xvals = srcvals[0,:]
bie.surf_vtk_plot_scalar(norders,ixyzs,iptype,srccoefs,srcvals,xvals,
  'surfx.vtk','a')

# plot surface info with real vector field. In this example we will
# use the normals which are stored in srcvals[9:12,0:npts]

normals = srcvals[9:12,0:npts]
bie.surf_vtk_plot_vec(norders,ixyzs,iptype,srccoefs,srcvals,normals,
   'surf_normals.vtk','a')


# plot surface info with a complex vector field. In this example
# we will use normals + 1j*[1,0,0]
ztest = np.zeros(np.shape(normals),dtype='complex')
ztest[0,:] = 1j
ztest = ztest + normals
bie.surf_vtk_plot_zvec(norders,ixyzs,iptype,srccoefs,srcvals,ztest,
   'surf_znormals.vtk','a')

# plot surface info with one of the tangential currents
# which is the output soln of the subroutine em_solver_wrap
# In this example we create a dummy solution whose real part is
# the tangential projection of [1,0,0], and whose imaginary part 
# is the tangential projection of [0,0.5,0] for the electric current
# and the magnetic current is the tangential projection of 
# [0,1,0] + 1j[0,0,0.5]
# 
#
# In order to compute the tangential projections, we need
# to obtain the orthonormal frame used in the construction
# of the solution. This can be computed via the subroutine
# orthonormalize_all
#
xu = srcvals[3:6,0:npts]
[ru,rv] = bie.orthonormalize_all(xu,normals)
soln = np.zeros(4*npts,dtype='complex')
soln[0:npts] = ru[0,:] + 1j*ru[1,:]/2
soln[npts:2*npts] = rv[0,:] + 1j*rv[1,:]/2
soln[2*npts:3*npts] = ru[1,:] + 1j*ru[2,:]/2
soln[3*npts:4*npts] = rv[1,:] + 1j*rv[2,:]/2

# first plot electric current
ifjk = 1
mw.em_plot_surf_current_vtk(norders,ixyzs,iptype,srccoefs,srcvals,soln,
  ifjk,'jvals.vtk','a')


# Compare jvals to exact vector field by removing normal projection
# from [1,0,0] + 1j*[0,0.5,0]

jvals_ex = np.zeros((3,npts),dtype='complex')
jvals_ex[0,:] = 1
jvals_ex[1,:] = 0.5j
rdot = jvals_ex[0,:]*normals[0,:] + jvals_ex[1,:]*normals[1,:] + jvals_ex[2,:]*normals[2,:]
jvals_ex[0:3,:] = jvals_ex[0:3,:] - rdot*normals[0:3,:]
bie.surf_vtk_plot_zvec(norders,ixyzs,iptype,srccoefs,srcvals,jvals_ex,
   'jvals_ex.vtk','a')

ifjk = 2
mw.em_plot_surf_current_vtk(norders,ixyzs,iptype,srccoefs,srcvals,soln,
  ifjk,'kvals.vtk','a')

