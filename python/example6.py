import fmm3dbie
import fmm3dpy
import numpy as np
import mwaspbie as mw
import pickle


# Define various components of the geometry and problem 
# specification
#

# 
# Refinement level of geometry used
# iref should be 0 or 1
iref = 0

#
#  iomegafrac is a paramater
#  that defines the wavelength in the
#  simulation
#  
#
#
iomegafrac = 2

if(iomegafrac == 1):
    omega = np.pi*2.0/600.0

if(iomegafrac == 2):
    omega = np.pi*2.0/550.0

if(iomegafrac == 3):
    omega = np.pi*2.0/450.0

if(iomegafrac == 4):
    omega = np.pi*2.0/350.0


#
# icase = 2. Denotes the fact that 
# we are solving a scattering problem 
# and not an analytic test
#
icase = 2

#
# idir. Flag for determining 
# what the direction of the incident field was
#
# currently supported value idir = 1
# for which the incident field is along 
# smallest principal component of the rhabdom for which
# direction[0:1] = []
# 
#
#
idir = 1
direction = np.zeros(2)
if (idir == 1):
    direction[0] = -1.530844785708483 
    direction[1] = 1.063728846738649 


#
#
#  ipol. Flag for determining the in plane
#  polarization directions.
#  Options for ipol = 1 or 2. 
#
#
ipol = 1
pol = np.zeros(2,dtype='complex')
if (ipol == 1):
    pol[0] = 1.0

if (ipol == 2):
    pol[1] = 1.0



# Define geometry and number of components
string1 = '../geometries/lens_r0'+str(iref)+'.go3?'
string1 = string1+'../geometries/con_r0'+str(iref)+'.go3?'
string1 = string1+'../geometries/rhabdom_r0'+str(iref)+'.go3?'

string2 = '../geometries/lens_r0'+str(iref)+'.msh?'
string2 = string1+'../geometries/con_r0'+str(iref)+'.msh?'
string2 = string1+'../geometries/rhabdom_r0'+str(iref)+'.msh?'



n_components = 3



fsol = '../test/data/soln_iref'+str(iref)+'_iomega'+str(iomegafrac)
fsol = fsol+'_icase'+str(icase)+'_idir'+str(idir)+'_ipol'+str(ipol)+'.dat'

ftarg = '../test/data/targ_iref'+str(iref)+'_iomega'+str(iomegafrac)
ftarg = ftarg+'_icase'+str(icase)+'_idir'+str(idir)+'_ipol'+str(ipol)+'.dat'

fpickle = 'pickle_iref'+str(iref)+'_iomega'+str(iomegafrac)
fpickle = fpickle+'_icase'+str(icase)+'_idir'+str(idir)+'_ipol'+str(ipol)+'.pickle'


# compute number of patches and points
npatches,npts = mw.em_solver_wrap_mem(string1,n_components)

# Set translation and scaling of each component
dP = np.zeros((4,n_components),order="F")
dP[3,:] = 1.0

# set wave number of problem, should be consistent with
# the units of prescribed geometry

# set material parameters on either side of each component
contrast_matrix = np.ones((4,n_components),order="F",dtype="complex")
contrast_matrix[2,0] = 1.452
contrast_matrix[2,1] = 1.348
contrast_matrix[2,2] = 1.363

# set tolerance for quadrature
eps = 1.0e-3

# read solution
xsol = np.loadtxt(fsol)
soln = xsol[:,0] + 1j*xsol[:,1]

# Get geometry info
[npatches_vect,npts_vect,norders,ixyzs,iptype,srcvals,srccoefs,wts,
sorted_vector,exposed_surfaces] = mw.em_solver_open_geom(string1,dP,npatches,
  npts,eps)

#
#  Set points per wavelength in each direction 
#  for target grid
#
ppw = 5


#
# create target grid
#

xmin = np.min(srcvals[0,:])
xmax = np.max(srcvals[0,:])

ymin = np.min(srcvals[1,:])
ymax = np.max(srcvals[1,:])

zmin = np.min(srcvals[2,:])
zmax = np.max(srcvals[2,:])

dx = xmax-xmin
dy = ymax-ymin
dz = zmax-zmin

# determine grid spacing to get the correct 
# resolution in each dimension of target grid
nx = int(np.ceil(1.5*(xmax-xmin)*omega/2/np.pi*ppw))
ny = int(np.ceil(1.5*(ymax-ymin)*omega/2/np.pi*ppw))
nz = int(np.ceil(1.5*(zmax-zmin)*omega/2/np.pi*ppw))

xs = np.linspace(xmin-0.25*dx,xmax+0.25*dx,nx)
ys = np.linspace(ymin-0.25*dy,ymax+0.25*dy,ny)
zs = np.linspace(zmin-0.25*dz,zmax+0.25*dz,nz)
xx,yy,zz = np.meshgrid(xs,ys,zs)

nt = nx*ny*nz
targs = np.zeros((3,nt),order="F")
targs[0,:] = xx.reshape(nt)
targs[1,:] = yy.reshape(nt)
targs[2,:] = zz.reshape(nt) 

#
# Compute fields at targets using computed surface currents
#
E,H = mw.em_solver_wrap_postproc(string1,dP,contrast_matrix,omega,eps,soln,
  targs)

write_obj = {'xgrid':xs,'ygrid':ys,'zgrid':zs,'E':E,'H':H}
with open(fpickle,'wb') as f:
    pickle.dump(write_obj,f)

