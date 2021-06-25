import fmm3dbie as bie
import fmm3dpy
import numpy as np
import mwaspbie as mw

# Define geometry and number of components
string1 = '../geometries/lens_r00.go3?../geometries/con_r00.go3?'
string2 = '../geometries/lens_r00.msh?../geometries/con_r00.msh?'
n_components = 2

# compute number of patches and points
npatches,npts = mw.em_solver_wrap_mem(string1,n_components)

# Set translation and scaling of each component
dP = np.zeros((4,n_components),order="F")
dP[3,:] = 1.0

eps = 1e-6

# Get geometry info
[npatches_vect,npts_vect,norders,ixyzs,iptype,srcvals,srccoefs,wts,
sorted_vector,exposed_surfaces] = mw.em_solver_open_geom(string1,dP,npatches,
  npts,eps)



#Estimate number of vertices and flat triangles in geometry
nverts,nel = mw.em_gen_plot_info_surf_mem(string2,n_components)

#Get the vertices and definition of the flat triangles in 
#the geometry
verts,elements = mw.em_gen_plot_info_surf(string2,n_components,nverts,nel)

# Transpose the list of elements, i.e. elements currently
# stores which vertices form the flat triangles, the transpose
# stores which elements meet at a given vertex and are stored
# in a sparse compressed format
#
# iver_el_list: stores a continuous list of element ids
# iverstart: iverstart[i] indicates the location in the iver_el_list array
#   where the list of elements for vertex i begin
# iverind: denotes the vertex number in element iver_el_list[j] 
#  corresponding to the vertex i of the vertex-element pair
#  [i,iver_el_list[j]]

iver_el_list,iverstart,iverind = mw.em_elem_trans(nverts,elements)

fvals = srcvals[0:3,:]
rscales = np.ones(nel)

fvalsout = mw.em_surf_fun_to_plot_fun(norders,ixyzs,iptype,fvals,iver_el_list,
  iverstart,iverind,rscales)

erra = np.linalg.norm(fvalsout-verts)/np.linalg.norm(verts)
print("error in estimating vertex locations="+str(erra))


