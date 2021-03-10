.. _pyt:

Python
=======

The Python interface has four callable subroutines:

*  `Memory estimator <python.html#mem>`__: Estimate number of points and patches in given geometry
*  `Solver <python.html#solver>`__: Solve the dielectric problem for given geometry 
*  `Post processor <python.html#postproc>`__: Using known solutions for the surface currents, compute the electric and magnetic fields at an arbitrary collection of points
*  `Geometry extractor <python.html#geom>`__: Extract the geometry information for using various routines from fmm3dbie
*  `Plotting routines <python.html#plot>`__: Plotting routines for generating vtk files


.. _mem:

Memory estimator
*******************

This subroutine estimates the total number of patches and points 
for the geometry. Recall that the geometry is defined as a string of go3 
files where each component is separted by a ``?``. For example, 
for a two component geometry with components ``a.go3`` and ``b.go3``, the
string should be set to ``a.go3?b.go3?``.

.. code:: python

   npatches,npts = em_solver_wrap_mem(s,ncomp)


Args:

-  s: string
      string defining the geometry
-  ncomp: integer
      number of components in the geometry

Returns:

-  npatches: integer
      total number of patches
-  npts: integer
      total number of discretization points

.. _solver:

Solver
*******************

This subroutine solves for the surface electric and magnetic currents.
In the muller representation, the electric field and magnetic field in
region $j$ denoted by $E_{j}$, $H_{j}$ respectively are represented as

.. math::

   E_{j}(x) &= \mu_{j} \nabla \times S_{k_{j},\Gamma_{j}}[a_{j}] - \frac{\nabla \times \nabla \times S_{k_{j},\Gamma_{j}}[b_{j}]}{i\omega} \\
   H_{j}(x) &= \varepsilon_{j} \nabla \times S_{k_{j},\Gamma_{j}}[b_{j}] + \frac{\nabla \times \nabla \times S_{k_{j},\Gamma_{j}}[a_{j}]}{i\omega}

where $\Gamma_{j}$ is the boundary of the region $j$, $a_{j},b_{j}$ are
the electric and magnetic currents on $\Gamma_{j}$, $\varepsilon_{j}$ 
is the electric permitivity in region $j$, and $\mu_{j}$ the magnetic
permeability, $k_{j} = \omega \sqrt{\varepsilon_{j}} \sqrt{\mu_{j}}$ is
the wave number in region $j$, and $S_{k,\Gamma}[\sigma]$ is the
single layer potential on the boundary $\Gamma$ with density $\sigma$,
and wave number $k$ given by

.. math::

   S_{k,\Gamma}[\sigma] = \int_{\Gamma} \frac{e^{ik|x-y|}}{|x-y|} \sigma(y) \, da (y)

The surface currents $a_{j},b_{j}$ are stored as their projections onto
the tangent vectors.

Note that the incident field can either be a plane wave, or an analytic
solution obtained by constructing a plane wave in interior regions and a
point electric and magnetic dipole in the exterior. Given a incident
direction of plane wave on $S^{2}$ denoted by $(\theta_{d},\phi_{d})$ and polarization vectors
$v_{\theta}$, $v_{\phi}$, the incident electric and magnetic fields are
given by

.. math::

   E_{j} &= \left(\hat{\theta} v_{\theta} + \hat{\phi} v_{\phi} \right)
   e^{ik x \cdot \hat{r}} \\
   H_{j} &= \frac{\sqrt{\mu_{j}}}{\sqrt{\varepsilon_{j}}}\left( -\hat{\theta} v_{\phi} + \hat{\phi} v_{\theta}  \right)
   e^{ik x \cdot \hat{r}} \, ,


where

.. math::

   \hat{r} = \begin{bmatrix} 
    \sin{\theta_{d}} \cos{\phi_{d}} \\
    \sin{\theta_{d}} \sin{\phi_{d}} \\
    \cos{\theta_{d}}
   \end{bmatrix} \, , \quad 
   \hat{\theta} = \begin{bmatrix}
   \cos{\theta_{d}} \cos{\phi_{d}} \\
   \cos{\theta_{d}} \sin{\phi_{d}} \\
   -\sin{\theta_{d}}
   \end{bmatrix} \, , \quad
   \hat{\phi} = \begin{bmatrix}
   -\sin{\phi_{d}} \\
   \cos{\phi_{d}} \\
   0
   \end{bmatrix} \, .


 
.. code:: python

   soln, err_est = em_solver_wrap(s,dP,contrast_matrix,npts,omega,icase,direction,pol,eps,eps_gmres)

Args:

-  s: character(len=*)
      string defining the go3 files
-  dP: double precision (4,n_components)
      | dP(1:3,i) translation for the ith component
      | dP(4,i) scaling
-  contrast_matrix: double complex (4,n_components)
      $\varepsilon$,$\mu$ in positive and negative normal direction of each component

         * contrast_matrix(1,i) = $\varepsilon$ on the positive normal side of $\Gamma_{i}$
         * contrast_matrix(2,i) = $\mu$ on the positive normal side of $\Gamma_{i}$
         * contrast_matrix(3,i) = $\varepsilon$ on the negative normal side of $\Gamma_{i}$
         * contrast_matrix(4,i) = $\mu$ on the negative normal side of $\Gamma_{i}$

-  omega: double precision
      wavenumber of problem (should be consistent with the units of the go3 file)
-  icase: integer
      solve either for analytic solution or incident plane wave
      (analytic solution in the interior regions is still the plane
      wave)
        
      *  icase = 1, analytic solution test. Exterior solution
         is constructed via dipole at centroid of geometry
         whose orientation is in the (1,1,1) direction
         and renormalized to size/omega (this makes sure exterior and
         interior solutions are commensurate)
      * icase = 2, incident plane wave computation
-  direction: double precision(2)

      * direction(1) - $\phi_d$ angle of the plane wave
      * direction(2) - $\theta_d$ angle of the plane wave
-  pol: double complex(2)

      * pol(1) - $v_{\phi}$, $\phi$ component of the plane wave
      * pol(2) - $v_{\theta}$, $\theta$ component of plane wave

-  eps: double precision
      tolerance of FMM and quadrature generation
-  eps_gmres: double precision
      relative residual tolerance for GMRES
    
Returns:

-  soln: double complex (4$\cdot$npts)
     
      The projection of surface currents $a$,$b$ onto tangent vectors

      *  soln(1:npts) is the projection of $a$ onto the first tangent vector
      *  soln(npts+1:2$\cdot$npts) is the projection of $a$ onto the second tangent vector
      *  soln(2$\cdot$npts:3$\cdot$npts) is the projection of $b$ onto the first tangent vector
      *  soln(3$\cdot$npts+1:4$\cdot$npts) is the projection of $b$ onto the second tangent vector

-  err_est: integer
      Estimated error. Only meaningful if icase = 1

.. _postproc:

Post processor
*******************

Given the computed surface currents $a$, $b$, this subroutine evaluates
the electric field and magnetic fields at a given collection of targets.


.. code:: python

   E, H = em_solver_wrap_postproc(s,dP,contrast_matrix,omega,eps,soln,targs)

Args:

-  s: character(len=*)
      string defining the go3 files
-  dP: double precision (4,n_components)
      | dP(1:3,i) translation for the ith component
      | dP(4,i) scaling
-  contrast_matrix: double complex (4,n_components)
      $\varepsilon$,$\mu$ in positive and negative normal direction of each component

         * contrast_matrix(1,i) = $\varepsilon$ on the positive normal side of $\Gamma_{i}$
         * contrast_matrix(2,i) = $\mu$ on the positive normal side of $\Gamma_{i}$
         * contrast_matrix(3,i) = $\varepsilon$ on the negative normal side of $\Gamma_{i}$
         * contrast_matrix(4,i) = $\mu$ on the negative normal side of $\Gamma_{i}$

-  omega: double precision
      wavenumber of problem (should be consistent with the units of the go3 file)
-  eps: double precision
      tolerance of FMM and quadrature generation

-  soln: double complex (4$\cdot$npts)
     
      The projection of surface currents $a$,$b$ onto tangent vectors

      *  soln(1:npts) is the projection of $a$ onto the first tangent vector
      *  soln(npts+1:2$\cdot$npts) is the projection of $a$ onto the second tangent vector
      *  soln(2$\cdot$npts:3$\cdot$npts) is the projection of $b$ onto the first tangent vector
      *  soln(3$\cdot$npts+1:4$\cdot$npts) is the projection of $b$ onto the second tangent vector

-  targs: double precision (3,ntarg)
      x,y,z coordinates of target locations where electric and magnetic
      fields are to be computed
    
Returns:

-  E: double complex (3,ntarg)
      x,y,z components of the electric field at the target locations
      
-  H: double complex (3,ntarg)
      x,y,z components of the magnetic field at the target locations
      

.. _geom:

Geometry extractor
*******************

This subroutine extracts the geometry information which can be later
used for various plotting routines.

.. code:: python

   npatches_v,npts_v,norders,ixyzs,iptype,srcvals,srccoefs,sorted_vector,exposed_surfaces = em_solver_open_geom(s,dP,npatches,npts,eps)
   
This subroutine extracts various geometric information which can be
utilized for the visualization and other fmm3dbie routines.  

Args:

-  s: character(len=*)
      string defining the go3 files
-  dP: double precision (4,n_components)
      | dP(1:3,i) translation for the ith component
      | dP(4,i) scaling
-  npatches: integer
      total number of patches 
-  npts: integer
      total number of discretization points 
-  eps: double precision
      tolerance of FMM and quadrature generation

Returns:

-  npatches_v: integer(ncomp)
     number of patches on each component

-  npts_v: integer(ncomp)
     number of discretization points on each component

-  norders: integer(npatches)
      order of discretization of each patch

-  ixyzs: integer(npatches+1)
      location in srcvals, and srccoefs where information for patches begin 

-  iptype: integer(npatches)
      patch type

-  srcvals: double precision(12,npts)
      surface samples of $\boldsymbol{x}^{j}, \partial_{u} \boldsymbol{x}^{j}, \partial_{v} \boldsymbol{x}^{j},$
      and $\boldsymbol{n}^{j}$, where
    
    .. math::
   
        \boldsymbol{n}^{j} = \frac{\partial_{u} \boldsymbol{x}^{j} \times 
        \partial_{v} \boldsymbol{x}^{j}}{|\partial_{u} \boldsymbol{x}^{j}
        \times \partial_{v} \boldsymbol{x}^{j}|}

-  srccoefs: double precision(9,npts)
      Orthogonal polynomial expansions of 
      $\boldsymbol{x}^{j}, \partial_{u} \boldsymbol{x}^{j}$, 
      and $\partial_{v} \boldsymbol{x}_{j}$

-  wts: double precision(npts)
      quadrature weights for integrating smooth functions

-  sorted_vector: integer (ncomp+1)
      sorting of components which respects the partial ordering
      induced by set inclusion

-  exposed_surfaces: logical(ncomp)
      true is the component shares boundary with unbounded component,
      false otherwise.
      
