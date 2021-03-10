.. FMM3D documentation master file, created by
   sphinx-quickstart on Wed Nov  1 16:19:13 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Optical simulators for a miniwasp eye   (mwaspbie)
==================================================

.. image:: mwasp.pdf
    :width: 90%
    :align: center
	    
`miniwaspbie <https://github.com/mrachh/miniwasp>`_ 
is a package for solving time harmonic Maxwell equations dielectric
problems specialized to the geometries arising in the eyes
of miniwasps.
The electric permitivity in each component of the eye is assumed
to be a piecewise constant, and has dependencies on 
`FMM3D <https://github.com/flatironinsitute/FMM3D>`_  and
`fmm3dbie <https://gitlab.com/fastalgorithms/fmm3dbie>`_.
The library is written in Fortran and has wrappers for Python.
As an example, suppose that we are 
given a component of the eye denoted by 
$\Omega_{-} \in \mathbb{R}^{3}$ 
whose boundary is $\Gamma$,
an incoming electric field $E^{\textrm{inc}}$, and incoming 
magnetic field $H^{\textrm{inc}}$.
Let $\Omega_{+} = \mathbb{R}^{3} \setminus \Omega_{-}$, denote
the exterior of $\Omega_{-}$.
Let $E_{-},H_{-}$ denote the electric and magnetic fields inside
$\Omega$, and $E^{+}, H_{+}$ denote the electric and magnetic
fields outisde $\Omega$. Let $\varepsilon_{\pm}$ denote the 
electric permitivities in the interior and exterior domains 
respectively, and $\mu_{\pm}$ denote the magnetic permeabilities
in the interior and exterior domains respectively. 
Let $\omega$ denote the wave number of the incoming wave. Then
the electric and magnetic fields satisfy the interface problem 
for time harmonic
Maxwell equations given by


.. math:: 
   
   \nabla \times E_{\pm} &= i\omega \mu_{\pm} H_{\pm} \quad \mbox{ in } \Omega_{\pm} \, , \\
   \nabla \times H_{pm} &= -i \omega \varepsilon_{\pm} E_{\pm} \quad
   \mbox{ in } \Omega_{\pm} \, , \\
   n \times (E_{+} - E_{-}) &= -n \times E^{\textrm{inc}} \quad \mbox {
   on } \Gamma \, , \\
   n \times (H_{+} - H_{-}) &= -n \times H^{\textrm{inc}} \quad \mbox {
   on } \Gamma \, , \\

The library currently supports high order triangulations of surfaces
stored in the .go3 format. In this setup each map from the standard 
right triangle is stored at order $p$ Vioreanu-Rokhlin nodes. 
The input format currently assumes that each patch is discretized using
the same order nodes. Upcoming support will be provided for
triangulations stored in .gmsh
format, and .step format, and quadrilaterizations in all of the above
formats.

.. toctree::
   :maxdepth: 2
	   
   install
   python
   
   
