!   This file has the following user-callable routines
!
!   em_solver_wrap_mem - estimate the number of patches and points
!     in the geometry prescription defined through a collection
!     of go3 files encoded in a string, with each separate
!     component separated by a '?'
!
!   em_solver_open_geom - get the surface discretization info
!     for the geometry defined above
!
!   em_solver_wrap - solve the dielectric problem using the
!   muller representation. 
!    * Geometry to be defined via a string  which encodes a 
!      collection of go3 files
!    * Material parameters must be passed in as a 4,n_components
!      vector, where for each component, \ep, \mu
!      is defined for either side of the normal
!    * Boundary condition is either a plane wave on the exposed
!      surface with a prescribed direction and polarization
!           OR
!      an analytical test where the field in each region
!      is defined via the plane wave, and the field in the
!      exterior region is defined by a point dipole
!      located at the centroid of the first component
!      with orientation vector (1,2,3) scaled by 1e-3
!    * On output, returns the u,v components of the surface
!      currents, and an estimated error if the analytical 
!      test is being done
!
!   em_solver_wrap_postproc - evaluate the solution
!   at an arbitrary collection of targets where the u,v
!   projection of surface currents are know
!    * Geometry as defined above
!    * Note: routine might currently fail if there are too
!      many targets in the volume close to a particular triangle
!      in the discretization. One temporary fix is to
!      change the following lines in
!      (fmm3dbie download folder)/src/quadratures/ggq-quads.f
!
!      line 215, 623: (Set ntrimax to a higher value)
!
!      Make sure to reinstall the library after making the
!      changes
!    
!
!    em_sol_exact - computes the exact solution used in the analytic
!    test in em_solver_wrap. 
!
!
!    em_plot_surf_current_vtk - plots the surface current 
!    along with the mesh


      subroutine em_solver_wrap_mem(string0,n_components,npatches,npts)
!
!f2py intent(in) string0,n_components
!f2py intent(out) npatches,npts
!
!  This subroutine estimates the number of points 
!  for the given string defining the geometry
!
!  Input arguments:
!    - string0: character(len=*)
!      string defining the go3 files
!    - n_components: integer
!      number of components
!    
!  Output arguments:
!    - npatches: integer
!        number of patches
!    - npts: integer
!        number of points in the discretization
!
      implicit none
      character (len=*), intent(in) :: string0
      character (len=2000) :: string1
      integer ll
      integer, intent(in) :: n_components
      integer, intent(out) :: npts,npatches
      integer, allocatable :: iwords(:)
      integer i,n0,n1
      character (len=1000) :: fname1
      
      ll = len(string0)
      string1(1:ll) = string0
      allocate(iwords(n_components+1))
      call text_process(string1,n_components,iwords)
      
      npts = 0
      npatches = 0
      do i=1,n_components
        fname1 = trim(string1(iwords(i)+1:iwords(i+1)-1))
        n0 = 0
        n1 = 0
        call open_gov3_geometry_mem(fname1,n1,n0)
        npts = npts + n0
        npatches = npatches + n1
      enddo

      return
      end subroutine em_solver_wrap_mem
!
!
!
!
!
      subroutine em_solver_open_geom(string0,n_components,dP, &
        npatches,npts,eps,npatches_vect,npts_vect,norders,ixyzs,iptype, &
        srcvals,srccoefs,wts,sorted_vector,exposed_surfaces)
!
!f2py intent(in) string0,n_componnts,dP,npatches,npts,eps
!f2py intent(out) npatches_vect,npts_vect,norders,ixyzs,iptype
!f2py intent(out) srcvals,srccoefs,wts,sorted_vector,exposed_surfaces
!
!  This subroutine takes in the string format for a collection
!  of go3 files and returns relevant surface information
!  to be subsequently used by the solver and other routines in
!  fmm3dbie
!
!  Input arguments:
!    - string0: character(len=*)
!      string defining the go3 files
!    - n_components: integer
!      number of components
!    - dP: real *8 (4,n_components)
!        dP(1:3,i) translation for the ith component
!        dP(4,i) scaling
!    - npatches: integer
!        total number of patches
!    - npts: integer
!        total number of discretization points on surface
!    - eps: real *8
!        precision to be used for topological sorting
!   
!    Output arguments:
!      - npatches_vect: integer(n_components)
!          number of patches on each component of the geometry
!      - npts_vect: integer(n_components)
!          number of discretization points on each component
!      - norders: integer(npatches)
!          order of discretization of each patch
!      - ixyzs: integer(npatches+1)
!          starting location of points on patch i
!      - iptype: integer(npatches)
!          type of patch
!          iptype = 1, triangle discretized using RV nodes
!      - npts: integer
!          total number of points on the surface
!      - srccoefs: double precision (9,npts)
!          koornwinder expansion coefs of geometry info
!      - srcvals: double precision (12,npts)
!          xyz, dxyz/du,dxyz/dv, normals at all nodes
!      - sorted_vector: integer(n_components+1)
!          sorting of components to respect partial ordering
!          of inclusion
!       - exposed_surfaces: logical(n_components)
!           true if surface is exposed to ambient space
!
      implicit real *8 (a-h,o-z) 
      character (len=*), intent(in) :: string0
      integer, intent(in) :: n_components
      real *8, intent(in) :: dP(4,n_components),eps
      integer, intent(in) :: npts,npatches

      integer, intent(out) :: npatches_vect(n_components)
      integer, intent(out) :: npts_vect(n_components)
      integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(out) :: iptype(npatches)
      real *8, intent(out) :: srcvals(12,npts),srccoefs(9,npts),wts(npts)
      integer, intent(out) :: sorted_vector(n_components+1)
      logical, intent(out) :: exposed_surfaces(n_components)
      integer, allocatable :: iwords(:)

      character *2000, string1
      character *1000, fname1

      allocate(iwords(n_components+1))

      ll = len(string0)
      string1(1:ll) = string0
      call text_process(string1,n_components,iwords)

!
!  estimate number of points and patches on each component
!  

      do i=1,n_components
        fname1=trim(string1(iwords(i)+1:iwords(i+1)-1))
        call open_gov3_geometry_mem(fname1,npatches_vect(i), &
          npts_vect(i))
      enddo
!
!  compute source discretization
!  
!
      do i=1,n_components
        fname1=trim(string1(iwords(i)+1:iwords(i+1)-1))
        if (i.eq.1) then
          i1=1
          i2=npts_vect(1)
          j1=1
          j2=npatches_vect(1)
        else
          i1=i2+1
          i2=i1+npts_vect(i)-1
          j1=j2+1
          j2=j1+npatches_vect(i)-1
        endif
        call open_gov3_geometry_v2(fname1,npatches_vect(i), &
          norders(j1:j2),ixyzs(j1:j2+1),iptype(j1:j2), &
          npts_vect(i),srcvals(:,i1:i2),srccoefs(:,i1:i2), &
          wts(i1:i2),dP(1,i))
        if (i.ge.2) then
          do count2=j1,j2+1
            ixyzs(count2)=ixyzs(count2)+n_aux-1
          enddo
        endif
        n_aux=ixyzs(j2+1)
      enddo

      
      call topological_sorting(npatches,norders,ixyzs,iptype, &
       npts,srccoefs,srcvals,wts,npatches_vect,n_components, &
       sorted_vector,exposed_surfaces,eps)

      
      end subroutine em_solver_open_geom
!
!
!
!
!
      subroutine em_solver_wrap(string0,n_components,dP, &
        contrast_matrix,npts,omega,icase,direction,pol, &
        eps,eps_gmres,soln,err_est)
!
!f2py intent(in) string0,n_components,dP,contrast_matrix
!f2py intent(in) npts,omega,icase,direction,pol
!f2py intent(in) eps,eps_gmres
!f2py intent(out) soln,err_est
!
!  This subroutine solves the dielectric probleem either
!  for a known analytic solution constructed using the
!  point sources and plane waves/computes the solution
!  due to a scattering problem.
!
!
!  Input arguments:
!    - string0: character(len=*)
!      string defining the go3 files
!    - n_components: integer
!      number of components
!    - dP: real *8 (4,n_components)
!        dP(1:3,i) translation for the ith component
!        dP(4,i) scaling
!    - contrast_matrix: complex *16 (4,n_components)
!        \ep,\mu in positive and negative normal direction
!        of each component
!        * contrast_matrix(1,i) = \ep on the positive normal side
!        * contrast_matrix(2,i) = \mu on the positive normal side
!        * contrast_matrix(3,i) = \ep on the negative normal side
!        * contrast_matrix(4,i) = \mu on the negative normal side
!    - npts: integer
!        total number of discretization points on surface
!    - omega: real *8
!        wavenumber of problem (should be consistent with the units
!          of the go3 file)
!    - icase: integer
!        solve either for analytic solution or incident plane wave
!        (analytic solution in the interior regions is still the plane
!         wave)
!        * icase = 1, analytic solution test. Exterior solution
!          is constructed via dipole at centroid of geometry
!          whose orientation is in the (1,1,1) direction
!          and renormalized to rsc/omega (this makes sure exterior and
!          interior solutions are commensurate)
!        * icase = 2, plane wave computation
!    - direction: real *8(2)
!        direction(1) - phi angle of the plane wave
!        direction(2) - theta angle of the plane wave
!    - pol: complex *16(2)
!        pol(1) - phi vector of the plane wave
!        pol(2) - theta vector of plane wave
!    - eps: real *8
!        tolerance of FMM and quadrature generation
!    - eps_gmres: real *8
!        relative residual tolerance for GMRES
!    
!  Output arguments:
!    - soln: complex *16 (npts*4)
!        u,v projections of surface currents
!    - err_est: real *8 
!        estimated error. non zero only for icase=1
!        
!
!
      implicit real *8 (a-h,o-z) 
      character (len=*), intent(in) :: string0
      integer, intent(in) :: n_components
      real *8, intent(in) :: dP(4,n_components)
      complex *16, intent(in) :: contrast_matrix(4,n_components)
      integer, intent(in) :: npts
      real *8, intent(in) :: omega
      integer, intent(in) :: icase
      real *8, intent(in) :: direction(2)
      complex *16, intent(in) :: pol(2)
      real *8, intent(in) :: eps,eps_gmres

      complex *16, intent(out) :: soln(4*npts)
      real *8, intent(out) :: err_est

      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)
      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: rhs(:)

      complex *16 vf(3)


      real *8, allocatable :: errs(:)
      real *8 thet,phi
      complex * 16  zpars(5),omega_use

      integer numit,niter

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima,zk
      complex *16 alpha_rhs

      integer count1,count2,count3,icount,n_aux,i1,i2,j1,j2,nx,ny,nz


      integer ntarg 
      integer, allocatable :: npatches_vect(:),iwords(:),npts_vect(:)
      integer, allocatable :: sorted_vector(:),location_targs(:)

      character *1000 fname1    
      character *2000 string1
      
      real *8, allocatable :: srcvals_extended(:,:)

      real *8, allocatable :: targs(:,:)
      logical, allocatable :: exposed_surfaces(:)

      real *8 x_min,x_max,y_min,y_max,z_min,z_max,dx,dy,dz
      real *8 xa,ya,za
      real *8 pi,done

      data ima/(0.0d0,1.0d0)/


      done = 1
      pi = atan(done)*4

      call em_solver_wrap_mem(string0,n_components,npatches,npts0)

      allocate(npatches_vect(n_components),npts_vect(n_components))
      allocate(sorted_vector(n_components+1))
      allocate(exposed_surfaces(n_components))
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      allocate(srcvals_extended(20,npts))
      call em_solver_open_geom(string0,n_components,dP, &
        npatches,npts,eps,npatches_vect,npts_vect,norders,ixyzs,iptype, &
        srcvals,srccoefs,wts,sorted_vector,exposed_surfaces)
!
!  Rescale sources to make box unit size
!
      nt = 0
      call get_bsize(12,npts,srcvals,3,nt,tmp,bsize)

      rsc = 1.0d0/bsize
      
      do i=1,npts
        do j=1,9
          srcvals(j,i) = srcvals(j,i)*rsc
          srccoefs(j,i) = srccoefs(j,i)*rsc
        enddo
        wts(i) = wts(i)*rsc**2
      enddo

      omega_use = omega/rsc
      zpars(1) = omega_use


      allocate(rhs(4*npts))

!
! Compute centroid of points on first component
!

      xa = 0
      ya = 0
      za = 0
      ra = 0
      do i=1,npts_vect(1)
        xa = xa + srcvals(1,i)*wts(i)
        ya = ya + srcvals(2,i)*wts(i)
        za = za + srcvals(3,i)*wts(i)
        ra = ra + wts(i)
      enddo
      xa = xa/ra
      ya = ya/ra
      za = za/ra
      xyz_in(1) = xa    
      xyz_in(2) = ya
      xyz_in(3) = za


      vf(1)=1.0d0*1.0d-3
      vf(2)=2.0d0*1.0d-3
      vf(3)=3.0d0*1.0d-3

      if(icase.eq.1) then
        call get_rhs_em_muller_trans_testing(xyz_in,vf, &
          direction,pol,srcvals,omega_use,rhs, &
          n_components,npts_vect,contrast_matrix, &
          npts,exposed_surfaces)
       else
!
!  This alternative RHS is used to illuminate the dielectrics 
!  with an incoming plane wave
!  defined by the angles direction(1:2), 
!  that is phi and theta; and the polarization
!  Pol(1:2), that is, Phi vector and 
!  Theta vector of the plane wave.
!
        call get_rhs_em_muller_trans_PW(direction,Pol,srcvals, & 
          omega_use,rhs,n_components,npts_vect,contrast_matrix, &
          npts,exposed_surfaces)
      endif
!
!  build extended source and target values
!
      call build_extended_targ(n_components,srcvals_extended, &
        srcvals,npts_vect,contrast_matrix,npts)
  
      numit = 400
      niter = 0
      ifinout = 1

      allocate(errs(numit+1))
      do i=1,4*npts
        soln(i) = 0
      enddo


      call cpu_time(t1)
!C$      t1 = omp_get_wtime()    
      call em_muller_trans_v2_solver(npatches,norders,ixyzs,iptype, &
        npts,srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres, &
        niter,errs,rres,soln,contrast_matrix,npts_vect, &
        n_components,srcvals_extended)
      call cpu_time(t2)
!C$       t2 = omp_get_wtime()

!
!  if computing pw soln return 
! 
      if(icase.eq.2) return

!
!  Compute extremes of bounding box
!
!
      x_min=minval(srcvals(1,:))
      y_min=minval(srcvals(2,:))
      z_min=minval(srcvals(3,:))

      x_max=maxval(srcvals(1,:))
      y_max=maxval(srcvals(2,:))
      z_max=maxval(srcvals(3,:))
      
      dx=x_max-x_min
      dy=y_max-y_min
      dz=z_max-z_min

      x_min=x_min-dx/10.0d0
      y_min=y_min-dy/10.0d0
      z_min=z_min-dz/10.0d0

      x_max=x_max+dx/10.0d0
      y_max=y_max+dy/10.0d0
      z_max=z_max+dz/10.0d0

! this is for a 2D cut in the XY plane

      nx=20
      ny=20
      nz=1
      ntarg=nx*ny*nz
      allocate(targs(3,ntarg))
      allocate(location_targs(ntargs))

      icount=1
      do count1=1,nx
       do count2=1,ny
         targs(1,icount)=x_min+(x_max-x_min)*(count1-1)/(nx-1)
         targs(2,icount)=y_min+(y_max-y_min)*(count2-1)/(ny-1)
         targs(3,icount)=(z_max+z_min)/2.0d0
         icount=icount+1
       enddo
      enddo

      err_est = 0
      call test_accuracy_em_muller(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,wts,targs,ntarg,npatches_vect, &
        n_components,sorted_vector,contrast_matrix, &
        exposed_surfaces,eps,zpars,soln,xyz_in,vf, &
        direction,Pol,err_est)
      

      return
      end subroutine em_solver_wrap
!
!
!
!
!
!
      subroutine em_solver_wrap_postproc(string0,n_components,dP, &
        contrast_matrix,npts,omega,eps,soln,ntarg,targs,E,H)
!
!f2py intent(in) string0,n_components,dP,contrast_matrix
!f2py intent(in) npts,omega
!f2py intent(in) eps,soln,ntarg,targs
!f2py intent(out) E,H
!
!  This subroutine handles the postprocessing for the solution
!  of the dielectric problem, where the u,v projections of the
!  surface currents are already computed.
!
!  Input arguments:
!    - string0: character(len=*)
!      string defining the go3 files
!    - n_components: integer
!      number of components
!    - dP: real *8 (4,n_components)
!        dP(1:3,i) translation for the ith component
!        dP(4,i) scaling
!    - contrast_matrix: complex *16 (4,n_components)
!        \ep,\mu in positive and negative normal direction
!        of each component
!        * contrast_matrix(1,i) = \ep on the positive normal side
!        * contrast_matrix(2,i) = \mu on the positive normal side
!        * contrast_matrix(3,i) = \ep on the negative normal side
!        * contrast_matrix(4,i) = \mu on the negative normal side
!    - npts: integer
!        total number of discretization points on surface
!    - omega: real *8
!        wavenumber of problem (should be consistent with the units
!          of the go3 file)
!    - eps: real *8
!        tolerance of FMM and quadrature generation
!    - soln: complex *16 (npts*4)
!        u,v projections of surface currents
!    - ntarg: integer
!        number of targets
!    - targs: real *8 (3,ntarg)
!        location of targets
!    
!  Output arguments:
!    - E: complex *16(3,ntarg)
!        E field at target locations
!    - H: complex *16(3,ntarg)
!        H field at target locations
!        
!
!
      implicit real *8 (a-h,o-z) 
      character (len=*), intent(in) :: string0
      integer, intent(in) :: n_components
      real *8, intent(in) :: dP(4,n_components)
      complex *16, intent(in) :: contrast_matrix(4,n_components)
      integer, intent(in) :: npts
      real *8, intent(in) :: omega
      real *8, intent(in) :: eps
      complex *16, intent(in) :: soln(4*npts)
      integer, intent(in) :: ntarg
      real *8, intent(in) :: targs(3,ntarg)

      complex *16, intent(out) :: E(3,ntarg),H(3,ntarg)

      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      real *8, allocatable :: targ0(:,:)
      integer ipars(2)
      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)

      complex *16 vf(3)
      real *8, allocatable :: errs(:)
      real *8 thet,phi
      complex * 16  zpars(5),omega_use

      integer numit,niter

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima,zk
      complex *16 alpha_rhs

      integer count1,count2,count3,icount,n_aux,i1,i2,j1,j2,nx,ny,nz

      integer, allocatable :: npatches_vect(:),iwords(:),npts_vect(:)
      integer, allocatable :: sorted_vector(:),location_targs(:)

      character *1000 fname1    
      character *2000 string1
      
      real *8, allocatable :: srcvals_extended(:,:)

      logical, allocatable :: exposed_surfaces(:)

      real *8 x_min,x_max,y_min,y_max,z_min,z_max,dx,dy,dz
      real *8 xa,ya,za
      real *8 pi,done

      data ima/(0.0d0,1.0d0)/


      done = 1
      pi = atan(done)*4
      call em_solver_wrap_mem(string0,n_components,npatches,npts0)

      allocate(npatches_vect(n_components),npts_vect(n_components))
      allocate(sorted_vector(n_components+1))
      allocate(exposed_surfaces(n_components))

      allocate(srcvals_extended(20,npts))
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call em_solver_open_geom(string0,n_components,dP, &
        npatches,npts,eps,npatches_vect,npts_vect,norders,ixyzs,iptype, &
        srcvals,srccoefs,wts,sorted_vector,exposed_surfaces)


      nt = 0
      call get_bsize(12,npts,srcvals,3,nt,tmp,bsize)

      rsc = 1.0d0/bsize
      
      do i=1,npts
        do j=1,9
          srcvals(j,i) = srcvals(j,i)*rsc
          srccoefs(j,i) = srccoefs(j,i)*rsc
        enddo
        wts(i) = wts(i)*rsc**2
      enddo

      allocate(targ0(3,ntarg))

      do i=1,ntarg
        targ0(1,i) = targs(1,i)*rsc
        targ0(2,i) = targs(2,i)*rsc
        targ0(3,i) = targs(3,i)*rsc
      enddo

      omega_use = omega/rsc

      zpars(1) = omega_use

      call evaluate_field_muller(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,wts,targ0,ntarg,npatches_vect, &
        n_components,sorted_vector,contrast_matrix,exposed_surfaces, &
        eps,zpars,soln,E,H)

      

      return
      end subroutine em_solver_wrap_postproc
!
!
!
!
!
      subroutine em_sol_exact(string0,n_components,dP, &
        contrast_matrix,npts,omega,eps,direction,pol,ntarg,targs,E,H)
!
!
!f2py intent(in) string0,n_components,dP,contrast_matrix
!f2py intent(in) npts,omega,direction,pol
!f2py intent(in) eps,ntarg,targs
!f2py intent(out) E,H
!
!  This subroutine handles the postprocessing for the solution
!  This subroutine handles the postprocessing for the solution
!  of the dielectric problem, and computes the exact analytic
!  solution used at a collection of target points
!
!  Input arguments:
!    - string0: character(len=*)
!      string defining the go3 files
!    - n_components: integer
!      number of components
!    - dP: real *8 (4,n_components)
!        dP(1:3,i) translation for the ith component
!        dP(4,i) scaling
!    - contrast_matrix: complex *16 (4,n_components)
!        \ep,\mu in positive and negative normal direction
!        of each component
!        * contrast_matrix(1,i) = \ep on the positive normal side
!        * contrast_matrix(2,i) = \mu on the positive normal side
!        * contrast_matrix(3,i) = \ep on the negative normal side
!        * contrast_matrix(4,i) = \mu on the negative normal side
!    - npts: integer
!        total number of discretization points on surface
!    - omega: real *8
!        wavenumber of problem (should be consistent with the units
!          of the go3 file)
!    - eps: real *8
!        tolerance of FMM and quadrature generation
!    - direction: real *8(2)
!        direction(1) - phi angle of the plane wave
!        direction(2) - theta angle of the plane wave
!    - pol: complex *16(2)
!        pol(1) - phi vector of the plane wave
!        pol(2) - theta vector of plane wave
!    - ntarg: integer
!        number of targets
!    - targs: real *8 (3,ntarg)
!        location of targets
!    
!  Output arguments:
!    - E: complex *16(3,ntarg)
!        exact E field at target locations
!    - H: complex *16(3,ntarg)
!        exact H field at target locations
!        
!
!
      implicit real *8 (a-h,o-z) 
      character (len=*), intent(in) :: string0
      integer, intent(in) :: n_components
      real *8, intent(in) :: dP(4,n_components)
      complex *16, intent(in) :: contrast_matrix(4,n_components)
      integer, intent(in) :: npts
      real *8, intent(in) :: omega
      real *8, intent(in) :: eps
      real *8, intent(in) :: direction(2)
      complex *16, intent(in) :: pol(2)
      integer, intent(in) :: ntarg
      real *8, intent(in) :: targs(3,ntarg)

      complex *16, intent(out) :: E(3,ntarg),H(3,ntarg)

      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)
      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)

      complex *16 vf(3)
      real *8, allocatable :: errs(:)
      real *8 thet,phi
      complex * 16  zpars(5),omega_use

      integer numit,niter

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima,zk
      complex *16 alpha_rhs

      integer count1,count2,count3,icount,n_aux,i1,i2,j1,j2,nx,ny,nz


      integer, allocatable :: npatches_vect(:),iwords(:),npts_vect(:)
      integer, allocatable :: sorted_vector(:),location_targs(:)

      character *1000 fname1    
      character *2000 string1
      
      real *8, allocatable :: srcvals_extended(:,:)
      real *8, allocatable :: targ0(:,:)

      logical, allocatable :: exposed_surfaces(:)

      real *8 x_min,x_max,y_min,y_max,z_min,z_max,dx,dy,dz
      real *8 xa,ya,za
      real *8 pi,done

      data ima/(0.0d0,1.0d0)/


      done = 1
      pi = atan(done)*4

      call em_solver_wrap_mem(string0,n_components,npatches,npts0)
      allocate(npatches_vect(n_components),npts_vect(n_components))
      allocate(sorted_vector(n_components+1))
      allocate(exposed_surfaces(n_components))
      allocate(srcvals_extended(20,npts))
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call em_solver_open_geom(string0,n_components,dP, &
        npatches,npts,eps,npatches_vect,npts_vect,norders,ixyzs,iptype, &
        srcvals,srccoefs,wts,sorted_vector,exposed_surfaces)

      nt = 0
      call get_bsize(12,npts,srcvals,3,nt,tmp,bsize)

      rsc = 1.0d0/bsize
      
      do i=1,npts
        do j=1,9
          srcvals(j,i) = srcvals(j,i)*rsc
          srccoefs(j,i) = srccoefs(j,i)*rsc
        enddo
        wts(i) = wts(i)*rsc**2
      enddo

      allocate(targ0(3,ntarg))

      do i=1,ntarg
        targ0(1,i) = targs(1,i)*rsc
        targ0(2,i) = targs(2,i)*rsc
        targ0(3,i) = targs(3,i)*rsc
      enddo

      omega_use = omega/rsc

      zpars(1) = omega_use


      
!
! Compute centroid of points on first component
!

      xa = 0
      ya = 0
      za = 0
      ra = 0
      do i=1,npts_vect(1)
        xa = xa + srcvals(1,i)*wts(i)
        ya = ya + srcvals(2,i)*wts(i)
        za = za + srcvals(3,i)*wts(i)
        ra = ra + wts(i)
      enddo
      xa = xa/ra
      ya = ya/ra
      za = za/ra
      xyz_in(1) = xa    
      xyz_in(2) = ya
      xyz_in(3) = za


      vf(1)=1.0d0*1.0d-3
      vf(2)=2.0d0*1.0d-3
      vf(3)=3.0d0*1.0d-3

      
      call evaluate_field_muller_exact(npatches,norders,ixyzs,iptype, &
        npts,srccoefs,srcvals,wts,targ0,ntarg,npatches_vect, &
        n_components,sorted_vector,contrast_matrix,exposed_surfaces, &
        eps,zpars,xyz_in,vf,direction,pol,E,H)


      return
      end subroutine em_sol_exact
!
!
!
!
!
!
      subroutine em_plot_surf_current_vtk(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,soln,ifjk,fname,title)
!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs
!f2py intent(in) srcvals,soln,ifjk,fname,title
!
!   This subroutine writes a vtk to plot the surface along
!   with a tangential vector field prescribed through its projection
!   onto X_{u}, X_{v}. 
!   Currently only supports triangular patches
!
!
      implicit none
      integer, intent(in) :: npatches,norders(npatches),ixyzs(npatches+1),npts
      integer, intent(in) :: iptype(npatches),ifjk
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      complex *16, intent(in) :: soln(4*npts)
      character (len=*), intent(in) :: fname,title
      complex *16, allocatable :: sigma(:,:)
      real *8, allocatable :: u_vect_s(:,:),v_vect_s(:,:)
      integer i

      allocate(sigma(3,npts),u_vect_s(3,npts),v_vect_s(3,npts))
      call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect_s, &
         v_vect_s,npts)

      if(ifjk.eq.1) then
        do i=1,npts
          sigma(1:3,i) = soln(i)*u_vect_s(1:3,i) + &
            soln(i+npts)*v_vect_s(1:3,i)
        enddo
      else
        do i=1,npts
          sigma(1:3,i) = soln(i+2*npts)*u_vect_s(1:3,i) + &
            soln(i+3*npts)*v_vect_s(1:3,i)
        enddo
      endif

      call surf_vtk_plot_zvec(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,sigma,fname,title)


      end subroutine em_plot_surf_current_vtk
!
!
