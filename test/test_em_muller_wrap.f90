      implicit real *8 (a-h,o-z)
      character *2000 string1
      integer n_components
      complex *16, allocatable :: soln(:)
      complex *16, allocatable :: contrast_matrix(:,:)
      real *8, allocatable :: dP(:,:)
      real *8 direction(2),omega,eps,eps_gmres
      real *8, allocatable :: targs(:,:)
      complex *16, allocatable :: E(:,:),H(:,:),Eex(:,:),Hex(:,:)
      integer, allocatable :: npatches_vect(:),npts_vect(:)
      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),wts(:)
      integer, allocatable :: sorted_vector(:)
      logical *8, allocatable :: exposed_surfaces(:)
      complex *16 pol(2)
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
      
      call prini(6,13)

      done = 1
      pi = atan(done)*4

      n_components = 1
      allocate(contrast_matrix(4,n_components),dP(4,n_components))
      contrast_matrix(1,1)=1.1d0  
      contrast_matrix(2,1)=1.1d0  
      contrast_matrix(3,1)=1.2d0  
      contrast_matrix(4,1)=1.0d0

      dP(1,1) = 0
      dP(2,1) = 0
      dP(3,1) = 0
      dP(4,1) = 1.1d0

!      string1 = '../geometries/simplest_cube_quadratic_v4_o08_r02.go3?'
      string1 = '../geometries/lens_r00.go3?'

!       estimate number of discretization points      
      call em_solver_wrap_mem(string1,n_components,npatches,npts)
      call prinf('npts=*',npts,1)

      omega = 2.0d0*pi/1200
      icase = 1
      eps = 0.51d-3
      eps_gmres = 0.51d-6
      allocate(soln(4*npts))
!
!  run the solver
!
      direction(1) = 0.0d0
      direction(2) = pi/2.0d0
      pol(1) = 3.0d0
      pol(2) = 3.0d0
      err_est = 0
      call em_solver_wrap(string1,n_components,dP,contrast_matrix,npts,&
        omega,icase,direction,pol,eps,eps_gmres,soln,err_est)
      call prin2('estimated error=*',err_est,1)
      stop
!
!  Now do the same computation using the postprocessing routine
!  and getting the analytic solution
!
! For that we first need to extract geometry info
      
      allocate(npatches_vect(n_components),npts_vect(n_components))
      allocate(sorted_vector(n_components+1))
      allocate(exposed_surfaces(n_components))
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))
      call em_solver_open_geom(string1,n_components,dP, &
        npatches,npts,eps,npatches_vect,npts_vect,norders,ixyzs,iptype, &
        srcvals,srccoefs,wts,sorted_vector,exposed_surfaces)

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

      allocate(E(3,ntarg),H(3,ntarg))

      icount=1
      do count1=1,nx
       do count2=1,ny
         targs(1,icount)=x_min+(x_max-x_min)*(count1-1)/(nx-1)
         targs(2,icount)=y_min+(y_max-y_min)*(count2-1)/(ny-1)
         targs(3,icount)=(z_max+z_min)/2.0d0
         icount=icount+1
       enddo
      enddo

      do i=1,ntarg
        E(1,i) = 0
        E(2,i) = 0
        E(3,i) = 0

        H(1,i) = 0
        H(2,i) = 0
        H(3,i) = 0
      enddo
      call prin2('targs=*',targs,24)
      soln = 1

      call em_solver_wrap_postproc(string1,n_components,dP, &
        contrast_matrix,npts,omega,eps,soln,ntarg,targs,E,H)



      stop
      end
