      implicit real *8 (a-h,o-z)
      character *2000 string1
      character *300 fname1,fname2,fname3,fsol,ftarg,fname4
      integer n_components
      complex *16, allocatable :: soln(:)
      complex *16, allocatable :: contrast_matrix(:,:)
      real *8, allocatable :: dP(:,:)
      real *8 direction(2),omega,eps,eps_gmres
      real *8, allocatable :: targs(:,:)
      complex *16, allocatable :: E(:,:),H(:,:),E_ex(:,:),H_ex(:,:)
      integer, allocatable :: npatches_vect(:),npts_vect(:)
      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),wts(:)
      real *8, allocatable :: cms(:,:),rads(:)
      integer, allocatable :: sorted_vector(:)
      logical, allocatable :: exposed_surfaces(:)
      complex *16 pol(2)
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
      
      call prini(6,13)

      done = 1
      pi = atan(done)*4

      iref = 1 


      ncmax = 4
      n_components = 4
      allocate(contrast_matrix(4,ncmax),dP(4,ncmax))
      contrast_matrix(1,1)=1.34d0**2  
      contrast_matrix(2,1)=1.0d0  
      contrast_matrix(3,1)=1.452d0**2  
      contrast_matrix(4,1)=1.0d0

      contrast_matrix(1,2)=1.34d0**2  
      contrast_matrix(2,2)=1.0d0  
      contrast_matrix(3,2)=1.348d0**2  
      contrast_matrix(4,2)=1.0d0

      contrast_matrix(1,3)=1.34d0**2  
      contrast_matrix(2,3)=1.0d0  
      contrast_matrix(3,3)=1.363d0**2  
      contrast_matrix(4,3)=1.0d0

      contrast_matrix(1,4) = 1.0d0
      contrast_matrix(2,4) = 1.0d0
      contrast_matrix(3,4) = 1.34d0**2
      contrast_matrix(4,4) = 1.0d0


      dP(1,1) = 0
      dP(2,1) = 0
      dP(3,1) = 0
      dP(4,1) = 1.0d0


      dP(1,2) = 0
      dP(2,2) = 0
      dP(3,2) = 0
      dP(4,2) = 1.0d0

      dP(1,3) = 0
      dP(2,3) = 0
      dP(3,3) = 0
      dP(4,3) = 1.0d0

      dP(1,4) = 71563.516d0
      dP(2,4) = 14128.477d0 
      dP(3,4) = 40817.917d0
      dP(4,4) = 10206.791d0

      write(fname1,'(a,i2.2,a)') '../geometries/lens_r',iref,'.go3?'
      write(fname2,'(a,i2.2,a)') '../geometries/con_r',iref,'.go3?'
      write(fname3,'(a,i2.2,a)') '../geometries/rhabdom_r',iref,'.go3?'
      write(fname4,'(a,i2.2,a)') &
        '../geometries/sphere_r',iref+3,'_o01.go3?'


!      string1 = '../geometries/simplest_cube_quadratic_v4_o08_r02.go3?'
      string1 = trim(fname1)//trim(fname2)//trim(fname3)//trim(fname4)
!

      print *, trim(string1)


!       estimate number of discretization points      
      call em_solver_wrap_mem(string1,n_components,npatches,npts)

      iomegafrac = 6
!      omega = 2.0d0*pi/600/(iomegafrac+0.0d0)
      omega = 2.0d0*pi/350.0d0
      icase = 1
      eps = 0.51d-3
      eps_gmres = 0.51d-5
      allocate(soln(4*npts))
!
!  run the solver
!
!
      idircase = 1

!
!  default values
!
      direction(1) = 0
      direction(2) = pi/2

      if(idircase.eq.1) then
        direction(1) = -1.530844785708483D0 
        direction(2) = 1.063728846738649D0 
      endif


      ipolcase = 1
!
!  default values
!
      
      pol(1) = 3.0d0
      pol(2) = 3.0d0

      if(ipolcase.eq.1) then
        pol(1) = 1.0d0
        pol(2) = 0.0d0
      else if(ipolcase.eq.2) then
        pol(1) = 0.0d0
        pol(2) = 1.0d0
      endif
      err_est = 0


      write(fsol,'(a,i1,a,i1,a,i1,a,i1,a,i1,a)')  &
        'data/soln_iref',iref,'_iomega',iomegafrac, &
        '_icase',icase,'_idir',idircase,'_ipol',ipolcase, &
        '.dat'

      write(ftarg,'(a,i1,a,i1,a,i1,a,i1,a,i1,a)')  &
      'data/targ_iref',iref,'_iomega',iomegafrac, &
        '_icase',icase,'_idir',idircase,'_ipol',ipolcase, &
        '.dat'
      print *, trim(fsol)
      print *, trim(ftarg)

      call prin2('omega=*',omega,1)
      call prin2('pol=*',pol,4)
      call prin2('direction=*',direction,2)




      call em_solver_wrap(string1,n_components,dP,contrast_matrix,npts,&
        omega,icase,direction,pol,eps,eps_gmres,soln,err_est)

      open(unit=33,file=trim(fsol))
      do i=1,4*npts
        write(33,'(2(2x,d22.16))') real(soln(i)),imag(soln(i))
      enddo

      close(33)


      err_est = 0
      call prin2('estimated error=*',err_est,1)
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
      x_min=  69211.29d0
      y_min = 5963.045d0
      z_min = 35501.2415d0

      x_max = 73915.741d0
      y_max = 35501.2415d0
      z_max = 46134.592d0

      if(1.eq.0) then
        x_min = -1.5d0    
        x_max = 1.5d0

        y_min = -1.5d0    
        y_max = 1.5d0

        z_min = -1.5d0    
        z_max = 1.5d0
      endif
      dx=x_max-x_min
      dy=y_max-y_min
      dz=z_max-z_min

      x_min=x_min-dx/2.0d0
      y_min=y_min-dy/2.0d0
      z_min=z_min-dz/2.0d0

      x_max=x_max+dx/2.0d0
      y_max=y_max+dy/2.0d0
      z_max=z_max+dz/2.0d0

! this is for a 2D cut in the XY plane

      nx=20
      ny=20
      nz=1
      ntarg=nx*ny*nz
      allocate(targs(3,ntarg))

      allocate(E(3,ntarg),H(3,ntarg),E_ex(3,ntarg),H_ex(3,ntarg))

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

        E_ex(1,i) = 0
        E_ex(2,i) = 0
        E_ex(3,i) = 0

        H_ex(1,i) = 0
        H_ex(2,i) = 0
        H_ex(3,i) = 0
      enddo

      call em_solver_wrap_postproc(string1,n_components,dP, &
        contrast_matrix,npts,omega,eps,soln,ntarg,targs,E,H)

      open(unit=34,file=trim(ftarg))
      do i=1,ntarg
        write(34,'(15(2x,d22.16))') targs(1,i),targs(2,i),targs(3,i), &
          real(E(1,i)),imag(E(1,i)),real(E(2,i)),imag(E(2,i)), &
          real(E(3,i)),imag(E(3,i)), &
          real(H(1,i)),imag(H(1,i)),real(H(2,i)),imag(H(2,i)), &
          real(H(3,i)),imag(H(3,i))
        
      enddo
      close(34)



      call em_sol_exact(string1,n_components,dP, &
        contrast_matrix,npts,omega,eps,direction,pol,ntarg,targs, &
        E_ex,H_ex)

      erra = 0
      ra = 0
      do i=1,ntarg
        do j=1,3
          erra = erra + abs(E(j,i)-E_ex(j,i))**2
          erra = erra + abs(H(j,i)-H_ex(j,i))**2
          ra = ra + abs(E_ex(j,i))**2
          ra = ra + abs(H_ex(j,i))**2
        enddo
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in E and H fields=*',erra,1)


      stop
      end
