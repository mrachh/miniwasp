      implicit real *8 (a-h,o-z)
      character *2000 string1,string2
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
      real *8, allocatable :: verts(:,:)
      integer, allocatable :: elements(:,:)
      integer, allocatable :: sorted_vector(:)
      logical, allocatable :: exposed_surfaces(:)
      real *8, allocatable :: fvals(:,:),fvalsout(:,:)
      integer, allocatable :: iver_el_list(:),iverstart(:),iverind(:)
      real *8 cm(3),inertia_tensor(3,3),uinertia(3,3),vinertia(3,3)
      real *8 sing_inertia(3)
      real *8 avals(3,2),xtest(3)
      complex *16 pol(2)
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
      
      call prini(6,13)

      done = 1
      pi = atan(done)*4

      n_components = 1
      ncmax = 3
      allocate(contrast_matrix(4,ncmax),dP(4,ncmax))
      contrast_matrix(1,1)=1.1d0  
      contrast_matrix(2,1)=1.1d0  
      contrast_matrix(3,1)=1.2d0  
      contrast_matrix(4,1)=1.0d0

      contrast_matrix(1,2)=1.1d0  
      contrast_matrix(2,2)=1.1d0  
      contrast_matrix(3,2)=1.2d0  
      contrast_matrix(4,2)=1.0d0

      contrast_matrix(1,3)=1.1d0  
      contrast_matrix(2,3)=1.1d0  
      contrast_matrix(3,3)=1.2d0  
      contrast_matrix(4,3)=1.0d0

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

      

!      string1 = '../geometries/simplest_cube_quadratic_v4_o08_r02.go3?'
      string1 = '../geometries/rhabdom_r00.go3?'

      eps = 0.51d-3
      eps_gmres = 0.51d-5
!
!  Now do the same computation using the postprocessing routine
!  and getting the analytic solution
!
! For that we first need to extract geometry info
      
      call em_solver_wrap_mem(string1,n_components,npatches,npts)
      allocate(npatches_vect(n_components),npts_vect(n_components))
      allocate(sorted_vector(n_components+1))
      allocate(exposed_surfaces(n_components))
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))
      call em_solver_open_geom(string1,n_components,dP, &
        npatches,npts,eps,npatches_vect,npts_vect,norders,ixyzs,iptype, &
        srcvals,srccoefs,wts,sorted_vector,exposed_surfaces)

      cm(1) = 0
      cm(2) = 0
      cm(3) = 0
      ra = 0
      do i=1,npts
        ra = ra + wts(i)
        cm(1) = cm(1) + srcvals(1,i)*wts(i)
        cm(2) = cm(2) + srcvals(2,i)*wts(i)
        cm(3) = cm(3) + srcvals(3,i)*wts(i)
      enddo

      cm(1) = cm(1)/ra
      cm(2) = cm(2)/ra
      cm(3) = cm(3)/ra

      ra2 = 0
      do i=1,npts
        do l=1,3
          ra2 = ra2 + (srcvals(l,i)-cm(l))**2*wts(i)
        enddo
      enddo

      call prin2('ra2=*',ra2,1)


      do j=1,3
        do l=1,3
          inertia_tensor(j,l) = 0
          do i=1,npts
            inertia_tensor(j,l) = inertia_tensor(j,l) - &
              (srcvals(j,i)-cm(j))*(srcvals(l,i)-cm(l))*wts(i)
          enddo
        enddo
      enddo

      do l=1,3
        inertia_tensor(l,l) = inertia_tensor(l,l) + ra2
      enddo

      call prin2('inertia_tensor=*',inertia_tensor,9)

      call dsvd(3,3,inertia_tensor,uinertia,sing_inertia,vinertia)
      call prin2('sing_inertia=*',sing_inertia,3)
      call prin2('uinertia=*',uinertia,9)
      call prin2('vinertia=*',vinertia,9)


!
!  find bounding box
!
!
      xmin = minval(srcvals(1,:))
      xmax = maxval(srcvals(1,:))

      ymin = minval(srcvals(2,:))
      ymax = maxval(srcvals(2,:))

      zmin = minval(srcvals(3,:))
      zmax = maxval(srcvals(3,:))

      xsize = xmax - xmin
      ysize = ymax - ymin
      zsize = zmax - zmin

      bsize = max(xsize,ysize)
      bsize = max(bsize,zsize)

      
      avals(1:3,1) = cm(1:3) - bsize*uinertia(1:3,3)
      avals(1:3,2) = cm(1:3) + bsize*uinertia(1:3,3)

      call vtk_curv_plot(2,3,avals,'axis.vtk','a')


      call prin2('avals=*',avals,6)

      theta = 0
      phi = 0

      call cart2polar(uinertia(1:3,3),r,theta,phi)
      call prin2('theta=*',theta,1)
      call prin2('phi=*',phi,1)

      call prin2_long('theta=*',theta,1)
      call prin2_long('phi=*',phi,1)

      xtest(1) = r*sin(theta)*cos(phi)
      xtest(2) = r*sin(theta)*sin(phi)
      xtest(3) = r*cos(theta)

      erra = abs(xtest(1)-uinertia(1,3))+abs(xtest(2) - uinertia(2,3))+&
        abs(xtest(3)-uinertia(3,3))
      call prin2('erra=*',erra,1)


      call surf_quadratic_msh_vtk_plot(npatches,norders,ixyzs, iptype, &
        npts,srccoefs,srcvals,'eye.vtk','a')
      

      stop
      end
!
!
!
!
!


      subroutine vtk_curv_plot(n,nda,avals,fname,title)
!
! This subroutine writes a vtk file to plot a given a curve 
! as a collection of line segments
!
!
!  Input arguments:
!    - n: integer
!        number of points
!    - nda: integer
!        leading dimension of data array
!    - avals: real *8 (nda,n)
!        curve to be plotted, the first three components
!        must be xyz coordinates
!    - fname: character (len=*)
!        file name where vtk output should be written
!
      implicit none
      integer, intent(in) :: nda,n
      real *8, intent(in) :: avals(nda,n)
      character (len=*), intent(in) :: fname,title


      integer i,j,k,l,n0,npuv,ipatch,ipt,i1,m,norder,npols,iunit1

  
      iunit1 = 877
      open(unit = iunit1, file=trim(fname), status='replace')

      write(iunit1,'(a)') "# vtk DataFile Version 3.0"
      write(iunit1,'(a)') trim(title)
      write(iunit1,'(a)') "ASCII"
      write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
      write(iunit1,'(a,i9,a)') "POINTS ", n, " float"

      do i = 1,n
        write(iunit1,"(E11.5,2(2x,e11.5))") avals(1,i), avals(2,i), avals(3,i)
      end do

      write(iunit1,'(a,i9,i9)') "CELLS ", n, n*3

      do i=1,n
        if(i.ne.n)  write(iunit1,'(a,i9,i9)') "2 ", i-1, i
        if(i.eq.n)  write(iunit1,'(a,i9,i9)') "2 ", i-1, 0
      enddo

      write(iunit1,'(a,i9)') "CELL_TYPES ", n
      do ipatch = 1,n
        write(iunit1,'(a)') "3"
      end do

      close(iunit1)

      end subroutine vtk_curv_plot
!
!
!
!
!

