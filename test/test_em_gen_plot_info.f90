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
      real *8, allocatable :: rscales(:)
      complex *16 pol(2)
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
      
      call prini(6,13)

      done = 1
      pi = atan(done)*4

      n_components = 2
      allocate(contrast_matrix(4,n_components),dP(4,n_components))
      contrast_matrix(1,1)=1.1d0  
      contrast_matrix(2,1)=1.1d0  
      contrast_matrix(3,1)=1.2d0  
      contrast_matrix(4,1)=1.0d0

      contrast_matrix(1,2)=1.1d0  
      contrast_matrix(2,2)=1.1d0  
      contrast_matrix(3,2)=1.2d0  
      contrast_matrix(4,2)=1.0d0

      dP(1,1) = 0
      dP(2,1) = 0
      dP(3,1) = 0
      dP(4,1) = 1.0d0

      dP(1,2) = 0
      dP(2,2) = 0
      dP(3,2) = 0
      dP(4,2) = 1.0d0

!      string1 = '../geometries/simplest_cube_quadratic_v4_o08_r02.go3?'
      string1 = '../geometries/lens_r00.go3?../geometries/con_r00.go3?'
      string2 = '../geometries/lens_r00.msh?../geometries/con_r00.msh?'
!      string1 = '../geometries/sphere_r02_o03.go3?'

!       estimate number of discretization points      
      call em_gen_plot_info_surf_mem(string2,n_components,nverts,nel)
      call prinf('nverts=*',nverts,1)
      call prinf('nel=*',nel,1)

      allocate(verts(3,nverts),elements(3,nel))
      call em_gen_plot_info_surf(string2,n_components,nverts,verts, &
        nel,elements)
      call prin2('verts=*',verts,24)
      call prinf('elements=*',elements,20)

      allocate(iver_el_list(3*nel),iverstart(nverts+1),iverind(3*nel))

      call em_elem_trans(nel,nverts,elements,iver_el_list,iverstart, &
        iverind)
      call prinf('iverind=*',iverind,20)
      call prinf('iver_el_list=*',iver_el_list,20)
      call prinf('iverstart=*',iverstart,20)


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


      nd = 3
      allocate(fvals(nd,npts),fvalsout(nd,npts))

      do i=1,npts
        fvals(1,i) = srcvals(1,i)
        fvals(2,i) = srcvals(2,i)
        fvals(3,i) = srcvals(3,i)
        fvalsout(1,i) = 0
        fvalsout(2,i) = 0
        fvalsout(3,i) = 0
      enddo

      call prinf('npatches=*',npatches,1)

      allocate(rscales(nel))
      do i=1,nel
        rscales(i) = 1.0d0
      enddo
      
      call em_surf_fun_to_plot_fun(nd,npatches,norders,ixyzs,&
        iptype,npts,fvals, &
        nverts,nel,iver_el_list,iverstart,iverind,rscales,fvalsout)
!
!
!  estimate error by comparing to actual vertices
!
!
      erra = 0
      ra = 0
      do i=1,nverts
        do j=1,3
          erra = erra + (verts(j,i)-fvalsout(j,i))**2
          ra = ra + verts(j,i)**2
        enddo
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in vertices=*',erra,1)


      stop
      end
