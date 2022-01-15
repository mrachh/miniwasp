      implicit real *8 (a-h,o-z) 
      real *8 abc_lens(3),abc_cone(3),abc_rhabdom(3),abc_pig(3)
      real *8 xyz_lens(3),xyz_cone(3),xyz_rhabdom(3),xyz_pig(3)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:),rhs(:)
      complex *16, allocatable :: sigma2(:),rhs2(:)
      complex *16, allocatable :: currents(:,:)
      real *8, allocatable :: u_vect_s(:,:),v_vect_s(:,:),errp(:)
      real *8, allocatable :: errp_plot(:)
      complex *16 vf(3)


      real *8, allocatable :: errs(:)
      real *8 thet,phi
      complex * 16  zpars(5), omega, ep0, mu0, ep1, mu1

      integer numit,niter
      character *200 title,fname

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima,zk
      complex *16 alpha_rhs

      integer count1,count2,count3,icount,n_aux,i1,i2,j1,j2,nx,ny,nz

!!! You need to edit these variables to change the number of components:

      integer n_components,ntarg     
      integer, allocatable :: npatches_vect(:),iwords(:),npts_vect(:)
      integer, allocatable :: sorted_vector(:),location_targs(:)
      complex *16, allocatable :: contrast_matrix(:,:)  

      character *200 fname1    
      character *2000 string1
      
      real *8, allocatable :: srcvals_extended(:,:)

      real *8, allocatable :: dP(:,:),targs(:,:)
      logical, allocatable :: exposed_surfaces(:)

      real *8 direction(2)
      complex *16 Pol(2)

      real *8 x_min,x_max,y_min,y_max,z_min,z_max,dx,dy,dz
      real *8 xa,ya,za


      complex *16, allocatable :: E_far(:,:),H_far(:,:)

      real *8 pi,done

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


!

! n_components is the number of different dielectrics (surfaces with contrast)

!     READ***
     n_components=3
     n_components0=4

     allocate(npatches_vect(n_components0),npts_vect(n_components0))
     allocate(contrast_matrix(4,n_components0))
     allocate(dP(4,n_components0))
     allocate(sorted_vector(n_components0+1))
     allocate(exposed_surfaces(n_components0))


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


      iref = 0
      npatches = 0
      call get_miniwasp_ellip_params(abc_lens,abc_cone, &
        abc_rhabdom,abc_pig,xyz_lens,xyz_cone,xyz_rhabdom,xyz_pig)

      call prin2('abc_lens=*',abc_lens,3)
!      call prin2('abc_cone=*',abc_cone,3)
!      call prin2('abc_rhabdom=*',abc_rhabdom,3)
!      call prin2('abc_pig=*',abc_pig,3)
      call prin2('xyz_lens=*',xyz_lens,3)
!      call prin2('xyz_cone=*',xyz_cone,3)
!      call prin2('xyz_rhabdom=*',xyz_rhabdom,3)
!      call prin2('xyz_pig=*',xyz_pig,3)

      call get_miniwasp_ellip_mem(n_components,abc_lens,abc_cone, &
        abc_rhabdom,abc_pig,iref,npatches_vect,npatches)
      call prinf('npatches=*',npatches,1)

      norder = 5
      npols = (norder+1)*(norder+2)/2
      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts),wts(npts))
      call get_miniwasp_ellip(n_components,abc_lens,abc_cone,abc_rhabdom, &
        abc_pig,xyz_lens,xyz_cone,xyz_rhabdom,xyz_pig,iref,npatches, &
        npatches_vect,norder,srcvals,srccoefs,npts,npts_vect)

      allocate(iptype(npatches),ixyzs(npatches+1),norders(npatches))

      do i=1,npatches
        ixyzs(i) = (i-1)*npols+1
        norders(i) = norder
        iptype(i) = 1
      enddo
      ixyzs(npatches+1) = 1+npts

      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)




!  dP(:,n) defines an affine transform to be applied to the n'th .go3 geometry
!  (x,y,z) -> (x',y',z')
!
!  x'=dP(4,n)*x+dP(1,n)
!  y'=dP(4,n)*y+dP(2,n)
!  z'=dP(4,n)*z+dP(3,n)
!
   do count1=1,n_components
      dP(:,count1) = (/0.0d0, 0.0d0, 0.0d0, 1.0d0 /)   
   enddo  

    allocate(srcvals_extended(20,npts))

    call prinf('npatches=*',npatches,1)
    call prinf('npts=*',npts,1)
    print *, "rat=",npts/(npatches+0.0d0)

!
!	Electromagnetic parameters
!
!     READ***
!     zk = omega*sqrt(epsilon*mu)  -> 2*pi/lambda_0
      omega=1.12d0


      zpars(1) = omega

      write (*,*) 'omega: ',omega

      allocate(sigma(4*npts),rhs(4*npts))
      ifinout = 1

      thet = hkrand(0)*pi
      phi = hkrand(0)*2*pi

!
!	Set tolerances
!
      eps = 1d-5    ! <- this eps is used in the fmm, 
                    ! for the quadratures and in the
                    ! sorting algorithm as tolerance to evaluate 
                    ! the double layer

      eps_gmres=1d-6

      call topological_sorting(npatches,norders,ixyzs,iptype, &
        npts,srccoefs,srcvals,wts,npatches_vect,n_components, &
        sorted_vector,exposed_surfaces,eps)

      write (*,*) 'sorted vector: ', sorted_vector(1:2)

!
! Parameters of the electric and magnetic dipoles 
!
      xyz_in(1:3) = xyz_lens(1:3)
      xyz_in(1) = xyz_in(1) + 1.0d-3
      xyz_in(2) = xyz_in(2) - 2.0d-3
      xyz_in(3) = xyz_in(3) + 3.0d-3
      

! orientation of the electric and magnetic dipoles used to create the solution in the 
! exterior region and test the accuracy

      vf(1)=1.0d0
      vf(2)=2.0d0
      vf(3)=3.0d0

!
! Plane wave parameters
!

!     READ***
      direction(1)=0.0d0    !phi angle of the plane wave
      direction(2)=pi/2.0d0 !theta angle of the plane wave
      Pol(1)=3.0d0          !Phi vector of the plane wave
      Pol(2)=3.0d0          !Theta vector of the plane wave

!
!  This RHS is used to test accuracy
!

      call get_rhs_em_muller_trans_testing(xyz_in,vf,direction, &
        Pol,srcvals,omega,rhs,n_components,npts_vect, &
        contrast_matrix,npts,exposed_surfaces)
!
!  This alternative RHS is used to illuminate the dielectrics with an incoming plane wave
!  defined by the angles direction(1:2), that is phi and theta; and the polarization
!  Pol(1:2), that is, Phi vector and Theta vector of the plane wave.
!

!      call get_rhs_em_muller_trans_PW(direction,Pol,srcvals,omega,rhs,&
!       &n_components,npts_vect,contrast_matrix,npts,exposed_surfaces)


      call build_extended_targ(n_components,srcvals_extended, &
        srcvals,npts_vect,contrast_matrix,npts)
  
      numit = 400
      niter = 0
      allocate(errs(numit+1))

      sigma= 0
!
!	Specify
!
!      goto 2111
      call cpu_time(t1)
!C$      t1 = omp_get_wtime()    
      call em_muller_trans_v2_solver(npatches,norders,ixyzs,iptype, &
       npts,srccoefs,srcvals,eps,zpars,numit,ifinout,rhs, &
       eps_gmres,niter,errs,rres,sigma,contrast_matrix, &
       npts_vect,n_components,srcvals_extended)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)

      call cpu_time(t2)
!C$       t2 = omp_get_wtime()
      call prin2('analytic solve time=*',t2-t1,1)
 2111 continue    
      
      allocate(currents(6,npts))
      allocate(u_vect_s(3,npts),v_vect_s(3,npts))
      call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect_s, &
        v_vect_s,npts)
     
      print *, "done computing uv basis"
      do i=1,npts
        do j=1,3
          currents(j,i) = u_vect_s(j,i)*sigma(i) + v_vect_s(j,i)*sigma(i+npts)
          currents(j+3,i) = u_vect_s(j,i)*sigma(i+2*npts) + &
            v_vect_s(j,i)*sigma(i+3*npts)
        enddo
      enddo
      print *, "done computing currents"

      allocate(errp(npatches))
      errm = 0

      call surf_fun_error(12,npatches,norders,ixyzs,iptype,npts,currents, &
        wts,errp,errm)
      print *, "done estimatting error"
      print *, "errm=",errm
      
      call prin2('estimated error from tails=*',errm,1)
      allocate(errp_plot(npts))
      do i=1,npatches
        do j=ixyzs(i),ixyzs(i+1)-1
          errp_plot(j) = log(errp(i))/log(10.0d0)
        enddo
      enddo

      write(fname,'(a,i1,a,i1,a)') 'errs_iref',iref,'_norder',norders(1),'.vtk'

      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype, &
         npts,srccoefs,srcvals,errp_plot,trim(fname),'a')
      stop

 


! Here I define a set of target points to calculate the field E,H and
! plot/test accuracuy

!      x_min,x_max,y_min,y_max,z_min,z_max
!     READ***    OR USE THIS GRID
      x_min=minval(srcvals(1,:))
      y_min=minval(srcvals(2,:))
      z_min=minval(srcvals(3,:))

      x_max=maxval(srcvals(1,:))
      y_max=maxval(srcvals(2,:))
      z_max=maxval(srcvals(3,:))
      
      dxlam = (x_max-x_min)*omega/2/pi
      dylam = (y_max-y_min)*omega/2/pi
      dzlam = (z_max-z_min)*omega/2/pi
     ! 
      print *, "dxlam=",dxlam
      print *, "dylam=",dylam
      print *, "dzlam=",dzlam

      
      
      dx=x_max-x_min
      dy=y_max-y_min
      dz=z_max-z_min

      x_min=x_min-dx/10.0d0
      y_min=y_min-dy/10.0d0
      z_min=z_min-dz/10.0d0

      x_max=x_max+dx/10.0d0/sqrt(pi)
      y_max=y_max+dy/10.0d0/sqrt(pi)
      z_max=z_max+dz/10.0d0/sqrt(pi)


! this is for a volume equispaced sampling

      nx=20
      ny=20
      nz=1
      ntarg=nx*ny*nz
      allocate(targs(3,ntarg))
      allocate(location_targs(ntargs))
      allocate(E_far(3,ntarg),H_far(3,ntarg))

      icount=1
      do count1=1,nx
       do count2=1,ny
         targs(1,icount)=x_min+(x_max-x_min)*(count1-1)/(nx-1)
         targs(2,icount)=y_min+(y_max-y_min)*(count2-1)/(ny-1)
         targs(3,icount)=(z_max+z_min)/2.0d0
         icount=icount+1
       enddo
      enddo


!call find_inclusion_vect(npatches,norders,ixyzs,iptype,npts,srccoefs,&
!    &srcvals,wts,targs,ntarg,npatches_vect,n_components,sorted_vector,location_targs)

!      write (*,*) 'location_targs: ',location_targs(1),location_targs(2)

!     call evaluate_field_muller(npatches,norders,ixyzs,iptype,npts,srccoefs,&
!    &srcvals,wts,targs,ntarg,npatches_vect,n_components,sorted_vector,&
!    &contrast_matrix,exposed_surfaces,eps,zpars,sigma,E_far,H_far)
     


! Next we test the accuracy of the solution if the RHS has been produced by the
! subroutine get_rhs_em_muller_trans_testing

     call test_accuracy_em_muller(npatches,norders,ixyzs,iptype, &
      npts,srccoefs,srcvals,wts,targs,ntarg,npatches_vect, &
      n_components,sorted_vector,contrast_matrix,exposed_surfaces, &
      eps,zpars,sigma,xyz_in,vf,direction,Pol,err_est)

      print *, "err_est:",err_est

      stop
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars,& 
     &srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri,& 
       &triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,& 
     &ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,&
     &npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),&
     &nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,& 
     &iptr3,iptr4, norder,'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,&
     &npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end

      subroutine test_exterior_pt(npatches,norder,npts,srcvals,&
     &srccoefs,wts,xyzout,isout)

!
!  this subroutine tests whether the pt xyzin, is
!  in the exterior of a surface, and also estimates the error
!  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
!  centered at the interior point. Whether a point 
!  is in the interior or not is tested using Gauss' 
!  identity for the flux due to a point charge
!
!
!  input:
!    npatches - integer
!       number of patches
!    norder - integer
!       order of discretization
!    npts - integer
!       total number of discretization points on the surface
!    srccoefs - real *8 (9,npts)
!       koornwinder expansion coefficients of geometry info
!    xyzout -  real *8 (3)
!       point to be tested
!
!  output: 
!    isout - boolean
!      whether the target is in the interior or not
!

      implicit none
      integer npatches,norder,npts,npols
      real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer ipatch,j,i
      real *8 ra,ds
      logical isout

      done = 1
      pi = atan(done)*4

      npols = (norder+1)*(norder+2)/2

      zk = 0
      ra = 0

      do ipatch=1,npatches
        do j=1,npols
          i = (ipatch-1)*npols + j
          call h3d_sprime(xyzout,srcvals(1,i),dpars,zk,ipars,val)
          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-3) isout = .false.
      if(abs(ra).le.1.0d-3) isout = .true.

      return
      end
!
!
!
!
!
      subroutine get_miniwasp_ellip(n_components,abc_lens,abc_cone, &
        abc_rhabdom,abc_pig,xyz_lens,xyz_cone,xyz_rhabdom,xyz_pig, &
        iref,npatches,npatches_vect,norder,srcvals,srccoefs,npts, &
        npts_vect)
      implicit real *8 (a-h,o-z)
      real *8 abc_lens(3),abc_cone(3),abc_rhabdom(3),abc_pig(3)
      real *8 xyz_lens(3),xyz_cone(3),xyz_rhabdom(3),xyz_pig(3)
      real *8 srcvals(12,npts),srccoefs(9,npts)
      integer npatches_vect(4),npts_vect(4)
      character *100 fname
      
      istart = 1
      npols = (norder+1)*(norder+2)/2
      print *, "npols=",npols


      a = abc_lens(1)
      b = abc_lens(2)
      c = abc_lens(3)
      fname = 'lens.vtk'
      print *, istart
      print *, "norder=",norder
      print *, "iref=",iref
      print *, "istart=",istart
      call prin2('xyz_lens=*',xyz_lens,3)
      call get_ellipsoid_geom(a,b,c,xyz_lens,norder,npatches_vect(1), &
       iref,srcvals(1,istart),srccoefs(1,istart),ifplot,trim(fname))
      npts_vect(1) = npatches_vect(1)*npols
      print *, npatches_vect(1)
      if(n_components.eq.1) return

      istart = istart + npts_vect(1)

      a = abc_cone(1)
      b = abc_cone(2)
      c = abc_cone(3)
      fname = 'cone.vtk'
      print *, "istart cone=",istart
      call get_ellipsoid_geom(a,b,c,xyz_cone,norder,npatches_vect(2), &
       iref,srcvals(1,istart),srccoefs(1,istart),ifplot,trim(fname))
      npts_vect(2) = npatches_vect(2)*npols
      if(n_components.eq.2) return

      istart = istart + npts_vect(2)

      a = abc_rhabdom(1)
      b = abc_rhabdom(2)
      c = abc_rhabdom(3)
      fname = 'rhabdom.vtk'
      call get_ellipsoid_geom(a,b,c,xyz_rhabdom,norder,npatches_vect(3), &
       iref,srcvals(1,istart),srccoefs(1,istart),ifplot,trim(fname))
      npts_vect(3) = npatches_vect(3)*npols
      if(n_components.eq.3) return

      istart = istart + npts_vect(3)

      a = abc_pig(1)
      b = abc_pig(2)
      c = abc_pig(3)
      fname = 'pigment.vtk'
      call get_ellipsoid_geom(a,b,c,xyz_pig,norder,npatches_vect(4), &
       iref,srcvals(1,istart),srccoefs(1,istart),ifplot,trim(fname))
      npts_vect(4) = npatches_vect(4)*npols

      return
      end
!
!
!

      subroutine get_miniwasp_ellip_mem(n_components, &
        abc_lens,abc_cone,abc_rhabdom,abc_pig,iref,npatches_vect, &
        npatches)
      implicit real *8 (a-h,o-z)
      real *8 abc_lens(3),abc_cone(3),abc_rhabdom(3),abc_pig(3)
      integer npatches_vect(4)
      
      npatches = 0
      call get_rectparapiped_mem(abc_lens(1),abc_lens(2), &
        abc_lens(3),iref,npatches_vect(1))

      call get_rectparapiped_mem(abc_cone(1),abc_cone(2), &
        abc_cone(3),iref,npatches_vect(2))

      call get_rectparapiped_mem(abc_rhabdom(1),abc_rhabdom(2), &
        abc_rhabdom(3),iref,npatches_vect(3))

      call get_rectparapiped_mem(abc_pig(1),abc_pig(2), &
        abc_pig(3),iref,npatches_vect(4))
      do i=1,n_components
        npatches = npatches + npatches_vect(i)
      enddo

      return
      end
!
!
!
!
!

      subroutine get_miniwasp_ellip_params(abc_lens,abc_cone, &
        abc_rhabdom,abc_pig,xyz_lens,xyz_cone,xyz_rhabdom,xyz_pig)
      implicit real *8 (a-h,o-z)
      real *8 abc_lens(1:3),xyz_lens(1:3)
      real *8 abc_cone(1:3),xyz_cone(1:3)
      real *8 abc_rhabdom(1:3),xyz_rhabdom(1:3)
      real *8 abc_pig(1:3),xyz_pig(1:3)
      
      abc_lens(1) = 0.3d0
      abc_lens(2) = 0.2d0
      abc_lens(3) = 0.14d0

      xyz_lens(1:3) = 0.0d0
      xyz_lens(3) = -abc_lens(3)


      abc_cone(1) = 0.48d0
      abc_cone(2) = 0.4d0
      abc_cone(3) = 0.37d0

      xyz_cone(1:3) = 0
      xyz_cone(3) = xyz_lens(3)-abc_lens(3) - 0.2d0 - abc_cone(3)

      abc_rhabdom(1) = 0.2d0
      abc_rhabdom(2) = 0.25d0
      abc_rhabdom(3) = 5.1d0
      xyz_rhabdom(1:3) = 0
      xyz_rhabdom(3) = xyz_cone(3)-0.2d0-abc_cone(3)-abc_rhabdom(3)

      abc_pig(1:3) = 2*(abc_lens(3) + abc_cone(3) + abc_rhabdom(3) + & 
         0.8d0)/2
      xyz_pig(1:3) = 0
      xyz_pig(3) = -2*(abc_lens(3) + abc_cone(3) + abc_rhabdom(3) + &
        0.4d0)/2


      return
      end
!
!
!
!
!
      subroutine get_ellipsoid_geom(a,b,c,xyz0,norder,npatches,iref, &
        srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)

      real *8 v1(3),v2(3),v3(3),v4(3),xyz0(1:3)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      real *8, allocatable :: uvs(:,:),wts(:),umatr(:,:),vmatr(:,:)
      real *8 srcvals(12,*),srccoefs(9,*)
      character (len=*) fname
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      external xtri_ellipsoid_eval



      allocate(triaskel(3,3,npatches))

      ause = a/sqrt(3.0d0)
      buse = b/sqrt(3.0d0)
      cuse = c/sqrt(3.0d0)
      call get_rectparapiped(ause,buse,cuse,iref,npatches,triaskel)
      
      call xtri_vtk_flat(34, npatches, triaskel, 'a')

      p2(1) = a
      p2(2) = b
      p2(3) = c

      p3(1:3) = xyz0(1:3)
      
      ptr1 => triaskel(1,1,1)
      ptr2 => p2(1)
      ptr3 => p3(1)
      ptr4 => p4(1)

      norder = 5
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),wts(npols),umatr(npols,npols), &
        vmatr(npols,npols))
      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
      npts = npols*npatches
      call getgeominfo(npatches,xtri_ellipsoid_eval,ptr1,ptr2,ptr3, &
        ptr4,npols,uvs,umatr,srcvals,srccoefs)

      call xtri_vtk_surf(trim(fname),npatches,xtri_ellipsoid_eval, &
        ptr1,ptr2,ptr3,ptr4,norder,'Triangulated surface of ellipsoid')
       

      return
      end
!
!
!
!
!
!

      subroutine get_rectparapiped(a,b,c,iref,npatches,triaskel)
      implicit real *8 (a-h,o-z)

      real *8 triaskel(3,3,npatches),vs(3,4)

      real *8 vcube(3,8),xnorm(3)

      
      
      vcube(1,1) = -a
      vcube(2,1) = -b
      vcube(3,1) = -c

      vcube(1,2) = a
      vcube(2,2) = -b
      vcube(3,2) = -c

      vcube(1,3) = a
      vcube(2,3) = b
      vcube(3,3) = -c

      vcube(1,4) = -a
      vcube(2,4) = b
      vcube(3,4) = -c

      vcube(1,5) = -a
      vcube(2,5) = -b
      vcube(3,5) = c

      vcube(1,6) = a
      vcube(2,6) = -b
      vcube(3,6) = c

      vcube(1,7) = a
      vcube(2,7) = b
      vcube(3,7) = c

      vcube(1,8) = -a
      vcube(2,8) = b
      vcube(3,8) = c

      rmin = a
      if(b.lt.rmin) rmin = b
      if(c.lt.rmin) rmin = c

      na = nint(a/rmin)*2**iref*3
      nb = nint(b/rmin)*2**iref*3
      nc = nint(c/rmin)*2**iref



!       z = -c face      
      vs(1:3,1) = vcube(1:3,1)
      vs(1:3,2) = vcube(1:3,4)
      vs(1:3,3) = vcube(1:3,3)
      vs(1:3,4) = vcube(1:3,2)
      ntri = 0
      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nb,na, &
        npatches,triaskel(1,1,ntri+1))

      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call prin2('xnorm1 z=-c=*',xnorm,3)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)
      call prin2('xnorm2 z=-c=*',xnorm,3)

      ntri = ntri + 2*na*nb
      

!       z = c face      
      vs(1:3,1:4) = vcube(1:3,5:8)
      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),na,nb, &
        npatches,triaskel(1,1,ntri+1))

      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call prin2('xnorm1 z=c=*',xnorm,3)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)
      call prin2('xnorm2 z=c=*',xnorm,3)

      ntri = ntri + 2*na*nb

!      y = -b face
!
      vs(1:3,1) = vcube(1:3,1)
      vs(1:3,2) = vcube(1:3,2)
      vs(1:3,3) = vcube(1:3,6)
      vs(1:3,4) = vcube(1:3,5)

      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),na,nc, &
        npatches,triaskel(1,1,ntri+1))
      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call prin2('xnorm1 y=-b=*',xnorm,3)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)
      call prin2('xnorm2 y=-b=*',xnorm,3)


      ntri = ntri + 2*na*nc

!      y = b face
!
      vs(1:3,1) = vcube(1:3,4)
      vs(1:3,2) = vcube(1:3,8)
      vs(1:3,3) = vcube(1:3,7)
      vs(1:3,4) = vcube(1:3,3)

      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nc,na, &
        npatches,triaskel(1,1,ntri+1))
      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call prin2('xnorm1 y=b=*',xnorm,3)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)
      call prin2('xnorm2 y=b=*',xnorm,3)

      ntri = ntri + 2*na*nc



!      x = -a face
!
      vs(1:3,1) = vcube(1:3,1)
      vs(1:3,2) = vcube(1:3,5)
      vs(1:3,3) = vcube(1:3,8)
      vs(1:3,4) = vcube(1:3,4)

      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nc,nb, &
        npatches,triaskel(1,1,ntri+1))
      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call prin2('xnorm1 x=-a=*',xnorm,3)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)
      call prin2('xnorm2 x=-a=*',xnorm,3)


      ntri = ntri + 2*nb*nc

!      x = 2 face
!
      vs(1:3,1) = vcube(1:3,2)
      vs(1:3,2) = vcube(1:3,3)
      vs(1:3,3) = vcube(1:3,7)
      vs(1:3,4) = vcube(1:3,6)

      call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nb,nc, &
        npatches,triaskel(1,1,ntri+1))
      call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
      call prin2('xnorm1 x=a=*',xnorm,3)
      call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)
      call prin2('xnorm2 x=a=*',xnorm,3)

      ntri = ntri + 2*nb*nc


      return
      end


      subroutine get_norm_triaskel(tria,xnorm)
      implicit real *8 (a-h,o-z)
      real *8 tria(3,3),xnorm(3),xu(3),xv(3)
      

      xu(1:3) = tria(1:3,2) - tria(1:3,1)
      xv(1:3) = tria(1:3,3) - tria(1:3,1)

      xnorm(1:3) = 0
      call cross_prod3d(xu,xv,xnorm)


      return
      end
      


      subroutine get_rectparapiped_mem(a,b,c,iref,npatches)
      implicit real *8 (a-h,o-z)
      integer na,nb,nc

      rmin = a
      if(b.lt.rmin) rmin = b
      if(c.lt.rmin) rmin = c
      na = nint(a/rmin)*3
      nb = nint(b/rmin)*3
      nc = nint(c/rmin)
      print *, na,nb,nc

      npatches = 4*(na*nb + nb*nc + nc*na)*(4**iref)

      return
      end subroutine get_rectparapiped_mem



      subroutine xtri_rectmesh_3d(v1,v2,v3,v4,nu,nv,npatches,triaskel)
      implicit real *8 (a-h,o-z)
      real *8 triaskel(3,3,npatches),v1(3),v2(3),v3(3),v4(3)
      real *8 vl(3),vr(3),vb(3),vt(3)
      real *8 uvw1(3),uvw2(3),uvw3(3),uvw4(3)

      vl(1:3) = v4(1:3)-v1(1:3)
      vr(1:3) = v3(1:3)-v2(1:3)



      ntri = 0

      do i=1,nv
        uvw1(1:3) = v1(1:3) + (i-1)*vl(1:3)/(nv+ 0.0d0)
        uvw4(1:3) = uvw1(1:3) + vl(1:3)/(nv+0.0d0)

        vb(1:3) =  v2(1:3) + (i-1)*vr(1:3)/(nv+0.0d0) - uvw1(1:3)
        vt(1:3) = uvw1(1:3) + vb(1:3) + vr(1:3)/(nv+0.0d0) - uvw4(1:3)
        
        uvw2(1:3) = uvw1(1:3) + vb(1:3)/(nu+0.0d0)
        uvw3(1:3) = uvw4(1:3) + vt(1:3)/(nu+0.0d0)

        do j=1,nu
          call xtri_rectmesh0_3d(uvw1,uvw2,uvw3,uvw4, &
            triaskel(1,1,ntri+1))
          ntri = ntri + 2
          uvw1(1:3) = uvw2(1:3)
          uvw2(1:3) = uvw2(1:3) + vb(1:3)/(nu+0.0d0)
          uvw4(1:3) = uvw3(1:3)
          uvw3(1:3) = uvw3(1:3) + vt(1:3)/(nu+0.0d0)
        enddo
      enddo

      end subroutine xtri_rectmesh_3d







      subroutine xtri_rectmesh0_3d(v1,v2,v3,v4,triaskel)
      implicit real *8 (a-h,o-z)
      real *8 v1(3),v2(3),v3(3),v4(3),triaskel(3,3,2)

      do i=1,3
        triaskel(i,1,1) = v1(i)
        triaskel(i,2,1) = v2(i)
        triaskel(i,3,1) = v4(i)
        triaskel(i,1,2) = v3(i)
        triaskel(i,2,2) = v4(i)
        triaskel(i,3,2) = v2(i)
      enddo

      return
      end subroutine xtri_rectmesh0_3d



      subroutine xtri_vtk_flat(iw, ntri, xtri1s, title)
      implicit real *8 (a-h,o-z)
      real *8 :: xtri1s(3,3,ntri)
      character(len=*) :: title

      character(len=1024) :: filename, dataname, valsname, imgname
      character(len=1024) :: trisname, vecsname, centname
      character(len=12) :: fmt, fmt3, fmt4
      character(len=25) :: fmt2

!
!  This routine plots a sequence of FLAT TRIANGLES with surface
! color vals.
!
! Input:
!   iw - plot number, controls the filenames
!   ntri - number of flat triangles
!   xtri1s - full triangle information
!
! Output:
!   files which can be executed in matlab to plot the surface
!
!

      if (iw .lt. 10) then
        fmt = "(A4,I1,A4)"
        fmt3 = "(A8,I1,A4)"
        fmt4 = "(A5,I1,A4)"
      elseif (iw .lt. 100) then
        fmt = "(A4,I2,A4)"
        fmt3 = "(A8,I2,A4)"
        fmt4 = "(A5,I2,A4)"
      elseif (iw .lt. 1000) then
        fmt = "(A4,I3,A4)"
        fmt3 = "(A8,I3,A4)"
        fmt4 = "(A5,I3,A4)"
      elseif (iw .lt. 10000) then
        fmt = "(A4,I4,A4)"
        fmt3 = "(A8,I4,A4)"
        fmt4 = "(A5,I4,A4)"
      end if

      write(filename, fmt) 'plot', iw, '.vtk'

      !
      ! write the vtk plotting script
      !
      iunit1 = 877
      open(unit = iunit1, file=trim(filename), status='replace')
    
      write(iunit1,'(a)') "# vtk DataFile Version 3.0"
      write(iunit1,'(a)') "vtk output"
      write(iunit1,'(a)') "ASCII"
      !write(iunit1,'(a)') "DATASET POLYDATA"
      write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
      write(iunit1,'(a,i8,a)') "POINTS ", ntri*3, " float"

      fmt2 = "(E11.5,2X,E11.5,2X,E11.5)"
      do i = 1,ntri
        write(iunit1,fmt2) xtri1s(1,1,i), xtri1s(2,1,i), xtri1s(3,1,i)
        write(iunit1,fmt2) xtri1s(1,2,i), xtri1s(2,2,i), xtri1s(3,2,i)
        write(iunit1,fmt2) xtri1s(1,3,i), xtri1s(2,3,i), xtri1s(3,3,i)
      end do


      write(iunit1,'(a,i8,i8)') "CELLS ", ntri, ntri*4

      do i = 1,ntri
        i1 = 3*(i-1) + 1
        write(iunit1,'(a,i8,i8,i8)') "3 ", i1-1, i1, i1+1
      end do

      write(iunit1,'(a,i8)') "CELL_TYPES ", ntri
      do i = 1,ntri
        write(iunit1,'(a)') "5"
      end do

      write(iunit1,'(a,i8)') "POINT_DATA ", ntri*3
      write(iunit1,'(a)') "SCALARS scalars float 1"
      write(iunit1,'(a)') "LOOKUP_TABLE default"
      do i = 1,ntri
        do j = 1,3
          write(iunit1,'(E11.5)') xtri1s(3,j,i)
        end do
      end do



      write(iunit1,'(a)') ""
      write(iunit1,'(a)') ""
      write(iunit1,'(a,i8)') "CELL_DATA ", ntri
      write(iunit1,'(a)') "SCALARS scalars float 1"
      write(iunit1,'(a)') "LOOKUP_TABLE default"
      do i = 1,ntri
        write(iunit1,'(E13.5)') (xtri1s(3,1,i) + &
        xtri1s(3,2,i) + xtri1s(3,3,i))/3
      end do
    
      close(iunit1)

      return
      end subroutine xtri_vtk_flat





      subroutine xtri_ellipsoid_eval(itri, u, v, xyz, dxyzduv, & 
          triainfo,p2, p3, p4)
      implicit real *8 (a-h,o-z)
      real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*),p2(3),p3(3)


      !
      ! project the triangle itri in triainfo onto the sphere
      !
      !    Input:
      ! itri - triangle number to map
      ! u,v - local uv coordinates on triangle itri
      ! triainfo - flat skeleton triangle info
      ! p2,p3,p4 - dummy parameters
      !
      !    Output:
      ! xyz - point on the sphere
      ! dxyzduv - first derivative information
      !
      !

      x0=triainfo(1,1,itri)
      y0=triainfo(2,1,itri)
      z0=triainfo(3,1,itri)

      x1=triainfo(1,2,itri)
      y1=triainfo(2,2,itri)
      z1=triainfo(3,2,itri)

      x2=triainfo(1,3,itri)
      y2=triainfo(2,3,itri)
      z2=triainfo(3,3,itri)

      !
      ! ... process the geometry, return the point location on the sphere
      ! and the derivatives with respect to u and v
      !
      x=x0+u*(x1-x0)+v*(x2-x0)
      y=y0+u*(y1-y0)+v*(y2-y0)
      z=z0+u*(z1-z0)+v*(z2-z0)

      dxdu = x1-x0
      dydu = y1-y0
      dzdu = z1-z0
    
      dxdv = x2-x0
      dydv = y2-y0
      dzdv = z2-z0

      !
      ! second derivatives are zero...
      !

      !
      ! project onto the sphere
      !
      r=sqrt(x**2+y**2+z**2)
      xyz(1)=p2(1)*x/r + p3(1)
      xyz(2)=p2(2)*y/r + p3(2)
      xyz(3)=p2(3)*z/r + p3(3)

      a = x0*(x1-x0) + y0*(y1-y0) + z0*(z1-z0)
      b = (x1-x0)*(x2-x0) + (y1-y0)*(y2-y0) + (z1-z0)*(z2-z0)
      c = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0)

      drdu = (a + v*b + u*c)/r
      drdu2 = (r*c - r*drdu*drdu)/r/r

      e = x0*(x2-x0) + y0*(y2-y0) + z0*(z2-z0)
      f = b
      g = (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0)

      drdv = (e + u*f + v*g)/r
      drdv2 = (r*g - r*drdv*drdv)/r/r

      drduv = (r*b - r*drdu*drdv)/r/r

      ! du
      dxyzduv(1,1) = p2(1)*(r*dxdu-x*drdu)/r/r
      dxyzduv(2,1) = p2(2)*(r*dydu-y*drdu)/r/r
      dxyzduv(3,1) = p2(3)*(r*dzdu-z*drdu)/r/r

      ! dv
      dxyzduv(1,2) = p2(1)*(r*dxdv-x*drdv)/r/r
      dxyzduv(2,2) = p2(2)*(r*dydv-y*drdv)/r/r
      dxyzduv(3,2) = p2(3)*(r*dzdv-z*drdv)/r/r

      return
      end subroutine xtri_ellipsoid_eval

!
!
!
!
!

      subroutine est_ppw(ncomp,npatches,norders,ixyzs,iptype,npts, &
         srccoefs,srcvals,contrast_matrix,omega,ppwmin,ppwp)
      implicit real *8 (a-h,o-z)
      integer ncomp,npatches,norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches),npts
      real *8 srccoefs(9,npts),srcvals(12,npts)
      complex *16 contrast_matrix(4,ncomp),zk
      real *8 rzkmax
      real *8 omega
      real *8 ppwmin,ppwp(npatches)
      real *8, allocatable :: cms(:,:),rads(:)

      done = 1.0d0
      pi = atan(done)*4.0d0

      allocate(cms(3,npatches),rads(npatches))

      rzkmax = 0
      do i=1,ncomp
        zk = omega*sqrt(contrast_matrix(1,i))*   &
          sqrt(contrast_matrix(2,i))
        if(abs(zk).gt.rzkmax) rzkmax = abs(zk)

        zk = omega*sqrt(contrast_matrix(3,i))*   &
          sqrt(contrast_matrix(4,i))
        if(abs(zk).gt.rzkmax) rzkmax = abs(zk)
      enddo

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
          srccoefs,cms,rads)
      
      do i=1,npatches
        ppwp(i) = norders(i)*pi/(rads(i)*rzkmax)
      enddo

      ppwmin = minval(ppwp)

      return
      end subroutine est_ppw






