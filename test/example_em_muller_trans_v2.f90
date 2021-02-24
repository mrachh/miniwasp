      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:),rhs(:)
      complex *16, allocatable :: sigma2(:),rhs2(:)
	  

!      complex ( kind = 8 ), allocatable :: a_vect(:),RHS_vect(:)
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

!!!!!! You need to edit these variables to change the number of components:

      integer n_components,ntarg     !!! Number of different conected components (dielectrics)
      integer, allocatable :: npatches_vect(:),iwords(:),npts_vect(:)
      integer, allocatable :: sorted_vector(:),location_targs(:)
      complex *16, allocatable :: contrast_matrix(:,:)  !!! store ep0,mu0,ep1,mu1 on each side of each surface

      character *200 fname1    !!! name of different geometries, of each dielectric
      character *2000 string1
      
      real *8, allocatable :: srcvals_extended(:,:)

      real *8, allocatable :: dP(:,:),targs(:,:)
      logical *8, allocatable :: exposed_surfaces(:)

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
!   simulation for .go3 files of DFIE

!

! n_components is the number of different dielectrics (surfaces with contrast)

!     READ***
     n_components=1
     n_components0 = 4
	 
     allocate(npatches_vect(n_components0),npts_vect(n_components0))
     allocate(contrast_matrix(4,n_components0))
     allocate(iwords(n_components0+1))
     allocate(dP(4,n_components0))
     allocate(sorted_vector(n_components0+1))
     allocate(exposed_surfaces(n_components0))

! contrast on the two sides of surface 1.

!     READ***
     contrast_matrix(1,1)=1.1d0  !ep0
     contrast_matrix(2,1)=1.1d0  !mu0
     contrast_matrix(3,1)=1.2d0  !ep1
     contrast_matrix(4,1)=1.0d0  !mu1

! contrast on the two sides of surface 2.

     contrast_matrix(1,2)=1.1d0  !ep0
     contrast_matrix(2,2)=1.1d0  !mu0
     contrast_matrix(3,2)=1.4d0  !ep1
     contrast_matrix(4,2)=1.0d0  !mu1

! contrast on the two sides of surface 3.

     contrast_matrix(1,3)=1.2d0  !ep0
     contrast_matrix(2,3)=1.0d0  !mu0
     contrast_matrix(3,3)=1.0d0  !ep1
     contrast_matrix(4,3)=1.2d0  !mu1

! contrast on the two sides of surface 4.

     contrast_matrix(1,4)=1.0d0  !ep0
     contrast_matrix(2,4)=1.0d0  !mu0
     contrast_matrix(3,4)=1.1d0  !ep1
     contrast_matrix(4,4)=1.1d0  !mu1

     print *, contrast_matrix


!     READ***
      string1 = '../geometries/simplest_cube_quadratic_v4_o04_r01.go3?' 
!      string1 = '0_lens.msh-2.go3?' 
!      string1 = 'lens_r00.go3?' 

!      & 'simplest_cube_quadratic_v4_o04_r01.go3?' // &

!      & 'simplest_cube_quadratic_v4_o04_r01.go3?' // &

!     & '../../../../Geometries_go3/' // &
!	& 'simplest_cube_quadratic_v4_o04_r01.go3?' // &

!     & '../../../../Geometries_go3/' // &
!	& 'simplest_cube_quadratic_v4_o04_r01.go3?' // &

!     & '../../../../Geometries_go3/' // &
!	& 'simplest_cube_quadratic_v4_o04_r01.go3?' // &

!	& 'simplest_cube_quadratic_v4_o04_r01.go3?'

!  dP(:,n) defines an affine transform to be applied to the n'th .go3 geometry
!  (x,y,z) -> (x',y',z')
!
!  x'=dP(4,n)*x+dP(1,n)
!  y'=dP(4,n)*y+dP(2,n)
!  z'=dP(4,n)*z+dP(3,n)
!
!do count1=1,n_components
!   dP(:,count1) = (/0.0d0, 0.0d0, 0.0d0, 1.0d0 /)   !affine transform applied to the 1st geometry
!enddo  
!    dP(:,1) = (/100.0d0, 100.0d0, 100.0d0, 10.0d0 /)   !affine transform applied to the 1st geometry
    dP(:,1) = (/0.0d0, 0.0d0, 0.0d0, 1.0d0 /)   !affine transform applied to the 1st geometry
!    dP(:,2) = (/5.0d0, 0.0d0, -1.5d0, 1.0d0 /) !affine transform applied to the 2nd geometry
!    dP(:,3) = (/0.0d0, 0.0d0, -1.5d0*0.50d0, 0.50d0 /)  ! ...
!    dP(:,4) = (/2.5d0, 0.0d0, -1.50d0*4.0d0, 4.0d0 /)  ! ...

! I'm substracting -1.5d0*dP(4,i) in the z coordinate because the cube is originally centered at xyz=(0,0,1.5d0)
! and size 3.0d0 and I want them to be centered at z=0

    call text_process(string1,n_components,iwords)
	
	npts=0
	npatches=0
	do count1=1,n_components
		fname1=trim(string1(iwords(count1)+1:iwords(count1+1)-1))
	    call open_gov3_geometry_mem(fname1,npatches_vect(count1),npts_vect(count1))
		npts=npts+npts_vect(count1)
		npatches=npatches+npatches_vect(count1)
	        write (*,*) 'npts: ',npts
	        write (*,*) 'npatches: ',npatches
	enddo
	write (*,*) 'npts: ',npts
	write (*,*) 'npatches: ',npatches
!	write (*,*) 'npts_vect: ',npts_vect
!	read (*,*)
    allocate(srcvals_extended(20,npts))
    allocate(srcvals(12,npts),srccoefs(9,npts))
    allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
    allocate(wts(npts))

    call prinf('npatches=*',npatches,1)
    call prinf('npts=*',npts,1)
    print *, "rat=",npts/(npatches+0.0d0)
	
	do count1=1,n_components
	    fname1=trim(string1(iwords(count1)+1:iwords(count1+1)-1))
		if (count1.eq.1) then
		    i1=1
			i2=npts_vect(1)
			j1=1
			j2=npatches_vect(1)
		else
		    i1=i2+1
			i2=i1+npts_vect(count1)-1
			j1=j2+1
			j2=j1+npatches_vect(count1)-1
		endif
	    call open_gov3_geometry_v2(fname1,npatches_vect(count1),norders(j1:j2),ixyzs(j1:j2+1),&
         &iptype(j1:j2),npts_vect(count1),srcvals(:,i1:i2),srccoefs(:,i1:i2),wts(i1:i2),dP(:,count1))
		if (count1.ge.2) then
			do count2=j1,j2+1
				ixyzs(count2)=ixyzs(count2)+n_aux-1
			enddo
		endif
		n_aux=ixyzs(j2+1)
	enddo
	

!
!	Electromagnetic parameters
!
!     READ***
!     zk = omega*sqrt(epsilon*mu)  -> 2*pi/lambda_0
      omega=1.12d2*1.6d0/2
      omega=1.12d0


      zpars(1) = omega

      write (*,*) 'omega: ',omega
!      write (*,*) 'ep1: ',ep1


!      xyz_in(1) = 0.11d0
!      xyz_in(2) = 0.0d-5
!      xyz_in(3) = 0.537d0

!      xyz_out(1) = -3.5d0
!      xyz_out(2) = 3.1d0
!      xyz_out(3) = 20.1d0
      
!      fname = '../../geometries/sphere_192_o03.go3'

      allocate(sigma(4*npts),rhs(4*npts))
      ifinout = 1

      thet = hkrand(0)*pi
      phi = hkrand(0)*2*pi

!
!	Set tolerances
!
      eps = 1d-5    ! <- this eps is used in the fmm, for the quadratures and in the
                    ! sorting algorithm as tolerance to evaluate the double layer

      eps_gmres=1d-6

      call topological_sorting(npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals,&
       &wts,npatches_vect,n_components,sorted_vector,exposed_surfaces,eps)

      write (*,*) 'sorted vector: ', sorted_vector(1:2)

!
! Parameters of the electric and magnetic dipoles 
!

! Location

      xa = 0
      ya = 0
      za = 0
      ra = 0
      do i=1,npts
        xa = xa + srcvals(1,i)*wts(i)
        ya = ya + srcvals(2,i)*wts(i)
        za = za + srcvals(3,i)*wts(i)
        ra = ra + wts(i)
      enddo
      xa = xa/ra
      ya = ya/ra
      za = za/ra
      xyz_in(1) = xa    ! this point must be inside any region (not in the extarior space)
      xyz_in(2) = ya
      xyz_in(3) = za

! orientation of the electric and magnetic dipoles used to create the solution in the 
! exterior region and test the accuracy

      vf(1)=1.0d-3
      vf(2)=2.0d-3
      vf(3)=3.0d-3

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

      call get_rhs_em_muller_trans_testing(xyz_in,vf,direction,Pol,srcvals,&
       &omega,rhs,n_components,npts_vect,contrast_matrix,npts,exposed_surfaces)
!
!  This alternative RHS is used to illuminate the dielectrics with an incoming plane wave
!  defined by the angles direction(1:2), that is phi and theta; and the polarization
!  Pol(1:2), that is, Phi vector and Theta vector of the plane wave.
!

!      call get_rhs_em_muller_trans_PW(direction,Pol,srcvals,omega,rhs,&
!       &n_components,npts_vect,contrast_matrix,npts,exposed_surfaces)


      call build_extended_targ(n_components,srcvals_extended,srcvals,npts_vect,contrast_matrix,npts)
  
      numit = 400
      niter = 0
      allocate(errs(numit+1))

!
!	Specify
!
!      goto 2111
      call cpu_time(t1)
!C$      t1 = omp_get_wtime()    
      call em_muller_trans_v2_solver(npatches,norders,ixyzs,iptype,npts,&
     &srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres,niter,errs,&
     &rres,sigma,contrast_matrix,npts_vect,n_components,srcvals_extended)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)

      call cpu_time(t2)
!C$       t2 = omp_get_wtime()
      call prin2('analytic solve time=*',t2-t1,1)
 2111 continue    
!      sigma = 1
!      do i=1,4*npts
!        write(79,*) real(sigma(i)),imag(sigma(i))
!      enddo
      



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

      x_max=x_max+dx/10.0d0
      y_max=y_max+dy/10.0d0
      z_max=z_max+dz/10.0d0

      write (*,*) 'x_min,x_max,y_min,y_max,z_min,z_max',x_min,x_max,y_min,y_max,z_min,z_max

! this is for a volume equispaced sampling
!      nx=8
!      ny=8
!      nz=8
!      ntarg=nx*ny*nz
!      allocate(targs(3,ntarg))
!      allocate(location_targs(ntargs))
!      allocate(E_far(3,ntarg),H_far(3,ntarg))

!      icount=1
!      do count1=1,nx
!       do count2=1,ny
!        do count3=1,nz
!         targs(1,icount)=x_min+(x_max-x_min)*(count1-1)/(nx-1)
!         targs(2,icount)=y_min+(y_max-y_min)*(count2-1)/(ny-1)
!         targs(3,icount)=z_min+(z_max-z_min)*(count3-1)/(nz-1)
!         icount=icount+1
!        enddo
!       enddo
!      enddo

! this is for a 2D cut in the XY plane

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

      nx = 2
      ntarg = 1
      
      targs(1,1) = xa
      targs(2,1) = ya
      targs(3,1) = za
      targs(1,2) = srcvals(1,1) + 0.1*srcvals(10,1)
      targs(2,2) = srcvals(2,1) + 0.1*srcvals(11,1)
      targs(3,2) = srcvals(3,1) + 0.1*srcvals(12,1)
!      call prin2('srcvals=*',srcvals(1:3,1:6),18)
!      call prin2('targs=*',targs,3)



!call find_inclusion_vect(npatches,norders,ixyzs,iptype,npts,srccoefs,&
!    &srcvals,wts,targs,ntarg,npatches_vect,n_components,sorted_vector,location_targs)

!      write (*,*) 'location_targs: ',location_targs(1),location_targs(2)

     call prinf('ntarg=*',ntarg,1)
     call prin2('sigma=*',sigma,48)
!     call evaluate_field_muller(npatches,norders,ixyzs,iptype,npts,srccoefs,&
!    &srcvals,wts,targs,ntarg,npatches_vect,n_components,sorted_vector,&
!    &contrast_matrix,exposed_surfaces,eps,zpars,sigma,E_far,H_far)
     call prin2('sigma=*',sigma,48)
     

!    write (*,*) E_far,H_far
! do count1=1,ntarg
!   write (*,*) 'E_far',count1,E_far(:,count1)
! enddo
! do count1=1,ntarg
!   write (*,*) 'H_far',count1,H_far(:,count1)
! enddo


! Next we test the accuracy of the solution if the RHS has been produced by the
! subroutine get_rhs_em_muller_trans_testing

     call test_accuracy_em_muller(npatches,norders,ixyzs,iptype,npts,srccoefs,&
    &srcvals,wts,targs,ntarg,npatches_vect,n_components,sorted_vector,&
    &contrast_matrix,exposed_surfaces,eps,zpars,sigma,xyz_in,vf,direction,Pol,nx,ny)

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
