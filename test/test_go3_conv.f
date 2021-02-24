      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3),tmp(3)
      real *8 errconv(2)
      complex *16, allocatable :: sigma(:),rhs(:)
      complex *16, allocatable :: sigma2(:),rhs2(:)
      real *8, allocatable :: errs(:)
      real *8 thet,phi,eps_gmres
      complex * 16 zpars(3),zk
      integer numit,niter
      character *100 title,fname

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

c
c   set wavelength
c
c
      rlam = 10.0d0

c
c   define scaling parameter to rescale geometry. 
c

      rsc = 1.0d0
c
c
c   select geometry to test
c    icase=1, lens
c    icase=2, rhabdom
c    icase=3, cone
c
      icase = 3
c
c   test interior or exterior problem
c    ifinout = 0, interior problem
c    ifinout = 1, exterior problem
c 

      ifinout = 1

 1211 format(a,i1,a)
      do iref=0,1
        if(icase.eq.1) 
     1     write(fname,1211) '../geometries/lens_r0',iref,'.go3' 
        if(icase.eq.2) 
     1     write(fname,1211) '../geometries/rhabdom_r0',iref,'.go3' 
        if(icase.eq.3) 
     1     write(fname,1211) '../geometries/con_r0',iref,'.go3'

        print *, "fname=",trim(fname)
      
        call open_gov3_geometry_mem(fname,npatches,npts)

        allocate(srcvals(12,npts),srccoefs(9,npts))
        allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
        allocate(wts(npts))

        call prinf('npts=*',npts,1)

        call open_gov3_geometry(fname,npatches,norders,ixyzs,
     1     iptype,npts,srcvals,srccoefs,wts)
      

        allocate(sigma(npts),rhs(npts))

c
c  rescale geometry
c
        do i=1,npts
          do j=1,9
            srcvals(j,i) = srcvals(j,i)*rsc
            srccoefs(j,i) = srccoefs(j,i)*rsc
          enddo
          wts(i) = wts(i)*rsc**2
        enddo
c
c
c  compute boxsize of bounding box
c
        nt = 0 
        bsize = 0.0d0
        call get_bsize(12,npts,srcvals,3,nt,tmp,bsize)
        zk = rlam*2*pi/bsize
        call prin2('bsize=*',bsize,1)

        zpars(1) = zk
        call prin2('zk=*',zk,2)
        if(ifinout.eq.0) zpars(2) = zk*ima
        if(ifinout.eq.1) zpars(2) = -zk*ima
        zpars(3) = 1.0d0

c
c   define interior and exterior points
c
c
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

        xyz_in(1) = xa + 1.0d-3*bsize*(hkrand(0)-0.5d0) 
        xyz_in(2) = ya + 1.0d-3*bsize*(hkrand(0)-0.5d0)
        xyz_in(3) = za + 1.0d-3*bsize*(hkrand(0)-0.5d0)

        call prin2('xyz_in=*',xyz_in,3)


        xyz_out(1) = xa + 1.7d0*bsize
        xyz_out(2) = ya + 2.1d0*bsize
        xyz_out(3) = za - 3.1d0*bsize
        call prin2('xyz_out=*',xyz_out,3)

c
c  compute boundary data for analytic solution
c
c
        do i=1,npts
          if(ifinout.eq.0) 
     1       call h3d_slp(xyz_out,3,srcvals(1,i),0,dpars,1,zpars,0,
     2         ipars,rhs(i))
          if(ifinout.eq.1) 
     1       call h3d_slp(xyz_in,3,srcvals(1,i),0,dpars,1,zpars,0,
     1       ipars,rhs(i))
          rhs(i) = rhs(i)
          sigma(i) = 0
        enddo

        call prin2('rhs=*',rhs,24)
c
c
c  set gmres parameters
c

        numit = 200
        niter = 0
        allocate(errs(numit+1))

        eps = 0.51d-3
        eps_gmres = 0.51d-6

        call helm_comb_dir_solver(npatches,norders,ixyzs,iptype,npts,
     1    srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres,
     2    niter,errs,rres,sigma)

        call prinf('niter=*',niter,1)
        call prin2('rres=*',rres,1)
        call prin2('errs=*',errs,niter)
        call prin2('sigma=*',sigma,24)


c
c       test solution at appropriate point
c
        call h3d_slp(xyz_out,3,xyz_in,0,dpars,1,zpars,0,ipars,potex)
        pot = 0
          do i=1,npts
          if(ifinout.eq.0) 
     1       call h3d_comb(srcvals(1,i),3,xyz_in,0,dpars,3,zpars,1,
     2        ipars,ztmp)
          if(ifinout.eq.1) 
     1     call h3d_comb(srcvals(1,i),3,xyz_out,0,dpars,3,zpars,1,ipars,
     1     ztmp)
          pot = pot + sigma(i)*wts(i)*ztmp
        enddo

        call prin2('potex=*',potex,2)
        call prin2('pot=*',pot,2)
        erra = abs(pot-potex)/abs(potex)
        call prin2('relative error=*',erra,1)

        ra = 0
        do i=1,npts
          ra = ra + abs(sigma(i))**2*wts(i)
        enddo
        ra = sqrt(ra)
        erra2 = abs(pot-potex)/ra
        call prin2('relative error to l2 norm of sigma=*',erra2,1)
        errconv(iref+1) = erra2
        deallocate(srcvals,srccoefs,ixyzs,iptype,norders,wts)
        deallocate(errs,rhs,sigma)
      enddo
      call prin2('err conv study=*',errconv,2)

      stop
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
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
        call xtri_platonic(itype, ipars(1), npmax, ntri, 
     1      triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
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
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)

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
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         iptr3,iptr4, norder,
     2         'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end


      subroutine test_exterior_pt(npatches,norder,npts,srcvals,
     1   srccoefs,wts,xyzout,isout)
c
c
c  this subroutine tests whether the pt xyzin, is
c  in the exterior of a surface, and also estimates the error
c  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
c  centered at the interior point. Whether a point 
c  is in the interior or not is tested using Gauss' 
c  identity for the flux due to a point charge
c
c
c  input:
c    npatches - integer
c       number of patches
c    norder - integer
c       order of discretization
c    npts - integer
c       total number of discretization points on the surface
c    srccoefs - real *8 (9,npts)
c       koornwinder expansion coefficients of geometry info
c    xyzout -  real *8 (3)
c       point to be tested
c
c  output: 
c    isout - boolean
c      whether the target is in the interior or not
c

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
          call h3d_sprime(xyzout,3,srcvals(1,i),0,dpars,1,zk,0,ipars,
     1       val)
          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-3) isout = .false.
      if(abs(ra).le.1.0d-3) isout = .true.

      return
      end

   




