      implicit real *8 (a-h,o-z)
      
      real *8, allocatable :: srcvals_rhab(:,:),srccoefs_rhab(:,:)
      integer, allocatable :: iptype_rhab(:),norders_rhab(:), &
        ixyzs_rhab(:)

      real *8, allocatable :: srcvals_cone(:,:),srccoefs_cone(:,:)
      integer, allocatable :: iptype_cone(:),norders_cone(:), &
        ixyzs_cone(:)

      real *8, allocatable :: srcvals_lens(:,:),srccoefs_lens(:,:)
      integer, allocatable :: iptype_lens(:),norders_lens(:), &
        ixyzs_lens(:)
      real *8, allocatable :: srcvals_pigm(:,:),srccoefs_pigm(:,:)
      integer, allocatable :: iptype_pigm(:),norders_pigm(:), & 
        ixyzs_pigm(:)
      
      real *8, allocatable :: srcvals_extended(:,:)
      real *8 xyzc_pigm(3)
      integer ppw,ppw0

      complex *16 contrast_matrix(4,4)
      integer, allocatable :: sorted_vector(:)
      logical, allocatable :: exposed_surfaces(:)
      complex *16 pol(2),vf(3),omega,omega_use
      real *8 direction(2)

      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: dpot(:),dsigma(:)
      real *8, allocatable :: errp(:),errp_plot(:)
      real *8 errm
      real *8 dpars(3)
      real *8 xyz0(3)
      character *100 fname
      integer ipars(2),npatches_vect(4),npts_vect(4)
      integer iinc_include(4)
      complex *16, allocatable :: rhs(:),sigma(:),pot(:),potex(:)
      complex *16, allocatable :: Evals(:,:),Evalsex(:,:), &
        Hvals(:,:),Hvalsex(:,:)

      complex *16 zk,ima,zpars(5)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ipatch_id(:),inode_id(:),ipatch_id_in(:)
      real *8, allocatable :: uvs_targ(:,:),uvs_targ_in(:,:)
      real *8, allocatable :: xyz_in(:,:)
      complex *16, allocatable :: charges_in(:)
      real *8, allocatable :: rzbdry(:,:)
      real *8 xyz_out(3),errs(1000)
      complex *16 pottmp,zpig
      character (len=2) :: arg_comm
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1.0d0
      pi = atan(done)*4

      
      

      rpigm = 9094.9269d0*1.25d0
      zpigm = 8844.9269d0

      xyzc_pigm(1:3) = 0
      xyzc_pigm(3) = zpigm

      n_components = 4

! 
!   which interior point to use to generate exterior
!   field, must be between 1-150
!  
      isrcin = 34

! default values
      ppw = 5
      norder = 5
      ibc = 0
      irlam = 4
      iref = 1


      call get_command_argument(1,arg_comm)
      read(arg_comm,*) irlam

      call get_command_argument(2,arg_comm)
      read(arg_comm,*) iref

      call get_command_argument(3,arg_comm)
      read(arg_comm,*) norder

      call get_command_argument(4,arg_comm)
      read(arg_comm,*) ibc 

      if(irlam.eq.1) rlam = 550.0d0
      if(irlam.eq.2) rlam = 450.0d0
      if(irlam.eq.3) rlam = 350.0d0
      if(irlam.eq.4) rlam = 550.0d0*4

      omega = 2.0d0*pi/rlam

      idis = 0
      if(idis.eq.0) zpig = 1.34d0
      if(idis.eq.1) zpig = 1.34d0 + 0.05d0*ima
      if(idis.eq.1) zpig = 1.34d0 + 0.005d0*ima
!      zpig = 1.13d0

      eps = 1.0d-7
      eps_gmres = 1.0d-9
      numit = 400

      direction(1) = 0
      direction(2) = pi/2

      ipolcase = 1
      pol(1) = 3.0d0
      pol(2) = 3.0d0
      if(ipolcase.eq.1) then
        pol(1) = 1.0d0
        pol(2) = 0.0d0
      else if(ipolcase.eq.2) then
        pol(1) = 0.0d0
        pol(2) = 1.0d0
      endif

      rsc = 1.0d0
      iinc_include(1) = 1
      iinc_include(2) = 0
      iinc_include(3) = 0
      iinc_include(4) = 1

      n_components = iinc_include(1) + iinc_include(2) + &
        iinc_include(3) + iinc_include(4)




      
      contrast_matrix(1,1)=zpig**2
      contrast_matrix(2,1)=1.0d0  
      contrast_matrix(3,1)=1.452d0**2  
!      contrast_matrix(3,1)=1.1d0**2  
      contrast_matrix(4,1)=1.0d0

      contrast_matrix(1,2)=zpig**2
      contrast_matrix(2,2)=1.0d0  
      contrast_matrix(3,2)=1.348d0**2  
!      contrast_matrix(3,2)=1.09d0**2  
      contrast_matrix(4,2)=1.0d0

      contrast_matrix(1,3)=zpig**2
      contrast_matrix(2,3)=1.0d0  
      contrast_matrix(3,3)=1.363d0**2  
!      contrast_matrix(3,3)=1.11d0**2  
      contrast_matrix(4,3)=1.0d0

      contrast_matrix(1,4) = 1.0d0
      contrast_matrix(2,4) = 1.0d0
      contrast_matrix(3,4) = zpig**2 
      contrast_matrix(4,4) = 1.0d0

      rlam0 = rlam*1.452d0


      call get_axissym_miniwasp_geom_mem(rlam0,ppw,iref,norder, &
        npatches_rhab, &
        npts_rhab, npatches_cone, npts_cone, &
        npatches_lens, npts_lens)
      
      ppw0 = ppw*iref
      call get_polar_uvmem_nc(rpigm,rpigm,rpigm,rlam0,ppw0,norder, &
        npatches_pigm)
      npts_pigm = npatches_pigm*(norder+1)*(norder+2)/2


      allocate(srcvals_rhab(12,npts_rhab),srccoefs_rhab(9,npts_rhab))
      allocate(iptype_rhab(npatches_rhab),ixyzs_rhab(npatches_rhab+1))
      allocate(norders_rhab(npatches_rhab))

      allocate(srcvals_cone(12,npts_cone),srccoefs_cone(9,npts_cone))
      allocate(iptype_cone(npatches_cone),ixyzs_cone(npatches_cone+1))
      allocate(norders_cone(npatches_cone))

      allocate(srcvals_lens(12,npts_lens),srccoefs_lens(9,npts_lens))
      allocate(iptype_lens(npatches_lens),ixyzs_lens(npatches_lens+1))
      allocate(norders_lens(npatches_lens))

      allocate(srcvals_pigm(12,npts_pigm),srccoefs_pigm(9,npts_pigm))
      allocate(iptype_pigm(npatches_pigm),ixyzs_pigm(npatches_pigm+1))
      allocate(norders_pigm(npatches_pigm))
      
      ifplot = 0

      call get_axissym_miniwasp_geom(rlam0,ppw,iref,norder, &
        npatches_rhab, &
        norders_rhab, ixyzs_rhab, iptype_rhab, npts_rhab, &
        srcvals_rhab, srccoefs_rhab, npatches_cone, &
        norders_cone, ixyzs_cone, iptype_cone, npts_cone, &
        srcvals_cone, srccoefs_cone, npatches_lens, &
        norders_lens, ixyzs_lens, iptype_lens, npts_lens, &
        srcvals_lens, srccoefs_lens,ifplot)

      call get_ellipsoid_geom_polar_nc(rpigm,rpigm,rpigm,xyzc_pigm, &
        rlam0,ppw0,norder,npatches_pigm,norders_pigm,ixyzs_pigm, &
        iptype_pigm,npts_pigm,srcvals_pigm,srccoefs_pigm)

      if(ifplot.eq.1) then
        call surf_vtk_plot(npatches_pigm,norders_pigm,ixyzs_pigm, &
         iptype_pigm,npts_pigm,srccoefs_pigm, srcvals_pigm, &
         'pigm_axissym.vtk','a')
      endif

      npatches = 0
      if(iinc_include(1).eq.1) npatches = npatches + npatches_rhab
      if(iinc_include(2).eq.1) npatches = npatches + npatches_cone
      if(iinc_include(3).eq.1) npatches = npatches + npatches_lens
      if(iinc_include(4).eq.1) npatches = npatches + npatches_pigm


!     
!      npatches = npatches_rhab
      npols = (norder+1)*(norder+2)/2
      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(iptype(npatches),ixyzs(npatches+1),norders(npatches))

!
!   combine all four components of geometry
!
!
      istart = 0

      if(iinc_include(1).eq.1) then
        do i=1,npts_rhab
          srcvals(1:12,istart+i) = srcvals_rhab(1:12,i)
          srccoefs(1:9,istart+i) = srccoefs_rhab(1:9,i)
        enddo
        istart = istart + npts_rhab
      endif

      if(iinc_include(2).eq.1) then
        do i=1,npts_cone
          srcvals(1:12,istart+i) = srcvals_cone(1:12,i)
          srccoefs(1:9,istart+i) = srccoefs_cone(1:9,i)
        enddo
        istart = istart + npts_cone
      endif

      if(iinc_include(3).eq.1) then
        do i=1,npts_lens
          srcvals(1:12,istart+i) = srcvals_lens(1:12,i)
          srccoefs(1:9,istart+i) = srccoefs_lens(1:9,i)
        enddo
        istart = istart + npts_lens
      endif

      if(iinc_include(4).eq.1) then
        do i=1,npts_pigm
          srcvals(1:12,istart+i) = srcvals_pigm(1:12,i)
          srccoefs(1:9,istart+i) = srccoefs_pigm(1:9,i)
        enddo
      endif
 

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = (i-1)*npols+1
        iptype(i) = 1
      enddo
      ixyzs(npatches+1) = npts+1
      

      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)
      nt = 0
      call get_bsize(12,npts,srcvals,3,nt,tmp,bsize)
      rsc = 1.0d0/bsize
      rsc = 1.0d0

      do i=1,npts
        do j=1,9
          srcvals(j,i) = srcvals(j,i)*rsc
          srccoefs(j,i) = srccoefs(j,i)*rsc
        enddo
        wts(i) = wts(i)*rsc**2
      enddo
      omega_use = omega/rsc

      allocate(sorted_vector(n_components+1), &
        exposed_surfaces(n_components))
      allocate(srcvals_extended(20,npts))

      icomp = 1

      if(iinc_include(1).eq.1) then
        npatches_vect(icomp) = npatches_rhab
        npts_vect(icomp) = npts_rhab
        icomp = icomp + 1
      endif

      if(iinc_include(2).eq.1) then
        npatches_vect(icomp) = npatches_cone
        npts_vect(icomp) = npts_cone
        icomp = icomp + 1
      endif

      if(iinc_include(3).eq.1) then
        npatches_vect(icomp) = npatches_lens
        npts_vect(icomp) = npts_lens
        icomp = icomp + 1
      endif

      if(iinc_include(4).eq.1) then
        npatches_vect(icomp) = npatches_pigm
        npts_vect(icomp) = npts_pigm
        icomp = icomp + 1
      endif



      print *, "eps=",eps
      print *, "starting topological sort"

      print *, "n_components=",n_components


      call topological_sorting(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,wts,npatches_vect,n_components,sorted_vector,&
        exposed_surfaces,eps)
!      sorted_vector(1) = 1      
!      sorted_vector(2) = 2      
!      sorted_vector(3) = 3      
!      sorted_vector(4) = 4      
!      sorted_vector(5) = 0     
!      exposed_surfaces(1:4) = .false.
!      exposed_surfaces(4) = .true.
      print *, "done with topological sort"


      call build_extended_targ(n_components,srcvals_extended, &
        srcvals,npts_vect,contrast_matrix,npts)


      print *, "npatches_use=",npatches

      allocate(sigma(4*npts),rhs(4*npts))
      vf(1) = 1.0d-3
      vf(2) = 2.0d-3
      vf(3) = 3.0d-3
!
!  get boundary definition at various z slices
!
!

      nz = 5
      nskip = 5
      allocate(rzbdry(2,nz))
      call get_rzbdry(nz,nskip,rzbdry)
      call prin2('rzbdry=*',rzbdry,2*nz)

      nrin = 3
      ntin = 10


!
!  get points in the interior for analytic solution test
!
!
      nin= nz*nrin*ntin
      allocate(xyz_in(3,nin),charges_in(nin))
      allocate(dpot(nin),ipatch_id_in(nin),uvs_targ_in(2,nin))
      allocate(dsigma(npts))
      do i=1,npts
        dsigma(i) = 1.0d0
      enddo
      ii = 0
      do iz=1,nz
        z0 = rzbdry(2,iz)
        r0 = rzbdry(1,iz)
        do ir=1,nrin
          rr = (ir-0.5d0)/(nrin+0.0d0)*0.5d0*r0
          do itt = 1,ntin
            thet = (itt-1.0d0)/(ntin+0.0d0)*2*pi
            ii = ii +1
            xyz_in(1,ii) = rr*cos(thet)*rsc
            xyz_in(2,ii) = rr*sin(thet)*rsc
            xyz_in(3,ii) = z0*rsc
            charges_in(ii) = (hkrand(0)-0.5d0) + ima*(hkrand(0)-0.5d0)
            ipatch_id_in(ii) = -1
            uvs_targ_in(1:2,ii) = 0
          enddo
        enddo
      enddo
      print *, "isrcin=",isrcin
      print *, "xyz_in=",xyz_in(1:3,isrcin)
      print *, "n_components=",n_components
      print *, "omega_use=",omega_use
      print *, "omega=",omega
      print *, "rsc=",rsc

      print *, "vf =",vf
      print *, "direction=",direction
      print *, "pol=",pol

      rhs = 0
      print *, "npts_vect=",npts_vect(1:n_components)
      


      call get_rhs_em_muller_trans_testing(xyz_in(1,isrcin), &
        vf,direction,pol,srcvals,omega_use,rhs, n_components, npts_vect, &
        contrast_matrix, npts, exposed_surfaces)

      zpars(1) = omega_use
      call prin2('contrast_matrix=*',contrast_matrix,32)
      sigma = 0

!      goto 1111 
      call em_muller_trans_v2_solver(npatches,norders,ixyzs,iptype, &
        npts,srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres, &
        niter,errs,rres,sigma,contrast_matrix,npts_vect,n_components, &
        srcvals_extended)
      call prin2('errs=*',errs,niter)

 1111 continue 

!      do i=1,4*npts
!        sigma(i) = hkrand(0) + ima*hkrand(0)
!      enddo
 

      nt = 10
      nr = 10

      ntarg = nt*nz*nr
      allocate(targs(3,ntarg))
      allocate(pot(ntarg),potex(ntarg))
      allocate(Evals(3,ntarg),Evalsex(3,ntarg))
      allocate(Hvals(3,ntarg),Hvalsex(3,ntarg))

      ndtarg = 3
      ii = 0
      do iz=1,nz
        z0 = rzbdry(2,iz)
        r0 = rzbdry(1,iz)
        do ir=1,nr
          rexp = 0 + (ir-1.0d0)/(nr-1.0d0)*0.7
          r =(1.0d0+10**(rexp))
          rr = r0*r
          do itt = 1,nt
            thet = (itt-1.0d0)/(ntin+0.0d0)*2*pi
            ii = ii +1
            targs(1,ii) = rr*cos(thet)*rsc
            targs(2,ii) = rr*sin(thet)*rsc
            targs(3,ii) = z0*rsc
            pot(ii) = 0
            potex(ii) = 0
            Evals(1:3,ii) = 0
            Evalsex(1:3,ii) = 0
            Hvals(1:3,ii) = 0
            Hvalsex(1:3,ii) = 0
          enddo
        enddo
      enddo
      

      allocate(ipatch_id(ntarg),uvs_targ(2,ntarg))
      do i=1,ntarg
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo

      err_est = 0
      call test_accuracy_em_muller_wfield(npatches,norders,ixyzs, &
       iptype,npts,srccoefs,srcvals,wts,targs,ntarg,npatches_vect, &
       n_components, sorted_vector, contrast_matrix, &
       exposed_surfaces, eps, zpars, sigma, xyz_in(1,isrcin), vf, &
       direction, pol, Evals,Evalsex, Hvals, Hvalsex,err_est)
      print *, "err_est=",err_est

      open(unit=33,file='soln.dat')
      do i=1,4*npts
        write(33,*) real(sigma(i)),imag(sigma(i))
      enddo
      close(33)

      open(unit=33,file='rhs.dat')
      do i=1,4*npts
        write(33,*) real(rhs(i)),imag(rhs(i)) 
      enddo
      close(33)

      open(unit=33,file='targinfo.dat')
      do i=1,ntarg
        write(33,*) targs(1,i),targs(2,i),targs(3,i)
      enddo
      close(33)

      open(unit=33,file='Evals.dat')
      do i=1,ntarg
        write(33,'(6(2x,e11.5))') real(Evals(1,i)), &
          imag(Evals(1,i)),real(Evals(2,i)), imag(Evals(2,i)), &
          real(Evals(3,i)),imag(Evals(3,i))
      enddo
      close(33)
      
      open(unit=33,file='Hvals.dat')
      do i=1,ntarg
        write(33,'(6(2x,e11.5))') real(Hvals(1,i)), &
          imag(Hvals(1,i)),real(Hvals(2,i)), imag(Hvals(2,i)), &
          real(Hvals(3,i)),imag(Hvals(3,i))
      enddo
      close(33)
 
      open(unit=33,file='Evalsex.dat')
      do i=1,ntarg
        write(33,'(6(2x,e11.5))') real(Evalsex(1,i)), &
          imag(Evalsex(1,i)),real(Evalsex(2,i)), imag(Evalsex(2,i)), &
          real(Evalsex(3,i)),imag(Evalsex(3,i))
      enddo
      close(33)
      
      open(unit=33,file='Hvalsex.dat')
      do i=1,ntarg
        write(33,'(6(2x,e11.5))') real(Hvalsex(1,i)), &
          imag(Hvalsex(1,i)),real(Hvalsex(2,i)), imag(Hvalsex(2,i)), &
          real(Hvalsex(3,i)),imag(Hvalsex(3,i))
      enddo
      close(33)

      print *, "estimated error=",err_est


      stop
      end
!
!
!
!
      
      
      subroutine get_axissym_miniwasp_geom_mem(rlam,ppw,iref,norder, &
        npatches_rhab, &
        npts_rhab, npatches_cone, npts_cone, &
        npatches_lens, npts_lens)
      
      implicit real *8 (a-h,o-z)
      real *8 verts_rhab(2,4),verts_cone(2,6), verts_lens(2,5)
      real *8 widths_rhab(4),widths_cone(6),widths_lens(5)
      real *8, allocatable :: chunks_rhab(:,:,:),ders_rhab(:,:,:)
      real *8, allocatable :: ders2_rhab(:,:,:),hs_rhab(:)
      integer, allocatable :: adjs_rhab(:,:)

      real *8, allocatable :: chunks_cone(:,:,:),ders_cone(:,:,:)
      real *8, allocatable :: ders2_cone(:,:,:),hs_cone(:)
      integer, allocatable :: adjs_cone(:,:)

      real *8, allocatable :: chunks_lens(:,:,:),ders_lens(:,:,:)
      real *8, allocatable :: ders2_lens(:,:,:),hs_lens(:)
      integer, allocatable :: adjs_lens(:,:)

      real *8, allocatable :: rzvals_rhab(:,:,:),rzcoeffs_rhab(:,:,:)
      real *8, allocatable :: whts_rhab(:,:)

      real *8, allocatable :: rzvals_cone(:,:,:),rzcoeffs_cone(:,:,:)
      real *8, allocatable :: whts_cone(:,:)

      real *8, allocatable :: rzvals_lens(:,:,:),rzcoeffs_lens(:,:,:)
      real *8, allocatable :: whts_lens(:,:)
      integer ppw


      character (len=1000) fname1,fname2,fname3

      done = 1.0d0
      pi = atan(done)*4

      call prini(6,13)

      nchmax = 10000
      k = 16
      allocate(chunks_rhab(2,k,nchmax),ders_rhab(2,k,nchmax))
      allocate(ders2_rhab(2,k,nchmax),hs_rhab(nchmax))
      allocate(adjs_rhab(2,nchmax))

      allocate(chunks_cone(2,k,nchmax),ders_cone(2,k,nchmax))
      allocate(ders2_cone(2,k,nchmax),hs_cone(nchmax))
      allocate(adjs_cone(2,nchmax))

      allocate(chunks_lens(2,k,nchmax),ders_lens(2,k,nchmax))
      allocate(ders2_lens(2,k,nchmax),hs_lens(nchmax))
      allocate(adjs_lens(2,nchmax))

!
!
!  read in vertices
!
      

      nv_rhab = 4
      nv_cone = 6
      nv_lens = 5
      call get_verts_mwaxissym(verts_rhab,nv_rhab,verts_cone,nv_cone, &
        verts_lens, nv_lens)
      call prin2('verts_rhab=*',verts_rhab,2*nv_rhab)
      call prin2('verts_rhab=*',verts_cone,2*nv_cone)
      call prin2('verts_rhab=*',verts_lens,2*nv_lens)


!
!  Compute widths for smoothing
!
!
      rfac = 3.0d0
      widths_rhab(1) = 0
      widths_rhab(nv_rhab) = 0

      widths_cone(1) = 0
      widths_cone(nv_cone) = 0

      widths_lens(1) = 0
      widths_lens(nv_lens) = 0

      do i=2,nv_rhab-1
        r1 = sqrt((verts_rhab(1,i)-verts_rhab(1,i-1))**2 + &
                  (verts_rhab(2,i)-verts_rhab(2,i-1))**2)
        
        r2 = sqrt((verts_rhab(1,i)-verts_rhab(1,i+1))**2 + &
                  (verts_rhab(2,i)-verts_rhab(2,i+1))**2)
        ruse = r1
        if(r2.le.ruse) ruse = r2
        widths_rhab(i) = ruse/rfac
      enddo


      do i=2,nv_cone-1
        r1 = sqrt((verts_cone(1,i)-verts_cone(1,i-1))**2 + &
                  (verts_cone(2,i)-verts_cone(2,i-1))**2)
        
        r2 = sqrt((verts_cone(1,i)-verts_cone(1,i+1))**2 + &
                  (verts_cone(2,i)-verts_cone(2,i+1))**2)
        ruse = r1
        if(r2.le.ruse) ruse = r2
        widths_cone(i) = ruse/rfac
      enddo



      do i=2,nv_lens-1
        r1 = sqrt((verts_lens(1,i)-verts_lens(1,i-1))**2 + &
                  (verts_lens(2,i)-verts_lens(2,i-1))**2)
        
        r2 = sqrt((verts_lens(1,i)-verts_lens(1,i+1))**2 + &
                  (verts_lens(2,i)-verts_lens(2,i+1))**2)
        ruse = r1
        if(r2.le.ruse) ruse = r2
        widths_lens(i) = ruse/rfac
      enddo
!
!  Chunk axissymmetric parts of the curve
!
!

      nch_rhab = 0
      nch_cone = 0
      nch_lens = 0

      ibell = 1
      nover = 0
      
      eps = 1.0d-10
      i1 = 0
      i2 = 0
      p1 = 0.0d0
      p2 = 0.0d0
      ifclosed = 0


      ier_rhab = 0
      ier_cone = 0 
      ier_lens = 0

      call prin2('widths_rhab=*',widths_rhab,nv_rhab)
      call prin2('verts_rhab=*',verts_rhab,2*nv_rhab)

      call chunkpolysmooth(ier_rhab,eps,widths_rhab,ibell,p1,p2,i1,i2, &
        nv_rhab,verts_rhab,ifclosed,nover,k,nch_rhab,chunks_rhab, &
        adjs_rhab,ders_rhab,ders2_rhab,hs_rhab)


      call chunkpolysmooth(ier_cone,eps,widths_cone,ibell,p1,p2,i1,i2, &
        nv_cone,verts_cone,ifclosed,nover,k,nch_cone,chunks_cone, &
        adjs_cone,ders_cone,ders2_cone,hs_cone)


      call chunkpolysmooth(ier_lens,eps,widths_lens,ibell,p1,p2,i1,i2, &
        nv_lens,verts_lens,ifclosed,nover,k,nch_lens,chunks_lens, &
        adjs_lens,ders_lens,ders2_lens,hs_lens)

      allocate(rzvals_rhab(8,k,nch_rhab),rzcoeffs_rhab(6,k,nch_rhab))
      allocate(whts_rhab(k,nch_rhab))
      allocate(rzvals_cone(8,k,nch_cone),rzcoeffs_cone(6,k,nch_cone))
      allocate(whts_cone(k,nch_cone))
      allocate(rzvals_lens(8,k,nch_lens),rzcoeffs_lens(6,k,nch_lens))
      allocate(whts_lens(k,nch_lens))
!
!  convert chunks format to vals, coeffs format
!
!
      call chunksort(k,nch_rhab,chunks_rhab,adjs_rhab,ders_rhab, &
        ders2_rhab,hs_rhab)
      call chunks_to_srcinfo(k,nch_rhab,chunks_rhab,ders_rhab, &
        ders2_rhab,hs_rhab,rzvals_rhab,rzcoeffs_rhab,whts_rhab)

      call chunksort(k,nch_cone,chunks_cone,adjs_cone,ders_cone, &
        ders2_cone,hs_cone)
      call chunks_to_srcinfo(k,nch_cone,chunks_cone,ders_cone, &
        ders2_cone,hs_cone,rzvals_cone,rzcoeffs_cone,whts_cone)

      call chunksort(k,nch_lens,chunks_lens,adjs_lens,ders_lens, &
        ders2_lens,hs_lens)
      call chunks_to_srcinfo(k,nch_lens,chunks_lens,ders_lens, &
        ders2_lens,hs_lens,rzvals_lens,rzcoeffs_lens,whts_lens)

      npatches_rhab = 0
      npatches_cone = 0
      npatches_lens = 0

!
!  note extra factor of 1.4 to account for refractive index
!

      call get_axissym_uvmem(nch_rhab,k,rzcoeffs_rhab,rzvals_rhab, &
         rlam,ppw,iref,norder,npatches_rhab)

      call get_axissym_uvmem(nch_cone,k,rzcoeffs_cone,rzvals_cone, &
         rlam,ppw,iref,norder,npatches_cone)

      call get_axissym_uvmem(nch_lens,k,rzcoeffs_lens,rzvals_lens, &
         rlam,ppw,iref,norder,npatches_lens)
      

      print *, npatches_rhab,npatches_cone,npatches_lens

      npols = (norder+1)*(norder+2)/2
      npts_rhab = npatches_rhab*npols
      npts_cone = npatches_cone*npols
      npts_lens = npatches_lens*npols

      return
      end

      




      
      subroutine get_axissym_miniwasp_geom(rlam,ppw,iref,norder, &
        npatches_rhab, &
        norders_rhab, ixyzs_rhab, iptype_rhab, npts_rhab, &
        srcvals_rhab, srccoefs_rhab, npatches_cone, &
        norders_cone, ixyzs_cone, iptype_cone, npts_cone, &
        srcvals_cone, srccoefs_cone, npatches_lens, &
        norders_lens, ixyzs_lens, iptype_lens, npts_lens, &
        srcvals_lens, srccoefs_lens,ifplot)
      
      implicit real *8 (a-h,o-z)
      real *8 verts_rhab(2,4),verts_cone(2,6), verts_lens(2,5)
      real *8 widths_rhab(4),widths_cone(6),widths_lens(5)
      real *8, allocatable :: chunks_rhab(:,:,:),ders_rhab(:,:,:)
      real *8, allocatable :: ders2_rhab(:,:,:),hs_rhab(:)
      integer, allocatable :: adjs_rhab(:,:)

      real *8, allocatable :: chunks_cone(:,:,:),ders_cone(:,:,:)
      real *8, allocatable :: ders2_cone(:,:,:),hs_cone(:)
      integer, allocatable :: adjs_cone(:,:)

      real *8, allocatable :: chunks_lens(:,:,:),ders_lens(:,:,:)
      real *8, allocatable :: ders2_lens(:,:,:),hs_lens(:)
      integer, allocatable :: adjs_lens(:,:)

      real *8, allocatable :: rzvals_rhab(:,:,:),rzcoeffs_rhab(:,:,:)
      real *8, allocatable :: whts_rhab(:,:)

      real *8, allocatable :: rzvals_cone(:,:,:),rzcoeffs_cone(:,:,:)
      real *8, allocatable :: whts_cone(:,:)

      real *8, allocatable :: rzvals_lens(:,:,:),rzcoeffs_lens(:,:,:)
      real *8, allocatable :: whts_lens(:,:)


      real *8 srcvals_rhab(12,npts_rhab),srccoefs_rhab(9,npts_rhab)
      integer iptype_rhab(npatches_rhab),norders_rhab(npatches_rhab), &
        ixyzs_rhab(npatches_rhab+1)

      real *8 srcvals_cone(12,npts_cone),srccoefs_cone(9,npts_cone)
      integer iptype_cone(npatches_cone),norders_cone(npatches_cone), &
        ixyzs_cone(npatches_cone+1)

      real *8 srcvals_lens(12,npts_lens),srccoefs_lens(9,npts_lens)
      integer iptype_lens(npatches_lens),norders_lens(npatches_lens), &
        ixyzs_lens(npatches_lens+1)


      integer ppw


      character (len=1000) fname1,fname2,fname3

      done = 1.0d0
      pi = atan(done)*4

      call prini(6,13)

      nchmax = 10000
      k = 16
      allocate(chunks_rhab(2,k,nchmax),ders_rhab(2,k,nchmax))
      allocate(ders2_rhab(2,k,nchmax),hs_rhab(nchmax))
      allocate(adjs_rhab(2,nchmax))

      allocate(chunks_cone(2,k,nchmax),ders_cone(2,k,nchmax))
      allocate(ders2_cone(2,k,nchmax),hs_cone(nchmax))
      allocate(adjs_cone(2,nchmax))

      allocate(chunks_lens(2,k,nchmax),ders_lens(2,k,nchmax))
      allocate(ders2_lens(2,k,nchmax),hs_lens(nchmax))
      allocate(adjs_lens(2,nchmax))

!
!
!  read in vertices
!
      

      nv_rhab = 4
      nv_cone = 6
      nv_lens = 5


      call get_verts_mwaxissym(verts_rhab,nv_rhab,verts_cone,nv_cone, &
        verts_lens, nv_lens)
      call prin2('verts_rhab=*',verts_rhab,2*nv_rhab)
      call prin2('verts_rhab=*',verts_cone,2*nv_cone)
      call prin2('verts_rhab=*',verts_lens,2*nv_lens)

!
!  Compute widths for smoothing
!
!
      rfac = 3.0d0
      widths_rhab(1) = 0
      widths_rhab(nv_rhab) = 0

      widths_cone(1) = 0
      widths_cone(nv_cone) = 0

      widths_lens(1) = 0
      widths_lens(nv_lens) = 0

      do i=2,nv_rhab-1
        r1 = sqrt((verts_rhab(1,i)-verts_rhab(1,i-1))**2 + &
                  (verts_rhab(2,i)-verts_rhab(2,i-1))**2)
        
        r2 = sqrt((verts_rhab(1,i)-verts_rhab(1,i+1))**2 + &
                  (verts_rhab(2,i)-verts_rhab(2,i+1))**2)
        ruse = r1
        if(r2.le.ruse) ruse = r2
        widths_rhab(i) = ruse/rfac
      enddo


      do i=2,nv_cone-1
        r1 = sqrt((verts_cone(1,i)-verts_cone(1,i-1))**2 + &
                  (verts_cone(2,i)-verts_cone(2,i-1))**2)
        
        r2 = sqrt((verts_cone(1,i)-verts_cone(1,i+1))**2 + &
                  (verts_cone(2,i)-verts_cone(2,i+1))**2)
        ruse = r1
        if(r2.le.ruse) ruse = r2
        widths_cone(i) = ruse/rfac
      enddo



      do i=2,nv_lens-1
        r1 = sqrt((verts_lens(1,i)-verts_lens(1,i-1))**2 + &
                  (verts_lens(2,i)-verts_lens(2,i-1))**2)
        
        r2 = sqrt((verts_lens(1,i)-verts_lens(1,i+1))**2 + &
                  (verts_lens(2,i)-verts_lens(2,i+1))**2)
        ruse = r1
        if(r2.le.ruse) ruse = r2
        widths_lens(i) = ruse/rfac
      enddo
!
!  Chunk axissymmetric parts of the curve
!
!

      nch_rhab = 0
      nch_cone = 0
      nch_lens = 0

      ibell = 1
      nover = 0
      
      eps = 1.0d-10
      i1 = 0
      i2 = 0
      p1 = 0.0d0
      p2 = 0.0d0
      ifclosed = 0


      ier_rhab = 0
      ier_cone = 0 
      ier_lens = 0

      call prin2('widths_rhab=*',widths_rhab,nv_rhab)
      call prin2('verts_rhab=*',verts_rhab,2*nv_rhab)

      call chunkpolysmooth(ier_rhab,eps,widths_rhab,ibell,p1,p2,i1,i2, &
        nv_rhab,verts_rhab,ifclosed,nover,k,nch_rhab,chunks_rhab, &
        adjs_rhab,ders_rhab,ders2_rhab,hs_rhab)


      call chunkpolysmooth(ier_cone,eps,widths_cone,ibell,p1,p2,i1,i2, &
        nv_cone,verts_cone,ifclosed,nover,k,nch_cone,chunks_cone, &
        adjs_cone,ders_cone,ders2_cone,hs_cone)


      call chunkpolysmooth(ier_lens,eps,widths_lens,ibell,p1,p2,i1,i2, &
        nv_lens,verts_lens,ifclosed,nover,k,nch_lens,chunks_lens, &
        adjs_lens,ders_lens,ders2_lens,hs_lens)

      print *, nch_rhab,nch_cone,nch_lens
      if(ifplot.eq.1) then
        open(unit=33,file='plottmp.dat')
        do i=1,nch_rhab
          do j=1,k
            write(33,*) chunks_rhab(1,j,i),chunks_rhab(2,j,i)
          enddo
        enddo

        do i=1,nch_cone
          do j=1,k
            write(33,*) chunks_cone(1,j,i),chunks_cone(2,j,i)
          enddo
        enddo

        do i=1,nch_lens
          do j=1,k
            write(33,*) chunks_lens(1,j,i),chunks_lens(2,j,i)
          enddo
        enddo
        close(33)
      endif

      allocate(rzvals_rhab(8,k,nch_rhab),rzcoeffs_rhab(6,k,nch_rhab))
      allocate(whts_rhab(k,nch_rhab))
      allocate(rzvals_cone(8,k,nch_cone),rzcoeffs_cone(6,k,nch_cone))
      allocate(whts_cone(k,nch_cone))
      allocate(rzvals_lens(8,k,nch_lens),rzcoeffs_lens(6,k,nch_lens))
      allocate(whts_lens(k,nch_lens))
!
!  convert chunks format to vals, coeffs format
!
!
      call chunksort(k,nch_rhab,chunks_rhab,adjs_rhab,ders_rhab, &
        ders2_rhab,hs_rhab)
      call chunks_to_srcinfo(k,nch_rhab,chunks_rhab,ders_rhab, &
        ders2_rhab,hs_rhab,rzvals_rhab,rzcoeffs_rhab,whts_rhab)

      call chunksort(k,nch_cone,chunks_cone,adjs_cone,ders_cone, &
        ders2_cone,hs_cone)
      call chunks_to_srcinfo(k,nch_cone,chunks_cone,ders_cone, &
        ders2_cone,hs_cone,rzvals_cone,rzcoeffs_cone,whts_cone)

      call chunksort(k,nch_lens,chunks_lens,adjs_lens,ders_lens, &
        ders2_lens,hs_lens)
      call chunks_to_srcinfo(k,nch_lens,chunks_lens,ders_lens, &
        ders2_lens,hs_lens,rzvals_lens,rzcoeffs_lens,whts_lens)

      
      call get_axissym_geom(nch_rhab,k,rzcoeffs_rhab,rzvals_rhab, &
         rlam,ppw,iref,norder,npatches_rhab,norders_rhab,ixyzs_rhab, &
         iptype_rhab,npts_rhab,srcvals_rhab,srccoefs_rhab)

      call get_axissym_geom(nch_cone,k,rzcoeffs_cone,rzvals_cone, &
         rlam,ppw,iref,norder,npatches_cone,norders_cone,ixyzs_cone, &
         iptype_cone,npts_cone,srcvals_cone,srccoefs_cone)

      call get_axissym_geom(nch_lens,k,rzcoeffs_lens,rzvals_lens, &
         rlam,ppw,iref,norder,npatches_lens,norders_lens,ixyzs_lens, &
         iptype_lens,npts_lens,srcvals_lens,srccoefs_lens)
      print *, npatches_rhab,npatches_cone,npatches_lens
      print *, npts_rhab,npts_cone,npts_lens
      print *, npts_rhab+npts_cone+npts_lens

      if(ifplot.eq.1) then

       call surf_vtk_plot(npatches_rhab,norders_rhab,ixyzs_rhab, &
         iptype_rhab,npts_rhab,srccoefs_rhab,srcvals_rhab, &
         'rhabdom_axissym.vtk','a')
       call surf_quadratic_msh_vtk_plot(npatches_rhab,norders_rhab, &
         ixyzs_rhab,iptype_rhab,npts_rhab,srccoefs_rhab,srcvals_rhab, &
         'rhabdom_axissym_msh.vtk','a')


        call surf_vtk_plot(npatches_cone,norders_cone,ixyzs_cone, &
         iptype_cone,npts_cone,srccoefs_cone, srcvals_cone, &
         'cone_axissym.vtk','a')

        call surf_vtk_plot(npatches_lens,norders_lens,ixyzs_lens, &
         iptype_lens,npts_lens,srccoefs_lens, srcvals_lens, &
         'lens_axissym.vtk','a')
       endif

      return
      end

      




      subroutine get_axissym_uvmem(nch,k,rzcoefs,rzvals,rlam,ppw,iref, &
         norder,npatches)
!
!  This subroutine estimates the number of patches required to
!  mesh an axissymetric object whose generating curve is described
!  via a chunked representation given by rzcoefs, and rzvals.
!
!  Note for developers: This is using the old chunkie representation
!  and needs to be updated to the new representation before pushing
!  changes to main repo
!
!  Input arguments:
!  
!    - nch: integer
!        number of chunks describing the generating curve
!    - k: integer
!        order of discretization of panels on the generating curve
!    - rzcoefs: real *8 (6,k,nch)
!        legendre coefficient expansions of r(s),z(s),drds,dzds,d2rds, and d2zds
!        Note that each panel is assumed to be a map from [-1,1]
!    - rzvals: real *8 (8,k,nch)
!        discretization points for r,z,drds,dzds,d2rds,d2zds, and normals
!    - rlam: real *8
!        wavelength
!    - ppw: integer
!        points per wavelength
!    - norder: integer
!        order of discretization of meshed 3d geometry
!
!  Output arguments:
!  
!    - npatches: integer
!        total number of patches
!
!

      implicit real *8 (a-h,o-z)
      real *8 rzcoefs(6,k,nch),rzvals(8,k,nch)
      real *8, allocatable :: ts(:),ws(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: amatrint(:,:),work(:)
      real *8 xs(2),rs(2)
      integer ppw

      done = 1.0d0
      pi = atan(done)*4

      itype = 1
      xs(1) = -1.0d0
      xs(2) = 1.0d0
      allocate(ts(k),ws(k),umat(k,k),vmat(k,k),amatrint(2,k))

      lw = 4*k*k + k +200
      allocate(work(lw))
      call legeexps(itype,k,ts,umat,vmat,ws)

      call lematrin(k,2,xs,amatrint,ts,work)

      iasp = 5

      npatches = 0
      do ich=1,nch
!
!  find max r on each patch, and length of each patch
!
        rmax = 0.0d0
        rlen = 0.0d0
        rs(1:2) = 0
        do j=1,k
          rs(1:2) = rs(1:2) + amatrint(1:2,j)*rzvals(1,j,ich)
          if(rzvals(1,j,ich).ge.rmax) rmax = rzvals(1,j,ich)
          dsdt = sqrt(rzvals(3,j,ich)**2 + rzvals(4,j,ich)**2)
          rlen = rlen + dsdt*ws(j)
        enddo
        if(rs(1).ge.rmax) rmax = rs(1)
        if(rs(2).ge.rmax) rmax = rs(2)

        if(rlen.ge.iasp*2*pi*rmax) then
          nt = ceiling(rmax*2*pi*(ppw+0.0d0)/rlam/(norder+0.0d0))
          ns = ceiling(nt*rlen/2/pi/rmax/(iasp+0.0d0))
        else if(2*pi*rmax.ge.iasp*rlen) then
          ns = ceiling(rlen*(ppw+0.0d0)/rlam/(norder+0.0d0))
          nt = ceiling(ns*2*pi*rmax/rlen/(iasp+0.0d0))
        else
          ns = ceiling(rlen*(ppw+0.0d0)/rlam/(norder+0.0d0))
          nt = ceiling(rmax*2*pi*(ppw+0.0d0)/rlam/(norder+0.0d0))
        endif
        ns = ns*iref
        nt = nt*iref
        npatches = npatches + 2*ns*nt
      enddo


      return
      end





      subroutine get_axissym_geom(nch,k,rzcoefs,rzvals,rlam,ppw,iref, &
         norder,npatches,norders,ixyzs,iptype,npts,srcvals,srccoefs)
!
!  This subroutine estimates the number of patches required to
!  mesh an axissymetric object whose generating curve is described
!  via a chunked representation given by rzcoefs, and rzvals.
!
!  Note for developers: This is using the old chunkie representation
!  and needs to be updated to the new representation before pushing
!  changes to main repo
!
!  Input arguments:
!  
!    - nch: integer
!        number of chunks describing the generating curve
!    - k: integer
!        order of discretization of panels on the generating curve
!    - rzcoefs: real *8 (6,k,nch)
!        legendre coefficient expansions of r(s),z(s),drds,dzds,d2rds, and d2zds
!        Note that each panel is assumed to be a map from [-1,1]
!    - rzvals: real *8 (8,k,nch)
!        discretization points for r,z,drds,dzds,d2rds,d2zds, and normals
!    - rlam: real *8
!        wavelength
!    - ppw: integer
!        points per wavelength
!    - norder: integer
!        order of discretization of meshed 3d geometry
!    - npatches: integer
!        total number of patches
!
!  Output arguments:
!    - norders: integer(npatches)
!        order of discretization of patches on axissymetric surface
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!

      implicit real *8 (a-h,o-z)
      integer :: k,nch
      integer, target :: kuse
      integer norders(npatches),ixyzs(npatches+1),iptype(npatches)
      real *8 srcvals(12,npts),srccoefs(9,npts)
      real *8 rzvals(8,k,nch)
      real *8 :: rzcoefs(6,k,nch)
      real *8, allocatable, target :: rzcoefs_tmp(:,:,:)
      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, allocatable :: ts(:),ws(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: amatrint(:,:),work(:)
      real *8, allocatable, target :: triaskel(:,:,:)
      integer, allocatable, target :: ichuse(:)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)
      real *8 xs(2),rs(2)
      integer ppw
      external xtri_axissym_chunk

      done = 1.0d0
      pi = atan(done)*4

      itype = 1
      xs(1) = -1.0d0
      xs(2) = 1.0d0
      allocate(ts(k),ws(k),umat(k,k),vmat(k,k),amatrint(2,k))

      lw = 4*k*k + k +200
      allocate(work(lw))
      call legeexps(itype,k,ts,umat,vmat,ws)

      call lematrin(k,2,xs,amatrint,ts,work)

      umin = -1.0d0
      umax = 1.0d0
      
      vmin = 2*pi
      vmax = 0

      allocate(triaskel(3,3,npatches))
      allocate(ichuse(npatches))


      npatches0 = 0
      nover = 0
      istart = 1
      iasp = 5
      do ich=1,nch
!
!  find max r on each patch, and length of each patch
!
        rmax = 0.0d0
        rlen = 0.0d0
        rs(1:2) = 0

        do j=1,k
          rs(1:2) = rs(1:2) + amatrint(1:2,j)*rzvals(1,j,ich)
          if(rzvals(1,j,ich).ge.rmax) rmax = rzvals(1,j,ich)
          dsdt = sqrt(rzvals(3,j,ich)**2 + rzvals(4,j,ich)**2)
          rlen = rlen + dsdt*ws(j)
        enddo
        if(rs(1).ge.rmax) rmax = rs(1)
        if(rs(2).ge.rmax) rmax = rs(2)
        if(rlen.ge.iasp*2*pi*rmax) then
          nt = ceiling(rmax*2*pi*(ppw+0.0d0)/rlam/(norder+0.0d0))
          ns = ceiling(nt*rlen/2/pi/rmax/(iasp+0.0d0))
        else if(2*pi*rmax.ge.iasp*rlen) then
          ns = ceiling(rlen*(ppw+0.0d0)/rlam/(norder+0.0d0))
          nt = ceiling(ns*2*pi*rmax/rlen/(iasp+0.0d0))
        else
          ns = ceiling(rlen*(ppw+0.0d0)/rlam/(norder+0.0d0))
          nt = ceiling(rmax*2*pi*(ppw+0.0d0)/rlam/(norder+0.0d0))
        endif
        ns = ns*iref
        nt = nt*iref

        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ns,nt,nover, &
          npatches,npatches0,triaskel(1,1,istart))
        ichuse(istart:(istart+2*ns*nt-1)) = ich 
        istart = istart + 2*ns*nt
      enddo

      kuse = k
      print *, "kuse=",kuse
      print *, "norder=",norder


      allocate(rzcoefs_tmp(6,k,nch))
      rzcoefs_tmp = rzcoefs

      ptr1 => triaskel(1,1,1)
      ptr2 => rzcoefs_tmp(1,1,1)
      iptr3 => ichuse(1)
      iptr4 => kuse

      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),wts(npols),umatr(npols,npols), &
        vmatr(npols,npols))
      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
      call getgeominfo(npatches,xtri_axissym_chunk,ptr1,ptr2,iptr3, &
        iptr4,npols,uvs,umatr,srcvals,srccoefs)

      do i=1,npatches
        ixyzs(i) = (i-1)*npols + 1
        iptype(i) = 1
        norders(i) = norder
      enddo
      ixyzs(npatches+1) = npts+1


      return
      end







      subroutine xtri_axissym_chunk(itri, u, v, xyz, dxyzduv, & 
          triainfo,rzcoefs,ichuse, k)
      implicit real *8 (a-h,o-z)
      real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*)
      real *8 :: rzcoefs(6,k,*),pols(k)
      integer :: ichuse(*)


      !
      ! project the triangle itri in triainfo onto the sphere
      !
      !    Input:
      ! itri - triangle number to map
      ! u,v - local uv coordinates on triangle itri
      ! triainfo - flat skeleton triangle info
      ! rzcoefs - legendre coefficient expansions of rzvals
      ! ichuse(itri) - tells which chunk or (r,z) to use for patch itri
      !
      !    Output:
      ! xyz - point on the sphere
      ! dxyzduv - first derivative information
      !
      !


      ich = ichuse(itri)
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
      s=x0+u*(x1-x0)+v*(x2-x0)
      t=y0+u*(y1-y0)+v*(y2-y0)

      dsdu = x1-x0
      dtdu = y1-y0
    
      dsdv = x2-x0
      dtdv = y2-y0
!
! s is coordinate in parametrization of r,z space
! t is coorindate in theta space
!


!
!  Compute r,z,drds,dzds
! 
!
      r = 0
      z = 0
      drds = 0
      dzds = 0
      pols = 0
      call legepols(s,k-1,pols)
      do j=1,k
        r = r + rzcoefs(1,j,ich)*pols(j)
        z = z + rzcoefs(2,j,ich)*pols(j)
        drds = drds + rzcoefs(3,j,ich)*pols(j)
        dzds = dzds + rzcoefs(4,j,ich)*pols(j)
      enddo
      !
      !
      !
      xyz(1)= r*cos(t)
      xyz(2)= r*sin(t)
      xyz(3)= z


      ! du
      dxyzduv(1,1) = drds*cos(t)*dsdu - r*sin(t)*dtdu
      dxyzduv(2,1) = drds*sin(t)*dsdu + r*cos(t)*dtdu 
      dxyzduv(3,1) = dzds*dsdu

      ! dv
      dxyzduv(1,2) = drds*cos(t)*dsdv - r*sin(t)*dtdv
      dxyzduv(2,2) = drds*sin(t)*dsdv + r*cos(t)*dtdv 
      dxyzduv(3,2) = dzds*dsdv

      return
      end subroutine xtri_axissym_chunk










!
!
!
!
!
      subroutine get_polar_uvmem_nc(a,b,c,rlam,ppw,norder,npatches)
      implicit real *8 (a-h,o-z)
      integer ppw,norder,npatches

      done = 1.0d0
      pi = atan(done)*4

      nthet = ceiling(2*c*(ppw+0.0d0)/rlam/(norder+0.0d0))
      hthet = (pi+0.0d0)/(nthet+0.0d0)

      npatches = 0
      do ithet=1,nthet
        t0 = (ithet-1)*hthet
        t1 = (ithet)*hthet

        tuse = t0
        if(abs(t0-pi/2).ge.abs(t1-pi/2)) tuse = t1
        if((t0-pi/2)*(t1-pi/2).le.0) tuse = pi/2


        alpha = a*sin(tuse)
        beta = b*sin(tuse)
        hh = (alpha-beta)**2/(alpha+beta)**2 

        ellip_p = pi*(alpha + beta)* &
           (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
        nphi = ceiling(ellip_p*(ppw+0.0d0)/rlam/(norder+0.0d0))
        npatches = npatches + 2*nphi
      enddo
      

      return
      end
!
!
!
!

      subroutine get_ellipsoid_geom_polar_nc(a,b,c,xyz0,rlam,ppw, &
        norder,npatches,norders,ixyzs,iptype,npts,srcvals,srccoefs)
      implicit real *8 (a-h,o-z)

      real *8 v1(3),v2(3),v3(3),v4(3),xyz0(1:3)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      real *8, allocatable :: uvs(:,:),wts(:),umatr(:,:),vmatr(:,:)
      real *8 srcvals(12,npts),srccoefs(9,npts)
      integer iptype(npatches),ixyzs(npatches+1),norders(npatches)
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      integer ppw
      external xtri_ellipsoid_polar_eval

      npols = (norder+1)*(norder+2)/2
      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = (i-1)*npols + 1
        iptype(i) = 1
      enddo
      ixyzs(npatches+1) = npts+1

      allocate(triaskel(3,3,npatches))

      done = 1.0d0
      pi = atan(done)*4

      nthet = ceiling(2*c*(ppw+0.0d0)/rlam/(norder+0.0d0))
      hthet = (pi+0.0d0)/(nthet+0.0d0)

      nthet0 = 1
      istart = 1

      do ithet = 1,nthet
        n0 = 0

        umin = (ithet-1)*hthet
        umax = ithet*hthet
        vmin = 0
        vmax = 2*pi

        tuse = umin
        if(abs(umin-pi/2).ge.abs(umax-pi/2)) tuse = umax
        if((umin-pi/2)*(umax-pi/2).le.0) tuse = pi/2

        nthet0 = 1
        nover = 0
        alpha = a*sin(tuse)
        beta = b*sin(tuse)
        hh = (alpha-beta)**2/(alpha+beta)**2 

        ellip_p = pi*(alpha + beta)* &
           (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
        nphi = ceiling(ellip_p*(ppw+0.0d0)/rlam/(norder+0.0d0))
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,nthet0,nphi,nover, &
          npatches,n0,triaskel(1,1,istart))
        istart = istart + n0
      enddo
      print *, npatches,n0

      p2(1) = a
      p2(2) = b
      p2(3) = c

      p3(1:3) = xyz0(1:3)
      
      ptr1 => triaskel(1,1,1)
      ptr2 => p2(1)
      ptr3 => p3(1)
      ptr4 => p4(1)

      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),wts(npols),umatr(npols,npols), &
        vmatr(npols,npols))
      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
      call getgeominfo(npatches,xtri_ellipsoid_polar_eval,ptr1,ptr2,ptr3, &
        ptr4,npols,uvs,umatr,srcvals,srccoefs)


      return
      end
!
!
!
!
!
!


      subroutine xtri_ellipsoid_polar_eval(itri, u, v, xyz, dxyzduv, & 
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
      xyz(1)=p2(1)*sin(x)*cos(y) + p3(1)
      xyz(2)=p2(2)*sin(x)*sin(y) + p3(2)
      xyz(3)=p2(3)*cos(x) + p3(3)


      ! du
      dxyzduv(1,1) = p2(1)*(cos(x)*cos(y)*dxdu - sin(x)*sin(y)*dydu)
      dxyzduv(2,1) = p2(2)*(cos(x)*sin(y)*dxdu + sin(x)*cos(y)*dydu)
      dxyzduv(3,1) = p2(3)*(-sin(x)*dxdu)

      ! dv
      dxyzduv(1,2) = p2(1)*(cos(x)*cos(y)*dxdv - sin(x)*sin(y)*dydv)
      dxyzduv(2,2) = p2(2)*(cos(x)*sin(y)*dxdv + sin(x)*cos(y)*dydv)
      dxyzduv(3,2) = p2(3)*(-sin(x)*dxdv)

      return
      end subroutine xtri_ellipsoid_polar_eval








!
      subroutine get_rzbdry(nz,nskip,rzbdry)
      implicit real *8 (a-h,o-z)
      real *8 rzbdryall(2,17),rzbdry(2,nz)
      data rzbdryall/ &
       150.779295752848,1000,&
       211.307081823032,2000,&
       271.834867893216,3000,&
       332.3626539634,4000,&
       392.890440033584,5000,&
       453.418226103768,6000,&
       513.946012173952,7000,&
       574.473798244136,8000,&
       635.00158431432,9000,&
       695.529370384504,10000,&
       756.057156454688,11000,&
       816.584942524872,12000,&
       877.112728595056,13000,&
       907.155903632639,14000,&
       1452.26059023509,15000,&
       1997.36527683753,16000,&
       1685.75362372017,17000/

       do i=1,nz
         rzbdry(1,i) = rzbdryall(1,i+nskip)
         rzbdry(2,i) = rzbdryall(2,i+nskip)
       enddo


      end subroutine get_rzbdry









      subroutine get_verts_mwaxissym(verts_rhab,nv_rhab,verts_cone, &
        nv_cone, verts_lens, nv_lens)
      implicit real *8 (a-h,o-z)
      real *8 verts_rhab(2,nv_rhab),verts_cone(2,nv_cone)
      real *8 verts_lens(2,nv_lens)

      real *8 vrhab(2,4)
      real *8 vcone(2,6)
      real *8 vlens(2,5)

      data vrhab/ &
        0,0,&
        90.251509682664,0,&
        902.51509682664,13419.6811064943,&
        0,13419.6811064943/
      data vcone/ &
        0,13660.3517989814,&
        722.012077461312,13660.3517989814,&
        2142.61551613909,16266.4630168691,&
        1787.46465646965,16917.990821341,&
        599.396363018662,16266.4630168691,&
        0,16266.4630168691/
      data vlens/&
        0,16386.7983631126,&
        541.802017255196,16386.7983631126,&
        1729.87031070618,17023.6482247334,&
        541.802017255196,17689.8539720565,&
        0,17689.8539720565/
      
      verts_rhab(1:2,1:4) = vrhab
      verts_cone(1:2,1:6) = vcone
      verts_lens(1:2,1:5) = vlens


      end subroutine get_verts_mwaxissym
