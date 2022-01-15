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
      integer ppw


      rlam = 550.0d0/1.4d0
      ppw = 5
      norder = 5


      call get_axissym_miniwasp_geom_mem(rlam,ppw,norder, &
        npatches_rhab, &
        npts_rhab, npatches_cone, npts_cone, &
        npatches_lens, npts_lens)

      allocate(srcvals_rhab(12,npts_rhab),srccoefs_rhab(9,npts_rhab))
      allocate(iptype_rhab(npatches_rhab),ixyzs_rhab(npatches_rhab+1))
      allocate(norders_rhab(npatches_rhab))

      allocate(srcvals_cone(12,npts_cone),srccoefs_cone(9,npts_cone))
      allocate(iptype_cone(npatches_cone),ixyzs_cone(npatches_cone+1))
      allocate(norders_cone(npatches_cone))

      allocate(srcvals_lens(12,npts_lens),srccoefs_lens(9,npts_lens))
      allocate(iptype_lens(npatches_lens),ixyzs_lens(npatches_lens+1))
      allocate(norders_lens(npatches_lens))

      ifplot = 0

      call get_axissym_miniwasp_geom(rlam,ppw,norder,npatches_rhab, &
        norders_rhab, ixyzs_rhab, iptype_rhab, npts_rhab, &
        srcvals_rhab, srccoefs_rhab, npatches_cone, &
        norders_cone, ixyzys_cone, iptype_cone, npts_cone, &
        srcvals_cone, srccoefs_cone, npatches_lens, &
        norders_lens, ixyzs_lens, iptype_lens, npts_lens, &
        srcvals_lens, srccoefs_lens,ifplot)

      stop
      end
      
      
      subroutine get_axissym_miniwasp_geom_mem(rlam,ppw,norder, &
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

      fname1 = '../geometries/rhabdom_axissym_pts.dat'
      fname2 = '../geometries/cone_axissym_pts.dat'
      fname3 = '../geometries/lens_axissym_pts.dat'

      open(unit=33,file=trim(fname1))
      do i=1,nv_rhab
        read(33,*) verts_rhab(1,i),verts_rhab(2,i)
      enddo
      close(33)

      open(unit=33,file=trim(fname2))
      do i=1,nv_cone
        read(33,*) verts_cone(1,i),verts_cone(2,i)
      enddo
      close(33)

      open(unit=33,file=trim(fname3))
      do i=1,nv_lens
        read(33,*) verts_lens(1,i),verts_lens(2,i)
      enddo
      close(33)
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
         rlam,ppw,norder,npatches_rhab)

      call get_axissym_uvmem(nch_cone,k,rzcoeffs_cone,rzvals_cone, &
         rlam,ppw,norder,npatches_cone)

      call get_axissym_uvmem(nch_lens,k,rzcoeffs_lens,rzvals_lens, &
         rlam,ppw,norder,npatches_lens)
      

      print *, npatches_rhab,npatches_cone,npatches_lens

      npols = (norder+1)*(norder+2)/2
      npts_rhab = npatches_rhab*npols
      npts_cone = npatches_cone*npols
      npts_lens = npatches_lens*npols

      return
      end

      




      
      subroutine get_axissym_miniwasp_geom(rlam,ppw,norder, &
        npatches_rhab, &
        norders_rhab, ixyzs_rhab, iptype_rhab, npts_rhab, &
        srcvals_rhab, srccoefs_rhab, npatches_cone, &
        norders_cone, ixyzys_cone, iptype_cone, npts_cone, &
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

      fname1 = '../geometries/rhabdom_axissym_pts.dat'
      fname2 = '../geometries/cone_axissym_pts.dat'
      fname3 = '../geometries/lens_axissym_pts.dat'

      open(unit=33,file=trim(fname1))
      do i=1,nv_rhab
        read(33,*) verts_rhab(1,i),verts_rhab(2,i)
      enddo
      close(33)

      open(unit=33,file=trim(fname2))
      do i=1,nv_cone
        read(33,*) verts_cone(1,i),verts_cone(2,i)
      enddo
      close(33)

      open(unit=33,file=trim(fname3))
      do i=1,nv_lens
        read(33,*) verts_lens(1,i),verts_lens(2,i)
      enddo
      close(33)
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
         rlam,ppw,norder,npatches_rhab,norders_rhab,ixyzs_rhab, &
         iptype_rhab,npts_rhab,srcvals_rhab,srccoefs_rhab)

      call get_axissym_geom(nch_cone,k,rzcoeffs_cone,rzvals_cone, &
         rlam,ppw,norder,npatches_cone,norders_cone,ixyzs_cone, &
         iptype_cone,npts_cone,srcvals_cone,srccoefs_cone)

      call get_axissym_geom(nch_lens,k,rzcoeffs_lens,rzvals_lens, &
         rlam,ppw,norder,npatches_lens,norders_lens,ixyzs_lens, &
         iptype_lens,npts_lens,srcvals_lens,srccoefs_lens)
      print *, npatches_rhab,npatches_cone,npatches_lens
      print *, npts_rhab,npts_cone,npts_lens
      print *, npts_rhab+npts_cone+npts_lens

      if(ifplot.eq.1) then

       call surf_vtk_plot(npatches_rhab,norders_rhab,ixyzs_rhab, &
         iptype_rhab,npts_rhab,srccoefs_rhab,srcvals_rhab, &
         'rhabdom_axissym.vtk','a')

        call surf_vtk_plot(npatches_cone,norders_cone,ixyzs_cone, &
         iptype_cone,npts_cone,srccoefs_cone, srcvals_cone, &
         'cone_axissym.vtk','a')

        call surf_vtk_plot(npatches_lens,norders_lens,ixyzs_lens, &
         iptype_lens,npts_lens,srccoefs_lens, srcvals_lens, &
         'lens_axissym.vtk','a')
       endif

      return
      end

      




      subroutine get_axissym_uvmem(nch,k,rzcoefs,rzvals,rlam,ppw, &
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

        ns = ceiling(rlen*(ppw+0.0d0)/rlam/(norder+0.0d0))
        nt = ceiling(rmax*2*pi*(ppw+0.0d0)/rlam/(norder+0.0d0))
        npatches = npatches + 2*ns*nt
      enddo


      return
      end





      subroutine get_axissym_geom(nch,k,rzcoefs,rzvals,rlam,ppw, &
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
      
      vmin = 0
      vmax = 2*pi

      allocate(triaskel(3,3,npatches))
      allocate(ichuse(npatches))


      npatches0 = 0
      nover = 0
      istart = 1
      do ich=1,nch
!
!  find max r on each patch, and length of each patch
!
        rmax = 0.0d0
        rlen = 0.0d0
        rs(1:2) = 0

        if(ich.eq.56) then
          call prin2('rzvals=*',rzvals(1,1:k,ich),k)
          call prin2('rzvals=*',rzvals(3,1:k,ich),k)
          call prin2('rzvals=*',rzvals(4,1:k,ich),k)
        endif
        do j=1,k
          rs(1:2) = rs(1:2) + amatrint(1:2,j)*rzvals(1,j,ich)
          if(rzvals(1,j,ich).ge.rmax) rmax = rzvals(1,j,ich)
          dsdt = sqrt(rzvals(3,j,ich)**2 + rzvals(4,j,ich)**2)
          rlen = rlen + dsdt*ws(j)
        enddo
        if(rs(1).ge.rmax) rmax = rs(1)
        if(rs(2).ge.rmax) rmax = rs(2)

        ns = ceiling(rlen*(ppw+0.0d0)/rlam/(norder+0.0d0))
        nt = ceiling(rmax*2*pi*(ppw+0.0d0)/rlam/(norder+0.0d0))
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ns,nt,nover, &
          npatches,npatches0,triaskel(1,1,istart))
        ichuse(istart:istart+2*ns*nt-1) = ich 
        istart = istart + 2*ns*nt
      enddo

      kuse = k
      print *, "kuse=",kuse
      print *, "norder=",norder

      open(unit=35,file='tmp1.dat')
      do i=1,npatches
        write(35,*) i,ichuse(i)

      enddo

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










