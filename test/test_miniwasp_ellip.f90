      implicit real *8 (a-h,o-z)
      real *8 abc_lens(3),abc_cone(3),abc_rhabdom(3),abc_pig(3)
      real *8 xyz_lens(3),xyz_cone(3),xyz_rhabdom(3),xyz_pig(3)
      integer npatches_vect(4),npts_vect(4) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      real *8 errs(6),ts(2),xyz0(3)
      real *8, allocatable :: rfacs(:,:)
      character *100 fname
      integer ipars(2)
      integer, allocatable :: row_ptr(:),col_ind(:)
      integer, allocatable :: iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      complex *16, allocatable :: uval(:),dudnval(:)
      complex *16, allocatable :: sigmaover(:),slp_near(:),dlp_near(:)
      complex *16, allocatable :: pot(:),potslp(:),potdlp(:)
      complex *16, allocatable :: potslp2(:)

      complex *16 zk

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ixyzso(:),nfars(:)

      integer, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:)
      complex * 16 zpars(3)

      call prini(6,13)

      iref = 1
      npatches = 0
      call get_miniwasp_ellip_params(abc_lens,abc_cone, &
        abc_rhabdom,abc_pig,xyz_lens,xyz_cone,xyz_rhabdom,xyz_pig)
      call get_miniwasp_ellip_mem(abc_lens,abc_cone, &
        abc_rhabdom,abc_pig,iref,npatches)
      call prinf('npatches=*',npatches,1)

      norder = 5
      npols = (norder+1)*(norder+2)/2
      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      call get_miniwasp_ellip(abc_lens,abc_cone,abc_rhabdom, &
        abc_pig,xyz_lens,xyz_cone,xyz_rhabdom,xyz_pig,iref,npatches, &
        npatches_vect,norder,srcvals,srccoefs,npts,npts_vect)
      stop


      allocate(iptype(npatches),ixyzs(npatches+1),norders(npatches))

      do i=1,npatches
        ixyzs(i) = (i-1)*npols+1
        norders(i) = norder
        iptype(i) = 1
      enddo
      ixyzs(npatches+1) = 1+npts

      xyz_in(1) = 0.01d0 + xyz0(1)
      xyz_in(2) = 0.32d0 + xyz0(2)
      xyz_in(3) = -0.13d0 + xyz0(3)

      
      xyz_out(1) = 7.2d0 + xyz0(1)
      xyz_out(2) = -6.3d0 + xyz0(2)
      xyz_out(3) = 5.4d0 + xyz0(3)

      zk  = 1.1d0
      allocate(targs(3,npts))
      allocate(ixyzso(npatches+1),nfars(npatches))
      
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)



      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      allocate(pot(npts),potslp(npts),potdlp(npts))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
          srccoefs,cms,rads)

      allocate(sigma(npts),uval(npts),dudnval(npts))

      do i=1,npts
        call h3d_slp(xyz_out,3,srcvals(1,i),0,dpars,1,zk,0, &
           ipars,uval(i))
        call h3d_sprime(xyz_out,12,srcvals(1,i),0,dpars,1,zk,0, &
            ipars,dudnval(i))
      enddo

      ndtarg = 3
     
      do i=1,npts
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo

      allocate(ipatch_id(npts),uvs_targ(2,npts))
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, & 
              ipatch_id,uvs_targ)

 
!
!    find near field
!
      iptype = 1
      call get_rfacs(norder,iptype,rfac,rfac0)
      do i=1,npatches 
        rad_near(i) = rads(i)*rfac
      enddo
      

      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr, &
             col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind, &
             iquad)

      nquad = iquad(nnz+1)-1
      allocate(slp_near(nquad),dlp_near(nquad))


      ndtarg = 3

      eps = 0.50001d-9

      ikerorder = -1



      call cpu_time(t1)
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms, &
         rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zk, &
         nnz,row_ptr,col_ind,rfac,nfars,ixyzso)
      call cpu_time(t2)
      tfar = t2-t1


      npts_over = ixyzso(npatches+1)-1

      print *, "npts_over=",npts_over


      allocate(srcover(12,npts_over),sigmaover(npts_over), &
              wover(npts_over))

          
      call oversample_geom(npatches,norders,ixyzs,iptype,npts, & 
         srccoefs,srcvals,nfars,ixyzso,npts_over,srcover)

      call get_qwts(npatches,nfars,ixyzso,iptype,npts_over, &
             srcover,wover)


      do i=1,nquad
        slp_near(i) = 0
        dlp_near(i) = 0
      enddo



      call cpu_time(t1)

      zpars(1) = zk
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0

      iquadtype = 1


      call getnearquad_helm_comb_dir(npatches,norders, &
           ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs, &
           ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind, &
           iquad,rfac0,nquad,slp_near)

      
      zpars(2) = 0.0d0
      zpars(3) = 1.0d0
      call getnearquad_helm_comb_dir(npatches,norders, &
           ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs, &
           ipatch_id,uvs_targ,eps,zpars,iquadtype, &
           nnz,row_ptr,col_ind,iquad,rfac0,nquad,dlp_near)
      
      call cpu_time(t2)
      tquadgen = t2-t1



      ifinout = 1     

      zpars(2) = 1.0d0
      zpars(3) = 0.0d0


      call cpu_time(t1)

      call lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs, &
       iptype,npts,srccoefs,srcvals,ndtarg,npts,targs, &
       eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,slp_near, &
       dudnval,nfars,npts_over,ixyzso,srcover,wover,potslp)


      zpars(2) = 0.0d0
      zpars(3) = 1.0d0


      call lpcomp_helm_comb_dir_setsub(npatches,norders,ixyzs, &
       iptype,npts,srccoefs,srcvals,ndtarg,npts,targs, &
       eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,dlp_near, &
       uval,nfars,npts_over,ixyzso,srcover,wover,potdlp)


      call cpu_time(t2)
      tlpcomp = t2-t1


!
!
!      compute error
!
      errl2 = 0
      rl2 = 0
      do i=1,npts
        pot(i) = (potslp(i) - potdlp(i))*2
        errl2 = errl2 + abs(uval(i)-pot(i))**2*wts(i)
        rl2 = rl2 + abs(uval(i))**2*wts(i)
      enddo


      err = sqrt(errl2/rl2)

      call prin2('error in greens identity=*',err,1)


      stop
      end
!
!
!
!
!
      subroutine get_miniwasp_ellip(abc_lens,abc_cone,abc_rhabdom, &
        abc_pig,xyz_lens,xyz_cone,xyz_rhabdom,xyz_pig,iref,npatches, &
        npatches_vect,norder,srcvals,srccoefs,npts,npts_vect)
      implicit real *8 (a-h,o-z)
      real *8 abc_lens(3),abc_cone(3),abc_rhabdom(3),abc_pig(3)
      real *8 xyz_lens(3),xyz_cone(3),xyz_rhabdom(3),xyz_pig(3)
      real *8 srcvals(12,npts),srccoefs(9,npts)
      integer npatches_vect(4),npts_vect(4)
      character *100 fname
      
      istart = 1
      npols = (norder+1)*(norder+2)/2

      a = abc_lens(1)
      b = abc_lens(2)
      c = abc_lens(3)
      fname = 'lens.vtk'
      call get_ellipsoid_geom(a,b,c,xyz_lens,norder,npatches_vect(1), &
       iref,srcvals(1,istart),srccoefs(1,istart),ifplot,trim(fname))
      npts_vect(1) = npatches_vect(1)*npols

      istart = istart + npts_vect(1)

      a = abc_cone(1)
      b = abc_cone(2)
      c = abc_cone(3)
      fname = 'cone.vtk'
      call get_ellipsoid_geom(a,b,c,xyz_cone,norder,npatches_vect(2), &
       iref,srcvals(1,istart),srccoefs(1,istart),ifplot,trim(fname))
      npts_vect(2) = npatches_vect(2)*npols

      istart = istart + npts_vect(2)

      a = abc_rhabdom(1)
      b = abc_rhabdom(2)
      c = abc_rhabdom(3)
      fname = 'lens.vtk'
      call get_ellipsoid_geom(a,b,c,xyz_rhabdom,norder,npatches_vect(3), &
       iref,srcvals(1,istart),srccoefs(1,istart),ifplot,trim(fname))
      npts_vect(3) = npatches_vect(3)*npols

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

      subroutine get_miniwasp_ellip_mem(abc_lens,abc_cone, &
        abc_rhabdom,abc_pig,iref,npatches)
      implicit real *8 (a-h,o-z)
      real *8 abc_lens(3),abc_cone(3),abc_rhabdom(3),abc_pig(3)
      
      npatches = 0
      call get_rectparapiped_mem(abc_lens(1),abc_lens(2), &
        abc_lens(3),iref,npatches0)
      npatches = npatches + npatches0

      call get_rectparapiped_mem(abc_cone(1),abc_cone(2), &
        abc_cone(3),iref,npatches0)
      npatches = npatches + npatches0

      call get_rectparapiped_mem(abc_rhabdom(1),abc_rhabdom(2), &
        abc_rhabdom(3),iref,npatches0)
      npatches = npatches + npatches0

      call get_rectparapiped_mem(abc_pig(1),abc_pig(2), &
        abc_pig(3),iref,npatches0)
      npatches = npatches + npatches0

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

      na = nint(a/rmin)*2**iref
      nb = nint(b/rmin)*2**iref
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
      na = nint(a/rmin)
      nb = nint(b/rmin)
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








