      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 dpars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: rhs(:)
      complex *16, allocatable :: psinm(:,:),phinm(:,:),dfuv(:,:)
      complex *16, allocatable :: vynm(:,:)

      complex *16, allocatable :: vynm_targ(:,:)
      complex *16, allocatable :: phinm_targ(:,:)
      complex *16, allocatable :: psinm_targ(:,:)
      
      complex *16, allocatable :: wnear(:,:)
      real *8, allocatable :: targs(:,:)

      complex *16 contrast_matrix(4,1)
      real *8, allocatable :: srcvals_extended(:,:)
      complex *16 omega,ep,mu,ep1,mu1,ep0,mu0
      complex *16, allocatable :: rhs_muller(:)
      complex *16, allocatable :: soln_muller(:)
      real *8, allocatable :: dxyzdu(:,:),dxyzdv(:,:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:),errs(:)
      real *8, allocatable :: srcover(:,:),wover(:)

      integer, allocatable :: col_ptr(:),row_ind(:)
      integer, allocatable :: ixyzso(:),novers(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)

      complex *16, allocatable :: ecomp(:,:),eex(:,:),hcomp(:,:)
      complex *16, allocatable :: h_ex(:,:)
      complex *16, allocatable :: pot1(:),pot2(:),pot3(:)

      complex *16 zalpha,zbeta,zgamma,zdelta,zeta,zteta,zk,ztetap
      complex *16 zk0
      complex *16 ztetam
      complex *16 fjvals(0:100),fhvals(0:100),fjder(0:100),fhder(0:100)
      complex *16 fjvals0(0:100),fhvals0(0:100)
      complex *16 fjder0(0:100),fhder0(0:100)
      complex *16 fjvalst(0:100),fhvalst(0:100)
      complex *16 fjdert(0:100),fhdert(0:100)
      complex *16 z1,z2,z3,z4
      complex *16 zvec1(3),zvec2(3),zvec3(3)
      real *8 dvec1(3),dvec2(3),dvec3(3)

      integer sorted_vector(100)
      logical exposed_surfaces(100)


      real *8 thet,phi,eps_gmres
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,dirname
      character *300 fname

      real *8, allocatable :: w(:,:)

      logical isout0,isout1

      complex *16 ztmp,ima
      procedure (), pointer :: fker
      external h3d_sgradx, h3d_sgrady, h3d_sgradz

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


      igeomtype = 1
      ipars(1) = 3 
      npatches=12*(4**ipars(1))

      norder = 3 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 0

      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      print *, 'npts=',npts

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


c
c       define rhs to be one of the ynm's
c
      nn = 3
      mm = 1
      nmax = nn
      allocate(w(0:nmax,0:nmax))
      allocate(rhs(npts))
      call l3getsph(nmax,mm,nn,12,srcvals,rhs,npts,w)
c
c
c   set material parameters
c
      contrast_matrix(1,1)=1.1d0  
      contrast_matrix(2,1)=1.1d0  
      contrast_matrix(3,1)=1.2d0  
      contrast_matrix(4,1)=1.0d0

      omega = 1.0d0
      ep0 = contrast_matrix(1,1)
      mu0 = contrast_matrix(2,1)
      ep = contrast_matrix(3,1)
      mu = contrast_matrix(4,1)
c
c  interior wave number
c
      zk = omega*sqrt(ep)*sqrt(mu)
      zk0 = omega*sqrt(ep0)*sqrt(mu0)

      call prin2('zk=*',zk,2)
      call prin2('zk0=*',zk0,2)


      njh = 30
      ifder = 1
      rscale = 1.0d0
      call prin2('zk=*',zk,2)
      call besseljs3d(njh,zk,rscale,fjvals,ifder,fjder)
      call h3dall(njh,zk,rscale,fhvals,ifder,fhder)

      call besseljs3d(njh,zk0,rscale,fjvals0,ifder,fjder0)
      call h3dall(njh,zk0,rscale,fhvals0,ifder,fhder0)
      call prin2('fjvals=*',fjvals,2*njh)
      call prin2('fjvals=*',fjvals0,2*njh)



      allocate(dfuv(2,npts))

      call get_surf_grad(2,npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,rhs,dfuv)


      allocate(psinm(3,npts),phinm(3,npts),vynm(3,npts))
      do i=1,npts
        psinm(1:3,i) = dfuv(1,i)*srcvals(4:6,i) + 
     1     dfuv(2,i)*srcvals(7:9,i) 
        call dzcross_prod3d(srcvals(10,i),psinm(1,i),phinm(1,i))
      enddo
cc      call l3getsph_vec(mm,nn,12,npts,srcvals,vynm,psinm,
cc     1   phinm)

      call prin2('psinm=*',psinm,24)
      call prin2('phinm=*',phinm,24)

      allocate(dxyzdu(3,npts),dxyzdv(3,npts))
      do i=1,npts
        dvec1(1) = srcvals(4,i)
        dvec1(2) = srcvals(5,i)
        dvec1(3) = srcvals(6,i)

        dvec2(1) = srcvals(10,i)
        dvec2(2) = srcvals(11,i)
        dvec2(3) = srcvals(12,i)

        dxyzdu(1:3,i) = 0
        dxyzdv(1:3,i) = 0

        call orthonormalize(dvec1,dvec2,dxyzdu(1,i),dxyzdv(1,i))

      enddo
      call prin2('dxyzdu=*',dxyzdu,6)
      call prin2('dxyzdv=*',dxyzdv,6)

c
c  get boundary data
c
c
      allocate(rhs_muller(4*npts),soln_muller(4*npts))
      call prin2('omega=*',omega,2)
      call prin2('ep=*',ep,2)
      call prin2('ep0=*',ep0,2)
      call prin2('mu=*',mu,2)
      call prin2('mu0=*',mu0,2)
      call prin2('zk=*',zk,2)
      call prin2('zk0=*',zk0,2)
      call prinf('nn=*',nn,1)
      z1 = ima*mu0*(fjvals0(nn) + zk0*fjder0(nn))*zk0*fhvals0(nn) 
      z1 = z1-ima*mu*fjvals(nn)*zk*(fhvals(nn)+fhder(nn)*zk)
      z1 = z1/(mu+mu0)

      z2=zk0*(fjvals0(nn)+fjder0(nn)*zk0)*(fhvals0(nn)+fhder0(nn)*zk0)
      z2=z2-zk*(fjvals(nn)+fjder(nn)*zk)*(fhvals(nn)+fhder(nn)*zk)
      z2 = z2/omega/(ep+ep0)
      print *, "z1=",z1
      print *, "z2=",z2
      do i=1,npts

        rhs_muller(i) = z1*(dxyzdu(1,i)*psinm(1,i) + 
     1      dxyzdu(2,i)*psinm(2,i) + dxyzdu(3,i)*psinm(3,i))


        rhs_muller(i+npts) = z1*(dxyzdv(1,i)*psinm(1,i) + 
     1      dxyzdv(2,i)*psinm(2,i) + dxyzdv(3,i)*psinm(3,i))

        rhs_muller(i+2*npts) = z2*(dxyzdu(1,i)*phinm(1,i) + 
     1      dxyzdu(2,i)*phinm(2,i) + dxyzdu(3,i)*phinm(3,i))
        rhs_muller(i+3*npts) = z2*(dxyzdv(1,i)*phinm(1,i) + 
     1      dxyzdv(2,i)*phinm(2,i) + dxyzdv(3,i)*phinm(3,i))
      enddo

      ra = 0
      do i=1,npts
        ra = ra + wts(i)
      enddo
      call prin2('error in surfce area=*',ra-4*pi,1)
c
c  find topological sorting
c
c
      eps = 1.0d-7
      n_components = 1
      call topological_sorting(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,wts,npatches,n_components,sorted_vector,
     2  exposed_surfaces,eps)
      print *, exposed_surfaces(1)
      call prinf('sorted_vector=*',sorted_vector,2)

      allocate(srcvals_extended(20,npts))

      call build_extended_targ(n_components,srcvals_extended,
     1  srcvals,npts,contrast_matrix,npts)
      call prin2('srcvals_extended=*',srcvals_extended,20)
      call prin2('srcvals=*',srcvals,20)

c
c  test self near quadrature to ensure correct boundary
c  data is generated
c
c
      allocate(ipatch_id(npts),uvs_targ(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)

      iptype_avg = 1
      norder_avg = norder
      call get_rfacs(norder_avg,iptype,rfac,rfac0)
      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,cms,rads)
      
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
c        rad_near(i) = 3.0d0
      enddo

      call findnearmem(cms,npatches,rad_near,12,srcvals,npts,nnz)
      allocate(row_ptr(npts+1),col_ind(nnz))
      call findnear(cms,npatches,rad_near,12,srcvals,npts,row_ptr,
     1  col_ind)
      allocate(novers(npatches),ixyzso(npatches+1))
      ikerorder = 0
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1 rads,npts,srccoefs,12,npts,srcvals,ikerorder,omega,
     2 nnz,row_ptr,col_ind,rfac,novers,ixyzso)
cc      do i=1,npatches
cc        novers(i) = norder+2
cc        npols0 = (novers(i)+1)*(novers(i)+2)/2
cc        ixyzso(i) = (i-1)*npols0+1
cc      enddo
cc      ixyzso(npatches+1) = npatches*npols0+1 

      npts_over = ixyzso(npatches+1)-1
      allocate(srcover(12,npts_over),wover(npts_over))
      call oversample_geom(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,novers,ixyzso,npts_over,srcover)
      
      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1 srcover,wover)

      allocate(iquad(nnz+1))
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind, 
     1  iquad)
      nquad = iquad(nnz+1)-1
      
      allocate(wnear(nquad,16))
      iquadtype = 1
      wnear = 0
      goto 1111
      call getnearquad_em_muller_trans_v2(npatches,norders,
     1  ixyzs,iptype,npts,srccoefs,srcvals,20,npts,srcvals_extended,
     2  ipatch_id,uvs_targ,eps,omega,iquadtype,nnz,row_ptr,col_ind,
     3  iquad,rfac0,nquad,wnear)
 1111 continue     
      call prinf('nnz=*',nnz,1)
      call prin2('eps=*',eps,1)
      call prinf('row_ptr=*',row_ptr,20)
      call prin2('wnear=*',wnear(1,4),24)

c
c   set boundary data
c 
      allocate(pot1(4*npts),pot2(4*npts))
      do i=1,4*npts
        pot1(i) = 0
        pot2(i) = 0
      enddo

      do i=1,npts
        pot1(i) = psinm(1,i)*dxyzdu(1,i) + psinm(2,i)*dxyzdu(2,i) + 
     1      psinm(3,i)*dxyzdu(3,i)
        pot1(i+npts) = psinm(1,i)*dxyzdv(1,i)+psinm(2,i)*dxyzdv(2,i) + 
     1      psinm(3,i)*dxyzdv(3,i)
      enddo
      call prin2('pot1=*',pot1,24)


      call lpcomp_em_muller_trans_v2_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,12,npts,srcvals,eps,omega,
     2   nnz,row_ptr,col_ind,iquad,nquad,pot1,novers,npts_over,
     3   ixyzso,srcover,wover,pot2,wnear,n_components,
     4   contrast_matrix,npts)
      

      call prin2('pot2=*',pot2(2*npts+1),24)
c
c
c  test accuracy of pot2, first test eletric field
c
      erra = 0
      ra = 0
      do j=1,4
        istart = (j-1)*npts
        do ii=1,npts
          i = istart+ii
          pot2(i) = pot2(i) + 0.5d0*pot1(i)
          erra = erra + abs(rhs_muller(i)-pot2(i))**2*wts(ii)
          ra = ra + abs(rhs_muller(i))**2*wts(ii)
        enddo
      enddo

      call prin2('pot2=*',pot2,24)
      call prin2('rhs_muller=*',rhs_muller,24)
      erra = sqrt(erra/ra)
      call prin2('erra in e dot ru=*',erra,1)




      numit = 200
      allocate(errs(numit+1))
      eps_gmres = 1.0d-10

      do i=1,4*npts
        soln_muller(i)=0
      enddo
      goto 1211
      call em_muller_trans_v2_solver(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,eps,omega,numit,ifinout,rhs_muller,
     2  eps_gmres,niter,errs,rres,soln_muller,contrast_matrix,
     3  npts,n_components,srcvals_extended)
      
      call prin2('soln_muller=*',soln_muller,24)
      call prin2('soln_muller2=*',soln_muller(npts+1),24)
      call prin2('soln_muller3=*',soln_muller(2*npts+1),24)
      call prin2('soln_muller4=*',soln_muller(3*npts+1),24)
      erra = 0
      ra = 0
      do i=1,npts
        zvec1(1:3) = (soln_muller(i)*dxyzdu(1:3,i) + 
     1     soln_muller(i+npts)*dxyzdv(1:3,i))
        erra = erra + abs(psinm(1,i) - zvec1(1))**2*wts(i)
        erra = erra + abs(psinm(2,i) - zvec1(2))**2*wts(i)
        erra = erra + abs(psinm(3,i) - zvec1(3))**2*wts(i)
        ra = ra + abs(psinm(1,i))**2*wts(i)

        if(i.lt.5) then
          call prin2('zvec=*',zvec1,6)
          call prin2('psinm=*',psinm(1,i),6)
          print *, psinm(1,i)/zvec1(1)
        endif
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in soln=*',erra,1)

c
c  now test target evaluator
c
c
 1211 continue
      
      if(1.eq.1) then   
        do i=1,4*npts
         soln_muller(i) = 0
        enddo

        do i=1,npts
          soln_muller(i) = psinm(1,i)*dxyzdu(1,i) + 
     1                     psinm(2,i)*dxyzdu(2,i) + 
     1                     psinm(3,i)*dxyzdu(3,i)  

          soln_muller(i+npts) = psinm(1,i)*dxyzdv(1,i) + 
     1                     psinm(2,i)*dxyzdv(2,i) + 
     1                     psinm(3,i)*dxyzdv(3,i)  

        enddo
      endif
      ntarg = 20
      allocate(targs(3,ntarg))
      do i=1,ntarg
        r = hkrand(0)*0.8d0
        thet = hkrand(0)*pi
        phi = hkrand(0)*2*pi

        r = 0.9999d0

        targs(1,i) = r*sin(thet)*cos(phi)
        targs(2,i) = r*sin(thet)*sin(phi)
        targs(3,i) = r*cos(thet)

      enddo

      allocate(eex(3,ntarg),ecomp(3,ntarg),hcomp(3,ntarg))
      allocate(h_ex(3,ntarg))
      call prin2('soln_muller=*',soln_muller,24)

      call evaluate_field_muller(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,wts,targs,ntarg,npatches,n_components,
     2  sorted_vector,contrast_matrix,exposed_surfaces,eps,omega,
     3  soln_muller,ecomp,hcomp)

      allocate(vynm_targ(3,ntarg),phinm_targ(3,ntarg))
      allocate(psinm_targ(3,ntarg))


      call l3getsph_vec(mm,nn,3,ntarg,targs,vynm_targ,psinm_targ,
     1   phinm_targ)

      erra = 0
      ra = 0

      erra2 = 0
      ra2 = 0
      do i=1,ntarg
        r = sqrt(targs(1,i)**2 + targs(2,i)**2 + targs(3,i)**2)
        z1 = zk*r
        call besseljs3d(njh,z1,rscale,fjvalst,ifder,fjdert)
        z2 = -ima*fjvalst(nn)*zk*(zk*fhder(nn)+fhvals(nn))*mu
        eex(1:3,i) = phinm_targ(1:3,i)*z2

        z1 = zk*r
        call besseljs3d(njh,z1,rscale,fjvalst,ifder,fjdert)
        z2 = ima*nn*(nn+1.0d0)*zk*fjvalst(nn)/r*(fhvals(nn) +
     1       zk*fhder(nn))/ima/omega
        z3 = ima*(zk**2*fjdert(nn) + fjvalst(nn)*zk/r)*(fhvals(nn) + 
     1       zk*fhder(nn))/ima/omega
        h_ex(1:3,i) = psinm_targ(1:3,i)*z3 + vynm_targ(1:3,i)*z2
        erra = erra + abs(h_ex(1,i)-hcomp(1,i))**2
        erra = erra + abs(h_ex(2,i)-hcomp(2,i))**2
        erra = erra + abs(h_ex(3,i)-hcomp(3,i))**2
        ra = ra + abs(h_ex(1,i))**2
        ra = ra + abs(h_ex(2,i))**2
        ra = ra + abs(h_ex(3,i))**2

        erra2 = erra2 + abs(eex(1,i)-ecomp(1,i))**2
        erra2 = erra2 + abs(eex(2,i)-ecomp(2,i))**2
        erra2 = erra2 + abs(eex(3,i)-ecomp(3,i))**2
        ra2 = ra2 + abs(eex(1,i))**2
        ra2 = ra2 + abs(eex(2,i))**2
        ra2 = ra2 + abs(eex(3,i))**2
      enddo
      call prin2('ecomp=*',ecomp,6)
      call prin2('eex=*',eex,6)

      erra = sqrt(erra/ra)
      erra2 = sqrt(erra2/ra2)
      call prin2('erra=*',erra,1)
      call prin2('ra=*',sqrt(ra),1)
      call prin2('error in H=*',erra,1)

      call prin2('error in E=*',erra2,1)

c
      


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
        vmin = 0
        vmax = 2*pi
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

   



      subroutine l3getsph(nmax,mm,nn,ndx,xyzs,ynms,npts,ynm)
      implicit real *8 (a-h,o-z)
      real *8 :: xyzs(ndx,npts)
      complex *16 ynms(npts),ima
      real *8 rat1(10000),rat2(10000)
      real *8 ynm(0:nmax,0:nmax)
      data ima/(0.0d0,1.0d0)/
  
      call ylgndrini(nmax, rat1, rat2)
  
      do i=1,npts
        x=xyzs(1,i)
        y=xyzs(2,i)
        z=xyzs(3,i)
        r=sqrt(x**2+y**2+z**2)
        call cart2polar(xyzs(1,i),r,theta,phi)
        ctheta = cos(theta)
        call ylgndrf(nmax, ctheta, ynm, rat1, rat2)
        ynms(i) = ynm(nn,abs(mm))*exp(ima*mm*phi)        
      enddo
       
      return
      end
c
c
c
c
c
      subroutine l3getsph_vec(mm,nn,ndx,npts,xyzs,vynm,psinm,phinm)
      implicit real *8 (a-h,o-z)
      real *8 :: xyzs(ndx,npts)
      complex *16 ima
      complex *16 vynm(3,npts),psinm(3,npts),phinm(3,npts)
      real *8 vtmp(3)
      real *8, allocatable :: wlege(:),ynm(:,:),ynmd(:,:)
      complex *16 zr,zt,zp
      data ima/(0.0d0,1.0d0)/

      nmax = nn+1
  
      nlege = nmax + 10
      lw = (nlege+1)**2*4
      allocate(wlege(lw),ynm(0:nmax,0:nmax),ynmd(0:nmax,0:nmax))
      call ylgndrfwini(nlege,wlege,lw,lused)
      
  
      do i=1,npts
        x=xyzs(1,i)
        y=xyzs(2,i)
        z=xyzs(3,i)
        call cart2polar(xyzs(1,i),r,thet,phi)
        ctheta = cos(thet)
        rx = sin(thet)*cos(phi)
        ry = sin(thet)*sin(phi)
        rz = cos(thet)

        thetx = cos(thet)*cos(phi)
        thety = cos(thet)*sin(phi)
        thetz = -sin(thet)

        phix = -sin(phi)
        phiy = cos(phi)
        phiz = 0

        call ylgndr2sfw(nmax,ctheta,ynm,ynmd,wlege,nlege)

        vtmp(1) = x/r
        vtmp(2) = y/r
        vtmp(3) = z/r

        if(mm.eq.0) then
          vynm(1,i) = ynm(nn,0)*rx 
          vynm(2,i) = ynm(nn,0)*ry 
          vynm(3,i) = ynm(nn,0)*rz

          psinm(1,i) = -sin(thet)*ynmd(nn,0)*thetx
          psinm(2,i) = -sin(thet)*ynmd(nn,0)*thety
          psinm(3,i) = -sin(thet)*ynmd(nn,0)*thetz
        else
          zr = ynm(nn,abs(mm))*sin(thet)*exp(ima*mm*phi)
          vynm(1,i) = zr*rx
          vynm(2,i) = zr*ry 
          vynm(3,i) = zr*rz

          zt = -ynmd(nn,abs(mm))*exp(ima*mm*phi)
          zp = ima*mm*ynm(nn,abs(mm))*exp(ima*mm*phi)

          psinm(1,i) = zt*thetx + zp*phix
          psinm(2,i) = zt*thety + zp*phiy
          psinm(3,i) = zt*thetz + zp*phiz
        endif
        call dzcross_prod3d(vtmp,psinm(1,i),phinm(1,i))
      enddo
       
      return
      end



