      subroutine getnearquad_em_muller_trans_v2(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the integral DFIE:
!
!  Note: the 4 \pi scaling is NOT!! included as the output of the FMM
!  has been rescaled.
!
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a chunk centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!
!        
!
!  input:
!    npatches - integer
!      number of patches
!
!    norders - integer(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(npatches+1)
!      starting location of data on patch i
!  
!    iptype - integer(npatches)
!      type of patch
!      iptype = 1 -> triangular patch discretized with RV nodes
!
!    npts - integer
!      total number of discretization points on the boundary
!
!    srccoefs - real *8 (9,npts)
!      koornwinder expansion coefficients of xyz, dxyz/du,
!      and dxyz/dv on each patch. 
!      For each point srccoefs(1:3,i) is xyz info
!        srccoefs(4:6,i) is dxyz/du info
!        srccoefs(7:9,i) is dxyz/dv info
!
!    srcvals - real *8 (12,npts)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
! 
!    ndtarg - integer
!      leading dimension of target array
!        
!    ntarg - integer
!      number of targets
!
!    targs - real *8 (ndtarg,ntarg)
!      target information
!
!    ipatch_id - integer(ntarg)
!      id of patch of target i, id = -1, if target is off-surface
!
!         uvs_targ - real *8 (2,ntarg)
!            local uv coordinates on patch if on surface, otherwise
!            set to 0 by default
!          (maybe better to find closest uv on patch using
!            newton)
!            
!          eps - real *8
!             precision requested
!
!          zpars - complex *16(3)
!              kernel parameters
!              zpars(1) = omega 
!              zpars(2) = ep0
!              zpars(3) = mu0
!              zpars(4) = ep1
!              zpars(5) = mu1
!
!           iquadtype - integer
!              quadrature type
!              iquadtype = 1, use ggq for self + adaptive integration
!                 for rest
! 
!
!           nnz - integer
!             number of source patch-> target interactions in the near
!             field
! 
!           row_ptr - integer(ntarg+1)
!              row_ptr(i) is the pointer
!              to col_ind array where list of relevant source patches
!              for target i start
!
!           col_ind - integer (nnz)
!               list of source patches relevant for all targets, sorted
!               by the target number
!
!           iquad - integer(nnz+1)
!               location in wnear array where quadrature for col_ind(i)
!               starts
!
!           rfac0 - integer
!               radius parameter for near field
!
!           nquad - integer
!               number of entries in wnear
!
!        output
!           wnear - complex *16(36*nquad)
!               the desired near field quadrature
!

      implicit none 
      integer npatches,norders(npatches),npts,nquad
      integer ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      real *8 srcvals_extended(20,npts)
      integer ndtarg,ntarg
      integer iquadtype
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(5)
      integer nnz,ipars(2)
      real *8 dpars(1)
      integer row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(16*nquad)
	  
      integer ipatch_id(ntarg)
      real *8 uvs_targ(2,ntarg)

      complex *16 alpha,beta
      integer i,j,ndi,ndd,ndz,count1,count2,icount

      integer ipv
      integer :: t1, t2,clock_rate, clock_max

      procedure (), pointer :: fker
      external  fker_em_muller_trans_v2
      external  em_muller_trans_v2

     ndz=5
	   ndd=1
	   ndi=2
	   ipv=1

!!      fker => fker_em_muller_trans_v2
	  fker => em_muller_trans_v2
!    call system_clock ( t1, clock_rate, clock_max )

!	  icount=0
!	  do count1=1,4
!	    do count2=1,4
!		  ipars(1)=count1
!      ipars(2)=count2
!      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
!		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
!		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
!		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
!		   &wnear(icount*nquad+1:(icount+1)*nquad))
!		   icount=icount+1
!      enddo
!	  enddo

		  ipars(1)=1
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(1:nquad))
       print *, "done with kernel 1"

		  ipars(1)=1
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(nquad+1:2*nquad))
       print *, "done with kernel 2"

		  ipars(1)=1
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(2*nquad+1:3*nquad))
       print *, "done with kernel 3"

		  ipars(1)=1
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(3*nquad+1:4*nquad))
       print *, "done with kernel 4"


		  ipars(1)=2
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(4*nquad+1:5*nquad))
       print *, "done with kernel 5"

		  ipars(1)=2
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(5*nquad+1:6*nquad))

		  ipars(1)=2
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(6*nquad+1:7*nquad))
       print *, "done with kernel 6"

		  ipars(1)=2
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(7*nquad+1:8*nquad))
       print *, "done with kernel 7"

		  ipars(1)=3
      ipars(2)=1
       wnear(8*nquad+1:9*nquad)=-wnear(2*nquad+1:3*nquad)

		  ipars(1)=3
      ipars(2)=2
       wnear(9*nquad+1:10*nquad)=-wnear(3*nquad+1:4*nquad)

		  ipars(1)=3
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(10*nquad+1:11*nquad))

		  ipars(1)=3
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(11*nquad+1:12*nquad))

		  ipars(1)=4
      ipars(2)=1
       wnear(12*nquad+1:13*nquad)=-wnear(6*nquad+1:7*nquad)

		  ipars(1)=4
      ipars(2)=2
       wnear(13*nquad+1:14*nquad)=-wnear(7*nquad+1:8*nquad)

		  ipars(1)=4
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(14*nquad+1:15*nquad))

		  ipars(1)=4
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(15*nquad+1:16*nquad))

!call system_clock ( t2, clock_rate, clock_max )
!write (*,*) 'time quadratures: ',real(t2-t1)/real(clock_rate)
!stop
    return
    end subroutine getnearquad_em_muller_trans_v2


      subroutine lpcomp_em_muller_trans_v2_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,wnear,&
     &n_components,contrast_matrix,npts_vect)

!
!  This subroutine evaluates the layer potential for
!  the DFIE boundary integral equation:
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  Note the 4\pi scaling is NOT included as the FMM output was scaled
!  appropriately
!
!  Note: the identities are not included as the gmres takes care of that
!
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
!
!
!  input:
!    npatches - integer
!      number of patches
!
!    norders- integer(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(npatches+1)
!      ixyzs(i) denotes the starting location in srccoefs,
!      and srcvals array corresponding to patch i
!   
!    iptype - integer(npatches)
!      type of patch
!      iptype = 1, triangular patch discretized using RV nodes
!
!    npts - integer
!      total number of discretization points on the boundary
! 
!    srccoefs - real *8 (9,npts)
!      koornwinder expansion coefficients of xyz, dxyz/du,
!      and dxyz/dv on each patch. 
!      For each point srccoefs(1:3,i) is xyz info
!        srccoefs(4:6,i) is dxyz/du info
!        srccoefs(7:9,i) is dxyz/dv info
!
!    srcvals - real *8 (12,npts)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
! 
!    ndtarg - integer
!      leading dimension of target array
!        
!    ntarg - integer
!      number of targets
!
!    targs - real *8 (ndtarg,ntarg)
!      target information
!
!    eps - real *8
!      precision requested
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = omega 
!      zpars(2) = ep0
!      zpars(3) = mu0
!      zpars(4) = ep1
!      zpars(5) = mu1
!
!    nnz - integer *8
!      number of source patch-> target interactions in the near field
! 
!    row_ptr - integer(ntarg+1)
!      row_ptr(i) is the pointer
!      to col_ind array where list of relevant source patches
!      for target i start
!
!    col_ind - integer (nnz)
!      list of source patches relevant for all targets, sorted
!      by the target number
!
!    iquad - integer(nnz+1)
!      location in wnear array where quadrature for col_ind(i)
!      starts
!
!    nquad - integer
!      number of entries in wnear
!
!    wnear  - complex *16(36*nquad)
!      near field precomputed quadrature
!
!    sigma - complex *16(2*ns)
!      induced charge and current on the surface
!      sigma(1:ns) - first component of 'a' along
!        the srcvals(4:6,i) direction
!      sigma(ns+1:2*ns) - second component of 'a' along
!        the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      sigma(2*ns+1:3*ns) - scalar sigma on the surface
!      sigma(3*ns+1:4*ns) - first component of 'b' along
!        the srcvals(4:6,i) direction
!      sigma(4*ns+1:5*ns) - second component of 'b' along
!        the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      sigma(5*ns+1:6*ns) - scalar rho on the surface
!
!    novers - integer(npatches)
!      order of discretization for oversampled sources and density
!
!    ixyzso - integer(npatches+1)
!      ixyzso(i) denotes the starting location in srcover,
!      corresponding to patch i
!   
!    nptso - integer
!      total number of oversampled points
!
!    srcover - real *8 (12,nptso)
!      oversampled set of source information
!
!    whtsover - real *8 (nptso)
!      smooth quadrature weights at oversampled nodes
!

      implicit none
      integer npatches,norder,npols,npts
      integer ndtarg,ntarg
      integer norders(npatches),ixyzs(npatches+1)
      integer ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(5),zpars_aux(5)
      integer nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      complex *16 sigma(4*npts),sigma2(npts)
      integer n_components
      integer npts_vect(n_components)
      complex *16 contrast_matrix(4,n_components)
	  
      complex *16 wnear(16*nquad)
	  
      integer novers(npatches+1)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(4*ntarg)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(4)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart,count1,count2,istart,ifinish


      integer ifaddsub,ifdir

      integer ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp2_s(:),dtmp2(:,:)
      complex *16, allocatable :: ctmp2_u(:),ctmp2_v(:)
      complex *16, allocatable :: ctmp2_a_u(:),ctmp2_a_v(:)
      complex *16, allocatable :: ctmp2_b_u(:),ctmp2_b_v(:)

	    complex *16, allocatable :: pot_aux(:)

      real *8 radexp,epsfmm

      integer ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      integer nss,ii,l,npover
      complex *16 ima
	    complex *16 omega,ep0,mu0,ep1,mu1

      integer nd,ntarg0
      integer icount

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)

      ns = nptso
      done = 1
      pi = atan(done)*4
	    ima=(0.0d0,1.0d0)

           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(3,ns),targvals(3,ntarg))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(4*ns))
	    allocate(pot_aux(4*ntarg))

! 
!       oversample density

    call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,&
	&npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))
	       
    call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,&
	&npts,sigma(npts+1:2*npts),novers,ixyzso,ns,sigmaover(ns+1:2*ns))

	call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
	&npts,sigma(2*npts+1:3*npts),novers,ixyzso,ns,sigmaover(2*ns+1:3*ns))

	call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
	&npts,sigma(3*npts+1:4*npts),novers,ixyzso,ns,sigmaover(3*ns+1:4*ns))

      ra = 0

!
!       fmm
!

		call get_fmm_thresh(12,ns,srcover,ndtarg,ntarg,targs,thresh)
       ifdir=0
	  
    istart=1
    ifinish=npts_vect(1)

    do count1=1,n_components
      zpars_aux(1)=zpars(1)
      zpars_aux(2)=contrast_matrix(1,count1)
      zpars_aux(3)=contrast_matrix(2,count1)
      zpars_aux(4)=contrast_matrix(3,count1)
      zpars_aux(5)=contrast_matrix(4,count1)

      !Calculate the far_field with FMM		
      call em_muller_trans_FMM(eps,zpars_aux,ns,npts_vect(count1),&
      &srcover,targs(:,istart:ifinish),whtsover,&
      &sigmaover(1:ns),sigmaover(ns+1:2*ns),sigmaover(2*ns+1:3*ns),&
      &sigmaover(3*ns+1:4*ns),pot_aux(istart:ifinish),&
      &pot_aux(istart+ntarg:ifinish+ntarg),&
      &pot_aux(istart+2*ntarg:ifinish+2*ntarg),&
      &pot_aux(istart+3*ntarg:ifinish+3*ntarg),thresh,ifdir)

      if (count1<n_components) then
        istart=ifinish+1
        ifinish=istart+npts_vect(count1+1)-1
      endif
    enddo

		call cpu_time(t1)
!C$      t1 = omp_get_wtime()
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
!C$OMP$PRIVATE(jstart,pottmp,npols)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols	
		        do count1=0,3
			       do count2=0,3
			         pot_aux(i+count1*npts) = pot_aux(i+count1*npts) + &
				      &wnear((count1*4+count2)*nquad+jquadstart+l-1)*&
				      &sigma(jstart+l-1+npts*count2)
			       enddo
			      enddo
          enddo
        enddo
      enddo

!C$OMP END PARALLEL DO
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
!		ipars(1)=1
!		ipars(2)=1
	ifdir=1
!      do i=1,ntarg
i=0
      do count1=1,n_components
       zpars_aux(1)=zpars(1)
       zpars_aux(2)=contrast_matrix(1,count1)
       zpars_aux(3)=contrast_matrix(2,count1)
       zpars_aux(4)=contrast_matrix(3,count1)
       zpars_aux(5)=contrast_matrix(4,count1)
       do count2=1,npts_vect(count1)
       i=i+1
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),wtmp2(nss))
        allocate(ctmp2_a_u(nss),ctmp2_a_v(nss))
	      allocate(ctmp2_b_u(nss),ctmp2_b_v(nss))

        rmin = 1.0d6
        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
			      ii = ii+1
			      srctmp2(:,ii) = srcover(:,jstart+l)
			      ctmp2_a_u(ii)=sigmaover(jstart+l)
			      ctmp2_a_v(ii)=sigmaover(jstart+l+ns)			
			      ctmp2_b_u(ii)=sigmaover(jstart+l+2*ns)
			      ctmp2_b_v(ii)=sigmaover(jstart+l+3*ns)
			      wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo
        call em_muller_trans_FMM(eps,zpars_aux,nss,ntarg0,&
        &srctmp2,targs(:,i),wtmp2,ctmp2_a_u,ctmp2_a_v,&
        &ctmp2_b_u,ctmp2_b_v,E(1),E(2),E(3),E(4),thresh,ifdir)
	 
        do j=0,3
          pot_aux(i+j*ntarg) = pot_aux(i+j*ntarg) - E(j+1)
        enddo
		
        deallocate(srctmp2,wtmp2)
        deallocate(ctmp2_a_u,ctmp2_a_v,ctmp2_b_u,ctmp2_b_v)

       enddo
      enddo

    omega = zpars(1)
	  ep0 = zpars(2)
    mu0 = zpars(3)
	  ep1 = zpars(4)
	  mu1 = zpars(5)

    i=0
      do count1=1,n_components
       zpars_aux(1)=zpars(1)
       ep0=contrast_matrix(1,count1)
       mu0=contrast_matrix(2,count1)
       ep1=contrast_matrix(3,count1)
       mu1=contrast_matrix(4,count1)
       do count2=1,npts_vect(count1)
       i=i+1
		   pot_aux(i)=pot_aux(i)/(mu0+mu1)
		   pot_aux(i+ntarg)=pot_aux(i+ntarg)/(mu0+mu1)
		   pot_aux(i+2*ntarg)=pot_aux(i+2*ntarg)/(ep0+ep1)
		   pot_aux(i+3*ntarg)=pot_aux(i+3*ntarg)/(ep0+ep1)
       enddo
      enddo
	  
	  do i=1,ntarg
		  pot(i)=pot_aux(i)
		  pot(i+ntarg)=pot_aux(i+ntarg)
		  pot(i+2*ntarg)=pot_aux(i+2*ntarg)
		  pot(i+3*ntarg)=pot_aux(i+3*ntarg)
	  enddo
	  
	  return
    end subroutine lpcomp_em_muller_trans_v2_addsub



subroutine em_muller_trans_v2_solver(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,&
     &rhs,eps_gmres,niter,errs,rres,soln,contrast_matrix,&
     &npts_vect,n_components,srcvals_extended)
!
!  This subroutine solves the Scattering Maxwell homogeneous dielectric
!  problem.
!  The the equations are:
!
!    curlE0=ik0H0; curlH0 =-ik0E0   (exterior region)
!    curlE1=ik1H1; curlH1 =-ik1E1   (interior region)
!
!  parameters: k0=omega*sqrt(ep0*mu0); k1=omega*sqrt(ep1*mu1)
!
!  Representation: (1)
!
!    E0 = mu0curlS_{k0}[a] - mu0S_{k0}[n·sigma] + mu0ep0S_{k0}[b] + 
!       + gradS_{k0}[rho]
!    E1 = mu1curlS_{k1}[a] - mu1S_{k1}[n·sigma] + mu1ep1S_{k1}[b] + 
!       + gradS_{k1}[rho]
!
!  Boundary conditions: (2)
!
!                    nxE0-nxE1 = -nxE_inc
!                  divE0-divE1 = 0
!    nxcurlE0/mu0-nxcurlE1/mu1 = -nxcurlE_inc/mu0
!              n·(ep0E0-ep1E1) = -n·ep0·E_inc
!
!
!  The incoming fields must be 'compatible' 
!  (solutions to Maxwell's equations in free space)
!
!  Boundary integral equation:
!
!    equation (86)of paper https://doi.org/10.1080/03605302.2018.1446447
!    is the result of hitting the representation (1) with the boundary 
!    conditions (2) and rescale the the equations with 1/(mu0+mu1) 
!    and 1/(ep0+ep1) to get 1/2 in the diagonal
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of 1e-15 is reached
!
!  input:
!    npatches - integer
!      number of patches
!
!    norders- integer(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(npatches+1)
!      ixyzs(i) denotes the starting location in srccoefs,
!      and srcvals array corresponding to patch i
!   
!    iptype - integer(npatches)
!      type of patch
!      iptype = 1, triangular patch discretized using RV nodes
!
!    npts - integer
!      total number of discretization points on the boundary
! 
!    srccoefs - real *8 (9,npts)
!      koornwinder expansion coefficients of xyz, dxyz/du,
!      and dxyz/dv on each patch. 
!      For each point srccoefs(1:3,i) is xyz info
!        srccoefs(4:6,i) is dxyz/du info
!        srccoefs(7:9,i) is dxyz/dv info
!
!    srcvals - real *8 (12,npts)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
! 
!    eps - real *8
!      precision requested for computing quadrature and fmm
!      tolerance
!
!    zpars - complex *16 (5)
!      kernel parameters 
!      zpars(1) = omega 
!      zpars(2) = ep0
!      zpars(3) = mu0
!      zpars(4) = ep1
!      zpars(5) = mu1
!
!    ifinout - integer
!      flag for interior or exterior problems (normals assumed to 
!        be pointing in exterior of region)
!      ifinout = 0, interior problem
!      ifinout = 1, exterior problem
!
!    rhs - complex *16(npts)
!      right hand side
!
!    eps_gmres - real *8
!      gmres tolerance requested
!
!    numit - integer
!      max number of gmres iterations
!
!    output
!      niter - integer
!      number of gmres iterations required for relative residual
!      to converge to 1e-15
!          
!    errs(1:iter) - relative residual as a function of iteration
!      number
! 
!    rres - real *8
!      relative residual for computed solution
!              
!    soln - complex *16(3*npts)
!      soln(1:npts) component of the tangent induced current a on the 
!        surface along srcvals(4:6,i) direction
!      soln(npts+1:2*npts) component of the tangent induced current a 
!        on the surface along (srcvals(10:12,i) x srcvals(4:6,i)) 
!        direction
!      soln(2*npts+1:3*npts) rho on the surface 
!      soln(3*npts+1:4*npts) component of the tangent induced current
!        a on the surface along srcvals(4:6,i) direction
!      soln(4*npts+1:5*npts) component of the tangent induced current
!        a on the surface along (srcvals(10:12,i) x srcvals(4:6,i))
!        direction
!      soln(5*npts+1:6*npts) rho on the surface 
!

      implicit none

      integer, intent(in) :: n_components
      integer, intent(in) :: npts_vect(n_components)
      complex *16, intent(in) :: contrast_matrix(4,n_components)
      real *8, intent(in) :: srcvals_extended(20,npts)

      integer npatches,norder,npols,npts
      integer ifinout
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(5)
      complex *16 rhs(4*npts)
      complex *16 soln(4*npts)

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg

      real *8 errs(numit+1)
      real *8 rres,eps2
      integer niter


      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)

      complex *16, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars(2)
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over
	    integer n_var

!
!
!       gmres variables
!
      complex *16 zid,ztmp
      real *8 rb,wnrm2
      integer numit,it,iind,it1,k,l,count1
      real *8 rmyerr
      complex *16 temp
      complex *16, allocatable :: vmat(:,:),hmat(:,:)
      complex *16, allocatable :: cs(:),sn(:)
      complex *16, allocatable :: svec(:),yvec(:),wtmp(:)
      complex *16 ima

      ima=(0.0d0,1.0d0)
	  
!
!   n_var is the number of unknowns in the linear system.
!   as we have tow vector unknown a,b and two scalars rho,sigma
!   we need n_var=6*npts
!

      n_var=4*npts

      allocate(vmat(n_var,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(n_var),svec(numit+1),yvec(numit+1))


      done = 1
      pi = atan(done)*4


!
!        setup targets as on surface discretization points
! 
      ndtarg = 12
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))
!C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(:,i)=srcvals(:,i)
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!C$OMP END PARALLEL DO   


!
!    initialize patch_id and uv_targ for on surface targets
!
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,&
     &ipatch_id,uvs_targ)

!
!
!        this might need fixing
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,& 
     &srccoefs,cms,rads)

!C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!C$OMP END PARALLEL DO      

!
!    find near quadrature correction interactions
!
      print *, "entering find near mem"
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr,&
     &col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,&
     &iquad)

      ikerorder = 0

!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,&
     &rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars(1),&
     &nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1
      print *, "npts_over=",npts_over

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts,&
     &srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,&
     &srcover,wover)


!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      allocate(wnear(16*nquad))
	  
!C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,16*nquad
	    wnear(i)=0
      enddo
!C$OMP END PARALLEL DO    

      iquadtype = 1

!!      eps2 = 1.0d-8

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      
ndtarg=20
      call getnearquad_em_muller_trans_v2(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals_extended,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear)
      ndtarg=12
      call cpu_time(t2)
!C$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time=*',t2-t1,1)
      
      print *, "done generating near quadrature, now starting gmres"

!
!
!     start gmres code here
!
!     NOTE: matrix equation should be of the form (z*I + K)x = y
!       the identity scaling (z) is defined via zid below,
!       and K represents the action of the principal value 
!       part of the matvec
!

	  zid=0.5d0


      niter=0

!
!      compute norm of right hand side and initialize v
! 
      rb = 0

      do i=1,numit
        cs(i) = 0
        sn(i) = 0
      enddo


!
      do i=1,n_var
        rb = rb + abs(rhs(i))**2
      enddo
      rb = sqrt(rb)

      do i=1,n_var
        vmat(i,1) = rhs(i)/rb
      enddo

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

!
!        NOTE:
!        replace this routine by appropriate layer potential
!        evaluation routine  
!

!write (*,*) 'about to multiply'
        call lpcomp_em_muller_trans_v2_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
     &vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp,wnear,&
     &n_components,contrast_matrix,npts_vect)
!write (*,*) 'done first multiplication'

        do k=1,it
          hmat(k,it) = 0
          do j=1,n_var      
            hmat(k,it) = hmat(k,it) + wtmp(j)*conjg(vmat(j,k))
          enddo

          do j=1,n_var
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
        wnrm2 = 0
        do j=1,n_var
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
        wnrm2 = sqrt(wnrm2)

        do j=1,n_var
          vmat(j,it1) = wtmp(j)/wnrm2
        enddo

        do k=1,it-1
          temp = cs(k)*hmat(k,it)+sn(k)*hmat(k+1,it)
          hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
          hmat(k,it) = temp
        enddo

        ztmp = wnrm2

        call zrotmat_gmres(hmat(it,it),ztmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+sn(it)*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)

        if(rmyerr.le.eps_gmres.or.it.eq.numit) then

!
!            solve the linear system corresponding to
!            upper triangular part of hmat to obtain yvec
!
!            y = triu(H(1:it,1:it))\s(1:it);
!
          do j=1,it
            iind = it-j+1
            yvec(iind) = svec(iind)
            do l=iind+1,it
              yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
            enddo
            yvec(iind) = yvec(iind)/hmat(iind,iind)
          enddo

!
!          estimate x
!
          do j=1,n_var
            soln(j) = 0
            do i=1,it
              soln(j) = soln(j) + yvec(i)*vmat(j,i)
            enddo
          enddo


          rres = 0
          do i=1,n_var
            wtmp(i) = 0
          enddo
!
!        NOTE:
!        replace this routine by appropriate layer potential
!        evaluation routine  
!


          call lpcomp_em_muller_trans_v2_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
     &soln,novers,npts_over,ixyzso,srcover,wover,wtmp,wnear,&
     &n_components,contrast_matrix,npts_vect)

            
          do i=1,npts
            rres = rres + abs(zid*soln(i) + wtmp(i)-rhs(i))**2
          enddo
          rres = sqrt(rres)/rb
          niter = it
          return
        endif
      enddo
!
      return
      end subroutine em_muller_trans_v2_solver




subroutine fker_em_muller_trans_v2(srcinfo, ndt,targinfo,ndd, dpars,ndz,&
 &zpars,ndi,ipars,E_val)
implicit none

!  THIS FUNCTION IS OBSOLETE, NOT USED, subroutine em_dfie_trans is ued
!  instead (much faster), but very hard to read..
!  this function provides the near field kernel that will use 
!  zgetnearquad_ggq_guru 
!  through getnearquad_DFIE

    !List of calling arguments
	integer, intent(in) :: ndt,ndd,ndz,ndi
	real ( kind = 8 ), intent(in) :: srcinfo(12)
	real ( kind = 8 ), intent(in) :: targinfo(ndt)
	integer, intent(in) :: ipars(ndi)
	real ( kind = 8 ), intent(in) :: dpars(ndd)
	complex ( kind = 8 ), intent(in) :: zpars(ndz)
	complex ( kind = 8 ), intent(out) :: E_val
	
	!List of local variables
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3),aux_real
	complex ( kind = 8 ) nxcurlSk0a(2,2),nxcurlSk1a(2,2),nxcurlcurlSk0a(2,2)
	complex ( kind = 8 ) nxcurlcurlSk1a(2,2),ncurlSk1a(1,2),ncurlSk0a(1,2)
	complex ( kind = 8 ) nxSk0b(2,2),nxSk1b(2,2),divSk0b(1,2),divSk1b(1,2)
	complex ( kind = 8 ) nxSk0nrho(2,1),nxSk1nrho(2,1),nxcurlSk0nrho(2,1)
	complex ( kind = 8 ) nxcurlSk1nrho(2,1),Dk0rho(1,1),Dk1rho(1,1)
	complex ( kind = 8 ) nxgradSk0lambda(2,1),nxgradSk1lambda(2,1),Sk0lambda(1,1)
	complex ( kind = 8 ) Sk1lambda(1,1),ngradSk0lambda(1,1),nSk1nrho(1,1)
	complex ( kind = 8 ) ngradSk1lambda(1,1),nSk0b(1,2),nSk1b(1,2),nSk0nrho(1,1)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima,my_exp_0,my_exp_1,omega,ep0	
	complex ( kind = 8 ) mu0,ep1,mu1,zk0,zk1
	real ( kind = 8 ) pi
	integer count1,count2
	complex ( kind = 8 )  E_mat(6,6)

	
	pi=3.1415926535897932384626433832795028841971d0
	ima=(0.0d0,1.0d0)
	omega=zpars(1)
	
	sour(1)=srcinfo(1)
	sour(2)=srcinfo(2)
	sour(3)=srcinfo(3)
	
	n_s(1)=srcinfo(10)
	n_s(2)=srcinfo(11)
	n_s(3)=srcinfo(12)	

	targ(1)=targinfo(1)
	targ(2)=targinfo(2)
	targ(3)=targinfo(3)

	n_t(1)=targinfo(10)
	n_t(2)=targinfo(11)
	n_t(3)=targinfo(12)

	dr(1)=targ(1)-sour(1)
	dr(2)=targ(2)-sour(2)
	dr(3)=targ(3)-sour(3)

    ep0=targinfo(13)+ima*targinfo(14)
    mu0=targinfo(15)+ima*targinfo(16)
	ep1=targinfo(17)+ima*targinfo(18)
    mu1=targinfo(19)+ima*targinfo(20)
	
	r=sqrt((dr(1))**2+(dr(2))**2+(dr(3))**2)
	zk0=omega*sqrt(ep0*mu0)
	zk1=omega*sqrt(ep1*mu1)

	R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
	R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*exp(ima*zk0*r)/&
	 &(4.0d0*pi)
	my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
	R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
	R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*exp(ima*zk1*r)/&
	 &(4.0d0*pi)
	my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	
	call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
	call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)

	call get_nxcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,nxcurlSk0a)		
	call get_nxcurlcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,R2_0,zk0,&
	 &my_exp_0,r,nxcurlcurlSk0a)
		
	call get_nxcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,nxcurlSk1a)		
	call get_nxcurlcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,R2_1,zk1,&
	 &my_exp_1,r,nxcurlcurlSk1a)


		E_mat(1,1)=mu0*nxcurlSk0a(1,1)-mu1*nxcurlSk1a(1,1)
		E_mat(1,2)=mu0*nxcurlSk0a(1,2)-mu1*nxcurlSk1a(1,2)
		E_mat(2,1)=mu0*nxcurlSk0a(2,1)-mu1*nxcurlSk1a(2,1)
		E_mat(2,2)=mu0*nxcurlSk0a(2,2)-mu1*nxcurlSk1a(2,2)

		E_mat(3,3)=ep0*nxcurlSk0a(1,1)-ep1*nxcurlSk1a(1,1)
		E_mat(3,4)=ep0*nxcurlSk0a(1,2)-ep1*nxcurlSk1a(1,2)
		E_mat(4,3)=ep0*nxcurlSk0a(2,1)-ep1*nxcurlSk1a(2,1)
		E_mat(4,4)=ep0*nxcurlSk0a(2,2)-ep1*nxcurlSk1a(2,2)

		E_mat(1,3)=(nxcurlcurlSk0a(1,1)-nxcurlcurlSk1a(1,1))/(-ima*omega)
		E_mat(1,4)=(nxcurlcurlSk0a(1,2)-nxcurlcurlSk1a(1,2))/(-ima*omega)
		E_mat(2,3)=(nxcurlcurlSk0a(2,1)-nxcurlcurlSk1a(2,1))/(-ima*omega)
		E_mat(2,4)=(nxcurlcurlSk0a(2,2)-nxcurlcurlSk1a(2,2))/(-ima*omega)
		
		E_mat(3,1)=-E_mat(1,3)
		E_mat(3,2)=-E_mat(1,4)
		E_mat(4,1)=-E_mat(2,3)
		E_mat(4,2)=-E_mat(2,4)
		
	E_val=E_mat(ipars(1),ipars(2))


return
end subroutine fker_em_muller_trans_v2



subroutine em_muller_trans_FMM(eps,zpars,ns,nt,srcvals,targvals,wts,&
 &a_u,a_v,b_u,b_v,AA_u,AA_v,BB_u,BB_v,thresh,ifdir)
implicit none

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zpars(5)
    integer, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns)
    real ( kind = 8 ), intent(in) :: targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns)
    complex ( kind = 8 ), intent(in) :: b_u(ns),b_v(ns)
    complex ( kind = 8 ), intent(out) :: AA_u(nt),AA_v(nt)
    complex ( kind = 8 ), intent(out) :: BB_u(nt),BB_v(nt)
    real ( kind = 8 ), intent(in) :: thresh
    integer, intent(in) :: ifdir 


    !List of local variables
	  real ( kind = 8 ), allocatable :: u_vect_s(:,:),v_vect_s(:,:)
	  real ( kind = 8 ), allocatable :: source(:,:),n_vect_s(:,:)

	  real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:)
	  real ( kind = 8 ), allocatable :: targets(:,:),v_vect_t(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:),rho(:),b_vect_t(:,:)
    complex ( kind = 8 ), allocatable :: E(:,:),curlE(:,:),divE(:)
    complex ( kind = 8 ) ima,zk0,zk1

    complex ( kind = 8 ) omega,ep0,mu0,ep1,mu1

    integer count1,count2
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

	  ima=(0.0d0,1.0d0)

	  omega=zpars(1)
	  ep0=zpars(2)
	  mu0=zpars(3)
	  ep1=zpars(4)
	  mu1=zpars(5)
	
	  zk0=omega*sqrt(ep0*mu0)
	  zk1=omega*sqrt(ep1*mu1)


    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(b_vect_t(3,nt))
    allocate(lambda(ns))
	  allocate(rho(ns))
    allocate(E(3,nt))
    allocate(curlE(3,nt))
    allocate(divE(nt))
	  allocate(n_vect_s(3,ns))
	  allocate(n_vect_t(3,nt))
	  allocate(u_vect_s(3,ns))
	  allocate(v_vect_s(3,ns))
	  allocate(u_vect_t(3,nt))
	  allocate(v_vect_t(3,nt))
	  allocate(source(3,ns))
	  allocate(targets(3,nt))

	  do count1=1,ns
      n_vect_s(:,count1)=srcvals(10:12,count1)
      source(:,count1)=srcvals(1:3,count1)
	  enddo
	  call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect_s,&
	   &v_vect_s,ns)
	
	  do count1=1,nt
      n_vect_t(:,count1)=targvals(10:12,count1)
      targets(:,count1)=targvals(1:3,count1)
	  enddo
	  call orthonormalize_all(targvals(4:6,:),targvals(10:12,:),u_vect_t,&
	   &v_vect_t,nt)
	


    do count1=1,ns
		  a_vect(1,count1)=(b_u(count1)*u_vect_s(1,count1)+b_v(count1)*v_vect_s(1,count1))/(-ima*omega)
		  a_vect(2,count1)=(b_u(count1)*u_vect_s(2,count1)+b_v(count1)*v_vect_s(2,count1))/(-ima*omega)
		  a_vect(3,count1)=(b_u(count1)*u_vect_s(3,count1)+b_v(count1)*v_vect_s(3,count1))/(-ima*omega)
				
      b_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*v_vect_s(1,count1))*mu0
      b_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*v_vect_s(2,count1))*mu0
      b_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*v_vect_s(3,count1))*mu0
	  enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

	  call Vector_Helmholtz_targ2(eps,zk0,ns,source,wts,ifa_vect,a_vect,&
	   &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	   &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

!	call Vector_Helmholtz_targ(eps,izk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
!	 &b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,curlE,ifdivE,divE,nt,targets)

    do count1=1,nt
      b_vect_t(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-n_vect_t(3,count1)*curlE(2,count1)
      b_vect_t(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-n_vect_t(1,count1)*curlE(3,count1)
      b_vect_t(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
      AA_u(count1)=b_vect_t(1,count1)*u_vect_t(1,count1)+b_vect_t(2,count1)*u_vect_t(2,count1)&
		   &+b_vect_t(3,count1)*u_vect_t(3,count1)
      AA_v(count1)=b_vect_t(1,count1)*v_vect_t(1,count1)+b_vect_t(2,count1)*v_vect_t(2,count1)&
		   &+b_vect_t(3,count1)*v_vect_t(3,count1)
    enddo

	  do count1=1,ns
		  b_vect(1,count1)=(b_u(count1)*u_vect_s(1,count1)+b_v(count1)*v_vect_s(1,count1))*ep0
		  b_vect(2,count1)=(b_u(count1)*u_vect_s(2,count1)+b_v(count1)*v_vect_s(2,count1))*ep0
		  b_vect(3,count1)=(b_u(count1)*u_vect_s(3,count1)+b_v(count1)*v_vect_s(3,count1))*ep0
				
      a_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*v_vect_s(1,count1))/(ima*omega)
      a_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*v_vect_s(2,count1))/(ima*omega)
      a_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*v_vect_s(3,count1))/(ima*omega)
	  enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

	  call Vector_Helmholtz_targ2(eps,zk0,ns,source,wts,ifa_vect,a_vect,&
	   &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	   &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

    do count1=1,nt
      b_vect_t(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-n_vect_t(3,count1)*curlE(2,count1)
      b_vect_t(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-n_vect_t(1,count1)*curlE(3,count1)
      b_vect_t(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
      BB_u(count1)=b_vect_t(1,count1)*u_vect_t(1,count1)+b_vect_t(2,count1)*u_vect_t(2,count1)&
       &+b_vect_t(3,count1)*u_vect_t(3,count1)
      BB_v(count1)=b_vect_t(1,count1)*v_vect_t(1,count1)+b_vect_t(2,count1)*v_vect_t(2,count1)&
       &+b_vect_t(3,count1)*v_vect_t(3,count1)
    enddo


!!! now zk1


    do count1=1,ns
		  a_vect(1,count1)=(b_u(count1)*u_vect_s(1,count1)+b_v(count1)*v_vect_s(1,count1))/(-ima*omega)
		  a_vect(2,count1)=(b_u(count1)*u_vect_s(2,count1)+b_v(count1)*v_vect_s(2,count1))/(-ima*omega)
		  a_vect(3,count1)=(b_u(count1)*u_vect_s(3,count1)+b_v(count1)*v_vect_s(3,count1))/(-ima*omega)
				
      b_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*v_vect_s(1,count1))*mu1
      b_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*v_vect_s(2,count1))*mu1
      b_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*v_vect_s(3,count1))*mu1
    enddo

    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

	  call Vector_Helmholtz_targ2(eps,zk1,ns,source,wts,ifa_vect,a_vect,&
	   &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	   &curlE,ifdivE,divE,nt,targets,thresh,ifdir)


    do count1=1,nt
        b_vect_t(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-n_vect_t(3,count1)*curlE(2,count1)
        b_vect_t(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-n_vect_t(1,count1)*curlE(3,count1)
        b_vect_t(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
        AA_u(count1)=AA_u(count1)-(b_vect_t(1,count1)*u_vect_t(1,count1)+b_vect_t(2,count1)*u_vect_t(2,count1)&
		     &+b_vect_t(3,count1)*u_vect_t(3,count1))
        AA_v(count1)=AA_v(count1)-(b_vect_t(1,count1)*v_vect_t(1,count1)+b_vect_t(2,count1)*v_vect_t(2,count1)&
		     &+b_vect_t(3,count1)*v_vect_t(3,count1))
    enddo


	  do count1=1,ns
		    b_vect(1,count1)=(b_u(count1)*u_vect_s(1,count1)+b_v(count1)*v_vect_s(1,count1))*ep1
		    b_vect(2,count1)=(b_u(count1)*u_vect_s(2,count1)+b_v(count1)*v_vect_s(2,count1))*ep1
		    b_vect(3,count1)=(b_u(count1)*u_vect_s(3,count1)+b_v(count1)*v_vect_s(3,count1))*ep1
				
        a_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*v_vect_s(1,count1))/(ima*omega)
        a_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*v_vect_s(2,count1))/(ima*omega)
        a_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*v_vect_s(3,count1))/(ima*omega)
	  enddo


    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

	  call Vector_Helmholtz_targ2(eps,zk1,ns,source,wts,ifa_vect,a_vect,&
	   &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	   &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

    do count1=1,nt
        b_vect_t(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-n_vect_t(3,count1)*curlE(2,count1)
        b_vect_t(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-n_vect_t(1,count1)*curlE(3,count1)
        b_vect_t(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
        BB_u(count1)=BB_u(count1)-(b_vect_t(1,count1)*u_vect_t(1,count1)+b_vect_t(2,count1)*u_vect_t(2,count1)&
		     &+b_vect_t(3,count1)*u_vect_t(3,count1))
        BB_v(count1)=BB_v(count1)-(b_vect_t(1,count1)*v_vect_t(1,count1)+b_vect_t(2,count1)*v_vect_t(2,count1)&
		     &+b_vect_t(3,count1)*v_vect_t(3,count1))
    enddo

	
	deallocate(a_vect)
	deallocate(b_vect)
  deallocate(b_vect_t)
	deallocate(lambda)
	deallocate(rho)
	deallocate(E)
	deallocate(curlE)
	deallocate(divE)
	deallocate(u_vect_s)
	deallocate(v_vect_s)
	deallocate(n_vect_s)
	deallocate(source)
	deallocate(u_vect_t)
	deallocate(v_vect_t)
	deallocate(n_vect_t)
	deallocate(targets)

return
end subroutine em_muller_trans_FMM


subroutine 	get_rhs_em_muller_trans_v2(P0,vf,srcvals,omega,RHS,n_components,npts_vect,contrast_matrix,ns)
implicit none

	!List of calling arguments
  integer, intent(in) :: n_components,ns
  integer, intent(in) :: npts_vect(n_components)
  complex (kind = 8 ), intent(in) :: contrast_matrix(4,n_components)
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
  complex ( kind = 8 ), intent(in) :: omega
	complex ( kind = 8 ), intent(out) :: RHS(4*ns)
	
	!List of local variables
  complex ( kind = 8 ) ep0,mu0,ep1,mu1
	complex ( kind = 8 ), allocatable :: E0(:,:), H0(:,:),E(:,:), H(:,:)
	integer count1,count2,icount
	real ( kind = 8 ), allocatable :: xyz_aux(:,:)
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)

 	allocate(E(3,ns), H(3,ns))
	allocate(E0(3,ns), H0(3,ns))
  allocate(xyz_aux(3,ns))

  ep0=contrast_matrix(1,1)
  mu0=contrast_matrix(2,1)	
	call fieldsEDomega(omega,ep0,mu0,P0,srcvals,ns,E0,H0,vf,0)
	call fieldsMDomega(omega,ep0,mu0,P0,srcvals,ns,E0,H0,vf,1)

!	call fieldsEDomega(omega,ep1,mu1,Pt,srcvals,ns,E,H,vf_minus,1)
!	call fieldsMDomega(omega,ep1,mu1,Pt,srcvals,ns,E,H,vf_minus,1)
!	read (*,*)

  do count1=1,ns
    xyz_aux(1,count1)=srcvals(1,count1)
    xyz_aux(2,count1)=srcvals(2,count1)
    xyz_aux(3,count1)=srcvals(3,count1)
  enddo

  icount=1
  do count2=1,n_components
    ep0=contrast_matrix(1,count2)
    mu0=contrast_matrix(2,count2)
    ep1=contrast_matrix(3,count2)
    mu1=contrast_matrix(4,count2)
    call fieldsPWomega(omega,ep1,mu1,xyz_aux,ns,E,H)
!    E = 0
!    H = 0
	  do count1=1,npts_vect(count2)	
		  call orthonormalize(srcvals(4:6,icount),srcvals(10:12,icount),ru,rv)		
		  RHS(icount)=-DOT_PRODUCT(rv,E0(:,icount)-E(:,icount))
		  RHS(ns+icount)=DOT_PRODUCT(ru,E0(:,icount)-E(:,icount))	
		  RHS(2*ns+icount)=-DOT_PRODUCT(rv,H0(:,icount)-H(:,icount))
		  RHS(3*ns+icount)=DOT_PRODUCT(ru,H0(:,icount)-H(:,icount))
    
      RHS(icount)=RHS(icount)/(mu0+mu1)
      RHS(icount+ns)=RHS(icount+ns)/(mu0+mu1)
      RHS(icount+2*ns)=RHS(icount+2*ns)/(ep0+ep1)
      RHS(icount+3*ns)=RHS(icount+3*ns)/(ep0+ep1)
      icount=icount+1
	  enddo
  enddo

return
end subroutine get_rhs_em_muller_trans_v2




subroutine 	get_rhs_em_muller_trans_testing(P0,vf,direction,Pol,srcvals,omega,&
 &RHS,n_components,npts_vect,contrast_matrix,ns,exposed_surfaces)
implicit none

	!List of calling arguments
  integer, intent(in) :: n_components,ns
  integer, intent(in) :: npts_vect(n_components)
  complex (kind = 8 ), intent(in) :: contrast_matrix(4,n_components)
	real ( kind = 8 ), intent(in) :: direction(2),P0(3)
	complex ( kind = 8 ), intent(in) :: Pol(2),vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
  complex ( kind = 8 ), intent(in) :: omega
	complex ( kind = 8 ), intent(out) :: RHS(4*ns)
  logical *8, intent(in) :: exposed_surfaces(n_components)
	
	!List of local variables
  complex ( kind = 8 ) ep0,mu0,ep1,mu1
	complex ( kind = 8 ), allocatable :: E0(:,:), H0(:,:),E(:,:), H(:,:),E2(:,:), H2(:,:)
	integer count1,count2,icount
  real ( kind = 8 ), allocatable :: xyz_aux(:,:)
	
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)

 	allocate(E(3,ns), H(3,ns),E2(3,ns), H2(3,ns))
	allocate(E0(3,ns), H0(3,ns))
  allocate(xyz_aux(3,ns))
		
  do count1=1,n_components
    if (exposed_surfaces(count1)) then
        ep0=contrast_matrix(1,count1)
        mu0=contrast_matrix(2,count1)
        exit
    endif
  enddo

  do count1=1,ns
    xyz_aux(1,count1)=srcvals(1,count1)
    xyz_aux(2,count1)=srcvals(2,count1)
    xyz_aux(3,count1)=srcvals(3,count1)
  enddo

! call fieldsPWomega(omega,ep0,mu0,direction,Pol,srcvals,ns,E0,H0)
	call fieldsEDomega(omega,ep0,mu0,P0,srcvals,ns,E0,H0,vf,0)   
	call fieldsMDomega(omega,ep0,mu0,P0,srcvals,ns,E0,H0,vf,1)

icount=1
do count2=1,n_components
  ep0=contrast_matrix(1,count2)
  mu0=contrast_matrix(2,count2)
  ep1=contrast_matrix(3,count2)
  mu1=contrast_matrix(4,count2)
  call fieldsPWomega(omega,ep1,mu1,xyz_aux,ns,E,H,direction,Pol)
  call fieldsPWomega(omega,ep0,mu0,xyz_aux,ns,E2,H2,direction,Pol)
	do count1=1,npts_vect(count2)	
		call orthonormalize(srcvals(4:6,icount),srcvals(10:12,icount),ru,rv)
    if (exposed_surfaces(count2)) then
		  RHS(icount)=-DOT_PRODUCT(rv,E0(:,icount)-E(:,icount))
		  RHS(ns+icount)=DOT_PRODUCT(ru,E0(:,icount)-E(:,icount))	
		  RHS(2*ns+icount)=-DOT_PRODUCT(rv,H0(:,icount)-H(:,icount))
		  RHS(3*ns+icount)=DOT_PRODUCT(ru,H0(:,icount)-H(:,icount))
    else
		  RHS(icount)=-DOT_PRODUCT(rv,E2(:,icount)-E(:,icount))
		  RHS(ns+icount)=DOT_PRODUCT(ru,E2(:,icount)-E(:,icount))	
		  RHS(2*ns+icount)=-DOT_PRODUCT(rv,H2(:,icount)-H(:,icount))
		  RHS(3*ns+icount)=DOT_PRODUCT(ru,H2(:,icount)-H(:,icount))
    endif
    RHS(icount)=RHS(icount)/(mu0+mu1)
    RHS(icount+ns)=RHS(icount+ns)/(mu0+mu1)
    RHS(icount+2*ns)=RHS(icount+2*ns)/(ep0+ep1)
    RHS(icount+3*ns)=RHS(icount+3*ns)/(ep0+ep1)
    icount=icount+1
	enddo
enddo

return
end subroutine get_rhs_em_muller_trans_testing




subroutine 	get_rhs_em_muller_trans_PW(direction,Pol,srcvals,omega,RHS,&
 &n_components,npts_vect,contrast_matrix,ns,exposed_surfaces)
implicit none

	!List of calling arguments
  integer, intent(in) :: n_components,ns
  integer, intent(in) :: npts_vect(n_components)
  complex (kind = 8 ), intent(in) :: contrast_matrix(4,n_components)
	real ( kind = 8 ), intent(in) :: direction(2)
	complex ( kind = 8 ), intent(in) :: Pol(2)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
  complex ( kind = 8 ), intent(in) :: omega
	complex ( kind = 8 ), intent(out) :: RHS(4*ns)
  logical *8, intent(in) :: exposed_surfaces(n_components)
	
	!List of local variables
  complex ( kind = 8 ) ep0,mu0,ep1,mu1
	complex ( kind = 8 ), allocatable :: E0(:,:), H0(:,:),E(:,:), H(:,:)
	integer count1,count2,icount
  real ( kind = 8 ), allocatable :: xyz_aux(:,:)
	
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)

 	allocate(E(3,ns), H(3,ns))
	allocate(E0(3,ns), H0(3,ns))
  allocate(xyz_aux(3,ns))
		
  do count1=1,n_components
    if (exposed_surfaces(count1)) then
        ep0=contrast_matrix(1,count1)
        mu0=contrast_matrix(2,count1)
        exit
    endif
  enddo

  do count1=1,ns
    xyz_aux(1,count1)=srcvals(1,count1)
    xyz_aux(2,count1)=srcvals(2,count1)
    xyz_aux(3,count1)=srcvals(3,count1)
  enddo
  call fieldsPWomega(omega,ep0,mu0,xyz_aux,ns,E0,H0,direction,Pol)
!	call fieldsEDomega(omega,ep1,mu1,Pt,srcvals,ns,E,H,vf_minus,1)
!	call fieldsMDomega(omega,ep1,mu1,Pt,srcvals,ns,E,H,vf_minus,1)
!	read (*,*)

icount=1
do count2=1,n_components
  ep0=contrast_matrix(1,count2)
  mu0=contrast_matrix(2,count2)
  ep1=contrast_matrix(3,count2)
  mu1=contrast_matrix(4,count2)
	do count1=1,npts_vect(count2)	
		call orthonormalize(srcvals(4:6,icount),srcvals(10:12,icount),ru,rv)
    if (exposed_surfaces(count2)) then
      RHS(icount)=-DOT_PRODUCT(rv,E0(:,icount))
      RHS(ns+icount)=DOT_PRODUCT(ru,E0(:,icount))	
      RHS(2*ns+icount)=-DOT_PRODUCT(rv,H0(:,icount))
      RHS(3*ns+icount)=DOT_PRODUCT(ru,H0(:,icount))
    else
      RHS(icount)=0.0d0
      RHS(ns+icount)=0.0d0	
      RHS(2*ns+icount)=0.0d0
      RHS(3*ns+icount)=0.0d0
    endif
    RHS(icount)=RHS(icount)/(mu0+mu1)
    RHS(icount+ns)=RHS(icount+ns)/(mu0+mu1)
    RHS(icount+2*ns)=RHS(icount+2*ns)/(ep0+ep1)
    RHS(icount+3*ns)=RHS(icount+3*ns)/(ep0+ep1)
    icount=icount+1
	enddo
enddo

return
end subroutine get_rhs_em_muller_trans_PW




subroutine test_accuracy_em_muller_trans_v2(eps_FMM,sol,omega,Pts_test,&
 &em_constants_test,n_components,ns,wts,srcvals,P0,vf,direction, Pol)
implicit none

    !List of calling arguments
    complex ( kind = 8 ), intent(in) :: omega
    integer, intent(in) :: ns,n_components
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),direction(2)
	  complex ( kind = 8 ), intent(in) :: sol(4*ns),vf(3),Pol(2)
    complex ( kind = 8 ), intent(in) :: em_constants_test(2,n_components+1)
    real ( kind = 8 ), intent(in) :: Pts_test(3,n_components+1)
	
    !List of local variables
    complex ( kind = 8 ) ep0,mu0,ep1,mu1,ep,mu
	  complex ( kind = 8 ) ima, R1, R2
    complex ( kind = 8 ) E(3),H(3)
    complex ( kind = 8 ) E0(3),H0(3)

	  real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
	  real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),error_E,rel_err_E,error_H,rel_err_H
	  real ( kind = 8 ) pi
	  integer count1


		ima=(0.0d0,1.0d0)
		pi=3.1415926535897932384626433832795028841971d0

    do count1=1,n_components+1
      ep = em_constants_test(1,count1)
      mu = em_constants_test(2,count1)
      call em_muller_trans_FMM_targ(eps_FMM,omega,ep,mu,ns,srcvals,1,Pts_test(:,count1),wts,&
      &sol(1:ns),sol(ns+1:2*ns),sol(2*ns+1:3*ns),sol(3*ns+1:4*ns),E,H)
 
      if (count1.eq.1) then      
        call fieldsEDomega(omega,ep,mu,P0,Pts_test(:,count1),1,E0,H0,vf,0)
        call fieldsMDomega(omega,ep,mu,P0,Pts_test(:,count1),1,E0,H0,vf,1)
      else
        call fieldsPWomega(omega,ep,mu,Pts_test(:,count1),1,E0,H0,direction,Pol)
      endif


      write (*,*) 'Errors in region:',count1
      write (*,*) 'targ point: ',Pts_test(:,count1)
      write (*,*) 'E0: ',E0
      write (*,*) 'E: ',E
      write (*,*) 'H0: ',H0
      write (*,*) 'H: ',H
      
      error_E=sqrt(abs(E0(1)-E(1))**2+abs(E0(2)-E(2))**2+abs(E0(3)-E(3))**2)
  !		write (*,*) 'Error E: ', error_E
      write (*,*) 'Relative Error in E: ', error_E/sqrt(abs(E0(1))**2+abs(E0(2))**2+abs(E0(3))**2)
      
      error_H=sqrt(abs(H0(1)-H(1))**2+abs(H0(2)-H(2))**2+abs(H0(3)-H(3))**2)
  !		write (*,*) 'Error H: ', error_H
      write (*,*) 'Relative Error in H: ', error_H/sqrt(abs(H0(1))**2+abs(H0(2))**2+abs(H0(3))**2)

    enddo

return
end subroutine test_accuracy_em_muller_trans_v2


subroutine em_muller_trans_FMM_targ(eps,omega,ep,mu,ns,srcvals,nt,targ,wts,a_u,a_v,b_u,b_v,E,H)
implicit none

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: omega,ep,mu
    integer, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),b_u(ns),b_v(ns)
    complex ( kind = 8 ), intent(out) :: E(3,nt),H(3,nt)

    !List of local variables
	real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:),v_vect(:,:),source(:,:)
    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:),lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: curlE(:,:),divE(:)
	complex ( kind = 8 ) ima,zk

    integer count1,count2
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

	ima=(0.0d0,1.0d0)
	zk=omega*sqrt(ep*mu)

    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
	allocate(rho(ns))
    allocate(curlE(3,nt))
	allocate(divE(nt))
	
	allocate(n_vect(3,ns))
	allocate(u_vect(3,ns))
	allocate(v_vect(3,ns))
	allocate(source(3,ns))

	do count1=1,ns
		n_vect(:,count1)=srcvals(10:12,count1)
		source(:,count1)=srcvals(1:3,count1)
	enddo
	call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect,v_vect,ns)

    do count1=1,ns
		a_vect(1,count1)=(b_u(count1)*u_vect(1,count1)+b_v(count1)*v_vect(1,count1))/(-ima*omega)
		a_vect(2,count1)=(b_u(count1)*u_vect(2,count1)+b_v(count1)*v_vect(2,count1))/(-ima*omega)
		a_vect(3,count1)=(b_u(count1)*u_vect(3,count1)+b_v(count1)*v_vect(3,count1))/(-ima*omega)
				
        b_vect(1,count1)=(a_u(count1)*u_vect(1,count1)+a_v(count1)*v_vect(1,count1))*mu
        b_vect(2,count1)=(a_u(count1)*u_vect(2,count1)+a_v(count1)*v_vect(2,count1))*mu
        b_vect(3,count1)=(a_u(count1)*u_vect(3,count1)+a_v(count1)*v_vect(3,count1))*mu
	enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
 &b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,curlE,ifcurlE,E,ifdivE,divE,nt,targ)	!watch out!! the E and curlE are fliped

    do count1=1,ns
		a_vect(1,count1)=(a_u(count1)*u_vect(1,count1)+a_v(count1)*v_vect(1,count1))/(ima*omega)
		a_vect(2,count1)=(a_u(count1)*u_vect(2,count1)+a_v(count1)*v_vect(2,count1))/(ima*omega)
		a_vect(3,count1)=(a_u(count1)*u_vect(3,count1)+a_v(count1)*v_vect(3,count1))/(ima*omega)
				
        b_vect(1,count1)=(b_u(count1)*u_vect(1,count1)+b_v(count1)*v_vect(1,count1))*ep
        b_vect(2,count1)=(b_u(count1)*u_vect(2,count1)+b_v(count1)*v_vect(2,count1))*ep
        b_vect(3,count1)=(b_u(count1)*u_vect(3,count1)+b_v(count1)*v_vect(3,count1))*ep
	enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
 &b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,curlE,ifcurlE,H,ifdivE,divE,nt,targ)	!watch out!! the E and curlE are fliped


	deallocate(a_vect)
	deallocate(b_vect)
	deallocate(lambda)
	deallocate(curlE)
	deallocate(rho)
	
	deallocate(u_vect)
	deallocate(v_vect)
	deallocate(n_vect)
	deallocate(source)
	deallocate(divE)

return
end subroutine em_muller_trans_FMM_targ



subroutine em_dfie_trans(srcinfo, ndt,targinfo,ndd, dpars,ndz,zpars,&
 &ndi,ipars,E_val)
implicit none

! this function provides the near field kernel that will use
! zgetnearquad_ggq_guru 
! through getnearquad_DFIE

    !List of calling arguments
	integer, intent(in) :: ndt,ndd,ndz,ndi
	real ( kind = 8 ), intent(in) :: srcinfo(12)
	real ( kind = 8 ), intent(in) :: targinfo(ndt)
	integer, intent(in) :: ipars(ndi)
	real ( kind = 8 ), intent(in) :: dpars(ndd)
	complex ( kind = 8 ), intent(in) :: zpars(ndz)
	complex ( kind = 8 ), intent(out) :: E_val
	
	!List of local variables
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3),aux_real
  real ( kind = 8 ) c_aux,c_aux_A,c_aux_B,c_aux_C,c_aux_D
  real ( kind = 8 ) xprod_aux3(3),xprod_aux4(3)	
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)
	complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima	
	complex ( kind = 8 ) my_exp_0,my_exp_1,omega,ep0	
	complex ( kind = 8 ) mu0,ep1,mu1,zk0,zk1,z_aux_0,z_aux_1
	real ( kind = 8 ) pi
	integer count1,count2
	
	pi=3.1415926535897932384626433832795028841971d0
	ima=(0.0d0,1.0d0)
	omega=zpars(1)
	ep0=zpars(2)
	mu0=zpars(3)
	ep1=zpars(4)
	mu1=zpars(5)
	
	sour(1)=srcinfo(1)
	sour(2)=srcinfo(2)
	sour(3)=srcinfo(3)
	
	n_s(1)=srcinfo(10)
	n_s(2)=srcinfo(11)
	n_s(3)=srcinfo(12)	

	targ(1)=targinfo(1)
	targ(2)=targinfo(2)
	targ(3)=targinfo(3)

	n_t(1)=targinfo(10)
	n_t(2)=targinfo(11)
	n_t(3)=targinfo(12)

	dr(1)=targ(1)-sour(1)
	dr(2)=targ(2)-sour(2)
	dr(3)=targ(3)-sour(3)
	
	r=sqrt((dr(1))**2+(dr(2))**2+(dr(3))**2)
	
	zk0=omega*sqrt(ep0*mu0)
	zk1=omega*sqrt(ep1*mu1)

!R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
!R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
	
	
!R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*exp(ima*zk0*r)/&
!&(4.0d0*pi)
!R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*exp(ima*zk1*r)/&
!&(4.0d0*pi)
	 
!my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
!my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	
	call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
	call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)

    if (ipars(1).eq.1) then
      if (ipars(2).eq.1) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
	  	call my_cross_v2(dr,ru_s,xprod_aux1)		
		c_aux=-DOT_PRODUCT(xprod_aux1,rv_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(1,1)=mu0*nxcurlSk0a(1,1)-mu1*nxcurlSk1a(1,1)      OK
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
    	call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=-DOT_PRODUCT(xprod_aux2,rv_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(1,2)=mu0*nxcurlSk0a(1,2)-mu1*nxcurlSk1a(1,2)       OK
	  elseif (ipars(2).eq.3) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	    call my_cross_v2(n_t,n_s,xprod_aux1)
	    c_aux=DOT_PRODUCT(ru_t,xprod_aux1)
	    E_val=c_aux/r*(-mu0*my_exp_0+mu1*my_exp_1)	  
!E_mat(1,3)=-mu0*nxSk0nrho(1,1)+mu1*nxSk1nrho(1,1)
	  elseif (ipars(2).eq.4) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(rv_t,ru_s)
        E_val=-c_aux/r*(ep0*mu0*my_exp_0-ep1*mu1*my_exp_1)
!E_mat(1,4)=ep0*mu0*nxSk0b(1,1)-ep1*mu1*nxSk1b(1,1)
	  elseif (ipars(2).eq.5) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(rv_t,rv_s)
        e_val=-c_aux/r*(ep0*mu0*my_exp_0-ep1*mu1*my_exp_1)
!E_mat(1,5)=ep0*mu0*nxSk0b(1,2)-ep1*mu1*nxSk1b(1,2)
	  else
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(rv_t,dr)
        E_val=-c_aux*(R1_0-R1_1)
!E_mat(1,6)=nxgradSk0lambda(1,1)-nxgradSk1lambda(1,1)
	  endif
	elseif (ipars(1).eq.2) then
      if (ipars(2).eq.1) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,ru_s,xprod_aux1)
	  	c_aux=DOT_PRODUCT(xprod_aux1,ru_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(2,1)=mu0*nxcurlSk0a(2,1)-mu1*nxcurlSk1a(2,1)     OK
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=DOT_PRODUCT(xprod_aux2,ru_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)	  
!E_mat(2,2)=mu0*nxcurlSk0a(2,2)-mu1*nxcurlSk1a(2,2)     OK
	  elseif (ipars(2).eq.3) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	  	call my_cross_v2(n_t,n_s,xprod_aux1)
	    c_aux=DOT_PRODUCT(rv_t,xprod_aux1)
	    E_val=c_aux/r*(-mu0*my_exp_0+mu1*my_exp_1)
!E_mat(2,3)=-mu0*nxSk0nrho(2,1)+mu1*nxSk1nrho(2,1)
	  elseif (ipars(2).eq.4) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(ru_t,ru_s)
	    E_val=c_aux/r*(ep0*mu0*my_exp_0-ep1*mu1*my_exp_1)
!E_mat(2,4)=ep0*mu0*nxSk0b(2,1)-ep1*mu1*nxSk1b(2,1)
	  elseif (ipars(2).eq.5) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	  	c_aux=DOT_PRODUCT(ru_t,rv_s)
		E_val=c_aux/r*(ep0*mu0*my_exp_0-ep1*mu1*my_exp_1)
!E_mat(2,5)=ep0*mu0*nxSk0b(2,2)-ep1*mu1*nxSk1b(2,2)
	  else
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(ru_t,dr)
        E_val=c_aux*(R1_0-R1_1)
!E_mat(2,6)=nxgradSk0lambda(2,1)-nxgradSk1lambda(2,1)
	  endif
	elseif (ipars(1).eq.3) then
      if (ipars(2).eq.1) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	    c_aux_A=DOT_PRODUCT(rv_t,ru_s)
	    c_aux_B=DOT_PRODUCT(rv_t,ru_s)
	    c_aux_C=DOT_PRODUCT(rv_t,dr)
	    c_aux_D=DOT_PRODUCT(ru_s,-dr)
	    z_aux_0=zk0**2*(-c_aux_A*my_exp_0/r)+(-c_aux_B*R1_0)-&
		 &(-c_aux_C*c_aux_D*R2_0)
	    z_aux_1=zk1**2*(-c_aux_A*my_exp_1/r)+(-c_aux_B*R1_1)-&
		 &(-c_aux_C*c_aux_D*R2_1)
	    E_val=z_aux_0-z_aux_1
!E_mat(3,1)=(nxcurlcurlSk0a(1,1)-nxcurlcurlSk1a(1,1))
	  elseif (ipars(2).eq.2) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	  	c_aux_A=DOT_PRODUCT(rv_t,rv_s)
	    c_aux_B=DOT_PRODUCT(rv_t,rv_s)
	    c_aux_C=DOT_PRODUCT(rv_t,dr)
	    c_aux_D=DOT_PRODUCT(rv_s,-dr)
	    z_aux_0=zk0**2*(-c_aux_A*my_exp_0/r)+(-c_aux_B*R1_0)-&
		 &(-c_aux_C*c_aux_D*R2_0)
        z_aux_1=zk1**2*(-c_aux_A*my_exp_1/r)+(-c_aux_B*R1_1)-&
		 &(-c_aux_C*c_aux_D*R2_1)
        E_val=z_aux_0-z_aux_1
!E_mat(3,2)=(nxcurlcurlSk0a(1,2)-nxcurlcurlSk1a(1,2))
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        call my_cross_v2(dr, n_s, xprod_aux1)
        c_aux=DOT_PRODUCT(rv_t,xprod_aux1)
        E_val=-c_aux*(-R1_0+R1_1)	
!E_mat(3,3)=-nxcurlSk0nrho(1,1)+nxcurlSk1nrho(1,1)
	  elseif (ipars(2).eq.4) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,ru_s,xprod_aux1)
		c_aux=-DOT_PRODUCT(xprod_aux1,rv_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(3,4)=ep0*nxcurlSk0a(1,1)-ep1*nxcurlSk1a(1,1)
	  elseif (ipars(2).eq.5) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=-DOT_PRODUCT(xprod_aux2,rv_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(3,5)=ep0*nxcurlSk0a(1,2)-ep1*nxcurlSk1a(1,2)
	  else
	    E_val=0.0d0
	  endif
	elseif (ipars(1).eq.4) then
      if (ipars(2).eq.1) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
        c_aux_A=DOT_PRODUCT(ru_t,ru_s)
        c_aux_B=DOT_PRODUCT(ru_t,ru_s)
        c_aux_C=DOT_PRODUCT(ru_t,dr)
        c_aux_D=DOT_PRODUCT(ru_s,-dr)
        z_aux_0=zk0**2*(c_aux_A*my_exp_0/r)+(+c_aux_B*R1_0)-&
		 &(+c_aux_C*c_aux_D*R2_0)
        z_aux_1=zk1**2*(c_aux_A*my_exp_1/r)+(+c_aux_B*R1_1)-&
		 &(+c_aux_C*c_aux_D*R2_1)
        E_val=z_aux_0-z_aux_1
!E_mat(4,1)=(nxcurlcurlSk0a(2,1)-nxcurlcurlSk1a(2,1))
	  elseif (ipars(2).eq.2) then	
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	  	c_aux_A=DOT_PRODUCT(ru_t,rv_s)
	    c_aux_B=DOT_PRODUCT(ru_t,rv_s)
	    c_aux_C=DOT_PRODUCT(ru_t,dr)
	    c_aux_D=DOT_PRODUCT(rv_s,-dr)
	    z_aux_0=zk0**2*(c_aux_A*my_exp_0/r)+(+c_aux_B*R1_0)-&
		 &(+c_aux_C*c_aux_D*R2_0)
	    z_aux_1=zk1**2*(c_aux_A*my_exp_1/r)+(+c_aux_B*R1_1)-&
		 &(+c_aux_C*c_aux_D*R2_1)
        E_val=z_aux_0-z_aux_1
!E_mat(4,2)=(nxcurlcurlSk0a(2,2)-nxcurlcurlSk1a(2,2))
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        call my_cross_v2(dr, n_s, xprod_aux1)
        c_aux=DOT_PRODUCT(ru_t,xprod_aux1)
        E_val=c_aux*(-R1_0+R1_1)
!E_mat(4,3)=-nxcurlSk0nrho(2,1)+nxcurlSk1nrho(2,1)
	  elseif (ipars(2).eq.4) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,ru_s,xprod_aux1)
	  	c_aux=DOT_PRODUCT(xprod_aux1,ru_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(4,4)=ep0*nxcurlSk0a(2,1)-ep1*nxcurlSk1a(2,1)
	  elseif (ipars(2).eq.5) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=DOT_PRODUCT(xprod_aux2,ru_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)	  
!E_mat(4,5)=ep0*nxcurlSk0a(2,2)-ep1*nxcurlSk1a(2,2)
	  else
	    E_val=0.0d0
	  endif
	elseif (ipars(1).eq.5) then
      if (ipars(2).eq.1) then
	    E_val=0.0d0
	  elseif (ipars(2).eq.2) then
        E_val=0.0d0
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(n_s,dr)
        E_val=-c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(5,3)=mu0*Dk0rho(1,1)-mu1*Dk1rho(1,1)
      elseif (ipars(2).eq.4) then	  
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(dr,ru_s)
        E_val=c_aux*(ep0*mu0*R1_0-ep1*mu1*R1_1)
!E_mat(5,4)=ep0*mu0*divSk0b(1,1)-ep1*mu1*divSk1b(1,1)
	  elseif (ipars(2).eq.5) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(dr,rv_s)
        E_val=c_aux*(ep0*mu0*R1_0-ep1*mu1*R1_1)
!E_mat(5,5)=ep0*mu0*divSk0b(1,2)-ep1*mu1*divSk1b(1,2)
	  else
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	    E_val=(-zk0**2*my_exp_0+zk1**2*my_exp_1)/r
!E_mat(5,6)=-zk0**2*Sk0lambda(1,1)+zk1**2*Sk1lambda(1,1)
	  endif
	else
      if (ipars(2).eq.1) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        call my_cross_v2(dr, ru_s, xprod_aux1)
        c_aux=DOT_PRODUCT(n_t,xprod_aux1)
        E_val=c_aux*(ep0*mu0*R1_0-ep1*mu1*R1_1)	  
!E_mat(6,1)=ep0*mu0*ncurlSk0a(1,1)-ep1*mu1*ncurlSk1a(1,1)
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        call my_cross_v2(dr, rv_s, xprod_aux2)
        c_aux=DOT_PRODUCT(n_t,xprod_aux2)
        E_val=c_aux*(ep0*mu0*R1_0-ep1*mu1*R1_1)
!E_mat(6,2)=ep0*mu0*ncurlSk0a(1,2)-ep1*mu1*ncurlSk1a(1,2)
	  elseif (ipars(2).eq.3) then
        c_aux=DOT_PRODUCT(n_t,n_s)
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        E_val=c_aux/r*(-ep0*mu0*my_exp_0+ep1*mu1*my_exp_1)
!E_mat(6,3)=-ep0*mu0*nSk0nrho(1,1)+ep1*mu1*nSk1nrho(1,1)
	  elseif (ipars(2).eq.4) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(n_t,ru_s)
        E_val=c_aux/r*(mu0*ep0**2*my_exp_0-mu1*ep1**2*my_exp_1)
!E_mat(6,4)=mu0*ep0**2*nSk0b(1,1)-mu1*ep1**2*nSk1b(1,1)
	  elseif (ipars(2).eq.5) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(n_t,rv_s)
        E_val=c_aux/r*(mu0*ep0**2*my_exp_0-mu1*ep1**2*my_exp_1)
!E_mat(6,5)=mu0*ep0**2*nSk0b(1,2)-mu1*ep1**2*nSk1b(1,2)
	  else
        c_aux=DOT_PRODUCT(n_t,dr)
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(6,6)=ep0*ngradSk0lambda(1,1)-ep1*ngradSk1lambda(1,1)
	  endif
	endif

return
end subroutine em_dfie_trans






subroutine em_muller_trans_v2(srcinfo, ndt,targinfo,ndd, dpars,ndz,zpars,&
 &ndi,ipars,E_val)
implicit none

! this function provides the near field kernel that will use
! zgetnearquad_ggq_guru 
! through getnearquad_DFIE

    !List of calling arguments
	integer, intent(in) :: ndt,ndd,ndz,ndi
	real ( kind = 8 ), intent(in) :: srcinfo(12)
	real ( kind = 8 ), intent(in) :: targinfo(ndt)
	integer, intent(in) :: ipars(ndi)
	real ( kind = 8 ), intent(in) :: dpars(ndd)
	complex ( kind = 8 ), intent(in) :: zpars(ndz)
	complex ( kind = 8 ), intent(out) :: E_val
	
	!List of local variables
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3),aux_real
  real ( kind = 8 ) c_aux,c_aux_A,c_aux_B,c_aux_C,c_aux_D
  real ( kind = 8 ) xprod_aux3(3),xprod_aux4(3)	
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)
	complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima	
	complex ( kind = 8 ) my_exp_0,my_exp_1,omega,ep0	
	complex ( kind = 8 ) mu0,ep1,mu1,zk0,zk1,z_aux_0,z_aux_1
	real ( kind = 8 ) pi
	integer count1,count2
	
	pi=3.1415926535897932384626433832795028841971d0
	ima=(0.0d0,1.0d0)
	omega=zpars(1)
	ep0=zpars(2)
	mu0=zpars(3)
	ep1=zpars(4)
	mu1=zpars(5)
	
	sour(1)=srcinfo(1)
	sour(2)=srcinfo(2)
	sour(3)=srcinfo(3)
	
	n_s(1)=srcinfo(10)
	n_s(2)=srcinfo(11)
	n_s(3)=srcinfo(12)	

	targ(1)=targinfo(1)
	targ(2)=targinfo(2)
	targ(3)=targinfo(3)

	n_t(1)=targinfo(10)
	n_t(2)=targinfo(11)
	n_t(3)=targinfo(12)

	dr(1)=targ(1)-sour(1)
	dr(2)=targ(2)-sour(2)
	dr(3)=targ(3)-sour(3)
	
	r=sqrt((dr(1))**2+(dr(2))**2+(dr(3))**2)
	
  ep0=targinfo(13)+ima*targinfo(14)
  mu0=targinfo(15)+ima*targinfo(16)
	ep1=targinfo(17)+ima*targinfo(18)
  mu1=targinfo(19)+ima*targinfo(20)


	zk0=omega*sqrt(ep0*mu0)
	zk1=omega*sqrt(ep1*mu1)

!R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
!R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
	
	
!R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*exp(ima*zk0*r)/&
!&(4.0d0*pi)
!R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*exp(ima*zk1*r)/&
!&(4.0d0*pi)
	 
!my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
!my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	
	call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
	call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)

    if (ipars(1).eq.1) then
      if (ipars(2).eq.1) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
	  	call my_cross_v2(dr,ru_s,xprod_aux1)		
		c_aux=-DOT_PRODUCT(xprod_aux1,rv_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(1,1)=mu0*nxcurlSk0a(1,1)-mu1*nxcurlSk1a(1,1)      OK
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
    	call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=-DOT_PRODUCT(xprod_aux2,rv_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(1,2)=mu0*nxcurlSk0a(1,2)-mu1*nxcurlSk1a(1,2)       OK
	  elseif (ipars(2).eq.3) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	    c_aux_A=DOT_PRODUCT(rv_t,ru_s)
	    c_aux_B=DOT_PRODUCT(rv_t,ru_s)
	    c_aux_C=DOT_PRODUCT(rv_t,dr)
	    c_aux_D=DOT_PRODUCT(ru_s,-dr)
	    z_aux_0=zk0**2*(-c_aux_A*my_exp_0/r)+(-c_aux_B*R1_0)-&
		 &(-c_aux_C*c_aux_D*R2_0)
	    z_aux_1=zk1**2*(-c_aux_A*my_exp_1/r)+(-c_aux_B*R1_1)-&
		 &(-c_aux_C*c_aux_D*R2_1)
	    E_val=(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(1,3)=(nxcurlcurlSk0a(1,1)-nxcurlcurlSk1a(1,1))/(-ima*omega)     OK
	  else
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	  	c_aux_A=DOT_PRODUCT(rv_t,rv_s)
	    c_aux_B=DOT_PRODUCT(rv_t,rv_s)
	    c_aux_C=DOT_PRODUCT(rv_t,dr)
	    c_aux_D=DOT_PRODUCT(rv_s,-dr)
	    z_aux_0=zk0**2*(-c_aux_A*my_exp_0/r)+(-c_aux_B*R1_0)-&
		 &(-c_aux_C*c_aux_D*R2_0)
        z_aux_1=zk1**2*(-c_aux_A*my_exp_1/r)+(-c_aux_B*R1_1)-&
		 &(-c_aux_C*c_aux_D*R2_1)
        E_val=(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(1,4)=(nxcurlcurlSk0a(1,2)-nxcurlcurlSk1a(1,2))/(-ima*omega)   OK
	  endif
	elseif (ipars(1).eq.2) then
      if (ipars(2).eq.1) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,ru_s,xprod_aux1)
	  	c_aux=DOT_PRODUCT(xprod_aux1,ru_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(2,1)=mu0*nxcurlSk0a(2,1)-mu1*nxcurlSk1a(2,1)     OK
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=DOT_PRODUCT(xprod_aux2,ru_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)	  
!E_mat(2,2)=mu0*nxcurlSk0a(2,2)-mu1*nxcurlSk1a(2,2)     OK
	  elseif (ipars(2).eq.3) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
        c_aux_A=DOT_PRODUCT(ru_t,ru_s)
        c_aux_B=DOT_PRODUCT(ru_t,ru_s)
        c_aux_C=DOT_PRODUCT(ru_t,dr)
        c_aux_D=DOT_PRODUCT(ru_s,-dr)
        z_aux_0=zk0**2*(c_aux_A*my_exp_0/r)+(+c_aux_B*R1_0)-&
		 &(+c_aux_C*c_aux_D*R2_0)
        z_aux_1=zk1**2*(c_aux_A*my_exp_1/r)+(+c_aux_B*R1_1)-&
		 &(+c_aux_C*c_aux_D*R2_1)
        E_val=(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(2,3)=(nxcurlcurlSk0a(2,1)-nxcurlcurlSk1a(2,1))/(-ima*omega)   OK
	  else
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	  	c_aux_A=DOT_PRODUCT(ru_t,rv_s)
	    c_aux_B=DOT_PRODUCT(ru_t,rv_s)
	    c_aux_C=DOT_PRODUCT(ru_t,dr)
	    c_aux_D=DOT_PRODUCT(rv_s,-dr)
	    z_aux_0=zk0**2*(c_aux_A*my_exp_0/r)+(+c_aux_B*R1_0)-&
		 &(+c_aux_C*c_aux_D*R2_0)
	    z_aux_1=zk1**2*(c_aux_A*my_exp_1/r)+(+c_aux_B*R1_1)-&
		 &(+c_aux_C*c_aux_D*R2_1)
        E_val=(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(2,4)=(nxcurlcurlSk0a(2,2)-nxcurlcurlSk1a(2,2))/(-ima*omega)
	  endif
	elseif (ipars(1).eq.3) then
      if (ipars(2).eq.1) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	    c_aux_A=DOT_PRODUCT(rv_t,ru_s)
	    c_aux_B=DOT_PRODUCT(rv_t,ru_s)
	    c_aux_C=DOT_PRODUCT(rv_t,dr)
	    c_aux_D=DOT_PRODUCT(ru_s,-dr)
	    z_aux_0=zk0**2*(-c_aux_A*my_exp_0/r)+(-c_aux_B*R1_0)-&
		 &(-c_aux_C*c_aux_D*R2_0)
	    z_aux_1=zk1**2*(-c_aux_A*my_exp_1/r)+(-c_aux_B*R1_1)-&
		 &(-c_aux_C*c_aux_D*R2_1)
	    E_val=-(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(3,1)=-(nxcurlcurlSk0a(1,1)-nxcurlcurlSk1a(1,1))/(-ima*omega)     OK
	  elseif (ipars(2).eq.2) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	  	c_aux_A=DOT_PRODUCT(rv_t,rv_s)
	    c_aux_B=DOT_PRODUCT(rv_t,rv_s)
	    c_aux_C=DOT_PRODUCT(rv_t,dr)
	    c_aux_D=DOT_PRODUCT(rv_s,-dr)
	    z_aux_0=zk0**2*(-c_aux_A*my_exp_0/r)+(-c_aux_B*R1_0)-&
		 &(-c_aux_C*c_aux_D*R2_0)
        z_aux_1=zk1**2*(-c_aux_A*my_exp_1/r)+(-c_aux_B*R1_1)-&
		 &(-c_aux_C*c_aux_D*R2_1)
        E_val=-(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(3,2)=-(nxcurlcurlSk0a(1,2)-nxcurlcurlSk1a(1,2))/(-ima*omega)   OK
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		    call my_cross_v2(dr,ru_s,xprod_aux1)
		    c_aux=-DOT_PRODUCT(xprod_aux1,rv_t)
		    E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(3,3)=ep0*nxcurlSk0a(1,1)-ep1*nxcurlSk1a(1,1)
	  else
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		    call my_cross_v2(dr,rv_s,xprod_aux2)
	      c_aux=-DOT_PRODUCT(xprod_aux2,rv_t)
		    E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(3,4)=ep0*nxcurlSk0a(1,2)-ep1*nxcurlSk1a(1,2)
	  endif
	elseif (ipars(1).eq.4) then
      if (ipars(2).eq.1) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
        c_aux_A=DOT_PRODUCT(ru_t,ru_s)
        c_aux_B=DOT_PRODUCT(ru_t,ru_s)
        c_aux_C=DOT_PRODUCT(ru_t,dr)
        c_aux_D=DOT_PRODUCT(ru_s,-dr)
        z_aux_0=zk0**2*(c_aux_A*my_exp_0/r)+(+c_aux_B*R1_0)-&
		 &(+c_aux_C*c_aux_D*R2_0)
        z_aux_1=zk1**2*(c_aux_A*my_exp_1/r)+(+c_aux_B*R1_1)-&
		 &(+c_aux_C*c_aux_D*R2_1)
        E_val=-(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(4,1)=-(nxcurlcurlSk0a(2,1)-nxcurlcurlSk1a(2,1))/(-ima*omega)   OK
	  elseif (ipars(2).eq.2) then	
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	  	c_aux_A=DOT_PRODUCT(ru_t,rv_s)
	    c_aux_B=DOT_PRODUCT(ru_t,rv_s)
	    c_aux_C=DOT_PRODUCT(ru_t,dr)
	    c_aux_D=DOT_PRODUCT(rv_s,-dr)
	    z_aux_0=zk0**2*(c_aux_A*my_exp_0/r)+(+c_aux_B*R1_0)-&
		 &(+c_aux_C*c_aux_D*R2_0)
	    z_aux_1=zk1**2*(c_aux_A*my_exp_1/r)+(+c_aux_B*R1_1)-&
		 &(+c_aux_C*c_aux_D*R2_1)
        E_val=-(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(4,2)=-(nxcurlcurlSk0a(2,2)-nxcurlcurlSk1a(2,2))/(-ima*omega)
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,ru_s,xprod_aux1)
	  	c_aux=DOT_PRODUCT(xprod_aux1,ru_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(4,3)=ep0*nxcurlSk0a(2,1)-ep1*nxcurlSk1a(2,1)
	  else
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,rv_s,xprod_aux2)
	  c_aux=DOT_PRODUCT(xprod_aux2,ru_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)	  
!E_mat(4,4)=ep0*nxcurlSk0a(2,2)-ep1*nxcurlSk1a(2,2)
	  endif
	endif
return
end subroutine em_muller_trans_v2


subroutine build_extended_targ(n_components,srcvals_extended,srcvals,&
 &npts_vect,contrast_matrix,npts)
implicit none

  !List of calling arguments
  integer, intent(in) :: npts,n_components
  integer, intent(in) :: npts_vect(n_components)
	real ( kind = 8 ), intent(in) :: srcvals(12,npts)
  complex ( kind = 8 ), intent(in) :: contrast_matrix(4,n_components)
	real ( kind = 8 ), intent(out) :: srcvals_extended(20,npts)

	!List of local variables
  integer count1,count2,icount

    icount=1
    do count2=1,n_components
      do count1=1,npts_vect(count2)
        srcvals_extended(1:12,icount)=srcvals(1:12,icount)
        srcvals_extended(13,icount)=real(contrast_matrix(1,count2))
        srcvals_extended(14,icount)=aimag(contrast_matrix(1,count2))
        srcvals_extended(15,icount)=real(contrast_matrix(2,count2))
        srcvals_extended(16,icount)=aimag(contrast_matrix(2,count2))
        srcvals_extended(17,icount)=real(contrast_matrix(3,count2))
        srcvals_extended(18,icount)=aimag(contrast_matrix(3,count2))
        srcvals_extended(19,icount)=real(contrast_matrix(4,count2))
        srcvals_extended(20,icount)=aimag(contrast_matrix(4,count2))
        icount=icount+1
      enddo
    enddo

  return
  end subroutine build_extended_targ

subroutine fieldsPWzomega(omega,ep,mu,xyz,n,E,H)
! a plane wave; that is: E=exp(ima*zk*z) ux
! H=exp(ima*zk*z) uy

	integer, intent(in) :: n
	real ( kind = 8 ), intent(in) :: xyz(12,n)
	complex ( kind = 8 ), intent(in) :: omega,ep,mu
	complex *16 E(3,n),H(3,n)
	complex *16 ima,zk,eta
	integer i
	data ima/(0.0d0,1.0d0)/

	  zk=omega*sqrt(ep*mu)
	  eta=sqrt(mu/ep)
      do i=1,n
        H(1,i)=0
        H(2,i)=exp(ima*zk*xyz(3,i))/eta
        H(3,i)=0
        E(1,i)=exp(ima*zk*xyz(3,i))
        E(2,i)=0
        E(3,i)=0
      enddo

return
end subroutine fieldsPWzomega


subroutine fieldsPWomega(omega,ep,mu,xyz,n,E,H,direction,Pol)
! a plane wave; that is: E=exp(ima*zk*z) ux
! H=exp(ima*zk*z) uy

	integer, intent(in) :: n
	real ( kind = 8 ), intent(in) :: xyz(3,n),direction(2)
	complex ( kind = 8 ), intent(in) :: omega,ep,mu,Pol(2)
	complex *16 E(3,n),H(3,n)
	complex *16 ima,zk,eta,k_vect(3),E_phi,E_theta,H_phi,H_theta,eprop
  real ( kind = 8 ) phi,theta
	integer i
	data ima/(0.0d0,1.0d0)/

    phi=direction(1)
    theta=direction(2)

	  zk=omega*sqrt(ep*mu)
    k_vect(1)=zk*cos(phi)*sin(theta)
    k_vect(2)=zk*sin(phi)*sin(theta)
    k_vect(3)=zk*cos(theta)
    E_phi=Pol(1)
    E_theta=Pol(2)
	  eta=sqrt(mu/ep)
    do i=1,n
      eprop=exp(ima*(k_vect(1)*xyz(1,i)+k_vect(2)*xyz(2,i)+k_vect(3)*xyz(3,i)))
      E(1,i)=(E_theta*cos(theta)*cos(phi)-E_phi*sin(phi))*eprop;
      E(2,i)=(E_theta*cos(theta)*sin(phi)+E_phi*cos(phi))*eprop;
      E(3,i)=(-E_theta*sin(theta))*eprop;
      H_phi=E_theta/eta
      H_theta=-E_phi/eta       
      H(1,i)=(H_theta*cos(theta)*cos(phi)-H_phi*sin(phi))*eprop;
      H(2,i)=(H_theta*cos(theta)*sin(phi)+H_phi*cos(phi))*eprop;
      H(3,i)=(-H_theta*sin(theta))*eprop;
    enddo

return
end subroutine fieldsPWomega


subroutine open_gov3_geometry_v2(filename,npatches,norders,ixyzs, &
  iptype,npoints,srcvals,srccoefs,wts,dP)
  implicit none

  !List of calling arguments
  character (len=*), intent(in) :: filename
	integer ( kind = 4 ), intent(in) :: npoints
	integer ( kind = 4 ), intent(in) :: npatches
  integer ( kind = 4), intent(out) :: norders(npatches)
  integer ( kind = 4), intent(out) :: ixyzs(npatches+1)
  integer ( kind = 4), intent(out) :: iptype(npatches)
	real ( kind = 8 ), intent(out) :: srcvals(12,npoints)
	real ( kind = 8 ), intent(out) :: srccoefs(9,npoints)
	real ( kind = 8 ), intent(out) :: wts(npoints)
	real ( kind = 8 ), intent(in) :: dP(4)

  !List of local variables
  integer ( kind = 4 ) umio,count1,count2,flag,aux,npols,icount
  integer ( kind = 4) norder
  integer :: ierror
	real ( kind = 8 ) aux_real,aux_vect(3)
	real ( kind = 8 ), allocatable :: h_points(:),h_coefs(:)
	real ( kind = 8 ), allocatable :: uv(:,:),umatr(:,:),vmatr(:,:),w(:)
  integer i

    open(UNIT=18, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)

        read(18,*) norder
        read(18,*) count1

        do count1=1,npoints
            read(18,*) srcvals(1,count1) !points(1,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(2,count1) !points(2,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(3,count1) !points(3,count1)
        enddo

        do count1=1,npoints
            read(18,*) srcvals(4,count1) !du(1,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(5,count1) !du(2,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(6,count1) !du(3,count1)
        enddo


        do count1=1,npoints
            read(18,*) srcvals(7,count1) !dv(1,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(8,count1) !dv(2,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(9,count1) !dv(3,count1)
        enddo

        do count1=1,npoints
            read(18,*) srcvals(10,count1) !normals(1,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(11,count1) !normals(2,count1)
        enddo
        do count1=1,npoints
            read(18,*) srcvals(12,count1) !normals(3,count1)
        enddo

        close (18)

        do count1=1,npoints
          srcvals(1,count1)=dP(4)*srcvals(1,count1)+dP(1)
          srcvals(2,count1)=dP(4)*srcvals(2,count1)+dP(2)
!          srcvals(3,count1)=dP(4)*(srcvals(3,count1)-1.5d0)+dP(3)
          srcvals(3,count1)=dP(4)*(srcvals(3,count1))+dP(3)
		  
		      srcvals(4,count1)=dP(4)*srcvals(4,count1)
          srcvals(5,count1)=dP(4)*srcvals(5,count1)
          srcvals(6,count1)=dP(4)*srcvals(6,count1)

		      srcvals(7,count1)=dP(4)*srcvals(7,count1)
          srcvals(8,count1)=dP(4)*srcvals(8,count1)
          srcvals(9,count1)=dP(4)*srcvals(9,count1)
        enddo

		    npols = (norder+1)*(norder+2)/2
        do i = 1,npatches
          norders(i) = norder
          iptype(i) = 1
          ixyzs(i) = (i-1)*npols + 1
        enddo
        ixyzs(npatches+1) = npoints+1

		allocate(h_points(npols))
		allocate(h_coefs(npols))
		allocate(uv(2,npols),umatr(npols,npols),vmatr(npols,npols),w(npols))
		call vioreanu_simplex_quad(norder,npols,uv,umatr,vmatr,w)
        
        call get_qwts(npatches,norders,ixyzs,iptype,npoints,srcvals,wts)


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(count1,h_points,h_coefs)
		do count1=1,npatches
			h_points(:)=srcvals(1,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(1,(count1-1)*npols+1:count1*npols)=h_coefs(:)					
			
			h_points(:)=srcvals(2,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(2,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

			h_points(:)=srcvals(3,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(3,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

			h_points(:)=srcvals(4,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(4,(count1-1)*npols+1:count1*npols)=h_coefs(:)					
			
			h_points(:)=srcvals(5,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(5,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

			h_points(:)=srcvals(6,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(6,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

			h_points(:)=srcvals(7,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(7,(count1-1)*npols+1:count1*npols)=h_coefs(:)					
			
			h_points(:)=srcvals(8,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(8,(count1-1)*npols+1:count1*npols)=h_coefs(:)		

			h_points(:)=srcvals(9,(count1-1)*npols+1:count1*npols)
			call dmatvec(npols,npols,umatr,h_points,h_coefs)		
			srccoefs(9,(count1-1)*npols+1:count1*npols)=h_coefs(:)		
		enddo
!$OMP END PARALLEL DO        

return
end subroutine open_gov3_geometry_v2


subroutine evaluate_field_muller(npatches,norders,ixyzs,iptype,npts,srccoefs,&
    &srcvals,wts,targ,ntarg,npatches_vect,n_components,sorted_vector,&
    &contrast_matrix,exposed_surfaces,eps,zpars,sigma,E_far,H_far)
implicit none

!
!  This function calculates the scattered/total fields E,H at a set of target
!  points in space without specifying the ep/mu on each point or the region
!  to which it belongs
!  uses FMM and 3dbie off surface quadratures
!
!  input:
!    npatches - integer
!      number of patches
!
!    norders - integer(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(npatches+1)
!      starting location of data on patch i
!  
!    iptype - integer(npatches)
!      type of patch
!      iptype = 1 -> triangular patch discretized with RV nodes
!
!    npts - integer
!      total number of discretization points on the boundary
!
!    srccoefs - real *8 (9,npts)
!      koornwinder expansion coefficients of xyz, dxyz/du,
!      and dxyz/dv on each patch. 
!      For each point srccoefs(1:3,i) is xyz info
!        srccoefs(4:6,i) is dxyz/du info
!        srccoefs(7:9,i) is dxyz/dv info
!
!    srcvals - real *8 (12,npts)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    wts - real *8 (npts)
!      integration weights for each point on the surface srcvals(1:3,:)
!      it allows integrate smooth functions on the surface
!
!    targ - real *8 (3,ntarg)
!      xyz coordinates of the target points
!
!    ntarg - integer
!      number of target points
!
!    n_components - integer
!      number of different connected components of the geometry
!      that is, number of different interfaces in the transmission
!      problem 
!
!    npatches_vect - integer(n_components)
!      number of patches that forms if the i'th component
!       (the i'th interface)
!      the first interface is formed by patches 1,...,npatcehs_vect(1)
!      the second interface is formed by patches:
!        npatcehs_vect(1)+1,....,npatcehs_vect(1)+npatcehs_vect(2)
!      the third interface is formed by patches:
!        npatcehs_vect(1)+npatcehs_vect(2)+1,...
!        ...,npatcehs_vect(1)+npatches_vect(2)+npatches_vect(3)
!
!    sorted_vector - integer(n_components+1)
!     resulting ordered components (interfaces)
!     sorted_vector(1) is always a surface that does not contain any
!     other surface in the interior
!     sorted_vector(n_components+1) = 0 always (meaning a big sphere
!     at infinity that contains everything)
!
!    contrast_matrix - complex*16(4,n_components)
!      ep and mu at each side of each component (or interface)
!      contrast_matrix(1,i)=ep0  exterior of the surface
!      contrast_matrix(2,i)=mu0  exterior of the surface
!      contrast_matrix(3,i)=ep   interior of the surface
!      contrast_matrix(4,i)=mu   interior of the surface
!
!    exposed_surfaces - logical(n_components)
!     exposed_surfaces(i)=.true. if the i'th interface is exposed to the 
!     ambient space, that is, is not contained
!     in the interior of any other surface.
!     this is useful and used to compute the RHS for an incoming plane wave.
!     exposed_surfaces(i)=.false. otherwise

!
!    eps - real*8
!     accuracty in the calculation of the fmm and near quadratures
!
!
!  output
!    location_targs integer(ntarg)
!      location of each target point, that is, index of the smaller 
!      interface that contains the target point
!


 !List of calling arguments
 integer, intent(in) :: npatches,npts,n_components,ntarg
 integer, intent(in) :: norders(npatches),npatches_vect(n_components),ixyzs(npatches+1),iptype(npatches)
 real ( kind = 8 ), intent(in) :: srcvals(12,npts), srccoefs(9,npts),targ(3,ntarg),wts(npts)
 integer, intent(in) :: sorted_vector(n_components+1)
 complex ( kind = 8 ), intent(in) :: contrast_matrix(4,n_components),zpars(3),sigma(4*npts)
 logical *8 exposed_surfaces(n_components)
 real ( kind = 8 ), intent(in) :: eps
 complex ( kind = 8 ), intent(out) :: E_far(3,ntarg),H_far(3,ntarg)

 !List of local variables
 integer count1,count2,icount,icount2,icount3,n_aux,i1,i2,j1,j2,npatches_aux,npts_aux,x,ntarg_vect(n_components+1)
 integer ndtarg,n_regions
 integer, allocatable :: location_targs(:),permutation_targs(:)
 integer, allocatable :: ipatch_id(:)
 real *8, allocatable :: uvs_targ(:,:)
 complex ( kind = 8 ), allocatable :: E_far_aux(:,:),H_far_aux(:,:)
 real *8, allocatable :: targ_sort(:,:)
 complex ( kind = 8 ) ep0,mu0
 

 allocate(location_targs(ntarg),permutation_targs(ntarg))
 allocate(ipatch_id(ntarg))
 allocate(uvs_targ(2,ntarg))
 allocate(targ_sort(7,ntarg))
 allocate(E_far_aux(3,ntarg),H_far_aux(3,ntarg))

  do count1=1,n_components+1
    ntarg_vect(count1)=0
  enddo

  call find_inclusion_vect(npatches,norders,ixyzs,iptype,npts,srccoefs,&
    &srcvals,wts,targ,ntarg,npatches_vect,n_components,sorted_vector,&
    &location_targs,eps)
  print *, "done computing inclusion vect"
  do count1=1,n_components
    if (exposed_surfaces(count1)) then
       ep0=contrast_matrix(1,count1)
       mu0=contrast_matrix(2,count1)
       exit
    endif
  enddo

 icount=1
 n_regions=0
 do count1=0,n_components
   icount2=0
   do count2=1,ntarg
    if (location_targs(count2).eq.count1) then
      targ_sort(1:3,icount)=targ(:,count2)
      permutation_targs(count2)=icount
      if (count1.eq.0) then
        targ_sort(4,icount)=real(ep0)
        targ_sort(5,icount)=aimag(ep0)
        targ_sort(6,icount)=real(mu0)
        targ_sort(7,icount)=aimag(mu0)
      else
        targ_sort(4,icount)=real(contrast_matrix(3,count1))
        targ_sort(5,icount)=aimag(contrast_matrix(3,count1))
        targ_sort(6,icount)=real(contrast_matrix(4,count1))
        targ_sort(7,icount)=aimag(contrast_matrix(4,count1))
      endif
      icount=icount+1
      icount2=icount2+1
    endif
   enddo
  if (icount2>0) then    
    n_regions=n_regions+1
    ntarg_vect(n_regions)=icount2
  endif
 enddo
 ndtarg=7
 do count1=1,ntarg
   ipatch_id(count1) = -1
   uvs_targ(1,count1) = 0
   uvs_targ(2,count1) = 0
 enddo
! do count1=1,ntarg
!   write (*,*) 'targ_sort: ',count1,targ_sort(:,count1)
! enddo
!  do count1=1,ntarg
!   write (*,*) 'targ: ',count1,targ(:,count1)
! enddo
!write (*,*) 'ntarg_vect: ',ntarg_vect,n_regions
 
 print *, "before lpcomp_em_muller_far_dir"

 call lpcomp_em_muller_far_dir(npatches,norders,ixyzs,&
  &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targ_sort,&
  &ipatch_id,uvs_targ,eps,zpars,sigma,E_far_aux,H_far_aux,n_regions,ntarg_vect)


 do count1=1,ntarg
   E_far(:,count1)=E_far_aux(:,permutation_targs(count1))
   H_far(:,count1)=H_far_aux(:,permutation_targs(count1))   
 enddo
return
end subroutine evaluate_field_muller


      subroutine lpcomp_em_muller_far_dir(npatches,norders,ixyzs,&
       &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       &ipatch_id,uvs_targ,eps,zpars,sigma,E_far,H_far,n_regions,ntarg_vect)
!c
!cf2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
!cf2py intent(in) ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps,zpars
!cf2py intent(in) sigma
!cf2py intent(out) pot
!c
!c
!c------------------------------
!c  This subroutine evaluates the layer potential for the representation 
!c
!c
!c  .. math ::
!c  
!c      u = (\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
!c
!c  Note: For targets on the boundary, this routine only computes
!c  the principal value part, the identity term corresponding to the jump
!c  in the layer potential is not included in the layer potential.
!c
!c
!c  Input arguments:
!c
!c    - npatches: integer
!c        number of patches
!c    - norders: integer(npatches)
!c        order of discretization on each patch 
!c    - ixyzs: integer(npatches+1)
!c        ixyzs(i) denotes the starting location in srccoefs,
!c        and srcvals array where information for patch i begins
!c    - iptype: integer(npatches)
!c        type of patch
!c    - npts: integer
!c        total number of discretization points on the boundary
!c    - srccoefs: double precision (9,npts)
!c        koornwinder expansion coefficients of x, $\partial_{u} x$,
!c        and $\partial_{v} x$. 
!c    - srcvals: double precision (12,npts)
!c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
!c        discretization nodes
!c    - ndtarg: integer
!c        leading dimension of target array
!c    - ntarg: integer
!c        number of targets
!c    - targs: double precision (ndtarg,ntarg)
!c        target information
!c    - ipatch_id: integer(ntarg)
!c        id of patch of target i, id = -1, if target is off-surface
!c    - uvs_targ: double precision (2,ntarg)
!c        local uv coordinates on patch if on surface, otherwise
!c        set to 0 by default
!c    - eps: double precision
!c        precision requested
!c    - zpars: double complex (3)
!c        kernel parameters (Referring to formula above)
!c        zpars(1) = k 
!c        zpars(2) = $\alpha$
!c        zpars(3) = $\beta$
!c     - sigma: double complex(npts)
!c         density for layer potential
!c
!c  Output arguments
!c    - pot: double complex(ntarg)
!c        layer potential evaluated at the target points
!c
!c-----------------------------------
!c
      implicit none
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      complex *16, intent(in) :: zpars(3)
      complex *16, intent(in) :: sigma(4*npts)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      integer, intent(in) :: n_regions
      integer, intent(in) :: ntarg_vect(n_regions+1)

      complex *16, intent(out) :: E_far(3,ntarg),H_far(3,ntarg)


      integer nptso,nnz,nquad


      integer nover,npolso
      integer norder,npols
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      real *8 over4pi
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over
      data over4pi/0.07957747154594767d0/


!c
!c
!c        this might need fixing
!c
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,& 
       &srccoefs,cms,rads)

!C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!C$OMP END PARALLEL DO      

!c
!c    find near quadrature correction interactions
!c
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,ntarg,nnz)

      allocate(row_ptr(ntarg+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,row_ptr,& 
       &col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind,&
       &iquad)

!      ikerorder = -1
!      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 0

      ikerorder = 0

!c
!c    estimate oversampling for far-field, and oversample geometry
!c

      allocate(novers(npatches),ixyzso(npatches+1))

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,&
       &rads,npts,srccoefs,ndtarg,ntarg,targs,ikerorder,zpars(1),&
       &nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts,& 
       &srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,&
       &srcover,wover)

!c
!c   compute near quadrature correction
!c
      nquad = iquad(nnz+1)-1
      write (*,*) 'nquad: ',nquad
      allocate(wnear(12*nquad))
      
!C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,12*nquad
        wnear(i) = 0
      enddo
!C$OMP END PARALLEL DO    


      iquadtype = 1
      call prinf('ntarg=*',ntarg,1)
!      goto 1111

      call getnearquad_muller_far_dir(npatches,norders,&
       &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
       &iquad,rfac0,nquad,wnear)
 1111  continue       
        print *, "done computing near quadrature"

        print *, "entering layer potential evaluator"

!c
!c
!c   compute layer potential
!c
      call lpcomp_em_muller_far_addsub(npatches,norders,ixyzs,&
       &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
       &npts_over,ixyzso,srcover,wover,E_far,H_far,wnear,n_regions,ntarg_vect)


      return
      end
!c


      subroutine getnearquad_muller_far_dir(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the integral DFIE:
!
!  Note: the 4 \pi scaling is NOT!! included as the output of the FMM
!  has been rescaled.
!
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a chunk centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!
!        
!
!  input:
!    npatches - integer
!      number of patches
!
!    norders - integer(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(npatches+1)
!      starting location of data on patch i
!  
!    iptype - integer(npatches)
!      type of patch
!      iptype = 1 -> triangular patch discretized with RV nodes
!
!    npts - integer
!      total number of discretization points on the boundary
!
!    srccoefs - real *8 (9,npts)
!      koornwinder expansion coefficients of xyz, dxyz/du,
!      and dxyz/dv on each patch. 
!      For each point srccoefs(1:3,i) is xyz info
!        srccoefs(4:6,i) is dxyz/du info
!        srccoefs(7:9,i) is dxyz/dv info
!
!    srcvals - real *8 (12,npts)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
! 
!    ndtarg - integer
!      leading dimension of target array
!        
!    ntarg - integer
!      number of targets
!
!    targs - real *8 (ndtarg,ntarg)
!      target information
!
!    ipatch_id - integer(ntarg)
!      id of patch of target i, id = -1, if target is off-surface
!
!         uvs_targ - real *8 (2,ntarg)
!            local uv coordinates on patch if on surface, otherwise
!            set to 0 by default
!          (maybe better to find closest uv on patch using
!            newton)
!            
!          eps - real *8
!             precision requested
!
!          zpars - complex *16(3)
!              kernel parameters
!              zpars(1) = omega 
!              zpars(2) = ep0
!              zpars(3) = mu0
!              zpars(4) = ep1
!              zpars(5) = mu1
!
!           iquadtype - integer
!              quadrature type
!              iquadtype = 1, use ggq for self + adaptive integration
!                 for rest
! 
!
!           nnz - integer
!             number of source patch-> target interactions in the near
!             field
! 
!           row_ptr - integer(ntarg+1)
!              row_ptr(i) is the pointer
!              to col_ind array where list of relevant source patches
!              for target i start
!
!           col_ind - integer (nnz)
!               list of source patches relevant for all targets, sorted
!               by the target number
!
!           iquad - integer(nnz+1)
!               location in wnear array where quadrature for col_ind(i)
!               starts
!
!           rfac0 - integer
!               radius parameter for near field
!
!           nquad - integer
!               number of entries in wnear
!
!        output
!           wnear - complex *16(36*nquad)
!               the desired near field quadrature
!

      implicit none 
      integer npatches,norders(npatches),npts,nquad
      integer ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      real *8 srcvals_extended(20,npts)
      integer ndtarg,ntarg
      integer iquadtype
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(3)
      integer nnz,ipars(2)
      real *8 dpars(1)
      integer row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(12*nquad)
	  
      integer ipatch_id(ntarg)
      real *8 uvs_targ(2,ntarg)

      complex *16 alpha,beta
      integer i,j,ndi,ndd,ndz,count1,count2,icount

      integer ipv
      integer :: t1, t2,clock_rate, clock_max

      procedure (), pointer :: fker
      external  fker_em_muller_far
      external  em_muller_far

     ndz=3
	   ndd=1
	   ndi=2
	   ipv=1

!!      fker => fker_em_muller_trans_v2
!!	  fker => fker_em_muller_far
    fker => em_muller_far

!    call system_clock ( t1, clock_rate, clock_max )

	  icount=0
	  do count1=1,3
	    do count2=1,4
		   ipars(1)=count1
       ipars(2)=count2
       call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		    &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		    &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		    &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		    &wnear(icount*nquad+1:(icount+1)*nquad))
		   icount=icount+1
      enddo
	  enddo


!call system_clock ( t2, clock_rate, clock_max )
!write (*,*) 'time quadratures: ',real(t2-t1)/real(clock_rate)
!stop
    return
    end subroutine getnearquad_muller_far_dir



subroutine fker_em_muller_far(srcinfo, ndt,targinfo,ndd, dpars,ndz,&
 &zpars,ndi,ipars,E_val)
implicit none

!  THIS FUNCTION IS OBSOLETE, NOT USED, subroutine em_dfie_trans is ued
!  instead (much faster), but very hard to read..
!  this function provides the near field kernel that will use 
!  zgetnearquad_ggq_guru 
!  through getnearquad_DFIE

    !List of calling arguments
	integer, intent(in) :: ndt,ndd,ndz,ndi
	real ( kind = 8 ), intent(in) :: srcinfo(12)
	real ( kind = 8 ), intent(in) :: targinfo(ndt)
	integer, intent(in) :: ipars(ndi)
	real ( kind = 8 ), intent(in) :: dpars(ndd)
	complex ( kind = 8 ), intent(in) :: zpars(ndz)
	complex ( kind = 8 ), intent(out) :: E_val
	
	!List of local variables
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(7)
	real ( kind = 8 ) r, dr(3),aux_real
	complex ( kind = 8 ) curlSka(3,2),curlcurlSka(3,2)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1,R2,ima,my_exp,omega	
	complex ( kind = 8 ) mu,ep,zk
	real ( kind = 8 ) pi
	integer count1,count2
	complex ( kind = 8 )  E_mat(6,4)

	
	pi=3.1415926535897932384626433832795028841971d0
	ima=(0.0d0,1.0d0)
	omega=zpars(1)
	
	sour(1)=srcinfo(1)
	sour(2)=srcinfo(2)
	sour(3)=srcinfo(3)
	
	n_s(1)=srcinfo(10)
	n_s(2)=srcinfo(11)
	n_s(3)=srcinfo(12)	

	targ(1)=targinfo(1)
	targ(2)=targinfo(2)
	targ(3)=targinfo(3)

  ep=targinfo(4)+ima*targinfo(5)
  mu=targinfo(6)+ima*targinfo(7)

	dr(1)=targ(1)-sour(1)
	dr(2)=targ(2)-sour(2)
	dr(3)=targ(3)-sour(3)

	
	r=sqrt((dr(1))**2+(dr(2))**2+(dr(3))**2)
	zk=omega*sqrt(ep*mu)

	R1=(ima*zk*r-1.0d0)/r**3*exp(ima*zk*r)/(4.0d0*pi)
	R2=((ima*zk)**2/r**3-3.0d0*ima*zk/r**4+3.0d0/r**5)*exp(ima*zk*r)/&
	 &(4.0d0*pi)
	my_exp=exp(ima*zk*r)/(4.0d0*pi)
	
	call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)

	call get_curlSka(ru_s,rv_s,n_s,dr,R1,my_exp,r,curlSka)		
	call get_curlcurlSka(ru_s,rv_s,n_s,dr,R1,R2,zk,&
	 &my_exp,r,curlcurlSka)
		
		E_mat(1,1)=curlSka(1,1)
    E_mat(2,1)=curlSka(2,1)
		E_mat(3,1)=curlSka(3,1)

		E_mat(1,2)=curlSka(1,2)
    E_mat(2,2)=curlSka(2,2)
		E_mat(3,2)=curlSka(3,2)

		E_mat(1,3)=curlcurlSka(1,1)/(-ima*omega)
		E_mat(2,3)=curlcurlSka(2,1)/(-ima*omega)
		E_mat(3,3)=curlcurlSka(3,1)/(-ima*omega)

		E_mat(1,4)=curlcurlSka(1,2)/(-ima*omega)
		E_mat(2,4)=curlcurlSka(2,2)/(-ima*omega)
		E_mat(3,4)=curlcurlSka(3,2)/(-ima*omega)


!		E_mat(4,3)=ep*curlSka(1,1)
!    E_mat(5,3)=ep*curlSka(2,1)
!		E_mat(6,3)=ep*curlSka(3,1)

!		E_mat(4,4)=ep*curlSka(1,2)
!    E_mat(5,4)=ep*curlSka(2,2)
!		E_mat(6,4)=ep*curlSka(3,2)

!		E_mat(4,1)=curlcurlSka(1,1)/(ima*omega)
!		E_mat(5,1)=curlcurlSka(2,1)/(ima*omega)
!		E_mat(6,1)=curlcurlSka(3,1)/(ima*omega)

!		E_mat(4,2)=curlcurlSka(1,2)/(ima*omega)
!		E_mat(5,2)=curlcurlSka(2,2)/(ima*omega)
!		E_mat(6,2)=curlcurlSka(3,2)/(ima*omega)


	E_val=E_mat(ipars(1),ipars(2))


return
end subroutine fker_em_muller_far



subroutine em_muller_far(srcinfo, ndt,targinfo,ndd, dpars,ndz,&
 &zpars,ndi,ipars,E_val)
implicit none

!  THIS FUNCTION IS OBSOLETE, NOT USED, subroutine em_dfie_trans is ued
!  instead (much faster), but very hard to read..
!  this function provides the near field kernel that will use 
!  zgetnearquad_ggq_guru 
!  through getnearquad_DFIE

    !List of calling arguments
	integer, intent(in) :: ndt,ndd,ndz,ndi
	real ( kind = 8 ), intent(in) :: srcinfo(12)
	real ( kind = 8 ), intent(in) :: targinfo(ndt)
	integer, intent(in) :: ipars(ndi)
	real ( kind = 8 ), intent(in) :: dpars(ndd)
	complex ( kind = 8 ), intent(in) :: zpars(ndz)
	complex ( kind = 8 ), intent(out) :: E_val
	
	!List of local variables
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(7)
	real ( kind = 8 ) r, dr(3),aux_real
	complex ( kind = 8 ) curlSka(3,2),curlcurlSka(3,2)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1,R2,ima,my_exp,omega	
	complex ( kind = 8 ) mu,ep,zk
	real ( kind = 8 ) pi
	integer count1,count2
	complex ( kind = 8 )  E_mat(6,4)

	
	pi=3.1415926535897932384626433832795028841971d0
	ima=(0.0d0,1.0d0)
	omega=zpars(1)
	
	sour(1)=srcinfo(1)
	sour(2)=srcinfo(2)
	sour(3)=srcinfo(3)
	
	n_s(1)=srcinfo(10)
	n_s(2)=srcinfo(11)
	n_s(3)=srcinfo(12)	

	targ(1)=targinfo(1)
	targ(2)=targinfo(2)
	targ(3)=targinfo(3)

  ep=targinfo(4)+ima*targinfo(5)
  mu=targinfo(6)+ima*targinfo(7)

	dr(1)=targ(1)-sour(1)
	dr(2)=targ(2)-sour(2)
	dr(3)=targ(3)-sour(3)

	
	r=sqrt((dr(1))**2+(dr(2))**2+(dr(3))**2)
	zk=omega*sqrt(ep*mu)

	R1=(ima*zk*r-1.0d0)/r**3*exp(ima*zk*r)/(4.0d0*pi)
	R2=((ima*zk)**2/r**3-3.0d0*ima*zk/r**4+3.0d0/r**5)*exp(ima*zk*r)/&
	 &(4.0d0*pi)
	my_exp=exp(ima*zk*r)/(4.0d0*pi)
	
	call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)


    if (ipars(1).eq.1) then
      if (ipars(2).eq.1) then
	      call my_cross_v2(dr,ru_s,xprod_aux1)
        E_val=xprod_aux1(1)*R1
      elseif (ipars(2).eq.2) then
	      call my_cross_v2(dr,rv_s,xprod_aux2)
        E_val=xprod_aux2(1)*R1
      elseif (ipars(2).eq.3) then
        E_val=(zk**2*ru_s(1)*my_exp/r+R1*ru_s(1)-dr(1)*DOT_PRODUCT(ru_s,-dr)*R2)/(-ima*omega)
      elseif (ipars(2).eq.4) then
        E_val=(zk**2*rv_s(1)*my_exp/r+R1*rv_s(1)-dr(1)*DOT_PRODUCT(rv_s,-dr)*R2)/(-ima*omega)
      endif
    elseif (ipars(1).eq.2) then
      if (ipars(2).eq.1) then
	      call my_cross_v2(dr,ru_s,xprod_aux1)
        E_val=xprod_aux1(2)*R1
      elseif (ipars(2).eq.2) then
	      call my_cross_v2(dr,rv_s,xprod_aux2)
        E_val=xprod_aux2(2)*R1
      elseif (ipars(2).eq.3) then
        E_val=(zk**2*ru_s(2)*my_exp/r+R1*ru_s(2)-dr(2)*DOT_PRODUCT(ru_s,-dr)*R2)/(-ima*omega)
      elseif (ipars(2).eq.4) then
        E_val=(zk**2*rv_s(2)*my_exp/r+R1*rv_s(2)-dr(2)*DOT_PRODUCT(rv_s,-dr)*R2)/(-ima*omega)
      endif
    elseif (ipars(1).eq.3) then 
      if (ipars(2).eq.1) then
	      call my_cross_v2(dr,ru_s,xprod_aux1)
        E_val=xprod_aux1(3)*R1
      elseif (ipars(2).eq.2) then
	      call my_cross_v2(dr,rv_s,xprod_aux2)
        E_val=xprod_aux2(3)*R1
      elseif (ipars(2).eq.3) then
        E_val=(zk**2*ru_s(3)*my_exp/r+R1*ru_s(3)-dr(3)*DOT_PRODUCT(ru_s,-dr)*R2)/(-ima*omega)
      elseif (ipars(2).eq.4) then
        E_val=(zk**2*rv_s(3)*my_exp/r+R1*rv_s(3)-dr(3)*DOT_PRODUCT(rv_s,-dr)*R2)/(-ima*omega)
      endif
    endif


!	curlcurlSka(1,1)=zk**2*ru_s(1)*my_exp/r+R1*ru_s(1)-dr(1)*DOT_PRODUCT(ru_s,-dr)*R2
!	curlcurlSka(2,1)=zk**2*ru_s(2)*my_exp/r+R1*ru_s(2)-dr(2)*DOT_PRODUCT(ru_s,-dr)*R2
!	curlcurlSka(3,1)=zk**2*ru_s(3)*my_exp/r+R1*ru_s(3)-dr(3)*DOT_PRODUCT(ru_s,-dr)*R2

!	curlcurlSka(1,2)=zk**2*rv_s(1)*my_exp/r+R1*rv_s(1)-dr(1)*DOT_PRODUCT(rv_s,-dr)*R2
!	curlcurlSka(2,2)=zk**2*rv_s(2)*my_exp/r+R1*rv_s(2)-dr(2)*DOT_PRODUCT(rv_s,-dr)*R2
!	curlcurlSka(3,2)=zk**2*rv_s(3)*my_exp/r+R1*rv_s(3)-dr(3)*DOT_PRODUCT(rv_s,-dr)*R2

return
end subroutine em_muller_far




subroutine get_curlSka(ru_s,rv_s,n_s,dr,R1,my_exp,r,curlSka)
implicit none

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,my_exp
	complex ( kind = 8 ), intent(out) :: curlSka(3,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)

	call my_cross_v2(dr,ru_s,xprod_aux1)
	call my_cross_v2(dr,rv_s,xprod_aux2)
	
	curlSka(1,1)=xprod_aux1(1)*R1
  curlSka(2,1)=xprod_aux1(2)*R1
	curlSka(3,1)=xprod_aux1(3)*R1


	curlSka(1,2)=xprod_aux2(1)*R1
  curlSka(2,2)=xprod_aux2(2)*R1
	curlSka(3,2)=xprod_aux2(3)*R1

return
end subroutine get_curlSka




subroutine get_curlcurlSka(ru_s,rv_s,n_s,dr,R1,R2,zk,my_exp,r,curlcurlSka)
implicit none

	!List of calling arguments
	real ( kind = 8 ), intent(in) :: ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ), intent(in) :: dr(3),r
	complex ( kind = 8 ), intent(in) :: R1,R2,zk,my_exp
	complex ( kind = 8 ), intent(out) :: curlcurlSka(3,2)
	
 	!List of local variables
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)
	complex ( kind = 8 ) Skb(3,2)	
	
	Skb(1,1)=ru_s(1)*my_exp/r
	Skb(2,1)=ru_s(2)*my_exp/r
	Skb(3,1)=ru_s(3)*my_exp/r

	Skb(1,2)=rv_s(1)*my_exp/r
	Skb(2,2)=rv_s(2)*my_exp/r
	Skb(3,2)=rv_s(3)*my_exp/r

	curlcurlSka(1,1)=zk**2*Skb(1,1)+R1*ru_s(1)-dr(1)*DOT_PRODUCT(ru_s,-dr)*R2
	curlcurlSka(2,1)=zk**2*Skb(2,1)+R1*ru_s(2)-dr(2)*DOT_PRODUCT(ru_s,-dr)*R2
	curlcurlSka(3,1)=zk**2*Skb(3,1)+R1*ru_s(3)-dr(3)*DOT_PRODUCT(ru_s,-dr)*R2

	curlcurlSka(1,2)=zk**2*Skb(1,2)+R1*rv_s(1)-dr(1)*DOT_PRODUCT(rv_s,-dr)*R2
	curlcurlSka(2,2)=zk**2*Skb(2,2)+R1*rv_s(2)-dr(2)*DOT_PRODUCT(rv_s,-dr)*R2
	curlcurlSka(3,2)=zk**2*Skb(3,2)+R1*rv_s(3)-dr(3)*DOT_PRODUCT(rv_s,-dr)*R2

!	nxcurlcurlSka(1,1)=nxcurlcurlSka(1,1)+(-DOT_PRODUCT(rv_t,ru_s)*R1)
!	nxcurlcurlSka(1,2)=nxcurlcurlSka(1,2)+(-DOT_PRODUCT(rv_t,rv_s)*R1)
!	nxcurlcurlSka(2,1)=nxcurlcurlSka(2,1)+(+DOT_PRODUCT(ru_t,ru_s)*R1)
!	nxcurlcurlSka(2,2)=nxcurlcurlSka(2,2)+(+DOT_PRODUCT(ru_t,rv_s)*R1)
	
!	nxcurlcurlSka(1,1)=nxcurlcurlSka(1,1)-(-DOT_PRODUCT(rv_t,dr)*DOT_PRODUCT(ru_s,-dr)*R2)
!	nxcurlcurlSka(1,2)=nxcurlcurlSka(1,2)-(-DOT_PRODUCT(rv_t,dr)*DOT_PRODUCT(rv_s,-dr)*R2)
!	nxcurlcurlSka(2,1)=nxcurlcurlSka(2,1)-(+DOT_PRODUCT(ru_t,dr)*DOT_PRODUCT(ru_s,-dr)*R2)
!	nxcurlcurlSka(2,2)=nxcurlcurlSka(2,2)-(+DOT_PRODUCT(ru_t,dr)*DOT_PRODUCT(rv_s,-dr)*R2)

return
end subroutine get_curlcurlSka





      subroutine lpcomp_em_muller_far_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,E_far,H_far,wnear,&
     &n_regions,ntarg_vect)

!
!  This subroutine evaluates the layer potential for
!  the DFIE boundary integral equation:
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  Note the 4\pi scaling is NOT included as the FMM output was scaled
!  appropriately
!
!  Note: the identities are not included as the gmres takes care of that
!
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
!
!
!  input:
!    npatches - integer
!      number of patches
!
!    norders- integer(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(npatches+1)
!      ixyzs(i) denotes the starting location in srccoefs,
!      and srcvals array corresponding to patch i
!   
!    iptype - integer(npatches)
!      type of patch
!      iptype = 1, triangular patch discretized using RV nodes
!
!    npts - integer
!      total number of discretization points on the boundary
! 
!    srccoefs - real *8 (9,npts)
!      koornwinder expansion coefficients of xyz, dxyz/du,
!      and dxyz/dv on each patch. 
!      For each point srccoefs(1:3,i) is xyz info
!        srccoefs(4:6,i) is dxyz/du info
!        srccoefs(7:9,i) is dxyz/dv info
!
!    srcvals - real *8 (12,npts)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
! 
!    ndtarg - integer
!      leading dimension of target array
!        
!    ntarg - integer
!      number of targets
!
!    targs - real *8 (ndtarg,ntarg)
!      target information
!
!    eps - real *8
!      precision requested
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = omega 
!      zpars(2) = ep0
!      zpars(3) = mu0
!      zpars(4) = ep1
!      zpars(5) = mu1
!
!    nnz - integer *8
!      number of source patch-> target interactions in the near field
! 
!    row_ptr - integer(ntarg+1)
!      row_ptr(i) is the pointer
!      to col_ind array where list of relevant source patches
!      for target i start
!
!    col_ind - integer (nnz)
!      list of source patches relevant for all targets, sorted
!      by the target number
!
!    iquad - integer(nnz+1)
!      location in wnear array where quadrature for col_ind(i)
!      starts
!
!    nquad - integer
!      number of entries in wnear
!
!    wnear  - complex *16(36*nquad)
!      near field precomputed quadrature
!
!    sigma - complex *16(2*ns)
!      induced charge and current on the surface
!      sigma(1:ns) - first component of 'a' along
!        the srcvals(4:6,i) direction
!      sigma(ns+1:2*ns) - second component of 'a' along
!        the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      sigma(2*ns+1:3*ns) - scalar sigma on the surface
!      sigma(3*ns+1:4*ns) - first component of 'b' along
!        the srcvals(4:6,i) direction
!      sigma(4*ns+1:5*ns) - second component of 'b' along
!        the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      sigma(5*ns+1:6*ns) - scalar rho on the surface
!
!    novers - integer(npatches)
!      order of discretization for oversampled sources and density
!
!    ixyzso - integer(npatches+1)
!      ixyzso(i) denotes the starting location in srcover,
!      corresponding to patch i
!   
!    nptso - integer
!      total number of oversampled points
!
!    srcover - real *8 (12,nptso)
!      oversampled set of source information
!
!    whtsover - real *8 (nptso)
!      smooth quadrature weights at oversampled nodes
!

      implicit none
      integer npatches,norder,npols,npts
      integer ndtarg,ntarg
      integer norders(npatches),ixyzs(npatches+1)
      integer ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(3),zpars_aux(5)
      integer nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      complex *16 sigma(4*npts),sigma2(npts)
      integer n_regions
      integer ntarg_vect(n_regions+1)
	  
      complex *16 wnear(12*nquad)
	  
      integer novers(npatches+1)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 E_far(3,ntarg),H_far(3,ntarg)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(6)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart,count1,count2,istart,ifinish


      integer ifaddsub,ifdir

      integer ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp2_s(:),dtmp2(:,:)
      complex *16, allocatable :: ctmp2_u(:),ctmp2_v(:)
      complex *16, allocatable :: ctmp2_a_u(:),ctmp2_a_v(:)
      complex *16, allocatable :: ctmp2_b_u(:),ctmp2_b_v(:)

	    complex *16, allocatable :: pot_aux(:)

      real *8 radexp,epsfmm

      integer ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      integer nss,ii,l,npover
      complex *16 ima
	    complex *16 omega,ep,mu

      integer nd,ntarg0
      integer icount

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)

      ns = nptso
      done = 1
      pi = atan(done)*4
	    ima=(0.0d0,1.0d0)

           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(3,ns),targvals(3,ntarg))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(4*ns))
	    allocate(pot_aux(4*ntarg))

! 
!       oversample density

    call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,&
	&npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))
	       
    call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,&
	&npts,sigma(npts+1:2*npts),novers,ixyzso,ns,sigmaover(ns+1:2*ns))

	call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
	&npts,sigma(2*npts+1:3*npts),novers,ixyzso,ns,sigmaover(2*ns+1:3*ns))

	call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
	&npts,sigma(3*npts+1:4*npts),novers,ixyzso,ns,sigmaover(3*ns+1:4*ns))

    call prin2('sigma=*',sigma,48)

      ra = 0

!
!       fmm
!

		call get_fmm_thresh(12,ns,srcover,ndtarg,ntarg,targs,thresh)
       ifdir=0
	  
    istart=1
    ifinish=ntarg_vect(1)

    print *, "here"
    write (*,*) 'n_regions: ', n_regions
    write (*,*) 'thresh: ', thresh

    do count1=1,n_regions
      zpars_aux(1)=zpars(1)
      zpars_aux(2)=targs(4,istart)+ima*targs(5,istart)
      zpars_aux(3)=targs(6,istart)+ima*targs(7,istart)

     ! write (*,*) 'contador regiones fmm: ',count1,zpars_aux(1:3)
     ! write (*,*) targs(:,istart)
      !Calculate the far_field with FMM

       call prin2('zpars_aux=*',zpars_aux,6)
       call prin2('eps=*',eps,1)
       call prinf('ns=*',ns,1)
       print *, "here2"
       call prinf('ntarg_vect=*',ntarg_vect(count1),1)

       call em_muller_far_FMM(eps,zpars_aux,ns,ntarg_vect(count1),srcover,ndtarg,targs(:,istart:ifinish),whtsover,&
      &sigmaover(1:ns),sigmaover(ns+1:2*ns),sigmaover(2*ns+1:3*ns),&
      &sigmaover(3*ns+1:4*ns),E_far(:,istart:ifinish),H_far(:,istart:ifinish),&
      &thresh,ifdir)

      if (count1<n_regions) then
        istart=ifinish+1
        ifinish=istart+ntarg_vect(count1+1)-1
      endif
    enddo
!!!!return
!write (*,*) 'aquí no llega'
!read (*,*)

      call prin2('E_far=*',E_far,6)
      call prin2('H_far=*',H_far,6)
      read * ,i
		call cpu_time(t1)
!C$      t1 = omp_get_wtime()
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
!C$OMP$PRIVATE(jstart,pottmp,npols,ep,mu,l)
      do i=1,ntarg
        ep=targs(4,i)+ima*targs(5,i)
        mu=targs(6,i)+ima*targs(7,i)
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols

              E_far(1,i) = E_far(1,i) + mu*wnear(jquadstart+l-1)*sigma(jstart+l-1)
              E_far(2,i) = E_far(2,i) + mu*wnear(4*nquad+jquadstart+l-1)*sigma(jstart+l-1)
              E_far(3,i) = E_far(3,i) + mu*wnear(8*nquad+jquadstart+l-1)*sigma(jstart+l-1)

              E_far(1,i) = E_far(1,i) + mu*wnear(nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)
              E_far(2,i) = E_far(2,i) + mu*wnear(5*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)
              E_far(3,i) = E_far(3,i) + mu*wnear(9*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)


              H_far(1,i) = H_far(1,i) + ep*wnear(jquadstart+l-1)*sigma(jstart+l-1+2*npts)
              H_far(2,i) = H_far(2,i) + ep*wnear(4*nquad+jquadstart+l-1)*sigma(jstart+l-1+2*npts)
              H_far(3,i) = H_far(3,i) + ep*wnear(8*nquad+jquadstart+l-1)*sigma(jstart+l-1+2*npts)

              H_far(1,i) = H_far(1,i) + ep*wnear(nquad+jquadstart+l-1)*sigma(jstart+l-1+3*npts)
              H_far(2,i) = H_far(2,i) + ep*wnear(5*nquad+jquadstart+l-1)*sigma(jstart+l-1+3*npts)
              H_far(3,i) = H_far(3,i) + ep*wnear(9*nquad+jquadstart+l-1)*sigma(jstart+l-1+3*npts)

              E_far(1,i) = E_far(1,i) + wnear(2*nquad+jquadstart+l-1)*sigma(jstart+l-1+2*npts)
              E_far(2,i) = E_far(2,i) + wnear(6*nquad+jquadstart+l-1)*sigma(jstart+l-1+2*npts)
              E_far(3,i) = E_far(3,i) + wnear(10*nquad+jquadstart+l-1)*sigma(jstart+l-1+2*npts)

              E_far(1,i) = E_far(1,i) + wnear(3*nquad+jquadstart+l-1)*sigma(jstart+l-1+3*npts)
              E_far(2,i) = E_far(2,i) + wnear(7*nquad+jquadstart+l-1)*sigma(jstart+l-1+3*npts)
              E_far(3,i) = E_far(3,i) + wnear(11*nquad+jquadstart+l-1)*sigma(jstart+l-1+3*npts)

              H_far(1,i) = H_far(1,i) - wnear(2*nquad+jquadstart+l-1)*sigma(jstart+l-1)
              H_far(2,i) = H_far(2,i) - wnear(6*nquad+jquadstart+l-1)*sigma(jstart+l-1)
              H_far(3,i) = H_far(3,i) - wnear(10*nquad+jquadstart+l-1)*sigma(jstart+l-1)

              H_far(1,i) = H_far(1,i) - wnear(3*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)
              H_far(2,i) = H_far(2,i) - wnear(7*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)
              H_far(3,i) = H_far(3,i) - wnear(11*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)

		        !do count1=0,2
			      ! do count2=0,3
            !  E_far(count1+1,i) = E_far(count1+1,i) + &
				    !   &wnear((count1*4+count2)*nquad+jquadstart+l-1)*&
				    !   &sigma(jstart+l-1+npts*count2)
            !  H_far(count1+1,i) = H_far(count1+1,i) + &
				    !   &wnear((count1*4+count2+12)*nquad+jquadstart+l-1)*&
				    !   &sigma(jstart+l-1+npts*count2)
			      ! enddo
			      !enddo

          enddo
        enddo
      enddo

!!C$OMP END PARALLEL DO
!!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
!		ipars(1)=1
!		ipars(2)=1
	ifdir=1
!      do i=1,ntarg
      i=0
      do count1=1,n_regions
       zpars_aux(1)=zpars(1)
       zpars_aux(2)=targs(4,i+1)+ima*targs(5,i+1)
       zpars_aux(3)=targs(6,i+1)+ima*targs(7,i+1)
       do count2=1,ntarg_vect(count1)
        i=i+1
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),wtmp2(nss))
        allocate(ctmp2_a_u(nss),ctmp2_a_v(nss))
        allocate(ctmp2_b_u(nss),ctmp2_b_v(nss))

        rmin = 1.0d6
        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
			      ii = ii+1
			      srctmp2(:,ii) = srcover(:,jstart+l)
			      ctmp2_a_u(ii)=sigmaover(jstart+l)
			      ctmp2_a_v(ii)=sigmaover(jstart+l+ns)			
			      ctmp2_b_u(ii)=sigmaover(jstart+l+2*ns)
			      ctmp2_b_v(ii)=sigmaover(jstart+l+3*ns)
			      wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo
        call em_muller_far_FMM(eps,zpars_aux,nss,ntarg0,srctmp2,ndtarg,targs(:,i),&
        &wtmp2,ctmp2_a_u,ctmp2_a_v,ctmp2_b_u,ctmp2_b_v,&
        &E(1:3),E(4:6),thresh,ifdir)

          E_far(1,i)=E_far(1,i)-E(1)
          E_far(2,i)=E_far(2,i)-E(2)
          E_far(3,i)=E_far(3,i)-E(3)

          H_far(1,i)=H_far(1,i)-E(4)
          H_far(2,i)=H_far(2,i)-E(5)
          H_far(3,i)=H_far(3,i)-E(6)

        deallocate(srctmp2,wtmp2)
        deallocate(ctmp2_a_u,ctmp2_a_v,ctmp2_b_u,ctmp2_b_v)

       enddo
      enddo

	  
	  return
    end subroutine lpcomp_em_muller_far_addsub



subroutine em_muller_far_FMM(eps,zpars,ns,nt,srcvals,ndtarg,targvals,wts,&
 &a_u,a_v,b_u,b_v,E_far,H_far,thresh,ifdir)
implicit none

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zpars(3)
    integer, intent(in) :: ns,nt,ndtarg
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns)
    real ( kind = 8 ), intent(in) :: targvals(ndtarg,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns)
    complex ( kind = 8 ), intent(in) :: b_u(ns),b_v(ns)
    complex ( kind = 8 ), intent(out) :: E_far(3,nt),H_far(3,nt)
    real ( kind = 8 ), intent(in) :: thresh
    integer, intent(in) :: ifdir 


    !List of local variables
	  real ( kind = 8 ), allocatable :: u_vect_s(:,:),v_vect_s(:,:)
	  real ( kind = 8 ), allocatable :: source(:,:),n_vect_s(:,:)

	  real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:)
	  real ( kind = 8 ), allocatable :: targets(:,:),v_vect_t(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: E(:,:),curlE(:,:),divE(:)
    complex ( kind = 8 ) ima,zk

    complex ( kind = 8 ) omega,ep,mu

    integer count1,count2
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

	  ima=(0.0d0,1.0d0)

	  omega=zpars(1)
	  ep=zpars(2)
	  mu=zpars(3)
	 
	  zk=omega*sqrt(ep*mu)



    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
    allocate(rho(ns))
    allocate(E(3,nt))
    allocate(curlE(3,nt))
    allocate(divE(nt))
    allocate(n_vect_s(3,ns))
    allocate(n_vect_t(3,nt))
    allocate(u_vect_s(3,ns))
    allocate(v_vect_s(3,ns))
    allocate(u_vect_t(3,nt))
    allocate(v_vect_t(3,nt))
    allocate(source(3,ns))
    allocate(targets(3,nt))

    do count1=1,ns
      n_vect_s(:,count1)=srcvals(10:12,count1)
      source(:,count1)=srcvals(1:3,count1)
   enddo
   do count1=1,nt
      targets(:,count1)=targvals(1:3,count1)
   enddo
   call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect_s,&
      &v_vect_s,ns)

    do count1=1,ns
      a_vect(1,count1)=(b_u(count1)*u_vect_s(1,count1)+b_v(count1)*v_vect_s(1,count1))/(-ima*omega)
      a_vect(2,count1)=(b_u(count1)*u_vect_s(2,count1)+b_v(count1)*v_vect_s(2,count1))/(-ima*omega)
      a_vect(3,count1)=(b_u(count1)*u_vect_s(3,count1)+b_v(count1)*v_vect_s(3,count1))/(-ima*omega)

      b_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*v_vect_s(1,count1))*mu
      b_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*v_vect_s(2,count1))*mu
      b_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*v_vect_s(3,count1))*mu
    enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0
    call prin2('zk=*',zk,2)
    call prin2('omega=*',omega,1)
    call prin2('a_vect=*',a_vect,24)
    call prin2('b_vect=*',b_vect,24)
    call prin2('wts=*',wts,24)
    call prinf('ifdir=*',ifdir,1)
    print *, "Here"

     call Vector_Helmholtz_targ2(eps,zk,ns,source,wts,ifa_vect,a_vect,&
      &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
      &E_far,ifdivE,divE,nt,targets,thresh,ifdir)

     call prin2('E_far=*',E_far,6)

!	call Vector_Helmholtz_targ(eps,izk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
!	 &b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,curlE,ifdivE,divE,nt,targets)


	  do count1=1,ns
		    b_vect(1,count1)=(b_u(count1)*u_vect_s(1,count1)+b_v(count1)*v_vect_s(1,count1))*ep
		    b_vect(2,count1)=(b_u(count1)*u_vect_s(2,count1)+b_v(count1)*v_vect_s(2,count1))*ep
		    b_vect(3,count1)=(b_u(count1)*u_vect_s(3,count1)+b_v(count1)*v_vect_s(3,count1))*ep
				
        a_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*v_vect_s(1,count1))/(ima*omega)
        a_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*v_vect_s(2,count1))/(ima*omega)
        a_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*v_vect_s(3,count1))/(ima*omega)
	  enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

	  call Vector_Helmholtz_targ2(eps,zk,ns,source,wts,ifa_vect,a_vect,&
	   &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	   &H_far,ifdivE,divE,nt,targets,thresh,ifdir)

	
	deallocate(a_vect)
	deallocate(b_vect)
	deallocate(lambda)
	deallocate(rho)
	deallocate(E)
	deallocate(curlE)
	deallocate(divE)
	deallocate(u_vect_s)
	deallocate(v_vect_s)
	deallocate(n_vect_s)
	deallocate(source)
	deallocate(u_vect_t)
	deallocate(v_vect_t)
	deallocate(n_vect_t)
	deallocate(targets)

return
end subroutine em_muller_far_FMM








subroutine test_accuracy_em_muller(npatches,norders,ixyzs,iptype,npts,srccoefs,&
    &srcvals,wts,targ,ntarg,npatches_vect,n_components,sorted_vector,&
    &contrast_matrix,exposed_surfaces,eps,zpars,sigma,P0,vf,direction, Pol,nx,ny)
implicit none

 !List of calling arguments
 integer, intent(in) :: npatches,npts,n_components,ntarg,nx,ny
 integer, intent(in) :: norders(npatches),npatches_vect(n_components),ixyzs(npatches+1),iptype(npatches)
 real ( kind = 8 ), intent(in) :: srcvals(12,npts), srccoefs(9,npts),targ(3,ntarg),wts(npts)
 integer, intent(in) :: sorted_vector(n_components+1)
 complex ( kind = 8 ), intent(in) :: contrast_matrix(4,n_components),zpars(3),sigma(4*npts)
 logical *8 exposed_surfaces(n_components)
 real ( kind = 8 ), intent(in) :: eps
 real ( kind = 8 ), intent(in) :: P0(3), direction(2)
 complex ( kind = 8 ), intent(in) :: vf(3),Pol(2)

 !List of local variables
 integer count1,count2,icount,icount2,n_aux,i1,i2,j1,j2,npatches_aux,npts_aux,x,ntarg_vect(n_components+1)
 integer ndtarg,n_regions
 integer, allocatable :: location_targs(:),permutation_targs(:)
 integer, allocatable :: ipatch_id(:)
 real *8, allocatable :: uvs_targ(:,:)

 real *8, allocatable :: targ_sort(:,:),error_E(:),error_H(:),error_rel_E(:),error_rel_H(:)
 complex ( kind = 8 ) ep0,mu0,ep,mu
 complex ( kind = 8 ), allocatable :: E_far(:,:), H_far(:,:),E_0(:,:), H_0(:,:)
 complex ( kind = 8 ), allocatable :: E_far_aux(:,:),H_far_aux(:,:)
 character (len=100) nombre_plot
 integer ( kind = 8 ) M_plot,N_plot


 allocate(location_targs(ntarg),permutation_targs(ntarg))
 allocate(ipatch_id(ntarg))
 allocate(uvs_targ(2,ntarg))
 allocate(targ_sort(7,ntarg))
 allocate(E_far(3,ntarg),H_far(3,ntarg))
 allocate(E_far_aux(3,ntarg),H_far_aux(3,ntarg))
 allocate(E_0(3,ntarg),H_0(3,ntarg))
 allocate(error_E(ntarg),error_H(ntarg))
 allocate(error_rel_E(ntarg),error_rel_H(ntarg))

  do count1=1,n_components+1
    ntarg_vect(count1)=0
  enddo

  call find_inclusion_vect(npatches,norders,ixyzs,iptype,npts,srccoefs,&
    &srcvals,wts,targ,ntarg,npatches_vect,n_components,sorted_vector,&
    &location_targs,eps)
  do count1=1,n_components
    if (exposed_surfaces(count1)) then
       ep0=contrast_matrix(1,count1)
       mu0=contrast_matrix(2,count1)
       exit
    endif
  enddo
 icount=1
 n_regions=0
 do count1=0,n_components
   print *, "n_regions=",n_regions
   print *, "count1=",count1
   icount2=0
   do count2=1,ntarg
    if (location_targs(count2).eq.count1) then
      permutation_targs(count2)=icount
      targ_sort(1:3,icount)=targ(:,count2)
      if (count1.eq.0) then
        targ_sort(4,icount)=real(ep0)
        targ_sort(5,icount)=aimag(ep0)
        targ_sort(6,icount)=real(mu0)
        targ_sort(7,icount)=aimag(mu0)
      else
        targ_sort(4,icount)=real(contrast_matrix(3,count1))
        targ_sort(5,icount)=aimag(contrast_matrix(3,count1))
        targ_sort(6,icount)=real(contrast_matrix(4,count1))
        targ_sort(7,icount)=aimag(contrast_matrix(4,count1))
      endif
      icount=icount+1
      icount2=icount2+1
    endif
   enddo
  if (icount2>0) then
    
    n_regions=n_regions+1
    ntarg_vect(n_regions)=icount2
  endif
 enddo
 ndtarg=7
 do count1=1,ntarg
   ipatch_id(count1) = -1
   uvs_targ(1,count1) = 0
   uvs_targ(2,count1) = 0
 enddo
! do count1=1,ntarg
!   write (*,*) 'targ_sort: ',count1,targ_sort(:,count1)
! enddo
!  do count1=1,ntarg
!   write (*,*) 'targ: ',count1,targ(:,count1)
!  enddo

 E_far_aux = 0
 H_far_aux = 0
 print *, "Entering routine in testing code"
 call prinf('ndtarg=*',ndtarg,1)
 call prinf('ntarg=*',ntarg,1)
 call prinf('ntarg_vect=*',ntarg_vect,n_regions)
 call lpcomp_em_muller_far_dir(npatches,norders,ixyzs,&
  &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targ_sort,&
  &ipatch_id,uvs_targ,eps,zpars,sigma,E_far_aux,H_far_aux,n_regions,ntarg_vect)
 

 call prin2('E_far_aux=*',E_far_aux,6)
 call prin2('H_far_aux=*',H_far_aux,6)
 
 do count1=1,ntarg
   E_far(:,count1)=E_far_aux(:,permutation_targs(count1))
   H_far(:,count1)=H_far_aux(:,permutation_targs(count1))   
 enddo

 do count1=1,ntarg
  if (location_targs(count1).eq.0) then
!    write (*,*) 'datum dipoles: ', zpars(1),ep0,mu0,P0,targ(1:3,count1),vf
!    read (*,*)
    call fieldsEDomega(zpars(1),ep0,mu0,P0,targ(1:3,count1),1,&
     &E_0(:,count1),H_0(:,count1),vf,0)
    call fieldsMDomega(zpars(1),ep0,mu0,P0,targ(1:3,count1),1,&
     &E_0(:,count1),H_0(:,count1),vf,1)
  else
    ep=contrast_matrix(3,location_targs(count1))
    mu=contrast_matrix(4,location_targs(count1))
    call fieldsPWomega(zpars(1),ep,mu,targ(1:3,count1),1,&
     &E_0(:,count1),H_0(:,count1),direction,Pol)
  endif
    error_E(count1)=sqrt(abs(E_0(1,count1)-E_far(1,count1))**2&
     &+abs(E_0(2,count1)-E_far(2,count1))**2+abs(E_0(3,count1)-E_far(3,count1))**2)
    error_rel_E(count1)=error_E(count1)/sqrt(abs(E_0(1,count1))**2&
     &+abs(E_0(2,count1))**2+abs(E_0(3,count1))**2)

    error_H(count1)=sqrt(abs(H_0(1,count1)-H_far(1,count1))**2+&
     &abs(H_0(2,count1)-H_far(2,count1))**2+abs(H_0(3,count1)-H_far(3,count1))**2)
    error_rel_H(count1)=error_H(count1)/sqrt(abs(H_0(1,count1))**2+&
     &abs(H_0(2,count1))**2+abs(H_0(3,count1))**2)
 enddo
 call prin2('E0=*',E_0,6)
 call prin2('E_far=*',E_far,6)
 call prin2('E error=*',E_0-E_far,6)
 call prin2('H0=*',H_0,6)
 call prin2('H_far=*',H_far,6)
 call prin2('H error=*',H_0-H_far,6)

      write (*,*) 'Relative Error in E: ', error_rel_E
      write (*,*) 'Relative Error in H: ', error_rel_H

      write (*,*) 'Relative Error in E (maxval): ', maxval(error_rel_E)
      write (*,*) 'Relative Error in H (maxval): ', maxval(error_rel_H)

! do count1=1,ntarg
!   write (*,*) 'targ : ',count1,targ(:,count1)
!   write (*,*) 'E_far: ',count1,E_far(:,count1)
!   write (*,*) 'E_0  : ',count1,E_0(:,count1)
!   write (*,*) ' '
! enddo
! write (*,*) ' '
! do count1=1,ntarg
!   write (*,*) 'targ : ',count1,targ(:,count1)
!   write (*,*) 'H_far: ',count1,H_far(:,count1)
!   write (*,*) 'H_0  : ',count1,H_0(:,count1)
!   write (*,*) ' '
! enddo

 nombre_plot='plot_err_E_fortran.txt'
 M_plot=int(nx,8)
 N_plot=int(ny,8)
 call plot2D(targ(1,:),targ(2,:),error_rel_E,M_plot,N_plot,nombre_plot)

 nombre_plot='plot_err_H_fortran.txt'
 call plot2D(targ(1,:),targ(2,:),error_rel_H,M_plot,N_plot,nombre_plot)

 nombre_plot='plot_E_fortran.txt'
 call plot2D(targ(1,:),targ(2,:),real(E_far(3,:)),M_plot,N_plot,nombre_plot)

 nombre_plot='plot_E0_fortran.txt'
 call plot2D(targ(1,:),targ(2,:),real(E_0(3,:)),M_plot,N_plot,nombre_plot)

 nombre_plot='plot_H_fortran.txt'
 call plot2D(targ(1,:),targ(2,:),real(H_far(3,:)),M_plot,N_plot,nombre_plot)

 nombre_plot='plot_H0_fortran.txt'
 call plot2D(targ(1,:),targ(2,:),real(H_0(3,:)),M_plot,N_plot,nombre_plot)


return
end subroutine test_accuracy_em_muller

