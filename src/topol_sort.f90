subroutine text_process(string1,n,iwords)
implicit none
!  This function locates the position of n characters ? on a string of 
!  size 2000 where n is known
!
!  input:
!    string1 - character *2000
!      e.g. 1234?67?9?
!
!    n - integer
!      number of characters '?' contained in the input string
!
!    iwords - integer(n+1)
!      location of the characters '?'
!      in the example 1234?67?9? iwords=(/0,5,8,10/)
!      notice that iwords(1)=0 always
!  

 !List of calling arguments
 character *2000 string1 
 integer, intent(in) :: n
 integer, intent(out) :: iwords(n+1)

 !List of local variables
 integer count1,count2,icount,n_aux,i1,i2,j1,j2,npatches_aux,npts_aux

  iwords(1)=0
  i1=0
  do count1=1,n
     iwords(count1+1) = index(string1(iwords(count1)+1:2000), '?') + &
         iwords(count1)
  enddo
  
return
end subroutine text_process


subroutine topological_sorting(npatches,norders,ixyzs,iptype,npts, &
 srccoefs,srcvals,wts,npatches_vect,n_components,sorted_vector, &
 exposed_surfaces,eps)
implicit none
!
!  This function sorts a set of non intersecting connected surfaces in 
!  a way that is compatible with the partial ordering defined by the 
!  inclusion relation 
!  Also provides which surfaces are exposed to the ambient space 
!  (that is, surfaces that are not in the interior of any other surface)
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
!    n_components - integer
!      number of different connected components of the geometry
!      that is, number of different interfaces in the transmission
!      problem 
!
!    npatches_vect - integer(n_components)
!      number of patches that forms if the i'th component 
!        (the i'th interface)
!      the first interface is formed by patches 1,...,npatcehs_vect(1)
!      the second interface is formed by patches:
!        npatcehs_vect(1)+1,....,npatcehs_vect(1)+npatcehs_vect(2)
!      the third interface is formed by patches:
!        npatcehs_vect(1)+npatcehs_vect(2)+1,..
!        ...,npatcehs_vect(1)+npatches_vect(2)+npatches_vect(3)
!
!    eps real*8
!      accuracy for the evaluation of the double layer
!
!  output
!    sorted_vector - integer(n_components+1)
!     resulting ordered components (interfaces)
!     sorted_vector(1) is always a surface that does not contain 
!     any other surface in the interior
!     sorted_vector(n_components+1) = 0 always (meaning a big sphere at 
!     infinity that contains everything)
!
!    exposed_surfaces - logical(n_components)
!     exposed_surfaces(i)=.true. if the i'th interface is exposed to the 
!     ambient space, that is, is not contained
!     in the interior of any other surface.
!     this is useful and used to compute the RHS for an incoming plane wave.
!     exposed_surfaces(i)=.false. otherwise
!

 !List of calling arguments
 integer, intent(in) :: npatches,n_components,npts
 integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
 integer, intent(in) :: iptype(npatches)
 integer, intent(in) :: npatches_vect(n_components)
 real ( kind = 8 ), intent(in) :: srcvals(12,npts),eps
 real ( kind = 8 ), intent(in) :: srccoefs(9,npts),wts(npts)
 integer, intent(out) :: sorted_vector(n_components+1)
 logical *8 exposed_surfaces(n_components)

 !List of local variables
 integer count1,count2,icount,n_aux
 real ( kind = 8 ) pt(3)
 logical *8 poset_matrix(n_components+1,n_components+1)
 logical *8 mascara(n_components+1)

  poset_matrix(:,:)=.false.
  n_aux=0
  do count1=1,n_components
    n_aux=n_aux+npatches_vect(count1)
    pt=srcvals(1:3,ixyzs(n_aux+1)-1)
    call find_inclusion(npatches,norders,ixyzs,iptype,npts,srccoefs,&
    &srcvals,wts,pt,npatches_vect,n_components,count1,&
    &poset_matrix(count1+1,:),eps)
  enddo
  do count1=1,n_components+1
    poset_matrix(count1,count1)=.false.
  enddo
  exposed_surfaces(:)=.false.
  do count1=1,n_components
    if (all(.not.poset_matrix(count1+1,2:n_components+1))) then
       exposed_surfaces(count1)=.true.
    endif
  enddo

  poset_matrix(:,:)=.not.poset_matrix(:,:)
  mascara(:)=.false.
  icount=1
  do count2=1,n_components+1
    do count1=1,n_components+1
      if (.not.mascara(count1)) then
        if (all(poset_matrix(:,count1).or.mascara)) then
          sorted_vector(icount)=count1-1
          mascara(count1)=.true.
          icount=icount+1
          exit
        endif
      endif
    enddo
  enddo
  
return
end subroutine topological_sorting



subroutine find_inclusion(npatches,norders,ixyzs,iptype,npts,srccoefs, &
   srcvals,wts,pt,npatches_vect,n_components,i_component, &
   vector_location,eps)
implicit none
!
!  This function finds the location of a point pt relative to a set of
!  surfaces vector_location is a logical output array that is 
!  vector_location(i)=.true. if the point is inside the i'th component
!  (interface) of the geometry. The calculation is done by direct sum,
!  not fmm. This function is used by topological_sorting to find the
!  partial ordering of the set of components (or interfaces) in the
!  geometry
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
!    pt - real *8 (3)
!      xyz coordinates of the point to locate
!
!    n_components - integer
!      number of different connected components of the geometry
!      that is, number of different interfaces in the transmission
!      problem 
!
!    npatches_vect - integer(n_components)
!      number of patches that forms if the i'th component 
!        (the i'th interface)
!      the first interface is formed by patches 1,...,npatcehs_vect(1)
!      the second interface is formed by patches:
!        npatcehs_vect(1)+1,....,npatcehs_vect(1)+npatcehs_vect(2)
!      the third interface is formed by patches:
!        npatcehs_vect(1)+npatcehs_vect(2)+1,...
!        ....,npatcehs_vect(1)+npatches_vect(2)+npatches_vect(3)
!
!    i_component integer
!      component of the geometry to which pt belongs (that way we 
!      avoid to locate a surface with respect to itself)
!
!
!  output
!    vector_location logical(n_components+1)
!      relative location of the point pt with respect to each 
!       component (or interface)
!      vector_location(i+1)=.true. if pt is inside the i'th interface
!      vectro-location(1)=.true. always as pt is always inside a 
!       big sphere at infinity
!

 !List of calling arguments
 integer, intent(in) :: npatches,npts,n_components,i_component
 integer, intent(in) :: norders(npatches),npatches_vect(n_components)
 integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
 real ( kind = 8 ), intent(in) :: srcvals(12,npts)
 real ( kind = 8 ), intent(in) :: srccoefs(9,npts),pt(3),wts(npts),eps

 logical *8, intent(out) :: vector_location(n_components+1)

 !List of local variables
 integer count1,count2,icount,n_aux,i1,i2,j1,j2,npatches_aux
 integer npts_aux,ipatch_id_aux
 integer, allocatable :: ixyzs_aux(:)
 real *8, allocatable :: sigma(:)
 
 integer uvs_targ_aux(2,1)
 real *8 pot,dpars(2)
 real *8 ttot,done,pi
 done = 1
 pi = atan(done)*4

 allocate(sigma(npts))
 ipatch_id_aux=-1
 sigma(:)=-1.0d0
 uvs_targ_aux(1,1)=0
 uvs_targ_aux(2,1)=0
 dpars(1)=0.0d0
 dpars(2)=1.0d0

 vector_location(1)=.true.

 i1=1
 n_aux=0
 do count1=1,n_components
  i1=n_aux+1
  n_aux=n_aux+npatches_vect(count1)
  npatches_aux=npatches_vect(count1)
  i2=n_aux
  if (count1.ne.i_component) then
    allocate(ixyzs_aux(npatches_aux+1))
    ixyzs_aux=ixyzs(i1:(i2+1))-ixyzs(i1)+1
    npts_aux=ixyzs_aux(npatches_aux+1)-1
    
!
!  call laplace double layer potential evaluator where 
!  far field is computed directly
!
!

    call lpcomp_lap_comb_dir_nonfmm(npatches_aux,norders(i1:i2),ixyzs_aux,&
    &iptype(i1:i2),npts_aux,srccoefs(1:9,ixyzs(i1):(ixyzs(i2+1)-1)),&
    &srcvals(1:12,ixyzs(i1):(ixyzs(i2+1)-1)),3,1,pt,&
    &ipatch_id_aux,uvs_targ_aux,eps,dpars,&
    &sigma(ixyzs(i1):(ixyzs(i2+1)-1)),pot)
!
!  
!

    if (pot.le.0.5d0) then 
      vector_location(count1+1)=.false.
    else
      vector_location(count1+1)=.true.
    endif
    deallocate(ixyzs_aux)
  endif
 enddo
 deallocate(sigma)

return
end subroutine find_inclusion


subroutine point_location(n,vector,target_pt,sol)
implicit none
!
!  This function finds the location of a point whos relative position
!  with respect to each interface is given by target_pt
!  Sol is an integer that indicates the smallest interface (component)
!  in which the point is. This is used to find ep,mu associated to
!  that point
!
!  input:
!    n - integer
!      number of interfaces
!
!    vector - integer (n)
!      topological sorted list of interfaces (linear ordering that 
!      respects the partial ordering)
!
!    target_pt - logical (n)
!      indicates the relative position of the point with respect to
!      each interface
!
!

 !List of calling arguments
 integer, intent(in) :: n
 integer, intent(in) :: vector(n)
 logical *8, intent(in) :: target_pt(n)
 integer, intent(out) :: sol

 !List of local variables
 integer count1

  if (all(.not.target_pt)) then
   sol=0
  else
   do count1=1,n
    if (target_pt(vector(count1))) then
      sol=vector(count1)
      exit
    endif
   enddo
  endif

return 
end subroutine point_location





subroutine find_inclusion_vect(npatches,norders,ixyzs,iptype,npts,&
 &srccoefs,srcvals,wts,targ,ntarg,npatches_vect,n_components,&
 &sorted_vector,location_targs,eps)
implicit none
!
!  This function finds the location of a set of points targ relative 
!  to a set of surfaces
!  Same as the subroutine find_inclusion but on a set of targets.
!  therefore the calculation is done by fmm.
!  This function is used by the subroutine evaluate_field_muller to
!  find the ep,mu associated to each point and its scattered field.
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
!    eps - real*8
!     accuracty in the calculation of the double layer
!
!
!  output
!    location_targs integer(ntarg)
!      location of each target point, that is, index of the smaller 
!      interface that contains the target point
!

 !List of calling arguments
 integer, intent(in) :: npatches,npts,n_components,ntarg
 integer, intent(in) :: norders(npatches),npatches_vect(n_components)
 integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
 real ( kind = 8 ), intent(in) :: srcvals(12,npts), srccoefs(9,npts)
 real ( kind = 8 ), intent(in) :: targ(3,ntarg),wts(npts) 
 integer, intent(in) :: sorted_vector(n_components+1)
 integer, intent(out) :: location_targs(ntarg)
 real ( kind = 8 ), intent(in) :: eps

 !List of local variables
 integer count1,count2,icount,n_aux,i1,i2,j1,j2,npatches_aux,npts_aux,x
 integer, allocatable :: ixyzs_aux(:),ipatch_id(:),uvs_targ(:,:)
 real *8, allocatable :: sigma(:),pot(:)
 
 real *8 dpars(2)
 real *8 ttot,done,pi
 logical *8 d(n_components)

 done = 1
 pi = atan(done)*4

 allocate(sigma(npts))
 allocate(ipatch_id(ntarg))
 allocate(uvs_targ(2,ntarg))
 allocate(pot(ntarg))
 ipatch_id(:)=-1
 sigma(:)=-1.0d0
 uvs_targ(1,:)=0
 uvs_targ(2,:)=0
 dpars(1)=0.0d0
 dpars(2)=1.0d0


 i1=1
 n_aux=0
 do count1=1,n_components
  i1=n_aux+1
  n_aux=n_aux+npatches_vect(count1)
  npatches_aux=npatches_vect(count1)
  i2=n_aux
  sigma(ixyzs(i1):(ixyzs(i2+1)-1))=-1.0d0*2.0d0**(count1-1)
 enddo

!
!  This is the double layer evaluator, for a 2D version you should
!  replace that part
!

    call lpcomp_lap_comb_dir(npatches,norders,ixyzs,&
    &iptype,npts,srccoefs,&
    &srcvals,3,ntarg,targ,&
    &ipatch_id,uvs_targ,eps,dpars,&
    &sigma,pot)
!
!
!
  do count1=1,ntarg
    x=nint(pot(count1))
    call digitsinbinary(x,n_components,d)
    call point_location(n_components,sorted_vector,d,&
     &location_targs(count1))
  enddo

 deallocate(sigma)
 deallocate(ipatch_id)
 deallocate(uvs_targ)
 deallocate(pot)
return
end subroutine find_inclusion_vect


subroutine digitsinbinary(x,n,d)
implicit none
!
!  This function provides the n least significant bits of the integer x
!  
!  input:
!    x - integer
!      input number to move to binary
!
!    n - integer
!      number of bits (from the least significant) desired  
!
!  output
!    d - logical(n)
!      n least significant bits of the integer x
!

 !List of calling arguments
 integer, intent(in) :: x,n
 logical *8, intent(out) :: d(n)

 !List of local variables
 integer count1,count2,dd,rr,xx

 xx=x
 dd=2
 d(:)=.false.
 do count1=1,n
  rr=mod(xx,dd)
  if (rr.eq.1) then
    d(count1)=.true.
    else
    d(count1)=.false.
  endif
  xx=xx/2
 enddo

return
end subroutine digitsinbinary






subroutine lpcomp_lap_comb_dir_nonfmm(npatches,norders,ixyzs,&
 &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
 &ipatch_id,uvs_targ,eps,dpars,sigma,pot)
!
!
!------------------------------
!  This subroutine evaluates the layer potential for the representation 
!
!
!  .. math ::
!  
!      u = (\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
!
!  Note: For targets on the boundary, this routine only computes
!  the principal value part, the identity term corresponding to the jump
!  in the layer potential is not included in the layer potential.
!
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array where information for patch i begins
!    - iptype: integer(npatches)
!        type of patch
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefficients of x, $\partial_{u} x$,
!        and $\partial_{v} x$. 
!    - srcvals: double precision (12,npts)
!        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
!        discretization nodes
!    - ndtarg: integer
!        leading dimension of target array
!    - ntarg: integer
!        number of targets
!    - targs: double precision (ndtarg,ntarg)
!        target information
!    - ipatch_id: integer(ntarg)
!        id of patch of target i, id = -1, if target is off-surface
!    - uvs_targ: double precision (2,ntarg)
!        local uv coordinates on patch if on surface, otherwise
!        set to 0 by default
!    - eps: double precision
!        precision requested
!    - dpars: double complex (2)
!        kernel parameters (Referring to formula above)
!        dpars(1) = $\alpha$
!        dpars(2) = $\beta$
!     - sigma: double precision(npts)
!         density for layer potential
!
!  Output arguments
!    - pot: double precision(ntarg)
!        layer potential evaluated at the target points
!
!-----------------------------------
!
      implicit none
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      real *8, intent(in) :: dpars(2)
      real *8, intent(in) :: sigma(npts)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)

      real *8, intent(out) :: pot(ntarg)


      integer nptso,nnz,nquad


      integer nover,npolso
      integer norder,npols
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      complex *16 zpars
      integer ipars
      real *8 timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

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
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,&
       &ntarg,nnz)

      allocate(row_ptr(ntarg+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,&
       &row_ptr,col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind,&
       &iquad)

      ikerorder = -1
      if(abs(dpars(2)).gt.1.0d-16) ikerorder = 0


!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches),ixyzso(npatches+1))

      zpars = 0
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,&
       &rads,npts,srccoefs,ndtarg,ntarg,targs,ikerorder,zpars,&
       &nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts,&
       &srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,&
       &srcover,wover)


!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad))
      
!C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(i) = 0
      enddo
!C$OMP END PARALLEL DO    


      iquadtype = 1

      call getnearquad_lap_comb_dir(npatches,norders,&
       &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       &ipatch_id,uvs_targ,eps,dpars,iquadtype,nnz,row_ptr,col_ind,&
       &iquad,rfac0,nquad,wnear)


!
!
!   compute layer potential
!
      call lpcomp_lap_comb_dir_addsub_nonfmm(npatches,norders,ixyzs,&
       &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       &eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,&
       &sigma,novers,npts_over,ixyzso,srcover,wover,pot)

return
end subroutine lpcomp_lap_comb_dir_nonfmm

subroutine lpcomp_lap_comb_dir_addsub_nonfmm(npatches,norders,ixyzs,&
 &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
 &eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,novers,&
 &nptso,ixyzso,srcover,whtsover,pot)
!
!
!      this subroutine evaluates the layer potential for
!      the representation u = (\alpha S_{0} + \beta D_{0}) 
!      where the near field is precomputed and stored
!      in the row sparse compressed format.
!
!     The fmm is used to accelerate the far-field and 
!     near-field interactions are handled via precomputed quadrature
!
!
!     Using add and subtract - no need to call tree and set fmm parameters
!      can directly call existing fmm library
!
!
!       input:
!         npatches - integer
!            number of patches
!
!         norders- integer(npatches)
!            order of discretization on each patch 
!
!         ixyzs - integer(npatches+1)
!            ixyzs(i) denotes the starting location in srccoefs,
!               and srcvals array corresponding to patch i
!   
!         iptype - integer(npatches)
!            type of patch
!             iptype = 1, triangular patch discretized using RV nodes
!
!         npts - integer
!            total number of discretization points on the boundary
! 
!         srccoefs - real *8 (9,npts)
!            koornwinder expansion coefficients of xyz, dxyz/du,
!            and dxyz/dv on each patch. 
!            For each point srccoefs(1:3,i) is xyz info
!                           srccoefs(4:6,i) is dxyz/du info
!                           srccoefs(7:9,i) is dxyz/dv info
!
!         srcvals - real *8 (12,npts)
!             xyz(u,v) and derivative info sampled at the 
!             discretization nodes on the surface
!             srcvals(1:3,i) - xyz info
!             srcvals(4:6,i) - dxyz/du info
!             srcvals(7:9,i) - dxyz/dv info
!             srcvals(10:12,i) - normals info
! 
!         ndtarg - integer
!            leading dimension of target array
!        
!         ntarg - integer
!            number of targets
!
!         targs - real *8 (ndtarg,ntarg)
!            target information
!
!          eps - real *8
!             precision requested
!
!          dpars - real *8 (2)
!              kernel parameters (Referring to formula (1))
!              dpars(1) = alpha
!              dpars(2) = beta
!
!           nnz - integer *8
!             number of source patch-> target interactions in the near field
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
!           nquad - integer
!               number of entries in wnear
!
!           wnear - real *8(nquad)
!               the near field quadrature correction
!
!           sigma - real *8(npts)
!               density for layer potential
!
!           novers - integer(npatches)
!              order of discretization for oversampled sources and
!               density
!
!         ixyzso - integer(npatches+1)
!            ixyzso(i) denotes the starting location in srcover,
!               corresponding to patch i
!   
!           nptso - integer
!              total number of oversampled points
!
!           srcover - real *8 (12,nptso)
!              oversampled set of source information
!
!           whtsover - real *8 (nptso)
!             smooth quadrature weights at oversampled nodes
!
!
!         output
!           pot - real *8(npts)
!              layer potential evaluated at the target points
!
!           
!               
!
      implicit none
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: ixyzso(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      real *8, intent(in) :: dpars(2)
      integer, intent(in) :: nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: wnear(nquad),sigma(npts)
      integer, intent(in) :: novers(npatches+1)
      integer, intent(in) :: nptso
      real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
      real *8, intent(out) :: pot(ntarg)

      integer norder,npols,nover,npolso
      real *8, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:)
      real *8, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer ns,nt
      real *8 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 tmp(10),val
      real *8 over4pi

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart


      integer ifaddsub

      integer ntj
      
      real *8 ddot,pottmp
      real *8, allocatable :: ctmp2(:),dtmp2(:,:)
      real *8 radexp,epsfmm

      integer ipars,nmax
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      integer nss,ii,l,npover

      integer nd,ntarg0
      integer ier,iper

      integer count1,count2

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)
      data over4pi/0.07957747154594767d0/

      ns = nptso
      done = 1
      pi = atan(done)*4

!
!    estimate max number of sources in neear field of 
!    any target
!
      nmax = 0
      call get_near_corr_max(ntarg,row_ptr,nnz,col_ind,npatches,&
       &ixyzso,nmax)
      allocate(srctmp2(3,nmax),ctmp2(nmax),dtmp2(3,nmax))
           
           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(3,ns),targvals(3,ntarg))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(ns))

! 
!       oversample density
!

      call oversample_fun_surf(1,npatches,norders,ixyzs,iptype,&
       &npts,sigma,novers,ixyzso,ns,sigmaover)


!
!       set relevatn parameters for the fmm
!
      alpha = dpars(1)*over4pi
      beta = dpars(2)*over4pi
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)

        charges(i) = sigmaover(i)*whtsover(i)*alpha
        dipvec(1,i) = sigmaover(i)*whtsover(i)*srcover(10,i)*beta
        dipvec(2,i) = sigmaover(i)*whtsover(i)*srcover(11,i)*beta
        dipvec(3,i) = sigmaover(i)*whtsover(i)*srcover(12,i)*beta
      enddo
!C$OMP END PARALLEL DO      

!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntarg
        targvals(1,i) = targs(1,i)
        targvals(2,i) = targs(2,i)
        targvals(3,i) = targs(3,i)
      enddo
!C$OMP END PARALLEL DO      

      ifcharge = 1
      ifdipole = 1

      if(alpha.eq.0) ifcharge = 0
      if(beta.eq.0) ifdipole = 0

      iper = 0
      ier = 0

!
!
!       call the fmm
!

      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      



!
!        compute threshold for ignoring local computation
!
      
      call get_fmm_thresh(3,ns,sources,3,ntarg,targvals,thresh)
      

      do count1=1,ntarg
        pot(count1)=0.0d0
      enddo
      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        call l3ddirectcp(nd,sources,charges,&
        &ns,targvals,ntarg,pot,thresh)
      endif

      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        call l3ddirectdp(nd,sources,dipvec,&
         &ns,targvals,ntarg,pot,thresh)
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        call l3ddirectcdp(nd,sources,charges,dipvec,&
          &ns,targvals,ntarg,pot,thresh)
      endif

!      call lfmm3d(nd,eps,ns,sources,ifcharge,charges,&
!       &ifdipole,dipvec,iper,ifpgh,tmp,tmp,tmp,ntarg,targvals,ifpghtarg,&
!       &pot,tmp,tmp,ier)
      call cpu_time(t2)
!C$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1



!
!
!       add in precomputed quadrature
!

      call cpu_time(t1)
!C$      t1 = omp_get_wtime()

!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
!C$OMP$PRIVATE(jstart,pottmp,npols,l)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + wnear(jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!C$OMP END PARALLEL DO


!

!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,ii,val,npover)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss + 1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            if(ifcharge.eq.1) ctmp2(nss) = charges(l)
            if(ifdipole.eq.1) then
              dtmp2(1,nss) = dipvec(1,l)
              dtmp2(2,nss) = dipvec(2,l)
              dtmp2(3,nss) = dipvec(3,l)
            endif
          enddo
        enddo

        val = 0
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          call l3ddirectcp(nd,srctmp2,ctmp2,&
           &nss,targvals(1,i),ntarg0,val,thresh)
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          call l3ddirectdp(nd,srctmp2,dtmp2,&
           &nss,targvals(1,i),ntarg0,val,thresh)
        endif

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          call l3ddirectcdp(nd,srctmp2,ctmp2,dtmp2,&
           &nss,targvals(1,i),ntarg0,val,thresh)
        endif
        pot(i) = pot(i) - val
      enddo
      
      call cpu_time(t2)
!C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1


!cc      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)

      
return
end subroutine lpcomp_lap_comb_dir_addsub_nonfmm
