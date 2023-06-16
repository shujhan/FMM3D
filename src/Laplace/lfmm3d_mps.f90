!
!--------------------------------------------------------------------
!
! A fast multi-particle scattering code, based on the code in
! hfmm3d.f of the Flatiron Institute FMM3D library.
!
! Original skeleton code by Manas Rachh, Leslie Greengard, etc.
! FMPS re-write by Mike O'Neil, 2019
! oneil@cims.nyu.edu
!
! The input is assumed to be a collection of multipole expansions,
! each of arbitrary order, and the output is a collection of local
! expansions at the same locations (and of the same order) which take
! into account the potentials from all other multipole expansions.
!
! It is assume that all multipole expansions are well-separated, so
! that even those in LIST 1 can be translated.
!
! We use exp(ikr)/r for the Green's function., without the 1/4\pi
! scaling.
!
!--------------------------------------------------------------------
!

subroutine lfmm3d_mps(nd, eps, nmpole, cmpole, rmpole, mterms, &
    mpole, impole, local, ier)
  !, &
  !  ntarg, targ
  !-----------------------------------------------------------------------
  !   INPUT PARAMETERS:
  !
  !   nd:    in: integer
  !             number of densities
  !   
  !   eps:   in: double precision
  !             requested precision
  !              
  !
  !   nsource in: integer  
  !                number of sources
  !
  !   source  in: double precision (3,nsource)
  !                source(k,j) is the kth component of the jth
  !                source locations
  !
  !
  !     nmpole:  in: integer
  !              number of multipole expansion centers
  !
  !     cmpole:  in: double precision (3,nmpole)
  !              multipole expansion centers
  !
  !     rmpole:  in: double precision (nmpole)
  !              scaling factors for each multipole expansion
  !
  !     mterms:  in: integer (nmpole)
  !              order of the multipole expansions, each expansion
  !              can be of a different order
  !
  !     mpole:   in: double precision (nd,*)
  !              coefficients in the multipole expansions
  !
  !     impole:  in: integer (nmpole)
  !              indexing array for mpole, the ith expansion is at
  !              location mpole(1,impole(i)) and is of order mterms(i)
  !
  !     OUTPUT parameters:
  !
  !     local:   out: double precision ()
  !              local expansions at each center, due to all incoming
  !              multipole expansions (self is ignored). The orders
  !              are the same as for the incoming mpole.
  !     ier:     out: integer
  !              error flag
  !              ier = 0, successful execution
  !              ier = 4, failed to allocate multipole/local expansion
  !              ier = 8, failed to allocate plane wave expansion
  !
  !------------------------------------------------------------------
  
  implicit none

  integer nd


  double precision eps

  integer :: nmpole, mterms(nmpole), impole(nmpole)
  double precision :: cmpole(3,nmpole), rmpole(nmpole)
  double complex :: mpole(*)
  double complex :: local(*)
  

  ! Tree variables
  integer idivflag,ndiv,isep,nboxes,nbmax,nlevels
  integer *8 ltree
  integer nlmax,nlmin,iper,ifunif
  integer ntarg  
  integer *8 ipointer(8)
  integer, allocatable :: itree(:)
  double precision :: targ(3)
  double precision, allocatable :: treecenters(:,:),boxsize(:)
  integer, allocatable :: isrcse(:,:),itargse(:,:),iexpcse(:,:)
  integer, allocatable :: isrc(:)
  integer iexpc,itarg

  !
  ! temporary sorted arrays
  !
  integer :: lmpole, mt, ilen

  integer, allocatable :: mtermssort(:), impolesort(:)
  double precision, allocatable :: cmpolesort(:,:)
  double precision, allocatable :: rmpolesort(:)
  double complex, allocatable :: mpolesort(:)
  double complex, allocatable :: localsort(:)
  
  !
  !  temporary fmm arrays
  !
  double precision epsfmm
  integer, allocatable :: nterms(:)
  integer *8, allocatable :: iaddr(:,:)
  double precision, allocatable :: scales(:)
  double precision, allocatable :: rmlexp(:)

  integer lmptemp,nmax
  integer *8 lmptot
  double precision, allocatable :: mptemp(:),mptemp2(:)

  !
  !       temporary variables not used in particle code
  !
  double precision expc(3),scjsort(1),radexp
  double complex texpssort(100)
  double precision expcsort(3),radssort(1)
  integer ntj,nexpc,nadd, npts, perm, ptr, ifunsort

  !
  !        other temporary variables
  !
  integer :: i, j, l, ijk, iert,ifprint,ilev,idim,ier
  integer :: nlege, lw7, lused7
  double precision :: wlege(40000)
  double precision time1,time2,omp_get_wtime,second


  !
  !
  ! ifprint is an internal information printing flag.  Suppressed if
  ! ifprint=0.  Prints timing breakdown and other things if ifprint=1.
  !      
  ifprint=1
  print *, 'ndiv still needs to be optimized'
  ndiv = 1

  if(ifprint.ge.1) print *, "ndiv =",ndiv



  nexpc = 0
  radexp = 0
  nadd = 0
  ntj = 0

  !
  ! set tree flags
  !
  isep = 1
  nlmax = 200
  nlevels = 0
  nboxes = 0
  ltree = 0

  idivflag = 0




  !
  ! memory management code for constructing level restricted tree
  !
  iert = 0

  ntarg = 0
  targ(1) = 0
  targ(2) = 0
  targ(3) = 0
!
!!      set tree flags
! 
  nlmax = 51
  nlevels = 0
  nboxes = 0
  ltree = 0
  nlmin = 0
  ifunif = 0
  iper = 0

!
!!     memory management code for contructing level restricted tree
  call pts_tree_mem(cmpole,nmpole,targ,ntarg,idivflag,ndiv, &
     nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree)
      

  allocate(itree(ltree))
  allocate(boxsize(0:nlevels))
  allocate(treecenters(3,nboxes))

  call pts_tree_build(cmpole,nmpole,targ,ntarg,idivflag,ndiv, &
     nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree,itree,ipointer, &
     treecenters,boxsize)
      

  allocate(isrcse(2,nboxes),itargse(2,nboxes),iexpcse(2,nboxes))
  allocate(isrc(nmpole))

  call pts_tree_sort(nmpole,cmpole,itree,ltree,nboxes,nlevels, &
    ipointer,treecenters,isrc,isrcse)


!
!   End of tree build
!

  if(ifprint.ge.1) print *, "nlevels=",nlevels
  if(ifprint.ge.1) print *, ltree/1.0d9


  !
  !     Allocate sorted source and target arrays      
  !
  
  allocate(scales(0:nlevels),nterms(0:nlevels))
  do ilev = 0,nlevels
    scales(ilev) = boxsize(ilev)
  enddo

  !
  !compute length of expansions at each level      
  !
  nmax = 0
  do i=0,nlevels
    call l3dterms(eps,nterms(i))
    if(nterms(i).gt.nmax) nmax = nterms(i)
  enddo


  !       
  ! Multipole and local expansions will be held in workspace in
  ! locations pointed to by array iaddr(2,nboxes).
  !
  ! iiaddr is pointer to iaddr array, itself contained in workspace.
  ! imptemp is pointer for single expansion (dimensioned by nmax)
  !
  ! ... allocate iaddr and temporary arrays
  !

  allocate(iaddr(2,nboxes))
  lmptemp = (nmax+1)*(2*nmax+1)*2*nd
  allocate(mptemp(lmptemp),mptemp2(lmptemp))


  !
  ! reorder multipole expansions, their centers, and rscales
  !
  allocate(cmpolesort(3,nmpole))
  allocate(rmpolesort(nmpole))
  allocate(impolesort(nmpole))
  allocate(mtermssort(nmpole))

  lmpole = 0
  do i = 1,nmpole
    lmpole = lmpole + (mterms(i)+1)*(2*mterms(i)+1)
  end do
  lmpole = nd*lmpole
  
  allocate(mpolesort(lmpole) )

  call dreorderf(3, nmpole, cmpole, cmpolesort, isrc)
  call dreorderf(1, nmpole, rmpole, rmpolesort, isrc)
  call ireorderf(1, nmpole, mterms, mtermssort, isrc)


  impolesort(1) = 1  
  do i = 1,nmpole

    mt = mtermssort(i)
    ilen = (mt+1)*(2*mt+1)

    ijk = 1
    do j = 1,ilen
      do l = 1,nd
        mpolesort(impolesort(i)+ijk-1) = &
            mpole(impole(isrc(i))+ijk-1)
        ijk = ijk + 1
      end do
    end do

    if (i .lt. nmpole) impolesort(i+1) = impolesort(i) + nd*ilen
   
  end do


  !
  ! allocate memory need by multipole, local expansions at all
  ! levels
  !
  ! irmlexp is pointer for workspace need by various fmm routines
  !
  call mpalloc(nd,itree(ipointer(1)),iaddr,nlevels,lmptot,nterms)
  if(ifprint.ge. 1) print *, "lmptot =",lmptot/1.0d9
  ier = 0
  allocate(rmlexp(lmptot),stat=iert)
  if(iert.ne.0) then
    print *, "Cannot allocate mpole expansion workspace"
    print *, "lmptot=", lmptot
    ier = 4
    return
  endif

  allocate( localsort(lmpole) )
  

  !
  ! Memory allocation is complete. 
  ! Call main fmm routine
  !
  call cpu_time(time1)
  !$ time1=omp_get_wtime()
  call lfmm3dmain_mps(nd, eps, &
      nmpole, cmpolesort, rmpolesort, mtermssort, mpolesort, &
      impolesort, localsort, &
      iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp, &
      itree,ltree,ipointer,ndiv,nlevels, &
      nboxes,iper,boxsize, treecenters, isrcse, &
      scales,itree(ipointer(1)),nterms,ier )
  if(ier.ne.0) return

  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  if( ifprint .eq. 1 ) call prin2('time in fmm main=*', &
      time2-time1,1)

  !
  ! now unsort the local expansions
  !
  do i = 1,nmpole

    mt = mtermssort(i)
    ilen = (mt+1)*(2*mt+1)

    ijk = 1
    do j = 1,ilen
      do l = 1,nd
        local(impole(isrc(i))+ijk-1) = &
            localsort(impolesort(i)+ijk-1) 
        ijk = ijk + 1
      end do
    end do
  end do


  return
end subroutine lfmm3d_mps





subroutine lfmm3dmain_mps(nd, eps, &
    nmpole, cmpolesort, rmpolesort, mtermssort, mpolesort, &
    impolesort, localsort, &
    iaddr, rmlexp, lmptot, mptemp, mptemp2, lmptemp, &
    itree, ltree, ipointer, ndiv, nlevels, &
    nboxes, iper, boxsize, centers, isrcse, &
    rscales, laddr, nterms, ier )
  implicit none

  !
  ! INPUT variables
  !
  integer :: nd, ndiv,nlevels
  double precision :: eps
  !double complex :: zk,zk2

  ! input multipole stuff
  integer :: nmpole, mtermssort(nmpole)
  double precision :: cmpolesort(3,nmpole), rmpolesort(nmpole)
  double complex :: mpolesort(*)
  integer :: impolesort(nmpole)

  ! storage stuff for tree and multipole expansions
  integer :: lmptemp
  integer *8 :: iaddr(2,nboxes), lmptot
  double precision :: rmlexp(lmptot)
  double precision :: mptemp(lmptemp)
  double precision :: mptemp2(lmptemp)

  ! tree variables
  integer :: isep,iper
  integer *8 :: ltree
  integer :: laddr(2,0:nlevels)
  integer :: nterms(0:nlevels)
  integer *8 :: ipointer(8)
  integer :: itree(ltree)
  integer :: nboxes
  integer :: mnbors,mnlist1, mnlist2,mnlist3,mnlist4
  integer :: isrcse(2,nmpole)
  integer, allocatable :: nlist1(:),list1(:,:)
  integer, allocatable :: nlist2(:),list2(:,:)
  integer, allocatable :: nlist3(:),list3(:,:)
  integer, allocatable :: nlist4(:),list4(:,:)
  double precision :: rscales(0:nlevels)
  double precision :: boxsize(0:nlevels)
  double precision :: centers(3,nboxes)

  !
  ! OUTPUT variables
  !
  double complex :: localsort(*)


  !
  ! LOCAL variables
  !

  
  integer nuall,ndall,nnall,nsall,neall,nwall
  integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
  integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
  integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8
  integer uall(200),dall(200),nall(120),sall(120),eall(72),wall(72)
  integer u1234(36),d5678(36),n1256(24),s3478(24)
  integer e1357(16),w2468(16),n12(20),n56(20),s34(20),s78(20)
  integer e13(20),e57(20),w24(20),w68(20)
  integer e1(20),e3(5),e5(5),e7(5),w2(5),w4(5),w6(5),w8(5)



  ! pw stuff
  integer nexpmax, nlams, nmax, nthmax, nphmax,nmax2,nmaxt
  integer lca
  double precision, allocatable :: carray(:,:), dc(:,:)
  double precision, allocatable :: cs(:,:),fact(:),rdplus(:,:,:)
  double precision, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
  double precision, allocatable :: rdmsq3(:,:,:)
  
  double precision, allocatable :: rlams(:),whts(:)

  double precision, allocatable :: rlsc(:,:,:)
  integer, allocatable :: nfourier(:), nphysical(:)
  integer nexptot, nexptotp
  double complex, allocatable :: xshift(:,:)
  double complex, allocatable :: yshift(:,:)
  double precision, allocatable :: zshift(:,:)

  double complex, allocatable :: fexpe(:),fexpo(:),fexpback(:)
  double complex, allocatable :: mexp(:,:,:,:)
  double complex, allocatable :: mexpf1(:,:),mexpf2(:,:)
  double complex, allocatable ::
         mexpp1(:,:),mexpp2(:,:),mexppall(:,:,:)

  double complex, allocatable :: tmp(:,:,:)
  double precision, allocatable :: mptmp(:,:)


  ! temp variables


   integer ix,iy,iz,ictr
      double precision rtmp
      double complex zmul

      integer nlege, lw7, lused7, itype
      double precision wlege(40000)
      integer nterms_eval(4,0:nlevels)
       double complex eye, ztmp





















! list 4 variables
  integer cntlist4
  integer, allocatable :: list4ct(:),ilist4(:),nlist4tmp(:)
      double complex, allocatable :: gboxmexp(:,:,:)
      double complex, allocatable :: gboxwexp(:,:,:,:,:)
      double complex, allocatable :: pgboxwexp(:,:,:,:)
      double precision, allocatable :: gboxsubcenters(:,:,:)
      double precision, allocatable :: gboxsort(:,:,:)
      integer, allocatable :: gboxind(:,:)
      integer, allocatable :: gboxfl(:,:,:)
      double precision, allocatable :: gboxcgsort(:,:,:)
      double precision, allocatable :: gboxdpsort(:,:,:,:)
! end of list 4 variables


!   Hessian variables
!
  double precision, allocatable :: scarray(:,:)

  integer *8 bigint
  integer iert

  double precision d,time1,time2,omp_get_wtime
  double precision timeinfo(10)
  double precision rtmp
  double precision wlege(40000)
  double precision thresh
  double precision alphaj
  double precision radius
  double precision pi
  double precision :: rtmp1,rtmp2,rtmp3,rtmp4, done
  double precision, allocatable :: xnodes(:),wts(:)
  double complex zmul
  double complex eye, ztmp,zmult
  double complex :: ima, cd, cd1(10), cd2(10), work(100000)

  data ima/(0.0d0,1.0d0)/



  mnlist1 = 0
  mnlist2 = 0
  mnlist3 = 0
  mnlist4 = 0
  mnbors = 27

  isep = 1
      
  call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize, &
    centers,itree(ipointer(3)),itree(ipointer(4)), &
    itree(ipointer(5)),isep,itree(ipointer(6)),mnbors, &
    itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      
  allocate(list1(mnlist1,nboxes),nlist1(nboxes))
  allocate(list2(mnlist2,nboxes),nlist2(nboxes))
  allocate(list3(mnlist3,nboxes),nlist3(nboxes))
  allocate(list4(mnlist4,nboxes),nlist4(nboxes))

  call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize, &
    centers,itree(ipointer(3)),itree(ipointer(4)), &
    itree(ipointer(5)),isep,itree(ipointer(6)),mnbors, &
    itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2, &
    mnlist2,list2,nlist3,mnlist3,list3,nlist4,mnlist4,list4)


!     Initialize routines for plane wave mp loc translation
 
      if(isep.eq.1) then
         if(eps.ge.0.5d-2) nlams = 12
         if(eps.lt.0.5d-2.and.eps.ge.0.5d-3) nlams = 12
         if(eps.lt.0.5d-3.and.eps.ge.0.5d-6) nlams = 20
         if(eps.lt.0.5d-6.and.eps.ge.0.5d-9) nlams = 29
         if(eps.lt.0.5d-9) nlams = 37
      endif
      if(isep.eq.2) then
         if(eps.ge.0.5d-3) nlams = 9
         if(eps.lt.0.5d-3.and.eps.ge.0.5d-6) nlams = 15
         if(eps.lt.0.5d-6.and.eps.ge.0.5d-9) nlams = 22
         if(eps.lt.0.5d-9) nlams = 29
      endif

      allocate(rlams(nlams),whts(nlams))
      allocate(nphysical(nlams),nfourier(nlams))



  cntlist4 = 0
  allocate(list4ct(nboxes))
  allocate(ilist4(nboxes),nlist4tmp(nboxes))

  !$omp parallel do default(shared) private(i)
  do i=1,nboxes
    list4ct(i) = 0
    ilist4(i) = 0
    nlist4tmp(i) = 0
  enddo
  !$omp end parallel do

  done = 1
  pi = 4.0d0*atan(done)




  nmax = 0
  do i=0,nlevels
    if(nmax.lt.nterms(i)) nmax = nterms(i)
  enddo
  allocate(rscpow(0:nmax))
  allocate(carray(4*nmax+1,4*nmax+1))
  allocate(rsc(0:nmax))
  allocate(dc(0:4*nmax,0:4*nmax))
  allocate(rdplus(0:nmax,0:nmax,-nmax:nmax))
  allocate(rdminus(0:nmax,0:nmax,-nmax:nmax))
  allocate(rdsq3(0:nmax,0:nmax,-nmax:nmax))
  allocate(rdmsq3(0:nmax,0:nmax,-nmax:nmax))
  allocate(rlsc(0:nmax,0:nmax,nlams))



!     generate rotation matrices and carray
      call getpwrotmat(nmax,carray,rdplus,rdminus,rdsq3,rdmsq3,dc)


!     generate rlams and weights (these are the nodes
!     and weights for the lambda integral)

      call vwts(rlams,whts,nlams)


!     generate the number of fourier modes required to represent the
!     moment function in fourier space

      call numthetahalf(nfourier,nlams)
 
!     generate the number of fourier modes in physical space
!     required for the exponential representation
      call numthetafour(nphysical,nlams)

!     Generate powers of lambda for the exponential basis
      call rlscini(rlsc,nlams,rlams,nmax)


!     nn = 10*(nmax+2)**2
      allocate(scarray(nn,0:nlevels))

      do ilev=0,nlevels
        call l3dtaevalhessdini(nterms(ilev),scarray(1,ilev))
      enddo


!     Compute total number of plane waves
      nexptotp = 0
      nexptot = 0
      nthmax = 0
      nphmax = 0
      nthd = 1
      nn = 0
      do i=1,nlams
         nexptot = nexptot + nfourier(i)
         nexptotp = nexptotp + nphysical(i)
         if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
         if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
         nn = nn + nphysical(i)*nfourier(i)
      enddo

      allocate(fexpe(nn),fexpo(nn),fexpback(nn))
      allocate(tmp(nd,0:nmax,-nmax:nmax))
      allocate(mptmp(lmptemp,nthd))

      allocate(xshift(-5:5,nexptotp))
      allocate(yshift(-5:5,nexptotp))
      allocate(zshift(5,nexptotp))

      allocate(mexpf1(nd,nexptot),mexpf2(nd,nexptot),
        mexpp1(nd,nexptotp))
      allocate(mexpp2(nd,nexptotp),mexppall(nd,nexptotp,16))

!
!!      NOTE: there can be some memory savings here
!
      bigint = 0
      bigint = nboxes
      bigint = bigint*6
      bigint = bigint*nexptotp*nd

      if(ifprint.ge.1) print *, "mexp memory=",bigint/1.0d9

      allocate(mexp(nd,nexptotp,nboxes,6),stat=iert)
      if(iert.ne.0) then
        print *, "Cannot allocate pw expansion workspace"
        print *, "bigint=", bigint
        ier = 8
        return
      endif

      allocate(list4ct(nboxes))
      allocate(ilist4(nboxes))
      do i=1,nboxes
        list4ct(i)=0
        ilist4(i)=0
      enddo
      cntlist4=0

!     Precompute table for shifting exponential coefficients in 
!     physical domain
      call mkexps(rlams,nlams,nphysical,nexptotp,xshift,yshift,zshift)

!     Precompute table of exponentials for mapping from
!     fourier to physical domain
      call mkfexp(nlams,nfourier,nphysical,fexpe,fexpo,fexpback)
      

!!    compute array of factorials

     
      nmax2 = 2*nmax
      allocate(fact(0:nmax2),cs(0:nmax,-nmax:nmax))
      
      d = 1
      fact(0) = d
      do i=1,nmax2
        d=d*sqrt(i+0.0d0)
        fact(i) = d
      enddo

      cs(0,0) = 1.0d0
      do l=1,nmax
        do m=0,l
          cs(l,m) = ((-1)**l)/(fact(l-m)*fact(l+m))
          cs(l,-m) = cs(l,m)
        enddo
      enddo
      
      if(ifprint.ge.1) 
     1   call prin2('end of generating plane wave info*',i,0)



































  !
  ! threshold for computing interactions:
  !
  ! interactions will be ignored for all pairs of sources and targets
  ! which satisfy |r| < thresh where r is the disance between them
  !
  thresh = 2.0d0**(-51)*boxsize(0)

  


















  !
  ! ifprint is an internal information printing flag. 
  ! Suppressed if ifprint=0.
  ! Prints timing breakdown and other things if ifprint=1.
  ! Prints timing breakdown, list information, and other things if ifprint=2.
  !       
  ifprint=1
  do i=1,6
    timeinfo(i)=0
  enddo






  !
  ! ... set all multipole and local expansions to zero
  !
  do ilev = 0,nlevels
    !$omp parallel do default(shared) private(ibox)
    do ibox = laddr(1,ilev),laddr(2,ilev)
      call mpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
      call mpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
    enddo
    !$omp end parallel do          
  enddo

  ltot = 0
  !$omp parallel do default(shared) private(l,mt)  &
  !$omp reduction(+:ltot)
  do l = 1,nmpole
    mt = mtermssort(l)
    ltot = ltot + (mt+1)*(2*mt+1)
  end do
  !$omp end parallel do

  
  !$omp parallel do default (shared) private(l)
  do l = 1,ltot*nd
      localsort(l) = 0
  end do
  !$omp end parallel do
  


  ! initialize legendre function evaluation routines
  nlege = 100
  lw7 = 40000
  call ylgndrfwini(nlege,wlege,lw7,lused7)


  lca = 4*nmax
      do ilev=1,nlevels-1
         do ibox=laddr(1,ilev),laddr(2,ilev)
            if(nlist3(ibox).gt.0) then
              cntlist4=cntlist4+1
              list4ct(ibox)=cntlist4
              ilist4(cntlist4)=ibox
            endif
         enddo
      enddo
      if(ifprint.ge.1) print *,"nboxes:",nboxes,"cntlist4:",cntlist4
      allocate(pgboxwexp(nd,nexptotp,cntlist4,6))
      allocate(gboxmexp(nd*(nterms(nlevels)+1)*
                        (2*nterms(nlevels)+1),8,cntlist4))



      allocate(gboxsubcenters(3,8,nthd))
      allocate(gboxfl(2,8,nthd))

      nmaxt = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts)
!$OMP$REDUCTION(max:nmaxt)
      do ibox=1,nboxes
        if(list4ct(ibox).gt.0) then
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = iend-istart+1
          if(npts.gt.nmaxt) nmaxt = npts
        endif
      enddo
!$OMP END PARALLEL DO

      allocate(gboxind(nmaxt,nthd))
      allocate(gboxsort(3,nmaxt,nthd))
      allocate(gboxwexp(nd,nexptotp,6,8,nthd))
      allocate(gboxcgsort(nd,nmaxt,nthd))
      allocate(gboxdpsort(nd,3,nmaxt,nthd))

!   note gboxmexp is an array not scalar
      pgboxwexp=0d0
      gboxmexp=0d0

!     form mexp for all list4 type box at first ghost box center
      do ilev=1,nlevels-1

         rscpow(0) = 1.0d0/boxsize(ilev+1)
         rtmp = rscales(ilev+1)/boxsize(ilev+1)
         do i=1,nterms(ilev+1)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo

!$OMP PARALLEL DO DEFAULT(SHARED)
!$OMP$PRIVATE(ibox,istart,iend,jbox,jstart,jend,npts,npts0,i)
!$OMP$PRIVATE(ithd)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
!$          ithd=omp_get_thread_num()
            ithd = ithd + 1
            if(list4ct(ibox).gt.0) then
              istart=isrcse(1,ibox)
              iend=isrcse(2,ibox)
              npts = iend-istart+1

              if(npts.gt.0) then
                call subdividebox(sourcesort(1,istart),npts,
                    centers(1,ibox),boxsize(ilev+1),
                    gboxind(1,ithd),gboxfl(1,1,ithd),
                    gboxsubcenters(1,1,ithd))
                call dreorderf(3,npts,sourcesort(1,istart),
                    gboxsort(1,1,ithd),gboxind(1,ithd))
                if(ifcharge.eq.1) then
                  call dreorderf(nd,npts,chargesort(1,istart),
                      gboxcgsort(1,1,ithd),gboxind(1,ithd))
                endif
                if(ifdipole.eq.1) then
                  call dreorderf(3*nd,npts,dipvecsort(1,1,istart),
                      gboxdpsort(1,1,1,ithd),gboxind(1,ithd))
                endif
                do i=1,8
                  if(gboxfl(1,i,ithd).gt.0) then
                    jstart=gboxfl(1,i,ithd)
                    jend=gboxfl(2,i,ithd)
                    npts0=jend-jstart+1
                    jbox=list4ct(ibox)
                    if(ifcharge.eq.1.and.ifdipole.eq.0) then
                      call l3dformmpc(nd,rscales(ilev+1),
                        gboxsort(1,jstart,ithd),
                        gboxcgsort(1,jstart,ithd),
                        npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1),
                        gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.0.and.ifdipole.eq.1) then
                      call l3dformmpd(nd,rscales(ilev+1),
                        gboxsort(1,jstart,ithd),
                        gboxdpsort(1,1,jstart,ithd),
                        npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1),
                        gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.1.and.ifdipole.eq.1) then
                      call l3dformmpcd(nd,rscales(ilev+1),
                        gboxsort(1,jstart,ithd),
                        gboxcgsort(1,jstart,ithd),
                        gboxdpsort(1,1,jstart,ithd),
                        npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1),
                        gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    call l3dmpmp(nd,rscales(ilev+1),
                        gboxsubcenters(1,i,ithd),gboxmexp(1,i,jbox),
                        nterms(ilev+1),rscales(ilev),centers(1,ibox),
                        rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
     
                    call mpscale(nd,nterms(ilev+1),gboxmexp(1,i,jbox),
                        rscpow,tmp(1,0,-nmax,ithd))
!
!!               process up down for current box
!
                    call mpoletoexp(nd,tmp(1,0,-nmax,ithd),
                        nterms(ilev+1),nlams,
                        nfourier,nexptot,mexpf1(1,1,ithd),
                        mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd),
                        nlams,rlams,nfourier,
                        nphysical,nthmax,gboxwexp(1,1,1,i,ithd),
                        fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd),
                        nlams,rlams,nfourier,
                        nphysical,nthmax,gboxwexp(1,1,2,i,ithd),
                        fexpe,fexpo)

                    call processgboxudexp(nd,gboxwexp(1,1,1,i,ithd),
                        gboxwexp(1,1,2,i,ithd),i,nexptotp,
                        pgboxwexp(1,1,jbox,1),pgboxwexp(1,1,jbox,2),
                        xshift,yshift,zshift)
!
!!                process north-south for current box
!
                    call rotztoy(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd),
                        mptmp(1,ithd),rdminus)
                    call mpoletoexp(nd,mptmp(1,ithd),
                        nterms(ilev+1),nlams,
                        nfourier,nexptot,mexpf1(1,1,ithd),
                        mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd),
                        nlams,rlams,nfourier,
                        nphysical,nthmax,gboxwexp(1,1,3,i,ithd),
                        fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd),
                        nlams,rlams,nfourier,
                        nphysical,nthmax,gboxwexp(1,1,4,i,ithd),
                        fexpe,fexpo)

                    call processgboxnsexp(nd,gboxwexp(1,1,3,i,ithd),
                        gboxwexp(1,1,4,i,ithd),i,nexptotp,
                        pgboxwexp(1,1,jbox,3),pgboxwexp(1,1,jbox,4),
                        xshift,yshift,zshift)
!
!!               process east-west for current box
!
                    call rotztox(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd),
                        mptmp(1,ithd),rdplus)
                    call mpoletoexp(nd,mptmp(1,ithd),
                        nterms(ilev+1),nlams,
                        nfourier,nexptot,mexpf1(1,1,ithd),
                        mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd),
                        nlams,rlams,nfourier,
                        nphysical,nthmax,gboxwexp(1,1,5,i,ithd),
                        fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd),
                        nlams,rlams,nfourier,
                        nphysical,nthmax,gboxwexp(1,1,6,i,ithd),
                        fexpe,fexpo)
                
                    call processgboxewexp(nd,gboxwexp(1,1,5,i,ithd),
                        gboxwexp(1,1,6,i,ithd),i,nexptotp,
                        pgboxwexp(1,1,jbox,5),pgboxwexp(1,1,jbox,6),
                        xshift,yshift,zshift)
                  endif
                enddo
              endif
            endif
         enddo
!$OMP END PARALLEL DO
      enddo
      deallocate(gboxfl,gboxsubcenters,gboxwexp,gboxcgsort)
      deallocate(gboxdpsort,gboxind,gboxsort)

      call cpu_time(time2)
!$    time2=omp_get_wtime()
      if(ifprint.ge.1) print *,"mexp list4 time:",time2-time1
      timeinfo(3)=time2-time1
!     end of count number of boxes are in list4

















  !
  ! ----- Step 1: Shift incoming multipole expansions to the center
  ! of each leaf-node box -----
  !
  if(ifprint .ge. 1) call prinf('=== STEP 1 (shift mp) ====*',i,0)

  call cpu_time(time1)
  !$ time1=omp_get_wtime()

  do ilev=2,nlevels

    nquad2 = nterms(ilev)*2.5
    nquad2 = max(6,nquad2)
    ifinit2 = 1
    call legewhts(nquad2,xnodes,wts,ifinit2)

    !!!!!!!!!
    ! this radius has a fudge factor in it, debug in future
    !!!!!!!!!
    radius = boxsize(ilev)/2*sqrt(3.0d0)*1.5d0
    
      !$omp parallel do default(shared) &
      !$omp   private(ibox,npts,istart,iend,nchild,i)
      do ibox = laddr(1,ilev),laddr(2,ilev)

        istart = isrcse(1,ibox) 
        iend = isrcse(2,ibox) 
        npts = iend-istart+1
        
        nchild = itree(ipointer(4)+ibox-1)

        if((npts.gt.0) .and. (nchild.eq.0)) then
          
          do i = istart,iend
           ! call h3dmpmp(nd, zk, rmpolesort(i), cmpolesort(1,i), &
            !    mpolesort(impolesort(i)), mtermssort(i), &
            !    rscales(ilev), centers(1,ibox), &
            !    rmlexp(iaddr(1,ibox)), nterms(ilev), &
            !    radius, xnodes, wts, nquad2)

            

            call l3dmpmp(nd, rmpolesort(i), cmpolesort(1,i), &
                mpolesort(impolesort(i)), mtermssort(i), &
                rscales(ilev), centers(1,ibox), &
                rmlexp(iaddr(1,ibox)), nterms(ilev), &
                dc,lca)


          end do
        
        endif
      enddo
      !$omp end parallel do            

  enddo

  call cpu_time(time2)
  !$ time2 = omp_get_wtime()
  timeinfo(1)=time2-time1



  !
  ! ----- Step 2: List 3 interactions, non-adjacent nearfield of all
  ! boxes, do these translations directly right now
  !

  if(ifprint.ge.1) &
      call prinf('=== STEP 2 (form lo) ===*',i,0)
  call cpu_time(time1)
  !$ time1=omp_get_wtime()

  
  do ilev=2,nlevels

    nquad2 = nterms(ilev)*2.5
    nquad2 = max(6,nquad2)
    ifinit2 = 1
    call legewhts(nquad2,xnodes,wts,ifinit2)
    radius = boxsize(ilev)/2*sqrt(3.0d0)

    !$omp parallel do default(shared) &
    !$omp   private(ibox,jbox,istart,iend,npts,i,j) &
    !$omp   schedule(dynamic)
    do ibox=laddr(1,ilev),laddr(2,ilev)
      do i=1,nlist4(ibox)
        jbox = list4(i,ibox) 

        !
        ! Form local expansion for all boxes in list3 of the current
        ! box
        !
        istart = isrcse(1,jbox) 
        iend = isrcse(2,jbox) 
        npts = iend-istart+1
        if(npts.gt.0) then
          do j = istart,iend
            !call h3dmploc(nd, zk, rmpolesort(j), cmpolesort(1,j), &
            !    mpolesort(impolesort(j)), mtermssort(j), &
             !   rscales(ilev), centers(1,ibox), &
             !   rmlexp(iaddr(2,ibox)), nterms(ilev), &
             !   radius, xnodes, wts, nquad2)

            call l3dmploc(nd, rmpolesort(j), cmpolesort(1,j), &
                mpolesort(impolesort(j)), mtermssort(j), &
                rscales(ilev), centers(1,ibox), &
                rmlexp(iaddr(2,ibox)), nterms(ilev), &
                dc,lca)


          end do
        endif
      enddo
    enddo
    !$omp end parallel do
  enddo

  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  timeinfo(2)=time2-time1

  !
  ! A further optimization could be done, which translations initial
  ! multipole expansions directly to the center of the parent box,
  ! eliminating Step 1
  !


  !
  ! Step 3: Upward pass, multipole-to-multipole merges
  !
  
  if(ifprint .ge. 1) call prinf('=== STEP 3 (merge mp) ====*',i,0)

  call cpu_time(time1)
  !$ time1=omp_get_wtime()

  do ilev=nlevels-1,0,-1
    nquad2 = nterms(ilev)*2.5
    nquad2 = max(6,nquad2)
    ifinit2 = 1
    call legewhts(nquad2,xnodes,wts,ifinit2)
    radius = boxsize(ilev)/2*sqrt(3.0d0)

    !$omp parallel do default(shared) &
    !$omp     private(ibox,i,jbox,istart,iend,npts)
    do ibox = laddr(1,ilev),laddr(2,ilev)
      do i=1,8
        jbox = itree(ipointer(5)+8*(ibox-1)+i-1)

        if(jbox.gt.0) then
          istart = isrcse(1,jbox) 
          iend = isrcse(2,jbox) 
          npts = iend-istart+1

          if(npts .gt. 0) then
          !  call h3dmpmp(nd,zk,rscales(ilev+1), &
          !      centers(1,jbox),rmlexp(iaddr(1,jbox)), &
          !      nterms(ilev+1),rscales(ilev),centers(1,ibox), &
          !      rmlexp(iaddr(1,ibox)),nterms(ilev), &
          !      radius,xnodes,wts,nquad2)



            call l3dmpmp(nd,rscales(ilev+1), &
                centers(1,jbox),rmlexp(iaddr(1,jbox)), &
                nterms(ilev+1),rscales(ilev),centers(1,ibox), &
                rmlexp(iaddr(1,ibox)),nterms(ilev), &
                dc,lca)

          endif
        endif
      enddo
    enddo
    !$omp end parallel do          
  enddo

  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  timeinfo(3)=time2-time1

  !
  ! ----- Step 4: Crossward pass, multipole-to-local -----
  ! Note: This is generally the most expensive component of the FMM
  !
  if (ifprint .ge. 1) &
      call prinf('=== Step 4 (mp to loc) ===*',i,0)

  call cpu_time(time1)
  !$ time1=omp_get_wtime()


  !
  ! Set the below flag depending on whether you want planewave
  ! translations or standard point-and-shoot translations
  !
  ifpw = 1
  ifmp = 0

  if ( (ifprint .ge. 1) .and. (ifmp .eq. 1) ) &
      call prinf('=== doing point and shoots... ===*',i,0)

  if ( (ifprint .ge. 1) .and. (ifmp .eq. 0) ) &
      call prinf('=== doing plane waves... ===*',i,0)

  
  if (ifpw .eq. 1) then

    !if (ifprint .ge. 1) &
    !    call prinf('=== doing plane waves... ===*',i,0)
    
    do ilev = 2,nlevels
        rscpow(0) = 1.0d0/boxsize(ilev)
        rtmp = rscales(ilev)/boxsize(ilev)
        do i=1,nterms(ilev)
          rscpow(i) = rscpow(i-1)*rtmp
        enddo


 



































































































        !
        ! create multipole to plane wave expansion for all boxes at this
        ! level
        !
        !$omp parallel do default (shared) &
        !$omp    private(ibox,istart,iend,npts,tmp,mexpf1,mexpf2,tmp2)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          istart = isrcse(1,ibox) 
          iend = isrcse(2,ibox) 
          npts = iend - istart+1
          if(npts.gt.0) then

            ! rescale multipole expansion
            call mpscale(nd,nterms(ilev),rmlexp(iaddr(1,ibox)), &
                rscpow,tmp)

            call mpoletoexp(nd,tmp,nterms(ilev), &
                nlams,nfourier,nexptot,mexpf1,mexpf2,rlsc) 

            call ftophys(nd,mexpf1,nlams,rlams,nfourier,nphysical,nthmax &
                mexp(1,1,ibox,1),fexp,fexpo)           

            call ftophys(nd,mexpf2,nlams,rlams,nfourier,nphysical, nphmax&
                mexp(1,1,ibox,2),fexp,fexpo)


            ! form mexpnorth, mexpsouth for current box

            ! Rotate mpole for computing mexpnorth and
            ! mexpsouth
            call rotztoy(nd,nterms(ilev),tmp, &
                mptmp(1,ithd),rdminus)

            call mpoletoexp(nd,mptmp(1,ithd),nterms(ilev),nlams, &
                nfourier,nexptot,mexpf1,mexpf2,rlsc)

            call ftophys(nd,mexpf1,nlams,rlams，nfourier, &
                nphysical,nthmax,mexp(1,1,ibox,3),fexpe,fexpo)           

            call ftophys(nd,mexpf2,nlams,,rlams,nfourier, &
                nphysical,nthmax,mexp(1,1,ibox,4),fexpe,fexpo)   


            ! Rotate mpole for computing mexpeast, mexpwest
            call rotztox(nd,nterms(ilev),tmp, &
                mptmp(1,ithd),rdplus)
            call mpoletoexp(nd,mptmp(1,ithd),nterms(ilev),nlams, &
                nfourier,nexptot,mexpf1,mexpf2,rlsc)

            call ftophys(nd,mexpf1,nlams,rlams，nfourier, &
                nphysical,nthmax,mexp(1,1,ibox,5),fexpe,fexpo)

            call ftophys(nd,mexpf2,nlams,rlams，nfourier, &
                nphysical,mexp(1,1,ibox,6),fexpe,fexpo)           

          endif
        enddo
        !$omp end parallel do       

        !
        ! Loop over parent boxes and ship plane wave expansions to the
        ! first child of parent boxes.
        !
        ! The codes are now written from a gathering perspective so the
        ! first child of the parent is the one recieving all the local
        ! expansions coming from all the lists
        !

        !$omp parallel do default (shared) &
        !$omp private(ibox,istart,iend,npts,nchild) &
        !$omp private(mexpf1,mexpf2,mexpp1,mexpp2,mexppall) &
        !$omp private(nuall,uall,ndall,dall,nnall,nall,nsall,sall) &
        !$omp private(neall,eall,nwall,wall,nu1234,u1234,nd5678,d5678)&
        !$omp private(nn1256,n1256,ns3478,s3478,ne1357,e1357,nw2468,w2468)&
        !$omp private(nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,e13,ne57,e57)&
        !$omp private(nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,ne7,e7)&
        !$omp private(nw2,w2,nw4,w4,nw6,w6,nw8,w8)

         rscpow(0) = 1.0d0
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo



        do ibox = laddr(1,ilev-1),laddr(2,ilev-1)

          npts = 0
          nchild = itree(ipointer(4)+ibox-1)
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox) 
          npts = npts + iend-istart+1

          if ((npts.gt.0) .and. (nchild.gt.0)) then

            call getpwlistall(ibox,boxsize(ilev),nboxes, &
                itree(ipointer(6)+ibox-1),itree(ipointer(7)+ &
                mnbors*(ibox-1)),nchild,itree(ipointer(5)),centers, &
                isep,nuall,uall,ndall,dall,nnall,nall,nsall,sall, &
                neall,eall,nwall,wall,nu1234,u1234,nd5678,d5678, &
                nn1256,n1256,ns3478,s3478,ne1357,e1357,nw2468,w2468, &
                nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,e13,ne57, &
                e57,nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,ne7,e7, &
                nw2,w2,nw4,w4,nw6,w6,nw8,w8)


            call processudexp(nd,ibox,ilev,nboxes,centers, &
                itree(ipointer(5)),rscales(ilev),boxsize(ilev), &
                nterms(ilev), &
                iaddr,rmlexp,rlams,whts, &
                nlams,nfourier,nphysical,nthmax,nexptot,nexptotp, &
                mexp, &
                nuall,uall,nu1234,u1234,ndall,dall,nd5678,d5678, &
                mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1), &
                mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4), &
                xshift,yshift,zshift,fexpback,rlsc,rscpow,pgboxwexp,cntlist4, &
                list4ct,nlist4tmp,list4,mnlist4)


            call processnsexp(nd,ibox,ilev,nboxes,centers,&
                itree(ipointer(5)),rscales(ilev),boxsize(ilev), &
                nterms(ilev),&
                iaddr,rmlexp,rlams,whts,&
                nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,&
                nnall,nall,nn1256,n1256,nn12,n12,nn56,n56,nsall,sall,&
                ns3478,s3478,ns34,s34,ns78,s78,&
                mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),&
                mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),&
                mexppall(1,1,5),mexppall(1,1,6),mexppall(1,1,7),&
                mexppall(1,1,8),rdplus,xshift,yshift,zshift, &
                fexpback,rlsc,rscpow, pgboxwexp,cntlist4, list4ct, &
                nlist4tmp,list4,mnlist4)

            call processewexp(nd,ibox,ilev,nboxes,centers,&
                itree(ipointer(5)),rscales(ilev),boxsize(ilev),&
                nterms(ilev),&
                iaddr,rmlexp,rlams,whts,&
                nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,&
                neall,eall,ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,&
                ne3,e3,ne5,e5,ne7,e7,nwall,wall,&
                nw2468,w2468,nw24,w24,nw68,w68,&
                nw2,w2,nw4,w4,nw6,w6,nw8,w8,&
                mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),&
                mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),&
                mexppall(1,1,5),mexppall(1,1,6),&
                mexppall(1,1,7),mexppall(1,1,8),mexppall(1,1,9),&
                mexppall(1,1,10),mexppall(1,1,11),mexppall(1,1,12),&
                mexppall(1,1,13),mexppall(1,1,14),mexppall(1,1,15),&
                mexppall(1,1,16),rdminus,xshift,yshift,zshift,&
                fexpback,rlsc,rscpow, pgboxwexp,cntlist4,list4ct, &
                nlist4tmp,list4,mnlist4)
          endif
        enddo

         !deallocate(iboxlexp) 


        !$omp end parallel do        

    end do
      deallocate(iboxsrcind,iboxsrc,iboxpot,iboxgrad,iboxhess)
      deallocate(iboxsubcenters,iboxfl)
      deallocate(uall,dall,nall,sall,eall,wall)
      deallocate(u1234,d5678,n1256,s3478)
      deallocate(e1357,w2468,n12,n56,s34,s78)
      deallocate(e13,e57,w24,w68)
      deallocate(e1,e3,e5,e7,w2,w4,w6,w8)
      deallocate(tmp,mptmp)

  end if

  



  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  timeinfo(4) = time2-time1




  !
  ! ----- Step 5: Downward pass, local-to-local -----
  !
  if(ifprint.ge.1) call prinf('=== Step 5 (split loc) ===*',i,0)

  call cpu_time(time1)
  !$ time1=omp_get_wtime()
  do ilev = 2,nlevels-1



    !$omp parallel do default(shared) &
    !$omp     private(ibox,i,jbox,istart,iend,npts)
    do ibox = laddr(1,ilev),laddr(2,ilev)

      npts = 0

      istart = isrcse(1,ibox) 
      iend = isrcse(2,ibox)
      npts = npts + iend-istart+1

      if (npts .gt. 0) then
        do i=1,8
          jbox = itree(ipointer(5)+8*(ibox-1)+i-1)
          if(jbox.gt.0) then
            call l3dlocloc(nd,rscales(ilev), &
                centers(1,ibox),rmlexp(iaddr(2,ibox)), &
                nterms(ilev),rscales(ilev+1),centers(1,jbox), &
                rmlexp(iaddr(2,jbox)),nterms(ilev+1), &
                dc,lca)
          endif
        enddo
      endif
    enddo
    !$omp end parallel do         
  enddo
  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  timeinfo(5) = time2-time1


  
  !
  ! Step 6: Ship multipole expansions to local expansion in List 4
  !
  if(ifprint.ge.1) call prinf('=== step 6 (mp eval) ===*',i,0)
  call cpu_time(time1)
  !$ time1=omp_get_wtime()

  do ilev=1,nlevels



    !$omp parallel do default(shared) &
    !$omp   private(ibox,istart,iend,npts,i,jbox,j) &
    !$omp   schedule(dynamic)
    do ibox=laddr(1,ilev),laddr(2,ilev)
      istart = isrcse(1,ibox) 
      iend = isrcse(2,ibox) 

      npts = iend-istart+1

      do i=1,nlist3(ibox)
        jbox = list3(i,ibox) 
        do j = istart,iend
          call l3dmploc(nd, rscales(ilev+1), centers(1,jbox), &
              rmlexp(iaddr(1,jbox)),nterms(ilev+1), &
              rmpolesort(j), cmpolesort(1,j), &
              localsort(impolesort(j)), mtermssort(j), &
              dc,lca)
        end do
      enddo
    enddo
    !$omp end parallel do          

  enddo

  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  timeinfo(6) = time2-time1




  !
  ! Step 7: Shift local expansions in leaf nodes to the correct
  ! respective multipole expansion centers
  !
  if(ifprint.ge.1) &
      call prinf('=== Step 7 (LOC to CEN) ===*',i,0)
  print *, 'dont forget to check box radius for translation...'

  call cpu_time(time1)
  !$ time1=omp_get_wtime()

  do ilev = 0,nlevels


    !$omp parallel do default(shared) &
    !$omp   private(ibox,nchild,istart,iend,npts,i) schedule(dynamic)
    do ibox = laddr(1,ilev),laddr(2,ilev)

      nchild=itree(ipointer(4)+ibox-1)

      if(nchild.eq.0) then 
        istart = isrcse(1,ibox) 
        iend = isrcse(2,ibox) 
        npts = iend-istart+1

        do i = istart,iend
          call l3dlocloc(nd, rscales(ilev), &
              centers(1,ibox), rmlexp(iaddr(2,ibox)), &
              nterms(ilev), rmpolesort(i), cmpolesort(1,i), &
              localsort(impolesort(i)), mtermssort(i), &
              dc,lca)
        end do
        
      endif

    enddo
    !$omp end parallel do

  enddo

  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  timeinfo(7) = time2 - time1



  
  !
  ! Step 8: Direct multipole to local translations for nearfield
  !
  if(ifprint .ge. 1) call prinf('=== STEP 8 (direct) =====*',i,0)

  call cpu_time(time1)
  !$ time1=omp_get_wtime()

  do ilev=0,nlevels

    nquad2 = nterms(ilev)*2
    nquad2 = max(6,nquad2)
    ifinit2 = 1
    call legewhts(nquad2, xnodes, wts, ifinit2)
    radius = boxsize(ilev)/2*sqrt(3.0d0)/2/2

    !
    ! do final mploc translations for nearfield, which are assumed to
    ! be valid translations, despite them being in the nearfield
    !

    !$omp parallel do default(shared) &
    !$omp   private(ibox,istarts,iends,npts0,i) &
    !$omp   private(jbox,jstart,jend,npts,d,j,iloc)
    do ibox = laddr(1,ilev),laddr(2,ilev)
      istarts = isrcse(1,ibox) 
      iends = isrcse(2,ibox)
      npts0 = iends-istarts+1

      do iloc = istarts,iends

        do i=1,nlist1(ibox)
          jbox = list1(i,ibox) 
          jstart = isrcse(1,jbox) 
          jend = isrcse(2,jbox)
          npts = jend-jstart+1

          do j = jstart,jend

            d = (cmpolesort(1,j)-cmpolesort(1,iloc))**2 &
                + (cmpolesort(2,j)-cmpolesort(2,iloc))**2 &
                + (cmpolesort(3,j)-cmpolesort(3,iloc))**2
            d = sqrt(d)

            if (d .gt. thresh) then
              call l3dmploc(nd, rmpolesort(j),&
                  cmpolesort(1,j), &
                  mpolesort(impolesort(j)), mtermssort(j), &
                  rmpolesort(iloc), cmpolesort(1,iloc), &
                  localsort(impolesort(iloc)), mtermssort(iloc), &
                  dc,lca)
            else
              if (j .ne. iloc) then
                print *, 'two MP centers closer than thresh... '
                print *, 'thresh = ', thresh
                print *, 'bombing code!!'
                stop
              end if
            end if

          end do


        enddo
      end do


    enddo
    !$omp end parallel do            

  enddo

  
  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  timeinfo(8) = time2-time1

  
  d = 0
  do i = 1,8
    d = d + timeinfo(i)
  enddo

  if (ifprint .eq. 1) then
    print *
    print *
    write(6,'(a)') '                             Time    Perct'
    write(6,'(a,f6.3,a,f6.2,a)') 'Step 1: SHIFT MPs           ',&
        timeinfo(1), ' ', timeinfo(1)/d*100, '%'
    write(6,'(a,f6.3,a,f6.2,a)') 'Step 2: FORM LOCAL (LIST 3) ',&
        timeinfo(2), ' ', timeinfo(2)/d*100, '%'
    write(6,'(a,f6.3,a,f6.2,a)') 'Step 3: MERGE MPs           ',&
        timeinfo(3), ' ', timeinfo(3)/d*100, '%'
    write(6,'(a,f6.3,a,f6.2,a)') 'Step 4: MP to LOCAL         ',&
        timeinfo(4), ' ', timeinfo(4)/d*100, '%'
    write(6,'(a,f6.3,a,f6.2,a)') 'Step 5: SPLIT LOCAL         ',&
        timeinfo(5), ' ', timeinfo(5)/d*100, '%'
    write(6,'(a,f6.3,a,f6.2,a)') 'Step 6: MP EVAL (LIST 4)    ',&
        timeinfo(6), ' ', timeinfo(6)/d*100, '%'
    write(6,'(a,f6.3,a,f6.2,a)') 'Step 7: EVAL LOCALS         ',&
        timeinfo(7), ' ', timeinfo(7)/d*100, '%'
    write(6,'(a,f6.3,a,f6.2,a)') 'Step 8: DIRECT EVAL         ',&
        timeinfo(8), ' ', timeinfo(8)/d*100, '%'
    write(6,'(a,f6.3)') 'Total time required         ', d
    print *
    print *
  end if


  return
end subroutine hfmm3dmain_mps
