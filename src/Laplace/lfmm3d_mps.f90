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
  integer :: i, j, l, ijk, iert,ifprint,ilev,idim
  integer:: ier
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
  ier = 0

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


  call pts_tree_sort(nexpc,expc,itree,ltree,nboxes,nlevels, &
     ipointer,treecenters,iexpc,iexpcse)


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
      impolesort, localsort, nexpc, ntj, &
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
    impolesort,nexpc,ntj, localsort, &
    iaddr, rmlexp, lmptot, mptemp, mptemp2, lmptemp, &
    itree, ltree, ipointer, ndiv, nlevels, &
    nboxes, iper, boxsize, centers, isrcse, &
    rscales, laddr, nterms, ier)
      implicit none



      double precision sourcesort(3,nmpole)

      integer ntj,nexpc
      integer ifnear
      double precision expcsort(3,nexpc)
      double complex tsort(nd,0:ntj,-ntj:ntj,nexpc)
      double precision scjsort(nexpc)

    !
    ! INPUT variables
    !
    integer :: nd, ndiv,nlevels
    double precision :: eps
    double complex :: zk

    ! input multipole stuff
    integer :: nmpole, mtermssort(nmpole)
    double precision :: cmpolesort(3,nmpole), rmpolesort(nmpole)
    double complex :: mpolesort(*)
    integer :: impolesort(nmpole)


      ! storage stuff for tree and multipole expansions
      integer *8 iaddr(2,nboxes), lmptot
      integer lmptemp
      double precision rmlexp(lmptot)
      double precision mptemp(lmptemp)
      double precision mptemp2(lmptemp)

  !
  ! OUTPUT variables
  !
    double complex :: localsort(*)


      integer ifcharge,ifdipole,ifpgh
      integer ifinit2,nquad2

      double precision thresh
       
      double precision timeinfo(6)
      double precision centers(3,nboxes)

      integer isep,iper,ier
      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer *8 ipointer(8),ltree
      integer itree(ltree)
      integer nboxes
      double precision rscales(0:nlevels)
      double precision boxsize(0:nlevels)
      integer isrcse(2,nboxes),itargse(2,nboxes),iexpcse(2,nboxes)
      integer, allocatable :: nlist1(:),list1(:,:)
      integer, allocatable :: nlist2(:),list2(:,:)
      integer, allocatable :: nlist3(:),list3(:,:)
      integer, allocatable :: nlist4(:),list4(:,:)

      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      integer, allocatable :: uall(:,:),dall(:,:),nall(:,:)
      integer, allocatable :: sall(:,:),eall(:,:),wall(:,:)
      integer, allocatable :: u1234(:,:),d5678(:,:)
      integer, allocatable :: n1256(:,:),s3478(:,:)
      integer, allocatable :: e1357(:,:),w2468(:,:)
      integer, allocatable :: n12(:,:),n56(:,:),s34(:,:),s78(:,:)
      integer, allocatable :: e13(:,:),e57(:,:),w24(:,:),w68(:,:)
      integer, allocatable :: e1(:,:),e3(:,:),e5(:,:),e7(:,:)
      integer, allocatable :: w2(:,:),w4(:,:),w6(:,:),w8(:,:)

!     temp variables
      integer i,j,k,l,ii,jj,kk,ll,m,idim,igbox
      integer ibox,jbox,ilev,npts,npts0,kbox,dir
      integer nchild

      integer istart,iend,istarts,iends,iloc
      integer istartt,iendt,istarte,iende
      integer isstart,isend,jsstart,jsend
      integer jstart,jend

      integer ifprint

      double precision d,time1,time2,second,omp_get_wtime
      double precision pottmp,fldtmp(3),hesstmp(3)

!     PW variables
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
      double complex, allocatable :: mexpf1(:,:,:),mexpf2(:,:,:)
      double complex, allocatable :: mexpp1(:,:,:),mexpp2(:,:,:),mexppall(:,:,:,:)

      double complex, allocatable :: tmp(:,:,:,:)
      double precision, allocatable :: mptmp(:,:)

      double precision sourcetmp(3)
      double complex chargetmp

      integer ix,iy,iz,ictr
      double precision rtmp
      double complex zmul

      integer nlege, lw7, lused7, itype
      double precision wlege(40000)
      integer nterms_eval(4,0:nlevels)

      double precision radius
      double precision, allocatable :: xnodes(:),wts(:)

      integer max_nodes


      integer mnlist1, mnlist2,mnlist3,mnlist4,mnbors
      double complex eye, ztmp
      double precision alphaj
      integer ctr,nn,iptr1,iptr2
      double precision, allocatable :: rscpow(:)
      double precision pi,errtmp
      double complex ima

      double precision ctmp(3)

!     list 3 variables
      double complex, allocatable :: iboxlexp(:,:,:)
      double precision, allocatable :: iboxsubcenters(:,:,:)
      double precision, allocatable :: iboxpot(:,:,:)
      double precision, allocatable :: iboxgrad(:,:,:,:)
      double precision, allocatable :: iboxhess(:,:,:,:)
      double precision, allocatable :: iboxsrc(:,:,:)
      integer, allocatable :: iboxsrcind(:,:)
      integer, allocatable :: iboxfl(:,:,:)
!     end of list 3 variables
!     list 4 variables
      integer cntlist4
      integer, allocatable :: list4ct(:),ilist4(:)
      double complex, allocatable :: gboxmexp(:,:,:)
      double complex, allocatable :: gboxwexp(:,:,:,:,:)
      double complex, allocatable :: pgboxwexp(:,:,:,:)
      double precision, allocatable :: gboxsubcenters(:,:,:)
      double precision, allocatable :: gboxsort(:,:,:)
      integer, allocatable :: gboxind(:,:)
      integer, allocatable :: gboxfl(:,:,:)
      double precision, allocatable :: gboxcgsort(:,:,:)
      double precision, allocatable :: gboxdpsort(:,:,:,:)
!     end of list 4 variables

!
!   hessian variables
!
      double precision, allocatable :: scarray(:,:)

      integer *8 bigint
      integer iert
      data ima/(0.0d0,1.0d0)/

      integer nthd,ithd
      integer omp_get_max_threads,omp_get_thread_num
      nthd = 1
!C$    nthd=omp_get_max_threads()

      pi = 4.0d0*atan(1.0d0)

      thresh = 2.0d0**(-51)*boxsize(0)


!     ifprint is an internal information printing flag. 
!     Suppressed if ifprint=0.
!     Prints timing breakdown and other things if ifprint=1.
!     Prints timing breakdown, list information, 
!     and other things if ifprint=2.
!       
      ifprint=0


      ifcharge = 1
      ifdipole = 0
      ifpgh = 1
!
!   initialize various tree lists
!
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
      mnbors = 27

      isep = 1
      
      call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize, &
       centers,itree(ipointer(3)),itree(ipointer(4)),&
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
       mnlist2,list2,nlist3,mnlist3,list3, &
       nlist4,mnlist4,list4)
      

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

      nmax = 0
      do i=0,nlevels
         if(nmax.lt.nterms(i)) nmax = nterms(i)
      enddo
      allocate(rscpow(0:nmax))
      allocate(carray(4*nmax+1,4*nmax+1))
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


      nn = 10*(nmax+2)**2
      allocate(scarray(nn,0:nlevels))

      do ilev=0,nlevels
        call l3dtaevalhessdini(nterms(ilev),scarray(1,ilev))
      enddo

!     Compute total number of plane waves
      nexptotp = 0
      nexptot = 0
      nthmax = 0
      nphmax = 0
      nn = 0
      do i=1,nlams
         nexptot = nexptot + nfourier(i)
         nexptotp = nexptotp + nphysical(i)
         if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
         if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
         nn = nn + nphysical(i)*nfourier(i)
      enddo

      allocate(fexpe(nn),fexpo(nn),fexpback(nn))
      allocate(tmp(nd,0:nmax,-nmax:nmax,nthd))
      allocate(mptmp(lmptemp,nthd))

      allocate(xshift(-5:5,nexptotp))
      allocate(yshift(-5:5,nexptotp))
      allocate(zshift(5,nexptotp))

      allocate(mexpf1(nd,nexptot,nthd),mexpf2(nd,nexptot,nthd), mexpp1(nd,nexptotp,nthd))
      allocate(mexpp2(nd,nexptotp,nthd),mexppall(nd,nexptotp,16,nthd))

!
!      NOTE: there can be some memory savings here
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
      

!    compute array of factorials

     
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





  max_nodes = 10000
  allocate(xnodes(max_nodes))
  allocate(wts(max_nodes))



      
      if(ifprint.ge.1) call prin2('end of generating plane wave info*',i,0)
!
!
!     ... set the expansion coefficients to zero
!
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,idim)
      do i=1,nexpc
        do k=-ntj,ntj
          do j = 0,ntj
            do idim=1,nd
              tsort(idim,j,k,i)=0
            enddo
          enddo
        enddo
      enddo
!C$OMP END PARALLEL DO

       
      do i=1,6
        timeinfo(i)=0
      enddo


!       ... set all multipole and local expansions to zero


      do ilev = 0,nlevels
!C$OMP PARALLEL DO DEFAULT(SHARED)
!C$OMP$PRIVATE(ibox)
        do ibox=laddr(1,ilev),laddr(2,ilev)
          call mpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
          call mpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
        enddo
!C$OMP END PARALLEL DO        
      enddo


!      set scjsort

      do ilev=0,nlevels
!C$OMP PARALLEL DO DEFAULT(SHARED)
!C$OMP$PRIVATE(ibox,nchild,istart,iend,i)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            nchild = itree(ipointer(4)+ibox-1)
            if(nchild.gt.0) then
               istart = iexpcse(1,ibox)
               iend = iexpcse(2,ibox) 
               do i=istart,iend
                  scjsort(i) = rscales(ilev)
               enddo
            endif
         enddo
!C$OMP END PARALLEL DO
      enddo


!    initialize legendre function evaluation routines
      nlege = 100
      lw7 = 40000
      call ylgndrfwini(nlege,wlege,lw7,lused7)


!     count number of boxes are in list4
      lca = 4*nmax
      if(ifprint.ge.1)  call prinf('=== STEP 0 list4===*',i,0)
      call cpu_time(time1)
!C$    time1=omp_get_wtime()
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
      allocate(gboxmexp(nd*(nterms(nlevels)+1)* (2*nterms(nlevels)+1),8,cntlist4))



      allocate(gboxsubcenters(3,8,nthd))
      allocate(gboxfl(2,8,nthd))

      nmaxt = 0
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts)
!C$OMP$REDUCTION(max:nmaxt)
      do ibox=1,nboxes
        if(list4ct(ibox).gt.0) then
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = iend-istart+1
          if(npts.gt.nmaxt) nmaxt = npts
        endif
      enddo
!C$OMP END PARALLEL DO

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

!C$OMP PARALLEL DO DEFAULT(SHARED)
!C$OMP$PRIVATE(ibox,istart,iend,jbox,jstart,jend,npts,npts0,i)
!C$OMP$PRIVATE(ithd)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
!C$          ithd=omp_get_thread_num()
            ithd = ithd + 1
            if(list4ct(ibox).gt.0) then
              istart=isrcse(1,ibox)
              iend=isrcse(2,ibox)
              npts = iend-istart+1

              if(npts.gt.0) then
                call subdividebox(sourcesort(1,istart),npts, &
                    centers(1,ibox),boxsize(ilev+1),&
                    gboxind(1,ithd),gboxfl(1,1,ithd),&
                    gboxsubcenters(1,1,ithd))
                call dreorderf(3,npts,sourcesort(1,istart),&
                    gboxsort(1,1,ithd),gboxind(1,ithd))

                do i=1,8
                  if(gboxfl(1,i,ithd).gt.0) then
                    jstart=gboxfl(1,i,ithd)
                    jend=gboxfl(2,i,ithd)
                    npts0=jend-jstart+1
                    jbox=list4ct(ibox)
                    if(ifcharge.eq.1.and.ifdipole.eq.0) then
                      call l3dformmpc(nd,rscales(ilev+1),&
                        gboxsort(1,jstart,ithd),&
                        gboxcgsort(1,jstart,ithd),&
                        npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1),&
                        gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.0.and.ifdipole.eq.1) then
                      call l3dformmpd(nd,rscales(ilev+1),&
                        gboxsort(1,jstart,ithd),&
                        gboxdpsort(1,1,jstart,ithd),&
                        npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1),&
                        gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    if(ifcharge.eq.1.and.ifdipole.eq.1) then
                      call l3dformmpcd(nd,rscales(ilev+1),&
                        gboxsort(1,jstart,ithd),&
                        gboxcgsort(1,jstart,ithd),&
                        gboxdpsort(1,1,jstart,ithd),&
                        npts0,gboxsubcenters(1,i,ithd),nterms(ilev+1),&
                        gboxmexp(1,i,jbox),wlege,nlege)          
                    endif
                    call l3dmpmp(nd,rscales(ilev+1),&
                        gboxsubcenters(1,i,ithd),gboxmexp(1,i,jbox),&
                        nterms(ilev+1),rscales(ilev),centers(1,ibox),&
                        rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
     
                    call mpscale(nd,nterms(ilev+1),gboxmexp(1,i,jbox),&
                        rscpow,tmp(1,0,-nmax,ithd))

!                process up down for current box

                    call mpoletoexp(nd,tmp(1,0,-nmax,ithd),&
                        nterms(ilev+1),nlams,&
                        nfourier,nexptot,mexpf1(1,1,ithd),&
                        mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd),&
                        nlams,rlams,nfourier,&
                        nphysical,nthmax,gboxwexp(1,1,1,i,ithd),&
                        fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd),&
                        nlams,rlams,nfourier,&
                        nphysical,nthmax,gboxwexp(1,1,2,i,ithd),&
                        fexpe,fexpo)

                    call processgboxudexp(nd,gboxwexp(1,1,1,i,ithd),&
                        gboxwexp(1,1,2,i,ithd),i,nexptotp,&
                        pgboxwexp(1,1,jbox,1),pgboxwexp(1,1,jbox,2),&
                        xshift,yshift,zshift)

!                process north-south for current box

                    call rotztoy(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd),&
                        mptmp(1,ithd),rdminus)
                    call mpoletoexp(nd,mptmp(1,ithd),&
                        nterms(ilev+1),nlams,&
                        nfourier,nexptot,mexpf1(1,1,ithd),&
                        mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd),&
                        nlams,rlams,nfourier,&
                        nphysical,nthmax,gboxwexp(1,1,3,i,ithd),&
                        fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd),&
                        nlams,rlams,nfourier,&
                        nphysical,nthmax,gboxwexp(1,1,4,i,ithd),&
                        fexpe,fexpo)

                    call processgboxnsexp(nd,gboxwexp(1,1,3,i,ithd),&
                        gboxwexp(1,1,4,i,ithd),i,nexptotp,&
                        pgboxwexp(1,1,jbox,3),pgboxwexp(1,1,jbox,4),&
                        xshift,yshift,zshift)

!                process east-west for current box

                    call rotztox(nd,nterms(ilev+1),tmp(1,0,-nmax,ithd),&
                        mptmp(1,ithd),rdplus)
                    call mpoletoexp(nd,mptmp(1,ithd),&
                        nterms(ilev+1),nlams,&
                        nfourier,nexptot,mexpf1(1,1,ithd),&
                        mexpf2(1,1,ithd),rlsc)

                    call ftophys(nd,mexpf1(1,1,ithd),&
                        nlams,rlams,nfourier,&
                        nphysical,nthmax,gboxwexp(1,1,5,i,ithd),&
                        fexpe,fexpo)

                    call ftophys(nd,mexpf2(1,1,ithd),&
                        nlams,rlams,nfourier,&
                        nphysical,nthmax,gboxwexp(1,1,6,i,ithd),&
                        fexpe,fexpo)
                
                    call processgboxewexp(nd,gboxwexp(1,1,5,i,ithd),&
                        gboxwexp(1,1,6,i,ithd),i,nexptotp,&
                        pgboxwexp(1,1,jbox,5),pgboxwexp(1,1,jbox,6),&
                        xshift,yshift,zshift)
                  endif
                enddo
              endif
            endif
         enddo
!C$OMP END PARALLEL DO
      enddo
      deallocate(gboxfl,gboxsubcenters,gboxwexp,gboxcgsort)
      deallocate(gboxdpsort,gboxind,gboxsort)

      call cpu_time(time2)
!C$    time2=omp_get_wtime()
      if(ifprint.ge.1) print *,"mexp list4 time:",time2-time1
      timeinfo(3)=time2-time1
!     end of count number of boxes are in list4




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

       
      if(ifprint .ge. 1)call prinf('=== STEP 2 (merge mp) ====*',i,0)

      call cpu_time(time1)
!C$    time1=omp_get_wtime()

      do ilev=nlevels-1,0,-1
!C$OMP PARALLEL DO DEFAULT(SHARED)
!C$OMP$PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            do i=1,8
               jbox = itree(ipointer(5)+8*(ibox-1)+i-1)
               if(jbox.gt.0) then
                  istart = isrcse(1,jbox)
                  iend = isrcse(2,jbox)
                  npts = iend-istart+1
                  if(npts.gt.0) then
                     call l3dmpmp(nd,rscales(ilev+1),&
                    centers(1,jbox),rmlexp(iaddr(1,jbox)),&
                    nterms(ilev+1),rscales(ilev),centers(1,ibox),&
                    rmlexp(iaddr(1,ibox)),nterms(ilev),dc,lca)
                  endif
               endif
            enddo
         enddo
!C$OMP END PARALLEL DO         
      enddo

      call cpu_time(time2)
!C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

      if(ifprint.ge.1) call prinf('=== Step 3 (mp to loc+formta+mpeval) ===*',i,0)

!      ... step 3, convert multipole expansions into local
!       expansions

      call cpu_time(time1)
!C$        time1=omp_get_wtime()


!     zero out mexp
 

!C$OMP PARALLEL DO DEFAULT(SHARED)
!C$OMP$PRIVATE(i,j,k,idim)
      do k=1,6
        do i=1,nboxes
          do j=1,nexptotp
            do idim=1,nd
              mexp(idim,j,i,k) = 0.0d0
            enddo
          enddo
        enddo
      enddo
!C$OMP END PARALLEL DO      

!     init uall,dall,...,etc arrays
      allocate(uall(200,nthd),dall(200,nthd),nall(120,nthd))
      allocate(sall(120,nthd),eall(72,nthd),wall(72,nthd))
      allocate(u1234(36,nthd),d5678(36,nthd),n1256(24,nthd))
      allocate(s3478(24,nthd))
      allocate(e1357(16,nthd),w2468(16,nthd),n12(20,nthd))
      allocate(n56(20,nthd),s34(20,nthd),s78(20,nthd))
      allocate(e13(20,nthd),e57(20,nthd),w24(20,nthd),w68(20,nthd))
      allocate(e1(20,nthd),e3(5,nthd),e5(5,nthd),e7(5,nthd))
      allocate(w2(5,nthd),w4(5,nthd),w6(5,nthd),w8(5,nthd))
      allocate(iboxsubcenters(3,8,nthd))
      allocate(iboxfl(2,8,nthd))

!  figure out allocations needed for iboxsrc,iboxsrcind,iboxpot
!  and so on

      nmaxt = 0
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,istart,iend,npts)
!C$OMP$REDUCTION(max:nmaxt)
      do ibox=1,nboxes
        if(nlist3(ibox).gt.0) then
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = iend-istart+1
          if(npts.gt.nmaxt) nmaxt = npts

          istart = itargse(1,ibox)
          iend = itargse(2,ibox)
          npts = iend - istart + 1
          if(npts.gt.nmaxt) nmaxt = npts
        endif
      enddo
!C$OMP END PARALLEL DO

      allocate(iboxsrcind(nmaxt,nthd))
      allocate(iboxsrc(3,nmaxt,nthd))
      allocate(iboxpot(nd,nmaxt,nthd))
      allocate(iboxgrad(nd,3,nmaxt,nthd))
      allocate(iboxhess(nd,6,nmaxt,nthd))

      do ilev=2,nlevels
        allocate(iboxlexp(nd*(nterms(ilev)+1)*(2*nterms(ilev)+1),8,nthd))
         rscpow(0) = 1.0d0/boxsize(ilev)
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo



!C$OMP PARALLEL DO DEFAULT (SHARED)
!C$OMP$PRIVATE(ibox,istart,iend,npts)
!C$OMP$PRIVATE(ithd)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            ithd = 0
!C$          ithd=omp_get_thread_num()
            ithd = ithd + 1
            istart = isrcse(1,ibox) 
            iend = isrcse(2,ibox)

            npts = iend-istart+1

            if(npts.gt.0) then
            !rescale the multipole expansion

                call mpscale(nd,nterms(ilev),rmlexp(iaddr(1,ibox)),&
                      rscpow,tmp(1,0,-nmax,ithd))

!                process up down for current box

                call mpoletoexp(nd,tmp(1,0,-nmax,ithd),nterms(ilev),&
                   nlams,nfourier,&
                   nexptot,mexpf1(1,1,ithd),mexpf2(1,1,ithd),rlsc)

                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier,&
               nphysical,nthmax,mexp(1,1,ibox,1),fexpe,fexpo)

                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier,&
               nphysical,nthmax,mexp(1,1,ibox,2),fexpe,fexpo)



!                process north-south for current box

                call rotztoy(nd,nterms(ilev),tmp(1,0,-nmax,ithd),&
                   mptmp(1,ithd),rdminus)
                call mpoletoexp(nd,mptmp(1,ithd),nterms(ilev),&
                   nlams,nfourier,&
                   nexptot,mexpf1(1,1,ithd),mexpf2(1,1,ithd),rlsc)

                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier,&
               nphysical,nthmax,mexp(1,1,ibox,3),fexpe,fexpo)

                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier,&
               nphysical,nthmax,mexp(1,1,ibox,4),fexpe,fexpo)


!                process east-west for current box

                call rotztox(nd,nterms(ilev),tmp(1,0,-nmax,ithd),&
                   mptmp(1,ithd),rdplus)
                call mpoletoexp(nd,mptmp(1,ithd),&
                   nterms(ilev),nlams,nfourier,&
                   nexptot,mexpf1(1,1,ithd),&
                   mexpf2(1,1,ithd),rlsc)

                call ftophys(nd,mexpf1(1,1,ithd),nlams,rlams,nfourier,&
               nphysical,nthmax,mexp(1,1,ibox,5),fexpe,fexpo)


                call ftophys(nd,mexpf2(1,1,ithd),nlams,rlams,nfourier,&
               nphysical,nthmax,mexp(1,1,ibox,6),fexpe,fexpo)

            endif

         enddo
!C$OMP END PARALLEL DO         


!         loop over parent boxes and ship plane wave
!          expansions to the first child of parent 
!          boxes. 
!          The codes are now written from a gathering perspective
!
!          so the first child of the parent is the one
!          recieving all the local expansions
!          coming from all the lists
!

         rscpow(0) = 1.0d0
         rtmp = rscales(ilev)/boxsize(ilev)
         do i=1,nterms(ilev)
            rscpow(i) = rscpow(i-1)*rtmp
         enddo
!C$OMP PARALLEL DO DEFAULT (SHARED)
!C$OMP$PRIVATE(ibox,istart,iend,npts,nchild)
!C$OMP$PRIVATE(nuall,ndall,nnall,nsall)
!C$OMP$PRIVATE(neall,nwall,nu1234,nd5678)
!C$OMP$PRIVATE(nn1256,ns3478,ne1357,nw2468)
!C$OMP$PRIVATE(nn12,nn56,ns34,ns78,ne13,ne57)
!C$OMP$PRIVATE(nw24,nw68,ne1,ne3,ne5,ne7)
!C$OMP$PRIVATE(nw2,nw4,nw6,nw8)
!C$OMP$PRIVATE(npts0,ctmp,jstart,jend,i)
!C$OMP$PRIVATE(ithd)
         do ibox = laddr(1,ilev-1),laddr(2,ilev-1)
           ithd = 0
!C$         ithd=omp_get_thread_num()
           ithd = ithd + 1
           npts = 0

           istart = iexpcse(1,ibox) 
           iend = iexpcse(2,ibox) 
           npts = npts + iend-istart+1

           nchild = itree(ipointer(4)+ibox-1)

           if(ifpgh.gt.0) then
             istart = isrcse(1,ibox) 
             iend = isrcse(2,ibox) 
             npts = npts + iend-istart+1
           endif


           if(npts.gt.0.and.nchild.gt.0) then

               call getpwlistall(ibox,boxsize(ilev),nboxes,&
              itree(ipointer(6)+ibox-1),itree(ipointer(7)+mnbors*(ibox-1)),nchild,itree(ipointer(5)),centers,&
              isep,nuall,uall(1,ithd),ndall,dall(1,ithd),&
              nnall,nall(1,ithd),nsall,sall(1,ithd),&
              neall,eall(1,ithd),nwall,wall(1,ithd),&
              nu1234,u1234(1,ithd),nd5678,d5678(1,ithd),&
              nn1256,n1256(1,ithd),ns3478,s3478(1,ithd),&
              ne1357,e1357(1,ithd),nw2468,w2468(1,ithd),&
              nn12,n12(1,ithd),nn56,n56(1,ithd),ns34,s34(1,ithd),&
              ns78,s78(1,ithd),ne13,e13(1,ithd),ne57,e57(1,ithd),&
             nw24,w24(1,ithd),nw68,w68(1,ithd),ne1,e1(1,ithd),&
              ne3,e3(1,ithd),ne5,e5(1,ithd),ne7,e7(1,ithd),&
              nw2,w2(1,ithd),nw4,w4(1,ithd),nw6,w6(1,ithd),&
              nw8,w8(1,ithd))

               call processudexp(nd,ibox,ilev,nboxes,centers,&
              itree(ipointer(5)),rscales(ilev),boxsize(ilev),&
              nterms(ilev),&
              iaddr,rmlexp,rlams,whts,&
              nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,&
              nuall,uall(1,ithd),nu1234,u1234(1,ithd),&
              ndall,dall(1,ithd),nd5678,d5678(1,ithd),&
             mexpf1(1,1,ithd),mexpf2(1,1,ithd),&
              mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd),&
              mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),&
              mexppall(1,1,4,ithd),xshift,&
              yshift,zshift,fexpback,rlsc,rscpow,&
              pgboxwexp,cntlist4,list4ct,&
              nlist4,list4,mnlist4)
               
               call processnsexp(nd,ibox,ilev,nboxes,centers,&
              itree(ipointer(5)),rscales(ilev),boxsize(ilev),&
              nterms(ilev),&
              iaddr,rmlexp,rlams,whts,&
              nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,&
              nnall,nall(1,ithd),nn1256,n1256(1,ithd),&
              nn12,n12(1,ithd),nn56,n56(1,ithd),nsall,sall(1,ithd),&
              ns3478,s3478(1,ithd),ns34,s34(1,ithd),ns78,s78(1,ithd),&
              mexpf1(1,1,ithd),mexpf2(1,1,ithd),&
              mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd),&
              mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),&
              mexppall(1,1,4,ithd),&
              mexppall(1,1,5,ithd),mexppall(1,1,6,ithd),&
              mexppall(1,1,7,ithd),&
              mexppall(1,1,8,ithd),rdplus,xshift,yshift,zshift,&
              fexpback,rlsc,rscpow,&
              pgboxwexp,cntlist4,list4ct,&
              nlist4,list4,mnlist4)

               
               call processewexp(nd,ibox,ilev,nboxes,centers,&
              itree(ipointer(5)),rscales(ilev),boxsize(ilev),&
              nterms(ilev),&
              iaddr,rmlexp,rlams,whts,&
             nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,&
              neall,eall(1,ithd),ne1357,e1357(1,ithd),&
              ne13,e13(1,ithd),ne57,e57(1,ithd),ne1,e1(1,ithd),&
              ne3,e3(1,ithd),ne5,e5(1,ithd),&
              ne7,e7(1,ithd),nwall,wall(1,ithd),&
              nw2468,w2468(1,ithd),&
              nw24,w24(1,ithd),nw68,w68(1,ithd),&
              nw2,w2(1,ithd),nw4,w4(1,ithd),nw6,w6(1,ithd),&
              nw8,w8(1,ithd),&
              mexpf1(1,1,ithd),mexpf2(1,1,ithd),&
              mexpp1(1,1,ithd),mexpp2(1,1,ithd),mexppall(1,1,1,ithd),&
              mexppall(1,1,2,ithd),mexppall(1,1,3,ithd),&
              mexppall(1,1,4,ithd),&
              mexppall(1,1,5,ithd),mexppall(1,1,6,ithd),&
              mexppall(1,1,7,ithd),mexppall(1,1,8,ithd),&
              mexppall(1,1,9,ithd),&
              mexppall(1,1,10,ithd),mexppall(1,1,11,ithd),&
              mexppall(1,1,12,ithd),&
              mexppall(1,1,13,ithd),mexppall(1,1,14,ithd),&
              mexppall(1,1,15,ithd),&
              mexppall(1,1,16,ithd),rdminus,xshift,yshift,zshift,&
              fexpback,rlsc,rscpow,&
              pgboxwexp,cntlist4,list4ct,nlist4,list4,mnlist4)


            endif

            if(nlist3(ibox).gt.0.and.npts.gt.0) then
              call getlist3pwlistall(ibox,boxsize(ilev),nboxes,&
                  nlist3(ibox),list3(1,ibox),isep,&
                  centers,nuall,uall(1,ithd),ndall,dall(1,ithd),&
                  nnall,nall(1,ithd),&
                  nsall,sall(1,ithd),neall,eall(1,ithd),&
                  nwall,wall(1,ithd))
              do i=1,8
                call mpzero(nd,iboxlexp(1,i,ithd),nterms(ilev))
              enddo

              call processlist3udexplong(nd,ibox,nboxes,centers,&
                  boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams,&
                  whts,nlams,nfourier,nphysical,nthmax,nexptot,&
                  nexptotp,mexp,nuall,uall(1,ithd),ndall,dall(1,ithd),&
                  mexpf1(1,1,ithd),mexpf2(1,1,ithd),&
                  mexpp1(1,1,ithd),mexpp2(1,1,ithd),&
                  mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),&
                  xshift,yshift,zshift,fexpback,rlsc,rscpow)

              call processlist3nsexplong(nd,ibox,nboxes,centers,&
                  boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams,&
                  whts,nlams,nfourier,nphysical,nthmax,nexptot,&
                  nexptotp,mexp,nnall,nall(1,ithd),nsall,sall(1,ithd),&
                  mexpf1(1,1,ithd),mexpf2(1,1,ithd),&
                  mexpp1(1,1,ithd),mexpp2(1,1,ithd),&
                 mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),rdplus,&
                  xshift,yshift,zshift,fexpback,rlsc,rscpow)

              call processlist3ewexplong(nd,ibox,nboxes,centers,&
                  boxsize(ilev),nterms(ilev),iboxlexp(1,1,ithd),rlams,&
                  whts,nlams,nfourier,nphysical,nthmax,nexptot,&
                  nexptotp,mexp,neall,eall(1,ithd),nwall,wall(1,ithd),&
                  mexpf1(1,1,ithd),mexpf2(1,1,ithd),&
                  mexpp1(1,1,ithd),mexpp2(1,1,ithd),&
                  mexppall(1,1,1,ithd),mexppall(1,1,2,ithd),rdminus,&
                  xshift,yshift,zshift,fexpback,rlsc,rscpow)

              if(ifpgh.eq.1) then
                istart = isrcse(1,ibox) 
                iend = isrcse(2,ibox) 
                npts = iend-istart+1
                if(npts.gt.0) then
                  call subdividebox(sourcesort(1,istart),npts,&
                         centers(1,ibox),boxsize(ilev),&
                        iboxsrcind(1,ithd),iboxfl(1,1,ithd),&
                         iboxsubcenters(1,1,ithd))
                  call dreorderf(3,npts,sourcesort(1,istart),&
                         iboxsrc(1,1,ithd),iboxsrcind(1,ithd))

                  do i=1,8
                    if(iboxfl(1,i,ithd).gt.0) then
                      jstart=iboxfl(1,i,ithd)
                      jend=iboxfl(2,i,ithd)
                      npts0=jend-jstart+1
                      if(npts0.gt.0) then
                        call l3dtaevalp(nd,rscales(ilev),&
                          iboxsubcenters(1,i,ithd),iboxlexp(1,i,ithd),&
                          nterms(ilev),iboxsrc(1,jstart,ithd),npts0,&
                          iboxpot(1,jstart,ithd),wlege,nlege)
                      endif
                    endif
                  enddo
                endif
              endif


            endif
         enddo
!C$OMP END PARALLEL DO        
        deallocate(iboxlexp)  
      enddo


      deallocate(iboxsrcind,iboxsrc,iboxpot,iboxgrad,iboxhess)
      deallocate(iboxsubcenters,iboxfl)
      deallocate(uall,dall,nall,sall,eall,wall)
      deallocate(u1234,d5678,n1256,s3478)
      deallocate(e1357,w2468,n12,n56,s34,s78)
      deallocate(e13,e57,w24,w68)
      deallocate(e1,e3,e5,e7,w2,w4,w6,w8)
      deallocate(tmp,mptmp)
      call cpu_time(time2)
!C$        time2=omp_get_wtime()
      timeinfo(3) = timeinfo(3) + time2-time1


      if(ifprint.ge.1)call prinf('=== Step 4 (split loc) ===*',i,0)


      call cpu_time(time1)
!C$        time1=omp_get_wtime()
      do ilev = 2,nlevels-1

!C$OMP PARALLEL DO DEFAULT(SHARED)
!C$OMP$PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = laddr(1,ilev),laddr(2,ilev)

            npts = 0
            istart = iexpcse(1,ibox) 
            iend = iexpcse(2,ibox) 
            npts = npts + iend-istart+1

            if(ifpgh.gt.0) then
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox) 
               npts = npts + iend-istart+1
            endif

            if(npts.gt.0) then
               do i=1,8
                  jbox = itree(ipointer(5)+8*(ibox-1)+i-1)
                  if(jbox.gt.0) then
                     call l3dlocloc(nd,rscales(ilev),&
                     centers(1,ibox),rmlexp(iaddr(2,ibox)),&
                     nterms(ilev),rscales(ilev+1),centers(1,jbox),&
                     rmlexp(iaddr(2,jbox)),nterms(ilev+1),dc,lca)
                  endif
               enddo
            endif
         enddo
!C$OMP END PARALLEL DO         
      enddo
      call cpu_time(time2)
!C$        time2=omp_get_wtime()
      timeinfo(4) = time2-time1


      if(ifprint.ge.1) call prinf('=== step 5 (eval lo) ===*',i,0)


!     ... step 6, evaluate all local expansions


      call cpu_time(time1)
!C$        time1=omp_get_wtime()



!       shift local expansion to local epxanion at expansion centers
!        (note: this part is not relevant for particle codes.
!        it is relevant only for qbx codes)

      do ilev = 0,nlevels
!C$OMP PARALLEL DO DEFAULT(SHARED)
!C$OMP$PRIVATE(ibox,nchild,istart,iend,i)
!C$OMP$SCHEDULE(DYNAMIC)      
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild=itree(ipointer(4)+ibox-1)
            if(nchild.eq.0) then 
               istart = iexpcse(1,ibox) 
               iend = iexpcse(2,ibox) 
               do i=istart,iend

                  call l3dlocloc(nd,rscales(ilev),&
                  centers(1,ibox),rmlexp(iaddr(2,ibox)),&
                  nterms(ilev),rscales(ilev),expcsort(1,i),&
                  tsort(1,0,-ntj,i),ntj,dc,lca)
               enddo
            endif
         enddo
!C$OMP END PARALLEL DO
      enddo

    
      call cpu_time(time2)
!C$        time2=omp_get_wtime()
      timeinfo(5) = time2 - time1


      if(ifprint .ge. 1) call prinf('=== STEP 6 (direct) =====*',i,0)

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


  !  !$omp parallel do default(shared) &
  !  !$omp   private(ibox,istarts,iends,npts0,i) &
  !  !$omp   private(jbox,jstart,jend,npts,d,j,iloc)
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
  !  !$omp end parallel do       

  enddo

  
  call cpu_time(time2)
  !$ time2=omp_get_wtime()
  timeinfo(8) = time2-time1

  
  d = 0
  do i = 1,8
    d = d + timeinfo(i)
  enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)

      return
      end
