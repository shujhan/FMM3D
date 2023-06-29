program test_lfmm3d_mps
  implicit none
  
  
  integer :: ns, nt, nc
  integer :: i,j,k,ntest,nd,idim,ier,ijk,ilen,isuccess
  integer :: ifcharge,ifdipole,ifpgh,ifpghtarg
  integer :: ipass(18),len1,ntests,isum
  integer, allocatable :: nterms(:), impole(:)
  integer :: interms,lca
  integer, allocatable :: tmp_vec(:)

  double precision :: dnorm, done, h, pi, rscale,sc, shift, thresh
  integer :: iper, lused, lw, n1, nlege, npts, ns1, ntarg, ntm,ntot

  
  double precision :: eps, err, hkrand, dnorms
  double precision, allocatable :: source(:,:), targ(:,:)
  double precision, allocatable :: centers(:,:)
  double precision, allocatable :: wlege(:), rscales(:)
  double precision, allocatable :: dc(:,:)


  double precision :: eye, zk, ima
  double precision, allocatable :: charge(:,:)
  double precision, allocatable :: dipvec(:,:,:)
  double precision, allocatable :: pot(:,:), pot2(:,:), pottarg(:,:)
  double precision, allocatable :: grad(:,:,:),gradtarg(:,:,:)
  double precision, allocatable :: hess(:,:,:),hesstarg(:,:,:)
  double complex, allocatable :: mpole(:), local(:)
  double complex, allocatable :: ilocal(:,:)


  data eye/(0.0d0,1.0d0)/
  ima = (0,1)
  done = 1
  pi = 4*atan(done)


  allocate(dc(0:50,0:50))
  call getsqrtbinomialcoeffs(50,dc)
  lca = 4*50



  !
  ! initialize printing routine
  !
  call prini(6,13)

  nd = 1


  n1 = 3
  ns = n1**3
  nc = ns

  call prinf('ns = *', ns, 1)
  call prinf('nc = *', nc, 1)
  
  nt = 19

  allocate(source(3,ns),targ(3,nt), centers(3,nc))
  allocate(charge(nd,ns),dipvec(nd,3,ns))
  allocate(pot(nd,ns), pot2(nd,ns))
  allocate(grad(nd,3,ns))
  allocate(hess(nd,6,ns))

  allocate(pottarg(nd,nt))
  allocate(gradtarg(nd,3,nt))
  allocate(hesstarg(nd,6,nt))
  eps = 0.5d-9

  write(*,*) "=========================================="
  write(*,*) "Testing suite for lfmm3d_mps"
  write(*,'(a,e12.5)') "Requested precision = ",eps



  !
  ! generate sources uniformly in the unit cube 
  !
  h = 1.0d0/(n1+1)
  ijk = 0
  do i = 1,n1
    do j = 1,n1
      do k = 1,n1
        ijk = ijk + 1
        source(1,ijk) = h*i
        source(2,ijk) = h*j
        source(3,ijk) = h*k
      enddo
    enddo
  enddo



  
  dnorm = 0
  do i=1,ns

    do idim=1,nd

      charge(idim,i) = hkrand(0) 
      dnorm = dnorm + abs(charge(idim,i))**2
      
      dipvec(idim,1,i) = hkrand(0) 
      dipvec(idim,2,i) = hkrand(0) 
      dipvec(idim,3,i) = hkrand(0) 

      pot(idim,i) = 0
      grad(idim,1,i) = 0
      grad(idim,2,i) = 0
      grad(idim,3,i) = 0
    enddo
  enddo

  dnorm = sqrt(dnorm)
  do i=1,ns
    do idim = 1,nd
      charge(idim,i) = charge(idim,i)/dnorm
    enddo
  enddo
  

  
  shift = h/1000
  call prin2('shift = *', shift, 1)
  do i = 1,ns
    centers(1,i) = source(1,i) + shift
    centers(2,i) = source(2,i)
    centers(3,i) = source(3,i)
  enddo



  !
  ! now form a multipole expansion at each center
  !
  allocate(nterms(nc), impole(nc))

  call l3dterms(eps, ntm)
  !ntm = 7
  ntot = 0
  do i = 1,nc
    nterms(i) = ntm 
    ntot = ntot + (nterms(i)+1)*(2*nterms(i)+1)
  enddo

  allocate(mpole(nd*ntot))

  impole(1) = 1
  do i = 1,nc-1
    ilen = (nterms(i)+1)*(2*nterms(i)+1)
    impole(i+1) = impole(i) + nd*ilen
  enddo

  
  nlege = 400
  lw = 5*(nlege+1)**2
  allocate( wlege(lw) )

  call prinf('before ylgndrfwini, lw = *', lw, 1)
  call ylgndrfwini(nlege, wlege, lw, lused)
  call prinf('after ylgndrfwini, lused = *', lused, 1)

  call zinitialize(nd*ntot*2, mpole)
  
  ns1 = 1
  rscale = 1
  sc = shift
  if (sc .lt. 1) rscale = sc
  call prin2('rscale = *', rscale, 1)
  
  allocate(rscales(nc))
  do i = 1,nc
    rscales(i) = rscale
    call l3dformmpc(nd,rscale, source(1,i), charge(1,i), &
        ns1, centers(1,i), nterms(i), mpole(impole(i)), &
        wlege, nlege)
  enddo



  
  !
  ! do the direct calculation
  !
  thresh = 1.0d-15
  ifcharge = 1
  ifdipole = 0
  ifpgh = 1
  ntarg = 0
  ifpghtarg = 0
  ier = 0
  call lfmm3d(nd, eps, ns, source, ifcharge, &
      charge, ifdipole, dipvec, iper, ifpgh, pot, grad, hess, ntarg, &
    targ, ifpghtarg, pottarg, gradtarg, hesstarg, ier)


 ! call lfmm3d_s_c_p(eps,ns,source,charge,pot,ier)
  
  call prin2('via fmm, potential = *', pot, 10)





  
  allocate(local(nd*ntot))
  local = 0
  
  !
  ! now test source to source, charge, 
  ! with potentials
  !
  print *
  print *
  print *
  write(6,*) 'testing multipoles to locals'
  write(6,*) 'input: multipole expansions'
  write(6,*) 'output: local expansions'
  write(6,*) 
  write(6,*) 

  
  call lfmm3d_mps(nd, eps,  &
      nc, centers, rscales, nterms, mpole, impole, local,ier)

! call prin2('from lfmm3d_mps, local expansions = *', local, 10)


       interms = (nterms(1)+1)*(2*nterms(1)+1)
       allocate(ilocal(nterms(1)+1,-nterms(1):nterms(1)))
       ilocal = 0
       allocate(tmp_vec(interms))

  call zinitialize(nd*nc, pot2)
  npts = 1
  do i = 1,nc
      call l3dtaevalp(nd, rscales(i), &
        centers(1,i), local(impole(i)), &
        nterms(i), source(1,i), npts, pot2(1,i), &
        wlege, nlege)

!   call l3dlocloc(nd, rscales(i),&
 !       centers(1,i), local(impole(i)),&
 !       nterms(i), rscales(i), centers(1,i),&
 !       ilocal, 7, &
 !       dc,lca)

 !   call l3dtaevalp(nd, rscales(i),&
 !        centers(1,i), ilocal,&
 !        7, source(1,i), npts, pot2(1,i),&
 !        wlege, nlege)
 !   tmp_vec = (/(k, k=impole(i),impole(i)+interms-1,1)/)
 !   ilocal = reshape(local(tmp_vec),(/nterms(1)+1,2*nterms(1)+1/))
  enddo

  call prin2('from lfmm3d_mps, potential = *', pot2, 10)




  err = 0
  dnorm = 0
  do j = 1,nc
    do i = 1,nd
      err = err + abs(pot(i,j)-pot2(i,j))**2
      dnorm = dnorm + abs(pot(i,j))**2
      pot2(i,j) = pot2(i,j) - pot(i,j)
    enddo
  enddo
  
  err = sqrt(err/dnorm)
  call prin2('l2 rel err=*',err,1)

  open(unit=33,file='print_testres.txt',access='append')
  isuccess = 0
  ntest = 1
  if(err.lt.eps) isuccess = 1

  write(33,'(a,i1,a,i1,a)') 'Successfully completed ', &
    isuccess,' out of ',ntest,' tests in lap3d_mps testing suite'
  close(33)

  stop
end program




       subroutine zinitialize(len, zs)
         implicit double precision (a-h,o-z)
         double precision :: zs(len)
      
         do i = 1,len
           zs(i) = 0
         end do
         return
       end subroutine zinitialize

