program test_lfmm3d_mp2loc
  implicit double precision (a-h,o-z)
 ! character(len=72) str1
 !  implicit none
      integer ns,nt,nd
      double precision, allocatable :: source(:,:),targ(:,:)
      double precision, allocatable :: charge(:)
      double precision, allocatable :: dipvec(:,:)
      double precision, allocatable :: pot(:),pottarg(:)
      double precision, allocatable :: grad(:,:),gradtarg(:,:)
      double precision, allocatable :: hess(:,:),hesstarg(:,:)

      double precision eps

      double complex eye
      integer i,j,k,ntest,ier
      integer ifcharge,ifdipole,ifpgh,ifpghtarg
      double precision err,hkrand
      integer ipass(27),len1,ntests,isum
      character(len=100) str1
      

      data eye/(0.0d0,1.0d0)/

!
!      initialize printing routine
!
      call prini(6,13)

      ns = 5000 
      nt = 4999

      ntest = 20

      allocate(source(3,ns),targ(3,nt))
      allocate(charge(ns),dipvec(3,ns))
      allocate(pot(ns))
      allocate(grad(3,ns))
      allocate(hess(6,ns))

      allocate(pottarg(nt))
      allocate(gradtarg(3,nt))
      allocate(hesstarg(6,nt))

      eps = 0.51d-3

      write(*,*) "=========================================="
      write(*,*) "Testing suite for lfmm3d"
      write(*,'(a,e11.4)') "Requested precision = ",eps

      open(unit=33,file='print_testres.txt',access='append')

      ntests = 27
      do i=1,ntests
        ipass(i) = 0
      enddo



!
!      generate sources uniformly in the unit cube 
!
!
      do i=1,ns
        source(1,i) = hkrand(0)**2
        source(2,i) = hkrand(0)**2
        source(3,i) = hkrand(0)**2

        charge(i) = hkrand(0) 

        dipvec(1,i) = hkrand(0) 
        dipvec(2,i) = hkrand(0)
        dipvec(3,i) = hkrand(0)

        pot(i) = 0
        grad(1,i) = 0
        grad(2,i) = 0
        grad(3,i) = 0
      enddo


!
!      generate targets uniformly in the unit cube
!
      do i=1,nt
        targ(1,i) = hkrand(0)
        targ(2,i) = hkrand(0)
        targ(3,i) = hkrand(0)

        pottarg(i) = 0
        gradtarg(1,i) = 0
        gradtarg(2,i) = 0
        gradtarg(3,i) = 0 
      enddo


!
!    now test source to source, charge, 
!      with potentials
!
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: charges'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 


  !
  ! do the direct calculation
  !
        ifcharge = 1
        ifdipole = 0
        ifpgh = 1
        ntarg = 0
        ifpghtarg = 0

       ier = 0
       nd = 1

     call lfmm3d(nd, eps, ns, source, ifcharge, &
      charge, ifdipole, dipvec, iper, ifpgh, pot, grad, hess, ntarg, &
      targ, ifpghtarg, pottarg, gradtarg, hesstarg, ier)
  
     call prin2('via fmm, potential = *', pot, 10)

     call lfmm3d_s_c_p(eps,ns,source,charge,pot,ier)


     call prin2('via fmm, potential = *', pot, 10)


  stop
end program

       
