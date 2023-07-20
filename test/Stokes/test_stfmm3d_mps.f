      program test_lfmm3d_mps
      implicit none

      integer :: i,j,k
      integer :: ns, nt, nc,ntarg,nd
      integer :: ifstoklet, ifstrslet, ifppreg, ifppregtarg, ier,iter
      integer :: npts, ntest,istress,iper
      double precision :: eye, ima, done, pi,eps,thresh
      double precision :: derr, derrg,drel,drelg,relerr,relerrg



      double precision, allocatable :: pot(:,:),pre(:)
      double precision, allocatable :: grad(:,:,:)
      double precision, allocatable :: pottarg(:,:),pretarg(:)
      double precision, allocatable :: gradtarg(:,:,:)

      double precision, allocatable :: pot2(:,:), pre2(:), grad2(:,:,:) 


      double precision, allocatable :: source(:,:), targ(:,:)


      double precision, allocatable :: stoklet(:,:), strslet(:,:)
      double precision, allocatable :: strsvec(:,:)




      data eye/(0.0d0,1.0d0)/
      ima = (0,1)
      done = 1
      pi = 4*atan(done)



      ns = 11*11*11
      nt = 0
      nc = ns
      ntarg = nt


      allocate(source(3,ns),targ(3,nt))
      allocate(stoklet(3,ns),strslet(3,ns),strsvec(3,ns))


      allocate(pot(3,ns),pre(ns),grad(3,3,ns))
      allocate(pottarg(3,nt),pretarg(nt),gradtarg(3,3,nt))  



      do i = 1,ns
         source(1,i) = cos(i*done)
         source(2,i) = cos(10*i*done)
         source(3,i) = cos(100*i*done)
         stoklet(1,i) = sin(2*i*done)
         stoklet(2,i) = sin(22*i*done)
         stoklet(3,i) = sin(222*i*done)
      enddo

      do i = 1,nt
         targ(1,i) = cos(7*i*done+1.2)
         targ(2,i) = cos(77*i*done+1.2)
         targ(3,i) = cos(777*i*done+1.2)
      enddo


      nd = 1
      eps = 1d-9
      ifstoklet = 1
      ifstrslet = 0
      ifppreg = 3
      ifppregtarg = 3
      ier = 0
      iper = 0
      npts = 1


      call stfmm3d_mps(nd,eps,ns,source,ifstoklet,stoklet,
     1     ifstrslet,strslet,strsvec,ifppreg,pot,pre,grad,
     2     nt,targ,ifppregtarg,pottarg,pretarg,gradtarg,ier)



c  test error for pot and grad we get 

      ntest = 1000
      
      allocate(pot2(3,ntest),pre2(ntest),grad2(3,3,ntest))

      do i = 1,ntest
         do j = 1,3
            pot2(j,i) = 0
            do k = 1,3
               grad2(k,j,i) = 0
            enddo
         enddo
         pre2(i) = 0
      enddo

      istress = 1
      thresh = 1d-15
      call st3ddirectstokg(nd,source,stoklet,
     1     ns,source,ntest,pot2,pre2,grad2,thresh)


      derr = 0
      drel = 0
      derrg = 0
      drelg = 0
      do i = 1,ntest
         derr = derr + (pot(1,i)-pot2(1,i))**2
         drel = drel + pot2(1,i)**2
         derr = derr + (pot(2,i)-pot2(2,i))**2
         drel = drel + pot2(2,i)**2
         derr = derr + (pot(3,i)-pot2(3,i))**2
         drel = drel + pot2(3,i)**2

         do j = 1,3
            do k = 1,3
               derrg = derrg + (grad(k,j,i)-grad2(k,j,i))**2
               drelg = drelg + grad2(k,j,i)**2
            enddo
         enddo
         
      enddo

      relerr = sqrt(derr/(ns*drel))
      relerrg = sqrt(derrg/(ns*drelg))
      

      call prin2('rel err pot srcs *',relerr,1)
      call prin2('rel err grad srcs *',relerrg,1)


      stop
      end program

















