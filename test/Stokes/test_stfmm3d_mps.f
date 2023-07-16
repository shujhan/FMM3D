      program test_lfmm3d_mps
      implicit none
  
  
      integer :: ns, nt, nc,ndl,l
      integer :: i,j,k,ntest,nd,idim,ier,ijk,ilen,isuccess,ii
      integer :: ifcharge,ifdipole,ifpgh,ifpghtarg,istress
      integer :: ipass(18),len1,ntests,isum
      integer, allocatable :: nterms(:), impole(:)
      integer :: interms,lca
      integer, allocatable :: tmp_vec(:)
      integer :: ifppreg,ifppregtarg,ifstoklet,ifstrslet,ifppreg1


      double precision :: dnorm, done, h, pi, rscale,sc, shift, thresh
      double precision :: dnormGrad,errGrad, errHess,dnormHess
      integer :: iper, lused, lw, n1, nlege, npts, ns1, ntarg, ntm,ntot
      integer :: npt

      double precision :: derr, derrg,drel,drelg,relerr,relerrg

  
      double precision :: eps, err, hkrand, dnorms
      double precision, allocatable :: source(:,:), targ(:,:)
      double precision, allocatable :: centers(:,:)
      double precision, allocatable :: wlege(:), rscales(:)
      double precision, allocatable :: dc(:,:)


      double precision :: eye, zk, ima
      double precision, allocatable :: charge(:,:)
      double precision, allocatable :: dipvec(:,:,:)
      double precision, allocatable :: pot(:,:),pre(:)
      double precision, allocatable :: pottarg(:,:)

      double precision, allocatable :: grad(:,:,:)
      double precision, allocatable :: gradtarg(:,:,:)
      double precision, allocatable :: hess(:,:,:),hesstarg(:,:,:)
      double complex, allocatable :: mpole(:), local(:)
      double complex, allocatable :: ilocal(:,:)


      double precision, allocatable :: stoklet(:,:), strslet(:,:)
      double precision, allocatable :: strsvec(:,:)
      double precision, allocatable :: potstokes(:,:), prestokes(:)
      double precision, allocatable :: gradstokes(:,:,:)
      double precision, allocatable :: potstokestarg(:,:)
      double precision, allocatable :: prestokestarg(:)
      double precision, allocatable :: gradstokestarg(:,:,:)



      double precision, allocatable :: charge1(:,:),charge2(:,:)
      double precision, allocatable :: charge3(:,:),charge4(:,:)


      double precision, allocatable :: potl1(:,:),potl2(:,:)
      double precision, allocatable :: potl3(:,:),potl4(:,:)

      double precision, allocatable :: gradl1(:,:,:),gradl2(:,:,:)
      double precision, allocatable :: gradl3(:,:,:),gradl4(:,:,:)

      double precision, allocatable :: hessl1(:,:,:),hessl2(:,:,:)
      double precision, allocatable :: hessl3(:,:,:),hessl4(:,:,:)

      double precision, allocatable :: pot2(:,:), pre2(:), grad2(:,:,:) 




      double complex, allocatable :: mpole1(:), mpole2(:)
      double complex, allocatable :: mpole3(:), mpole4(:)
      double complex, allocatable :: local1(:), local2(:)
      double complex, allocatable :: local3(:), local4(:)


      double precision :: pt(3), gl(3), hl(6), vel(3), velgrad(3,3)
      double precision :: press, pl, pv, dmu(3), dnu(3), sigma(3),pl4




      double complex, allocatable :: scarray(:)





      data eye/(0.0d0,1.0d0)/
      ima = (0,1)
      done = 1
      pi = 4*atan(done)


      allocate(dc(0:50,0:50))
      call getsqrtbinomialcoeffs(50,dc)
      lca = 4*50



      call prini(6,13)

      nd = 1


      ns = 12*12*12
      nt = 0
      nc = ns
      done = 1

      call prinf('ns = *', ns, 1)
      call prinf('nc = *', nc, 1)
  
 

      allocate(source(3,ns),targ(3,nt), centers(3,nc))
      allocate(charge(nd,ns),dipvec(nd,3,ns))


      allocate(pot(3,ns),pre(ns))
      allocate(grad(3,3,ns))
      allocate(hess(nd,6,ns))


      allocate(pottarg(1,nt))
      allocate(gradtarg(1,3,nt))
      allocate(hesstarg(1,6,nt))





      allocate(stoklet(3,ns),strslet(3,ns),strsvec(3,ns))
      allocate(potstokes(3,ns),prestokes(ns))
      allocate(gradstokes(3,3,ns))
      allocate(potstokestarg(3,ns),prestokestarg(ns))
      allocate(gradstokestarg(3,3,ns))

      allocate(charge1(1,ns),charge2(1,ns),
     1     charge3(1,ns),charge4(1,ns))



      allocate(potl1(nd,ns))
      allocate(potl2(nd,ns))
      allocate(potl3(nd,ns))
      allocate(potl4(nd,ns))

      allocate(gradl1(nd,3,ns),gradl2(nd,3,ns))
      allocate(gradl3(nd,3,ns),gradl4(nd,3,ns))

      allocate(hessl1(nd,6,ns),hessl2(nd,6,ns))
      allocate(hessl3(nd,6,ns),hessl4(nd,6,ns))


      eps = 0.5d-9

      do i = 1,ns
         source(1,i) = cos(i*done)
         source(2,i) = cos(10*i*done)
         source(3,i) = cos(100*i*done)
         charge(1,i) = sin(2*i*done)
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



c     set-up appropriate vector charge and dipole arrays

      do i = 1,ns
         charge1(1,i) = 0
         charge2(1,i) = 0
         charge3(1,i) = 0
         charge4(1,i) = 0

         sigma(1) = stoklet(1,i)
         sigma(2) = stoklet(2,i)
         sigma(3) = stoklet(3,i)

         charge1(1,i) = charge1(1,i) + sigma(1)/2
         charge2(1,i) = charge2(1,i) + sigma(2)/2
         charge3(1,i) = charge3(1,i) + sigma(3)/2


               pl4 = sigma(1)*source(1,i) + sigma(2)*source(2,i) +
     1              sigma(3)*source(3,i)
         charge4(1,i) = charge4(1,i) + pl4/2
      enddo





      h = 1.0d0/14
      shift = h/10000
      call prin2('shift = *', shift, 1)
      do i = 1,ns
        centers(1,i) = source(1,i) + shift
        centers(2,i) = source(2,i)
        centers(3,i) = source(3,i)
      enddo




      allocate(nterms(nc), impole(nc))

      call l3dterms(eps, ntm)
      ntot = 0
      do i = 1,nc
        nterms(i) = ntm 
        ntot = ntot + (nterms(i)+1)*(2*nterms(i)+1)
      enddo


      allocate(scarray(10*(ntm+2)**2))
      call l3dtaevalhessdini(ntm,scarray)





      allocate(mpole(nd*ntot))

      allocate(mpole1(nd*ntot))
      allocate(mpole2(nd*ntot))
      allocate(mpole3(nd*ntot))
      allocate(mpole4(nd*ntot))



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

      call zinitialize(nd*ntot*2, mpole1)
      call zinitialize(nd*ntot*2, mpole2)
      call zinitialize(nd*ntot*2, mpole3)
      call zinitialize(nd*ntot*2, mpole4)
      
      ns1 = 1
      rscale = 1
      sc = shift
      if (sc .lt. 1) rscale = sc
      call prin2('rscale = *', rscale, 1)
     


      allocate(local(nd*ntot))
      allocate(local1(nd*ntot))
      allocate(local2(nd*ntot))
      allocate(local3(nd*ntot))
      allocate(local4(nd*ntot))

      npts = 1

      call zinitialize(nd*nc, potl1)
      call zinitialize(nd*nc, potl2)
      call zinitialize(nd*nc, potl3)
      call zinitialize(nd*nc, potl4)

      thresh = 1.0d-15



c   First ndl
      allocate(rscales(nc))
      do i = 1,nc
        rscales(i) = rscale
      call l3dformmpc(nd,rscale, source(1,i), charge1(1,i),
     1   ns1, centers(1,i), nterms(i), mpole1(impole(i)),
     2   wlege, nlege)
      enddo



  
      call lfmm3d_mps(nd, eps,
     1   nc, centers, rscales, nterms, mpole1, impole, local1,ier)



c      do i = 1,nc
c      call l3dtaevalg(nd, rscales(i),
c     1   centers(1,i), local1(impole(i)),
c     2   nterms(i), source(1,i), npts, potl1(1,i),gradl1(1,1,i),
c     3   wlege, nlege)
c      enddo

      do i = 1,nc
      call l3dtaevalh(nd, rscales(i),
     1   centers(1,i), local1(impole(i)),
     2   nterms(i), source(1,i), npts, potl1(1,i),gradl1(1,1,i),
     3   hessl1(1,1,i),scarray)
      enddo




      call prin2('from lfmm3d_mps, potential = *', potl1, 10)
      call prin2('from lfmm3d_mps, grad = *', gradl1, 10)
      call prin2('from lfmm3d_mps, hess = *', hessl1, 10)









c   Second ndl
      do i = 1,nc
        rscales(i) = rscale
      call l3dformmpc(nd,rscale, source(1,i), charge2(1,i),
     1   ns1, centers(1,i), nterms(i), mpole2(impole(i)),
     2   wlege, nlege)
      enddo

  
      call lfmm3d_mps(nd, eps,
     1   nc, centers, rscales, nterms, mpole2, impole, local2,ier)

c     do i = 1,nc
c      call l3dtaevalg(nd, rscales(i),
c     1   centers(1,i), local2(impole(i)),
c     2   nterms(i), source(1,i), npts, potl2(1,i),gradl2(1,1,i),
c     3   wlege, nlege)
c      enddo

      do i = 1,nc
      call l3dtaevalh(nd, rscales(i),
     1   centers(1,i), local2(impole(i)),
     2   nterms(i), source(1,i), npts, potl2(1,i),gradl2(1,1,i),
     3   hessl2(1,1,i),scarray)
      enddo

      call prin2('from lfmm3d_mps, potential = *', potl2, 10)
      call prin2('from lfmm3d_mps, grad = *', gradl2, 10)







c   Third ndl
      do i = 1,nc
        rscales(i) = rscale
      call l3dformmpc(nd,rscale, source(1,i), charge3(1,i),
     1   ns1, centers(1,i), nterms(i), mpole3(impole(i)),
     2   wlege, nlege)
      enddo

  
      call lfmm3d_mps(nd, eps,
     1   nc, centers, rscales, nterms, mpole3, impole, local3,ier)

c      do i = 1,nc
c      call l3dtaevalg(nd, rscales(i),
c     1   centers(1,i), local3(impole(i)),
c     2   nterms(i), source(1,i), npts, potl3(1,i),gradl3(1,1,i),
c     3   wlege, nlege)
c      enddo


      do i = 1,nc
      call l3dtaevalh(nd, rscales(i),
     1   centers(1,i), local3(impole(i)),
     2   nterms(i), source(1,i), npts, potl3(1,i),gradl3(1,1,i),
     3   hessl3(1,1,i),scarray)
      enddo




      call prin2('from lfmm3d_mps, potential = *', potl3, 10)
      call prin2('from lfmm3d_mps, grad = *', gradl3, 10)







c   Fourth ndl
      do i = 1,nc
        rscales(i) = rscale
      call l3dformmpc(nd,rscale, source(1,i), charge4(1,i),
     1   ns1, centers(1,i), nterms(i), mpole4(impole(i)),
     2   wlege, nlege)
      enddo

  
      call lfmm3d_mps(nd, eps,
     1   nc, centers, rscales, nterms, mpole4, impole, local4,ier)

c      do i = 1,nc
c      call l3dtaevalg(nd, rscales(i),
c     1   centers(1,i), local4(impole(i)),
c     2   nterms(i), source(1,i), npts, potl4(1,i),gradl4(1,1,i),
c     3   wlege, nlege)
c      enddo


      do i = 1,nc
      call l3dtaevalh(nd, rscales(i),
     1   centers(1,i), local4(impole(i)),
     2   nterms(i), source(1,i), npts, potl4(1,i),gradl4(1,1,i),
     3   hessl4(1,1,i),scarray)
      enddo



      call prin2('from lfmm3d_mps, potential = *', potl4, 10)
      call prin2('from lfmm3d_mps, grad = *', gradl4, 10)



c      open(4, file = 'potl4.txt')  
c      do i=1,100  
c         write(4,*) potl4(1,i)
c      end do 
c      close(4) 


      open(4, file = 'gradl.txt')  
      do i=1,100  
         write(4,*) gradl1(1,1,i),gradl1(1,2,i),gradl2(1,1,i)
      end do 
      close(4) 



c     we only consider no targ condition now

      ntarg = 0
      npt = ntarg + ns

c     unpack stacked Laplace FMM calls

c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l,ifppreg1,ii)
c$OMP$ PRIVATE(pt,pl,gl,hl,vel,velgrad,press)      
      do i = 1,npt

         do j = 1,nd

            vel(1) = 0
            vel(2) = 0
            vel(3) = 0
            velgrad(1,1) = 0
            velgrad(2,1) = 0
            velgrad(3,1) = 0
            velgrad(1,2) = 0
            velgrad(2,2) = 0
            velgrad(3,2) = 0
            velgrad(1,3) = 0
            velgrad(2,3) = 0
            velgrad(3,3) = 0
            press = 0

c process l = 1
            if (i .gt. ntarg) then
               ifppreg1 = ifppreg
               ii = i-ntarg
               pt(1) = source(1,ii)
               pt(2) = source(2,ii)
               pt(3) = source(3,ii)
               pl = potl1(j,ii)
               gl(1) = gradl1(j,1,ii)
               gl(2) = gradl1(j,2,ii)
               gl(3) = gradl1(j,3,ii)
               hl(1) = hessl1(j,1,ii)
               hl(2) = hessl1(j,2,ii)
               hl(3) = hessl1(j,3,ii)
               hl(4) = hessl1(j,4,ii)
               hl(5) = hessl1(j,5,ii)
               hl(6) = hessl1(j,6,ii)                       
            endif


            vel(1) = vel(1) + pl
            vel(1) = vel(1) - pt(1)*gl(1)
            vel(2) = vel(2) - pt(1)*gl(2)
            vel(3) = vel(3) - pt(1)*gl(3)
            press = press - gl(1)*2
            if (ifppreg1 .eq. 3) then
               velgrad(1,1) =  velgrad(1,1) + gl(1)
               velgrad(2,1) =  velgrad(2,1) + gl(2)
               velgrad(3,1) =  velgrad(3,1) + gl(3)

               velgrad(1,1) = velgrad(1,1) - gl(1)
               velgrad(1,2) = velgrad(1,2) - gl(2)
               velgrad(1,3) = velgrad(1,3) - gl(3)                     

c     confirm hessian ordering convention...
               velgrad(1,1) = velgrad(1,1) - pt(1)*hl(1)
               velgrad(2,1) = velgrad(2,1) - pt(1)*hl(4)
               velgrad(3,1) = velgrad(3,1) - pt(1)*hl(5)
               velgrad(1,2) = velgrad(1,2) - pt(1)*hl(4)
               velgrad(2,2) = velgrad(2,2) - pt(1)*hl(2)
               velgrad(3,2) = velgrad(3,2) - pt(1)*hl(6)
               velgrad(1,3) = velgrad(1,3) - pt(1)*hl(5)
               velgrad(2,3) = velgrad(2,3) - pt(1)*hl(6)
               velgrad(3,3) = velgrad(3,3) - pt(1)*hl(3)
            endif



c  when l = 2
            if (i .gt. ntarg) then
               ifppreg1 = ifppreg
               ii = i-ntarg
               pt(1) = source(1,ii)
               pt(2) = source(2,ii)
               pt(3) = source(3,ii)
               pl = potl2(j,ii)
               gl(1) = gradl2(j,1,ii)
               gl(2) = gradl2(j,2,ii)
               gl(3) = gradl2(j,3,ii)
               hl(1) = hessl2(j,1,ii)
               hl(2) = hessl2(j,2,ii)
               hl(3) = hessl2(j,3,ii)
               hl(4) = hessl2(j,4,ii)
               hl(5) = hessl2(j,5,ii)
               hl(6) = hessl2(j,6,ii)                       
            endif


            vel(2) = vel(2) + pl
            vel(1) = vel(1) - pt(2)*gl(1)
            vel(2) = vel(2) - pt(2)*gl(2)
            vel(3) = vel(3) - pt(2)*gl(3)
            press = press - gl(2)*2
            if (ifppreg1 .eq. 3) then
               velgrad(1,2) =  velgrad(1,2) + gl(1)
               velgrad(2,2) =  velgrad(2,2) + gl(2)
               velgrad(3,2) =  velgrad(3,2) + gl(3)

               velgrad(2,1) = velgrad(2,1) - gl(1)
               velgrad(2,2) = velgrad(2,2) - gl(2)
               velgrad(2,3) = velgrad(2,3) - gl(3)                     

c     confirm hessian ordering convention...
               velgrad(1,1) = velgrad(1,1) - pt(2)*hl(1)
               velgrad(2,1) = velgrad(2,1) - pt(2)*hl(4)
               velgrad(3,1) = velgrad(3,1) - pt(2)*hl(5)
               velgrad(1,2) = velgrad(1,2) - pt(2)*hl(4)
               velgrad(2,2) = velgrad(2,2) - pt(2)*hl(2)
               velgrad(3,2) = velgrad(3,2) - pt(2)*hl(6)
               velgrad(1,3) = velgrad(1,3) - pt(2)*hl(5)
               velgrad(2,3) = velgrad(2,3) - pt(2)*hl(6)
               velgrad(3,3) = velgrad(3,3) - pt(2)*hl(3)
            endif



c l = 3

            if (i .gt. ntarg) then
               ifppreg1 = ifppreg
               ii = i-ntarg
               pt(1) = source(1,ii)
               pt(2) = source(2,ii)
               pt(3) = source(3,ii)
               pl = potl3(j,ii)
               gl(1) = gradl3(j,1,ii)
               gl(2) = gradl3(j,2,ii)
               gl(3) = gradl3(j,3,ii)
               hl(1) = hessl3(j,1,ii)
               hl(2) = hessl3(j,2,ii)
               hl(3) = hessl3(j,3,ii)
               hl(4) = hessl3(j,4,ii)
               hl(5) = hessl3(j,5,ii)
               hl(6) = hessl3(j,6,ii)                       
            endif


            vel(3) = vel(3) + pl
            vel(1) = vel(1) - pt(3)*gl(1)
            vel(2) = vel(2) - pt(3)*gl(2)
            vel(3) = vel(3) - pt(3)*gl(3)
            press = press - gl(3)*2
            if (ifppreg1 .eq. 3) then
               velgrad(1,3) =  velgrad(1,3) + gl(1)
               velgrad(2,3) =  velgrad(2,3) + gl(2)
               velgrad(3,3) =  velgrad(3,3) + gl(3)

               velgrad(3,1) = velgrad(3,1) - gl(1)
               velgrad(3,2) = velgrad(3,2) - gl(2)
               velgrad(3,3) = velgrad(3,3) - gl(3)                     

c     confirm hessian ordering convention...
               velgrad(1,1) = velgrad(1,1) - pt(3)*hl(1)
               velgrad(2,1) = velgrad(2,1) - pt(3)*hl(4)
               velgrad(3,1) = velgrad(3,1) - pt(3)*hl(5)
               velgrad(1,2) = velgrad(1,2) - pt(3)*hl(4)
               velgrad(2,2) = velgrad(2,2) - pt(3)*hl(2)
               velgrad(3,2) = velgrad(3,2) - pt(3)*hl(6)
               velgrad(1,3) = velgrad(1,3) - pt(3)*hl(5)
               velgrad(2,3) = velgrad(2,3) - pt(3)*hl(6)
               velgrad(3,3) = velgrad(3,3) - pt(3)*hl(3)
            endif




c     l = 4

            if (i .gt. ntarg) then
               ifppreg1 = ifppreg
               ii = i-ntarg
               pt(1) = source(1,ii)
               pt(2) = source(2,ii)
               pt(3) = source(3,ii)
               pl = potl4(j,ii)
               gl(1) = gradl4(j,1,ii)
               gl(2) = gradl4(j,2,ii)
               gl(3) = gradl4(j,3,ii)
               hl(1) = hessl4(j,1,ii)
               hl(2) = hessl4(j,2,ii)
               hl(3) = hessl4(j,3,ii)
               hl(4) = hessl4(j,4,ii)
               hl(5) = hessl4(j,5,ii)
               hl(6) = hessl4(j,6,ii)                       
            endif



               vel(1) = vel(1) + gl(1)
               vel(2) = vel(2) + gl(2)
               vel(3) = vel(3) + gl(3)                  
               if (ifppreg1 .eq. 3) then
c     confirm hessian ordering convention...
                  velgrad(1,1) = velgrad(1,1) + hl(1)
                  velgrad(2,1) = velgrad(2,1) + hl(4)
                  velgrad(3,1) = velgrad(3,1) + hl(5)
                  velgrad(1,2) = velgrad(1,2) + hl(4)
                  velgrad(2,2) = velgrad(2,2) + hl(2)
                  velgrad(3,2) = velgrad(3,2) + hl(6)
                  velgrad(1,3) = velgrad(1,3) + hl(5)
                  velgrad(2,3) = velgrad(2,3) + hl(6)
                  velgrad(3,3) = velgrad(3,3) + hl(3)
               endif



c    end part

            if (i .gt. ntarg) then
               if (ifppreg1 .ge. 1) then
                  pot(1,ii) = vel(1)
                  pot(2,ii) = vel(2)
                  pot(3,ii) = vel(3)
               endif
               if (ifppreg1 .ge. 2) then
                  pre(ii) = press
               endif
               if (ifppreg1 .ge. 3) then
                  grad(1,1,ii) = velgrad(1,1)
                  grad(2,1,ii) = velgrad(2,1)
                  grad(3,1,ii) = velgrad(3,1)
                  grad(1,2,ii) = velgrad(1,2)
                  grad(2,2,ii) = velgrad(2,2)
                  grad(3,2,ii) = velgrad(3,2)
                  grad(1,3,ii) = velgrad(1,3)
                  grad(2,3,ii) = velgrad(2,3)
                  grad(3,3,ii) = velgrad(3,3)
               endif
            endif               
         enddo
      enddo
c$OMP END PARALLEL DO



      open(2, file = 'potMPS.txt')  
      do i=1,100  
         write(2,*) pot(1,i),pot(2,i),pot(3,i)
      enddo 
      close(2) 












c  test error for pot and grad we get 


      ntest = 100
      
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
      if (relerr .lt. eps) ipass(1) = 1
      if (relerrg .lt. eps) ipass(2) = 1
      

      call prin2('rel err pot srcs *',relerr,1)
      call prin2('rel err grad srcs *',relerrg,1)






      open(unit=33,file='print_testres.txt',access='append')
      isuccess = 0
      if(err.lt.eps) isuccess = 1


      stop
      end program


