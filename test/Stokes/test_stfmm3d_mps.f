      program test_lfmm3d_mps
      implicit none
  
  
      integer :: ns, nt, nc
      integer :: i,j,k,ntest,nd,idim,ier,ijk,ilen,isuccess
      integer :: ifcharge,ifdipole,ifpgh,ifpghtarg
      integer :: ipass(18),len1,ntests,isum
      integer, allocatable :: nterms(:), impole(:)
      integer :: interms,lca
      integer, allocatable :: tmp_vec(:)
      integer :: ifppreg,ifppregtarg,ifstoklet,ifstrslet


      double precision :: dnorm, done, h, pi, rscale,sc, shift, thresh
      double precision :: dnormGrad,errGrad
      integer :: iper, lused, lw, n1, nlege, npts, ns1, ntarg, ntm,ntot

  
      double precision :: eps, err, hkrand, dnorms
      double precision, allocatable :: source(:,:), targ(:,:)
      double precision, allocatable :: centers(:,:)
      double precision, allocatable :: wlege(:), rscales(:)
      double precision, allocatable :: dc(:,:)


      double precision :: eye, zk, ima
      double precision, allocatable :: charge(:,:)
      double precision, allocatable :: dipvec(:,:,:)
      double precision, allocatable :: pot(:,:), pot2(:,:)
      double precision, allocatable :: pottarg(:,:)

      double precision, allocatable :: grad(:,:,:),grad2(:,:,:)
      double precision, allocatable :: gradtarg(:,:,:)
      double precision, allocatable :: hess(:,:,:),hesstarg(:,:,:)
      double complex, allocatable :: mpole(:), local(:)
      double complex, allocatable :: ilocal(:,:)


      double complex, allocatable :: stoklet(:,:), strslet(:,:)
      double complex, allocatable :: strsvec(:,:)
      double complex, allocatable :: potstokes(:,:), prestokes(:)
      double complex, allocatable :: gradstokes(:,:,:)
      double complex, allocatable :: potstokestarg(:,:)
      double complex, allocatable :: prestokestarg(:)
      double complex, allocatable :: gradstokestarg(:,:,:)



      double precision, allocatable :: charge1(:,:),charge2(:,:)
      double precision, allocatable :: charge3(:,:),charge4(:,:)


      double precision, allocatable :: potl1(:,:),potl2(:,:)
      double precision, allocatable :: potl3(:,:),potl4(:,:)

      double precision, allocatable :: gradl1(:,:,:),gradl2(:,:,:)
      double precision, allocatable :: gradl3(:,:,:),gradl4(:,:,:)


      double complex, allocatable :: mpole1(:), mpole2(:)
      double complex, allocatable :: mpole3(:), mpole4(:)
      double complex, allocatable :: local1(:), local2(:)
      double complex, allocatable :: local3(:), local4(:)


      double precision :: pt(3), gl(3), hl(6), vel(3), velgrad(3,3)
      double precision :: press, pl, pv, dmu(3), dnu(3), sigma(3),pl4








      data eye/(0.0d0,1.0d0)/
      ima = (0,1)
      done = 1
      pi = 4*atan(done)


      allocate(dc(0:50,0:50))
      call getsqrtbinomialcoeffs(50,dc)
      lca = 4*50



      call prini(6,13)

      nd = 1


      ns = 100
      nt = 120
      nc = ns
      done = 1

      call prinf('ns = *', ns, 1)
      call prinf('nc = *', nc, 1)
  
 

      allocate(source(3,ns),targ(3,nt), centers(3,nc))
      allocate(charge(nd,ns),dipvec(nd,3,ns))
      allocate(pot(nd,ns), pot2(nd,ns))
      allocate(grad(nd,3,ns),grad2(nd,3,ns))
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





      h = 1.0d0/5
      shift = h/100
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

      call zinitialize(nd*nc, pot2)
      npts = 1

      call zinitialize(nd*nc, potl1)
      call zinitialize(nd*nc, potl2)
      call zinitialize(nd*nc, potl3)
      call zinitialize(nd*nc, potl4)



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



      do i = 1,nc
      call l3dtaevalg(nd, rscales(i),
     1   centers(1,i), local1(impole(i)),
     2   nterms(i), source(1,i), npts, potl1(1,i),gradl1(1,1,i),
     3   wlege, nlege)
      enddo

      call prin2('from lfmm3d_mps, potential = *', potl1, 10)
      call prin2('from lfmm3d_mps, grad = *', gradl1, 10)









c   Second ndl
      do i = 1,nc
        rscales(i) = rscale
      call l3dformmpc(nd,rscale, source(1,i), charge2(1,i),
     1   ns1, centers(1,i), nterms(i), mpole2(impole(i)),
     2   wlege, nlege)
      enddo

  
      call lfmm3d_mps(nd, eps,
     1   nc, centers, rscales, nterms, mpole2, impole, local2,ier)

      do i = 1,nc
      call l3dtaevalg(nd, rscales(i),
     1   centers(1,i), local2(impole(i)),
     2   nterms(i), source(1,i), npts, potl2(1,i),gradl2(1,1,i),
     3   wlege, nlege)
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

      do i = 1,nc
      call l3dtaevalg(nd, rscales(i),
     1   centers(1,i), local3(impole(i)),
     2   nterms(i), source(1,i), npts, potl3(1,i),gradl3(1,1,i),
     3   wlege, nlege)
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

      do i = 1,nc
      call l3dtaevalg(nd, rscales(i),
     1   centers(1,i), local4(impole(i)),
     2   nterms(i), source(1,i), npts, potl4(1,i),gradl4(1,1,i),
     3   wlege, nlege)
      enddo

      call prin2('from lfmm3d_mps, potential = *', potl4, 10)
      call prin2('from lfmm3d_mps, grad = *', gradl4, 10)




      thresh = 1.0d-15
      ifcharge = 1
      ifdipole = 0
      ifpgh = 2
      ntarg = 0
      ifpghtarg = 0
      ier = 0
      call lfmm3d(nd, eps, ns, source, ifcharge,
     1   charge1, ifdipole, dipvec, iper, ifpgh, pot, grad, hess, ntarg,
     2   targ, ifpghtarg, pottarg, gradtarg, hesstarg, ier)
     
      call prin2('via fmm, potential = *', pot, 10)

      call prin2('via fmm, grad = *', grad, 10)



      err = 0
      dnorm = 0
      do j = 1,nc
        do i = 1,nd
          err = err + abs(pot(i,j)-potl1(i,j))**2
          dnorm = dnorm + abs(pot(i,j))**2
          potl1(i,j) = potl1(i,j) - pot(i,j)
        enddo
      enddo
  
      err = sqrt(err/dnorm)
      call prin2('l2 rel err=*',err,1)


      errGrad = 0
      dnormGrad = 0
      do j = 1,nc
         do i = 1,nd
            do k = 1,3
               errGrad = errGrad + abs(grad(i,k,j)-gradl1(i,k,j))**2
                dnormGrad = dnormGrad + abs(grad(i,k,j))**2
            enddo
         enddo
      enddo
  
      errGrad = sqrt(errGrad/dnormGrad)
      call prin2('Grad l2 rel err=*',errGrad,1)





      open(unit=33,file='print_testres.txt',access='append')
      isuccess = 0
      ntest = 1
      if(err.lt.eps) isuccess = 1


      stop
      end program


