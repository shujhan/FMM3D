      subroutine stfmm3d_mps(nd, eps, 
     $                 nsource, source,
     $                 ifstoklet, stoklet, ifstrslet, strslet, strsvec,
     $                 ifppreg, pot, pre, grad, ntarg, targ, 
     $                 ifppregtarg, pottarg, pretarg, gradtarg,ier)
cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source
cf2py  intent(in) ifstoklet,stoklet
cf2py  intent(in) ifstrslet,strslet,strsvec
cf2py  intent(in) ifppreg,ifppregtarg
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,pre,grad
cf2py  intent(out) pottarg,pretarg,gradtarg
cf2py  intent(out) ier
c
      implicit none
      integer :: nd, ifstoklet, ifstrslet, nsource, ntarg
      integer :: i,j,k,ii,l
      integer :: ns,nc,ntm,ntot,ilen,nlege,lw,lused,ns1
      integer :: npts, npt, ifppreg1, ifppreg,ifppregtarg
      integer :: ier, iper

      double precision :: h, shift,eps, rscale,sc


      double precision, allocatable :: centers(:,:)


      double precision, allocatable :: charge1(:,:),charge2(:,:)
      double precision, allocatable :: charge3(:,:),charge4(:,:)

      double precision, allocatable :: potl1(:,:),potl2(:,:)
      double precision, allocatable :: potl3(:,:),potl4(:,:)

      double precision, allocatable :: gradl1(:,:,:),gradl2(:,:,:)
      double precision, allocatable :: gradl3(:,:,:),gradl4(:,:,:)

      double precision, allocatable :: hessl1(:,:,:),hessl2(:,:,:)
      double precision, allocatable :: hessl3(:,:,:),hessl4(:,:,:)


      double complex, allocatable :: mpole1(:), mpole2(:)
      double complex, allocatable :: mpole3(:), mpole4(:)
      double complex, allocatable :: local1(:), local2(:)
      double complex, allocatable :: local3(:), local4(:)




      integer, allocatable :: nterms(:), impole(:)
      double complex, allocatable :: scarray(:)

      double precision, allocatable :: wlege(:), rscales(:)



      double precision :: pt(3), gl(3), hl(6), vel(3), velgrad(3,3)
      double precision :: press, pl, pv, dmu(3), dnu(3), sigma(3),pl4



      double precision pot(3,nsource),pre(nsource)
      double precision grad(3,3,nsource)
      double precision pottarg(3,ntarg),pretarg(ntarg)
      double precision gradtarg(3,3,ntarg)


      double precision source(3,nsource), targ(3,ntarg)
      

      double precision stoklet(3,nsource)
      double precision strslet(3,nsource)
      double precision strsvec(3,nsource)





c     allocate each ndl 

      ns = nsource
      nc = ns


      allocate(centers(3,nc))

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



      call zinitialize(nd*ntot*2, mpole1)
      call zinitialize(nd*ntot*2, mpole2)
      call zinitialize(nd*ntot*2, mpole3)
      call zinitialize(nd*ntot*2, mpole4)



      ns1 = 1
      rscale = 1
      sc = shift
      if (sc .lt. 1) rscale = sc
      call prin2('rscale = *', rscale, 1)


      allocate(local1(nd*ntot))
      allocate(local2(nd*ntot))
      allocate(local3(nd*ntot))
      allocate(local4(nd*ntot))


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


      do i = 1,nc
      call l3dtaevalh(nd, rscales(i),
     1   centers(1,i), local4(impole(i)),
     2   nterms(i), source(1,i), npts, potl4(1,i),gradl4(1,1,i),
     3   hessl4(1,1,i),scarray)
      enddo



      call prin2('from lfmm3d_mps, potential = *', potl4, 10)
      call prin2('from lfmm3d_mps, grad = *', gradl4, 10)







c     we only consider no targ condition now

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




      return
      end




      subroutine zinitialize(len, zs)
      implicit double precision (a-h,o-z)
      double precision :: zs(len)
      
      do i = 1,len
         zs(i) = 0
      end do
      return
      end subroutine zinitialize























