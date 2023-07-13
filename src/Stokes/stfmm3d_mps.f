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
c     Stokes FMM in R^{3}: evaluate all pairwise particle
c     interactions (ignoring self-interactions) and
c     interactions with targs.
c      
c     This routine computes sums of the form
c
c       u(x) = sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
c
c     where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
c     stresslet charge, and nu^{(m)} is the stresslet orientation
c     (note that each of these is a 3 vector per source point y^{(m)}).
c     For x a source point, the self-interaction in the sum is omitted. 
c
c     Optionally, the associated pressure p(x) and gradient grad u(x)
c     are returned
c
c       p(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
c          + sum_m T_{ijk}(x,y^{(m)}) PI_{jk} mu^{(m)}_j nu^{(m)}_k
c
c       grad u(x) = grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k]
c
c-----------------------------------------------------------------------
c     INPUT PARAMETERS:
c     
c   nd:    in: integer
c              number of densities
c   
c   eps:   in: double precision
c              requested precision
c
c   nsource in: integer  
c               number of sources
c
c   source  in: double precision (3,nsource)
c               source(k,j) is the kth component of the jth
c               source locations
c
c   ifstoklet  in: integer  
c               Stokeslet charge computation flag
c               ifstoklet = 1   =>  include Stokeslet contribution
c                                   otherwise do not
c 
c   stoklet in: double precision (nd,3,nsource) 
c               Stokeslet charge strengths (sigma vectors above)
c
c   ifstrslet in: integer
c               stresslet computation flag
c               ifstrslet = 1   =>  include standard stresslet
c                                   (type I)
c
c            NOT YET IMPLEMENTED
c      
c               ifstrslet = 2   =>  include symmetric stresslet
c                                   (type II)
c               ifstrslet = 3   =>  include rotlet
c               ifstrslet = 4   =>  include Stokes doublet
c                      otherwise do not include
c
c   strslet  in: double precision (nd,3,nsource) 
c               stresslet strengths (mu vectors above)
c
c   strsvec  in: double precision (nd,3,nsource)   
c               stresslet orientations (nu vectors above)
c
c     ifppreg    in: integer      
c               flag for evaluating potential, gradient, and pressure
c               at the sources
c               ifppreg = 1, only potential
c               ifppreg = 2, potential and pressure
c         GRADIENT NOT IMPLEMENTED
c               ifppreg = 3, potential, pressure, and gradient 
c      
c   ntarg   in: integer  
c              number of targs 
c
c   targ    in: double precision (3,ntarg)
c             targ(k,j) is the kth component of the jth
c             targ location
c      
c   ifppregtarg in: integer
c                flag for evaluating potential, gradient, and pressure
c                at the targets
c                ifppregtarg = 1, only potential
c                ifppregtarg = 2, potential and pressure
c                ifppregtarg = 3, potential, pressure, and gradient
c
c-----------------------------------------------------------------------
c
c   OUTPUT parameters:
c
c   pot   out: double precision(nd,3,nsource) 
c           velocity at the source locations
c      
c   pre   out: double precision(nd,nsource)
c           pressure at the source locations
c      
c         GRADIENT NOT IMPLEMENTED
c   grad   out: double precision(nd,3,3,nsource) 
c              gradient of velocity at the source locations
c              grad(l,i,j,k) is the ith component of the
c              gradient of the jth component of the velocity
c              for the lth density at the kth source location
c     
c   pottarg   out: double precision(nd,3,ntarg) 
c               velocity at the targets
c      
c   pretarg   out: double precision(nd,ntarg)
c               pressure at the targets
c      
c   gradtarg   out: double precision(nd,3,3,ntarg) 
c               gradient of velocity at the targets
c               gradtarg(l,i,j,k) is the ith component of the
c               gradient of the jth component of the velocity
c               for the lth density at the kth target
c     ier     out: integer
c               error flag
c
c     TODO: implement other stresslet options and gradient
c------------------------------------------------------------------
      implicit none
      integer nd, ifstoklet, ifstrslet, ntarg
      double precision eps
      integer nsource, ifppreg, ifppregtarg
      double precision source(3, nsource), targ(3, ntarg)
      double precision stoklet(nd, 3, nsource), strslet(nd, 3, nsource)
      double precision strsvec(nd, 3, nsource)
      double precision pot(nd, 3, nsource), pre(nd,nsource)
      double precision grad(nd, 3, 3, nsource)
      double precision pottarg(nd, 3, ntarg), pretarg(nd,ntarg),
     1     gradtarg(nd, 3, 3, ntarg) 

c     local
      double precision, allocatable :: charge(:,:,:), dipvec(:,:,:,:)
      double precision, allocatable :: potl(:,:,:), gradl(:,:,:,:),
     1     hessl(:,:,:,:), pottargl(:,:,:), gradtargl(:,:,:,:),
     2     hesstargl(:,:,:,:)
      double precision :: pt(3), gl(3), hl(6), vel(3), velgrad(3,3)
      double precision :: press, pl, pv, dmu(3), dnu(3), sigma(3)

      integer ndl, ifchargel, ifdipolel, ifpghl, ifpghtargl
      integer ndper

      integer i, j, ii, ifppreg1, l, npt, ier,iper









      double precision h, shift, rscale, lused,sc,npts
      integer :: nc, ntm, ntot, lw, nlege, ns1, ilen
      double precision, allocatable :: centers(:,:)
      double precision, allocatable :: wlege(:), rscales(:)
      double complex, allocatable :: mpole1(:), mpole2(:)
      double complex, allocatable :: mpole3(:), mpole4(:)
      double complex, allocatable :: local1(:), local2(:)
      double complex, allocatable :: local3(:), local4(:)
      integer, allocatable :: nterms(:), impole(:)


      double complex, allocatable :: local(:,:), mpole(:,:)


      double complex, allocatable :: potl1(:,:)




      double complex, allocatable :: test(:,:,:)
      integer :: init
      init = 1
      allocate(test(4,1,3))

      do i = 1,4
        do j = 1,1
          do ii = 1,3
            test(i,j,ii) = init
            init = init + 1
          enddo
        enddo
      enddo











      nc = nsource
      allocate(centers(3,nc))



      allocate(potl1(nd,nsource))




      

      ndper = 0
      
      if (ifstrslet .eq. 1 .or. ifstoklet .eq. 1) then
         ndper = 4
      endif

      ifdipolel = 0
      ifchargel = 0

      if (ifstoklet .eq. 1) ifchargel = 1
      if (ifstrslet .eq. 1) ifdipolel = 1
      
      ndper = 4
      ndl = ndper*nd

      ifpghl = 3
      ifpghtargl = 3

c     allocate necessary arrays
      
      allocate(charge(ndper,nd,nsource),dipvec(ndper,nd,3,nsource),
     1     potl(ndper,nd,nsource),pottargl(ndper,nd,ntarg),
     2     gradl(ndper,nd,3,nsource),gradtargl(ndper,nd,3,ntarg),
     3     hessl(ndper,nd,6,nsource),hesstargl(ndper,nd,6,ntarg),
     4     stat=ier)
      if(ier .ne. 0) then
         print *, "In stfmm3d: cannot allocate Laplace call storage"
         print *, "ndper =",ndper
         print *, "nd =",nd
         print *, "nsource =",nsource
         print *, "ntarg =",ntarg        
         stop
      endif

c     set-up appropriate vector charge and dipole arrays

c$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l,sigma,dmu,dnu)
c$OMP$ PRIVATE(pl,pv)      
      do i = 1,nsource

         do j = 1,nd
            do l = 1,ndper
               charge(l,j,i) = 0
               dipvec(l,j,1,i) = 0
               dipvec(l,j,2,i) = 0
               dipvec(l,j,3,i) = 0
            enddo
         enddo

         do j = 1,nd
            if (ifstoklet .eq. 1) then
               sigma(1) = stoklet(j,1,i)
               sigma(2) = stoklet(j,2,i)
               sigma(3) = stoklet(j,3,i)
            endif
            if (ifstrslet .ge. 1) then
               dmu(1) = strslet(j,1,i)
               dmu(2) = strslet(j,2,i)
               dmu(3) = strslet(j,3,i)
               dnu(1) = strsvec(j,1,i)
               dnu(2) = strsvec(j,2,i)
               dnu(3) = strsvec(j,3,i)
            endif

            do l = 1,3
               
               if (ifstoklet .eq. 1) then
                  charge(l,j,i) = charge(l,j,i) + sigma(l)/2
               endif
               if (ifstrslet .eq. 1) then
                  dipvec(l,j,1,i) = dipvec(l,j,1,i) - (dmu(l)*dnu(1) + 
     1                 dmu(1)*dnu(l))/2
                  dipvec(l,j,2,i) = dipvec(l,j,2,i) - (dmu(l)*dnu(2) + 
     1                 dmu(2)*dnu(l))/2
                  dipvec(l,j,3,i) = dipvec(l,j,3,i) - (dmu(l)*dnu(3) + 
     1                 dmu(3)*dnu(l))/2
               endif
            enddo
            
            l = 4
            
            if (ifstoklet .eq. 1) then
               pl = sigma(1)*source(1,i) + sigma(2)*source(2,i) +
     1              sigma(3)*source(3,i)
               charge(l,j,i) = charge(l,j,i) + pl/2
            endif
            if (ifstrslet .eq. 1) then
               pl = dmu(1)*source(1,i) + dmu(2)*source(2,i) +
     1              dmu(3)*source(3,i)
               pv = dnu(1)*source(1,i) + dnu(2)*source(2,i) +
     1              dnu(3)*source(3,i)
               
               dipvec(l,j,1,i) = dipvec(l,j,1,i) -
     1              (dmu(1)*pv + dnu(1)*pl)/2
               dipvec(l,j,2,i) = dipvec(l,j,2,i) -
     1              (dmu(2)*pv + dnu(2)*pl)/2
               dipvec(l,j,3,i) = dipvec(l,j,3,i) -
     1              (dmu(3)*pv + dnu(3)*pl)/2
            endif
            
         enddo
         
      enddo
c$OMP END PARALLEL DO      













      h = 1.0d0/5
      shift = h/1000
      call prin2('shift = *', shift, 1)

      do i = 1,nsource
         centers(1,i) = source(1,i) + shift
         centers(2,i) = source(2,i)
         centers(3,i) = source(3,i)
      enddo


      allocate(nterms(nc), impole(nc))

c      ntm = 10
      call l3dterms(eps, ntm)
      call prinf('ntm = *', ntm, 1)
      ntot = 0
      do i = 1,nc
         nterms(i) = ntm 
         ntot = ntot + (nterms(i)+1)*(2*nterms(i)+1)
      enddo
      call prinf('ntot = *', ntot, 1)



      allocate(mpole(ndper,nd*ntot))
      allocate(mpole1(nd*ntot))
      allocate(mpole2(nd*ntot))
      allocate(mpole3(nd*ntot))
      allocate(mpole4(nd*ntot))

      impole(1) = 1
      do i = 1,nc-1
         ilen = (nterms(i)+1)*(2*nterms(i)+1)
         impole(i+1) = impole(i) + nd*ilen
      enddo
      call prinf('ilen = *', ilen, 1)
      call prin2('impole = *', impole, 10)
  
      nlege = 300
      lw = 5*(nlege+1)**2
      allocate( wlege(lw) )

      call prinf('before ylgndrfwini, lw = *', lw, 1)
      call ylgndrfwini(nlege, wlege, lw, lused)
      call prinf('after ylgndrfwini, lused = *', lused, 1)


      call zinitialize(ndper*nd*ntot*2, mpole)
      call zinitialize(nd*ntot*2, mpole1)
      call zinitialize(nd*ntot*2, mpole2)
      call zinitialize(nd*ntot*2, mpole3)
      call zinitialize(nd*ntot*2, mpole4)
      
      ns1 = 1
      rscale = 1
      sc = shift
      if (sc .lt. 1) rscale = sc
      call prin2('rscale = *', rscale, 1)
      
      allocate(rscales(nc))



      allocate(local(ndper,nd*ntot))
      allocate(local1(nd*ntot))
      allocate(local2(nd*ntot))
      allocate(local3(nd*ntot))
      allocate(local4(nd*ntot))

      npts = 1
      iper = 0
      ier = 0







c   First ndl
      do i = 1,nc
         rscales(i) = rscale
      call l3dformmpc(nd,rscale, source(1,i), charge(1,nd,i),
     1     ns1, centers(1,i), nterms(i), mpole1(impole(i)),
     2     wlege, nlege)
      enddo


      call prin2('source = *', source, 10)
      call prin2('charge = *', charge, 10)
      call prin2('centers = *', centers, 10)


      call lfmm3d_mps(nd, eps, nc, centers, rscales, nterms,
     1     mpole1, impole, local1, ier)

 


      call zinitialize(nd*nc*ndper, potl)
      call zinitialize(nd*nc*ndper*3, gradl)


      call zinitialize(nd*nc, potl1)



c      do i = 1,nc
c      call l3dtaevalg(nd, rscales(i),
c     1     centers(1,i), local1(impole(i)), 
c     2     nterms(i), source(1,i), npts, potl(1,1,i),gradl(1,1,1,i), 
c     3     wlege, nlege)
c      enddo


      do i = 1,nc
      call l3dtaevalp(nd, rscales(i),
     1     centers(1,i), local1(impole(i)),
     2     nterms(i), source(1,i), npts, potl1(1,i),
     3     wlege, nlege)
      enddo



      call prin2('via mps, potl = *', potl, 10)
      call prin2('via mps, gradl = *', gradl, 10)


c   Second ndl
      do i = 1,nc
         rscales(i) = rscale
      call l3dformmpc(ndl,rscale, source(1,i), charge(2,1,i),
     1     ns1, centers(1,i), nterms(i), mpole2(impole(i)),
     2     wlege, nlege)
      enddo

      call lfmm3d_mps(nd, eps, nc, centers, rscales, nterms,
     1     mpole2, impole, local2, ier)


      call prin2('local2 = *', local2, 10)


      do i = 1,nc
      call l3dtaevalg(ndl, rscales(i),
     1     centers(1,i), local2(impole(i)), 
     2     nterms(i), source(1,i), npts, potl(2,1,i),gradl(2,1,1,i), 
     3     wlege, nlege)
      enddo
      call prin2('via mps, potl = *', potl, 10)
      call prin2('via mps, gradl = *', gradl, 10)




c   third ndl
      do i = 1,nc
         rscales(i) = rscale
      call l3dformmpc(ndl,rscale, source(1,i), charge(3,1,i),
     1     ns1, centers(1,i), nterms(i), mpole3(impole(i)),
     2     wlege, nlege)
      enddo

      call lfmm3d_mps(nd, eps, nc, centers, rscales, nterms,
     1     mpole3, impole, local3, ier)


      call prin2('local3 = *', local3, 10)

      do i = 1,nc
      call l3dtaevalg(ndl, rscales(i),
     1     centers(1,i), local3(impole(i)), 
     2     nterms(i), source(1,i), npts, potl(3,1,i),gradl(3,1,1,i), 
     3     wlege, nlege)
      enddo
      call prin2('via mps, potl = *', potl, 10)
      call prin2('via mps, gradl = *', gradl, 10)






c   fourth ndl
      do i = 1,nc
         rscales(i) = rscale
      call l3dformmpc(ndl,rscale, source(1,i), charge(4,1,i),
     1     ns1, centers(1,i), nterms(i), mpole4(impole(i)),
     2     wlege, nlege)
      enddo

      call lfmm3d_mps(nd, eps, nc, centers, rscales, nterms,
     1     mpole4, impole, local4, ier)

      call prin2('local4 = *', local4, 10)

      do i = 1,nc
      call l3dtaevalg(ndl, rscales(i),
     1     centers(1,i), local4(impole(i)), 
     2     nterms(i), source(1,i), npts, potl(4,1,i),gradl(4,1,1,i), 
     3     wlege, nlege)
      enddo
      call prin2('via mps, potl = *', potl, 10)
      call prin2('via mps, gradl = *', gradl, 10)






c     Try to put together 

      do j = 1, ndper
        do i = 1,nc
         rscales(i) = rscale
      call l3dformmpc(ndl,rscale, source(1,i), charge(j,1,i),
     1     ns1, centers(1,i), nterms(i), mpole(j,impole(i)),
     2     wlege, nlege)
        enddo
      enddo

      call lfmm3d_mps(ndl,eps,nc,centers,rscales,nterms,
     1     mpole,impole,local,ier)

      do j = 1, ndper
        do i = 1,nc
      call l3dtaevalg(ndl, rscales(i),
     1     centers(1,i), local(j,impole(i)), 
     2     nterms(i), source(1,i), npts, potl(j,1,i),gradl(j,1,1,i), 
     3     wlege, nlege)
        enddo
      enddo
      call prin2('via mps, potl = *', potl, 10)
      call prin2('via mps, gradl = *', gradl, 10)





      call lfmm3d(ndl,eps,nsource,source,ifchargel,charge,
     1     ifdipolel,dipvec,iper,ifpghl,potl,gradl,hessl,ntarg,
     2     targ,ifpghtargl,pottargl,gradtargl,hesstargl,ier)


      call prin2('via fmm, potl = *', potl, 10)
      call prin2('via fmm, gradl = *', gradl, 10)











      npt = ntarg + nsource

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
            
            do l = 1,ndper
               
               if (i .gt. ntarg) then
                  ifppreg1 = ifppreg
                  ii = i-ntarg
                  pt(1) = source(1,ii)
                  pt(2) = source(2,ii)
                  pt(3) = source(3,ii)
                  pl = potl(l,j,ii)
                  gl(1) = gradl(l,j,1,ii)
                  gl(2) = gradl(l,j,2,ii)
                  gl(3) = gradl(l,j,3,ii)
                  hl(1) = hessl(l,j,1,ii)
                  hl(2) = hessl(l,j,2,ii)
                  hl(3) = hessl(l,j,3,ii)
                  hl(4) = hessl(l,j,4,ii)
                  hl(5) = hessl(l,j,5,ii)
                  hl(6) = hessl(l,j,6,ii)            
               else
                  ifppreg1 = ifppregtarg                  
                  ii = i
                  pt(1) = targ(1,ii)
                  pt(2) = targ(2,ii)
                  pt(3) = targ(3,ii)
                  pl = pottargl(l,j,ii)
                  gl(1) = gradtargl(l,j,1,ii)
                  gl(2) = gradtargl(l,j,2,ii)
                  gl(3) = gradtargl(l,j,3,ii)
                  hl(1) = hesstargl(l,j,1,ii)
                  hl(2) = hesstargl(l,j,2,ii)
                  hl(3) = hesstargl(l,j,3,ii)
                  hl(4) = hesstargl(l,j,4,ii)
                  hl(5) = hesstargl(l,j,5,ii)
                  hl(6) = hesstargl(l,j,6,ii)            
               endif

               
               if (l .ge. 1 .and. l .le. 3) then
                  
                  vel(l) = vel(l) + pl
                  vel(1) = vel(1) - pt(l)*gl(1)
                  vel(2) = vel(2) - pt(l)*gl(2)
                  vel(3) = vel(3) - pt(l)*gl(3)
                  press = press - gl(l)*2
                  if (ifppreg1 .eq. 3) then
                     velgrad(1,l) =  velgrad(1,l) + gl(1)
                     velgrad(2,l) =  velgrad(2,l) + gl(2)
                     velgrad(3,l) =  velgrad(3,l) + gl(3)

                     velgrad(l,1) = velgrad(l,1) - gl(1)
                     velgrad(l,2) = velgrad(l,2) - gl(2)
                     velgrad(l,3) = velgrad(l,3) - gl(3)                     

c     confirm hessian ordering convention...
                     velgrad(1,1) = velgrad(1,1) - pt(l)*hl(1)
                     velgrad(2,1) = velgrad(2,1) - pt(l)*hl(4)
                     velgrad(3,1) = velgrad(3,1) - pt(l)*hl(5)
                     velgrad(1,2) = velgrad(1,2) - pt(l)*hl(4)
                     velgrad(2,2) = velgrad(2,2) - pt(l)*hl(2)
                     velgrad(3,2) = velgrad(3,2) - pt(l)*hl(6)
                     velgrad(1,3) = velgrad(1,3) - pt(l)*hl(5)
                     velgrad(2,3) = velgrad(2,3) - pt(l)*hl(6)
                     velgrad(3,3) = velgrad(3,3) - pt(l)*hl(3)
                  endif

               else if (l .eq. 4) then

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
                  
               endif
            enddo

            if (i .gt. ntarg) then
               if (ifppreg1 .ge. 1) then
                  pot(j,1,ii) = vel(1)
                  pot(j,2,ii) = vel(2)
                  pot(j,3,ii) = vel(3)
               endif
               if (ifppreg1 .ge. 2) then
                  pre(j,ii) = press
               endif
               if (ifppreg1 .ge. 3) then
                  grad(j,1,1,ii) = velgrad(1,1)
                  grad(j,2,1,ii) = velgrad(2,1)
                  grad(j,3,1,ii) = velgrad(3,1)
                  grad(j,1,2,ii) = velgrad(1,2)
                  grad(j,2,2,ii) = velgrad(2,2)
                  grad(j,3,2,ii) = velgrad(3,2)
                  grad(j,1,3,ii) = velgrad(1,3)
                  grad(j,2,3,ii) = velgrad(2,3)
                  grad(j,3,3,ii) = velgrad(3,3)
               endif
            else
               if (ifppreg1 .ge. 1) then
                  pottarg(j,1,ii) = vel(1)
                  pottarg(j,2,ii) = vel(2)
                  pottarg(j,3,ii) = vel(3)
               endif
               if (ifppreg1 .ge. 2) then
                  pretarg(j,ii) = press
               endif
               if (ifppreg1 .ge. 3) then
                  gradtarg(j,1,1,ii) = velgrad(1,1)
                  gradtarg(j,2,1,ii) = velgrad(2,1)
                  gradtarg(j,3,1,ii) = velgrad(3,1)
                  gradtarg(j,1,2,ii) = velgrad(1,2)
                  gradtarg(j,2,2,ii) = velgrad(2,2)
                  gradtarg(j,3,2,ii) = velgrad(3,2)
                  gradtarg(j,1,3,ii) = velgrad(1,3)
                  gradtarg(j,2,3,ii) = velgrad(2,3)
                  gradtarg(j,3,3,ii) = velgrad(3,3)
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

