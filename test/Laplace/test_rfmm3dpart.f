      implicit none
      integer ns,nt
      double precision, allocatable :: source(:,:),targ(:,:)
      double precision, allocatable :: charge(:)
      double precision, allocatable :: dipvec(:,:)
      double precision, allocatable :: pot(:),pottarg(:)
      double precision, allocatable :: grad(:,:),gradtarg(:,:)

      double precision eps
      double complex eye
      integer i,j,k,ntest
      integer ifcharge,ifdipole,ifpgh,ifpghtarg
      double precision err,hkrand
      

      data eye/(0.0d0,1.0d0)/

c
cc      initialize printing routine
c
      call prini(6,13)

      ns = 2000
      nt = 1999

      ntest = 10

      allocate(source(3,ns),targ(3,nt))
      allocate(charge(ns),dipvec(3,ns))
      allocate(pot(ns))
      allocate(grad(3,ns))

      allocate(pottarg(nt))
      allocate(gradtarg(3,nt))


c
cc      generate sources uniformly in the unit cube 
c
c
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


c
cc      generate targets uniformly in the unit cube
c
      do i=1,nt
        targ(1,i) = hkrand(0)
        targ(2,i) = hkrand(0)
        targ(3,i) = hkrand(0)

        pottarg(i) = 0
        gradtarg(1,i) = 0
        gradtarg(2,i) = 0
        gradtarg(3,i) = 0 
      enddo

      eps = 0.5d-4

c
cc     now test source to source, charge, 
c      with potentials
c
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: charges'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstoscp(eps,ns,source,charge,pot)


       ifcharge = 1
       ifdipole = 0
       ifpgh = 1
       ifpghtarg = 0

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'


c

c
cc     now test source to source, charge, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: charges'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstoscg(eps,ns,source,charge,
     1      pot,grad)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 2
       ifpghtarg = 0

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      


c
cc     now test source to source, dipole, 
c      with potentials
c
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstosdp(eps,ns,source,dipvec,
     1      pot)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 0

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'


c
cc     now test source to source, dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstosdg(eps,ns,source,dipvec,
     1      pot,grad)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 0

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'

c
cc     now test source to source, charge + dipole, 
c      with potentials
c
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstoscdp(eps,ns,source,charge,dipvec,
     1      pot)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 0 

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'


c
cc     now test source to source, charge + dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstoscdg(eps,ns,source,charge,dipvec,
     1      pot,grad)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 0

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'



c
cc     now test source to target, charge, 
c      with potentials
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: charges'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstotcp(eps,ns,source,charge,
     1      nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 0
       ifpghtarg = 1

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'

       stop


c
cc     now test source to target, charge, 
c      with potentials and gradients
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: charges'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstotcg(eps,ns,source,charge,
     1      nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 0
       ifpghtarg = 2

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      


c
cc     now test source to target, dipole, 
c      with potentials
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstotdp(eps,ns,source,dipvec,
     1      nt,targ,pottarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 1

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'


c
cc     now test source to target, dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstotdg(eps,ns,source,dipvec,
     1      nt,targ,pottarg,gradtarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 2

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'

c
cc     now test source to target, charge + dipole, 
c      with potentials
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstotcdp(eps,ns,source,charge,dipvec,
     1      nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 1

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'


c
cc     now test source to target, charge + dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to target'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstotcdg(eps,ns,source,charge,dipvec,
     1      nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 0
       ifpghtarg = 2

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'

c
cc     now test source to source + target, charge, 
c      with potentials
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: charges'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstostcp(eps,ns,source,charge,
     1      pot,nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 1
       ifpghtarg = 1

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'


c
cc     now test source to source + target, charge, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstostcg(eps,ns,source,charge,
     1      pot,grad,nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 0
       ifpgh = 2
       ifpghtarg = 2

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'
      


c
cc     now test source to source + target, dipole, 
c      with potentials
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstostdp(eps,ns,source,dipvec,
     1      pot,nt,targ,pottarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 1

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'


c
cc     now test source to source + target, dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstostdg(eps,ns,source,dipvec,
     1      pot,grad,nt,targ,pottarg,gradtarg)

       ifcharge = 0
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 2

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'

c
cc     now test source to source + target, charge + dipole, 
c      with potentials
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstostcdp(eps,ns,source,charge,dipvec,
     1      pot,nt,targ,pottarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 1
       ifpghtarg = 1

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'


c
cc     now test source to source + target, charge + dipole, 
c      with potentials and gradients
c
       write(6,*) 'testing source to source and target'
       write(6,*) 'interaction: charges + dipoles'
       write(6,*) 'output: potentials + gradients'
       write(6,*) 
       write(6,*) 

       call rfmm3dpartstostcdg(eps,ns,source,charge,dipvec,
     1      pot,grad,nt,targ,pottarg,gradtarg)

       ifcharge = 1
       ifdipole = 1
       ifpgh = 2
       ifpghtarg = 2

       
       call comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

       call prin2('l2 rel err=*',err,1)
       write(6,*)
       write(6,*)
       write(6,*) '================'


      stop
      end
c----------------------------------------------------------
c
cc
c
c
c
c
      subroutine comperr(ns,source,ifcharge,charge,ifdipole,
     1   dipvec,ifpgh,pot,grad,nt,targ,ifpghtarg,pottarg,gradtarg,
     2   ntest,err)

      implicit none
      double complex zk
      integer ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg
      
      double precision source(3,*),targ(3,*)
      double precision dipvec(3,*)
      double precision charge(*)

      double precision pot(*),pottarg(*),grad(3,*),gradtarg(3,*)

      integer i,j,ntest,nd

      double precision err,ra
      
      double precision potex(ntest),gradex(3,ntest),pottargex(ntest),
     1                  gradtargex(3,ntest)

      double precision thresh

      nd = 1
      err = 0 
      do i=1,ntest
        potex(i) = 0
        pottargex(i) = 0

        gradex(1,i) = 0
        gradex(2,i) = 0
        gradex(3,i) = 0

        gradtargex(1,i) = 0
        gradtargex(2,i) = 0
        gradtargex(3,i) = 0
      enddo

      thresh = 1.0d-16

      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        if(ifpgh.eq.1) then
          call l3ddirectcp(nd,source,charge,ns,source,ntest,
     1       potex,thresh)
        endif

        if(ifpgh.eq.2) then
          call l3ddirectcg(nd,source,charge,ns,source,ntest,
     1       potex,gradex,thresh)
        endif

        if(ifpghtarg.eq.1) then
          call l3ddirectcp(nd,source,charge,ns,targ,ntest,
     1       pottargex,thresh)
        endif

        if(ifpghtarg.eq.2) then
          call l3ddirectcg(nd,source,charge,ns,targ,ntest,
     1       pottargex,gradtargex,thresh)
        endif
      endif

      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        if(ifpgh.eq.1) then
          call l3ddirectdp(nd,source,dipvec,
     1      ns,source,ntest,potex,thresh)
        endif

        if(ifpgh.eq.2) then
          call l3ddirectdg(nd,source,dipvec,
     1       ns,source,ntest,potex,gradex,thresh)
        endif

        if(ifpghtarg.eq.1) then
          call l3ddirectdp(nd,source,dipvec,
     1       ns,targ,ntest,pottargex,thresh)
        endif

        if(ifpghtarg.eq.2) then
          call l3ddirectdg(nd,source,dipvec,
     1       ns,targ,ntest,pottargex,gradtargex,thresh)
        endif
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        if(ifpgh.eq.1) then
          call l3ddirectcdp(nd,source,charge,dipvec,
     1      ns,source,ntest,potex,thresh)
        endif

        if(ifpgh.eq.2) then
          call l3ddirectcdg(nd,source,charge,dipvec,
     1       ns,source,ntest,potex,gradex,thresh)
        endif

        if(ifpghtarg.eq.1) then
          call l3ddirectcdp(nd,source,charge,dipvec,
     1       ns,targ,ntest,pottargex,thresh)
        endif

        if(ifpghtarg.eq.2) then
          call l3ddirectcdg(nd,source,charge,dipvec,
     1       ns,targ,ntest,pottargex,gradtargex,thresh)
        endif
      endif

      err = 0
      ra = 0

      if(ifpgh.eq.1) then
        do i=1,ntest
          ra = ra + abs(potex(i))**2
          err = err + abs(pot(i)-potex(i))**2
        enddo
      endif

      if(ifpgh.eq.2) then
        do i=1,ntest
          ra = ra + abs(potex(i))**2
          ra = ra + abs(gradex(1,i))**2
          ra = ra + abs(gradex(2,i))**2
          ra = ra + abs(gradex(3,i))**2

          err = err + abs(pot(i)-potex(i))**2
          err = err + abs(grad(1,i)-gradex(1,i))**2
          err = err + abs(grad(2,i)-gradex(2,i))**2
          err = err + abs(grad(3,i)-gradex(3,i))**2
        enddo
      endif


      if(ifpghtarg.eq.1) then
        do i=1,ntest
          ra = ra + abs(pottargex(i))**2
          err = err + abs(pottarg(i)-pottargex(i))**2
        enddo
      endif

      if(ifpghtarg.eq.2) then
        do i=1,ntest
          ra = ra + abs(pottargex(i))**2
          ra = ra + abs(gradtargex(1,i))**2
          ra = ra + abs(gradtargex(2,i))**2
          ra = ra + abs(gradtargex(3,i))**2

          err = err + abs(pottarg(i)-pottargex(i))**2
          err = err + abs(gradtarg(1,i)-gradtargex(1,i))**2
          err = err + abs(gradtarg(2,i)-gradtargex(2,i))**2
          err = err + abs(gradtarg(3,i)-gradtargex(3,i))**2
        enddo
      endif

      err = sqrt(err/ra)
      return
      end
      
