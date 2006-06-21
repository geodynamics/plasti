C#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  ne = total number of elements
c  nn = total number of nodes
c  coord(2,nn) - x,z, coordinates
c  node(j,ne)  - j=1,2,3,4,5,6 node numbers for ith element
c                j=7     code number for viscosity
c                j=8     code number for solid compressibility
c                j=9     code number for density
c  velx(nn)    - x velocity
c  vely(nn)    - y velocity
c  sbar(nvert) - mean stress
c  den(ne)    - density
c
c
c       boundary conditions:
c
c  bvel        - velocities of boundary nodes
c  nvnd        - node numbers of boundary nodes, component codes
c  bp          - pressure at boundary node
c  npnd        - node numbers of pressure boundary node
C#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Copyright (C) 1995-2006 Sean Willett, Chris Fuller

C This program is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation; either version 2 of the License, or (at
C your option) any later version.

C This program is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.

C You should have received a copy of the GNU General Public License
C along with this program; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA


C Copyright (C) 1995 Philippe Fullsack
C Permission is hereby granted to use, reproduce, prepare derivative
C works, and to redistribute to others, so long as this original
C copyright notice is retained.

c####################################################################
c define arrays that will be dynamically allocated in subroutines
c .....................................................................
c arrays for mech model
      module dyn_arrays_mech
      real(kind=8),allocatable::den(:),phi(:),coh(:),coord(:,:),
     *bvel(:),bp(:),basvel(:),unvel(:),zeq(:),bside(:),
     *zinit(:),tpoint(:,:),wdinit(:),f1prev(:),f2prev(:),
     *xp1(:),yp1(:),dyinit1(:),xp2(:),yp2(:),dyinit2(:),
     *fnode1(:),fnode2(:),slen1(:),slen2(:),basinfill(:),
     *peakchop(:),wd1prev(:),wd2prev(:),bvelt(:),vely(:),velx(:),
     *sbar(:,:),visc(:,:),bulkmod(:,:),rhs(:),abd(:,:),soln(:),
     *vbound(:),xsur(:,:),stress(:,:),srate(:,:),sprev(:),ziso(:),
     *zinc(:,:),vpower(:,:),theta(:),veros(:,:),temptc(:),toldc(:),
     *cbase(:),vsur(:,:),rsur(:,:),vdiff(:,:),rdiff(:,:),xsurold(:,:),
     *exhum(:),bastrk(:,:),vmin(:),q(:),prex(:),expn(:)
      integer,allocatable::nvnd(:,:),npnd(:),nbase(:),
     *nsnd(:,:),ieletp(:),node(:,:),nbasinfill(:),npeakchop(:),
     *nvtnd(:),ipflag(:,:),ip(:),ibastrk(:,:),ieletpb(:)
      end module dyn_arrays_mech
c .....................................................................
c arrays for thermal model
      module dyn_arrays_therm
      real(kind=8),allocatable::hprod(:),spheat(:),trho(:),
     *tcond(:,:),btem(:),a(:,:),asf(:,:),bsf(:,:),area(:),
     *rhst(:),flux(:)
      integer,allocatable::ntbnd(:),ipt(:),neflux(:,:)
      end module dyn_arrays_therm
c .....................................................................
c arrays for both
      module dyn_arrays
      real (kind=8),allocatable::coordt(:,:),temp(:),vx(:),vz(:),
     *tempt(:),told(:)
      integer, allocatable::nodet(:,:),output_flags(:)
      end module dyn_arrays
c end of definitions
c####################################################################
c####################################################################
      use dyn_arrays
      use dyn_arrays_mech
      use dyn_arrays_therm
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer,allocatable::ldf(:)
      character date*10,time*10,time2*10,time3*10
      character trans*1
      integer ::count,ioutpt=0

C    Random number generator seed
      idum=-2
      call date_and_time(date,time)
      print*,'Real Time Starting Plasti is:  ',time

c  input the bulk of the data
      call input(nn,ne,lbw,numvbn,numpbn,nout,ntsts,ncol,nrow,
     *ndf,minter,lda,miter,toler,nrowt,ncolt,nnt,net,rhof,rhoman,ncom,
     *nsing,numsid,vrig,delt,nbn,npass,npoint,convel,
     *epsinv,nout_t,nlrow,upveln,erosl,erosr,peros,rpow,iunflag,
     *iunbeg,ntherm,dtherm,w_depth,beta,prig,rrig,sload,smomen,xadd,
     *ctoler,wtoler,np1,np2,npad,nsing1,ndom,ntbn,lbwt,ldat,nlcol,
     *slpmax,tmax,sealev,dy_flex_init1,dy_flex_init2,wdepth,nsthick,
     *plthick,numvetbn,leqflag,iplasflg,iblay,iblayt,isedl,isedr,
     *ibasflg,intmrkb,nbastrk,nbastary,nbastind,ninbas,ipkfill,ibasfill,
     *sedmax)
 
c calculate the initial load on each node for flexure calculation in remesh
c	also store the initial water depth at each node for flex. calc.
      call calc_init_force(np1,np2,nrow,npad,nsing,ncol,nsing1,
     *rhoavinitl,rhoavinitr,dy_flex_init1,dy_flex_init2,sealev,nn)
 
c  identify degrees of freedom of each node
      allocate(ldf(nn))
      do in=1,nn
      	call gdf(in,nrow,ldf(in),nod)
      end do

c  Send to thermal to get initial temp
      itst=0
      deltt=delt
      cbase=0
      call thermal(deltt,itst,nnt,net,nout,nrow,nrowt,ncolt,
     *ldat,lbwt,ntbn,ioutpt)
      do i=1,nnt
      	told(i)=tempt(i)
      end do

c  Make a temp array for just the crust
      ncrustbeg=nrowt-nrow
      l=0
      count=0
      do i=ncrustbeg,nnt,nrowt
	  	l=(i-ncrustbeg)/nrowt+1
        	do j=1,nrow
        		k=i+j
        		count=((l-1)*nrow)+j
        		toldc(count)=told(k)
        		temptc(count)=tempt(k)
      		end do
      end do

c calculate initial viscosities/rheologic parameters
c     rheology parameters
      do 35 ie=1,ne
      	bl=dexp(q(ie)/(8.3144d0*tmax))/(vmin(ie))
      	tele=(temptc(node(1,ie))+temptc(node(2,ie))+temptc(node(3,ie))
     *	+temptc(node(4,ie)))/4.
c     POWER-LAW VISCOSITY
      	if(leqflag.ne.1) then
      		vpow2=(prex(ie)*dexp(-q(ie)/(8.3144d0*tele)))**(-1./expn(ie))
      		vpow=vpow2*(epsinv**(1./expn(ie)-1.))
      	else	
c     LINEAR VISCOUS
     		vpow=dexp(q(ie)/(8.3144d0*tele))/bl
     		vpow2=vpow
     	endif	
      if(vpow.gt.vrig)vpow=vrig                                         
      if(vpow.lt.vmin(ie))vpow=vmin(ie)
      vpower(1,ie)=vpow                                                 
      vpower(2,ie)=vpow2                                                
c     OVERWRITE FOR PLASTIC CASE
      if(iplasflg.eq.1) then
      	vpower(1,ie)=vrig
      	vpower(2,ie)=vrig
      endif

   35 continue
      do i=1,ne
      	do k=1,4
      		visc(i,k)=vpower(1,i)
      		bulkmod(i,k)=1.0/beta
      		ipflag(i,k)=0
      	end do
      end do	
      nb=ne/(2*nrow-2)
      nb=nb+nrow-1
 
c initialize time 
      ttime=0.0
 
c initialize ridge and valley profiles
      do i= 1,nn/nrow
      	rsur(1,i)=coord(1,i*nrow-1)
      	vsur(1,i)=coord(1,i*nrow)
      	rsur(2,i)=coord(2,i*nrow)
      	vsur(2,i)=coord(2,i*nrow)
      	rdiff(1,i)=0.0;
      	vdiff(1,i)=0.0;
      	rdiff(2,i)=0.0;
      	vdiff(2,i)=0.0;
      end do	
 
c initialize stress field
      call sinit(nrow,ne)
 
c loop over time steps
      do 500 itst=1,ntsts
      call date_and_time(date,time2)
      print*,'Real Time in Plasti Loop is:  ',time2

      ttime=ttime+delt
 
c iterate for nonlinearity 
      do 400 iter=1,miter
 
c assemble global stiffness matrix 
      call globe(ne,nn,lbw,delt,lda,ndf,nrow,ldf,c,beta,
     *vrig,sigav,epsav,itst)
 
c  determine the region for underplating to occur
      call unplate(nrow,ncol,itst,nsing,ibegup,ibegmx)

c apply boundary conditions
      call bc(numvbn,numpbn,ndf,lbw,lda,nrow,numsid,
     *nbn,upveln,nsing,ibegup,delt,rhoman,numvetbn,ncol)
 
c apply underplating velocity in z-dir to thermal vel field
      call unplate_therm(nbn,nrowt,nrow)

c LAPACK routine
      call dgbtrf(ndf,ndf,lbw,lbw,abd,lda,ip,info)
      if(info.ne.0) then
	  	print*,'#####  ERROR IN FACTORIZATION, PLASTI DGBTRF'
	  	print*,'info =',info
	  	stop
	  endif
      write(6,603)itst,iter
  603 format('       tstep ',i5,' iteration ',i5)
 
c  solve system of equations
	  trans='N'
	  call dgbtrs(trans,ndf,lbw,lbw,1,abd,lda,ip,rhs,ndf,info)
	  if(info.ne.0) then
	  	print*,'#####  ERROR IN FACTORIZATION, PLASTI DGBTRS'
	  	print*,'info =',info
	  	stop
	  endif
      do idf=1,ndf
      	soln(idf)=rhs(idf)
      end do
 
c  check for convergence
      call conver(velx,vely,soln,toler,icflag,nrow,ncol
     *,ndf,nn,coord)

      if(iter.lt.minter)icflag=0
      if(miter.eq.1)icflag=1
c
c   filter pressure field
c
      call pfilt(ne,nrow,npass)
c
c   calculate stresses and strain rates
c
      call ss(ne,nn,lbw,delt,ndf,nrow,c,beta,vrig)
c
c  branch off after convergence
c
      if(icflag.eq.1)go to 450
  400 continue
  
c #############################
c  failed to converge, so exit
c #############################
      print*,' timestep failed to converge, iterations = ',miter
      nout=1

c   filter pressure field
      call pfilt(ne,nrow,npass)
      call output (nn,ne,itst,iter,nout,ttime,nout_t,nrow,nbn,
     *vrig,tstart,npoint,convel,ntsts,delt,nlrow,sealev,w_depth,
     *nbastrk,ibasflg,ninbas,ioutpt)
      deltt=delt
      call thermal(deltt,itst,nnt,net,nout,nrow,nrowt,ncolt,
     *ldat,lbwt,ntbn,ioutpt)
      stop

c ###########
c converged
c ###########
  450 continue

      upmass=(coord(1,(nsing-1)*nrow+1)-coord(1,ibegup*nrow+1))
     **upveln*ttime/1000000. 
      acmass=ttime*convel*(coord(2,nrow)-coord(2,1))/1000000.

      deltol=delt

c  track material points
      print*,'Entering Lagrangian Mesh Tracking Routine'
      call date_and_time(date,time2)
c      print*,'Time before L-mesh tracki:  ',time2
      call tracki(nn,ne,nrow,npoint,delt,itst)
      call date_and_time(date,time2)
c      print*,'Time after L-mesh tracki:  ',time2
      if(ibasflg.eq.1) then
c      	print*,'Entering Basin Tracking Routine'
      	call date_and_time(date,time2)
c      	print*,'Time before Basin track:  ',time2
      	call track_basin(nn,ne,nrow,ninbas,delt,itst)
      	call date_and_time(date,time2)
c      	print*,'Time after Basin track:  ',time2
      endif	

c   filter pressure field
      call pfilt(ne,nrow,npass)
      print*,'Entering output'
      call output (nn,ne,itst,iter,nout,ttime,nout_t,nrow,nbn,
     *vrig,tstart,npoint,convel,ntsts,delt,nlrow,sealev,w_depth,
     *nbastrk,ibasflg,ninbas,ioutpt)

      dum=delt
      write(6,602)dum
  602 format('  next time step = ',e12.6,' my')
 
c  CALL THERMAL
c      print*, 'time step before thermal is:', itst
      deltt=delt
      call thermal(deltt,itst,nnt,net,nout,nrow,nrowt,ncolt,
     *ldat,lbwt,ntbn,ioutpt)
 
c   remesh
      print*,'Entering remesh'
      call remesh(deltol,nn,nrow,rhoman,rhof,ncom,ne,vrig,npoint,
     *nnt,nrowt,net,nsing,convel,itst,erosl,erosr,
     *peros,rpow,sealev,slpmax,npad,np1,np2,prig,rrig,nsing1,ctoler,
     *wtoler,tmax,rhoavinitl,rhoavinitr,dy_flex_init1,dy_flex_init2,
     *wdepth,nsthick,plthick,leqflag,iplasflg,iblay,iblayt,isedl,isedr,
     *ibasflg,intmrkb,nbastrk,nbastary,nbastind,ninbas,ipkfill,ibasfill,
     *sedmax)
c
c reset rigid viscosity to power law viscosity
      do i=1,ne
      	do k=1,4
      		if(ipflag(i,k).eq.0)then
      			visc(i,k)=vpower(1,i)
      		endif
      	end do
      	
      end do

C Allow thermal conditions to evolve for sub zone w/o collision
      if(ntherm.gt.0) then
      	if(itst.eq.1) then
      		print*,'thermal runup'
      		call therm_runup(deltt2,itst2,nnt,net,nout,nrow,nrowt,
     *		ncolt,ldat,lbwt,ntbn,dtherm,ntherm)
      	endif
      endif

  500 continue
      stop
      end
c#CCCCCCCCCCCCCCCCCCCCCCCC
c                        C
c  END OF MAIN PROGRAM   C
c                        C
c#CCCCCCCCCCCCCCCCCCCCCCCC

c#########################################################
c thermal runup 
c#########################################################
      subroutine therm_runup(deltt2,itst2,nnt,net,nout,nrow,nrowt,
     *ncolt,ldat,lbwt,ntbn,dtherm,ntherm)

      use dyn_arrays
      use dyn_arrays_therm
      use dyn_arrays_mech
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real(kind=8),allocatable::vx2(:),vz2(:)

        allocate(vx2(net),vz2(net))
        itst2=99999

c set crustal velocities to zero for thermal equilibation        
      	nlthick=nrowt-nrow
      	vx2=vx
      	vy2=vy
      	do itimes=2*nlthick,net,2*(nrowt-1)
      		do jtimes=1,2*(nrow-1)
      			vx(itimes+jtimes)=0.0
      			vz(itimes+jtimes)=0.0
            end do
       	end do

      	do itherm=1,ntherm
      		if(itherm.gt.1) then
      			do inodes=1,nnt
        			told(inodes)=tempt(inodes)
                end do       			
      		endif
        	deltt2=dtherm
        	print*,'Thermal Setup loop:',itherm
        	call thermal(deltt2,itst2,nnt,net,nout,nrow,nrowt,
     *   	ncolt,ldat,lbwt,ntbn,ioutpt)
      	end do
      	vx=vx2
      	vy=vy2
      	deallocate(vx2,vz2)
       
       end
      subroutine calc_flex_remesh(nrow,ncol,np1,np2,prig,rrig,rhom,
     *npad,nsing,nsing1,ctoler,sload,smomen,wtoler,rhof,g,rhoavinitl,
     *rhoavinitr,dy_flex_init1,dy_flex_init2,wdepth,sealev,nn,itst)

      use dyn_arrays_mech

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  NOTE: there is no rhof in the def of alpha.
c 	see project notes for description of method of calculating flexure and
c	why the rhof is left out. In brief, it is left out because the forces
c	acting on the imaginary plate are calculated as the load from the crust
c	and the load from the water.  Another way to do this problem would be 
c	to use rhof in the eqn and calculate just the loads from the crust.
c	In this case the force from any portion of a colm of crust that is below a 
c	defined sea level is (rhoc-rhof)*g*h and the force form the portion of 
c	the same colm above sea level (if there is a sub aerial portion) is 
c	rhoc*g*h', where h' is the height of the colm above sea level
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c calculate flexural parameters
      alpha1=(4.0*prig/((rhom)*g))**0.25
      alpha2=(4.0*rrig/((rhom)*g))**0.25
      plam1=1.0/alpha1
      plam2=1.0/alpha2
      fk=rhom*g

c####################################################################
c caculate the change in force on each node from previous timestep
c####################################################################
c initial average density of entire model domain is
c	 used for the density of model in the padded regions
      rhoavtl=rhoavinitl
      rhoavtr=rhoavinitr

c total force from crustal loads
c 	plate 1
      call calc_force_p1(slen1,xp1,nrow,den,np1,fnode1,
     *npad,rhoavtl,nsing,g,coord,ncol,nsing1,dy_flex_init1)
c 	plate 2
      call calc_force_p2(slen2,xp2,nrow,den,np2,fnode2,
     *npad,rhoavtr,nsing,g,coord,ncol,dy_flex_init2)
c change in force
      fnode1=fnode1-f1prev
      fnode2=fnode2-f2prev

c calculate defection of plates from change in distributed load
c 	and coupling load (sub momment and load are not reapplied)
      call deflect(np1,np2,xp1,xp2,yp1,yp2,fnode1,fnode2,plam1,plam2,
     *fk,ctoler,nsing1,sload,smomen,xbase,nsing,npad)
c store force used in flexure calculation for det. change in force
c	at the next timestep.
      f1prev=fnode1+f1prev
      f2prev=fnode2+f2prev     

c calculate the deflection from the load of overlying water
      if(wdepth.gt.0.0) then
     	call deflectw(wdepth,xp1,xp2,slen1,slen2,wtoler,fnode1,
     *	fnode2,plam1,plam2,fk,rhof,np1,np2,npad,g,nsing1,ctoler,
     *	dyinit1,dyinit2,yp2,yp1,sealev,nrow,wd1prev,wd2prev,
     *	nsing,nn,coord)
      	dif=0.0
      endif
      end
c#########################################################
C calclate the additional deflection of the coupled plates
c	from the overlying load of water
c#########################################################
      subroutine deflectw(wdepth,xp1,xp2,slen1,slen2,wtoler,fnode1
     *,fnode2,plam1,plam2,fk,rhof,np1,np2,npad,g,nsing1,ctoler,
     *dyinit1,dyinit2,yp2,yp1,sealev,nrow,wd1prev,wd2prev,nsing,nn,
     *coord)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 xp1(*),xp2(*),yp1(np1),yp2(np2),fnode1(np1),fnode2(np2),
     *slen1(*),slen2(*),dyinit1(*),dyinit2(*),coord(2,*),wd1prev(np1),
     *wd2prev(np2)
      real(kind=8),allocatable::dloc1(:),dloc2(:),yp1pre(:),yp2pre(:)
      real(kind=8),allocatable::wd1cur(:),wd2cur(:)

      ychange=100.0*wtoler
      icount=0
      allocate(dloc1(np1),dloc2(np2),yp1pre(np1),yp2pre(np2))
      allocate(wd1cur(np1),wd2cur(np2))
      dloc1=0.0
      yp1pre=yp1
      fnode1=0.0
      wd1cur=0.0
      dloc2=0.0
      yp2pre=yp2
      fnode2=0.0
      wd2cur=0.0

c initerate for convergence on water flexure      
      do while(ychange.gt.wtoler)
      	icount=icount+1
c 	calculate curent water depths
      	call calc_wd(nsing1,npad,np1,sealev,nrow,wd1cur,
     *	wd2cur,nsing,np2,nn,coord,yp1,yp2)
c	calculate force from change in water depth
c	plate 1      		
      	do i=1,np1
c			change in water depth
      		deltad=wd1cur(i)-wd1prev(i)
      		if(wd1cur(i).le.0.0) then
      			if(wd1prev(i).le.0.0) then
      				fnode1(i)=0.0
      			else
      				fnode1(i)=slen1(i)*g*rhof*(-wd1prev(i))
      			endif
      		else
      			fnode1(i)=slen1(i)*g*rhof*deltad
      		endif
      	end	do
c	plate 2      		
      	do i=1,np2
c			change in water depth
      		deltad=wd2cur(i)-wd2prev(i)
      		if(wd2cur(i).le.0.0) then
      			if(wd2prev(i).le.0.0) then
      				fnode2(i)=0.0
      			else
      				fnode2(i)=slen2(i)*g*rhof*(-wd2prev(i))
      			endif
      		else
      			fnode2(i)=slen2(i)*g*rhof*deltad
      		endif
      	end	do

c############################
c calculate plate deflection
c############################
c calculate deflection,moment,shear force at the desired break point for 
c two infinite plates
c     	plate 1
      	amom1=0.0
      	ashear1=0.0
      	call calc_dms(xp1,amom1,ashear1,np1,fnode1,plam1,yp1,fk)
c     	plate 2
      	amom2=0.0
      	ashear2=0.0
      	call calc_dms(xp2,amom2,ashear2,np2,fnode2,plam2,yp2,fk)
c calculate deflection of semi-infinite plates using 
c	the end cond forces and subduction load/moment
c     	plate 1
      	call deflect2(np1,plam1,fk,xp1,sload,smomen,amom1,ashear1,yp1)
c     	plate 2
      	call deflect2(np2,plam2,fk,xp2,sload,smomen,amom2,ashear2,yp2)
c calculate the coupling load
      	ido_again=1
      	jcount=0
      	do while(ido_again==1) 
      		jcount=jcount+1
                call couple(yp1,yp2,nsing,plam1,plam2,fk,xp1,xp2,np1,np2
     $               ,nsing1)
      		if(abs(yp1(nsing1)-yp2(1)).le.ctoler) then
      			ido_again=0
      		else if(jcount.gt.100) then
      			ido_again=0
      			print*,'########################################'
      			print*,'## coupling iteration exceeded 100    ##'
      			print*,'##     inside water loop              ##'
      			print*,'########################################'
      			call profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
      		endif	
      	end do
      	dif=0.0
      	do j=1,np1
      		dif=dif+abs(yp1(j)-yp1pre(j))
      		yp1pre(j)=yp1(j)
      	end do
      	do j=1,np2
      		dif=dif+abs(yp2(j)-yp2pre(j))
      		yp2pre(j)=yp2(j)
      	end do	
      	ychange=dif
      	if(icount.gt.100) then
      		print*,'########################################'
      		print*,'## water depth iteration exceeded     ##'
      		print*,'##	iterations:',icount
      		print*,'##	diff. in plate position at spoint (m):',
     *		abs(yp1(nsing1)-yp2(1))
      		print*,'##	change in base elev:',dif
      		print*,'########################################'
      		call profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
      	endif	
c      	store water depths for next iteration/timestep      	
      	wd1prev=wd1cur
      	wd2prev=wd2cur
      end do
      print*,'Water depth coupling'
      print*,'  iterations:',icount
      print*,'  diff. in plate position at spoint (m):',
     *abs(yp1(nsing1)-yp2(1))
      print*,'  change in base elev:',dif
      deallocate(dloc1,dloc2,yp1pre,yp2pre,wd1cur,wd2cur)
      end

c########################################################
c output the flexural profiles when then code dumps
c########################################################
      subroutine profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
      implicit integer (i-n)
      implicit real (a-h,o-z)
      real(kind=8) xbase(*),yp1(*),yp2(*)

      open(21,file='profiles/pro_plate_dump')
      open(22,file='profiles/retro_plate_dump')
      do k=1,np1
      	ip=np1-k+1
	  	write(21,198)xbase(k)/1000.0,-yp1(ip)/1000.0
	  end do
	  do k=1,np2
      	write(22,198)xbase(nsing-1+k+npad)/1000.0,-yp2(k)/1000.0
      end do
  198 format(2e17.8)
	  close(21)
	  close(22)
      stop
      end 
      
c########################################################
C calculate the deflection of two semi-infinite plates 
c	coupled together at the s-point from a distributed
c	load as stored in fnode
c########################################################
      subroutine deflect(np1,np2,xp1,xp2,yp1,yp2,fnode1,fnode2
     *,plam1,plam2,fk,ctoler,nsing1,sload,smomen,xbase,nsing,npad)

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 xp1(*),xp2(*),yp1(*),yp2(*),fnode1(*),fnode2(*),xbase(*)

c calculate deflection,moment,shear force at the desired break point for 
c two infinite plates
c     plate 1
      amom1=0.0
      ashear1=0.0
      call calc_dms(xp1,amom1,ashear1,np1,fnode1,plam1,yp1,fk)
c     plate 2
      amom2=0.0
      ashear2=0.0
      call calc_dms(xp2,amom2,ashear2,np2,fnode2,plam2,yp2,fk)

c calculate deflection of semi-infinite plates using 
c	the end cond forces and subduction load/moment
c     plate 1
c      print*,'PLATE 1',yp1(1),yp1(100)
      call deflect2(np1,plam1,fk,xp1,sload,smomen,amom1,ashear1,yp1)
c     plate 2
      call deflect2(np2,plam2,fk,xp2,sload,smomen,amom2,ashear2,yp2)
c calculate the coupling load
      ido_again=1
      icount=0
      do while(ido_again==1) 
      	icount=icount+1
		call couple(yp1,yp2,nsing,plam1,plam2,fk,xp1,xp2,np1,np2,nsing1)
      	if(abs(yp1(nsing1)-yp2(1)).le.ctoler) then
      		ido_again=0
      	else if(icount.gt.100) then
      		ido_again=0
      		print*,'########################################'
      		print*,'## coupling iteration exceeded 100    ##'
      		print*,'##         first loop                 ##'
      		print*,'########################################'
      		call profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
      	endif	
      end do
      print*,'Plate coupling:'
      print*,'  iterations:',icount
      print*,'  diff. at s-point: ',abs(yp1(nsing1)-yp2(1))
      end

C########################################################
c calculate the plate ocupling load
c########################################################
      subroutine couple(yp1,yp2,nsing,plam1,plam2,fk,xp1,xp2,
     *np1,np2,nsing1)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 yp1(*),yp2(*),xp1(*),xp2(*)

c difference in deflection between the two plates at s-point
	  G0=yp1(nsing1)-yp2(1)
c calculations for plate 1	  
c deflection of infinite beam with the coupling load at s-point

      G1=plam1/(2.0*fk)*exp(-plam1*xp1(nsing1))*(cos(plam1*xp1(nsing1))
     *+sin(plam1*xp1(nsing1)))
c moment at the plate end from coupling load
      G2=1.0/(4.0*plam1)*exp(-plam1*xp1(nsing1))*(cos(plam1*xp1(nsing1))
     *-sin(plam1*xp1(nsing1)))
c shear force at the plate end from coupling load
      G3=0.5*exp(-plam1*xp1(nsing1))*cos(plam1*xp1(nsing1))
c end conditioning load
      G4=4.0*plam1*G2+4.0*G3
c end conditioning moment
      G5=-4.0*G2-2.0*G3/plam1
c deflection from end conditioning load
      G6=G4*plam1/(2.0*fk)*exp(-plam1*xp1(nsing1))*
     *(cos(plam1*xp1(nsing1))+sin(plam1*xp1(nsing1)))
c deflection from end conditioning moment
      G7=G5*plam1**2/fk*exp(-plam1*xp1(nsing1))*sin(plam1*xp1(nsing1))
      fcouple=G0/((2.0*plam2)/fk+G1+G6+G7)
c calculate deflection from coupling load
      do i=1,np1
      	yp1(i)=yp1(i)
     *	-(2.0*fcouple*plam1/fk*exp(-plam1*xp1(i))*cos(plam1*xp1(i)))
      end do
      do i=1,np2
      	yp2(i)=yp2(i)
     *	+(2.0*fcouple*plam2/fk*exp(-plam2*xp2(i))*cos(plam2*xp2(i)))
      end do
      end

c########################################################
c calculate deflection of semi-infinite plates using 
c	the end cond forces and subduction load/moment
c########################################################
      subroutine deflect2(np,plam,fk,xp,sload,smomen,amom,ashear,yp)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 xp(*),yp(*)
c calculate end conditioning forces
      fmo=-4.0*amom-2.0*ashear/plam
      fpo=4.0*(plam*amom+ashear)
c calculate deflection      
      do i=1,np
      	ypo=fpo*plam/(2.0*fk)*exp(-plam*xp(i))*(cos(plam*xp(i))
     *  +sin(plam*xp(i)))
      	ymo=(fmo*plam**2)/fk*exp(-plam*xp(i))*sin(plam*xp(i))
      	ysload=2.0*sload*plam/fk*exp(-plam*xp(i))*cos(plam*xp(i))
      	ysmom=-2.0*smomen*plam**2/fk*exp(-plam*xp(i))*(cos(plam*xp(i))
     *	-sin(plam*xp(i)))
        yp(i)=ypo+ymo+ysload+ymom+yp(i)
      end do  
      end

c########################################################
c  calculate the moment and shear in an infinite plate 
c########################################################
      subroutine calc_dms(xp,amom,ashear,np,fnode,plam,yp,fk)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 xp(*),fnode(*),yp(*)
      
      do i=1,np
      	do j=1,np
      		dist=abs(xp(i)-xp(j))
      		yp(j)=yp(j)+fnode(i)*plam/(2.0*fk)*exp(-plam*dist)*
     *		(cos(plam*dist)
     *		+sin(plam*dist))
     	end do
      	dist=xp(i)
      	amom=amom+fnode(i)/(4.0*plam)*exp(-plam*dist)*(cos(plam*dist)
     *	-sin(plam*dist))
        ashear=ashear+fnode(i)/2.0*exp(-plam*dist)*cos(plam*dist)
      end do  
      end

c####################################################################
c calculate the initial forces and water depths for flexure problem 
c####################################################################
      subroutine calc_init_force(np1,np2,nrow,npad,nsing,ncol,nsing1,
     *rhoavinitl,rhoavinitr,dy_flex_init1,dy_flex_init2,sealev,nn)

      use dyn_arrays_mech
      implicit real(kind=8) (a-h,o-z)
      implicit integer (i-n)
      integer i,j,irow,icol,np1,np2,nrow,npad,nsing,ncol,nsing1
      real*8 g,alpha1,alpha2,plam1,plam2,fk,height,base,heightl,
     *heightr,area,areat,arealoc,areacheck,rhoavt,areat2

      allocate(f1prev(np1),f2prev(np2),fnode1(np1),fnode2(np2))
      allocate(wd1prev(np1),wd2prev(np2))
      allocate(slen1(np1),slen2(np2))
      f1prev=0.0
      f2prev=0.0
      fnode1=0.0
      fnode2=0.0
      wd1prev=0.0
      wd2prev=0.0
      slen1=0.0
      slen2=0.0

c calculate flexural parameters
      g=9.8
      alpha1=(4.0*prig/((rhom)*g))**0.25
      alpha2=(4.0*rrig/((rhom)*g))**0.25
      plam1=1.0/alpha1
      plam2=1.0/alpha2
      fk=rhom*g
c########################################################
c caculate force on each node from mech. model thickness
c########################################################
c average density of entire model domain
c	 used for the density of model in the padded regions
c	calculate the area of L and R edge colms
      areatl=0.0
      areatr=0.0
c     lower left triangle
      base=abs(coord(2,nrow)-coord(2,1))
      height=abs(coord(1,nrow+1)-coord(1,1))
      areatl=0.5*base*height
c     upper left triangle
      base=abs(coord(2,nrow*2)-coord(2,nrow+1))
      areatl=0.5*base*height+areatl
c     lower right triangle
      base=abs(coord(2,nrow*(ncol-1))-coord(2,nrow*(ncol-2)+1))
      height=abs(coord(1,nrow*(ncol-1)+1)-coord(1,nrow*(ncol-2)+1))
      areatr=0.5*base*height
c     upper right triangle
      base=abs(coord(2,nrow*ncol)-coord(2,nrow*(ncol-1)+1))
      areatr=0.5*base*height+areatr

c determine the average density of R and L edge colm
      areatl2=0.0
      areatr2=0.0
      rhoavtl=0.0
      rhoavtr=0.0
c	left edge
      height=abs(coord(1,2*nrow)-coord(1,nrow))
      do irow=1,nrow-1
c     	lower triangle
      	base=abs(coord(2,1+irow)-coord(2,irow))
      	arealoc=base*height*.5
c     	upper triangle
      	base=abs(coord(2,nrow+irow+1)-coord(2,nrow+irow))
      	arealoc=arealoc+base*height*.5
      	areatl2=areatl2+arealoc
c     	calculate desity contribution      		
      	rhoavtl=rhoavtl+den(irow)*arealoc/areatl
      end do

c	right edge
      height=abs(coord(1,ncol*nrow)-coord(1,(ncol-1)*nrow))
      do irow=1,nrow-1
c     	lower triangle
      	base=abs(coord(2,(ncol-2)*nrow+1+irow)
     *	-coord(2,(ncol-2)*nrow+irow))
      	arealoc=base*height*.5
c     	upper triangle
      	base=abs(coord(2,(ncol-1)*nrow+irow+1)
     *	-coord(2,(ncol-1)*nrow+irow))
      	arealoc=arealoc+base*height*.5
      	areatr2=areatr2+arealoc
c     	calculate desity contribution      		
      	rhoavtr=rhoavtr+den((ncol-2)*(nrow-1)+irow)*arealoc/areatr
      end do

c store initial average denisty of model for use in flexure calculation
c	it is used in the padded regions
      rhoavinitl=rhoavtl
      rhoavinitr=rhoavtr
c      print*,'############'
c      print*,'model density and model area of L colm',rhoavtl,areatl
c      print*,'area check',areatl2
c      print*,'model density and model area of R colm',rhoavtr,areatr
c      print*,'area check',areatr2
c      print*,'############'

c plate 1
      call calc_force_p1(slen1,xp1,nrow,den,np1,fnode1,
     *npad,rhoavtl,nsing,g,coord,ncol,nsing1,dy_flex_init1)
      f1prev=fnode1
c plate 2
      call calc_force_p2(slen2,xp2,nrow,den,np2,fnode2,
     *npad,rhoavtr,nsing,g,coord,ncol,dy_flex_init2)
      f2prev=fnode2

c calculate the initial water depth
      call calc_wd(nsing1,npad,np1,sealev,nrow,wd1prev,
     *wd2prev,nsing,np2,nn,coord,yp1,yp2)

      end
c##########################################################
c calculate the water depth above the flexural profles 
c##########################################################
      subroutine calc_wd(nsing1,npad,np1,sealev,nrow,wd1,
     *wd2,nsing,np2,nn,coord,yp1,yp2)

      implicit real *8 (a-h,o-z)
      implicit integer (i-n)
      dimension wd1(*),wd2(*),coord(2,*),yp1(*),yp2(*)

c calculate initial water depths
      yshift=coord(2,1)+yp1(np1-npad)
c	Plate 1:
c     extended sub plate region
      do i=1,nsing1-1
      	wd1(i)=0.0
      end do	
c     region in mech model
      icount=0
      do i=nsing1,np1-npad
      	icount=icount+1
      	icol=nsing+1-icount
      	wd1(i)=sealev-coord(2,icol*nrow)
      	if(wd1(i).lt.0.0) wd1(i)=0.0
      end do	
c	  padded region
      dy=coord(2,nrow)-coord(2,1)
      do i=np1-npad+1,np1
      	wd1(i)=sealev-(yshift-yp1(i)+dy)
      	if(wd1(i).lt.0.0) wd1(i)=0.0
      end do
c	Plate 2:
c     region in mech model
      do i=1,np2-npad
      	wd2(i)=sealev-coord(2,nsing*nrow+(i-1)*nrow)
      	if(wd2(i).lt.0.0) wd2(i)=0.0
      end do
c     padded region
      dy=coord(2,nn)-coord(2,nn-nrow+1)
      do i=np2-npad+1,np2
      	wd2(i)=sealev-(yshift-yp2(i)+dy)
      	if(wd2(i).lt.0.0) wd2(i)=0.0
      end do	
      end
      

c########################################################
c calculate the force from the thickness of the mech model
c	for calculating the flexure
c########################################################
      subroutine calc_force_p2(slen,xp,nrow,den,np,fnode,
     *npad,rhoavt,nsing,g,coord,ncol,dy_flex_init2)      

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 slen(*),xp(*),den(*),fnode(*),coord(2,*)

      ifstrow=(nrow-1)*(nsing-1)
      ilstrow=(ncol-2)*(nrow-1)

c first node
      slen(1)=(xp(2)-xp(1))/2.0
      rhoav=0.0
      areacheck=0.0
c     calculate the area of the first colm
c     lower triangle
      base=abs(coord(2,nsing*nrow)-coord(2,(nsing-1)*nrow+1))
      height=abs(coord(1,(nsing+1)*nrow)-coord(1,nsing*nrow))
      areacol=.5*base*height
c     upper triangle
      base=abs(coord(2,(nsing+1)*nrow)-coord(2,nsing*nrow+1))
      areacol=areacol+.5*base*height
c     calculate area and density contribution of each elem in colm      
      do irow=1,nrow-1
c     	lower triangle      	
      	base=abs(coord(2,(nsing-1)*nrow+1+irow)
     *	-coord(2,(nsing-1)*nrow+irow))
      	areaelm=base*height*.5
c     	upper triangle
      	base=abs(coord(2,nsing*nrow+irow+1)-coord(2,nsing*nrow+irow))
      	areaelm=areaelm+base*height*.5
c     	contribution to average density
      	rhoav=rhoav+den(irow+ifstrow)*areaelm/areacol
      	areacheck=areacheck+areaelm
      end do
c      print*,'first colm',areacol,areacheck,rhoav
      fnode(1)=slen(1)*g*rhoav*(coord(2,nsing*nrow)
     *-coord(2,(nsing-1)*nrow+1))

c last node
      slen(np)=(xp(np)-xp(np-1))/2.0
      if(npad.eq.0) then
c if there is no padding on model edges      
      	rhoav=0.0
      	areacheck=0.0
c     	calculate area of the last colm
c     	lower triangle
      	base=abs(coord(2,(ncol-1)*nrow)-coord(2,(ncol-2)*nrow+1))
      	height=abs(coord(1,ncol*nrow)-coord(1,(ncol-1)*nrow))
      	areacol=.5*base*height
c      	upper triangle
      	base=abs(coord(2,ncol*nrow)-coord(2,(ncol-1)*nrow+1))
      	areacol=areacol+.5*base*height
c     	calculate area and density contribution of each elem in colm      	
      	do irow=1,nrow-1
c     		lower triangle      	
      		base=abs(coord(2,(ncol-2)*nrow+1+irow)
     *		-coord(2,(ncol-2)*nrow+irow))
      		areaelm=base*height*.5
c     		upper triangle
     		base=abs(coord(2,(ncol-1)*nrow+irow+1)
     *		-coord(2,(ncol-1)*nrow+irow))
      		areaelm=areaelm+base*height*.5
c     		contribution to average density
      		rhoav=rhoav+den(irow+ilstrow)*areaelm/areacol
      		areacheck=areacheck+areaelm
      	end do
      	dy=coord(2,nn)-coord(2,nn+1-nrow)
      	print*,'##################################'
      	print*,'##################################'
      	print*,'last colm',areacol,areacheck,rhoav
      else
c if there is padding on the model edges, use model average density for 
c	the padded regions
      	rhoav=rhoavt
      	dy=dy_flex_init2
      endif	
      fnode(np)=slen(np)*g*rhoav*dy

c all other nodes
      do icol=2,np-1
      	slen(icol)=(xp(icol+1)-xp(icol-1))/2.0
c       add catch for padded edges of model where density is not defined
      	if(icol.ge.np-npad) then
     		rhoav=rhoavt
     		dy=dy_flex_init2
      	else	
	      	rhoav=0.0
	      	areacheck=0.0
	      	icol2=icol+nsing-1
	      	dy=coord(2,icol2*nrow)-coord(2,(icol2-1)*nrow+1)
c      		now there are two colms to calculate areas for	      	
c	     	lower left triangle
      		heightl=abs(coord(1,icol2*nrow)-coord(1,(icol2-1)*nrow))
      		base=abs(coord(2,(icol2-1)*nrow)-coord(2,(icol2-2)*nrow+1))
      		areacol=.5*base*heightl
c	      	upper left triangle
      		base=abs(coord(2,icol2*nrow)-coord(2,(icol2-1)*nrow+1))
      		areacol=areacol+.5*base*heightl
c	     	lower right triangle
      		heightr=abs(coord(1,(icol2+1)*nrow)-coord(1,icol2*nrow))
c     		the base is the same for lright and uleft
      		areacol=areacol+.5*base*heightr
c	      	upper right triangle
      		base=abs(coord(2,(icol2+1)*nrow)-coord(2,icol2*nrow+1))
      		areacol=areacol+.5*base*heightr
      		do irow=1,nrow-1
c	     		lower left triangle
      			base=abs(coord(2,(icol2-2)*nrow+1+irow)
     *			-coord(2,(icol2-2)*nrow+irow))
      			areaelml=base*heightl*.5
c	      		upper left triangle
      			base=abs(coord(2,(icol2-1)*nrow+1+irow)
     *			-coord(2,(icol2-1)*nrow+irow))
      			areaelml=areaelml+.5*base*heightl
c	     		lower right triangle
      			areaelmr=base*heightr*.5
c	      		upper right triangle
      			base=abs(coord(2,icol2*nrow+1+irow)
     *			-coord(2,icol2*nrow+irow))
      			areaelmr=areaelmr+.5*base*heightr
c      			contribution to average density
      			rhoav=rhoav+den((icol2-2)*(nrow-1)+irow)*areaelml
     *			/areacol+den((icol2-1)*(nrow-1)+irow)*areaelmr/areacol
     			areacheck=areacheck+areaelml+areaelmr
      		end do
c      	    print*,'plate 2',abs(areacol-areacheck),rhoav
      	end	if
      	fnode(icol)=slen(icol)*g*rhoav*dy
      end do	
      end


c########################################################
c calculate the force from the thickness of the mech model
c	for calculating the flexure
c########################################################
      subroutine calc_force_p1(slen,xp,nrow,den,np,fnode,
     *npad,rhoavt,nsing,g,coord,ncol,nsing1,dy_flex_init1)

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 slen(*),xp(*),den(*),fnode(*),coord(2,*)

      ifstrow=(nerowm-1)*(np-1)
      ilstrow=0

c if there is an extended area (sub. plate), set all forces there to zero
      slen(1)=(xp(2)-xp(1))/2.0
      fnode(1)=0.0
      do i=2,nsing1-1
      	slen(i)=(xp(i+1)-xp(i-1))/2.0
      	fnode(i)=0.0
      end do	

c calculte the force at the spoint as if it was the first node
      slen(nsing1)=(xp(nsing1+1)-xp(nsing1))/2.0
      rhoav=0.0
      areacheck=0.0
c     calculate area of the colm
c     lower triangle
      base=abs(coord(2,(nsing-1)*nrow)-coord(2,(nsing-2)*nrow+1))
      height=abs(coord(1,nsing*nrow)-coord(1,(nsing-1)*nrow))
      areacol=.5*base*height
c     upper triangle
      base=abs(coord(2,nsing*nrow)-coord(2,(nsing-1)*nrow+1))
      areacol=areacol+.5*base*height
c     calculate area and density contribution of each elem in colm       	
      do irow=1,nrow-1
c     	lower triangle
      	base=abs(coord(2,(nsing-2)*nrow+irow+1)-
     *	coord(2,(nsing-2)*nrow+irow))
      	areaelm=base*height*.5
c      	upper triangle
      	base=abs(coord(2,(nsing-1)*nrow+irow+1)-
     *	coord(2,(nsing-1)*nrow+irow))
      	areaelm=areaelm+base*height*.5
c      	contribution to average density
      	rhoav=rhoav+den((nsing-2)*nrow+irow)*areaelm/areacol
      	areacheck=areacheck+areaelm
      end do
      fnode(nsing1)=slen(nsing1)*g*rhoav*abs(coord(2,nsing*nrow)
     *-coord(2,(nsing-1)*nrow+1))
c      print*,'spoint colm plate 1',areacol,areacheck,rhoav
c last node
      slen(np)=(xp(np)-xp(np-1))/2.0
      if(npad.eq.0) then
c if there is no padding on model edges      
      	rhoav=0.0
      	areacheck=0.0
c      	calculate area of last colm
c     	lower triangle
      	base=abs(coord(2,nrow)-coord(2,1))
      	height=abs(coord(1,nrow+1)-coord(1,1))
      	areacol=.5*base*height
c     	upper triangle
      	base=abs(coord(2,nrow*2)-coord(2,nrow+1))
      	areacol=areacol+.5*base*height
c     	calculate area and density contribution of each elem in colm      	
      	do irow=1,nrow-1
c     		lower triangle
     		base=abs(coord(2,irow+1)-coord(2,irow))
     		areaelm=.5*base*height
c     		upper triangle
      		base=abs(coord(2,nrow+1+i)-coord(2,nrow+i))
      		areaelm=areaelm+.5*base*height
c     		contribution to average density
      		rhoav=rhoav+den(irow)*areaelm/areacol
      		areacheck=areacheck+areaelm
      	end do
      	dy=coord(2,nrow)-coord(2,1)
      else
      	rhoav=rhoavt
      	dy=dy_flex_init1
      endif	
      fnode(np)=slen(np)*g*rhoav*dy


c all other nodes
      index=0
      do icol=nsing1+1,np-1
      	slen(icol)=(xp(icol+1)-xp(icol-1))/2.0
c       add catch for padded edges of model where density is not defined
      	if(icol.ge.np-npad) then
     		rhoav=rhoavt
     		dy=dy_flex_init1
      	else	
      		rhoav=0.0
      		areacheck=0.0
      		index=index+1
      		icol2=nsing-index
      		dy=coord(2,icol2*nrow)-coord(2,(icol2-1)*nrow+1)
c      		now there are two colms to calculate areas for	      	
c	     	lower left triangle
      		heightl=abs(coord(1,icol2*nrow)-coord(1,(icol2-1)*nrow))
      		base=abs(coord(2,(icol2-1)*nrow)-coord(2,(icol2-2)*nrow+1))
      		areacol=.5*base*heightl
c	      	upper left triangle
      		base=abs(coord(2,icol2*nrow)-coord(2,(icol2-1)*nrow+1))
      		areacol=areacol+.5*base*heightl
c	     	lower right triangle
      		heightr=abs(coord(1,(icol2+1)*nrow)-coord(1,icol2*nrow))
c     		the base is the same for lright and uleft
      		areacol=areacol+.5*base*heightr
c	      	upper right triangle
      		base=abs(coord(2,(icol2+1)*nrow)-coord(2,icol2*nrow+1))
      		areacol=areacol+.5*base*heightr
      		do irow=1,nrow-1
c	     		lower left triangle
      			base=abs(coord(2,(icol2-2)*nrow+1+irow)
     *			-coord(2,(icol2-2)*nrow+irow))
      			areaelml=base*heightl*.5
c	      		upper left triangle
      			base=abs(coord(2,(icol2-1)*nrow+1+irow)
     *			-coord(2,(icol2-1)*nrow+irow))
      			areaelml=areaelml+.5*base*heightl
c	     		lower right triangle
      			areaelmr=base*heightr*.5
c	      		upper right triangle
      			base=abs(coord(2,icol2*nrow+1+irow)
     *			-coord(2,icol2*nrow+irow))
      			areaelmr=areaelmr+.5*base*heightr
c      			contribution to average density
      			rhoav=rhoav+den((icol2-2)*(nrow-1)+irow)*areaelml
     *			/areacol+den((icol2-1)*(nrow-1)+irow)*areaelmr/areacol
      			areacheck=areacheck+areaelml+areaelmr
      		end do
      	end	if
      	fnode(icol)=slen(icol)*g*rhoav*dy
      end do	

      end



c***********************************************************************
c*                                                                     *
c*    routine to input bulk of data                                    *
c*                                                                     *
c***********************************************************************

      subroutine input(nn,ne,lbw,numvbn,numpbn,nout,ntsts,ncol,nrow,
     *ndf,minter,lda,miter,toler,nrowt,ncolt,nnt,net,rhof,rhoman,ncom,
     *nsing,numsid,vrig,delt,nbn,npass,npoint,convel,
     *epsinv,nout_t,nlrow,upveln,erosl,erosr,peros,rpow,iunflag,
     *iunbeg,ntherm,dtherm,w_depth,beta,prig,rrig,sload,smomen,xadd,
     *ctoler,wtoler,np1,np2,npad,nsing1,ndom,ntbn,lbwt,ldat,nlcol,
     *slpmax,tmax,sealev,dy_flex_init1,dy_flex_init2,wdepth,nsthick,
     *plthick,numvetbn,leqflag,iplasflg,iblay,iblayt,isedl,isedr,
     *ibasflg,intmrkb,nbastrk,nbastary,nbastind,ninbas,ipkfill,ibasfill,
     *sedmax)

      use dyn_arrays
      use dyn_arrays_mech
      use dyn_arrays_therm
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      open(3,file='input/mesh',position='rewind')
      open(2,file='input/connections.dat',position='rewind')

c mech model (# nodes, # elements, l-mesh style, dof and bandwidth)
      read(3,102) nn,ne
      read(3,102) nrow,ncol
      read(3,103)plscale,rlscale,blscale,dfact
      write(6,*)'Mech Mesh'
      write(6,*)' number of nodes=',nn,'number of elems=',ne
      write(6,*) ' nrow=',nrow,'ncol=',ncol
      print*,'Lagrangian mesh'
      print*,' Pro-side stretching factor =',plscale
      print*,' Retro-side strectching factor =',rlscale
      print*,' Base stretchig factor =',blscale
      print*,' Density of mesh compared to eulerian =',dfact
      lbw=2*nrow+3
      ndf=nn*2
      write(6,*)'degrees of freedom=',ndf,'bandwidth=',lbw
c calculate the # rows in abd (stiff. matrix stored in banded form 
c	for use in LAPAK solvers).  Assumes that stiff. is symmetric
      lda=(3*lbw+1)
      print*,'number of rows in abd = ',lda
c input parameters for thermal mesh
      read(3,102)nnt,net
      read(3,102) nrowt,ncolt,nsthick
      write(6,*)'Parameters for Thermal Mesh'
      write(6,*) ' number of nodes = ',nnt,' number of elements = ',net
      write(6,*) ' nrow = ',nrowt,' ncol = ',ncolt      
      write(6,*) ' Sub. slab thickness in nodes =',nsthick
c reference thickness of lith, used in thermal remeshing
      read(3,103)plthick
      write(6,*) ' Sub. slab thickness in meters =',plthick
c spoint node
      read(3,102)nsing
      print*,'Spoint= ',nsing
c convergence velocity and underplating parameters
      read(3,119)convel,upveln
      read(3,107)iunflag,iunbeg
      write(6,*)'Convergence Velocity = ',convel
      print*,'Underplating normal velocity =',upveln
      if(iunflag.eq.0)
     *print*,' Using node location criteria ibeg= ',iunbeg
c rigid visc
      read(3,103)vrig       
      print*,'Rigid Viscosity =',vrig
c compressibility, epsinv
      read(3,103)beta,epsinv,tmax
      print*,'Compressibility =',beta
      print*,'epsinv =',epsinv
      print*,'tmax =',tmax
c flag for linear or non-linear eqns
      read(3,102)leqflag
      if(leqflag.eq.1) print*,'using linear visc. eqns'
      if(leqflag.ne.1) print*,'using non-linear visc. eqns'
c flag for allowing purely plastic deformation (no visc)
      read(3,102)iplasflg
      if(iplasflg.eq.1) then
      	print*,'purely plastic formulation (visc=vrig)'
      else
      	print*,'visco-plasti formulation'
      endif	
c overlying fluid (ocean) and mantle density      
      read(3,103)rhof,rhoman
      print*,'Fluid density	=',rhof
      print*,'Mantle density=',rhoman
c flexural/isostacy parameters
      read(3,102)ncom
      read(3,103)prig,rrig,sload,smomen
      read(3,103)xadd,ctoler,wdepth,wtoler
      read(3,102)np1,np2,npad,nsing1
      if(ncom.eq.0) then
      	print*,'Local isostacy'
      else if(ncom.eq.1) then
      	print*,'One plate flex. Compensation'
      	print*,'D					=',prig
      else if(ncom.eq.2) then
      	print*,'Two plate Flex. Compensation'
      	print*,'  Pro D	=',prig
      	print*,'  Retro D=',rrig
      	print*,'  Sub. load	=',sload
      	print*,'  Sub. moment=',smomen
      	print*,'  Extension of plate 1=',xadd
      	print*,'  Coupling toler.=',ctoler
      	print*,'  Water iter. toler	=',wtoler
      	print*,'  Nodes in plate 1=',np1
      	print*,'  Nodes in plate 2=',np2
      	print*,'  Node padding at edges	=',npad
      	print*,'  Spoint in plate 1 ref	=',nsing1
      endif	
c ref. water depth
      print*,'Water Depth=',wdepth
c Mechanical model boundary conditions
      read(3,102)numvbn,numvetbn,numpbn,numsid,nbn
      print*,'Number of x,y edge velocity boundary nodes =',numvbn
      print*,'Number of tangent edge velocity boundary nodes ='
     *,numvetbn
      print*,'Number of pressure boundary nodes =',numpbn
      print*,'Number of loaded sides =',numsid
      print*,'Number of basal tangent velocity nodes =',nbn
c Time stepping and iteration parameters      
      read(3,104)ntsts,delt
      read(3,102)nout,nout_t
      read(3,109)minter,miter,npass,toler
      print*,'Number of timesteps=',ntsts
      print*,'Timestep length=',delt
      print*,'Output interval for all=',nout
      print*,'Output interval for L-temp=',nout_t
      print*,'Min number of iterations=',minter
      print*,'Max number of iterations =',miter
      print*,'Convergence tolerance	=',toler
      print*,'Number of filtering passes=',npass
c output dt for plasti output for dx
      open(7,file='output/dt_out')
      write(7,103)dble(nout*delt)
      close(7)
      print*,'Plasti output every (my):',dble(nout*delt)
c Erosion parameters
      read(3,103)erosl,erosr,peros,rpow
      print*,'Pro side erosion coef.=',erosl
      print*,'Retro side erosion coef.=',erosr
      print*,'Ridge erosion coef.=',peros
      print*,'Ridge erosion power=',rpow
c sedimentation parameters
      read(3,112)ipkfill,ibasfill,isedl,isedr,sedmax
      if(ipkfill.eq.1) print*,'Sedimentation between peaks'
      if(ibasfill.eq.1) then
      	print*,'Sedimentation in bounding basins, max fill',sedmax
      endif	
      if(ipkfill.eq.1) print*,'Sedimentation bounds =',isedl,isedr
c basin tracking parameters
      read(3,101)ibasflg,intmrkb,nbastary,nbastind
      if(ibasflg.eq.1) then
      	print*,'Basin tracking is on'
      	print*,'	mark every',intmrkb,'time steps'
      else
      	print*,'Basin Tracking is off'
      endif	
      allocate(bastrk(4,nbastary),ibastrk(2,nbastind),ieletpb(nbastary))
      nbastrk=0
      ninbas=0
      bastrk=0.0
      ibastrk=0
      ieletpb=100
c maximum surface slope
      read(3,103)slpmax
      print*,'Maximum Surface Slope =',slpmax
c thermal runup parameters
      read(3,104)ntherm,dtherm
      print*,'Num of therm runup steps =',ntherm
      print*,'Length of runup steps	=',dtherm
c number of boundary layers
      read(3,102)iblay,iblayt
      print*,'Number of boundary layers'
      print*,'	base:',iblay
      print*,'	top:',iblayt
c read in arrays of density, int. angle of friction and cohesion
c	for the mech model
      allocate(den(ne),phi(ne),coh(ne),vmin(ne),q(ne),prex(ne),expn(ne))
      do i=1,ne
      	read(3,113)den(i),phi(i),coh(i),vmin(i),q(i),prex(i),expn(i)
      end do	
c read in nodal coordinates and divide into therm and mech parts
      allocate(coordt(2,nnt),coord(2,nn))
      do i=1,nnt
		read(3,113)coordt(1,i),coordt(2,i)
	  end do	
	  icount=0
	  do i=1,nnt,nrowt
	  	do j=1,nrow
	  		k=i+(j-1)+(nrowt-nrow)
	  		icount=((i-1)/nrowt)*nrow+j
        	coord(1,icount)=coordt(1,k)
        	coord(2,icount)=coordt(2,k)
        end do
      end do  
c read in connections for mech model
      allocate(node(4,ne))
      do i=1,ne
      	read(3,102)node(1,i),node(2,i),node(3,i),node(4,i)
      end do	
c  set up initial thickness vector for isostacy problem
      allocate(zinit(ncol))
      do icol=1,ncol
      	itop=icol*(nrow)
      	ibot=itop-nrow+1
      	zinit(icol)=coord(2,itop)-coord(2,ibot)
      end do
c set initital position of base as equilibrium position for isostacy
      ncount=0
      allocate(zeq(ncol))
      do i=1,nn,nrow
  	  	ncount=(i-1)/nrow +1
      	zeq(ncount)=coord(2,i)
      end do
c define sealevel and make array of water depths
      allocate(wdinit(ncol))
      sealev=wdepth+coord(2,nrow)
      do i=1,ncol
      	wdtemp=sealev-coord(2,i*nrow)
      	if(wdtemp.gt.0.0) then
      		wdinit(i)=wdtemp
      	else
      		wdinit(i)=0.0
      	endif
      end do
      
c output sealevel for DX ploting
      open(8,file='output/sea_level',position='rewind')
      if(wdepth.lt.1.0) then
      	write(8,113)wdepth,coord(1,nrow),coord(1,nn)
      else
		write(8,113)sealev,coord(1,nrow),coord(1,nn)
      end if		
      close(8)
c input velocity boundary nodes
      allocate(nvnd(2,numvbn),bvel(numvbn))
      do i=1,numvbn
      	read(3,114)nvnd(1,i),nvnd(2,i),bvel(i)
      end do
c input pressure boundary nodes
      allocate(npnd(numpbn),bp(numpbn))
      do i=1,numpbn
      	read(3,115)npnd(i),bp(i)
      end do
c  input edge tangential velocity BCs
      allocate(nvtnd(numvetbn),bvelt(numvetbn))
      do i=1,numvetbn
      	read(3,115) nvtnd(i),bvelt(i)
      end do	
c  input base tangential velocity BCs
      allocate(nbase(nbn),basvel(nbn),unvel(nbn))
      do i=1,nbn
      	read(3,115) nbase(i),basvel(i),unvel(i)
      	basvel(i)=basvel(i)
      end do	
c input loaded sides
      allocate(nsnd(4,numsid),bside(numsid))
      if(numsid.eq.0)go to 211
      do i=1,numsid
      	read(3,109)nsnd(1,i),nsnd(2,i),nsnd(3,i),nsnd(4,i),bside(i)
      end do
  211 continue

c####################
c thermal parameters 
c####################

c number of nodes, number of elements, number of domains, num temp BCs
      read(3,102)njunk,njunk1,ndom,ntbn
c node connections and bandwidth
      allocate(nodet(net,5))
      lbwt=0
      read(2,102)njunk
      do i=1,net
      	read(2,102)nodet(i,1),nodet(i,2),nodet(i,3)
      	nodet(i,4)=0
      	nodet(i,5)=0
      	id1=iabs(nodet(i,1)-nodet(i,2))
      	id2=iabs(nodet(i,1)-nodet(i,3))
      	id3=iabs(nodet(i,2)-nodet(i,3))
      	if(id1.gt.lbwt)lbwt=id1
      	if(id2.gt.lbwt)lbwt=id2
      	if(id3.gt.lbwt)lbwt=id3
      end do
      ldat=(3*lbwt+1)
c domain map
      do i=1,net
      	read(3,102)nodet(i,5)
      end do	
c thermal conductivity
      allocate(tcond(2,net))
      do i=1,net
      	read(3,113)tcond(1,i),tcond(2,i)
      end do	
c density in thermal calc
      allocate(trho(net))
      do i=1,net
      	read(3,113)trho(i)
      end do	
c spec. heat
      allocate(spheat(net))
      do i=1,net
      	read(3,113)spheat(i)
      end do	
c heat production
      allocate(hprod(net))
      do i=1,net
      	read(3,113)hprod(i)
      end do	
c input constant temp nodes
      allocate(ntbnd(ntbn),btem(ntbn))
      do i=1,ntbn
      	read(3,115)ntbnd(i),btem(i)
      end do	
c#######################
c flexure profiles
c#######################
      allocate(xp1(np1),yp1(np1),dyinit1(np1),
     *xp2(np2),yp2(np2),dyinit2(np2))
      do i=1,np1
      	read(3,113)xp1(i),yp1(i),dyinit1(i)
      end do
      do i=1,np2
      	read(3,113)xp2(i),yp2(i),dyinit2(i)
      end do	

c#######################
c output file flags
c#######################
      read(3,101)noutput
      allocate(output_flags(noutput))
      do i=1,noutput
      	read(3,101)output_flags(i)
      end do	
c##########################
c Allocate other arrays
c##########################

c	mech. model vel arrays
      allocate(vely(nn),velx(nn))
      do i=1,nn
      	vely(i)=0.0
      	velx(i)=basvel(1)
      end do	
c mean stress and others
      allocate(sbar(ne,3),stress(ne,4),srate(ne,4),sprev(ne))
c gauss point viscosity, plasti failure flag
      allocate(visc(ne,4),ipflag(ne,4),bulkmod(ne,4))
c linear system of eqns
      allocate(rhs(2*nn),abd(lda,2*nn),soln(2*nn),ip(2*nn)) 
c not sure what this is. used in frictg routine which is no longer used
      allocate(vbound(nbn),theta(nbn))
c surface profile
      allocate(xsur(2,ncol),xsurold(2,ncol))
c used in isostacy/flexure calc and remeshing routine
      allocate(ziso(ncol),zinc(2,ncol),cbase(ncol))
c power-law and other viscosities
      allocate(vpower(2,ne))
c mechanical model temps
      allocate(temp(ne),temptc(nn),toldc(nn))
c surface erosion, valley and ridge surfaces
      allocate(veros(2,ncol),vsur(2,ncol),rsur(2,ncol))
      allocate(vdiff(2,ncol),rdiff(2,ncol))
c thermal velocities and temps
      allocate(vx(net),vz(net),tempt(nnt),told(nnt))
      

c record the amount of closed basin filling and max slope cutting
      allocate(basinfill(ncol),nbasinfill(ncol))
      allocate(npeakchop(ncol),peakchop(ncol))
      basinfill=0.0
      nbasinfill=0
      npeakchop=0
      peakchop=0.0
c set initial thickness of the padded regions used in the flexure
c	calculation
      dy_flex_init2=coord(2,nn)-coord(2,nn+1-nrow)
      dy_flex_init1=coord(2,nrow)-coord(2,1)

c make lagrangian mesh
      call mk_lmesh(convel,ntsts,delt,plscale,ncol,rlscale,
     *upveln,blscale,nlrow,nrow,dfact,nn,nlcol,npoint)
      allocate(exhum(npoint))


  101 format(9i5)
  102 format(9i8)
  103 format(4e16.8,i5)
  104 format(i5,2e16.8)
  107 format(i2,i4,e16.8)  
  109 format(3i5,e16.8)
  112 format(4i5,e10.2)
  113 format(9e23.15)
  114 format(2i8,4e23.15)
  115 format(i8,4e23.15)
  119 format(2f8.1,e13.8)

      close(2)
      close(3)
      end  

c*******************************************************************
c*                                                                 *
c*    routine to make lagrangian mesh                              *
c*                                                                 *
c*******************************************************************
      subroutine mk_lmesh(convel,ntsts,delt,plscale,ncol,rlscale,
     *upveln,blscale,nlrow,nrow,dfact,nn,nlcol,npoint)

      use dyn_arrays
      use dyn_arrays_mech
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

c calculate the extent of the mesh past the pro/retro model edges
      xlmin=coord(1,1)-convel*dble(ntsts)*delt*plscale
      xlmax=coord(1,ncol*nrow)*rlscale
c calculate the factor for how far below the model base to extend mesh
c	allows for tracking of underplated particles
      if(upveln.gt.0.0) then
      	blscale=blscale
      else
      	blscale=1.0
      endif	
c number of rows,colms and nodes in  mesh
      nlrow=nrow*dfact
c	desired spacing      
      dxl=(coord(1,nn)-coord(1,1))/(dble(ncol)*dfact)
c	number of colms
      nlcol=floor((xlmax-xlmin)/dxl)+1
c 	number of nodes      
      npoint=nlcol*nlrow      
      allocate(tpoint(7,npoint),ieletp(npoint))

c mesh domain starting from retro side
c	find begining eulerian node 
      do i=1,ncol
      	if(coord(1,nrow*i).ge.xlmax) then
      		ibeg=i
      		exit
      	end if
      end do	

c 	mesh domain
c	first row
      index=0
	  ylbase=coord(2,(ibeg-1)*nrow+1)/blscale
	  xl=coord(1,ibeg*nrow)
	  dy=(coord(2,ibeg*nrow)-ylbase)/dble(nlrow-1)
	  do i=1,nlrow
	  	index=index+1
	  	tpoint(1,index)=xl
	  	tpoint(2,index)=ylbase+dy*dble(i-1)
	  	ieletp(index)=100
	  	do k=3,7
	  		tpoint(k,index)=0.0
	  	end do
	  end do	
c	all other nodes
      iflag=0
c	e nodes on the l and r if l node
      iel=ibeg-1
      ier=ibeg
      do j=2,nlcol
      	xl=xl-dxl
c	if l node is outside bounds of e nodes (iel and ier) find new e node bounds
      	if(xl.lt.coord(1,iel*nrow)) then
c			if at pro edge, use pro edge thickness for the rest of the mesh
      		if(iel.eq.1) then
      			iel=iel
      			ier=ier
      			iflag=1
      		else
      			do i=iel,1,-1
      				if(coord(1,i*nrow).lt.xl) then
      					iel=i
      					ier=i+1
      					exit
      				endif
      			end do	
      		end if
      	end if
c		if still in mech domain      	
      	if(iflag.eq.0) then
c 			slope at model surface and base
      		tslope=(coord(2,iel*nrow)-coord(2,ier*nrow))/
     *		(coord(1,iel*nrow)-coord(1,ier*nrow))
      		bslope=(coord(2,(iel-1)*nrow+1)-coord(2,(ier-1)*nrow+1))/
     * 		(coord(1,iel*nrow)-coord(1,ier*nrow))
c     		x dist of lnode past rhs e node
      		difx=xl-coord(1,ier*nrow)
c     		surface and base of lmesh at l node
      		ylbase=(coord(2,(ier-1)*nrow+1)+difx*bslope)/blscale
      		yltop=coord(2,ier*nrow)+difx*tslope
c			vert. spacing of l nodes
      		dy=(yltop-ylbase)/dble(nlrow-1)
      		do i=1,nlrow
      			index=index+1
      			tpoint(1,index)=xl
      			tpoint(2,index)=ylbase+dble(i-1)*dy
      			ieletp(index)=100
      			do k=3,7
      				tpoint(k,index)=0.0
      			end do
      		end do
c		if outside mech domain      		
      	else
      		yltop=coord(2,nrow)
      		ylbase=coord(2,1)/blscale
      		dy=(yltop-ylbase)/dble(nlrow-1)
      		do i=1,nlrow
      			index=index+1
      			tpoint(1,index)=xl
      			tpoint(2,index)=ylbase+dble(i-1)*dy
      			ieletp(index)=100
      			do k=3,7
      				tpoint(k,index)=0.0
      			end do
      		end do	
      	endif
      end do	
      end
      				

c**********************************************************************
c*
c  routine to filter checkerboard out of pressure field
c
c***********************************************************************
      subroutine pfilt(ne,nr,npass)
      use dyn_arrays_mech
      implicit real*8 (a-h,o-z)
      nrow=nr-1
      ncol=ne/(nrow)
      do 100 ipass=1,npass
c
c  filter corners
c  
      sprev(1)=.5*(sbar(1,1)+sbar(1+nrow,1))
      sprev(nrow)=.5*(sbar(nrow,1)+sbar(nrow+nrow,1))
      sprev(ne)=.5*(sbar(ne,1)+sbar(ne-nrow,1))
      sprev(ne-nrow+1)=.5*(sbar(ne-nrow+1,1)+sbar(ne-nrow+1-nrow,1))
c
c  filter edges
c  
      icol=1
c!OMP parallel
c!OMP do private(irow,iele)
      do 10 irow=2,nrow-1
      iele=(icol-1)*(nrow)+irow
      sprev(iele)=sbar(iele,1)/2.0d0+sbar(iele-1,1)/6.0d0
     *+sbar(iele+1,1)
     */6.0d0+sbar(iele+nrow,1)/6.0d0
   10 continue
c!OMP end do
      icol=ncol
c!OMP do private(irow,iele)
      do 20 irow=2,nrow-1
      iele=(icol-1)*(nrow)+irow
      sprev(iele)=sbar(iele,1)/2.0d0+sbar(iele-1,1)/6.0d0
     *+sbar(iele+1,1)
     */6.0d0+sbar(iele-nrow,1)/6.0d0
   20 continue
c!OMP end do
      irow=1
c!OMP do private(icol,iele)
      do 30 icol=2,ncol-1
      iele=(icol-1)*(nrow)+irow
      sprev(iele)=sbar(iele,1)/2.0d0
     *+sbar(iele+nrow,1)/4.0d0
     *+sbar(iele-nrow,1)/4.0d0
   30 continue
c!OMP end do
      irow=nrow
c!OMP do private(icol,iele)
      do 40 icol=2,ncol-1
      iele=(icol-1)*(nrow)+irow
      sprev(iele)=sbar(iele,1)/2.0d0
     *+sbar(iele+nrow,1)/4.0d0
     *+sbar(iele-nrow,1)/4.0d0
   40 continue
c!OMP end do
c
c  filter interior
c
c!OMP do private(icol,irow,iele)
      do 49 icol=2,ncol-1
      do 50 irow=2,nrow-1
      iele=(icol-1)*(nrow)+irow
      sprev(iele)=sbar(iele,1)/2.0d0+sbar(iele-1,1)/8.0d0
     *+sbar(iele+1,1)/8.0d0
     *+sbar(iele+nrow,1)/8.0d0
     *+sbar(iele-nrow,1)/8.0d0
   50 continue
   49 continue
c!OMP end do
c
c  put back into sbar
c
c!OMP do private(iele)
      do 80 iele=1,ne
      sbar(iele,1)=sprev(iele)
   80 continue
c!OMP end do
c!OMP end parallel
  100 continue
      return
      end
c***********************************************************************
c*                                                                     *
c* routine to calculate stresses, strain rates and viscosity           *
c*                                                                     *
c***********************************************************************
      subroutine ss(ne,nn,lbw,delt,ndf,nrow,c,beta,vrig)

	  use dyn_arrays_mech
	  use dyn_arrays
	  parameter(nstbis=21)
      implicit real*8 (a-h,o-z)
      double precision delt,c,beta,vrig
      dimension ix(9),lr(9),lz(9),lw(9)
      real*8 ul(2,4),d(10),xl(2,4),bulkl(1,9),viscl(1,9),p(2,1)
     *,sigavl(1,4),epsavl(1,4),kel(9,9),sbarl(1,3),pflagl(1,9),
     *deltl,xs(2,2),sx(2,2),erhs(3),cmpp1(6),shps(4),shpt(4)
	  logical*1 flg
      common /eltvar/ iele,d,irow,ul,xl,i,j,kel,deltl,
     c viscl,bulkl,pflagl,sbarl,p,sigavl,epsavl,epsinv,vpow,l,
     c ptot,volt,ptemp,lint,shpp,k,v,ddv,dv,epstra,epsdev,xvol,
     c xlam,xcom,xmu,xrho,devstre,rj2d,presl,ivmises,
     c sigy,cosphi,sinphi,radret,rj2de,lp,mp,k1,a1,a2,a3,
     c cmpp1,erhs,etemp,lloc,iax,ino,piv,ii,nn2,ipjp,dd,ijp,jip,
     c ij,fac,nstu,flg,inopredv,isw,
     c sg,tg,wg,shp,fp,dperm,pgg,ig,sx,xs,xsj,tp,
     c a4,a5,a6,b1,b2,j1,sum,cmpp,iia,cdpu,ix,stressl,d5,xdiv

      data lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
      data lw/4*25,4*40,64/
      data shps/-0.5,0.5,0.5,-0.5/,shpt/-0.5,-0.5,0.5,0.5/

c add declarations for elt03n vars
	dimension shp(3,9),sg(9),tg(9),wg(9),sig(6),eps(3),wd(2),
     *  v(2),dv(2,2),shpp(3),indx(nstbis),cmpp(6),
     *  cdpu(54),fp(3),ptot(3),devstre(4),epsdev(4),stressl(4)
      ndfe=2
      ndm=2
      nst=9
      nen=4
      nel=4
      kstep=2
      n=1
      maxn=1
c
c loop over each element
c
      do iele=1,ne
      	press=sbar(iele,1)
      	if(press.lt.0.0)press=0.00
      end do

	do 100 iele=1,ne
      inopredv=0
      isw=3
      d(1)=vpower(1,iele)
      d(2)=(1.0/beta)
      d(3)=(den(iele))
      d(4)=2
      irow=mod(iele,nrow-1)
      if(irow.eq.0)irow=nrow-1
      d(5)=phi(iele)
      d(6)=coh(iele)
      d(7)=1.0e14
      d(8)=0.0
      d(9)=0.0
      d(10)=0.0
      do 35 j=1,4
      ul(1,j)=(velx(node(j,iele)))
      ul(2,j)=(vely(node(j,iele)))
   35 continue
      do 45 j=1,4
      xl(1,j)=(coord(1,node(j,iele)))
      xl(2,j)=(coord(2,node(j,iele)))
   45 continue
      do 291 i=1,9
      do 290 j=1,9
      kel(i,j)=0.0
  290 continue
  291 continue
      deltl=(delt)
c     ix(1)=0.0
c     p(1,1)=0.0
      do 150 j=1,4
      viscl(1,j)=(visc(iele,j))
      bulkl(1,j)=(bulkmod(iele,j))
      pflagl(1,j)=dble(ipflag(iele,j))
  150 continue
      sbarl(1,1)=(sbar(iele,1))
c      call elt03n(inopredv,d,ul,xl,ix,kel,p,ndfe,ndm,nst,isw,deltl
c     *,nen,n,nel,viscl,bulkl,sbarl,pflagl,sigavl,epsavl,maxn,kstep)
*********  add in subroutine
	nelp=1
      l=d(4)
      nstu=ndfe*nen
c so dp is in pinc
      ptot(1)=sbarl(n,1)
      ptot(2)=sbarl(n,2)
      ptot(3)=sbarl(n,3)
c
c replace call with subroutine -- OK
c      call pgauss(l,lint,sg,tg,wg)
	pgg=1./dsqrt(3.0d0)
	lint=l*l
	do ig=1,4
	  sg(ig)=pgg*lr(ig)
	  tg(ig)=pgg*lz(ig)
	  wg(ig)=1.
	end do
c end of pgauss
      volt=0.
      ptemp=0.
      sigavl(n,1)=0.
      sigavl(n,2)=0.
      sigavl(n,3)=0.
      sigavl(n,4)=0.
      epsavl(n,1)=0.
      epsavl(n,2)=0.
      epsavl(n,3)=0
      epsavl(n,4)=0.
      do 65 l=1,lint
c replace with subroutine lines -- creates error in elimp
c      call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
	flg=.false.
      do 103 i=1,4
      shp(3,i)=(0.5+shps(i)*sg(l))*(0.5+shpt(i)*tg(l))
      shp(1,i)=shps(i)*(0.5+shpt(i)*tg(l))
      shp(2,i)=shpt(i)*(0.5+shps(i)*sg(l))
  103 continue
      if(nel.ge.4)goto 120
      do 110 i=1,3
      shp(i,3)=shp(i,3)+shp(i,4)
  110 continue
  120 if(nel.gt.4)call shap2(sg(l),tg(l),shp,ix,nel)
      do 132 i=1,ndm
      do 131 j=1,2
      xs(i,j)=0.0
      do 130 k=1,nel
      xs(i,j)=xs(i,j)+xl(i,k)*shp(j,k)
  130 continue
  131 continue
  132 continue
      xsj=xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if(flg) goto 141
      sx(1,1)=xs(2,2)/xsj
      sx(2,2)=xs(1,1)/xsj
      sx(1,2)=-xs(1,2)/xsj
      sx(2,1)=-xs(2,1)/xsj
      do 140 i=1,nel
      tp=shp(1,i)*sx(1,1)+shp(2,i)*sx(2,1)
      shp(2,i)=shp(1,i)*sx(1,2)+shp(2,i)*sx(2,2)
      shp(1,i)=tp
  140 continue
  141 continue
c end of shape
      shpp(1)=1.
      shpp(2)=sg(l)
      shpp(3)=tg(l)
c compute v at l
      do 38 i=1,2
      v(i)=0.
      do 31 k=1,nel
      v(i)=v(i)+shp(3,k)*ul(i,k)
   31 continue
c compute gradv at l
      do 37 j=1,2
      ddv=0.0
c FIX this loop -- don't bother, nel is small
      do 32 k=1,nel
      ddv=ddv+shp(j,k)*ul(i,k)
   32 continue
      dv(i,j)=ddv
   37 continue
   38 continue
c from dv every strain or spin rate ...
      epstra=(dv(1,1)+dv(2,2))/3.
c convention 1=xx 2=yy 3=xy(not 2*xy) 4=zz=out of plane
      epsdev(1)=dv(1,1)-epstra
      epsdev(2)=dv(2,2)-epstra
      epsdev(3)=(dv(1,2)+dv(2,1))/2.
c because this is the plane strain elmt
      epsdev(4)=0.
      epsavl(n,1)=epsavl(n,1)+dv(1,1)
      epsavl(n,2)=epsavl(n,2)+dv(2,2)
      epsavl(n,3)=epsavl(n,3)+0.5*(dv(1,2)+dv(2,1))
      epsavl(n,4)=epsavl(n,3)+0.5*(dv(1,2)-dv(2,1))
c     linear case   or no predictor
      if(inopredv.eq.1)then
c or restart
      if(kstep.eq.1)then
      viscl(n,l)=d(1)
      bulkl(n,l)=d(2)
                    endif
      xvol=1.
      xlam=deltl*bulkl(n,l)
      xcom=1.0/xlam
      xmu=viscl(n,l)
      xrho=d(3)
                      else
c     nonlinear case = nonlinear iteration technique .
c  in the general case use this to compute stress predictor and next stif
      xvol=1.
      xlam=deltl*bulkl(n,l)
      xcom=1.0/xlam
      xrho=d(3)
      xmu=viscl(n,l)
c strain stress law here viscous '!isotropic!'
      devstre(1)=2.*xmu*epsdev(1)
      devstre(2)=2.*xmu*epsdev(2)
      devstre(3)=2.*xmu*epsdev(3)
      devstre(4)=2.*xmu*epsdev(4)
c invariant
      rj2d=(devstre(1)*devstre(1)+devstre(2)*devstre(2))/2.0+
     *devstre(3)*devstre(3)
      rj2d=sqrt(rj2d)
c refind pressure at gauss point level
      presl=0.0
      do 1155 i=1,nelp
c dont forget to update press(n,i)
      presl=presl+(ptot(i))*shpp(i)
 1155 continue
      stressl(1)=-presl+devstre(1)
      stressl(2)=-presl+devstre(2)
      stressl(3)=devstre(3)
      stressl(4)=-presl+devstre(4)
      ptemp=ptemp+presl
      sigavl(n,1)=sigavl(n,1)+stressl(1)
      sigavl(n,2)=sigavl(n,2)+stressl(2)
      sigavl(n,3)=sigavl(n,3)+stressl(3)
      sigavl(n,4)=sigavl(n,4)+stressl(4)
c compute state variable control
      ivmises=0
      if(d(5).lt.0.)then
      ivmises=1
      sigy=-d(5)
                    else
      d5=3.14159*d(5)/180.
      cosphi=dcos(d5)
      sinphi=dsin(d5)
      coh2=d(6)
      if(presl.lt.0.0)then
      	sigy=coh2*cosphi
      else
      	sigy=presl*sinphi+coh2*cosphi
       if(sigy.lt.0.0)print*,d(5)
c      	if(sigy.lt.0.0)print*,presl,sinphi
c      	if(sigy.lt.0.0)print*,coh2,cosphi
      endif
                    endif
      if(sigy.gt.d(7))sigy=d(7)
      if(sigy.lt.0.0)write(*,*)'pos 2: sigy < 0 elt n= ',n
c
c     radial return
c
      radret=rj2d/sigy
         if(pflagl(n,l).gt.0.or.radret.gt.1.)then
c plastic flow
         pflagl(n,l)=1
c notice that the following computation is redundant.
         rj2de=(epsdev(1)**2+epsdev(2)**2)/2.0+epsdev(3)**2
         rj2de=dsqrt(rj2de)
         xmu=sigy/(2.0*rj2de)
         if(xmu.gt.d(1))then
                        xmu=d(1)
                        pflagl(n,l)=0
                             endif
                                             endif
c     update nonlinear rheology
c here in general the whole rheology is reparametrized
             viscl(n,l)=xmu
                      endif
c

      xvol=xvol*xsj*wg(l)
      xlam=xlam*xsj*wg(l)
      xcom=xcom*xsj*wg(l)
      xmu=xmu*xsj*wg(l)
      xrho=xrho*xsj*wg(l)
      volt=volt+xvol
c     write(2,*)'end control                        '
c
c
c     isotropic operator     : spp
c     (dev-is  coupling)
c
      do 400 lp=1,nelp
      do 401 mp=1,nelp
      kel(nstu+lp,nstu+mp)=kel(nstu+lp,nstu+mp)+
     1xcom*shpp(lp)*shpp(mp)
  401 continue
  400 continue
c     write(2,*)'end spp                            '
c
c
c
c
      if(isw.eq.6)goto 60
      k1=1
c nel = 4, so not worth parallelizing?
      do 34 k=1,nel
c add this line
c	k1=1+(k-1)*ndfe
      a1=xmu*shp(1,k)
      a2=xmu*shp(2,k)
      a3=xrho*(dv(1,1)*shp(3,k)+v(1)*shp(1,k)+v(2)*shp(2,k))
      a4=xrho*(dv(2,2)*shp(3,k)+v(1)*shp(1,k)+v(2)*shp(2,k))
      a5=xrho*dv(1,2)*shp(3,k)
      a6=xrho*dv(2,1)*shp(3,k)
c eliminate deviatoric part
c     b1=xlam*shp(1,k)
c     b2=xlam*shp(2,k)
      b1=0.
      b2=0.
      j1=1
      do 33 j=1,nel
c add this line
c	j1=1+(j-1)*ndfe
c
c
c     deviatoric operator    : suu
c     (dev-dev coupling)
c
c xj xk
      kel(j1,k1)=kel(j1,k1)+shp(1,j)*a1+shp(2,j)*a2
      kel(j1,k1)=kel(j1,k1)+(shp(1,j)*a1)/3.0
c xj yk
      kel(j1,k1+1)=kel(j1,k1+1)+0.
c     kel(j1,k1+1)=kel(j1,k1+1)+a1*shp(2,j)/3.0
      kel(j1,k1+1)=kel(j1,k1+1)-2.*a2*shp(1,j)/3.0+a1*shp(2,j)
c yj xk
      kel(j1+1,k1)=kel(j1+1,k1)+0.
c     kel(j1+1,k1)=kel(j1+1,k1)+a2*shp(1,j)/3.0
      kel(j1+1,k1)=kel(j1+1,k1)-2.*a1*shp(2,j)/3.0+a2*shp(1,j)
c yj yk
      kel(j1+1,k1+1)=kel(j1+1,k1+1)+shp(1,j)*a1+shp(2,j)*a2
      kel(j1+1,k1+1)=kel(j1+1,k1+1)+(shp(2,j)*a2)/3.0
c this if statement breaks the elegance of the code helas!
c     write(2,*)'end suu                            '
      if(k.eq.1)then
c
c
c     iso-dev   operator     : sup
c     (dev-is  coupling)
c
      do 333 mp=1,nelp
      kel(nstu+mp,j1)=kel(nstu+mp,j1)+xvol*shpp(mp)*shp(1,j)
      kel(nstu+mp,j1+1)=kel(nstu+mp,j1+1)+xvol*shpp(mp)*shp(2,j)
      kel(j1,nstu+mp)=kel(nstu+mp,j1)
      kel(j1+1,nstu+mp)=kel(nstu+mp,j1+1)
  333 continue
c	write(6,*) 'kel:',kel(1,1),kel(9,9)
      endif
      j1=j1+ndfe
   33 continue
      k1=k1+ndfe
   34 continue
c     write(2,*)'end sup                            '
c
c
c
c     solve iso-dev coupling at the element level :
c     elimination of internal dofs .here pressure.
c
c     if u-u convective term is not zero s is not symmetric
c     if u-p convective term is not zero s is not symmetric
c        u-p convective term arises from stress rate computations
      goto 65
c  force-computation
   60 continue
c it is very important to notice here that 3 modes can exist:
c a neutral mode: 1.no reaction computed  ,return.(ex a simple fluid)
c                 2.reaction computed but not fed to rhs
c                 3.standard mode compute and feed reactions.
c we lock here in mode 1 but macro could call other modes as well
      xdiv=(dv(1,1)+dv(2,2))*xlam
      do 67 k=1,nel
      do 64 j=1,2
      sum=xdiv*shp(j,k)
      do 63 i=1,2
      sum=sum+xmu*(dv(j,i)+dv(i,j))*shp(i,k)+
     1   xrho*v(i)*dv(j,i)*shp(3,k)
   63 continue
      p(j,k)=p(j,k)-sum
   64 continue
   67 continue
c     write(2,*)'end loop 65 l                      '
   65 continue
c     write(2,*)'loop 65 terminated'
c save cmpp and cdpu in file 3
      if(isw.eq.3.and.inopredv.eq.0)then
      if(nelp.eq.1)then
      cmpp(1)=kel(nstu+1,nstu+1)
c	write(6,*) kel(nstu+1,nstu+1)
      iia=0
      do 1111 j=1,nstu
c commented this line
c         iia=iia+1
c added this line
	 iia=j
         cdpu(j)=kel(nstu+1,j)
c	 write(6,*) cdpu(iia),iia,nstu,j,kel(nstu+1,j)
 1111 continue
                   endif
cc    nelps=nelp*nelp
cc    write(mswap1,*)(cmpp(i),i=1,nelps)
cc    nelpu=nelp*nstu
cc    write(mswap1,*)(cdpu(i),i=1,nelpu)
c solve cmpp*p+cdpu*u=cmpp*p0
cc dont  forget to update fp in compressible case
c	write(6,*) 'calling elimp',ndf,nelp,nstu
c
c replace with subroutine
c      call elimp(ndfe,fp,cmpp,cdpu,ptot,nelp,nstu,ul)
      lloc=0
      iax=0
      ino=1
      do 10 i=1,nelp
      etemp=-0.
c      etemp=-fp(i)
c	write(6,*) 'fp(1)',fp(1), etemp, nelp, nstu
      do 20 j=1,nstu
      iax=iax+1
      lloc=lloc + 1
      etemp=etemp - cdpu(lloc)*ul(iax,ino)
c	write(6,*) etemp,cdpu(lloc),lloc,ul(iax,ino),iax,ino
      if(iax.eq.2)then
                  ino=ino+1
                  iax=0
                  endif
   20 continue
      erhs(i)=etemp
   10 continue
      if (nelp.ne.1) go to 600
      piv=cmpp(1)
      if (piv.eq.0.00) then
c	 write(6,*) piv,cmpp(1),etemp,erhs(1)
          print*,'error  kpp is not invertible - stop in elimp'
          stop
      endif
      ptot(1)=erhs(1)/piv
	goto 630
c
  600 continue
c	write(6,*) 'nelp ne 1',nelp
c
c     move akpp to the working array akpp1
c
c      ii=0
c      do 114 i=1,nelp
c      do 105 j=i,nelp
c      ii=ii+1
c      cmpp1(ii)=cmpp(ii)
c  105 continue
c  114 continue
cc
c      nn2=nelp + 2
c      ipjp=-nelp
c      do 116 ip=1,nelp-1
c      ipjp=ipjp + nn2 - ip
c      piv=cmpp1(ipjp)
c      if (piv.eq.0.00) then
c     	print*,'error: kpp is not invertible - stop in elimp'
c          stop
c      endif
c      dd=1.0/piv
cc
c      ii=ipjp
c      ijp=ipjp
c      do 151 i=ip+1,nelp
c      ii=ii + nn2 - i
c      ijp=ijp + 1
c      fac=cmpp1(ijp)*dd
cc
c      jip=ijp - 1
c      ij=ii - 1
c      do 121 j=i,nelp
c      jip=jip + 1
c      ij=ij + 1
c      cmpp1(ij)=cmpp1(ij) - fac*cmpp1(jip)
c  121 continue
c      erhs(i)=erhs(i) - fac*erhs(ip)
c  151 continue
c  116 continue
cc
c      piv=cmpp1(ij)
c      if (piv.eq.0.00) then
c      	print*,'error: kpp is not invertible - stop in elimp'
c          stop
c      endif
c      ptot(nelp)=erhs(nelp)/piv
c      ii=ij
c      do 133 i=nelp-1,1,-1
c      temp=erhs(i)
c      ii=ii - nn2 + i + 1
c      ij=ii
c      do 142 j=i+1,nelp
c      ij=ij + 1
c      temp=temp - cmpp1(ij)*ptot(j)
c  142 continue
c      ptot(i)=temp/cmpp1(ii)
c  133 continue
  630 continue
c
c
c
c      end of subroutine

      sbarl(n,1)=ptot(1)
      sbarl(n,2)=ptot(2)
      sbarl(n,3)=ptot(3)
                   endif
      if(nelp.eq.1)then
c     sigav(n,1)=(sigav(n,1)+ptemp-ptot(1))/lint
c     sigav(n,2)=(sigav(n,2)+ptemp-ptot(1))/lint
c     sigav(n,4)=(sigav(n,4)+ptemp-ptot(1))/lint
c     sigav(n,3)=sigav(n,3)/lint
      sigavl(n,1)=(sigavl(n,1)+ptemp)/lint
      sigavl(n,2)=(sigavl(n,2)+ptemp)/lint
      sigavl(n,4)=(sigavl(n,4)+ptemp)/lint
      sigavl(n,3)=sigavl(n,3)/lint
                   endif
      epsavl(n,1)=epsavl(n,1)/lint
      epsavl(n,2)=epsavl(n,2)/lint
      epsavl(n,3)=epsavl(n,3)/lint
      epsavl(n,4)=epsavl(n,4)/lint
c
c try hand elimination for n/1 elt
      if(nelp.eq.1)then
c     compress=volt/d(2)
c compressibility included in spp for possible nonlinear iterations
      do 5555 i=1,nstu
      do 5556 j=1,nstu
c     s(i,j)=s(i,j)+s(i,9)*s(9,j)/compress
      kel(i,j)=kel(i,j)+kel(i,9)*kel(9,j)/kel(9,9)
 5556 continue
 5555 continue
                   endif
      if(nelp.eq.3)then
c not invoked in this case
	  	write(6,*) 'nelp = 3'
      	call ludcmp(kel,nst,nst,indx(1),dperm,nelp)
      endif
c**** end of added subr elt03n

      do 300 j=1,4
      stress(iele,j)=sigavl(1,j)
      srate(iele,j)=epsavl(1,j)
      visc(iele,j)=viscl(1,j)
      ipflag(iele,j)=int(pflagl(1,j))
  300 continue
      sbar(iele,1)=sbarl(1,1)
      epsinv=(srate(iele,1)**2+srate(iele,2)**2)/2.+srate(iele,3)*
     *srate(iele,3)
      if(epsinv.lt.0.0)epsinv=0.0
      epsinv=dsqrt(epsinv)
      vpow=vpower(2,iele)*(epsinv**(1./expn(iele)-1.))
      if(vpow.gt.vrig)vpow=vrig
      if(vpow.lt.vmin(iele))vpow=vmin(iele)
      vpower(1,iele)=vpow
      do 99 j=1,4
      if(ipflag(iele,j).eq.0)visc(iele,j)=vpower(1,iele)
   99 continue
  100 continue
c!OMP end parallel do
      return
      end
c***********************************************************************
c*                                                                     *
c* routine to assemble global stiffness matrix and rhs                 *
c*                                                                     *
c***********************************************************************
      subroutine globe(ne,nn,lbw,delt,lda,ndf,nrow,ldf,c,beta,
     *vrig,sigav,epsav,itst)

      use dyn_arrays_mech
      use dyn_arrays
      implicit real*8 (a-h,o-z)
	  parameter(nstbis=21)
      double precision delt,c,beta,vrig
      logical flg
      real*8 ul(2,4),d(10),xl(2,4),bulkl(1,9),viscl(1,9),p(2,1)
     *,sigavl(1,4),epsavl(1,4),kel(9,9),sbarl(1,3),pflagl(1,9),
     *amass(9,9),bel(8),sigav(ne,4),epsav(ne,4),deltl,xs(2,2),
     *sx(2,2),erhs(3),cmpp1(6)
      character date*10, time3*10
	  dimension shp(3,9),sg(9),tg(9),wg(9),sig(6),eps(3),wd(2),
     *v(2),dv(2,2),shpp(3),indx(nstbis),cmpp(6),cdpu(54),fp(3),
     *ptot(3),devstre(4),epsdev(4),stressl(4),ldf(*),ix(9),lr(9),
     *lz(9),lw(9),shps(4),shpt(4)

      data lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
      data lw/4*25,4*40,64/
      data shps/-0.5,0.5,0.5,-0.5/,shpt/-0.5,-0.5,0.5,0.5/

      ndfe=2
      ndm=2
      nst=9
      nen=4
      nel=4
      kstep=2
      n=1
      maxn=1
      g=-9.8d0
      flg=.false.
c
c
c   initialize stiffness matrix and rhs
c
      ibd=2*(nrow-1)
      mbw=2*lbw+1
      do 11 j=1,ndf
      rhs(j)=0.d0
      do 10 i=1,lda
      abd(i,j)=0.d0
   10 continue
   11 continue
c
c loop over each element
c
      do iele=1,ne
      	press=sbar(iele,1)
      	if(press.lt.0.0)press=0.0d0
      end do

c
c  calc element stiffness matrix
c

      do 100 iele=1,ne
      inopredv=1
      isw=3
      d(1)=vpower(1,iele)
      d(2)=(1.0/beta)
      d(3)=(den(iele))
      d(4)=2.0
      d(5)=0.0
      d(6)=0.0
      d(7)=0.0
      d(8)=0.0
      d(9)=0.0
      d(10)=0.0
      do 35 j=1,4
      ul(1,j)=(velx(node(j,iele)))
      ul(2,j)=(vely(node(j,iele)))
   35 continue
      do 45 j=1,4
      xl(1,j)=(coord(1,node(j,iele)))
      xl(2,j)=(coord(2,node(j,iele)))
   45 continue
      do 292 i=1,9
      do 290 j=1,9
      kel(i,j)=0.0
  290 continue
  292 continue
      deltl=(delt)
c     ix(1)=0.0
c     p(1,1)=0.0
      do 150 j=1,4
      viscl(1,j)=(visc(iele,j))
      bulkl(1,j)=(bulkmod(iele,j))
      pflagl(1,j)=dble(ipflag(iele,j))
  150 continue
      sbarl(1,1)=(sbar(iele,1))
c      call elt03n(inopredv,d,ul,xl,ix,kel,p,ndfe,ndm,nst,isw,deltl
c     *,nen,n,nel,viscl,bulkl,sbarl,pflagl,sigavl,epsavl,maxn,kstep)
c replace with subroutine lines
	  nelp=1
      l=d(4)
      nstu=ndfe*nen
c get pressure back from saved matrices and velocities
c
c      call pgauss(l,lint,sg,tg,wg)
c replace with subr
	pgg=1./dsqrt(3.0d0)
	lint=l*l
	do ig=1,4
	  sg(ig)=pgg*lr(ig)
	  tg(ig)=pgg*lz(ig)
	  wg(ig)=1.
	end do
c end of pgauss
      volt=0.
      ptemp=0.
      sigavl(n,1)=0.
      sigavl(n,2)=0.
      sigavl(n,3)=0.
      sigavl(n,4)=0.
      epsavl(n,1)=0.
      epsavl(n,2)=0.
      epsavl(n,3)=0
      epsavl(n,4)=0.
      do 65 l=1,lint
c replace with subr
c      call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
      do 103 i=1,4
      shp(3,i)=(0.5+shps(i)*sg(l))*(0.5+shpt(i)*tg(l))
      shp(1,i)=shps(i)*(0.5+shpt(i)*tg(l))
      shp(2,i)=shpt(i)*(0.5+shps(i)*sg(l))
  103 continue
      if(nel.ge.4)goto 120
      do 110 i=1,3
      shp(i,3)=shp(i,3)+shp(i,4)
  110 continue
  120 if(nel.gt.4)call shap2(sg(l),tg(l),shp,ix,nel)
      do 132 i=1,ndm
      do 131 j=1,2
      xs(i,j)=0.0
      do 130 k=1,nel
      xs(i,j)=xs(i,j)+xl(i,k)*shp(j,k)
  130 continue
  131 continue
  132 continue
      xsj=xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if(flg) goto 141
      sx(1,1)=xs(2,2)/xsj
      sx(2,2)=xs(1,1)/xsj
      sx(1,2)=-xs(1,2)/xsj
      sx(2,1)=-xs(2,1)/xsj
      do 140 i=1,nel
      tp=shp(1,i)*sx(1,1)+shp(2,i)*sx(2,1)
      shp(2,i)=shp(1,i)*sx(1,2)+shp(2,i)*sx(2,2)
      shp(1,i)=tp
  140 continue
  141 continue
cc end of shape
      shpp(1)=1.
      shpp(2)=sg(l)
      shpp(3)=tg(l)
c compute v at l
      do 38 i=1,2
      v(i)=0.
      do 31 k=1,nel
      v(i)=v(i)+shp(3,k)*ul(i,k)
   31 continue
c compute gradv at l
      do 37 j=1,2
      ddv=0.0
c FIX this loop -- don't bother, nel is small
      do 32 k=1,nel
      ddv=ddv+shp(j,k)*ul(i,k)
   32 continue
      dv(i,j)=ddv
   37 continue
   38 continue
c from dv every strain or spin rate ...
      epstra=(dv(1,1)+dv(2,2))/3.
c convention 1=xx 2=yy 3=xy(not 2*xy) 4=zz=out of plane
      epsdev(1)=dv(1,1)-epstra
      epsdev(2)=dv(2,2)-epstra
      epsdev(3)=(dv(1,2)+dv(2,1))/2.
c because this is the plane strain elmt
      epsdev(4)=0.
      epsavl(n,1)=epsavl(n,1)+dv(1,1)
      epsavl(n,2)=epsavl(n,2)+dv(2,2)
      epsavl(n,3)=epsavl(n,3)+0.5*(dv(1,2)+dv(2,1))
      epsavl(n,4)=epsavl(n,3)+0.5*(dv(1,2)-dv(2,1))
c     linear case   or no predictor
      if(inopredv.eq.1)then
c 		or restart
      	if(kstep.eq.1)then
      		viscl(n,l)=d(1)
      		bulkl(n,l)=d(2)
        endif
      	xvol=1.
      	xlam=deltl*bulkl(n,l)
      	xcom=1.0/xlam
      	xmu=viscl(n,l)
      	xrho=d(3)
      endif
c

      xvol=xvol*xsj*wg(l)
      xlam=xlam*xsj*wg(l)
      xcom=xcom*xsj*wg(l)
      xmu=xmu*xsj*wg(l)
      xrho=xrho*xsj*wg(l)
      volt=volt+xvol
c     write(2,*)'end control                        '
c
c
c     isotropic operator     : spp
c     (dev-is  coupling)
c
      do 400 lp=1,nelp
      do 401 mp=1,nelp
      kel(nstu+lp,nstu+mp)=kel(nstu+lp,nstu+mp)+
     1xcom*shpp(lp)*shpp(mp)
  401 continue
  400 continue
c     write(2,*)'end spp                            '
c
      k1=1
c nel = 4, so not worth parallelizing?
      do 34 k=1,nel
c add this line
c	k1=1+(k-1)*ndfe
      a1=xmu*shp(1,k)
      a2=xmu*shp(2,k)
      a3=xrho*(dv(1,1)*shp(3,k)+v(1)*shp(1,k)+v(2)*shp(2,k))
      a4=xrho*(dv(2,2)*shp(3,k)+v(1)*shp(1,k)+v(2)*shp(2,k))
      a5=xrho*dv(1,2)*shp(3,k)
      a6=xrho*dv(2,1)*shp(3,k)
c eliminate deviatoric part
c     b1=xlam*shp(1,k)
c     b2=xlam*shp(2,k)
      b1=0.
      b2=0.
      j1=1
      do 33 j=1,nel
c add this line
c	j1=1+(j-1)*ndfe
c
c
c     deviatoric operator    : suu
c     (dev-dev coupling)
c
c xj xk
c     *a1,shp(2,j),a2

      kel(j1,k1)=kel(j1,k1)+shp(1,j)*a1+shp(2,j)*a2
      kel(j1,k1)=kel(j1,k1)+(shp(1,j)*a1)/3.0
c xj yk
      kel(j1,k1+1)=kel(j1,k1+1)+0.
c     s(j1,k1+1)=s(j1,k1+1)+a1*shp(2,j)/3.0
      kel(j1,k1+1)=kel(j1,k1+1)-2.*a2*shp(1,j)/3.0+a1*shp(2,j)
c yj xk
      kel(j1+1,k1)=kel(j1+1,k1)+0.
c     s(j1+1,k1)=s(j1+1,k1)+a2*shp(1,j)/3.0
      kel(j1+1,k1)=kel(j1+1,k1)-2.*a1*shp(2,j)/3.0+a2*shp(1,j)
c yj yk
      kel(j1+1,k1+1)=kel(j1+1,k1+1)+shp(1,j)*a1+shp(2,j)*a2
      kel(j1+1,k1+1)=kel(j1+1,k1+1)+(shp(2,j)*a2)/3.0
c this if statement breaks the elegance of the code helas!
c     write(2,*)'end suu                            '
      if(k.eq.1)then
c
c
c     iso-dev   operator     : sup
c     (dev-is  coupling)
c
      do 333 mp=1,nelp
      kel(nstu+mp,j1)=kel(nstu+mp,j1)+xvol*shpp(mp)*shp(1,j)
      kel(nstu+mp,j1+1)=kel(nstu+mp,j1+1)+xvol*shpp(mp)*shp(2,j)
      kel(j1,nstu+mp)=kel(nstu+mp,j1)
      kel(j1+1,nstu+mp)=kel(nstu+mp,j1+1)
  333 continue

      	endif

      j1=j1+ndfe
   33 continue
      k1=k1+ndfe
   34 continue
c     write(2,*)'end sup                            '
c
c
c
c     solve iso-dev coupling at the element level :
c     elimination of internal dofs .here pressure.
c
c     if u-u convective term is not zero s is not symmetric
c     if u-p convective term is not zero s is not symmetric
c        u-p convective term arises from stress rate computations
   65 continue
c     write(2,*)'loop 65 terminated'

      if(nelp.eq.1)then
c     sigav(n,1)=(sigav(n,1)+ptemp-ptot(1))/lint
c     sigav(n,2)=(sigav(n,2)+ptemp-ptot(1))/lint
c     sigav(n,4)=(sigav(n,4)+ptemp-ptot(1))/lint
c     sigav(n,3)=sigav(n,3)/lint
      sigavl(n,1)=(sigavl(n,1)+ptemp)/lint
      sigavl(n,2)=(sigavl(n,2)+ptemp)/lint
      sigavl(n,4)=(sigavl(n,4)+ptemp)/lint
      sigavl(n,3)=sigavl(n,3)/lint
                   endif
      epsavl(n,1)=epsavl(n,1)/lint
      epsavl(n,2)=epsavl(n,2)/lint
      epsavl(n,3)=epsavl(n,3)/lint
      epsavl(n,4)=epsavl(n,4)/lint
c
c try hand elimination for n/1 elt
      if(nelp.eq.1)then
c     compress=volt/d(2)
c compressibility included in spp for possible nonlinear iterations
      	do i=1,nstu
      		do j=1,nstu
c     			s(i,j)=s(i,j)+s(i,9)*s(9,j)/compress
      			kel(i,j)=kel(i,j)+kel(i,9)*kel(9,j)/kel(9,9)
      		end do
      	end do
      endif


                   
cc      return
cc  end of elt03n
c
c call routine to calculate mass matrix entries to  rhs
c
      do 391 i=1,9
      do 390 j=1,9
      amass(i,j)=0.0
  390 continue
  391 continue
      deltl=(delt)
      isw=5
c replace with subr
c      call elt03n(inopredv,d,ul,xl,ix,amass,p,ndfe,ndm,nst,isw,deltl
c     *,nen,n,nel,viscl,bulkl,sbarl,pflagl,sigavl,epsavl,maxn,kstep)
      nelp=1
      l=d(4)
c      call pgauss(l,lint,sg,tg,wg)
c replace with subr
	pgg=1./dsqrt(3.0d0)
	lint=l*l
	do ig=1,4
	  sg(ig)=pgg*lr(ig)
	  tg(ig)=pgg*lz(ig)
	  wg(ig)=1.
	end do
c end of pgauss
      do 503 l=1,lint
c      call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
      do 2103 i=1,4
      shp(3,i)=(0.5+shps(i)*sg(l))*(0.5+shpt(i)*tg(l))
      shp(1,i)=shps(i)*(0.5+shpt(i)*tg(l))
      shp(2,i)=shpt(i)*(0.5+shps(i)*sg(l))
 2103 continue
      if(nel.ge.4)goto 2120
      do 2110 i=1,3
      shp(i,3)=shp(i,3)+shp(i,4)
 2110 continue
 2120 if(nel.gt.4)call shap2(sg(l),tg(l),shp,ix,nel)
      do 2132 i=1,ndm
      do 2131 j=1,2
      xs(i,j)=0.0
      do 2130 k=1,nel
      xs(i,j)=xs(i,j)+xl(i,k)*shp(j,k)
 2130 continue
 2131 continue
 2132 continue
      xsj=xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if(flg) goto 2141
      sx(1,1)=xs(2,2)/xsj
      sx(2,2)=xs(1,1)/xsj
      sx(1,2)=-xs(1,2)/xsj
      sx(2,1)=-xs(2,1)/xsj
      do 2140 i=1,nel
      tp=shp(1,i)*sx(1,1)+shp(2,i)*sx(2,1)
      shp(2,i)=shp(1,i)*sx(1,2)+shp(2,i)*sx(2,2)
      shp(1,i)=tp
 2140 continue
 2141 continue
c end of shape
c or any rho replacing d(3)!
      dvscal=wg(l)*xsj*d(3)
      j1=1
      do 500 j=1,nel
      w11=shp(3,j)*dvscal
      k1=j1
      do 510 k=j,nel
      amass(j1,k1)=amass(j1,k1)+shp(3,k)*w11
      k1=k1+ndfe
  510 continue
      j1=j1+ndfe
  500 continue
  503 continue
      nsl=nel*ndfe
      do 521 j=1,nsl,ndfe
      do 520 k=j,nsl,ndfe
      amass(j+1,k+1)=amass(j,k)
      amass(k,j)=amass(j,k)
      amass(k+1,j+1)=amass(j,k)
  520 continue
  521 continue
c end of elt03n
      do 303 i=1,8
      bel(i)=0.0
      do 300 j=2,8,2
      bel(i)=bel(i)+dble(amass(i,j))*g
  300 continue
  303 continue
c
c assemble global stiffness matrix and rhs
c
c
c  write elem s m
c
c     do 333 i=1,8
c     write(6,334)(kel(i,j),j=1,8)
c 334 format(8e10.4)
c 333 continue
      locrow=0
      do 60 l=1,4
      iglrow=ldf(node(l,iele))-1
      do 50 idf=1,2
      iglrow=iglrow+1
      locrow=locrow+1
      rhs(iglrow)=rhs(iglrow)+bel(locrow)
      loccol=0
      do 40 m=1,4
      iglcol=ldf(node(m,iele))-1
      do 30 jdf=1,2
      iglcol=iglcol+1
      loccol=loccol+1
      k=iglrow-iglcol+mbw
      abd(k,iglcol)=abd(k,iglcol)+dble(kel(locrow,loccol))
   30 continue
   40 continue
   50 continue
   60 continue
  100 continue
c!OMP end parallel do
      return
      end



c #################################################################
c ## dertermine the region of underplating                       ##
c #################################################################

      subroutine unplate(nrow,ncol,itst,nsing,ibegup,ibegmx)

      use dyn_arrays_mech
      use dyn_arrays
      implicit real*8 (a-h,o-z)

      ibegup=iunbeg
      if(itst.eq.1) ibegup=nsing-1 	
c      print*,'underplate=',ibegup,'   tstep=',itst
      return 
      end

C####################################################
C  Random number generator from numerical recipies ##
C####################################################

      function ran1(idum)
      integer idum,ia,im,iq,ir,ntab,ndiv
      real*8 ran1,am,eps,rnmx
      parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,
     *ntab=32,
     *eps=.00000012,
     *rnmx=1.-eps,
     *ndiv=1+(im-1)/ntab)

      integer j,k,iv(ntab),iy
      save iv,iy

      
      data iv /ntab*0/, iy /0/
      if(idum.le.o.or.iy.eq.0) then
      idum=max(-idum,1)
      do 11 j=ntab+8,1,-1
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if(idum.lt.0) idum=idum+im
      if(j.le.ntab) iv(j)=idum
  11  continue
      iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if(idum.lt.o) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)
      end


      
c***********************************************************************
c*                                                                     *
c*routine to apply frictional forces to global stiffness matrix and rhs*
c*                                                                     *
c***********************************************************************
      subroutine frictg(ne,lbw,lda,rforce2,ndf,nrow,nbn,delt)

      use dyn_arrays_mech
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      real*8 rforce2(2,*)
      mbw=2*lbw+1
c
c loop over boundary nodes
c
      do 320 ibn=1,nbn
      xforce=0.0
      yforce=0.0
c
c  identify node and degrees of freedom
c
      call gdf(nbase(ibn),nrow,ld,no)
      ld=ld+1
c
c  calculate vertical reaction force
c
      j1=max0(1,ld-lbw)
      j2=min0(ndf,ld+lbw)
      do 330 jb=j1,j2
      k=ld-jb+mbw
      yforce=yforce+abd(k,jb)*soln(jb)
  330 continue
      yforce=(yforce-rhs(ld))
      ld=ld-1
c
c  calculate horizontal reaction force
c
      j1=max0(1,ld-lbw)
      j2=min0(ndf,ld+lbw)
      do 340 jb=j1,j2
      k=ld-jb+mbw
      xforce=xforce+abd(k,jb)*soln(jb)
  340 continue
      xforce=(xforce-rhs(ld))
c
c calculate tangential force
c
      rforce2(1,ibn)=xforce*dcos(theta(ibn))+yforce*dsin(theta(ibn))
c
c calculate normal force
c
      rforce2(2,ibn)=-xforce*dsin(theta(ibn))+yforce*dcos(theta(ibn))
c
c  add term to global stiffness matrix for horizontal friction force
c
      abd(mbw,ld)=abd(mbw,ld)+vbound(ibn)/dcos(theta(ibn))
c
c  add term to rhs
c
      rhs(ld)=rhs(ld)+vbound(ibn)*basvel(ibn)
c
      write(6,666)ld,yforce,xforce,rforce2(2,ibn),rforce2(1,ibn)
666   format(i5,4e15.6)
c
  320 continue
      return 
      end
c*****************************************************************
c*                                                               *
c*  routine to determine degrees of freedom associated with node *
c*                                                               *
c*****************************************************************
      subroutine gdf(inode,nrow,ldf,nodf)
      ldf=2*inode-1
      nodf=2
      return
      end
c ***************************************************************
c *                                                             *
c *        routine to apply boundary conditions                 *
c *                                                             *
c ***************************************************************
c
      subroutine bc(numvbn,numpbn,ndf,lbw,lda,nrow,numsid,
     *nbn,upveln,nsing,ibegup,delt,rhoman,numvetbn,ncol)

      use dyn_arrays_mech
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      m=2*lbw+1
 
c  define a number large wrt stiffness components
      bv=0.0
      do i=1,ndf
      	if(abd(m,i).gt.bv)bv=abd(m,i)
      end do
      bv=bv*.1d5

c  apply constant pressure boundary condition to mass conser eqn
      if(numpbn.eq.0)go to 401
      do in=1,numpbn
      	inode=npnd(in)
      	call gdf(inode,nrow,ldf,nodf)
      	ldf=ldf+2
c  		set corresponding row of global stiffness matrix to 0
      	j1=max0(1,ldf-lbw)
      	j2=min0(ndf,ldf+lbw)
      	do jb=j1,j2
      		kb=ldf-jb+m
      		abd(kb,jb)=0.d0
        end do
c  		set principle diagonal component to large value
c       	and rhs to prescribed value
      	abd(m,ldf)=bv 
      	rhs(ldf)=bv*(bp(in))
      end do
  401 continue
 
c  apply boundary stresses to  equation of motion
      if(numsid.eq.0)go to 201
      print*,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
      print*,'ERROR: subroutine BC is not setup'
      print*,'           to handle loaded sides'
      print*,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
c      fac(1)=1./6.
c      fac(2)=2./3.
c      fac(3)=1./6.
c      do in=1,nbs
c      	slen=dsqrt((coord(1,nsnd(1,in))-coord(1,nsnd(3,in)))**2+
c     *	(coord(2,nsnd(1,in))-coord(2,nsnd(3,in)))**2)
c      	do isn=1,3
c      		inode=nsnd(isn,in)
c      		call gdf(inode,nrow,ldf,nodf)
c      		ldf=ldf+nsnd(4,in)-1
c      		rhs(ldf)=rhs(ldf)+(bside(in)*fac(isn)*slen)
c      	end do
c      end do
  201 continue

c  apply constant x,y velocity boundary conditions at edges
      if(numvbn.eq.0)go to 301
      do in=1,numvbn
      	inode=nvnd(1,in)
      	ldf=2*inode-1
      	nodf=2
      	ldf=ldf+nvnd(2,in)-1
c  		set principle diagonal component to large value
c       	and rhs to prescribed value
      	abd(m,ldf)=bv
      	rhs(ldf)=bv*(bvel(in))
      end do
  301 continue

c  apply constant tangent velocity boundary conditions at edges
      if(numvetbn.eq.0)go to 309
      do in=1,numvetbn
      	inode=nvtnd(in)
c		det tangent angle at base
      	if(mod(inode,nrow).eq.0) then
      		icol=inode/nrow
      	else
      		icol=floor(dble(inode)/dble(nrow))+1
      	endif	
      	if(icol.eq.1) then
      		dely=coord(2,(icol-1)*nrow+1)-coord(2,icol*nrow+1)
      		delx=coord(1,(icol-1)*nrow+1)-coord(1,icol*nrow+1)
      		delx2=coord(1,icol*nrow)-coord(1,(icol+1)*nrow)
      		dely2=coord(2,icol*nrow)-coord(2,(icol+1)*nrow)
      		ang1=atan(dely/delx)
      		ang2=atan(dely2/delx2)
      	else if(icol.eq.ncol) then
      		dely=coord(2,(icol-1)*nrow+1)-coord(2,(icol-2)*nrow+1)
      		delx=coord(1,(icol-1)*nrow+1)-coord(1,(icol-2)*nrow+1)
      		ang1=atan(dely/delx)
      	else 
      		print*,'####################################'
      		print*,'### ERROR: tangent edge vel BCs'
      		print*,'### 	are not being applied at the'
      		print*,'###     model edge. icol=',icol
      		print*,'####################################'
      		stop
      	endif	
c       velocity components at node      	
      	vytemp=sin(ang1)*bvelt(in)
      	vxtemp=cos(ang1)*bvelt(in)
c       apply velocities to stiffness matrix and rhs
c  			set principle diagonal component to large value
c       		and rhs to prescribed value
c       x vel      	
      	ldf=2*inode-1
      	abd(m,ldf)=bv
      	rhs(ldf)=bv*vxtemp
c      	y vel
      	ldf=ldf+1
      	abd(m,ldf)=bv
      	rhs(ldf)=bv*vytemp
      end do
  309 continue
  
 
c define constraint equation for basal surface
c
c  set tangential velocity condition on base
c		determine the x and y comp of vel from the basal
c		tangental velocity; unvel is an additional underplating
c		velocity added to the y vel
      if(nbn.eq.0)go to 501
      do ibn=1,nbn
      	if(ibn.eq.1)then
      		dely2=coord(2,nbase(ibn)+nrow)-coord(2,nbase(ibn))
      		delx2=coord(1,nbase(ibn)+nrow)-coord(1,nbase(ibn))
      		dely1=dely2
      		delx1=delx2
      	elseif(ibn.eq.nbn)then
      		dely1=coord(2,nbase(ibn))-coord(2,nbase(ibn)-nrow)
      		delx1=coord(1,nbase(ibn))-coord(1,nbase(ibn)-nrow)
      		dely2=dely1
      		delx2=delx1
      	else
      		dely2=coord(2,nbase(ibn)+nrow)-coord(2,nbase(ibn))
      		delx2=coord(1,nbase(ibn)+nrow)-coord(1,nbase(ibn))
      		dely1=coord(2,nbase(ibn))-coord(2,nbase(ibn)-nrow)
      		delx1=coord(1,nbase(ibn))-coord(1,nbase(ibn)-nrow)
      	endif
      	xlen1=dsqrt(delx1**2+dely1**2)
      	xlen2=dsqrt(delx2**2+dely2**2)
      	ang1=datan2(dely1,delx1)
      	ang2=datan2(dely2,delx2)
c 		xfnum is equivalent to dely1+dely2
c 		xfden == delx1+delx2
      	xfnum=xlen1*dsin(ang1)+xlen2*dsin(ang2)
      	xfden=xlen1*dcos(ang1)+xlen2*dcos(ang2)
      	thet=datan2(xfnum,xfden)
      	xlen3=(xlen2+xlen1)*.5
        if((ibn.ge.ibegup).and.(.not. (ibegup.eq.nsing-1)).and.
     $       (.not.(ibn .gt. nsing-1))) then
      		upvelx=upveln*dsin(-1.0*thet)
      		upvely=upveln*dcos(-1.0*thet)
      		basvelx=basvel(ibn)*dcos(thet)
      		basvely=basvel(ibn)*dsin(thet)
      		vxtemp=upvelx+basvelx
      		vytemp=upvely+basvely
c 			calcuate the additional vertical velocity needed in the mantle to 
c				have mass balance due to underplating
      		flx_nrm=upveln*xlen2
      		vz_therm=flx_nrm/delx2
      		unvel(ibn)=vz_therm
      	else 
  987 		continue      
      		unvel(ibn)=0
      		vxtemp=basvel(ibn)*dcos(thet)
      		vytemp=basvel(ibn)*dsin(thet)
      	endif	
      	call gdf(nbase(ibn),nrow,ldf,nodf)
c 		x component
c  		set principle diagonal component to large value
c       	and rhs to prescribed value
      	abd(m,ldf)=bv
      	rhs(ldf)=bv*vxtemp
c  		y component
c  		set principle diagonal component to large value
c       	and rhs to prescribed value
      	ldf=ldf+1
      	abd(m,ldf)=bv
      	rhs(ldf)=bv*vytemp
      end do
  501 continue
      return
      end

c********************************************************
c                                                       * 
c routine to add unvel to thermal vel field for lithos  *
c                                                       *
c********************************************************
      subroutine unplate_therm(nbn,nrowt,nrow)

      use dyn_arrays_mech
      use dyn_arrays
      implicit real*8 (a-h,o-z)

      do 210 i=1,nbn-1
      	if(unvel(i).eq.0.0) then
      		do 205 j=(i-1)*2*(nrowt-1)+1,i*2*(nrowt-1)-2*(nrow-1)
      			vz(j)=vz(j)
  205 		continue      			
      	else
      		do 206 j=(i-1)*2*(nrowt-1)+1,i*2*(nrowt-1)-2*(nrow-1)
      			vz(j)=unvel(i)/3.15578e13
  206 		continue      		
    	endif	
  210 continue
      end	
c*************************************************************
c*                                                           *
c*     routine to output results                             *
c*                                                           *
c*************************************************************
c
      subroutine output (nn,ne,itst,iter,nout,ttime,nout_t,nrow,
     *nbn,vrig,tstart,npoint,convel,ntsts,delt,nlrow,sealev,w_depth,
     *nbastrk,ibasflg,ninbas,ioutpt)

      use dyn_arrays
      use dyn_arrays_mech
      implicit real*8 (a-h,o-z)
      integer thdpl,fstpl,secpl,temp1
      character(30):: coord_op='coord_',vel_op='vel_',press_op='press_',
     *stress_xx_op='stress_xx_',stress_yy_op='stress_yy_',
     *stress_xy_op='stress_xy_',stress_zz_op='stress_zz_',
     *stress_secinv_op='stress_secinv_',stress_yield_op='stress_yield_',
     *stress_flag_op='stress_flag_',srate_xx_op='srate_xx_',
     *srate_yy_op='srate_yy_',srate_xy_op='srate_xy_',
     *srate_zz_op='srate_zz_',srate_dilt_op='srate_dilt_',
     *srate_secinv_op='srate_secinv_',lmesh_op='lmesh_',
     *temp_mech_op='temp_mech_',visc_elem_op='visc_elem_',
     *visc_gp_op='visc_gp_',erosion_op='erosion_',
     *temp_track_op='temp_track_',unvel_op='unvel_',exhum_op='exhum_',
     *sur_prof_op='sur_prof_',duc_flag_op='duc_flag_',
     *matp_phi_op='matp_phi_',matp_den_op='matp_den_',
     *matp_coh_op='matp_coh_',matp_prex_op='matp_prex_',
     *matp_vmin_op='matp_vmin_',matp_activ_op='matp_activ_',
     *matp_expon_op='matp_expon_',basinfill_op='basinfill_',
     *peakchop_op='peakchop_',basin_track_op='basin_track_',
     *l_temp_all_op='l_temp_all',dir,fextn
      character(10):: nums='0123456789'

      time=ttime
      write(6,107) itst,time,iter

c output directory
      dir='output/'

c
c output l-temp every nout_t timestep
c
      if(itst.eq.1) then
      	if(output_flags(37).eq.1) then
      		open(2,file=trim(dir)//trim(l_temp_all_op),position='rewind')
      		write(2,*)ntsts
      		write(2,*)nout
      		write(2,402)delt
      		write(2,*)nout_t
      		do i=1,npoint
      			write(2,402)tpoint(5,i)
      		end do
      		close(2)
      	endif	
      else
      	itest=mod(itst,nout_t)
      	if(itest.eq.0)then
      		if(output_flags(37).eq.1) then
      			open(2,file=dir//trim(l_temp_all_op),position='append')
      			do i=1,npoint
      				write(2,402)tpoint(5,i)
     			end do
     			close(2)
     		endif	
      	endif	
      endif	

      if(nout.eq.1)go to 5
      if(itst.eq.1)go to 5
      itest=mod(itst,nout)
      if(itest.ne.0)return
    5 continue

      write(6,108) itst,time

c determine extension for output file names
c 	track number of outputs      
      if(itst.eq.1) then
      	ioutpt=0
      endif	
      ioutpt=ioutpt+1
      if(ioutpt.lt.10) then
      	fextn=nums(ioutpt+1:ioutpt+1)
      elseif(ioutpt.lt.100) then
      	fstpl=(ioutpt)/10+1
      	secpl=(ioutpt-10*(fstpl-1))+1
      	fextn=nums(fstpl:fstpl)//nums(secpl:secpl)
      elseif(ioutpt.lt.1000) then
      	fstpl=(ioutpt)/100+1
      	temp1=(ioutpt-(ioutpt/100)*100)
      	secpl=temp1/10+1
      	thdpl=ioutpt-((fstpl-1)*100+(secpl-1)*10)+1
      	if(temp1.lt.10)secpl=1
      	fextn=nums(fstpl:fstpl)//nums(secpl:secpl)//
     *	nums(thdpl:thdpl)
      endif

c#############
c coord file
c#############
      if(output_flags(1).eq.1) then
        open(2,file=trim(dir)//trim(coord_op)//trim(fextn),
     *  position='rewind')
      	write(2,101)nn
      	do i=1,nn
      		write(2,102)coord(1,i),coord(2,i)
      	end do
      	close(2)
      endif    

c#############
c crustal velocity
c#############
      if(output_flags(2).eq.1) then
        open(3,file=trim(dir)//trim(vel_op)//trim(fextn),
     *  position='rewind')
      	write(3,101)nn
      	do i=1,nn
      		write(3,102)velx(i),vely(i)
      	end do
      	close(3)
      endif	

c#############
c stress files
c#############
c pressure
      if(output_flags(3).eq.1) then
        open(4,file=trim(dir)//trim(press_op)//trim(fextn),
     *  position='rewind')
      	write(4,101)ne
      	do i=1,ne
      		write(4,402)sbar(j,1)
      	end do
      	close(4)
      endif
c stress xx
      if(output_flags(4).eq.1) then
      	open(7,file=trim(dir)//trim(stress_xx_op)//trim(fextn),
     *	position='rewind')
      	write(7,101)ne
      	do i=1,ne
      		write(7,103)stress(i,1)
      	end do	
      	close(7)
      endif	
c stress yy
      if(output_flags(5).eq.1) then
      	open(8,file=trim(dir)//trim(stress_yy_op)//trim(fextn),
     *	position='rewind')
      	write(8,101)ne
      	do i=1,ne
      		write(8,103)stress(i,2)
      	end do	
      	close(8)
      endif      
c stress xy
      if(output_flags(6).eq.1) then
      	open(9,file=trim(dir)//trim(stress_xy_op)//trim(fextn),
     *	position='rewind')
      	write(9,101)ne
      	do i=1,ne
      		write(9,103)stress(i,3)
      	end do	
      	close(9)
      endif      
c stress zz
      if(output_flags(7).eq.1) then
      	open(10,file=trim(dir)//trim(stress_zz_op)//trim(fextn),
     *	position='rewind')
      	write(10,101)ne
      	do i=1,ne
      		write(10,103)stress(i,4)
      	end do	
      	close(10)
      endif      
c secinv
      if(output_flags(8).eq.1) then
      	open(11,file=trim(dir)//trim(stress_secinv_op)//trim(fextn),
     *	position='rewind')
      	write(11,101)ne
      	do i=1,ne
      		secinv=-stress(i,1)*stress(i,2)+stress(i,3)*stress(i,3)
      		if(secinv.lt.0.0)then
      			secinv=0.0
      		endif
      		write(11,103)dsqrt(secinv)
      	end do
      	close(11)
      endif      
c yield stress
      if(output_flags(9).eq.1) then
      	open(12,file=trim(dir)//trim(stress_yield_op)//trim(fextn),
     *	position='rewind')
      	write(12,101)ne
      	do i=1,ne
      		phi2=3.14159*phi(i)/180.
      		press=sbar(i,1)
      		if(press.lt.0.)press=0.0
      		cosphi=dcos(phi2)
      		sinphi=dsin(phi2)
      		yield=press*sinphi+coh(i)*cosphi
      		write(12,103)yield
      	end do
      	close(12)
      endif
c plasti failure flag
      if(output_flags(10).eq.1) then
      	open(13,file=trim(dir)//trim(stress_flag_op)//trim(fextn),
     *	position='rewind')
      	write(13,101)ne
      	write(13,101)(ipflag(i,1),i=1,ne)
      	close(13)
      endif	

c#############
c strain rates
c#############
c srate_xx
      if(output_flags(11).eq.1) then
      	open(14,file=trim(dir)//trim(srate_xx_op)//trim(fextn),
     *	position='rewind')
      	write(14,101)ne
      	do i=1,ne
      		write(14,444)srate(i,1)
      	end do	
      	close(14)
      endif
c srate_yy
      if(output_flags(12).eq.1) then
      	open(15,file=trim(dir)//trim(srate_yy_op)//trim(fextn),
     *	position='rewind')
      	write(15,101)ne
      	do i=1,ne
      		write(15,444)srate(i,2)
      	end do	
      	close(15)
      endif      
c srate_xy
      if(output_flags(13).eq.1) then
      	open(16,file=trim(dir)//trim(srate_xy_op)//trim(fextn),
     *	position='rewind')
      	write(16,101)ne
      	do i=1,ne
      		write(16,444)srate(i,3)
      	end do
      	close(16)
      endif      
c srate_zz
      if(output_flags(14).eq.1) then
      	open(17,file=trim(dir)//trim(srate_zz_op)//trim(fextn),
     *	position='rewind')
      	write(17,101)ne
      	do i=1,ne
      		write(17,444)srate(i,4)
      	end do	
      	close(17)
      endif       
c srate_dilt 
      if(output_flags(15).eq.1) then
      	open(18,file=trim(dir)//trim(srate_dilt_op)//trim(fextn),
     *	position='rewind')
      	write(18,101)ne
      	dilmax=0.0
      	dilav=0.0
      	dilav2=0.0
      	dilav3=0.0
      	dilav4=0.0
      	do i=1,ne
      		dilit=srate(i,1)+srate(i,2)
      		secdef=(srate(i,1)*srate(i,1)+srate(i,2)*srate(i,2))/2.
     *		+srate(i,3)*srate(i,3)
      		if(secdef.lt.0.0)then
      				secdef=0.0
      		endif
      		secdef=dsqrt(secdef)
      		ditest=dabs(dilit/secdef)
      		dilav=dilit+ditest
      		dilav2=dilav2+dabs(srate(i,1)/dilit)
      		dilav3=dilav3+dabs(srate(i,2)/dilit)
      		dilav4=dilav4+dabs(srate(i,3)/dilit)
      		if(ditest.gt.dilmax)then
      				dilmax=ditest
      				imax=i
      		endif
      		write(18,444)dilit
      	end do
      	close(18)
      endif	
      dilav=dilav/dble(ne)
      write(6,804)imax,dilmax,dilav
c srate_secinv
      if(output_flags(16).eq.1) then
      	open(19,file=trim(dir)//trim(srate_secinv_op)//trim(fextn),
     *	position='rewind')
      	write(19,101)ne
      	do i=1,ne
      		secdef=(srate(i,1)*srate(i,1)+srate(i,2)*srate(i,2))/2.
     *		+srate(i,3)*srate(i,3)
      		if(secdef.lt.0.0)then
      			secdef=0.0
      		endif
      		secdef=dsqrt(secdef)
      		write(19,444)secdef
      	end do
      	close(19)
      endif	

c#############
c lmesh coords
c#############
      if(output_flags(17).eq.1) then
      	open(20,file=trim(dir)//trim(lmesh_op)//trim(fextn),
     *	position='rewind')
      	write(20,101)nlrow
      	write(20,101)npoint
      	do i=1,npoint
      		write(20,102)tpoint(1,i),tpoint(2,i)
      	end do	
      	close(20)
      endif	

c#############
c crustal temps
c#############
      if(output_flags(18).eq.1) then
      	open(21,file=trim(dir)//trim(temp_mech_op)//trim(fextn),
     *	position='rewind')
      	write(21,101)ne
      	do i=1,ne
      		write(21,103)temptc(i)
      	end do	
      	close(21)
      endif	

c#############
c viscosity
c#############
c for element
      if(output_flags(19).eq.1) then
      	open(22,file=trim(dir)//trim(visc_elem_op)//trim(fextn),
     *	position='rewind')
      	write(22,101)ne
      	do i=1,ne
      		write(22,103)visc(i,1)
      	end do	
      	close(22)
      endif
c for gauss points
      if(output_flags(20).eq.1) then
      	open(23,file=trim(dir)//trim(visc_gp_op)//trim(fextn),
     *	position='rewind')
      	write(23,101)ne
      	do i=1,ne
      		write(23,103)(visc(i,j),j=1,4)
      	end do
      	close(23)
      endif	

c############
c surface erosion
c############
      if(output_flags(21).eq.1) then
      	open(24,file=trim(dir)//trim(erosion_op)//trim(fextn),
     *	position='rewind')
      	write(24,101)(nn/nrow)
      	nfree=0
      	do inc=nrow,nn,nrow
      		nfree=nfree+1
      		write(24,102)coord(1,inc),coord(2,inc),veros(1,nfree),
     *		veros(2,nfree)
      	end do
      	close(24)
      endif

c#############
c temp of lagrangian nodes, only at normal output interval (temp_track)
c#############
      if(output_flags(22).eq.1) then
      	open(25,file=trim(dir)//trim(temp_track_op)//trim(fextn),
     *	position='rewind')
      	write(25,101)npoint
      	write(25,402)(tpoint(5,j),j=1,npoint)
      	close(25)
      endif	

c#############
c underplating velocity
c#############
      if(output_flags(23).eq.1) then
      	open(26,file=trim(dir)//trim(unvel_op)//trim(fextn),
     *	position='rewind')
      	write(26,101) (nn/nrow)
      	nfree=0
      	do inc=1,nn,nrow
      		nfree=nfree+1
      		write(26,102)coord(1,inc),coord(2,inc),unvel(nfree)
      	end do
      	close(26)
      endif

c#############
c exhumation rate
c#############
      if(output_flags(24).eq.1) then
      	open(27,file=trim(dir)//trim(exhum_op)//trim(fextn),
     *	position='rewind')
		write(27,101)npoint
      	write(27,402)(exhum(j),j=1,npoint)
      	close(27)
      endif	

c#############
c surface profiles
c#############
      if(output_flags(25).eq.1) then
      	open(28,file=trim(dir)//trim(sur_prof_op)//trim(fextn),
     *	position='rewind')
      	write(28,101)(nn/nrow)
      	nfree=0
      	do inc=nrow,nn,nrow
      		nfree=nfree+1
      		write(28,102)xsur(1,nfree),vsur(2,nfree),xsur(2,nfree),
     *		rsur(2,nfree)
      	end do
      	close(28)
      endif

c#############
c ductile flag for lagrangian nodes
c#############
      if(output_flags(26).eq.1) then
      	open(29,file=trim(dir)//trim(duc_flag_op)//trim(fextn),
     *	position='rewind')
      	write(29,101)npoint
      	do j=1,npoint
      		write(29,403)int(tpoint(6,j)),int(tpoint(7,j))
      	end do
      	close(29)
      endif

c#############
c material props, mechanical
c#############
c phi
      if(output_flags(27).eq.1) then
      	open(30,file=trim(dir)//trim(matp_phi_op)//trim(fextn),
     *	position='rewind')
      	write(30,101)ne
      	do i=1,ne
      		write(30,104)phi(i)
      	end do
      	close(30)
      endif
c den
      if(output_flags(28).eq.1) then
      	open(31,file=trim(dir)//trim(matp_den_op)//trim(fextn),
     *	position='rewind')
      	write(31,101)ne
      	do i=1,ne
      		write(31,104)den(i)
      	end do
      	close(31)
      endif
c coh
      if(output_flags(29).eq.1) then
      	open(32,file=trim(dir)//trim(matp_coh_op)//trim(fextn),
     *	position='rewind')
      	write(32,101)ne
      	do i=1,ne
      		write(32,104)coh(i)
      	end do	
      	close(32)
      endif      
c pre exponential
      if(output_flags(30).eq.1) then
      	open(33,file=trim(dir)//trim(matp_prex_op)//trim(fextn),
     *	position='rewind')
      	write(33,101)ne
      	do i=1,ne
      		write(33,104)prex(i)
      	end do	
      	close(33)
      endif
c min viscosity
      if(output_flags(31).eq.1) then
      	open(34,file=trim(dir)//trim(matp_vmin_op)//trim(fextn),
     *	position='rewind')
      	write(34,101)ne
      	do i=1,ne
      		write(34,104)vmin(i)
      	end do	
      	close(34)
      endif
c activation energy
      if(output_flags(32).eq.1) then
      	open(35,file=trim(dir)//trim(matp_activ_op)//trim(fextn),
     *	position='rewind')
      	write(35,101)ne
      	do i=1,ne
      		write(35,104)q(i)
      	end do	
      	close(35)
      endif
c powerlaw exponent
      if(output_flags(33).eq.1) then
      	open(36,file=trim(dir)//trim(matp_expon_op)//trim(fextn),
     *	position='rewind')
      	write(36,101)ne
      	do i=1,ne
      		write(36,104)expn(i)
      	end do	
      	close(36)
      endif      

c#############
c basinfill
c#############
      if(output_flags(34).eq.1) then
      	open(37,file=trim(dir)//trim(basinfill_op)//trim(fextn),
     *	position='rewind')
      	do i=1,nn/nrow
      		write(37,129)nbasinfill(i),coord(1,i*nrow),basinfill(i)
      	end do
      	close(37)
      endif	

c#############
c peakchop
c#############
      if(output_flags(35).eq.1) then
      	open(38,file=trim(dir)//trim(peakchop_op)//trim(fextn),
     *	position='rewind')
      	do i=1,nn/nrow
      		write(38,129)npeakchop(i),coord(1,i*nrow),peakchop(i)
      	end do
      	close(38)
      endif	

c#############
c basin surfaces
c#############
      if(output_flags(36).eq.1) then
      	open(39,file=trim(dir)//trim(basin_track_op)//trim(fextn),
     *	position='rewind')
      	if(ibasflg.eq.1) then
      		write(39,101)nbastrk
      		if(nbastrk.gt.0) write(39,101)ninbas
      		do i=1,nbastrk
      			write(39,101)ibastrk(2,i)-ibastrk(1,i)+1
      			do j=ibastrk(1,i),ibastrk(2,i)
      				write(39,102)bastrk(1,j),bastrk(2,j)
      			end do	
      		end do	
      	else
      		write(39,101)0
      	endif
      	close(39)
      endif

c#############
c number of outputs
c#############
      open(3,file=trim(dir)//'n_out',position='rewind')
      write(3,101)ioutpt
      close(3)
      

  101 format(i9)
  444 format(6e26.16)
  102 format(SP,6e12.6,/5e12.6)
  402 format(e20.9)
  403 format(2i4)
  104 format(6e12.6)
  103 format(7e11.5,i3)
  107 format(/(2x,'time step ',i5,' at time',e12.5,
     *' converged after',i5,' iterations')/)
  108 format(//(5x,'*****   output at time step ',i5,' at time',e12.5
     *,'   *****')//)
  111 format(2f8.1)
  129 format(i5,2e15.7)
  804 format(' max dilitation ', i5,e15.6,' average dil ',e15.6)

      return
      end
c*****************************************************************
c*                                                               *
c*    routine to check for convergence of each timestep          *
c*                                                               *
c*****************************************************************
      subroutine conver(velx,vely,rhs,toler,icflag,nrow,ncol,ndf,nn,
     *coord)
      implicit real*8 (a-h,o-z)
      dimension velx(*),vely(*),coord(2,*)
      real*8 rhs(*)

      vmax=0.0
      mnode=1

      do 100 inode=1,nn
      	k=inode*2-1
      	test=dabs(velx(inode)-rhs(k))
      	if(test.gt.vmax)then
      		vmax=test
      		mnode=inode
      	endif
      	velx(inode)=rhs(k)
  100 continue

      ivf=1
      do 200 inode=1,nn
      	k=inode*2
      	test=dabs(vely(inode)-rhs(k))
      	if(test.gt.vmax)then
      		vmax=test
      		mnode=inode
      		ivf=2
      	endif
      	vely(inode)=rhs(k)
  200 continue

      icflag=0
      if(vmax.lt.toler)then
      	icflag=1
      endif

      if(ivf.eq.1)then
      	write(6,101)vmax,mnode,velx(mnode)
      	write(6,109)int(mnode/nrow),mnode-int(mnode/nrow)*nrow,
     *	coord(1,mnode) 
      else
      	write(6,102)vmax,mnode,vely(mnode)
      	write(6,109)int(mnode/nrow),mnode-int(mnode/nrow)*nrow,
     *	coord(1,mnode)
      endif


  101 format('              vel norm= ',e12.6,'  x vel at node: ',i5,
     *2e12.4)
  109 format('              colm=',i4,'  row=',i4,'  x-pos= ',e12.6)      
  102 format('              vel norm= ',e12.6,'  y vel at node: ',i5,
     *2e12.4)  
      return
      end
c**********************************************************************
c*                                                                    *
c*    routine to initialize stress field                              *
c*                                                                    *
c**********************************************************************
      subroutine sinit(nrow,ne)

      use dyn_arrays
      use dyn_arrays_mech
      implicit real*8 (a-h,o-z)
      g=9.8
      ninc=nrow-1
      do 100 ie=1,ne 
      itop=node(1,ie)-mod(node(1,ie),nrow)+nrow
      depth=(coord(2,itop)+coord(2,itop+nrow))/2.0d0-
     *(coord(2,node(1,ie))
     *+coord(2,node(2,ie))+coord(2,node(3,ie))+
     *coord(2,node(4,ie)))/4.0d0
c
c  set initial stress to den(ie)
c
c     sprev(inode)=den(ie)*g*depth
c     sprev(inode)=0.0d0
      sbar(ie,1)=den(ie)*g*depth
  100 continue
      return
      end
c***********************************************************************
c*                                                                     *
c*    routine to reposition free surface and update mesh               *
c*    reassign temperatures to new Eulerian positions                  *
c*                                                                     *
c***********************************************************************
c
      subroutine remesh(delt,nn,nrow,rhoman,rhof,ncom,ne,vrig,npoint,
     *nnt,nrowt,net,nsing,convel,itst,erosl,erosr,
     *peros,rpow,sealev,slpmax,npad,np1,np2,prig,rrig,nsing1,ctoler,
     *wtoler,tmax,rhoavinitl,rhoavinitr,dy_flex_init1,dy_flex_init2,
     *wdepth,nsthick,plthick,leqflag,iplasflg,iblay,iblayt,isedl,isedr,
     *ibasflg,intmrkb,nbastrk,nbastary,nbastind,ninbas,ipkfill,ibasfill,
     *sedmax)
      
      use dyn_arrays
      use dyn_arrays_mech
      implicit real*8 (a-h,o-z)
      integer count
      real(kind=8),allocatable::yp1prev(:),yp2prev(:),cnew(:),
     *fnode(:),vdiffnew(:),rdiffnew(:),xold(:,:)

      g=9.8
      ninc=nrow
      numfre=nn/ninc
      ncol=numfre

c allocate arrays
      allocate(cnew(ncol),fnode(ncol),vdiffnew(ncol),rdiffnew(ncol),
     *xold(2,ncol))
      vdiffnew=0.0
      rdiffnew=0.0
      fnode=0.0
      xold=0.0
      cnew=0.0

c set all n-1 temperatures to n temperature at current euler position
      do i=1,nnt
      	told(i)=tempt(i)
      end do
c  Make a temp array for just the crust
      ncrustbeg=nrowt-nrow
      l=0
      count=0
      do i=ncrustbeg,nnt,nrowt
     	l=(i-ncrustbeg)/nrowt+1
        do j=1,nrow
        	k=i+j
        	count=((l-1)*nrow)+j
        	toldc(count)=told(k)
        	temptc(count)=tempt(k)
      	end do
      end do

c interpolate vdiff and rdiff from previous tstep to the 
c	xcoords used in this time step, ie account for the
c	remesh of the upper surface
      do nfree=ninc,nn,ninc
      	i=(nfree-ninc)/ninc+1
      	xold(1,i)=coord(1,nfree)
      	xold(2,i)=coord(2,nfree)
      end do
      do i=1,numfre
      	do j=1,numfre
      		if(vdiff(1,j).gt.xold(1,i)) then
      			if(j.eq.1) then
      				vslope=(vdiff(2,j+1)-vdiff(2,j))/
     *				  (vdiff(1,j+1)-vdiff(1,j))
     				delxv=vdiff(1,j)-xold(1,i)
     				vdiffnew(i)=vdiff(2,j)-delxv*vslope
     				rslope=(rdiff(2,j+1)-rdiff(2,j))/
     *				  (rdiff(1,j+1)-rdiff(1,j))
     				delxr=rdiff(1,j)-xold(1,i)
     				rdiffnew(i)=rdiff(2,j)-delxr*rslope
     				goto 929
     			else
     				vslope=(vdiff(2,j)-vdiff(2,j-1))/
     *                (vdiff(1,j)-vdiff(1,j-1))
     				delxv=vdiff(1,j)-xold(1,i)
     				vdiffnew(i)=vdiff(2,j)-delxv*vslope
     				rslope=(rdiff(2,j)-rdiff(2,j-1))/
     *                (rdiff(1,j)-rdiff(1,j-1))
     				delxr=rdiff(1,j)-xold(1,i)
     				rdiffnew(i)=rdiff(2,j)-delxr*rslope
     				go to 929
     			end if
     		else if(vdiff(1,j).eq.xold(1,i)) then
     			vdiffnew(i)=vdiff(2,j)
     			rdiffnew(i)=rdiff(2,j)
     		end if	
     	end do
  929 continue     	
      end do	

      do i=1,numfre
      	vdiff(2,i)=vdiffnew(i)
      	rdiff(2,i)=rdiffnew(i)
      end do 	

      do nfree=ninc,nn,ninc
      	i=(nfree-ninc)/ninc+1
      	xsurold(1,i)=xsur(1,i)
      	xsurold(2,i)=xsur(2,i)
      	xsur(1,i)=coord(1,nfree)+delt*velx(nfree)
      	xsur(2,i)=coord(2,nfree)+delt*vely(nfree)
      	vsur(1,i)=xsur(1,i)
      	rsur(1,i)=xsur(1,i)
      	vsur(2,i)=(coord(2,nfree)-vdiff(2,i))+delt*vely(nfree)
      	rsur(2,i)=(coord(2,nfree)-rdiff(2,i))+delt*vely(nfree)
      end do
c
c call routine to calculate erosion
      call erosion(nn,nrow,ninc,numfre,delt,itst,erosl,erosr,peros,
     *rpow,sealev,w_depth,isedl,isedr,ibasflg,intmrkb,nbastrk,
     *nbastary,nbastind,ninbas,ipkfill,ibasfill,
     *sedmax)

c apply filter from left to right to check that no slopes exceed
c 	a maximum value (slpmax)
c  	  check negative slopes
      icount=0
   61 islope=0
      icount=icount+1
      do i=2,numfre
      	slope= (-xsur(2,i)+xsur(2,i-1))/(xsur(1,i)-xsur(1,i-1))
      	if(slope.gt.slpmax)then
      		islope=1
      		diff=-xsur(2,i)+xsur(2,i-1)
      		diffm=slpmax*(xsur(1,i)-xsur(1,i-1))
      		npeakchop(i)=npeakchop(i)+1
      		peakchop(i)=peakchop(i)+.5*(diff-diffm)
      		peakchop(i-1)=peakchop(i-1)-.5*(diff-diffm)
      		xsur(2,i)=xsur(2,i)+.5*(diff-diffm)
      		xsur(2,i-1)=xsur(2,i-1)-.5*(diff-diffm)
      	endif
      end do
      if(islope.eq.1.and.icount.lt.999)go to 61
c     check positive slopes
      icount=0
   63 islope=0
      icount=icount+1
      do i=2,numfre
      	slope= (xsur(2,i)-xsur(2,i-1))/(xsur(1,i)-xsur(1,i-1))
      	if(slope.gt.slpmax)then
      		islope=1
      		diff=xsur(2,i)-xsur(2,i-1)
      		diffm=slpmax*(xsur(1,i)-xsur(1,i-1))
      		npeakchop(i)=npeakchop(i)+1
      		peakchop(i)=peakchop(i)-.5*(diff-diffm)
      		peakchop(i+1)=peakchop(i+1)+.5*(diff-diffm)
      		xsur(2,i)=xsur(2,i)-.5*(diff-diffm)
      		xsur(2,i-1)=xsur(2,i-1)+.5*(diff-diffm)
      	endif
      end do
      if(islope.eq.1.and.icount.lt.999)go to 63

c  force boundary nodes to move verticaly
      xsur(1,1)=coord(1,nrow)
      xsur(1,numfre)=coord(1,nn)
c	added 5-24-03
c	also force boundary nodes to hold their y-position,
      xsur(2,1)=coord(2,nrow)
      xsur(2,numfre)=coord(2,nn)

c  loop over columns
      ymax=0.0
      do 60 jfree=ninc,nn,ninc
      ibase=jfree-ninc+1+iblay
c
c  interpolate surface coords to local vector cnew
c     xsur is surface from old surface + vel*time
c     interpolate coord(2,) to the new y pos in xsur
c
      do k=2,numfre
      	if(coord(1,jfree).ge.xsur(1,k-1).and.coord(1,jfree).le.
     *	xsur(1,k))then
      		cnew(ninc)=xsur(2,k-1)+((coord(1,jfree)-xsur(1,k-1))/
     *		(xsur(1,k)-xsur(1,k-1)))*(xsur(2,k)-xsur(2,k-1))
      		ytest=dabs(cnew(ninc)-coord(2,jfree))
      		ymax=max(ymax,ytest)
      		go to 26
      	endif
      end do
      write(6,*)' error in remesh: upper surface node not repositioned '
      cnew(ninc)=xsur(2,1)
   26 continue

c  interpolate internal coords in a column to local vector cnew
      do i=2+iblay,ninc-1-iblayt
      	cnew(i)=coord(2,ibase)+dble(float(i-1-iblay)/
     *	float(ninc-1-iblay-iblayt))
     *	*(cnew(ninc)-iblayt*(coord(2,jfree)-
     *	coord(2,jfree-1))-coord(2,ibase))
      end do

c  assign new coords to points in boundary layer
      do i=0,iblay
      	cnew(i+1)=coord(2,jfree-ninc+i+1)
      end do

c loop for new coords in upper boundary layer
      do i=0,iblayt
      	j=ninc-iblayt+i
      	cnew(j)=cnew(ninc)-iblayt*(coord(2,jfree)-coord(2,jfree-1))+
     *	i*(coord(2,jfree)-coord(2,jfree-1))
      end do

c interpolate temperatures
      do i=2+iblay,ninc-1
      	inode=jfree-ninc+i
c     find first node above  point of interest
      	do j=jfree-ninc+1,jfree
      		if(cnew(i).lt.coord(2,j))then
      			nup=j
      			nbelow=j-1
      			go to 43
      		endif
      	end do
      	write(6,*)'error in remesh temperature not interpolated'
      	write(6,*) inode
      	print*,cnew(i),coord(2,j),coord(2,jfree-ninc+1)
c     interpolate temperature
   43 	z1=coord(2,nup)-cnew(i)
      	z2=cnew(i)-coord(2,nbelow)
      	w1=z1/(z1+z2)
      	w2=z2/(z1+z2)
      	toldc(inode)=w2*temptc(nup)+w1*temptc(nbelow)
      end do

c  interpolate remaining nodal coords in the column
      do i=2+iblay,ninc
      	inode=jfree-ninc+i
      	coord(2,inode)=cnew(i)
      end do
   60 continue
      write(6,*)'  max displacement in remesh = ',ymax

c  add toldc (remeshed crustal told) to told (thermal told)
      l=0
      count=0
      do m=ncrustbeg,nnt,nrowt
	  	l=(m-ncrustbeg)/nrowt +1
        do j=1,nrow
          k=m+j
          count=((l-1)*nrow)+j
          told(k)=toldc(count)
        end do
      end do

c##############################
c Flexure/isostacy calculation
c##############################

c  offset each column by isostatic displacement
      if(ncom.eq.0)then
c     local airy isostacy
      	ncol=nn/nrow
      do 351 icol=1,ncol
c  		overwrite ziso with simple column height calculation
      	itop=icol*(nrow)
      	ibot=itop-nrow+1
      	zbase=coord(2,ibot)
      	stiso=(rhof)*((coord(2,itop)-coord(2,ibot)-zinit(icol)))
      	ziso(icol)=zeq(icol)-stiso/rhoman
      	do 350 irow=1,nrow
      		inode=(icol-1)*nrow+irow
      		coord(2,inode)=coord(2,inode)+ziso(icol)-zbase
  350 	continue
  351 continue
      elseif(ncom.eq.1) then
c     1 plate flexural compensation
      	ncol=nn/nrow
      	do icol=1,ncol
      		itop=icol*nrow
      		ibot=itop-nrow+1
c  			find length over which force operates
      		if(icol.eq.1)then
      			slen=(coord(1,ibot+nrow)-coord(1,ibot))/2.0d0
      		elseif(icol.eq.ncol)then
      			slen=(coord(1,ibot)-coord(1,ibot-nrow))/2.0d0
      		else
      			slen=(coord(1,ibot+nrow)-coord(1,ibot-nrow))/2.0d0
      		endif
      		fnode(icol)=slen*g*(rhof)*((coord(2,itop)-
     *		coord(2,ibot)-zinit(icol)))
      	end do

c define flexural parameter
      	alpha=((4.0d0*prig)/((rhoman-rhof)*g))**.25
      	alph2=(alpha**3)/(8.0d0*prig)
      	do icol=1,ncol
      		itop=icol*nrow
      		ibot=itop-nrow+1
      		w=0.0d0
      		do icol2=1,ncol
      			itop2=icol2*nrow
      			ibot2=itop2-nrow+1
      			dist=dabs(coord(1,ibot)-coord(1,ibot2))
      			dist=dist/alpha
      			w=w+fnode(icol2)*alph2*dexp(-dist)*(dcos(dist)
     *			+dsin(dist))
      		end do
      		zbase=coord(2,ibot)
      		do irow=1,nrow
      			inode=(icol-1)*(nrow)+irow
      			coord(2,inode)=coord(2,inode)+zeq(icol)-zbase-w
        	end do
c  			save incremental drops
      		zinc(1,icol)=coord(1,ibot)
      		zinc(2,icol)=coord(2,ibot)-zbase
      	end do
      elseif(ncom.eq.2) then
c     2 plate flexural compensation
c		since superimposing changes to current profile, the original
c		sub loads and moments used to get the subduction profile should
c		not be reapplied, so set them to zero here.
      	sload=0.0
      	smomen=0.0
      	allocate(yp1prev(np1),yp2prev(np2))
      	yp1prev=yp1
      	yp2prev=yp2
c		calculate flexure      	
      	call calc_flex_remesh(nrow,numfre,np1,np2,prig,rrig,rhoman,
     *	npad,nsing,nsing1,ctoler,sload,smomen,wtoler,rhof,g,rhoavinitl,
     *	rhoavinitr,dy_flex_init1,dy_flex_init2,wdepth,sealev,nn,itst)
c       shift eulerian coords
c       plate 1 
        do i=nsing1,np1-npad
        	icol=nsing-i+nsing1
        	zinc(1,icol)=coord(1,icol*nrow)
        	zinc(2,icol)=-(yp1(i)-yp1prev(i))
        	do j=1,nrow
        		coord(2,(icol-1)*nrow+j)=coord(2,(icol-1)*nrow+j)
     *  		+zinc(2,icol)
     		end do
     	end do	
c     	plate 2
      	do i=2,np2-npad
      		icol=nsing+i-1
      		zinc(1,icol)=coord(1,icol*nrow)
        	zinc(2,icol)=-(yp2(i)-yp2prev(i))
      		do j=1,nrow
      			coord(2,(icol-1)*nrow+j)=coord(2,(icol-1)*nrow+j)
     *			+zinc(2,icol)
      		end do
      	end do	
      	deallocate(yp1prev,yp2prev)
      endif

c subside lagrangian track points
      	do ipoint=1,npoint
      		do k=2,numfre
      			if(tpoint(1,ipoint).le.zinc(1,1))then
      				tpoint(2,ipoint)=tpoint(2,ipoint)+zinc(2,1)
      				go to 326
      			elseif(tpoint(1,ipoint).ge.zinc(1,k-1).and.
     *			tpoint(1,ipoint).le.zinc(1,k))then
      				tpoint(2,ipoint)=tpoint(2,ipoint)+zinc(2,k-1)+
     *				((tpoint(1,ipoint)-zinc(1,k-1))/(zinc(1,k)-
     *				zinc(1,k-1)))*(zinc(2,k)-zinc(2,k-1))
      				go to 326
      			endif
      		end do
  326 		continue
      	end do

c subside tracked basin surfaces
      	do ipoint=1,ninbas
      		do k=2,numfre
      			if(bastrk(1,ipoint).le.zinc(1,1))then
      				bastrk(2,ipoint)=bastrk(2,ipoint)+zinc(2,1)
      				go to 328
      			elseif(bastrk(1,ipoint).ge.zinc(1,k-1).and.
     *			bastrk(1,ipoint).le.zinc(1,k))then
      				bastrk(2,ipoint)=bastrk(2,ipoint)+zinc(2,k-1)+
     *				((bastrk(1,ipoint)-zinc(1,k-1))/(zinc(1,k)-
     *				zinc(1,k-1)))*(zinc(2,k)-zinc(2,k-1))
      				go to 328
      			endif
      		end do
  328 		continue
      	end do


c thermal remeshing routine
      call tmesh(nrow,ncol,nrowt,net,nnt,nn,nsing,convel,nsthick,
     *plthick,itst)

c  calculate new power law viscosity
      do ie=1,ne
      	bl=dexp(q(ie)/(8.3144d0*tmax))/(vmin(ie))
      	epsinv=(srate(ie,1)**2+srate(ie,2)**2)/2.+srate(ie,3)*
     *	srate(ie,3)
      	if(epsinv.lt.0.)epsinv=0.0d0
      	epsinv=dsqrt(epsinv)     
      	tele=(temptc(node(1,ie))+temptc(node(2,ie))+temptc(node(3,ie))
     *	+temptc(node(4,ie)))/4.
c 		POWER-LAW VISCOSITY
      	if(leqflag.ne.1) then
      		vpow2=(prex(ie)*dexp(-q(ie)/(8.3144d0*tele)))**(-1./expn(ie))
      		vpow=vpow2*(epsinv**(1./expn(ie)-1.))
      	else	
c		LINEAR VISCOUS
        	vpow=dexp(q(ie)/(8.3144d0*tele))/bl
        	vpow2=vpow
        endif	
      	if(vpow.gt.vrig)vpow=vrig
      	if(vpow.lt.vmin(ie))vpow=vmin(ie)
      	vpower(1,ie)=vpow
      	vpower(2,ie)=vpow2
c OVERWRITE FOR PLASTIC CASE      	
      	if(iplasflg.eq.1) then
      		vpower(1,ie)=vrig
      	 	vpower(2,ie)=vrig
      	endif	 
      end do

      deallocate(cnew,fnode,vdiffnew,rdiffnew,xold)

      return
      end
c**********************************************************************
c*   Routine to Remesh the Thermal Mesh                               *
c**********************************************************************
      subroutine tmesh(nrow,ncol,nrowt,net,nnt,nn,nsing,
     *convel,nsthick,plthick,itst)

      use dyn_arrays
      use dyn_arrays_mech
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer count,vcol,rcount,tricount
      real(kind=8),allocatable::cdown(:),cpresent(:),sstop(:,:),
     *ccbot(:,:),rlthick(:),tangle(:),cangle(:),fluxman(:),rwidth(:),
     *vman(:),sstopold(:),dysstop(:),ssbot(:),ssbotold(:),dyssbot(:),
     *wwidth(:),zadd(:),bvelm(:),flux(:),snewthick(:)

c allocate space for arrays
      allocate(cdown(ncol),cpresent(ncol),sstop(2,nsing+nrowt-nrow))
      allocate(ccbot(2,ncol),rlthick(ncol),tangle(ncol),cangle(ncol))
      allocate(vman(ncol),fluxman(ncol),rwidth(ncol),ssbot(nsing+nrowt
     *-nrow-nsthick),ssbotold(nsing+nrowt-nrow-nsthick))
      allocate(sstopold(nsing+nrowt-nrow),dysstop(nsing+nrowt-nrow))
      allocate(dyssbot(nsing+nrowt-nrow-nsthick))
      allocate(wwidth(nsing+nrowt-nrow-1),zadd(ncol),bvelm(ncol),
     *flux(ncol),snewthick(ncol))

      cdown=0.0
      cpresent=0.0
      sstop=0.0
      ccbot=0.0
      rlthick=0.0
      tangle=0.0
      cangle=0.0
      vman=0.0
      fluxman=0.0
      rwidth=0.0
      ssbot=0.0
      ssbotold=0.0
      sstopold=0.0
      dysstop=0.0
      dyssbot=0.0
      wwidth=0.0
      zadd=0.0
      bvelm=0.0
      flux=0.0
      snewthick=0.0
      
      nmesh=0

c store new mech model base position
      count=0
      do i=1,nn,nrow
      	count=(i-1)/nrow +1
      	cdown(count)=coord(2,i)
      end do
c store previous tst mech model base pos      
      count=0
      do j=(nrowt-nrow+1),nnt,nrowt
      	count=(j-(nrowt-nrow+1))/nrowt+1
      	cpresent(count)=coordt(2,j)
      end do
c change in mech base position over tst
      do k=1,ncol
      	cbase(k)=(cdown(k)-cpresent(k))
      end do
c store previous position of slab top
      do i=1,nsing
      	sstopold(i)=cpresent(i)
      end do	
c  Remesh Thermal Mesh within the mech model
      count=0
      k=0
      do i=(nrowt-nrow+1),nnt,nrowt
        do j=1,nrow
        	k=i+j-1
 	    	count=((i-(nrowt-nrow+1))/nrowt)*nrow+j
        	coordt(2,k)=coord(2,count)
        	nmesh=nmesh+1
      	end do
      end do

c remesh thermal mesh in retro lith and asthen
c	portion above slab
      do i=1,nrowt-nrow
      	icol=nsing+i
      	istart=(nsing+i)*nrowt-nrow-i+1
      	istop=icol*nrowt-nrow
c      	store previous pos of slab top
     	sstopold(icol)=coordt(2,istart)
      	do jnode=istart,istop
      		coordt(2,jnode)=coordt(2,jnode)+cbase(icol)
        	nmesh=nmesh+1
      	end do	
c     	store positions of slab top
      	sstop(1,icol)=coordt(1,istart)      	
      	sstop(2,icol)=coordt(2,istart)
c     	store thickness of retro lith and athen      	
      	rlthick(icol)=coordt(2,istop+1)-coordt(2,istart)
      end do
c	rest of retro lith and asthen
      do i=1,ncol-(nsing+nrowt-nrow)
      	icol=nsing+nrowt-nrow+i
      	istart=(icol-1)*nrowt+1
      	istop=icol*nrowt-nrow
      	do jnode=istart,istop
      		coordt(2,jnode)=coordt(2,jnode)+cbase(icol)
        	nmesh=nmesh+1
      	end do
c     	store thickness of retro lith and athen      	
      	rlthick(icol)=coordt(2,istop+1)-coordt(2,istart)
      end do	

c store remaining positions of slab top and beg of crust base
      do i=1,nsing
      	sstop(1,i)=coordt(1,i*nrowt-nrow+1)
      	sstop(2,i)=coordt(2,i*nrowt-nrow+1)
      	ccbot(1,i)=coordt(1,i*nrowt-nrow+1)
      	ccbot(2,i)=coordt(2,i*nrowt-nrow+1)
      end do	
c store the retro portion of crust base
      do i=nsing+1,ncol
      	ccbot(1,i)=coordt(1,i*nrowt-nrow+1)
      	ccbot(2,i)=coordt(2,i*nrowt-nrow+1)
      end do	

c calculate the change in slab top height
      do i=1,nsing+nrowt-nrow
      	dysstop(i)=sstop(2,i)-sstopold(i)
      end do	

c number of colms in the slab (from model lhs to base)
      nslabcol=nsing+nrowt-nrow

c  Calculate the tangent angle at each slabtop (sstop) 
c 	point from lhs to end of slab
      do j=2,nslabcol-1
        xdif1=abs(sstop(1,j-1)-sstop(1,j))
        zdif1=abs(sstop(2,j-1)-sstop(2,j))
        ang1=atan(zdif1/xdif1) 
        xdif2=abs(sstop(1,j)-sstop(1,j+1))
        zdif2=abs(sstop(2,j)-sstop(2,j+1))
        ang2=atan(zdif2/xdif2)
        tangle(j)=(ang1-ang2)/2+ang2
      end do
      xdif1=abs(sstop(1,1)-sstop(1,2))
      zdif1=abs(sstop(2,1)-sstop(2,2))
      tangle(1)=atan(zdif1/xdif1)
      xdif1=abs(sstop(1,nslabcol-1)-sstop(1,nslabcol))
      zdif1=abs(sstop(2,nslabcol-1)-sstop(2,nslabcol))
      tangle(nslabcol)=atan(zdif1/xdif1)
c calculate the tangent angle of crust base
      do j=2,ncol-1
      	xcdif1=abs(ccbot(1,j-1)-ccbot(1,j))
        zcdif1=abs(ccbot(2,j-1)-ccbot(2,j))
        cang1=atan(zcdif1/xcdif1)
        xcdif2=abs(ccbot(1,j)-ccbot(1,j+1))
        zcdif2=abs(ccbot(2,j)-ccbot(2,j+1))
        cang2=atan(zcdif2/xcdif2)
        cangle(j)=(cang1-cang2)/2+cang2
      end do
      cangle(1)=tangle(1)
      xcdif1=abs(ccbot(1,ncol-1)-ccbot(1,ncol))
      zcdif1=abs(ccbot(2,ncol-1)-ccbot(2,ncol))
      cangle(ncol)=atan(zcdif1/xcdif1)

c number of nodes in region up to the base of the slab
      nslength=nrowt*nslabcol

c Calculate the Velocity for the retro-lithosphere
c	left over from ablation code. since ablation is not implemented 
c	here, set rlvelx to zero. will need to chnage the setting of velocties
c	to zero below to allow for the application of vel in the reto lith
c   Also, I have commented out the portions of the code that remesh the lith
c	to take into account the change in lith thickness due to ablative 
c	velocities.  Will need to put these back in and redo them if ablation is
c	included.  also commented out where bvelm was adjusted for ablation
      rvelx=0.0
      m=0
      l=-1
      n=0
c      do 81 i=2*(nrowt-nrow),net,2*(nrowt-1)
c        m=m+1
c        if(m.lt.nsing) goto 81
c        l=l+2
c        if(l.gt.2*nsthick) l=2*nsthick
c        do j=1,l
c          vx(i-j+1)=rlvelx*cos(cangle(m))/3.15578e13
c          vz(i-j+1)=rlvelx*sin(cangle(m))/3.15578e13
c      	end do
c   81 continue
c      do i=1,ncol
c        vman(i)=vx((i-1)*2*(nrowt-1)+2*(nrowt-nrow))
c        if(i.ge.2) then
c          rwidth(i-1)=(rlthick(i-1)+rlthick(i))/2
c          fluxman(i-1)=vman(i-1)*rwidth(i-1)
c        endif
c        if(i.eq.ncol) fluxman(i)=fluxman(i-1)*(rlthick(i)/rlthick(i-1))
c      end do

c Remesh the Slab keeping flux through a vert colm constant
c	As written, it is assumed that the elements in the sub. lith. have
c	a uniform thickness in each colm.
c	  mesh the first row on lhs
      inodet=nrowt-nrow
      inodeb=nrowt-nrow-nsthick+1
      dx=sstop(1,2)-sstop(1,1)
      dy=sstop(2,2)-sstop(2,1)
      hyplen=(dx**2+dy**2)**0.5
c     new position for bot of lithos      
      ybase=-plthick*hyplen/dx+sstop(2,icol)
c     new heights of elms in lithos
      dyelm=(sstop(2,icol)-ybase)/dble(nsthick)
      icount=0
      ssbotold(1)=coordt(2,inodeb)
      do i=inodet,inodeb,-1
      	icount=icount+1
      	coordt(2,i)=coordt(2,i+1)-dyelm
        	nmesh=nmesh+1
      end do	
      ssbot(1)=coordt(2,inodeb)
c	mesh from lhs+1 to the spoint
      do icol=2,nsing
      	inodet=nrowt*icol-nrow
      	inodeb=inodet-nsthick+1
      	dx=sstop(1,icol+1)-sstop(1,icol-1)
      	dy=sstop(2,icol+1)-sstop(2,icol-1)
      	hyplen=(dx**2+dy**2)**0.5
      	ybase=-plthick*hyplen/dx+sstop(2,icol)
      	dyelm=(sstop(2,icol)-ybase)/dble(nsthick)
      	icount=0
      	ssbotold(icol)=coordt(2,inodeb)
      	do i=inodet,inodeb,-1
      		icount=icount+1
      		coordt(2,i)=coordt(2,i+1)-dyelm
        	nmesh=nmesh+1
      	end do	
      	ssbot(icol)=coordt(2,inodeb)
      end do
c mesh from spoint+1 to when base of slab hits model base
      jcount=0
      do icol=nsing+1,nsing+nrowt-nrow-nsthick
        jcount=jcount+1
      	inodet=nrowt*icol-nrow-jcount
      	inodeb=inodet-nsthick+1
      	dx=sstop(1,icol+1)-sstop(1,icol-1)
      	dy=sstop(2,icol+1)-sstop(2,icol-1)
      	hyplen=(dx**2+dy**2)**0.5
      	ybase=-plthick*hyplen/dx+sstop(2,icol)
      	ssbotold(icol)=coordt(2,inodeb)
      	dyelm=(sstop(2,icol)-ybase)/dble(nsthick)
      	icount=0
      	do i=inodet,inodeb,-1
      		icount=icount+1
      		coordt(2,i)=coordt(2,i+1)-dyelm
        	nmesh=nmesh+1
      	end do	
      	ssbot(icol)=coordt(2,inodeb)
      end do
c mesh from where bottom of slab hits model base to the end of the slab
c	In this region, just need to shift the nodes the same amount as the
c		slab top was shifted
      kcount=0
      do icol=nsing+1+nrowt-nrow-nsthick,nslabcol-1
      	jcount=jcount+1
      	kcount=kcount+1
      	inodet=nrowt*icol-nrow-jcount
      	inodeb=nrowt*icol-nrow-nsthick+1-jcount+kcount
      	icount=0
      	do i=inodeb,inodet
      		icount=icount+1
      		coordt(2,i)=coordt(2,i)+dysstop(icol)
        	nmesh=nmesh+1
      	end do	
      end do

c calculate the drop in the bottom of the slab
      do i=1,nsing+nrowt-nrow-nsthick
      	dyssbot(i)=ssbot(i)-ssbotold(i)
      end do	

c remesh the pro asthenosphere
c	from lhs to spoint
      do icol=1,nsing
      	inodet=nrowt*icol-nrow-nsthick
      	inodeb=nrowt*(icol-1)+1
      	do i=inodeb,inodet
      		coordt(2,i)=coordt(2,i)+dyssbot(icol)
        	nmesh=nmesh+1
      	end do
      end do	
c	from spoint+1 to where base of slab touches model base -1
      icount=0
      do icol=nsing+1,nsing+nrowt-nrow-nsthick-1
      	icount=icount+1
      	inodet=nrowt*icol-nrow-nsthick-icount
      	inodeb=nrowt*(icol-1)+1
      	do i=inodeb,inodet
      		coordt(2,i)=coordt(2,i)+dyssbot(icol)
        	nmesh=nmesh+1
      	end do	
      end do	

c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c compute velocity for thermal model
c CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c mech domain
      nlthick=nrowt-nrow
      tricount=0
      k=0
      l=0
      do i=2*nlthick,net,2*(nrowt-1) 
      	do j=1,2*(nrow-1),2
      		k=i+j
	  		tricount=(i-2*nlthick)/(2*(nrowt-1))*((nrow-1))+
     *		(j+1)/2
      		vx(k)=(velx(node(1,tricount))+velx(node(2,tricount))+
     *		velx(node(4,tricount)))/(3*3.15578e13)
      		vx(k+1)=(velx(node(4,tricount))+velx(node(2,tricount))+
     *		velx(node(3,tricount)))/(3*3.15578e13)
      		vz(k)=(vely(node(1,tricount))+vely(node(2,tricount))+
     *		vely(node(4,tricount)))/(3*3.15578e13)
      		vz(k+1)=(vely(node(4,tricount))+vely(node(2,tricount))+
     *		vely(node(3,tricount)))/(3*3.15578e13)
c  			Find Basal Velocity Values
      		if(j.eq.1) then
	  			l=(i-2*nlthick)/(2*(nrowt-1))+1
      			bvelm(l)=convel/3.15578e13
      		endif
      	end do
      end do

c	average thickness of the subducting slab
c     from lhs to when base of slab hits bottom of model -1
      do i=1,nsing+nrowt-nrow-nsthick-1
      	wwidth(i)=((sstop(2,i)-ssbot(i))+(sstop(2,i+1)-ssbot(i+1)))/2.0
      end do	
c     node at the point where slab base hits model bottom
      icol=nsing+nrowt-nrow-nsthick
      wwidth(icol)=((sstop(2,icol)-ssbot(icol))+(sstop(2,i+1)
     *-coordt(2,icol*nrowt+1)))/2.0
c     rest of slab
      do i=nsing+nrowt-nrow-nsthick+1,nslabcol-1
      	wwidth(i)=((sstop(2,i)-coordt(2,(i-1)*nrowt+1))
     *	+(sstop(2,i+1)-coordt(2,nrowt*i+1)))/2.0
      end do	

c subducting slab domain
      k=0
      l=0
      n=0
      m=0
      ll=0
      nc=0
      mcount=0
      nslabhit=nsing+nrowt-nrow-nsthick
      neslength=2*((nslength/nrowt)-1)*(nrowt-1)
      do i=2*(nrowt-nrow-nsthick),neslength,2*(nrowt-1)
        nc=nc+1
        mcount=mcount+1
        if(nc.ge.nsing) then
        	if(nc.lt.nslabhit) then
          		m=1
          		n=n+1
        	endif
        endif
        if(nc.ge.nsing) mcount=nsing-1
        if(nc.gt.nslabhit) then
        	l=l+1
        endif
        if(nc.ge.nslabhit) then
        	wwidth(nc)=wwidth(nc-1)
        	ll=1
        	m=0
        endif
c
c 		commneted out since not set up for ablation        
c
        
c        if(nc.gt.1) then
c        	fac=(flux(nc-1)/(wwidth(nc)*cos(-tangle(nc))*
c     *  	bvelm(mcount)))
c        	bvelm(mcount)=bvelm(mcount)*fac
c        endif

        do j=1,2*(nsthick-l)-ll
        	k=i+j+m-2*n
        	vx(k)=bvelm(mcount)*cos(-tangle(nc))
        	vz(k)=bvelm(mcount)*sin(-tangle(nc))
c
c			for ablation. commented out since not worried about 
c				ablation, yet. (5-21-03)
c
c        	if(j.eq.2*(nsthick-l)-ll) then
c        		flux(nc)=vx(k)*wwidth(nc)
c        		if(nc.eq.1) then
c        			snewthick(1)=(flux(1)+abs(fluxman(1)))
c     *				/(vx(2*(nrowt-nrow))*nsthick)
c        		 endif
c          		snewthick(nc+1)=((((flux(nc)+abs(fluxman(nc)))
c     *			/(vx(k)*nsthick))-snewthick(nc))*2)+snewthick(nc)
c        	endif
      	end do
      end do

c set all other velocities to zero (only need to do on 1st tst???)
      if(itst.eq.1) then
      	do i=1,net
c 		retro lith
  	      if(nodet(i,5).eq.3) then
  	        vx(i)=0.0
  	        vz(i)=0.0
  	      endif
c 		pro asthen
  	      if(nodet(i,5).eq.4) then
  	        vx(i)=0.0
  	        vz(i)=0.0
  	      endif
c 		retro asthen
  	      if(nodet(i,5).eq.5) then
  	        vx(i)=0.0
  	        vz(i)=0.0
  	      endif        
  	    end do
      endif


c  Redue slab thickness due to Ablative Subduction
c	for ablation. commented out since not worried about ablation yet 5-21-03
c      k=0
c      m=0
c      n=0
c      l=0
c      do 96 i=(nrowt-nrow+1),nslength,nrowt
c        m=m+1
c        if(m.gt.nsing) k=k+1
c        zadd(m)=abs(snewthick(m))
c        do 97 j=1,nsthick
c          if(j.gt.nsthick-n) then
c            goto 98
c          endif
c          l=i-j-k
c          coordt(2,l)=coordt(2,l+1)-(zadd(m))
c          if(l.le.((nrowt)*(m-1)+1)) then
c            n=n+1
c          endif
c   98 continue
c   97 continue
c   96 continue
c  Re-Calculate Slab bottom
c      m=0
c      k=0
c      do 44 i=(nrowt-nrow-nsthick+1),nslength,nrowt
c        m=m+1
c        if(m.ge.nsing+1)k=k+1
c        ssbot(m)=coordt(2,i-k)
c        if(i-k.eq.(nrowt*(m-1)+1)) goto 43
c   44 continue
c   43 continue
c  Remesh Asthenoshere below slab, after ablation adjustment
c      botdep=coordt(2,((ncount-1)*nrowt+1))
c      k=0
c      l=0
c      do i=1,ncount
c        if(i.gt.nsing) then
c           k=k+1
c        endif
c        nsbot=nrowt-nrow-nsthick-k+1
c        if(nsbot.le.0) nsbot=1
c        zup=coordt(2,nrowt*i)-ssbot(i)
c        zdown=coordt(2,nrowt*i)+abs(botdep)-zup
c        if(nsbot.eq.1) goto 99 
c        zspace=zdown/(nsbot-1)
c        do j=1,nsbot-1
c          l=(i-1)*nrowt+nsbot-j
c          coordt(2,l)=coordt(2,l+1)-zspace
c      	end do
c   99 continue
c      end do
      deallocate(cdown,cpresent,sstop,ccbot,rlthick,tangle,cangle,
     *vman,fluxman,rwidth,ssbot,ssbotold,sstopold,dysstop,dyssbot,
     *wwidth,zadd,bvelm,flux,snewthick)

      return
      end
c**********************************************************************
c*                                                                    *
c*    routine to calculate erosion flux as function of slope
c*                                                                    *
c**********************************************************************

c modified so that there is no erosion, only filling of basins below sealev
      subroutine erosion(nn,nrow,ninc,ncol,delt,itst,erosl,erosr,
     *peros,rpow,sealev,w_depth,isedl,isedr,ibasflg,intmrkb,
     *nbastrk,nbastary,nbastind,ninbas,ipkfill,ibasfill,
     *sedmax)

      use dyn_arrays_mech
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      real(kind=8),allocatable::xsur_old(:),peakloc(:),rst(:),
     *temp_array(:,:),vsur_old(:)
      integer,allocatable::locmax(:),locmin(:),itemp_array(:,:),
     *itemp_array1d(:)
      real*8 rightdrop,leftdrop

c allocate arrays
      allocate(xsur_old(ncol),peakloc(ncol),rst(ncol),locmax(ncol),
     *locmin(ncol),vsur_old(ncol))

      xsur_old=0.0
      peakloc=0.0
      rst=0.0
      locmax=0.0
      locmin=0.0
      istart=0
      iend=0
      do i=1,ncol
      	vsur_old(i)=vsur(2,i)
      end do	

c     indexes for type of fill
      ihbnd=0
      ihbas=0
      ihpk=0
      ihsea=0

c find peak in valley surface above sea level
      hpeak=0.
      ipeak=0
      do i=1,ncol
      	if(vsur(2,i).gt.hpeak)then
      		hpeak=vsur(2,i)
      		ipeak=i
      	endif
      end do
c      print*,'Peak',ipeak,vsur(1,ipeak)

c initialize arrays
      do i=1,ncol
      	veros(1,i)=0.0
      	veros(2,i)=0.0
      	vdiff(1,i)=0.0
      	vdiff(2,i)=0.0
      	rdiff(1,i)=0.0
      	rdiff(2,i)=0.0
      	xsur_old(i)=xsur(2,i)
      end do

c #####################################
c sedimentation 
c ##############
c bounds for filling

      ilbound=isedl
      irbound=isedr
      if(ipkfill.eq.1) then
c	Fill between peaks
  150 	continue
c     	find local peaks  
      	icount=0
      	do k=ilbound+1,irbound-1
      	  slopel=vsur(2,k)-vsur(2,k-1)
      	  sloper=vsur(2,k+1)-vsur(2,k)
      	  if(slopel.gt.0.and.sloper.lt.0) then
      	    icount=icount+1
      	    locmax(icount)=k
      	    peakloc(icount)=vsur(2,k)
      	  endif
      	end do
c 		if more than one peak, fill in basins
      	if(icount.gt.1) then
      		do j=1,icount-1
      			do i=locmax(j)+1,locmax(j+1)-1
      				if(peakloc(j).le.peakloc(j+1)) then
      					if(vsur(2,i).lt.peakloc(j)) then
      						nbasinfill(i)=nbasinfill(i)+1
      						basinfill(i)=basinfill(i)+(peakloc(j)
     *							-vsur(2,i))
      						vsur(2,i)=peakloc(j)
      					endif	
      				else
      					if(vsur(2,i).lt.peakloc(j+1)) then
      						nbasinfill(i)=nbasinfill(i)+1
      						basinfill(i)=basinfill(i)+(peakloc(j)
     *							-vsur(2,i))
      						vsur(2,i)=peakloc(j+1)
      					endif	
      				endif
       			end do 
      		end do
c    		loop until only one peak
      		goto 150
      	endif	
c	Fill in bounding basins
      	if(ibasfill.eq.1) then
      		do k=ilbound+1,irbound-1
        		slopel1=vsur(2,k)-vsur(2,k-1)
        		sloper1=vsur(2,k+1)-vsur(2,k)
        		if(slopel1.lt.0.0.and.sloper1.gt.0.0) then
c     				is the basin to the L or R of peak, then fill in basin
c					to the max fill height,height at sed bounds or height
c					of peak, which ever is less
      				if(vsur(1,k).lt.vsur(1,locmax(1))) then
      					hbnd=vsur(2,ilbound)
      					hbas=vsur(2,k)+sedmax
      					hpk=vsur(2,locmax(1))
      					if(hbnd.le.hbas.and.hbnd.le.hpk.
     *						and.hbnd.le.sealev) then
      						hfill=hbnd
      						ihbnd=ihbnd+1
      					elseif(hbas.le.hpk.and.hbas.le.hbnd.
     *						and.hbas.le.sealev) then
      						hfill=hbas
      						ihbas=ihbas+1
      					elseif(hpk.le.hbas.and.hpk.le.hbnd.
     *						and.hpk.le.sealev) then
      						hfill=hpk
      						ihpk=ihpk+1
      					elseif(sealev.le.hpk.and.sealev.le.hbas.
     *						and.sealev.le.hbnd) then
      						hfill=sealev
      						ihsea=ihsea+1
      					endif	
c      					print*,'fill:(itst,bnd,bas,pk,sea)\nPfill',
c     *					itst,ihbnd,ihbas,ihpk,ihsea
      					do ii=ilbound,locmax(1)
      						if(vsur(2,ii).lt.hfill) vsur(2,ii)=hfill
      					end do	
      				elseif(vsur(1,k).gt.vsur(1,locmax(1))) then
      					hbnd=vsur(2,irbound)
      					hbas=vsur(2,k)+sedmax
      					hpk=vsur(2,locmax(1))
      					if(hbnd.le.hbas.and.hbnd.le.hpk.
     *						and.hbnd.le.sealev) then
      						hfill=hbnd
      						ihbnd=ihbnd+1
      					elseif(hbas.le.hpk.and.hbas.le.hbnd.
     *						and.hbas.le.sealev) then
      						hfill=hbas
      						ihbas=ihbas+1
      					elseif(hpk.le.hbas.and.hpk.le.hbnd.
     *						and.hpk.le.sealev) then
      						hfill=hpk
      						ihpk=ihpk+1
      					elseif(sealev.le.hpk.and.sealev.le.hbas.
     *						and.sealev.le.hbnd) then
      						hfill=sealev
      						ihsea=ihsea+1
      					endif	
c      					print*,'fill:(itst,bnd,bas,pk,sea)\nRfill',
c     *					itst,ihbnd,ihbas,ihpk,ihsea
      					do ii=locmax(1),irbound
      						if(vsur(2,ii).lt.hfill) vsur(2,ii)=hfill
      					end do	
      				endif	
      			endif
      		end do
      	endif	

c 	check for basin fill tracking
      	if(ibasflg.eq.1) then
c     		check if a tracking timestep      	
      		if(itst.eq.1.or.mod(itst,intmrkb).eq.0) then
c				check that index array is big enough      		
      			if(nbastrk+ncol.gt.nbastind) then
      				print*,'Resize Basin Index Array'
      				allocate(itemp_array(2,nbastind))
      				itemp_array=ibastrk
      				deallocate(ibastrk)
      				allocate(ibastrk(2,(nbastrk+ncol)*2))
      				ibastrk=0
      				do i=1,2
      					do j=1,nbastind
      						ibastrk(i,j)=itemp_array(i,j)
      					end do
      				end do	
      				deallocate(itemp_array)
      				nbastind=(nbastrk+ncol)*2
      			endif	
c     			check that basin array is big enough      		
     			if(ninbas+ncol.gt.nbastary) then
c     				basin track array     			
     				print*,'Resize Basin Tracking Array'
     				allocate(temp_array(4,nbastary))
      				temp_array=bastrk
      				deallocate(bastrk)
      				allocate(bastrk(4,(ninbas+ncol)*2))
      				bastrk=0.0
      				do i=1,4
      					do j=1,nbastary
      						bastrk(i,j)=temp_array(i,j)
      					end do
      				end do	
      				deallocate(temp_array)
c     				eulerian elements for basin track array
      				allocate(itemp_array1d(nbastary))
      				itemp_array1d=ieletpb
      				deallocate(ieletpb)
      				allocate(ieletpb((ninbas+ncol)*2))
      				ieletpb=100
      				do j=1,nbastary
      					ieletpb(j)=itemp_array1d(j)
      				end do
      				deallocate(itemp_array1d)
      				nbastary=(ninbas+ncol)*2
      			endif
c     			store the basin surface location and the position in array 
      			istart=2
  398 			continue      		
      			do i=istart,ncol
      				if(vsur(2,i).gt.vsur_old(i)) then
       					nbastrk=nbastrk+1
      					ninbas=ninbas+1
       					ibastrk(1,nbastrk)=ninbas
       					bastrk(1,ninbas)=vsur(1,i-1)
       					bastrk(2,ninbas)=vsur(2,i-1)
       					do j=i,ncol
       						ninbas=ninbas+1
       						bastrk(1,ninbas)=vsur(1,j)
       						bastrk(2,ninbas)=vsur(2,j)
       						ibastrk(2,nbastrk)=ninbas
       						if(vsur(2,j).eq.vsur_old(j)) exit
       					end do
       					istart=j+1
       					goto 398
       				end if
       			end do	
      		endif
      	endif	
      endif	

c ##########
c End of sedimentation
c ###################################



c fluvial erosion from R basin to peak
      veros(1,ilbound)=0.0
      do 230 i=ilbound+1,ipeak
      delx=vsur(1,i)-vsur(1,i-1)
      delh=vsur(2,i-1)-vsur(2,i)
      veros(1,i)=erosl*(vsur(1,ipeak)-vsur(1,i))*
     *(delh+delt*veros(1,i-1))
     */(delx+erosl*(vsur(1,ipeak)-vsur(1,i))*delt)
  230 continue

c fluvial erosion from L basin to peak
c integrate from right basin to peak for erosion
c      implicit formula
      veros(1,irbound)=0.0
      do 240 ii=ipeak+1,irbound
      i=irbound-ii+ipeak
      delx=vsur(1,i+1)-vsur(1,i)
      delh=vsur(2,i+1)-vsur(2,i)
      veros(1,i)=erosr*(vsur(1,i)-vsur(1,ipeak))*
     *(delh+delt*veros(1,i+1))
     */(delx+erosr*(vsur(1,i)-vsur(1,ipeak))*delt)
  240 continue

c John's erosion at peak nodes
      leftdrop=-1.0*(peros*(vsur(2,ipeak)-vsur(2,ipeak-1)))**rpow
      rightdrop=-1.0*(peros*(vsur(2,ipeak)-vsur(2,ipeak+1)))**rpow
      veros(1,ipeak)=min(rightdrop,leftdrop)

C CC  bounds of i have changed to affect only the surface between the two 
C CC  peak bounding basins
      do 20 i=ilbound,irbound
      vtemp=vsur(2,i)
      vsur(2,i)=vsur(2,i)+delt*veros(1,i)
   20 continue

c jon, august. calculate range front seperately.
      do i=1,ncol
		rst(i)=rsur(2,i)
	  end do	

      do i = ilbound,irbound
c commented out on 5-17, I think this catch may allow the rsur to be
c		less than the vsur
        if(rsur(2,i)-vsur(2,i).lt.0.0) goto 353
      	rtemp = rsur(2,i) - delt*(peros*(rsur(2,i)-vsur(2,i)))**rpow
      	if(rtemp.lt.vsur(2,i)) then
      		veros(2,i)=(rsur(2,i)-vsur(2,i))/delt
	    	rsur(2,i) =vsur(2,i)
	    else 	
			rsur(2,i)=rsur(2,i) - delt*(peros*(rsur(2,i)-
     *				vsur(2,i)))**rpow
      		veros(2,i)=peros*(rsur(2,i)-vsur(2,i))
      	endif
  353 continue      	
      end do
      rsur(2,ipeak) = vsur(2,ipeak)
      veros(2,ipeak)=0.0
      do  i=ilbound,irbound
      if(rsur(2,i).lt.vsur(2,i)) then
      	rsur(2,i)=vsur(2,i)
      endif
      xsur(2,i)=(vsur(2,i)+rsur(2,i))/2.0
      vdiff(2,i)=xsur(2,i)-vsur(2,i)
      rdiff(2,i)=xsur(2,i)-rsur(2,i)
      rdiff(1,i)=xsur(1,i)
      vdiff(1,i)=xsur(1,i)
      end do

      if(ilbound.gt.1) then
      	do i=1,ilbound-1
      		vsur(2,i)=xsur(2,i)
      		rsur(2,i)=xsur(2,i)
      		rst(i)=xsur(2,i)
      		rdiff(2,i)=0.0
      		vdiff(2,i)=0.0
      	end do
      end if	
      if(irbound.lt.ncol) then
      	do i=irbound+1,ncol
      		rst(i)=xsur(2,i)
      		vsur(2,i)=xsur(2,i)
      		rsur(2,i)=xsur(2,i)
      		rdiff(2,i)=0.0
      		vdiff(2,i)=0.0
      	end do
      end if	
      do i=1,ncol
      	rdiff(1,i)=xsur(1,i)
      	vdiff(1,i)=xsur(1,i)
      end do

c redefine veros(2,i) to store the erosion rate as the change in xsur between 
c timesteps as opposed to the ridge erosion rate
      do i=1,ncol
      	veros(2,i)=(xsur_old(i)-xsur(2,i))/delt
      end do	

      deallocate(xsur_old,peakloc,rst,locmax,locmin,vsur_old)

      return
      end


c ############################################################
c Old erosion routine
c ############################################################
      subroutine erosion_old(nn,nrow,ninc,ncol,delt,itst,erosl,erosr,
     *peros,rpow,sealev,w_depth)

      use dyn_arrays_mech
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      real(kind=8),allocatable::xsur_old(:),peakloc(:),rst(:)
      integer,allocatable::locmax(:),locmin(:)
      real*8 rightdrop,leftdrop

c allocate arrays
      allocate(xsur_old(ncol),peakloc(ncol),rst(ncol),locmax(ncol),
     *locmin(ncol))
      xsur_old=0.0
      peakloc=0.0
      rst=0.0
      locmax=0.0
      locmin=0.0

c find peak in valley surface above sea level
      hpeak=0.
      ipeak=0
      do i=1,ncol
      	if(vsur(2,i).gt.hpeak)then
      		hpeak=vsur(2,i)
      		ipeak=i
      	endif
      end do
      if(hpeak.le.sealev) goto 333
      print*,'Peak',ipeak,vsur(1,ipeak)

c initialize arrays
      do i=1,ncol
      	veros(1,i)=0.0
      	veros(2,i)=0.0
      	vdiff(1,i)=0.0
      	vdiff(2,i)=0.0
      	rdiff(1,i)=0.0
      	rdiff(2,i)=0.0
      	xsur_old(i)=xsur(2,i)
      end do

c find the base level to the L and R of peak
      ilbound=-100
      irbound=-100
      do i=ipeak-1,1,-1
      	if(vsur(2,i).le.sealev) then
      		if(abs(vsur(2,i)-sealev).lt.abs(vsur(2,i+1)-sealev)) then
      			ilbound=i
      			exit
      		else
      			ilbound=i+1
      			if(ilbound.eq.ipeak) ilbound=ipeak-1
      			exit
      		endif
      	endif
      end do	
      do i=ipeak+1,ncol
      	if(vsur(2,i).le.sealev) then
      		if(abs(vsur(2,i)-sealev).lt.abs(vsur(2,i-1)-sealev)) then
      			irbound=i
      			exit
      		else
      			irbound=i-1
      			if(irbound.eq.ipeak) irbound=ipeak+1
      			exit
      		endif
      	endif
      end do	
      if(ilbound.eq.-100) then
      	print*,'@@@@ NOTE: L erosion bound not found'
      	print*,'     using model edge'
      	ilbound=1
      endif	
      if(irbound.eq.-100) then
      	print*,'@@@@ NOTE: R erosion bound not found'
      	print*,'     using model edge'
      	irbound=ncol
      endif
      print*,'#######  Erosion bounds:',ilbound,irbound
      print*,'       ',vsur(1,ilbound),vsur(1,irbound)

c fill in any basins between the base level on the L and R of the peak
  150 continue
c        finds local peaks  
      icount=0
      do k=ilbound+1,irbound-1
        slopel=vsur(2,k)-vsur(2,k-1)
        sloper=vsur(2,k+1)-vsur(2,k)
        if(slopel.gt.0.and.sloper.lt.0) then
          icount=icount+1
          locmax(icount)=k
          peakloc(icount)=vsur(2,k)
        endif
      end do

c if more than one peak, fill in basins
      if(icount.gt.1) then
      	do j=1,icount-1
      		do i=locmax(j)+1,locmax(j+1)-1
      			if(peakloc(j).le.peakloc(j+1)) then
      				if(vsur(2,i).lt.peakloc(j)) then
      					nbasinfill(i)=nbasinfill(i)+1
      					basinfill(i)=basinfill(i)+(peakloc(j)-vsur(2,i))
      					vsur(2,i)=peakloc(j)
      				endif	
      			else
      				if(vsur(2,i).lt.peakloc(j+1)) then
      					nbasinfill(i)=nbasinfill(i)+1
      					basinfill(i)=basinfill(i)+(peakloc(j)-vsur(2,i))
      					vsur(2,i)=peakloc(j+1)
      				endif	
      			endif
       		end do 
      	end do
c    	loop until only one peak
      	goto 150
      endif	




c fluvial erosion from R basin to peak
      veros(1,ilbound)=0.0
      do 230 i=ilbound+1,ipeak
      delx=vsur(1,i)-vsur(1,i-1)
      delh=vsur(2,i-1)-vsur(2,i)
      veros(1,i)=erosl*(vsur(1,ipeak)-vsur(1,i))*
     *(delh+delt*veros(1,i-1))
     */(delx+erosl*(vsur(1,ipeak)-vsur(1,i))*delt)
  230 continue

c fluvial erosion from L basin to peak
c integrate from right basin to peak for erosion
c      implicit formula
      veros(1,irbound)=0.0
      do 240 ii=ipeak+1,irbound
      i=irbound-ii+ipeak
      delx=vsur(1,i+1)-vsur(1,i)
      delh=vsur(2,i+1)-vsur(2,i)
      veros(1,i)=erosr*(vsur(1,i)-vsur(1,ipeak))*
     *(delh+delt*veros(1,i+1))
     */(delx+erosr*(vsur(1,i)-vsur(1,ipeak))*delt)
  240 continue

c John's erosion at peak nodes
      leftdrop=-1.0*(peros*(vsur(2,ipeak)-vsur(2,ipeak-1)))**rpow
      rightdrop=-1.0*(peros*(vsur(2,ipeak)-vsur(2,ipeak+1)))**rpow
      veros(1,ipeak)=min(rightdrop,leftdrop)

C CC  bounds of i have changed to affect only the surface between the two 
C CC  peak bounding basins 
      do 20 i=ilbound,irbound
      vtemp=vsur(2,i)
      vsur(2,i)=vsur(2,i)+delt*veros(1,i)
   20 continue

c jon, august. calculate range front seperately.
      do i=1,ncol
		rst(i)=rsur(2,i)
	  end do	

      do i = ilbound,irbound
c commented out on 5-17, I think this catch may allow the rsur to be
c		less than the vsur
        if(rsur(2,i)-vsur(2,i).lt.0.0) goto 353
      	rtemp = rsur(2,i) - delt*(peros*(rsur(2,i)-vsur(2,i)))**rpow
      	if(rtemp.lt.vsur(2,i)) then
      		veros(2,i)=(rsur(2,i)-vsur(2,i))/delt
	    	rsur(2,i) =vsur(2,i)
	    else 	
			rsur(2,i)=rsur(2,i) - delt*(peros*(rsur(2,i)-
     *				vsur(2,i)))**rpow
      		veros(2,i)=peros*(rsur(2,i)-vsur(2,i))
      	endif
  353 continue      	
      end do
      rsur(2,ipeak) = vsur(2,ipeak)
      veros(2,ipeak)=0.0
      do  i=ilbound,irbound
      if(rsur(2,i).lt.vsur(2,i)) then
      	rsur(2,i)=vsur(2,i)
      endif
      xsur(2,i)=(vsur(2,i)+rsur(2,i))/2.0
      vdiff(2,i)=xsur(2,i)-vsur(2,i)
      rdiff(2,i)=xsur(2,i)-rsur(2,i)
      rdiff(1,i)=xsur(1,i)
      vdiff(1,i)=xsur(1,i)
      end do

      if(ilbound.gt.1) then
      	do i=1,ilbound-1
      		vsur(2,i)=xsur(2,i)
      		rsur(2,i)=xsur(2,i)
      		rst(i)=xsur(2,i)
      		rdiff(2,i)=0.0
      		vdiff(2,i)=0.0
      	end do
      end if	
      if(irbound.lt.ncol) then
      	do i=irbound+1,ncol
      		rst(i)=xsur(2,i)
      		vsur(2,i)=xsur(2,i)
      		rsur(2,i)=xsur(2,i)
      		rdiff(2,i)=0.0
      		vdiff(2,i)=0.0
      	end do
      end if	
      do i=1,ncol
      	rdiff(1,i)=xsur(1,i)
      	vdiff(1,i)=xsur(1,i)
      end do

c redefine veros(2,i) to store the erosion rate as the change in xsur between 
c timesteps as opposed to the ridge erosion rate
      do i=1,ncol
      	veros(2,i)=(xsur_old(i)-xsur(2,i))/delt
      end do	

  333 continue
      deallocate(xsur_old,peakloc,rst,locmax,locmin)

      return
      end

c *****************************************************************
      subroutine elt03n(inopredv,d,ul,xl,ix,s,p,ndf,ndm,nst,isw,dt
     *,nen,n,nel,viscos,press,plasfl,sigav,epsav,maxnel,kstep)
c
c
c
c
c this is a modified version of elmt03 that should
c allow for     1.    4/1 or 9/3 v-p or u-p computations
c               2.    nonlinear viscous models
c                     both in shear and isotropic part
c
c
c
c
c
      use dyn_arrays_mech
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      parameter (nstbis=21)
      dimension viscos(maxnel,9),press(maxnel,3)
     1       ,plasfl(maxnel,9)
      dimension sigav(maxnel,4),epsav(maxnel,4)
      dimension  d(*),ul(ndf,*),xl(ndm,*),ix(*),s(nst,*)
     1  ,shp(3,9),sg(9),tg(9),wg(9),sig(6),eps(3),wd(2)
     2  ,p(ndf,*),v(2),dv(2,2),shpp(3),indx(nstbis)
c
      dimension cmpp(6),cdpu(54),fp(3),ptot(3)
      dimension devstre(4),epsdev(4),stressl(4)
c
c      data wd/4hress,4hrain/
c     write(*,*)inopredv,(d(j),j=1,10),(ul(1,j),j=1,4),(ul(2,j),j=1,4),
c    *(xl(1,j),j=1,4),ix(1),s(1,1),p(1,1),p(2,1),
c    *ndf,ndm,nst,isw,dt
c    *,nen,n,nel,(viscos(1,j),j=1,4),(bulkmod(1,j),j=1,4),
c    *press(1,1),(plasfl(1,j),j=1,4),
c    *(sigav(1,j),j=1,4),(epsav(1,j),j=1,4),maxnel,kstep
c     write(*,*)inopredv,ndf,ndm,nst,isw
c     write(*,*)nen,n,nel,kstep
c     write(*,*)dt
c     write(*,*)(d(j),j=1,10)
c     write(*,*)(xl(1,j),j=1,4)
c     write(*,*)(xl(2,j),j=1,4)
c     write(*,*)(viscos(1,j),j=1,4)
c     write(*,*)(bulkmod(1,j),j=1,4)
c     write(*,*)(press(1,j),j=1,1)
c     write(*,*)(plasfl(1,j),j=1,4)
      nelp=1
      if(isw.eq.1) goto 1
      if(isw.eq.2) goto 2
      if(isw.eq.3) goto 3
      if(isw.eq.4) goto 4
      if(isw.eq.5) goto 5
c      goto (1,2,3,4,5,3),isw
c notice mode 1 is used in the input section (and can be reentered
c by macro: mesh although without ''change of b.c. '' if we quote
c mozart 's author. (mister taylor).but this should not be capital
c in particular when the topology is constant it poses no problem.
c the problem is more a interpolation/extrapolation problem
c from grid to grid than a profile recomputation pb.)
c
c
c
c     notice a storage (see d(10,)) of 10 is allowed for ech material
c     set.
c
    1 read(1,1000)d(1),d(2),d(3),l,d(5),d(6),d(7)
c     if the run is a restart the read phase should be in
c     restart mode.
      d(4)=l
      lint=0
      return
    2 return
    3 l=d(4)
      nstu=ndf*nen
      if(isw.eq.6)then
c it is very important to notice here that 4 modes can exist:
c a neutral mode: 1.no reaction computed  ,return.(ex a simple fluid)
c                 2.reaction computed but not fed to rhs
c                 3.standard mode compute and feed reactions.
c                 4.half standard mode :
c
c
c                   compute isotropic reactions say isw=7
c                           deviatoric reactions    isw=8
c                dissipated deviatoric rections     isw=9
c                dissipated isotropic  rections     isw=10 (if bulk visc
c                or plast)
c                etc!(any computation that would differ from
c                standard total div(sig(eps(u))
c                useful!!!!!! particular case are the following:
c                4.1.
c we lock here in mode 1 but macro could call other modes as well
                  goto 70
                  endif
c get pressure back from saved matrices and velocities
      if(inopredv.eq.0.and.isw.eq.3)then
cc    nelps=nelp*nelp
cc    read(mswap1,*)(cmpp(i),i=1,nelps)
cc    nelpu=nelp*nstu
cc    read(mswap1,*)(cdpu(i),i=1,nelpu)
c solve cmpp*p+cdpu*u=cmpp*p0
cc    call elimp(ndf,fp,cmpp,cdpu,ptot,nelp,nstu,ul)
c so dp is in pinc
      ptot(1)=press(n,1)
      ptot(2)=press(n,2)
      ptot(3)=press(n,3)
                      endif
c
      call pgauss(l,lint,sg,tg,wg)
      volt=0.
      ptemp=0.
      sigav(n,1)=0.
      sigav(n,2)=0.
      sigav(n,3)=0.
      sigav(n,4)=0.
      epsav(n,1)=0.
      epsav(n,2)=0.
      epsav(n,3)=0
      epsav(n,4)=0.
      do 65 l=1,lint
      call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
      shpp(1)=1.
      shpp(2)=sg(l)
      shpp(3)=tg(l)
c compute v at l
      do 38 i=1,2
      v(i)=0.
      do 31 k=1,nel
      v(i)=v(i)+shp(3,k)*ul(i,k)
   31 continue
c compute gradv at l
      do 37 j=1,2
      ddv=0.0
      do 32 k=1,nel
      ddv=ddv+shp(j,k)*ul(i,k)
   32 continue
      dv(i,j)=ddv
   37 continue
   38 continue
c from dv every strain or spin rate ...
      epstra=(dv(1,1)+dv(2,2))/3.
c convention 1=xx 2=yy 3=xy(not 2*xy) 4=zz=out of plane
      epsdev(1)=dv(1,1)-epstra
      epsdev(2)=dv(2,2)-epstra
      epsdev(3)=(dv(1,2)+dv(2,1))/2.
c because this is the plane strain elmt
      epsdev(4)=0.
      epsav(n,1)=epsav(n,1)+dv(1,1)
      epsav(n,2)=epsav(n,2)+dv(2,2)
      epsav(n,3)=epsav(n,3)+0.5*(dv(1,2)+dv(2,1))
      epsav(n,4)=epsav(n,3)+0.5*(dv(1,2)-dv(2,1))
c     linear case   or no predictor
      if(inopredv.eq.1)then
c or restart
      if(kstep.eq.1)then
      viscos(n,l)=d(1)
      bulkmod(n,l)=d(2)
                    endif
      xvol=1.
      xlam=dt*bulkmod(n,l)
      xcom=1.0/xlam
      xmu=viscos(n,l)
      xrho=d(3)
                      else
c     nonlinear case = nonlinear iteration technique .
c  in the general case use this to compute stress predictor and next stif
      xvol=1.
      xlam=dt*bulkmod(n,l)
      xcom=1.0/xlam
      xrho=d(3)
      xmu=viscos(n,l)
c strain stress law here viscous '!isotropic!'
      devstre(1)=2.*xmu*epsdev(1)
      devstre(2)=2.*xmu*epsdev(2)
      devstre(3)=2.*xmu*epsdev(3)
      devstre(4)=2.*xmu*epsdev(4)
c invariant
      rj2d=(devstre(1)*devstre(1)+devstre(2)*devstre(2))/2.0+
     *devstre(3)*devstre(3)
      rj2d=sqrt(rj2d)
c refind pressure at gauss point level
      presl=0.0
      do 1155 i=1,nelp
c dont forget to update press(n,i)
      presl=presl+(ptot(i))*shpp(i)
 1155 continue
      stressl(1)=-presl+devstre(1)
      stressl(2)=-presl+devstre(2)
      stressl(3)=devstre(3)
      stressl(4)=-presl+devstre(4)
      ptemp=ptemp+presl
      sigav(n,1)=sigav(n,1)+stressl(1)
      sigav(n,2)=sigav(n,2)+stressl(2)
      sigav(n,3)=sigav(n,3)+stressl(3)
      sigav(n,4)=sigav(n,4)+stressl(4)
c compute state variable control
      ivmises=0
      if(d(5).lt.0.)then
      ivmises=1
      sigy=-d(5)
                    else
      d5=3.14159*d(5)/180.
      cosphi=dcos(d5)
      sinphi=dsin(d5)
      coh2=d(6)
      if(presl.lt.0.0)then
      sigy=coh2*cosphi
      else
      sigy=presl*sinphi+coh2*cosphi
      endif
                    endif
      if(sigy.gt.d(7))sigy=d(7)
      if(sigy.lt.0)write(*,*)'pos 1: sigy < 0 elt n= ',n
c
c     radial return
c
      radret=rj2d/sigy
         if(plasfl(n,l).gt.0.or.radret.gt.1.)then
c plastic flow
         plasfl(n,l)=1
c notice that the following computation is redundant.
         rj2de=(epsdev(1)**2+epsdev(2)**2)/2.0+epsdev(3)**2
         rj2de=dsqrt(rj2de)
         xmu=sigy/(2.0*rj2de)
         if(xmu.gt.d(1))then
                        xmu=d(1)
                        plasfl(n,l)=0
                             endif
                                             endif
c     update nonlinear rheology
c here in general the whole rheology is reparametrized
             viscos(n,l)=xmu
                      endif
c

      xvol=xvol*xsj*wg(l)
      xlam=xlam*xsj*wg(l)
      xcom=xcom*xsj*wg(l)
      xmu=xmu*xsj*wg(l)
      xrho=xrho*xsj*wg(l)
      volt=volt+xvol
c     write(2,*)'end control                        '
c
c
c     isotropic operator     : spp
c     (dev-is  coupling)
c
      do 400 lp=1,nelp
      do 401 mp=1,nelp
      s(nstu+lp,nstu+mp)=s(nstu+lp,nstu+mp)+
     1xcom*shpp(lp)*shpp(mp)
  401 continue
  400 continue
c
      if(isw.eq.6)goto 60
      k1=1
c nel = 4, so not worth parallelizing?
      do 34 k=1,nel
c add this line
	k1=1+(k-1)*ndf
      a1=xmu*shp(1,k)
      a2=xmu*shp(2,k)
      a3=xrho*(dv(1,1)*shp(3,k)+v(1)*shp(1,k)+v(2)*shp(2,k))
      a4=xrho*(dv(2,2)*shp(3,k)+v(1)*shp(1,k)+v(2)*shp(2,k))
      a5=xrho*dv(1,2)*shp(3,k)
      a6=xrho*dv(2,1)*shp(3,k)
c eliminate deviatoric part
c     b1=xlam*shp(1,k)
c     b2=xlam*shp(2,k)
      b1=0.
      b2=0.
      j1=1
      do 33 j=1,nel
c add this line
	j1=1+(j-1)*ndf
c
c
c     deviatoric operator    : suu
c     (dev-dev coupling)
c
c xj xk
      s(j1,k1)=s(j1,k1)+shp(1,j)*a1+shp(2,j)*a2
      s(j1,k1)=s(j1,k1)+(shp(1,j)*a1)/3.0
c xj yk
      s(j1,k1+1)=s(j1,k1+1)+0.
c     s(j1,k1+1)=s(j1,k1+1)+a1*shp(2,j)/3.0
      s(j1,k1+1)=s(j1,k1+1)-2.*a2*shp(1,j)/3.0+a1*shp(2,j)
c yj xk
      s(j1+1,k1)=s(j1+1,k1)+0.
c     s(j1+1,k1)=s(j1+1,k1)+a2*shp(1,j)/3.0
      s(j1+1,k1)=s(j1+1,k1)-2.*a1*shp(2,j)/3.0+a2*shp(1,j)
c yj yk
      s(j1+1,k1+1)=s(j1+1,k1+1)+shp(1,j)*a1+shp(2,j)*a2
      s(j1+1,k1+1)=s(j1+1,k1+1)+(shp(2,j)*a2)/3.0
c this if statement breaks the elegance of the code helas!
c     write(2,*)'end suu                            '
      if(k.eq.1)then
c
c
c     iso-dev   operator     : sup
c     (dev-is  coupling)
c
      do 333 mp=1,nelp
      s(nstu+mp,j1)=s(nstu+mp,j1)+xvol*shpp(mp)*shp(1,j)
      s(nstu+mp,j1+1)=s(nstu+mp,j1+1)+xvol*shpp(mp)*shp(2,j)
      s(j1,nstu+mp)=s(nstu+mp,j1)
      s(j1+1,nstu+mp)=s(nstu+mp,j1+1)
  333 continue
      endif
c      j1=j1+ndf
   33 continue
c      k1=k1+ndf
   34 continue
c
c     solve iso-dev coupling at the element level :
c     elimination of internal dofs .here pressure.
c
c     if u-u convective term is not zero s is not symmetric
c     if u-p convective term is not zero s is not symmetric
c        u-p convective term arises from stress rate computations
      goto 65
c  force-computation
   60 continue
c it is very important to notice here that 3 modes can exist:
c a neutral mode: 1.no reaction computed  ,return.(ex a simple fluid)
c                 2.reaction computed but not fed to rhs
c                 3.standard mode compute and feed reactions.
c we lock here in mode 1 but macro could call other modes as well
      xdiv=(dv(1,1)+dv(2,2))*xlam
      do 67 k=1,nel
      do 64 j=1,2
      sum=xdiv*shp(j,k)
      do 63 i=1,2
      sum=sum+xmu*(dv(j,i)+dv(i,j))*shp(i,k)+
     1   xrho*v(i)*dv(j,i)*shp(3,k)
   63 continue
      p(j,k)=p(j,k)-sum
   64 continue
   67 continue
c     write(2,*)'end loop 65 l                      '
   65 continue
c     write(2,*)'loop 65 terminated'
c save cmpp and cdpu in file 3
      if(isw.eq.3.and.inopredv.eq.0)then
      if(nelp.eq.1)then
      cmpp(1)=s(nstu+1,nstu+1)
      iia=0
c!OMP parallel do private(j,iia)
      do 1111 j=1,nstu
c commented this line
c         iia=iia+1
c added this line
	 iia=j
         cdpu(iia)=s(nstu+1,j)
 1111 continue
c!OMP end parallel do
                   endif
      if(nelp.eq.3)then
      cmpp(1)=s(nstu+1,nstu+1)
      cmpp(2)=s(nstu+1,nstu+2)
      cmpp(3)=s(nstu+1,nstu+3)
      cmpp(4)=s(nstu+2,nstu+2)
      cmpp(5)=s(nstu+2,nstu+3)
      cmpp(6)=s(nstu+3,nstu+3)
      iia=0
c!OMP parallel do private(i,j,iia)
      do 1112 i=1,3
         do 1113 j=1,nstu
c commented this line
c         iia=iia+1
c added this line
	  iia=(i-1)*3+j
         cdpu(iia)=s(nstu+i,j)
 1113 continue
 1112 continue
c!OMP end parallel do
                   endif
cc    nelps=nelp*nelp
cc    write(mswap1,*)(cmpp(i),i=1,nelps)
cc    nelpu=nelp*nstu
cc    write(mswap1,*)(cdpu(i),i=1,nelpu)
c solve cmpp*p+cdpu*u=cmpp*p0
cc dont  forget to update fp in compressible case
c	write(6,*) 'calling elimp',ndf,nelp,nstu
      call elimp(ndf,fp,cmpp,cdpu,ptot,nelp,nstu,ul)
      press(n,1)=ptot(1)
      press(n,2)=ptot(2)
      press(n,3)=ptot(3)
                   endif
      if(nelp.eq.1)then
c     sigav(n,1)=(sigav(n,1)+ptemp-ptot(1))/lint
c     sigav(n,2)=(sigav(n,2)+ptemp-ptot(1))/lint
c     sigav(n,4)=(sigav(n,4)+ptemp-ptot(1))/lint
c     sigav(n,3)=sigav(n,3)/lint
      sigav(n,1)=(sigav(n,1)+ptemp)/lint
      sigav(n,2)=(sigav(n,2)+ptemp)/lint
      sigav(n,4)=(sigav(n,4)+ptemp)/lint
      sigav(n,3)=sigav(n,3)/lint
                   endif
      epsav(n,1)=epsav(n,1)/lint
      epsav(n,2)=epsav(n,2)/lint
      epsav(n,3)=epsav(n,3)/lint
      epsav(n,4)=epsav(n,4)/lint
c
c try hand elimination for n/1 elt
      if(nelp.eq.1)then
c     compress=volt/d(2)
c compressibility included in spp for possible nonlinear iterations
      do 5555 i=1,nstu
      do 5556 j=1,nstu
c     s(i,j)=s(i,j)+s(i,9)*s(9,j)/compress
      s(i,j)=s(i,j)+s(i,9)*s(9,j)/s(9,9)
 5556 continue
 5555 continue
                   endif
      if(nelp.eq.3)then
c       write(2,*)'before ludcmp'
      call ludcmp(s,nst,nst,indx(1),dperm,nelp)
c     write(2,*)'after  ludcmp'
                   endif
c     write(*,*)'stifness matrix'
c     write(*,*)s(1,1),s(1,2),s(1,3),s(1,4)
c     write(*,*)s(1,5),s(1,6),s(1,7),s(1,8)
c     write(*,*)s(2,2),s(2,3),s(2,4),s(2,5)
c     write(*,*)s(2,6),s(2,7),s(2,8)
c     write(*,*)s(3,3),s(3,4),s(3,5),s(3,6)
c     write(*,*)s(3,7),s(3,8)
c     write(*,*)s(4,4),s(4,5),s(4,6),s(4,7)
c     write(*,*)s(4,8)
c     write(*,*)s(5,5),s(5,6),s(5,7),s(5,8)
c     write(*,*)s(6,6),s(6,7),s(6,8)
c     write(*,*)s(7,7),s(7,8)
c     write(*,*)s(8,8)
c     stop
   70 continue
      return
c stress and vel gradient computation
    4 l=d(4)
c an selective integration rule could easily be inplemented on
c stif and stress computation:
c would essentially yield v-p: 2 on dev and 1 on trace term
c     if(isw.eq.4)l=d(6)
      call pgauss(l,lint,sg,tg,wg)
      do 600 l=1,lint
      call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
c here eps represent epsdot tensor
      do 410 i=1,3
      eps(i)=0.0
  410 continue
      xx=0.0
      yy=0.0
c added these 3 lines
	eps1=eps(1)
	eps2=eps(2)
	eps3=eps(3)
      do 420 j=1,nel
      xx=xx+shp(3,j)*xl(1,j)
      yy=yy+shp(3,j)*xl(2,j)
c of course ul is a velocity
c commented these three lines
c      eps(1)=eps(1)+shp(1,j)*ul(1,j)
c      eps(3)=eps(3)+shp(2,j)*ul(2,j)
c      eps(2)=eps(2)+shp(1,j)*ul(2,j)+shp(2,j)*ul(1,j)
c added these 3 lines
      eps1=eps1+shp(1,j)*ul(1,j)
      eps3=eps3+shp(2,j)*ul(2,j)
      eps2=eps2+shp(1,j)*ul(2,j)+shp(2,j)*ul(1,j)
  420 continue
      sig(1)=d(1)*eps(1)+d(2)*eps(3)
      sig(3)=d(2)*eps(1)+d(1)*eps(3)
      sig(2)=d(3)*eps(2)
c     if(isw.eq.6) go to 620
      call pstres(sig,sig(4),sig(5),sig(6))
      mct=mct-2
      if(mct.le.0) mct=50
      goto 600
c 620 dv=xsj*wg(l)
c     j1=1
c     do 610 j=1,nel
c     p(j1)=p(j1)-(shp(1,j)*sig(1)+shp(2,j)*sig(2))*dv
c     if(isw.eq.6)then
c                 endif
c     p(j1+1)=p(j1+1)-(shp(1,j)*sig(2)+shp(2,j)*sig(3))*dv
c 610 j1=j1+ndf
  600 continue
      return
c mass matrix
    5 l=d(4)
      call pgauss(l,lint,sg,tg,wg)
      do 503 l=1,lint
      call shape(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
c or any rho replacing d(3)!
      dvscal=wg(l)*xsj*d(3)
      j1=1
      do 500 j=1,nel
c added this line
	j1=(j-1)*ndf+1
      w11=shp(3,j)*dvscal
      k1=j1
      do 510 k=j,nel
      s(j1,k1)=s(j1,k1)+shp(3,k)*w11
      k1=k1+ndf
  510 continue
c commented this line
c      j1=j1+ndf
  500 continue
  503 continue
      nsl=nel*ndf
      do 521 j=1,nsl,ndf
      do 520 k=j,nsl,ndf
      s(j+1,k+1)=s(j,k)
      s(k,j)=s(j,k)
      s(k+1,j+1)=s(j,k)
  520 continue
  521 continue
c     write(*,*)'mass matrix'
c     write(*,*)s(1,1),s(2,2),s(3,3),s(4,4)
c     write(*,*)s(5,5),s(6,6),s(7,7),s(8,8)
c     write(*,*)s(9,9)
 1000 format(3e9.2,i2,3e9.2)
      return
      end
c *****************************************************************************
      subroutine ludcmp(s,nst,ndims,indx,dperm,nelim)
      implicit real*8 (a-h,o-z)
      parameter (maxnel=1,nstbis=21)
      dimension s(nst,*),indx(nstbis)
c     1,scalro(nstbis)
      parameter (tiny=0.)
c we are not aiming here to time efficiency
c flip s to p-u order
	write(6,*) 'nst:',nst
      do 5000 i=1,nst
      do 5001 j=1,nst
      s(i,j)=s(nst+1-i,nst+1-j)
 5001 continue
 5000 continue
      dperm=1.
c scalro not needed if s symmetric
c     do 12 i=1,nst
c           smax=0.
c           do 11 j=1,nst
c           if(abs(s(i,j).gt.smax)smax=abs(s(i,j))
c  11 continue
c     if(smax.eq.0.)then
c     write(*,*)'s is singular in lu routine'
c                   stop
c                   endif
c     scalro(i)=1./smax
c  12 continue
c
c
c     only eliminate internal dofs :1 to nelim.
c          for pressure elimination,nelim=nelp
      do 19 j=1,nelim
      if(j.gt.1)then
                do 14 i=1,j-1
                sum=s(i,j)
                if(i.gt.1)then
                do 13 k=1,i-1
                   sum=sum-s(i,k)*s(k,j)
   13 continue
                s(i,j)=sum
                          endif
   14 continue
                endif
      smax=0.
      do 16 i=j,nelim
         sum=s(i,j)
         if(j.gt.1)then
            do 15 k=1,j-1
               sum=sum-s(i,k)*s(k,j)
   15 continue
         s(i,j)=sum
                   endif
c no merit for pivoting if matrix s is symmetric.
c         dum=scalro(i)*abs(sum)
c         if(dum.ge.smax)then
c                        imax=i
c                        smax=dum
c                        endif
   16 continue
c      if(j.ne.imax)then
c         do 17 k=1,nelim
c            dum=s(imax,k)
c            s(imax,k)=s(j,k)
c            s(j,k)=dum
c   17 continue
c      dperm=-dperm
c      scalro(imax)=scalro(j)
c                   endif
      if(j.ne.nelim)then
              if(s(j,j).eq.0)s(j,j)=tiny
      if(s(j,j).eq.0.)then
              write(*,*)'0 pivot in internal dof elimination'
                      stop
                      endif
              dum=1./s(j,j)
              do 18 i=j+1,nelim
                    s(i,j)=s(i,j)*dum
   18 continue
                    endif
   19 continue
      if(s(nelim,nelim).eq.0.)s(nelim,nelim)=tiny
c flip s to u-p order
      do 6000 i=1,nst
      do 6001 j=1,nst
      s(i,j)=s(nst+1-i,nst+1-j)
 6001 continue
 6000 continue
      return
      end
c ######################################################################      
      subroutine elimp(ndf,fp,cmpp,cdpu,pinc,nelp,nstu,ul)
      implicit real*8 (a-h,o-z)
      dimension ul(ndf,*),pinc(3),cmpp(6),cdpu(54),fp(3),cmpp1(6)
      dimension rhs(3)
      loc=0
      iax=0
      ino=1
      do 10 i=1,nelp
      temp=-fp(i)
      do 20 j=1,nstu
      iax=iax+1
      loc=loc + 1
      temp=temp - cdpu(loc)*ul(iax,ino)
      if(iax.eq.2)then
                  ino=ino+1
                  iax=0
                  endif
   20 continue
      rhs(i)=temp
   10 continue
      if (nelp.ne.1) go to 100
      piv=cmpp(1)
      if (piv.eq.0.00) then
          print*,'error: kpp is not invertible - stop in elimp '
          stop
      endif
      pinc(1)=rhs(1)/piv
      return
c
  100 continue
c
c     move akpp to the working array akpp1
c
      ii=0
      do 114 i=1,nelp
      do 105 j=i,nelp
      ii=ii+1
      cmpp1(ii)=cmpp(ii)
  105 continue
  114 continue
c
      nn2=nelp + 2
      ipjp=-nelp
      do 116 ip=1,nelp-1
      ipjp=ipjp + nn2 - ip
      piv=cmpp1(ipjp)
      if (piv.eq.0.00) then
          print*,'error: kpp is not invertible - stop in elimp'
          stop
      endif
      dd=1.0/piv
c
      ii=ipjp
      ijp=ipjp
      do 110 i=ip+1,nelp
      ii=ii + nn2 - i
      ijp=ijp + 1
      fac=cmpp1(ijp)*dd
c
      jip=ijp - 1
      ij=ii - 1
      do 120 j=i,nelp
      jip=jip + 1
      ij=ij + 1
      cmpp1(ij)=cmpp1(ij) - fac*cmpp1(jip)
  120 continue
      rhs(i)=rhs(i) - fac*rhs(ip)
  110 continue
  116 continue
c
      piv=cmpp1(ij)
      if (piv.eq.0.00) then
          print*,'error: kpp is not invertible - stop in elimp'
          stop
      endif
      pinc(nelp)=rhs(nelp)/piv
      ii=ij
      do 130 i=nelp-1,1,-1
      temp=rhs(i)
      ii=ii - nn2 + i + 1
      ij=ii
      do 140 j=i+1,nelp
      ij=ij + 1
      temp=temp - cmpp1(ij)*pinc(j)
  140 continue
      pinc(i)=temp/cmpp1(ii)
  130 continue
      return
      end
c ####################################################################
      subroutine pgauss(l,lint,r,z,w)
      implicit real*8 (a-h,o-z)
      dimension lr(9),lz(9),lw(9),r(*),z(*),w(*)
      data lr/-1,1,1,-1,0,1,0,-1,0/,lz/-1,-1,1,1,-1,0,1,0,0/
      data lw/4*25,4*40,64/
      lint=l*l
      if(l.eq.1) goto 1
      if(l.eq.2) goto 2
      if(l.eq.3) goto 3
c      go to(1,2,3),l
    1 r(1)=0.
      z(1)=0.
      w(1)=4.
      return
    2 g=1./dsqrt(3.0d0)
      do 21 i=1,4
      r(i)=g*lr(i)
      z(i)=g*lz(i)
      w(i)=1.
   21 continue
      return
    3 g=dsqrt(0.6d0)
      h=1./81.
      do 31 i=1,9
      r(i)=g*lr(i)
      z(i)=g*lz(i)
      w(i)=h*lw(i)
   31 continue
      return
      end
      subroutine pstres(sig,p1,p2,p3)
      implicit real*8 (a-h,o-z)
      dimension sig(3)
      xi1=(sig(1)+sig(3))/2.
      xi2=(sig(1)-sig(3))/2.
      rho=dsqrt(xi2*xi2+sig(2)*sig(2))
      p1=xi1+rho
      p1=xi1-rho
      p3=45.0
      if(xi2.ne.0.0)p3=22.5*atan2(sig(2),xi2)/atan(1.0)
      return
      end
      subroutine shape(ss,tt,x,shp,xsj,ndm,nel,ix,flg)
      implicit real*8 (a-h,o-z)
      logical flg
      dimension shp(3,9),x(ndm,*),s(4),t(4),xs(2,2),sx(2,2),ix(*)
      data s/-0.5,0.5,0.5,-0.5/,t/-0.5,-0.5,0.5,0.5/
      do 100 i=1,4
      shp(3,i)=(0.5+s(i)*ss)*(0.5+t(i)*tt)
      shp(1,i)=s(i)*(0.5+t(i)*tt)
      shp(2,i)=t(i)*(0.5+s(i)*ss)
  100 continue
      if(nel.ge.4)goto 120
      do 110 i=1,3
      shp(i,3)=shp(i,3)+shp(i,4)
  110 continue
  120 if(nel.gt.4)call shap2(ss,tt,shp,ix,nel)
      do 132 i=1,ndm
      do 131 j=1,2
      xs(i,j)=0.0
      do 130 k=1,nel
      xs(i,j)=xs(i,j)+x(i,k)*shp(j,k)
  130 continue
  131 continue
  132 continue
      xsj=xs(1,1)*xs(2,2)-xs(1,2)*xs(2,1)
      if(flg)return
      sx(1,1)=xs(2,2)/xsj
      sx(2,2)=xs(1,1)/xsj
      sx(1,2)=-xs(1,2)/xsj
      sx(2,1)=-xs(2,1)/xsj
      do 140 i=1,nel
      tp=shp(1,i)*sx(1,1)+shp(2,i)*sx(2,1)
      shp(2,i)=shp(1,i)*sx(1,2)+shp(2,i)*sx(2,2)
      shp(1,i)=tp
  140 continue
      return
      end
      subroutine shap2(s,t,shp,ix,nel)
      implicit real*8 (a-h,o-z)
      dimension ix(*),shp(3,*)
      s2=(1.-s*s)/2.
      t2=(1.-t*t)/2.
      do 99 i=5,9
      do 100 j=1,3
      shp(j,i)=0.0
  100 continue
   99 continue
      if(ix(5).eq.0)goto 101
      shp(1,5)=-s*(1.-t)
      shp(2,5)=-s2
      shp(3,5)=-s2*(1.-t)
  101 if(nel.lt.6)goto 107
      if(ix(6).eq.0)goto 102
      shp(1,6)=t2
      shp(2,6)=-t*(1.+s)
      shp(3,6)=t2*(1.+s)
  102 if(nel.lt.7)goto 107
      if(ix(7).eq.0)goto 103
      shp(1,7)=-s*(1.+t)
      shp(2,7)=s2
      shp(3,7)=s2*(1.+t)
  103 if(nel.lt.8)goto 107
      if(ix(8).eq.0)goto 104
      shp(1,8)=-t2
      shp(2,8)=-t*(1.-s)
      shp(3,8)=t2*(1.-s)
  104 if(nel.lt.9)goto 107
      if(ix(9).eq.0)goto 103
      shp(1,9)=-4.0*s*t2
      shp(2,9)=-4.0*t*s2
      shp(3,9)=-4.*s2*t2
      do 119 j=1,3
      do 105 i=1,4
      shp(j,i)=shp(j,i)-0.25*shp(j,9)
  105 continue
      do 106 i=5,8
      if(ix(i).ne.0)shp(j,i)=shp(j,i)-.5*shp(j,9)
  106 continue
  119 continue
  107 k=8
      do 109 i=1,4
      l=i+4
      do 108 j=1,3
      shp(j,i)=shp(j,i)-0.5*(shp(j,k)+shp(j,l))
  108 continue
      k=l
  109 continue
      return
      end
c********************************************************************
      subroutine tracki(nn,ne,nrow,npoint,delt,itst)

      use dyn_arrays
      use dyn_arrays_mech
      implicit real*8 (a-h,o-z)
      dimension node1(4),node2(4),node3(4)
	common /trackivar/ ipoint,xm1,ym1,vxt1,vyt1,xmp,ymp,toler,
     * iloop,ie,iele,in,x2,y2,x3,y3,scaler,x4,y4,x5,y5,x6,y6,
     * a1,a0,b0,b1,b2,b3,r,s,vx1,vx2,vx3,vx4,vy1,vy2,vy3,vy4,a2,a3
     * vxtemp,vytemp,xm,ym,t1,t2,t3,t4,t_track
      if(npoint.eq.0)return
      node1(1)=1
      node2(1)=2
      node3(1)=3
      node1(2)=2
      node2(2)=3
      node3(2)=4
      node1(3)=3
      node2(3)=4
      node3(3)=1
      node1(4)=4
      node2(4)=1
      node3(4)=2

      do 500 ipoint=1,npoint
c
c  locally store previous tstep position and velocity
c
      xm1=tpoint(1,ipoint)
      ym1=tpoint(2,ipoint)
      vxt1=tpoint(3,ipoint)
      vyt1=tpoint(4,ipoint)
c predictor:
      xm=tpoint(1,ipoint)+delt*vxt1
      ym=tpoint(2,ipoint)+delt*vyt1
      xmp=xm
      ymp=ym
c
      toler=.005*dsqrt((xm-xm1)**2+(ym-ym1)**2)
c
c loop to find implicitly velocity at new point
c
      do 400 iloop=1,10
c
c   find current element
c
      if(ieletp(ipoint).eq.0)go to 450
c
c  search vicinity of previous element first
c
      do 200 iele=1,ne+5
      ie=iele-5
      if(iele.eq.1)ie=ieletp(ipoint)
      if(iele.eq.2)ie=min0(ieletp(ipoint)+1,ne)
      if(iele.eq.3)ie=max0(ieletp(ipoint)-1,1)
      if(iele.eq.4)ie=min0(ieletp(ipoint)+nrow-1,ne)
      if(iele.eq.5)ie=max0(ieletp(ipoint)-nrow+1,1)
      if(ie.le.0) then
      	write(6,*)' error element =  ',ie
      	print*,'   Looking for ipoint= ',ipoint,iele
      endif	
      if(ie.gt.ne) then
      	write(6,*)' error element = ',ie
      	print*,'   Looking for ipoint= ',ipoint,iele
      endif 
c
c  calculate scaler products
c
c
c  loop over each side
c
      do 250 in=1,4
      if(mod(ie,(nrow-1)).eq.1.and.in.eq.1)go to 249
      if(mod(ie,(nrow-1)).eq.0.and.in.eq.3)go to 249
      if(ie.lt.(nrow).and.in.eq.4)go to 249
      if(ie.gt.(ne-nrow+1).and.in.eq.2)go to 249
c
c  find point on edge whose normal intersects other node
c                (x4,y4)
c
      x2=coord(1,node(node2(in),ie))-coord(1,node(node1(in),ie))
      y2=coord(2,node(node2(in),ie))-coord(2,node(node1(in),ie))
      x3=coord(1,node(node3(in),ie))-coord(1,node(node1(in),ie))
      y3=coord(2,node(node3(in),ie))-coord(2,node(node1(in),ie))
      scaler=(x2*x3+y2*y3)/(x2*x2+y2*y2)
      x4=coord(1,node(node1(in),ie))+x2*scaler
      y4=coord(2,node(node1(in),ie))+y2*scaler
c
c find scaler product of edge normal and tracked point
c
      x5=coord(1,node(node3(in),ie))-x4
      y5=coord(2,node(node3(in),ie))-y4
      x6=xm-x4
      y6=ym-y4
      scaler=x5*x6+y5*y6
      if(scaler.lt.0.0)go to 190
  249 continue
  250 continue
      ieletp(ipoint)=ie
      go to 100
  190 continue
  200 continue
  450 continue
c
c  point not in any element
c
      write(6,*)' warning: lost point: ',xm,ym
  100 continue
c
c  find current velocity at track point
c
      x1=coord(1,node(1,ieletp(ipoint)))
      x2=coord(1,node(2,ieletp(ipoint)))
      x3=coord(1,node(3,ieletp(ipoint)))
      x4=coord(1,node(4,ieletp(ipoint)))
      y1=coord(2,node(1,ieletp(ipoint)))
      y2=coord(2,node(2,ieletp(ipoint)))
      y3=coord(2,node(3,ieletp(ipoint)))
      y4=coord(2,node(4,ieletp(ipoint)))
c
c  find r and s
c
c     find coefs for space interpolation
c
      a0=.5*(x1+x2)
      a1=.5*(x1-x2)
      b0=.25*(y1+y2+y3+y4)
      b1=.25*(y1-y2-y3+y4)
      b2=.25*(y1+y2-y3-y4)
      b3=.25*(y1-y2+y3-y4)
c
      r=(xm-a0)/a1
      s=(ym-b0-b1*r)/(b2+b3*r)
      if(r.gt.1.0)r=1.0
      if(r.lt.-1.0)r=-1.0
      if(s.gt.1.0)s=1.0
      if(s.lt.-1.0)s=-1.0
c
c  interpolate velocity
c
      vx1=velx(node(1,ieletp(ipoint)))
      vx2=velx(node(2,ieletp(ipoint)))
      vx3=velx(node(3,ieletp(ipoint)))
      vx4=velx(node(4,ieletp(ipoint)))
      vy1=vely(node(1,ieletp(ipoint)))
      vy2=vely(node(2,ieletp(ipoint)))
      vy3=vely(node(3,ieletp(ipoint)))
      vy4=vely(node(4,ieletp(ipoint)))
c
c     find coefs for velocity interpolation
c
      a0=.25*(vx1+vx2+vx3+vx4)
      a1=.25*(vx1-vx2-vx3+vx4)
      a2=.25*(vx1+vx2-vx3-vx4)
      a3=.25*(vx1-vx2+vx3-vx4)
      b0=.25*(vy1+vy2+vy3+vy4)
      b1=.25*(vy1-vy2-vy3+vy4)
      b2=.25*(vy1+vy2-vy3-vy4)
      b3=.25*(vy1-vy2+vy3-vy4)
c
c  find velocity at updated track point
c
      vxtemp=a0+a1*r+a2*s+a3*r*s
      vytemp=b0+b1*r+b2*s+b3*r*s


c#################################################
c#####  track temp for lagrangian mesh ##########
c#################################################

c determine temperature at eulerian nodes surrounding element 
c containing lagrangian node       
      t1=temptc(node(1,ieletp(ipoint)))
      t2=temptc(node(2,ieletp(ipoint)))
      t3=temptc(node(3,ieletp(ipoint)))
      t4=temptc(node(4,ieletp(ipoint)))


c find coefs for temperature interpolation

      a0=.25*(t1+t2+t3+t4)
      a1=.25*(t1-t2-t3+t4)
      a2=.25*(t1+t2-t3-t4)
      a3=.25*(t1-t2+t3-t4)
      
c find temperature at updated track point

      t_track=a0+a1*r+a2*s+a3*r*s
c####################################################
c ###  track ductile/plastic deformation for lmesh ##
c####################################################
c after finding correct element, check if any gauss point
c	has ductile deformation
      rigmax=.004
      iduc=0
      do 376 igp=1,4
      	if(ipflag(ieletp(ipoint),igp).eq.0) then
c		check 2nd invar of strain rate: if above rigmax, behave ductile
c			if below rigmax, behave rigid (only want ductile)
c			record ductile behavior as tpoint(6,ipoint) as =1
c			in tpoint(7,ipoint) record the number of timesteps 
c			of ductile def.
      		secdef=(srate(ieletp(ipoint),1)*srate(ieletp(ipoint),1)+
     *		srate(ieletp(ipoint),2)*srate(ieletp(ipoint),2))/2.
     *		+srate(ieletp(ipoint),3)*srate(ieletp(ipoint),3)
     		if(secdef.lt.0.0)then
      			secdef=0.0
      		endif
      		secdef=dsqrt(secdef)
      		if(secdef.ge.rigmax) then
     			tpoint(6,ipoint)=1.0
      			tpoint(7,ipoint)=tpoint(7,ipoint)+1.0
      			goto 377
      		endif	
      	endif
  376 continue	
  377 continue      
c
c  update position of material point
c
c    check for convergence in velocity-position loop
c
c  predictor
c
c
c corrector
c

c#### ensure that all lmesh nodes to the lhs of the domain do not
c####	have a vert component of displacement.  this is important 
c####	with the edge bcs adjusted so that they are tangential to the
c####	base of the model. if there is a vert component to the velocity 
c####	at the model edge, this causes the mesh outside the domain to
c####	advect downwards

      if(xm.lt.coord(1,1)) then
c		no vert displacement outside the model      	
      	ym=ym1
      else
c		vert displacement inside the model      
		ym=ym1+.5*(vyt1+vytemp)*delt
	  end if	

      xm=xm1+.5*(vxt1+vxtemp)*delt
c comment this out since no vert disp outside model domain catch is above
c      ym=ym1+.5*(vyt1+vytemp)*delt
      ctest=dsqrt((xm-xmp)**2+(ym-ymp)**2)
      if(ctest.gt.toler)go to 390
      go to 410
  390 continue
      xmp=xm
      ymp=ym
  400 continue
  410 continue
      tpoint(1,ipoint)=xm
      tpoint(2,ipoint)=ym
      tpoint(3,ipoint)=vxtemp
      tpoint(4,ipoint)=vytemp
      tpoint(5,ipoint)=t_track
c
c  routine to calculate the exhumation rate for each l-node
c
      call exhum_rate(xm,ym,xm1,ym1,nn,nrow,exhum,ipoint,
     *xsur,xsurold,itst,delt)

  500 continue
      return
      end
c ################################################################
c material trcking for basin surfaces
c ################################################################
      subroutine track_basin(nn,ne,nrow,npoint,delt,itst)

      use dyn_arrays
      use dyn_arrays_mech
      implicit real*8 (a-h,o-z)
      dimension node1(4),node2(4),node3(4)
	common /trackivar/ ipoint,xm1,ym1,vxt1,vyt1,xmp,ymp,toler,
     * iloop,ie,iele,in,x2,y2,x3,y3,scaler,x4,y4,x5,y5,x6,y6,
     * a1,a0,b0,b1,b2,b3,r,s,vx1,vx2,vx3,vx4,vy1,vy2,vy3,vy4,a2,a3
     * vxtemp,vytemp,xm,ym,t1,t2,t3,t4,t_track
      if(npoint.eq.0) return
      node1(1)=1
      node2(1)=2
      node3(1)=3
      node1(2)=2
      node2(2)=3
      node3(2)=4
      node1(3)=3
      node2(3)=4
      node3(3)=1
      node1(4)=4
      node2(4)=1
      node3(4)=2

      do 500 ipoint=1,npoint
c
c  locally store previous tstep position and velocity
c
      xm1=bastrk(1,ipoint)
      ym1=bastrk(2,ipoint)
      vxt1=bastrk(3,ipoint)
      vyt1=bastrk(4,ipoint)
c predictor:
      xm=bastrk(1,ipoint)+delt*vxt1
      ym=bastrk(2,ipoint)+delt*vyt1
      xmp=xm
      ymp=ym

      toler=.005*dsqrt((xm-xm1)**2+(ym-ym1)**2)

c loop to find implicitly velocity at new point
c
      do 400 iloop=1,10
c
c   find current element
c
      if(ieletpb(ipoint).eq.0)go to 450
c
c  search vicinity of previous element first
c
      do 200 iele=1,ne+5
      ie=iele-5
      if(iele.eq.1)ie=ieletpb(ipoint)
      if(iele.eq.2)ie=min0(ieletpb(ipoint)+1,ne)
      if(iele.eq.3)ie=max0(ieletpb(ipoint)-1,1)
      if(iele.eq.4)ie=min0(ieletpb(ipoint)+nrow-1,ne)
      if(iele.eq.5)ie=max0(ieletpb(ipoint)-nrow+1,1)
      if(ie.le.0) then
      	write(6,*)' error basin element =  ',ie
      	print*,'   Looking for ipoint= ',ipoint,iele
      endif	
      if(ie.gt.ne) then
      	write(6,*)' error basin element = ',ie
      	print*,'   Looking for ipoint= ',ipoint,iele
      endif 
c
c  calculate scaler products
c
c
c  loop over each side
c
      do 250 in=1,4
      if(mod(ie,(nrow-1)).eq.1.and.in.eq.1)go to 249
      if(mod(ie,(nrow-1)).eq.0.and.in.eq.3)go to 249
      if(ie.lt.(nrow).and.in.eq.4)go to 249
      if(ie.gt.(ne-nrow+1).and.in.eq.2)go to 249
c
c  find point on edge whose normal intersects other node
c                (x4,y4)
c
      x2=coord(1,node(node2(in),ie))-coord(1,node(node1(in),ie))
      y2=coord(2,node(node2(in),ie))-coord(2,node(node1(in),ie))
      x3=coord(1,node(node3(in),ie))-coord(1,node(node1(in),ie))
      y3=coord(2,node(node3(in),ie))-coord(2,node(node1(in),ie))
      scaler=(x2*x3+y2*y3)/(x2*x2+y2*y2)
      x4=coord(1,node(node1(in),ie))+x2*scaler
      y4=coord(2,node(node1(in),ie))+y2*scaler
c
c find scaler product of edge normal and tracked point
c
      x5=coord(1,node(node3(in),ie))-x4
      y5=coord(2,node(node3(in),ie))-y4
      x6=xm-x4
      y6=ym-y4
      scaler=x5*x6+y5*y6
      if(scaler.lt.0.0)go to 190
  249 continue
  250 continue
      ieletp(ipoint)=ie
      go to 100
  190 continue
  200 continue
  450 continue
c
c  point not in any element
c
      write(6,*)' warning: lost basin point: ',xm,ym
  100 continue
c
c  find current velocity at track point
c
      x1=coord(1,node(1,ieletp(ipoint)))
      x2=coord(1,node(2,ieletp(ipoint)))
      x3=coord(1,node(3,ieletp(ipoint)))
      x4=coord(1,node(4,ieletp(ipoint)))
      y1=coord(2,node(1,ieletp(ipoint)))
      y2=coord(2,node(2,ieletp(ipoint)))
      y3=coord(2,node(3,ieletp(ipoint)))
      y4=coord(2,node(4,ieletp(ipoint)))
c
c  find r and s
c
c     find coefs for space interpolation
c
      a0=.5*(x1+x2)
      a1=.5*(x1-x2)
      b0=.25*(y1+y2+y3+y4)
      b1=.25*(y1-y2-y3+y4)
      b2=.25*(y1+y2-y3-y4)
      b3=.25*(y1-y2+y3-y4)
c
      r=(xm-a0)/a1
      s=(ym-b0-b1*r)/(b2+b3*r)
      if(r.gt.1.0)r=1.0
      if(r.lt.-1.0)r=-1.0
      if(s.gt.1.0)s=1.0
      if(s.lt.-1.0)s=-1.0
c
c  interpolate velocity
c
      vx1=velx(node(1,ieletp(ipoint)))
      vx2=velx(node(2,ieletp(ipoint)))
      vx3=velx(node(3,ieletp(ipoint)))
      vx4=velx(node(4,ieletp(ipoint)))
      vy1=vely(node(1,ieletp(ipoint)))
      vy2=vely(node(2,ieletp(ipoint)))
      vy3=vely(node(3,ieletp(ipoint)))
      vy4=vely(node(4,ieletp(ipoint)))
c
c     find coefs for velocity interpolation
c
      a0=.25*(vx1+vx2+vx3+vx4)
      a1=.25*(vx1-vx2-vx3+vx4)
      a2=.25*(vx1+vx2-vx3-vx4)
      a3=.25*(vx1-vx2+vx3-vx4)
      b0=.25*(vy1+vy2+vy3+vy4)
      b1=.25*(vy1-vy2-vy3+vy4)
      b2=.25*(vy1+vy2-vy3-vy4)
      b3=.25*(vy1-vy2+vy3-vy4)
c
c  find velocity at updated track point
c
      vxtemp=a0+a1*r+a2*s+a3*r*s
      vytemp=b0+b1*r+b2*s+b3*r*s

      xm=xm1+.5*(vxt1+vxtemp)*delt
      ym=ym1+.5*(vyt1+vytemp)*delt
      ctest=dsqrt((xm-xmp)**2+(ym-ymp)**2)
      if(ctest.gt.toler)go to 390
      go to 410
  390 continue
      xmp=xm
      ymp=ym
  400 continue
  410 continue
      bastrk(1,ipoint)=xm
      bastrk(2,ipoint)=ym
      bastrk(3,ipoint)=vxtemp
      bastrk(4,ipoint)=vytemp
  500 continue
      return
      end



c#################################################################
      subroutine exhum_rate(xm,ym,xm1,ym1,nn,nrow,exhum,ipoint,
     *xsur,xsurold,itst,delt)
      implicit real*8 (a-h,o-z)
      dimension xsur(2,*),exhum(*),xsurold(2,*)

      if(itst.eq.1) then
      	exhum(ipoint)=0.0
      	return
      endif	
      if(itst.eq.2) then
      	exhum(ipoint)=0.0
      	return
      endif	

      ncolm=nn/nrow

c if l-node at prev is L of model edge, set =0 
      if(xm1.lt.0.0) then
      	exhum(ipoint)=0.0
      	return
      endif

c determine depth below surface (down is +) for previous tst
      do 200 i=1,ncolm
      	if(xm1.lt.xsurold(1,i)) then
      		dx=xsurold(1,i)-xsurold(1,i-1)
      		dy=xsurold(2,i)-xsurold(2,i-1)
      		xdiff=xm1-xsurold(1,i-1)
      		ydepth1=xsurold(2,i-1)+xdiff*dy/dx-ym1
      		if(ydepth1<0.0) then
      			exhum(ipoint)=0.0
      			return
      		endif	
      		goto 215
      	endif
  200 continue
  215 continue

c determine depth below surface (down is +) for current tst
      do 300 i=1,ncolm
      	if(xm.lt.xsur(1,i)) then
      		dx=xsur(1,i)-xsur(1,i-1)
      		dy=xsur(2,i)-xsur(2,i-1)
      		xdiff=xm-xsur(1,i-1)
      		ydepth=xsur(2,i-1)+xdiff*dy/dx-ym
      		goto 315
      	endif
  300 continue
  315 continue
  
      exhum(ipoint)=(ydepth1-ydepth)/(delt*1000.0)

      return
      end
