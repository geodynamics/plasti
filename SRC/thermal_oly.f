c   SUBROUTINE THERMAL (a derivative of the
c                 program heatran)

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

c#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c  ne = Total number of elements
c  nn = total number of nodes
c  n = total number of boundary nodes?
c
c  coordcoordcoordt,Z, coordinates
c  nodet(ne,j)  - J=1,2,3 Node numbers for the Ith element
c		 J=4 Propery Map Value
c                J=5 Velocity Domain
c
c  vx(ne)	- velocity in x-dir
c  vz(ne)	- velocity in z-dir
c  velx(nn)	- velocity in x-dir from plasti
c  vely(nn)	- velocity in z-dir from plasti
c  asf(ne,3)	- shape function A coefficents
c  bsf(ne,3)	- shape function B coefficents
c  area(ne)		- Areas of elements
c  tcond(2,ne)	- Rock thermal conductivities (anisotropic)
c  trho(ne) 	- Rock densities
c  spheat(ne)	- Rock specific heats
c  hprod(ne)	- Heat production
c  tempt(nn)	- Temperatures
c  temp(ne/2)  - Temperatures averaged over quad elements
c  told(nn)	- Temperatures at previous time step
c  ntbnd(n)	- Constant temperature nodes
c  btem(n)	- Constant temperature values
c  neflux(n,2)	- Constant heat flux element and nodes on that side
c			=(1,2, and/or 3)
c  flux(n)	- Constant heat flux value
c  iter         - Iteration number
c  itst         - Time step number
c#cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine thermal(deltt,itst,nn,ne,nout,nrowp,nrow,ncol,
     *lda,lbw,ntbn,ioutpt)

      use dyn_arrays
      use dyn_arrays_therm
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      character date*10,time*10,time2*10
      integer quadcount

c ALLOCATE LOCAL ARRAY STORAGE
	  allocate(a(lda,nn),asf(ne,3),bsf(ne,3),area(ne),ipt(nn))
	  allocate(rhst(nn))

      call date_and_time(date,time)
      print*,'Real Time Entering Thermal is:  ',time
      print*,'Time Step Number =',itst
      print*,'Delt = ',deltt
      deltt=(3.15578e13)*deltt

c first (0) interation sets up initial conditions
      if(itst.eq.0) then
       	vx=0.0
      	vz=0.0
      endif      

c SET UP SHAPE FUNCTION COEFFICIENTS
      call sfcoef(nn,ne,itst)

c ASSEMBLE GLOBAL STIFFNESS MATRIX FOR HEAT TRANSPORT PROBLEM
      call globet(ne,nn,lbw,lda,deltt,itst)

c APPLY BOUNDARY CONDITIONS
c	NOTE: flux bcs are not implemented in the platis meshg,
c			so set nfel == to 0 here
      nfel=0
      call bct(nn,ne,ntbn,nfel,lda,lbw)

c SOLVE SYSTEM OF EQUATIONS (LAPAK ROUTINES)	  
      call dgbtrf(nn,nn,lbw,lbw,a,lda,ipt,info)
      if(info.ne.0) then
	  	print*,'#####  ERROR IN FACTORIZATION, THERMAL DGBTRF'
	  	print*,'info from dgbtrf',info
       	call outputt(nn,ne,itst,nout,nrow,ncol,ioutpt)
	  	stop
	  endif
	  call dgbtrs('N',nn,lbw,lbw,1,a,lda,ipt,rhst,nn,info)
	  if(info.ne.0) then
	  	print*,'#####  ERROR IN FACTORIZATION, THERMAL DGBTRS'
	  	print*,'info from dgbtrs',info
	  	call outputt(nn,ne,itst,nout,nrow,ncol,ioutpt)
	  	stop
	  endif

c STORE NEW TEMPERATURES
      do i=1,nn
      	tempt(i)=rhst(i)
      end do

c  OUTPUT RESULTS AFTER CONVERGENCE OR SPECIFIED NUMBER OF ITERATIONS
       call outputt(nn,ne,itst,nout,nrow,ncol,ioutpt)  

c  AVERAGE TEMP OVER QUAD-ELEMENTS
      quadcount=0
      k=0
      do j=2*(nrow-nrowp),ne,2*(nrow-1)
      	do i=1,2*(nrowp-1),2
      		k=i+j
      		quadcount=quadcount+1
      		temp(quadcount)=(tempt(nodet(k,1))+tempt(nodet(k,2))+
     *		tempt(nodet((k+1),1))+tempt(nodet((k+1),3)))
      		temp(quadcount)=temp(quadcount)/4
      	end do
      end do

      call date_and_time(date,time2)
      print*,'Real Time Leaving Thermal is:  ',time2
      deallocate(a,asf,bsf,area,ipt,rhst)
      return
      end
c ######################################################################
c ######################################################################
c************************************************************
c*   SUBROUTINE TO OUTPUT RESULTS                           *
c************************************************************
      subroutine outputt(nn,ne,itst,nout,nrow,ncol,ioutpt)

      use dyn_arrays
      use dyn_arrays_therm
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer thdpl,fstpl,secpl,temp1
      character(30):: coordt_op='coordt_',veltherm_alt_op=
     *'velthermal_alt_',velthermal_op='velthermal_',temp_op='temp_',
     *matp_hprod_op='matp_hprod_',matp_tcond_y_op='matp_tcond_y_',
     *matp_spec_ht_op='matp_spec_ht_',dir,fextn
      character(10):: nums='0123456789'

      if(itst.eq.99999)return
      if(itst.eq.0)return
      if(itst.eq.1) goto 39
      if(nout.eq.1) goto 39
c catch for no output to equilibrate thermal model to subduction      
      itest=mod(itst,nout)
      if(itest.ne.0)return
   39 continue
c      print*,' Inside Thermal Output, Timestep=',itst

c output directory
      dir='output/'
c determine extension for output file names
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
c coords
c#############
      if(output_flags(38).eq.1) then
      	open(19,file=trim(dir)//trim(coordt_op)//trim(fextn),
     *	position='rewind')
      	write(19,101)nn
      	do i=1,nn
      		write(19,105)coordt(1,i),coordt(2,i)
      	end do
      	close(19)
      endif	

c#############
c reduced resolution thermal vel
c#############
      if(output_flags(39).eq.1) then
      	open(14,file=trim(dir)//trim(veltherm_alt_op)//trim(fextn),
     *	position='rewind')
      	write(14,104)ne/2,ncol,nrow
      	do i=1,ne,2
        	write(14,103)vx(i)*3.15578e13
      	end do
      	do i=1,ne,2
        	write(14,103)vz(i)*3.15578e13
        end do
        close(14)
      endif  

c############
c full resolution velocities
c############
      if(output_flags(40).eq.1) then
      	open(16,file=trim(dir)//trim(velthermal_op)//trim(fextn),
     *	position='rewind')
      	write(16,101)ne,ncol,nrow
      	do i=1,ne
      		write(16,103)vx(i)*3.15578e13
      	end do
      	do i=1,ne
      		write(16,103)vz(i)*3.15578e13

      	end do	
      	close(16)
      endif	
      
c#############
c temp
c#############
      if(output_flags(41).eq.1) then
      	open(18,file=trim(dir)//trim(temp_op)//trim(fextn),
     *	position='rewind')
      	write(18,101) nn
      	do i=1,nn,8
      		count=(nn-i)/8
      		if(count.ge.1) then
      			j=i+7
      		else
      			j=nn
      		endif
      		write(18,102)(tempt(k)-273,k=i,j)
      	end do
      	close(18)
      endif

c#############
c thermal material props
c#############
c heat production      	
      if(output_flags(42).eq.1) then	
      	open(21,file=trim(dir)//trim(matp_hprod_op)//trim(fextn),
     *	position='rewind')
      	write(21,101)ne
      	do i=1,ne
      		write(21,103)hprod(i)
      	end do	
      	close(21)
      endif	
c thermal cond y-dir   	
      if(output_flags(43).eq.1) then	
      	open(21,file=trim(dir)//trim(matp_tcond_y_op)//trim(fextn),
     *	position='rewind')
      	write(21,101)ne
      	do i=1,ne
      		write(21,103)tcond(2,i)
      	end do
      	close(21)
      endif      
c specific heat
      if(output_flags(44).eq.1) then	
      	open(21,file=trim(dir)//trim(matp_spec_ht_op)//trim(fextn),
     *	position='rewind')
      	write(21,101)ne
      	do i=1,ne
      		write(21,103)spheat(i)
      	end do	
      	close(21)
      endif      

  101 format(i6)
  102 format(13f12.2)
  103 format(2e15.9)
  104 format(3i6)
  105 format(SP,6e12.6,/5e12.6)
      return
      end

c********************************************************************
c*      SUBROUTINE TO APPLY BOUNDARY CONDITIONS                     *
c********************************************************************
      subroutine bct(nn,ne,ntbn,nfel,lda,lbw)

      use dyn_arrays
      use dyn_arrays_therm
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

c      print*,'nfel=',nfel,'ntbn=',ntbn
c      print*,'lda=',lda,'lbw=',lbw
      m=2*lbw+1

c   apply constant flux bc
      if(nfel.gt.0) then
      	do in=1,nfel
      		n1=neflux(in,1)
      		n2=neflux(in,2)
c  find length of element side
      		xl=(coordt(1,n1)-coordt(1,n2))**2
      		yl=(coordt(2,n1)-coordt(2,n2))**2
      		tl=dsqrt(xl+yl)
c  calculate nodal value and insert in rhs
      		fln=flux(in)*tl/2.000000
      		rhst(n1)=fln+rhst(n1)
      		rhst(n2)=fln+rhst(n2)
      	end do
      endif	

c  apply constant value boundary condition
      if(ntbn.eq.0) go to 301
      do in=1,ntbn
      	ib=ntbnd(in)
      	llb=ib-lbw
      	iub=ib+lbw
      	if(llb.lt.1)llb=1
      	if(iub.gt.nn)iub=nn
c  set corresponding row of global stiffness matrix to 0
      	do jb=llb,iub
      		kb=ib-jb+m
      		a(kb,jb)=0.0
      	end do
c  set principle diagonal component to 1.
c       and rhs to prescribed value
      	a(m,ib)=1.00000
      	rhst(ib)=btem(in)
      end do
  301 continue 

      return
      end
c**************************************************************
c*  ROUTINE TO ASSEMBLE GLOBAL STIFFNESS MATRIX AND RHS       *
c*         FOR HEAT TRANSPORT EQUATION                        *
c**************************************************************
      subroutine globet(ne,nn,lbw,lda,deltt,itst)

      use dyn_arrays_therm
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 massk,massb
      Dimension massb(3),massk(3,3),
     *t(3),s(3,3)

c initialize stiffness matrix and rhs

      mbw=2*lbw+1
      ma=(3*lbw+1)
      do j=1,nn
      	rhst(j)=0.0
      	do i=1,ma
      		a(i,j)=0.0
      	end do
      end do

c loop over each element
c calculate element stiffness matrix
      do 100 iele=1,ne
c conductivity for the element
c dispersion tensor
      dxx=tcond(1,iele)
      dzz=tcond(2,iele)
      dxz=0.0

c  assemble element stiffness matrix
      do jj=1,3
      	aj=asf(iele,jj)
      	bj=bsf(iele,jj)
      	massb(jj)=hprod(iele)*area(iele)/3.0d0
      	term2=trho(iele)*spheat(iele)*(vx(iele)*aj+
     *	vz(iele)*bj)/6.0000
      	do ii=1,3
      		ai=asf(iele,ii)
      		bi=bsf(iele,ii)
      		massk(ii,jj)=(dxx*ai*aj+dxz*(ai*bj+bi*aj)+dzz*bi*bj)/(4.*
     *		area(iele))+term2
      	end do
      end do

c Calculate Transient Component To E.S.M. and RHS

c  Define density and specific heat for the rock
	if(itst.eq.0) goto 31
      denr=trho(iele)
      spec=spheat(iele)

c  Calculate bulk volume specific heat for the element
      bhc=denr*spec

c  Define nodal temperatures
      do ii=1,3
      	t(ii)=told(nodet(iele,ii))
      end do

c  Define Transient Mass (s) Matrix

      do 15 ii=1,3
      do 110 jj=1,3
      s(ii,jj)=bhc*area(iele)/12
  110 continue
      s(ii,ii)=s(ii,ii)*2
   15 continue

c  Insert into ESM(element stiffness matrix) and RHS

      do 30 ii=1,3
      do 20 jj=1,3
      massb(ii)=massb(ii)+s(ii,jj)*t(jj)/deltt
      massk(ii,jj)=massk(ii,jj)+s(ii,jj)/deltt
   20 continue
   30 continue
   31	continue

c assemble global stiffness matrix and rhs

      do 41 l=1,3
      i=nodet(iele,l)
      rhst(i)=rhst(i)-massb(l)
      do 40 m=1,3
      j=nodet(iele,m)
      k=i-j+mbw
      a(k,j)=a(k,j)-massk(l,m)
   40 continue
   41 continue
  100 continue

      return 
      end



c****************************************************************
c*    SUBROUTINE TO SET UP SHAPE FUNCTION COEFFICIENTS          *
c*                 FOR EACH ELEMENT                             *
c****************************************************************
      subroutine sfcoef(nn,ne,itst)

      use dyn_arrays_therm
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      do ie=1,ne
      	asf(ie,1)=coordt(2,nodet(ie,2))-coordt(2,nodet(ie,3))
      	asf(ie,2)=coordt(2,nodet(ie,3))-coordt(2,nodet(ie,1))
      	asf(ie,3)=coordt(2,nodet(ie,1))-coordt(2,nodet(ie,2))
      	bsf(ie,1)=coordt(1,nodet(ie,3))-coordt(1,nodet(ie,2))
      	bsf(ie,2)=coordt(1,nodet(ie,1))-coordt(1,nodet(ie,3))
      	bsf(ie,3)=coordt(1,nodet(ie,2))-coordt(1,nodet(ie,1))
      	area(ie)=(asf(ie,1)*bsf(ie,2)-bsf(ie,1)*asf(ie,2))/2.0000
      	area(ie)=dabs(area(ie))
      end do
      return
      end

