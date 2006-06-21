c   mesh generator and parameter  input for thermal and mechanical mesh
c		used in plasti
c

c####################################################################
c define arrays that will be dynamically allocated in subroutines
      module dyn_arrays
      real(kind=8),allocatable::ypltop(:),xpl(:),yplbase(:),yrlbase(:),
     *xrlbase(:),xltemp(:),ranode(:,:),ymtop(:),yltemp(:),xmbase(:),
     *coh(:),slen1(:),slen2(:),fnode1(:),dyinit1(:),dyinit2(:),
     *fnode2(:),bc(:,:),dyinit(:),xp1(:),xp2(:),yp1(:),phi(:),
     *panode(:,:),rhoc(:),yp2(:),ymbase(:),rlnode(:,:),xbase(:),
     *cmnode(:,:),plnode(:,:),pos(:,:),therm_prop(:,:),thermbcs(:,:),
     *tm_prop(:,:),therm_cond(:,:),therm_rho(:),heat_prod(:),
     *spec_heat(:),dencol(:),therm_bc(:,:),vmin(:),q(:),prex(:),
     *expn(:)
      integer,allocatable::domain(:),ndomain(:),mecht_nodes(:),
     *plithb_nodes(:),plitht_nodes(:),mechb_nodes(:),rlithb_nodes(:),
     *output_flags(:)
      end module dyn_arrays
c end of definitions
c####################################################################

      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

c read in input file

      call input(ncol,nerowm,nstype,sing,
     *pmthick,plthick,athick,rlthick,ypmbase,yrmbase,wdepth,
     *wtoler,npad,xpad,prigp,rrigp,prigi,rrigi,sload,smomen,xadd,ctoler,
     *plvel,upvel,iunflag,iunbeg,xunbeg,vrig,
     *beta,epsinv,rhof,rhom,numvebn,numpbn,numsid,numvtbn,
     *ntst,intout,intoutl,delt,minint,maxint,npass,toler,erosl,
     *erosr,peros,rpow,ntt2,deltt2,iso,ntmchg,plscale,rlscale,blscale,
     *dfact,slpmax,tmax,numvetbn,ioutflag,inflag,dyc,linflag,sdip,
     *ipflag,itrench,iplasflg,iblay,iblayt,isedl,isedr,iexflg,
     *ibasflg,nbastary,nbastind,intmrkb,ipkfill,ibasfill,sedmax,iflgcl,
     *agecl,iflgblt,iflgbl,tblayt,tblay,noutput)

c find the spoint (used spoint is defined by position instead of node
      call find_node(nsing,nstype,sing,ncol,npad,xsing)
      print*,'s-point (desired, node, location): ',sing,nsing,xsing

C Make model outlines/boundaries 
      if(iso.ne.2) then
      	print*,'#########  STOP!!!   ############'
      	print*,'## Only have two plate case    ##'
      	print*,'## implemented, cannot do only ##'
      	print*,'## one plate.                  ##'
      	print*,'#################################'
      	stop
      endif
      if(iso.eq.2) then
c make x,y array for each plate
      	call mk_plates(ncol,npad,xsing,xadd,np1,nsing,nsing1,np2)
      	if(ipflag.eq.0) then
c  		calculate flexural profile 
      		call calc_flex(nerowm,ncol,np1,np2,prigi,rrigi,rhom,
     *		npad,nsing,nsing1,ctoler,sload,smomen,wdepth,wtoler,rhof,
     *		ypmbase,yrmbase,wheight,inflag,dyc)
    	elseif(ipflag.eq.1) then
			call arc_prof(nerowm,ncol,np1,np2,prigp,rrigp,rhom,
     *			npad,nsing,nsing1,wdepth,ypmbase,yrmbase,wheight,
     *			inflag,dyc,itrench,sdip)
    	endif
c output flexural profiles of just the two plates
      	call profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
c make arrays of flexure profile for model (mech and sub lithos)
      	call mech_bndry(ncol,nsing1,nsing,np1,npad,np2,plthick,
     *	athick,yshift)
      endif
c make lithosphere domain boundaries
      call lith_bndry(nsing1,nsing,plthick,athick,ncol,npltop,plvel,
     *nplbase,rlthick,irlbeg,nrlbase)
c output the domain boundaries for plotting     
      call bndry_output(wheight,ncol,npad,wdepth,np1,np2,yshift
     *,nsing,nsing1,npltop,nrlbase,nplbase)

c only make and output mesh if desired     
      if(ioutflag.eq.1) then
c 		make nodes in mech model
      	call mk_cmnodes(ncol,nerowm,iblay,iblayt,iflgblt,iflgbl,
     *		tblayt,tblay)
c 		make nodes in pro lithosphere
      	call mk_plnodes(npltop,nplbase,nplrow,base)
c 		make nodes in retro lithos      
      	call mk_rlnodes(nsing,nrlrow,irlbeg,ncol)
c 		make nodes in retro asthenosphere      
      	call mk_ranodes(nrlrow,nsing,npltop,nrarow,ncol,base,
     *	irlbeg,irabeg)
c 		make nodes in pro-asthenosphere
      	call mk_panodes(nplbase,nsing,nparow,base)
c 		make array of all nodes
      	call mk_node_array(nerowm,nrowl,nrowa,nplrow,nalrow,
     *	ncol,ntrow,nplbase,nsing,nparow,irabeg,npltop)
c 		make arrays of thermal properties and BCs
      	call mk_therm_para(ntrow,ncol,ntmchg,nrowl,nrowa,ntbcs,
     *	iflgcl,agecl,nplbase,npltop)
c 		output parameters and mesh
      	call output(ncol,nerowm,ntrow,nsing,plvel,upvel,
     *	iunflag,iunbeg,vrig,beta,epsinv,rhof,
     *	rhom,iso,prigi,rrigi,sload,smomen,xadd,ctoler,wdepth,wtoler,
     *	numvebn,numpbn,numsid,numvtbn,ntst,delt,intout,intoutl,minint,
     *	maxint,npass,toler,erosl,erosr,peros,rpow,ntt2,deltt2,np1,np2,
     *	nsing1,npad,nplbase,npltop,plscale,rlscale,blscale,dfact,slpmax,
     *	tmax,nrowl,plthick,numvetbn,linflag,iplasflg,iblay,iblayt,isedl,
     *	isedr,iexflg,ibasflg,nbastary,nbastind,intmrkb,ipkfill,ibasfill,
     *	sedmax,ntbcs,noutput)

      	deallocate(panode,xmbase,xbase,coh,phi,rhoc,bc,xp1,yp1,xp2,yp2,
     *	dyinit,yrlbase,xrlbase,xpl,ypltop,yplbase,ymtop,ymbase,xltemp,
     *	yltemp,rlnode,plnode,cmnode,ranode,pos,dyinit1,dyinit2,domain,
     *	ndomain,thermbcs,tm_prop,therm_cond,therm_rho,spec_heat,
     *	heat_prod,vmin,q,prex,expn)
      endif

      end
c########################################################
c############### END OF MAIN PROGRAM ####################
c########################################################


c########################################################
c make array of all nodes in mesh
c########################################################
      subroutine mk_node_array(nerowm,nrowl,nrowa,nplrow,
     *nalrow,ncol,ntrow,nplbase,nsing,nparow,irabeg,npltop)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer plitht(200)
c allocate space for the storage of node numbers of domain boundaries
c	used in output and used in plasti2dx to make boundaries 
c		for dx plotting
      allocate(mecht_nodes(ncol),plithb_nodes(nplbase),
     *plitht_nodes(npltop),mechb_nodes(ncol),
     *rlithb_nodes(ncol-nsing))
c initialize indicies for the above arrays
c     top of mech model
      imtnodes=0
c     base of pro lith      
      iplbnode=0
c     top of pro lith
      ipltnode=0
c     base of retro lith
      irlbnode=0
c     top of retro lith
      imbnode=0
      
c define number of rows/nodes in mesh
c     number of rows on pro lith side (-1 since top node is in mech)      
      nrowl=nplrow-1
c     number of rows in pro asthen (-1 since top node is in lith)
      nrowa=nparow-1
c     total number of rows
      ntrow=nerowm+nrowl+nrowa
      allocate(pos(ntrow*ncol,2))
c array for element domain defs      
      allocate(domain((ntrow)*(ncol)))
      allocate(ndomain((ntrow-1)*(ncol-1)*2))
c indicies for domains
      idomain=0
      index=0
      indexap=0
      indexlp=0
      indexm=0
c     these start at 1 since corner node is not counted in the domain      
      indexlr=1
      indexar=1

c##########################################################
c Determine which sub geometry this is
c	icase 1: pro-asthen ends at spoint
c	icase 2: pro-asthen ends before the spoint
c	icase  : pro-asthen ends after the spoint
c		  3: pro-asthen ends after retro-asthen begins
c		  4: pro-asthen ends when retro-asthen begins
c		  5: pro-asthen ends after retro-lith begins
c		  6: pro-asthen ends 1 colm before retro-asthen begins
c##########################################################
      icase=0
      if(nplbase.eq.nsing) then
      	icase=1
      else if(nsing.gt.nplbase) then
      	icase=2
      else if(nsing.lt.nplbase) then
      	if(irabeg.eq.nplbase) then
      		icase=4
      	else if(irabeg.eq.nplbase-1) then
      		icase=6
      	else if(irabeg.lt.nplbase-1) then
      		icase=3
      	else if(irabeg.gt.nplbase) then
      		icase=5
      	endif
      endif
      if(icase.eq.0) then
      	print*,'##################################'
      	print*,'## ERROR: Case not determined   ##'
      	print*,'##################################'
      endif	

c####################################################
c always begin with meshing to spoint (except case=1)
c	loop to/including the spoint
c####################################################

c     allow for a change in the number of rows (1=yes, 0=no)
      itog_ap=0; itog_lp=0; itog_lr=0; itog_ar=0;
c     starting value
      nshiftap=0; nshiftlp=0; nshiftar=0; nshiftlr=0;
c     include these domains (yes=1, no=0)
      iapon=1; ilpon=1; iaron=0; ilron=0;
      call mk_array(1,nrowa,iapon,1,nrowl,ilpon,1,nerowm,0,0,
     *iaron,1,0,ilron,1,nsing,nshiftap,nshiftlp,nshiftar,nshiftlr,
     *index,indexap,indexlp,indexm,indexar,indexlr,itog_ap,itog_lp,
     *itog_ar,itog_lr,idomain,imtnode,iplbnode,ipltnode,irlbnode,
     *imbnode)
      if(icase.eq.5) then
c#############
c Case 5
c#############
c loop from spoint+1 to final node of pro-astheno
      	nshiftap=0; nshiftlp=0; nshiftar=0; nshiftlr=0;
      	itog_ap=1; itog_lp=0; itog_ar=0; itog_lr=1;
      	iapon=1; ilpon=1; iaron=0; ilron=1;
      	index2=index
      	call mk_array(1,nrowa,iapon,1,nrowl,ilpon,1,nerowm,0,0,
     *  iaron,1,0,ilron,nsing+1,nplbase-1,nshiftap,nshiftlp,nshiftar,
     *  nshiftlr,
     *  index,indexap,indexlp,indexm,indexar,indexlr,itog_ap,itog_lp,
     *	itog_ar,itog_lr,idomain,imtnode,iplbnode,ipltnode,irlbnode,
     *	imbnode)
c loop till/including beg of retro-astheno
      	nshiftlr=nshiftlr; nshiftap=0; nshiftlp=0; nshiftar=0;
      	itog_ap=0; itog_lp=1; itog_lr=1; itog_ar=0;
      	iapon=0; ilpon=1; iaron=0; ilron=1;
      	call mk_array(0,0,iapon,1,nrowl+1,ilpon,1,nerowm,0,0,
     *	iaron,1,0,ilron,nplbase,irabeg,nshiftap,nshiftlp,nshiftar,
     *  nshiftlr,
     *	index,indexap,indexlp,indexm,indexar,indexlr,itog_ap,itog_lp,
     *	itog_ar,itog_lr,idomain,imtnode,iplbnode,ipltnode,irlbnode,
     *	imbnode)
c loop till/including end of pro-lith
      	nrowlr=nshiftlr; 
      	nshiftlr=0; nshiftap=0; nshiftlp=nshiftlp; nshiftar=0;
      	itog_ap=0; itog_lp=1; itog_lr=0; itog_ar=1;
      	iapon=0; ilpon=1; iaron=1; ilron=1;
      	call mk_array(0,0,iapon,1,nrowl+1,ilpon,1,nerowm,1,0,
     *	iaron,1,nrowlr,ilron,irabeg+1,npltop-1,nshiftap,nshiftlp,
     *  nshiftar,nshiftlr,
     *	index,indexap,indexlp,indexm,indexar,indexlr,itog_ap,itog_lp,
     *	itog_ar,itog_lr,idomain,imtnode,iplbnode,ipltnode,irlbnode,
     *	imbnode);
      	nrowar=nshiftar+1;
      else if(icase.eq.3) then	
c#############
c Case 3
c#############
c loop till the beg of retro-asthen      
      	nshiftlr=0; nshiftap=0; nshiftlp=0; nshiftar=0;
      	itog_ap=1; itog_lp=0; itog_lr=1; itog_ar=0;
      	iapon=1; ilpon=1; iaron=0; ilron=1;
      	call mk_array(1,nrowa,iapon,1,nrowl,ilpon,1,nerowm,0,0,
     *	iaron,1,0,ilron,nsing+1,irabeg,nshiftap,nshiftlp,nshiftar,
     *  nshiftlr,
     *	index,indexap,indexlp,indexm,indexar,indexlr,itog_ap,itog_lp,
     *	itog_ar,itog_lr,idomain,imtnode,iplbnode,ipltnode,irlbnode,
     *	imbnode);
c loop till end of pro-asthen
      	nrowlr=nshiftap;
      	nshiftlr=0; nshiftap=nshiftap; nshiftlp=0; nshiftar=0;
      	itog_ap=1; itog_lp=0; itog_lr=0; itog_ar=1;
      	iapon=1; ilpon=1; iaron=1; ilron=1;
      	call mk_array(1,nrowa,iapon,1,nrowl,ilpon,1,nerowm,1,0,
     *	iaron,1,nrowlr,ilron,irabeg+1,nplbase-1,nshiftap,nshiftlp,
     *  nshiftar,nshiftlr,
     *	index,indexap,indexlp,indexm,indexar,indexlr,itog_ap,itog_lp,
     *	itog_ar,itog_lr,idomain,imtnode,iplbnode,ipltnode,irlbnode,
     *	imbnode);
c loop till end of pro-lith
c     	temp # rows in retro-astthen
      	nartemp=nshiftar;
      	nshiftlr=0; nshiftap=0; nshiftlp=0; nshiftar=0;
      	itog_ap=0; itog_lp=1; itog_lr=0; itog_ar=1;
      	iapon=0; ilpon=1; iaron=1; ilron=1;
      	call mk_array(0,0,iapon,1,nrowl+1,ilpon,1,nerowm,1,nartemp,
     *	iaron,1,nrowlr,ilron,nplbase,npltop-1,nshiftap,nshiftlp,
     *  nshiftar,nshiftlr,
     *	index,indexap,indexlp,indexm,indexar,indexlr,itog_ap,itog_lp,
     *	itog_ar,itog_lr,idomain,imtnode,iplbnode,ipltnode,irlbnode,
     *	imbnode);
      	nrowar=nartemp+1+nshiftar;
      else if(icase.eq.1) then
      	print*,'####################################################'
      	print*,'####################################################'
      	print*,'CASE 1: this is not implemented. the pro-asthen     '
      	print*,'	ends where the retro-lith ends.  change the     '
      	print*,'	lith thickness or increase the mesh density     '
      	print*,'####################################################'
      else if(icase.eq.2) then
      	print*,'####################################################'
      	print*,'####################################################'
      	print*,'CASE 2: this is not implemented. the pro-asthen     '
      	print*,'	ends before the spoint. are you sure the        '
      	print*,'	orogen dimensions are reasonable                '
      	print*,'####################################################'
      else if(icase.eq.6) then
      	print*,'####################################################'
      	print*,'####################################################'
      	print*,'CASE 6: this is not implemented. see notes for     '
      	print*,'	details. change lith thickness or node density '
      	print*,'	to avoid this                                  '
      	print*,'####################################################'
      else if(icase.eq.4) then
      	print*,'####################################################'
      	print*,'####################################################'
      	print*,'CASE 4: this is not implemented. see notes for     '
      	print*,'	details. change lith thickness or node density '
      	print*,'	to avoid this                                  '
      	print*,'####################################################'
      endif	
c#############################################
c always end meshing from the end of the top 
c 	of the pro-lith to the model edge
c#############################################
      nshiftlr=0; nshiftap=0; nshiftlp=0; nshiftar=0;
      itog_ap=0; itog_lp=0; itog_lr=0; itog_ar=0;
      iapon=0; ilpon=0; iaron=1; ilron=1;
      call mk_array(0,0,iapon,0,0,ilpon,1,nerowm,1,nrowar,
     *iaron,1,nrowlr,ilron,npltop,ncol,nshiftap,nshiftlp,nshiftar,
     *nshiftlr,
     *index,indexap,indexlp,indexm,indexar,indexlr,itog_ap,itog_lp,
     *itog_ar,itog_lr,idomain,imtnode,iplbnode,ipltnode,irlbnode,
     *imbnode);

      index=0
      do j=1,ncol-1
      	do i=1,ntrow-1
      		index=index+1
      		ndomain(index)=domain(i+(j-1)*(ntrow-1))
      		index=index+1
      		ndomain(index)=domain(i+(j)*(ntrow-1))
      	end do
      end do	

c set final node in array of lith base boundary
      plitht_nodes(npltop)=plitht_nodes(npltop-1)+ntrow-1
      plithb_nodes(nplbase)=plithb_nodes(nplbase-1)+ntrow-1
      end


c########################################################
c combine nodes from differnt domains into one array
c########################################################
      subroutine mk_array(iapstart,iapstop,iapon,ilpstart,
     *ilpstop,ilpon,imstart,imstop,iarstart,iarstop,
     *iaron,ilrstart,ilrstop,ilron,icolstart,icolstop,
     *nshiftap,nshiftlp,nshiftar,nshiftlr,index,indexap,indexlp,
     *indexm,indexar,indexlr,itog_ap,itog_lp,itog_ar,itog_lr,idomain,
     *imtnode,iplbnode,ipltnode,irlbnode,imbnode)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

c loop over colms
      do icol=icolstart,icolstop
      	if(itog_ap.eq.1) then
      		nshiftap=nshiftap+1;
      	endif
      	if(itog_lp.eq.1) then
      		nshiftlp=nshiftlp+1;
      	endif
      	if(itog_ar.eq.1) then
      		nshiftar=nshiftar+1;
      	endif
      	if(itog_lr.eq.1) then
      		nshiftlr=nshiftlr+1;
      	endif
c pro-asthenosphere
      	if(iapon.eq.1) then
      		do irow=iapstart,iapstop-nshiftap
      			index=index+1;
      			indexap=indexap+1;
      			pos(index,1)=panode(indexap,1);
      			pos(index,2)=panode(indexap,2);
      			idomain=idomain+1
      			domain(idomain)=4
      		end do
      		indexap=indexap+1;
      		iplbnode=iplbnode+1
      		plithb_nodes(iplbnode)=index+1
      	endif	
c pro-lithosphere
      	if(ilpon.eq.1) then
      		do irow=ilpstart,ilpstop-nshiftlp
      			index=index+1;
      			indexlp=indexlp+1;
      			pos(index,1)=plnode(indexlp,1);
      			pos(index,2)=plnode(indexlp,2);
      			idomain=idomain+1
      			domain(idomain)=2
      		end do
      		indexlp=indexlp+1;
      		ipltnode=ipltnode+1
      		plitht_nodes(ipltnode)=index+1
      	endif
c retro-asthenosphere	
      	if(iaron.eq.1) then
      		do irow=iarstart,iarstop+nshiftar
      			index=index+1;
      			indexar=indexar+1;
      			pos(index,1)=ranode(indexar,1);
      			pos(index,2)=ranode(indexar,2);
      			idomain=idomain+1
      			domain(idomain)=5
      		end do
      		indexar=indexar+1;
      	endif
c retro-lith
      	if(ilron.eq.1) then
      		irlbnode=irlbnode+1
      		rlithb_nodes(irlbnode)=index+1
      		do irow=ilrstart,ilrstop+nshiftlr
      			index=index+1;
      			indexlr=indexlr+1;
      			pos(index,1)=rlnode(indexlr,1);
      			pos(index,2)=rlnode(indexlr,2);
      			idomain=idomain+1
      			domain(idomain)=3
      		end do
      		indexlr=indexlr+1;
      	endif
c mech
      	imbnode=imbnode+1
      	mechb_nodes(imbnode)=index+1
      	do irow=imstart,imstop
      		index=index+1;
      		indexm=indexm+1;
      		pos(index,1)=cmnode(indexm,1);
      		pos(index,2)=cmnode(indexm,2);
      			idomain=idomain+1
      			domain(idomain)=1
      	end do
      	idomain=idomain-1
      	imtnode=imtnode+1
      	mecht_nodes(imtnode)=index
      end do
      end

c########################################################
c make nodes in pro-asthenosphere
c########################################################
      subroutine mk_panodes(nplbase,nsing,nparow,base)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c det number of nodes in the domain      
      nparow=nplbase-nsing+1
      index=0
      npa=nsing*nparow
      do i=1,nparow-1
      	npa=npa+i
      end do
      allocate(panode(npa,2))
      do i=1,nsing
      	dy=(yplbase(i)-base)/dble(nparow-1)
      	do j=1,nparow
      		index=index+1
      		panode(index,1)=xmbase(i)
      		panode(index,2)=base+dble(j-1)*dy
      	end do
      end do	
      icount=0
      do i=nsing+1,nplbase-1
      	icount=icount+1
      	dy=(yplbase(i)-base)/dble(nparow-icount-1)
      	do j=1,nparow-icount
      		index=index+1
      		panode(index,1)=xmbase(i)
      		panode(index,2)=base+dble(j-1)*dy
      	end do
      end do
      index=index+1
      panode(index,1)=xmbase(nplbase)
      panode(index,2)=yplbase(nplbase)
c      do i=1,npa
c      	print*,panode(i,1),panode(i,2)
c      end do	
      end

c########################################################
c make nodes in retro asthenosphere
c########################################################
      subroutine mk_ranodes(nrlrow,nsing,npltop,nrarow,ncol,
     *base,irlbeg,irabeg)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c det number of rows and number of nodes in region      
      irabeg=nrlrow+nsing-1
      nrarow=npltop-irlbeg+1
      nra=(ncol-npltop+1)*nrarow
      do i=nrarow-1,1,-1
      	nra=nra+i
      end do	
      allocate(ranode(nra,2))
      ranode(1,1)=xmbase(irabeg)
      ranode(1,2)=ypltop(irabeg)
      ranode(2,1)=xmbase(irabeg+1)
      ranode(2,2)=ypltop(irabeg+1)
      ranode(3,1)=xmbase(irabeg+1)
      ranode(3,2)=yrlbase(2)
      index=3
      do i=2,npltop-irabeg
      	icount=i+1
      	dy=(yrlbase(i+1)-ypltop(irabeg+i))/dble(icount-1)
      	do j=1,icount
      		index=index+1
            ranode(index,1)=0.0
      		ranode(index,1)=xmbase(i+irabeg)
      		ranode(index,2)=ypltop(i+irabeg)+dy*dble(j-1)
      	end do
      end do	
      do i=2,ncol-npltop+1
      	dy=(yrlbase(npltop-irabeg+i)-base)/dble(nrarow-1)
      	do j=1,nrarow
      		index=index+1
      		ranode(index,1)=xrlbase(i+npltop-irabeg)
      		ranode(index,2)=base+dy*dble(j-1)
      	end do
      end do	
c      do i=1,nra
c      	print*,ranode(i,1),ranode(i,2)
c      end do	
      end

c########################################################
c make node in the retro lithosphere
c########################################################
      subroutine mk_rlnodes(nsing,nrlrow,irlbeg,ncol)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c det number of rows and number of nodes in region      
      nrlrow=irlbeg-nsing+1
      nrl=(ncol-irlbeg+1)*nrlrow
      do i=nrlrow-1,1,-1
      	nrl=nrl+i
      end do	
      allocate(rlnode(nrl,2))
      rlnode(1,1)=xmbase(nsing)
      rlnode(1,2)=ypltop(nsing)
      rlnode(2,1)=xmbase(nsing+1)
      rlnode(2,2)=ypltop(nsing+1)
      rlnode(3,1)=xmbase(nsing+1)
      rlnode(3,2)=ymbase(nsing+1)
      index=3
      do i=2,irlbeg-nsing
      	icount=i+1
      	dy=(ymbase(nsing+i)-ypltop(nsing+i))/dble(icount-1)
      	do j=1,icount
      		index=index+1
      		rlnode(index,1)=xmbase(i+nsing)
      		rlnode(index,2)=ypltop(i+nsing)+dble(j-1)*dy
      	end do
      end do
      do i=1,ncol-irlbeg
      	dy=(ymbase(irlbeg+i)-yrlbase(i+1))/dble(nrlrow-1)
      		do j=1,nrlrow
      			index=index+1
      			rlnode(index,1)=xmbase(i+irlbeg)
      			rlnode(index,2)=yrlbase(i+1)+dy*dble(j-1)
      		end do
      	end do 	
c      	do i=1,nrl
c      		print*,rlnode(i,1),rlnode(i,2)
c      	end do	
      end

c########################################################
c make nodes in pro-lithosphere
c########################################################
      subroutine mk_plnodes(npltop,nplbase,nplrow,base)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      base=ypltop(npltop)
c det # no rows in pro-lithosphere
      do i=1,npltop
      	if(xpl(i).ge.xpl(nplbase)) then
      		istart=i
      		exit
      	endif
      end do
      nplrow=npltop-istart+1
c det size of node array
      npl=nplrow*nplbase
      do i=nplrow-1,1,-1
      	npl=npl+i
      end do
      allocate(plnode(npl,2))
c det node locations
      index=0
      do i=1,nplbase
      	dy=(ypltop(i)-yplbase(i))/dble(nplrow-1)
      		do j=1,nplrow
      			index=index+1
      			plnode(index,1)=xmbase(i)
      			plnode(index,2)=yplbase(i)+dble(j-1)*dy
      		end do
      	end do	
      	icount=0
      	do i=nplbase+1,npltop-1
      		icount=icount+1
      		dy=(ypltop(i)-base)/dble(nplrow-1-icount)
      		do j=1,nplrow-icount
      			index=index+1
      			plnode(index,1)=xmbase(i)
      			plnode(index,2)=base+dble(j-1)*dy
      		end do
      	end do
      	index=index+1
      	plnode(index,1)=xmbase(npltop)
      	plnode(index,2)=base
c      	do i=1,npl
c      		print*,plnode(i,1),plnode(i,2)
c      	end do	

      end
      
c########################################################
c make nodes in the mech model 
c########################################################
      subroutine mk_cmnodes(ncol,nerowm,iblay,iblayt,iflgblt,iflgbl,
     *tblayt,tblay)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      allocate(cmnode(ncol*nerowm,2))
      index=0

c check that defined thicknesses of boundary layers do not exceed
c	the pro side crustal thickness
      if(tblayt+tblay.ge.ymtop(1)-ymbase(1)) then
      	print*,'######################################'
      	print*,'## ERROR: prescribed thickness of   ##'
      	print*,'##	boundary layers is greater than ##'
      	print*,'##	the crustal thickness.          ##'
      	print*,'##  tblay+tblayt=',tblay+tblayt
      	print*,'##  crustal thickness=',ymtop(1)-ymbase(1)
      	print*,'######################################'
      endif	
      
c define thickness for boundary layers so that their thickness remains 
c	constant
c     user defined thickness of bounadry layers
c     top layer      
      if(iflgblt.eq.1) then
      	dyblayt=tblayt/dble(iblayt)
      endif
c     bottom layer      
      if(iflgbl.eq.1) then
		dyblay=tblay/dble(iblay)
      endif	
c     mixed: user defined and automatic even spacing
c     top layer
      if(iflgbl.eq.1.and.iflgblt.eq.0) then
      	dyblayt=(ymtop(1)-ymbase(1)-tblay)/dble(nerowm-1-iblay)
      	tblayt=dble(iblayt)*dyblayt
      endif
c     bottom layer
      if(iflgblt.eq.1.and.iflgbl.eq.0) then
      	dyblay=(ymtop(1)-ymbase(1)-tblayt)/dble(nerowm-1-iblayt)
      	tblay=dble(iblay)*dyblay
      endif	
c     automatic even spacing of both layers      
      if(iflgbl.eq.0.and.iflgblt.eq.0) then
      	dyblay=(ymtop(1)-ymbase(1))/dble(nerowm-1)
      	dyblayt=dyblay
      	tblayt=dble(iblayt)*dyblayt
      	tblay=dble(iblay)*dyblay
      endif	

c make sure that no boundary layers -> boundary layer thickness = 0
      if(iblayt.eq.0) then
      	tblayt=0.0
      endif	
      if(iblay.eq.0) then
      	tblay=0.0
      endif	

c loop over all colms      
      do i=1,ncol
      	dy=(ymtop(i)-ymbase(i)-tblayt-tblay)
     *	/dble(nerowm-1-iblay-iblayt)
c	const thickness for lower boundary layers
      	do j=1,iblay
      		index=index+1
      		cmnode(index,1)=xmbase(i)
      		cmnode(index,2)=ymbase(i)+dble(j-1)*dyblay
      	end do	
c	fanning thickness for all other layers
      	do j=1,nerowm-iblayt-iblay
      		index=index+1
      		cmnode(index,1)=xmbase(i)
      		cmnode(index,2)=ymbase(i)+tblay+dble(j-1)*dy
      	end do
c	constant thickness for top boundary layers
      	base=cmnode(index,2)
      	do j=1,iblayt
      		index=index+1
      		cmnode(index,1)=xmbase(i)
      		cmnode(index,2)=base+dble(j)*dyblayt
      	end do	
      end do	
      end

c#########################################################
c make the boundaris of the lithosphere/mesh/asthenosphere
c 	on the pro side and lithosphere/asthenosphere on the
c	retro side
c#########################################################
      subroutine lith_bndry(nsing1,nsing,plthick,athick,ncol,npltop,
     *plvel,nplbase,rlthick,irlbeg,nrlbase)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real(kind=8),allocatable::ylbtemp(:)
c determine slope at the end of the loaded sub. lithosphere
c	used in continuing sub plate down to model boundary
      slp=(yltemp(nsing1)-yltemp(nsing1-1))/(xltemp(nsing1)
     *-xltemp(nsing1-1))
      print*,' Angle at end of plate (deg): ',atan(slp)*180.0/3.1416

c det number of nodes in top of sub lithos      
      base=ymbase(1)-plthick-athick
c if sub entension ends above model base      
      if(yltemp(nsing1).ge.base) then
      	icount=0
      	do i=nsing+nsing1-1,ncol
      		icount=icount+1
      		dx=abs(xmbase(nsing+nsing1-1)-xmbase(i))
      		ytest=yltemp(nsing1)+dx*slp
      		if(ytest.lt.base) then
				dif1=abs(ytest-base)
				dif2=abs(yltemp(nsing1)+(abs(xmbase(nsing+nsing1-1)
     *			-xmbase(i-1)))*slp-base)
     			if(dif1.ge.dif2) then
     				npltop=nsing1+icount-1
     				exit
     			else
     				npltop=nsing1+icount-2
     				exit
     			endif
     		end if	
      	end do	
c	make array for sub lithos top
      	allocate(ypltop(nsing+npltop-1))
      	allocate(xpl(nsing+npltop-1))
      	do i=1,nsing1
      		ypltop(i-1+nsing)=yltemp(i)
      		xpl(i-1+nsing)=xltemp(i)
      	end do
      	index=0
      	do i=nsing1+1,npltop-1
      		xpl(i-1+nsing)=xmbase(nsing+i-1)
      		dx=abs(xpl(i-1+nsing)-xmbase(nsing-1+nsing1))
      		ypltop(i-1+nsing)=yltemp(nsing1)+slp*dx
      	end do
      	ypltop(npltop+nsing-1)=base
      	xpl(npltop+nsing-1)=xmbase(npltop+nsing-1)
      	do i=1,nsing-1
      		xpl(i)=xmbase(i)
      		ypltop(i)=ymbase(i)
      	end do	
      	npltop=npltop+nsing-1
c if sub extension ends below model base
      else if (yltemp(nsing1).lt.base) then
     	icount=0
     	do i=1,nsing1
     		if(yltemp(i).lt.base) then
     			dif1=abs(base-yltemp(i-1))
     			dif2=abs(base-yltemp(i))
     			if(dif1.lt.dif2) then
     				npltop=i-1
     			else
     				npltop=i
     			endif	
     			exit
     		endif
     	end do
c	make arrays for sub lithos top     	
     	allocate(ypltop(nsing+npltop-1))
     	allocate(xpl(nsing+npltop-1))
     	do i=1,nsing-1
     		xpl(i)=xmbase(i)
     		ypltop(i)=ymbase(i)
     	end do
     	do i=1,npltop-2
     		xpl(i-1+nsing)=xltemp(i)
     		ypltop(i-1+nsing)=yltemp(i)
     	end do	
     	xpl(nsing+npltop-1)=xmbase(nsing+npltop-1)
     	ypltop(nsing+npltop-1)=base
     	npltop=npltop+nsing-1
      endif 	

c make array for lithos/asthenos boundary on pro side      
c allocate temp storage for the array that is as long as the array 
c	for the top.  when final length is determined, store in perm.
c	array
      allocate(ylbtemp(npltop))
c ref flux value for conserving mass in sub. plate
      qx=plvel*plthick
c thickness at lhs
      dx=xpl(2)-xpl(1)
      dy=ypltop(2)-ypltop(1)
      hyplen=(dx**2+dy**2)**0.5
      ylbtemp(1)=-qx/plvel*hyplen/dx+ypltop(1)
      do i=2,npltop-1
      	dx=xpl(i+1)-xpl(i-1)
      	dy=ypltop(i+1)-ypltop(i-1)
      	hyplen=(dx**2+dy**2)**0.5
      	test=-qx/plvel*hyplen/dx+ypltop(i)
      	if(test.lt.base) then
      		dif1=abs(test-base)
      		dif2=abs(ylbtemp(i-1)-base)
      		if(dif1.lt.dif2) then
      			ylbtemp(i)=base
      			nplbase=i
      		else
      			ylbtemp(i-1)=base
      			nplbase=i-1
      		endif
      		exit
      	endif
      	ylbtemp(i)=test
      end do	

      allocate(yplbase(nplbase))
      do i=1,nplbase
      	yplbase(i)=ylbtemp(i)
      end do
      deallocate(ylbtemp)

c make array for lithos/asthenos boundary on retro side
c determine hit point of boundary of sub lithos slab
      difmin=10.0e5
      imin=0
      do i=nsing+1,npltop
      	dif=abs(ymbase(i)-rlthick-ypltop(i))
      	if(dif.lt.difmin) then
      		difmin=dif
      		imin=i
      	endif
      end do	
c x node where base of retro lith touches sub lith
      irlbeg=imin
      nrlbase=ncol-irlbeg+1
      index=1
      allocate(yrlbase(nrlbase))
      allocate(xrlbase(nrlbase))
      yrlbase(1)=ypltop(irlbeg)
      xrlbase(1)=xmbase(irlbeg)
      do i=irlbeg+1,ncol
      	index=index+1
      	yrlbase(index)=ymbase(i)-rlthick
      	xrlbase(index)=xmbase(i)
      end do	
      end

c#########################################################
C calclate the additional deflection of the coupled plates
c	from the overlying load of water
c#########################################################
      subroutine deflectw(wdepth,xp1,xp2,slen1,slen2,wtoler,fnode1
     *,fnode2,plam1,plam2,fk,rhof,np1,np2,npad,g,nsing1,ctoler,
     *dyinit1,dyinit2,yp2,yp1,wheight)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 xp1(*),xp2(*),yp1(*),yp2(*),fnode1(*),fnode2(*),slen1(*),
     *slen2(*),dyinit1(*),dyinit2(*)
      real(kind=8),allocatable::dloc1(:),dloc2(:),yp1pre(:),yp2pre(:)

      ychange=100.0*wtoler
      icount=0
      allocate(dloc1(np1))
      allocate(dloc2(np2))
      allocate(yp1pre(np1))
      allocate(yp2pre(np2))
      do i=1,np1
      	dloc1(i)=0.0
      	yp1pre(i)=yp1(i)
      	fnode1(i)=0.0
      end do
      do i=1,np2
      	dloc2(i)=0.0
      	yp2pre(i)=yp2(i)
      	fnode2(i)=0.0
      end do	
      do while(ychange.gt.wtoler)
      	icount=icount+1
c define water height
      	wheight=-yp1(np1-npad)+dyinit1(np1-npad)+wdepth
c calculate loading due to water
      	if(icount.eq.1) then
c	first step      	
c	plate 1      		
c			local water depth and resulting force
      		do i=nsing1,np1
      			dloc1(i)=wheight-(-yp1(i)+dyinit1(i))
      			if(dloc1(i).gt.0.0) then
      				fnode1(i)=slen1(i)*g*rhof*dloc1(i)
      			else
      				fnode1(i)=0.0
      				dloc1(i)=0.0
      			endif
      		end do
c	plate 2
c			local water depth and resulting force
      		do i=1,np2
      			dloc2(i)=wheight-(-yp2(i)+dyinit2(i))
      			if(dloc2(i).gt.0.0) then
      				fnode2(i)=slen2(i)*g*rhof*dloc2(i)
      			else
      				fnode2(i)=0.0
      				dloc2(i)=0.0
      			endif
      		end do	
      	else
c all other steps
c	plate 1      		
      		do i=nsing1,np1
      			dtemp=wheight-(-yp1(i)+dyinit1(i))
c				change in water depth
      			deltad=dtemp-dloc1(i)
      			if(dtemp.le.0.0) then
      				dtemp=0.0
      				if(dloc1(i).le.0.0) then
      					fnode1(i)=0.0
      				else
      					fnode1(i)=slen1(i)*g*rhof*(-dloc1(i))
      				endif
      			else
      				fnode1(i)=slen1(i)*g*rhof*deltad
      			endif
      			dloc1(i)=dtemp
      		end	do
c	plate 2      		
      		do i=1,np2
      			dtemp=wheight-(-yp2(i)+dyinit2(i))
c				change in water depth
      			deltad=dtemp-dloc2(i)
      			if(dtemp.le.0.0) then
      				dtemp=0.0
      				if(dloc2(i).le.0.0) then
      					fnode2(i)=0.0
      				else
      					fnode2(i)=slen2(i)*g*rhof*(-dloc2(i))
      				endif
      			else
      				fnode2(i)=slen2(i)*g*rhof*deltad
      			endif
      			dloc2(i)=dtemp
      		end	do
      	endif	
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
      	call deflect2(np1,plam1,fk,xp1,0.0,0.0,amom1,ashear1,yp1)
c     	plate 2
      	call deflect2(np2,plam2,fk,xp2,0.0,0.0,amom2,ashear2,yp2)
c calculate the coupling load
      	ido_again=1
      	jcount=0
      	do while(ido_again==1) 
      		jcount=jcount+1
			call couple(yp1,yp2,nsing,plam1,plam2,fk,xp1,xp2,np1,np2,
     *		nsing1)
      		if(abs(yp1(nsing1)-yp2(1)).le.ctoler) then
      			ido_again=0
      			print*,'diff. at s-point: ',abs(yp1(nsing1)-yp2(1))
      		else if(jcount.gt.1000) then
      			ido_again=0
      			print*,'########################################'
      			print*,'## coupling iteration exceeded 1000   ##'
      			print*,'##      inside water loop             ##'
      			print*,'########################################'
      			call profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
      			stop
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
      	print*,'change in base due to water loading: ',dif
      	ychange=dif
      	if(icount.gt.1000) then
      		print*,'########################################'
      		print*,'## water depth iteration exceeded 1000##'
      		print*,'########################################'
      		call profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
      		stop
      	endif	
      end do
      deallocate(dloc1)
      deallocate(dloc2)
      deallocate(yp1pre)
      deallocate(yp2pre)
      end


c########################################################
C calculate the deflection of two semi-infinite plates 
c	coupled together at the s-point from a distributed
c	load as stored in fnode
c########################################################
      subroutine deflect(np1,np2,xp1,xp2,yp1,yp2,fnode1,fnode2
     *,plam1,plam2,fk,ctoler,nsing1,sload,smomen,xbase,nsing,npad,
     *dyc)

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 xp1(*),xp2(*),yp1(*),yp2(*),fnode1(*),fnode2(*),xbase(*)

c calculate deflection everywhere and moment,shear force at the 
c desired break point for two infinite plates
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
      call deflect2(np1,plam1,fk,xp1,sload,smomen,amom1,ashear1,yp1)
c     plate 2
      call deflect2(np2,plam2,fk,xp2,0.0,0.0,amom2,ashear2,yp2)
c      call profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
c      stop
c calculate the coupling load and coupled position of plates
      ido_again=1
      icount=0
      do while(ido_again==1) 
      	icount=icount+1
		call couple(yp1,yp2,nsing,plam1,plam2,fk,xp1,xp2,np1,
     *	np2,nsing1)
      	if(abs(yp1(nsing1)-yp2(1)).le.ctoler) then
      		ido_again=0
      		print*,'diff. at s-point: ',abs(yp1(nsing1)-yp2(1))
      	else if(icount.gt.1000) then
      		ido_again=0
      		print*,'########################################'
      		print*,'## coupling iteration exceeded 1000   ##'
      		print*,'##    in first calc                   ##'
      		print*,'########################################'
      		call profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
      		stop
      	endif	
      end do
c allow for shifting of the coupling point
      if(abs(dyc).gt.0.0) then
      	ido_again=1
      	icount=0
      	ycfinal=yp2(1)-dyc
		do while(ido_again==1) 
      		icount=icount+1
			call shift_couple(yp1,yp2,nsing,plam1,plam2,fk,xp1,xp2,np1,
     *		np2,nsing1,ycfinal)
      		if(abs(yp1(nsing1)-yp2(1)).le.ctoler) then
      			ido_again=0
      		else if(icount.gt.1000) then
      			ido_again=0
      			print*,'########################################'
      			print*,'## coupling iteration exceeded 1000   ##'
      			print*,'##    in couple point shift           ##'
      			print*,'########################################'
      			call profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
      			stop
      		endif	
      	end do
      endif
      end

C########################################################
c calculate the new position of the plates after the 
c	coupling point has been shifted
C########################################################
      subroutine shift_couple(yp1,yp2,nsing,plam1,plam2,fk,xp1,xp2,
     *np1,np2,nsing1,ycfinal)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 yp1(*),yp2(*),xp1(*),xp2(*)

c distances to shift both plates at the coupling point.
      G0p1=yp1(nsing1)-ycfinal
      G0p2=-yp2(1)+ycfinal
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
      fcouplep1=G0p1/(G1+G6+G7)
      fcouplep2=G0p2*fk/(2.0*plam2)
c calculate deflection from coupling load and shift in couple point
      do i=1,np1
      	yp1(i)=yp1(i)
     *	-(2.0*fcouplep1*plam1/fk*exp(-plam1*xp1(i))*cos(plam1*xp1(i)))
      end do
      do i=1,np2
      	yp2(i)=yp2(i)
     *	+(2.0*fcouplep2*plam2/fk*exp(-plam2*xp2(i))*cos(plam2*xp2(i)))
      end do
      end

C########################################################
c calculate the plate coupling load and couple the plates
c---Couples the plates so that both plates experience
c	equal and opposite forces to couple them
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
        yp(i)=ypo+ymo+ysload+ysmom+yp(i)
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

c########################################################
c calculate the force from the thickness of the mech model
c	for calculating the flexure
c########################################################
      subroutine calc_force_p2(slen,xp,nerowm,rhoc,np,fnode,dyinit,
     *npad,iplate,rhoavt,nsing,g,ncol,inflag,dencol)      

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 slen(*),xp(*),rhoc(*),fnode(*),dyinit(*),dencol(*)

      ifstrow=(nerowm-1)*(nsing-1)
      ilstrow=(ncol-2)*(nerowm-2)

c first node
      slen(1)=(xp(2)-xp(1))/2.0
c	average density of overlying colm, need for variable density models
c		!!!!!ONLY WORKS WHEN Y SPACING OF ELEMENTS IS CONSTANT!!!!!!
c		ALSO: the thickness of the plate 1 extension past nsing must be zero
c       prob. want to check this before implementing varying densities
      rhosum=0.0
      do i=1,nerowm-1
      	rhosum=rhosum+rhoc(i+ifstrow)
      end do
      rhoav=rhosum/dble(nerowm-1)
      if(inflag.eq.1) rhoav=dencol(nsing)
      fnode(1)=slen(1)*g*rhoav*dyinit(1)
c last node
      slen(np)=(xp(np)-xp(np-1))/2.0
      rhosum=0.0
      if(npad.eq.0) then
      	do i=1,nerowm-1
      		rhosum=rhosum+rhoc(i+ilstrow)
      	end do
      	rhoav=rhosum/dble(nerowm-1)
      else
      	rhoav=rhoavt
      endif	
      if(inflag.eq.1) rhoav=dencol(ncol)
      fnode(np)=slen(np)*g*rhoav*dyinit(np)
c all other nodes
      do icol=2,np-1
      	slen(icol)=(xp(icol+1)-xp(icol-1))/2.0
      	rhosum=0.0
c       add catch for padded edges of model where density is not defined
      	if(icol.ge.np-npad) then
     		rhoav=rhoavt
     		if(inflag.eq.1) rhoav=dencol(ncol)
      	else	
     		icol2=icol+nsing-1
     		dxr=xp(icol+1)-xp(icol)
      		dxl=xp(icol)-xp(icol-1)
      		sl=dxl/(dxl+dxr)
      		sr=dxr/(dxl+dxr)
      		do j=1,nerowm-1
      			rhosum=rhosum+sl*rhoc(j+(icol2-2)*(nerowm-1))
     *				+sr*rhoc(j+(icol2-1)*(nerowm-1))
      		end do
      		rhoav=rhosum/dble((nerowm-1)*1)
      		if(inflag.eq.1) rhoav=dencol(icol2)
      	end	if
      	fnode(icol)=slen(icol)*g*rhoav*dyinit(icol)
      end do	
      end

c########################################################
c calculate the force from the thickness of the mech model
c	for calculating the flexure
c########################################################
      subroutine calc_force_p1(slen,xp,nerowm,rhoc,np,fnode,dyinit,
     *npad,iplate,rhoavt,nsing,g,nsing1,inflag,dencol)      

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*8 slen(*),xp(*),rhoc(*),fnode(*),dyinit(*),dencol(*)

      ifstrow=(nerowm-1)*(np-1)
      ilstrow=0

c	average density of overlying colm, need for variable density models
c		!!!!!ONLY WORKS WHEN Y SPACING OF ELEMENTS IS CONSTANT!!!!!!
c		ALSO: the thickness of the plate 1 extension past nsing must be zero
c       prob. want to check this before implementing varying densities


c first node
c	if there is an extended plate (sub plate), set forces to zero
      slen(1)=(xp(2)-xp(1))/2.0
      fnode(1)=0.0
      do i=2,nsing1-1
      	slen(i)=(xp(i+1)-xp(i-1))/2.0
      	fnode(i)=0.0
      end do	
c calculate force at spoint as if it was the first colm      
      slen(nsing1)=(xp(nsing1+1)-xp(nsing1))/2.0
      rhosum=0.0
      do i=1,nerowm-1
      	rhosum=rhosum+rhoc((nsing-2)*(nerowm-1)+i)
      end do
      rhoav=rhosum/dble(nerowm-1)
c     catch for reading in colms of ave density     
      if(inflag.eq.1) rhoav=dencol(nsing)
      fnode(nsing1)=slen(nsing1)*g*rhoav*dyinit(nsing1)
c last node
      slen(np)=(xp(np)-xp(np-1))/2.0
      rhosum=0.0
      if(npad.eq.0) then
      	do i=1,nerowm-1
      		rhosum=rhosum+rhoc(i)
      	end do
      	rhoav=rhosum/dble(nerowm-1)
      else
      	rhoav=rhoavt
      endif	
      if(inflag.eq.1) rhoav=dencol(1)
      fnode(np)=slen(np)*g*rhoav*dyinit(np)

c all other nodes
      index=0
      do icol=nsing1+1,np-1
      	slen(icol)=(xp(icol+1)-xp(icol-1))/2.0
      	rhosum=0.0
c       add catch for padded edges of model where density is not defined
      	if(icol.ge.np-npad) then
     		rhoav=rhoavt
     		if(inflag.eq.1) rhoav=dencol(1)
      	else	
      		index=index+1
      		icol2=nsing-index
      		dxl=xp(icol+1)-xp(icol)
      		dxr=xp(icol)-xp(icol-1)
      		sl=dxl/(dxl+dxr)
      		sr=dxr/(dxl+dxr)
      		do j=1,nerowm-1
      			rhosum=rhosum+sl*rhoc((icol2-2)*(nerowm-1)+j)
     *				+sr*rhoc((icol2-1)*(nerowm-1)+j)
      		end do
      		rhoav=rhosum/dble((nerowm-1)*1)
      		if(inflag.eq.1) rhoav=dencol(icol2)
      	end	if
      	fnode(icol)=slen(icol)*g*rhoav*dyinit(icol)
      end do	

      end


c##########################################################################
c read in input file
c##########################################################################
      subroutine input(ncol,nerowm,nstype,sing,
     *pmthick,plthick,athick,rlthick,ypmbase,yrmbase,wdepth,wtoler,
     *npad,xpad,prigp,rrigp,prigi,rrigi,sload,smomen,xadd,ctoler,
     *plvel,upvel,iunflag,iunbeg,xunbeg,vrig,beta,epsinv,
     *rhof,rhom,numvebn,numpbn,numsid,numvtbn,ntst,intout,intoutl,
     *delt,minint,maxint,npass,toler,erosl,erosr,peros,rpow,ntt2,
     *deltt2,iso,ntmchg,plscale,rlscale,blscale,dfact,slpmax,tmax,
     *numvetbn,ioutflag,inflag,dyc,linflag,sdip,ipflag,itrench,
     *iplasflg,iblay,iblayt,isedl,isedr,iexflg,ibasflg,nbastary,
     *nbastind,intmrkb,ipkfill,ibasfill,sedmax,iflgcl,agecl,iflgblt,
     *iflgbl,tblayt,tblay,noutput)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      open(2,file='meshin_oly',position='rewind')

c style of output (1=all, 0= just flexural profiles       
      read(2,106)dummy
      read(2,101)ioutflag
c allow input files for xpositions, thickness and densities of mech model
c	(1=read input files, 0= just use meshin)
      read(2,106)dummy
      read(2,101)inflag
      if(inflag.eq.1) open(9,file='../data/flex_data',position='rewind')
c number of colms in model and lagrangian mesh style
      read(2,106)dummy
      read(2,101)ncol
c number of rows in mech      
      read(2,106)dummy
      read(2,106)dummy
      read(2,101)nerowm
c lagrangian mesh parameters      
c 	(extent past pro side, extent past retro side, extent past base, 
c 	node density compred to eulerian mesh)
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      read(2,103)plscale,rlscale,blscale,dfact
c s-point location
      read(2,106)dummy
      read(2,104)nstype,sing
c model thicknesses on the pro and retro side
      read(2,106)dummy
      read(2,103)pmthick,plthick,athick
      read(2,106)dummy
      read(2,106)dummy
      read(2,103)rlthick
c relative dif in initial mech. base height for pro and reto
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      read(2,103)ypmbase,yrmbase
c water depth and tolerance for change in water depth
      read(2,106)dummy
      read(2,106)dummy
      read(2,103)wdepth,wtoler
c x padding for aprox. an infinite plate
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      read(2,104)npad,xpad
c horizontal spacing
      read(2,106)dummy
      read(2,106)dummy
      read(2,101)ninc
      read(2,106)dummy
      allocate(xbase(ncol+npad*2))
      allocate(dyinit(ncol+npad*2))
      allocate(dencol(ncol))
c     make model array      
      index=0
      do i=1,ninc
		read(2,104)nnodes,xstrt,xstp
		do j=1,nnodes
			index=index+1
			xbase(index+npad)=xstrt+(xstp-xstrt)/float(nnodes-1)*(j-1)
		end do	
      end do
      if(index.ne.ncol) then
      	print*,'###########################################'
      	print*,'ERROR:'
      	print*,'Number of colms != num of nodes in x array'
      	print*,'###########################################'
      	stop
      endif
c     add on padding beyond model edges
      do i=1,npad
      	dx=abs(xbase(1+npad)-xpad)/dble(npad)*dble(i)
      	xbase(npad+1-i)=xbase(npad+1)-dx
      	xbase(npad+ncol+i)=xbase(ncol+npad)+dx
      end do	
c dimensions used when making model arrays from obs data
c	ie, when the profiles will be read in
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      read(2,103)t1,t2
      if(inflag.eq.1) then
      	if(t2-t1.ne.xbase(npad+ncol))then
      		print*,'##############################################'
      		print*,'ERROR: the width of model array',xbase(npad+ncol)
      		print*,'       does not equal the width'
      		print*,'       defined by coastline ref',t2-t1
      		print*,'##############################################'
      	endif
      endif	

c deviation in mechanical thickness from ref thickness defined above      
      read(2,106)dummy
      read(2,106)dummy
      read(2,101)ninc
      read(2,106)dummy
      do i=1,ncol+npad*2
      	dyinit(i)=pmthick
      end do	
      do j=1,ninc
      	read(2,104)ntyp,beg,slp
      	if(ntyp.eq.0) then
			ibeg=int(beg)
		else if(ntyp.eq.1) then
			do i=1+npad,ncol+npad
      			if(xbase(i).gt.beg) then
      				if(abs(beg-xbase(i)).ge.abs(beg-xbase(i-1))) then
      					ibeg=i-1-npad
      					exit
      				else
      					ibeg=i-npad
	    	  			exit
	    	  		endif	
      			else if(xbase(i).eq.beg) then
      				ibeg=i-npad
      				exit
      			endif
      		end do	
		else
			print*,'###########################################'
			print*,'ERROR: change in thickness must be defined '
			print*,'	on either a node or x position with the'
			print*,' 	flag set as 0 (node) or 1 (pos)        '
			print*,'###########################################'
			stop
		endif	
		do i=ibeg+1+npad,ncol+npad*2
			dyinit(i)=(xbase(i)-xbase(i-1))*slp+dyinit(i-1)
		end do
      end do
c if reading in x array and thickness from file, redo the above
      if(inflag.eq.1.) then
      	read(9,101)ncol2,nsing2
      	if(ncol2.ne.ncol.or.nsing2.ne.int(sing)) then
      		print*,'###############################################'
      		print*,'ERROR:'
      		print*,'ncol or nsing from flex_data dont match meshin'
      		print*,'###############################################'
      	endif
c		read in data      	
      	do icol=1,ncol
      		read(9,239)xbase(npad+icol),dyinit(npad+icol),dencol(icol)
      		xbase(icol+npad)=xbase(icol+npad)*1000.0
      		dyinit(icol+npad)=dyinit(icol+npad)*1000.0
      	end do	
c       shift x coords so they start at zero
        xshift=xbase(npad+1)
        do icol=1,ncol
        	xbase(npad+icol)=xbase(npad+icol)-xshift
        end do	
c 		set thickness in padded region to the thickness at the model edge
     	do i=1,npad
     		dyinit(i)=dyinit(npad+1)
     		dyinit(npad+ncol+i)=dyinit(npad+ncol)
     	end do	
      endif	
      
c type of isostatic comp.
      read(2,106)dummy
      read(2,101)iso
c initial geometry flag
      read(2,106)dummy
      read(2,101)ipflag
c slab dip for prescribed profile
      read(2,106)dummy
      read(2,106)dummy
      read(2,103)sdip
c begining of circular arc on pro side
      read(2,106)dummy
      read(2,106)dummy
      read(2,104)ntyp,beg
      if(ntyp.eq.0) then
		ibeg=int(beg)
	  else if(ntyp.eq.1) then
		do i=1+npad,ncol+npad
      		if(xbase(i).gt.beg) then
      			if(abs(beg-xbase(i)).ge.abs(beg-xbase(i-1))) then
      				ibeg=i-1-npad
      				exit
      			else
      				ibeg=i-npad
	    	  		exit
	    	  	endif	
      		else if(xbase(i).eq.beg) then
      			ibeg=i-npad
      			exit
      		endif
      	end do	
      endif
      itrench=ibeg
c flexural rigidity for plate profile calc
      read(2,106)dummy
      read(2,103)prigp,rrigp
c flexural rigidity for isotatic calc
      read(2,106)dummy
      read(2,103)prigi,rrigi
c subduction load
      read(2,106)dummy
      read(2,103)sload
c subduction moment
      read(2,106)dummy
      read(2,103)smomen
c shift in coupling point
      read(2,106)dummy
      read(2,103)dyc
c length of extra pro-plate for sub load     
      do i=1,5
      	read(2,106)dummy
      end do	
      read(2,103)xadd
c extension flag (=1 don't inlcude extension for plasti input)      
      do i=1,4
      	read(2,106)dummy
      end do
      read(2,101)iexflg
c tolerance for plate coupling position
      read(2,106)dummy
      read(2,105)ctoler
c convergence and undeplating velocity
      read(2,106)dummy
      read(2,103)plvel,upvel
c underplating flag and position
      read(2,106)dummy
      read(2,106)dummy
      read(2,107)iunflag,iunbeg,xunbeg
c find node for begining of up
      if(iunflag.eq.2) then
      	call find_node(iunbeg,1,xunbeg,ncol,npad,xuse)
      	print*,'beg underplating at (node, xposition)'
      	print*,iunbeg,xuse
      else
      	print*,'beg underplating at (node, xposition)'
      	print*,iunbeg,xbase(iunbeg+npad)
      endif	
c variable cohesion, int angle frict, density, min viscosity, activation
c	energy, pre-exponential for mech model
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      read(2,101)ninc
      read(2,106)dummy
      read(2,106)dummy
      allocate(coh((nerowm-1)*(ncol-1)))
      allocate(phi((nerowm-1)*(ncol-1)))
      allocate(rhoc((nerowm-1)*(ncol-1)))
      allocate(vmin((nerowm-1)*(ncol-1)))
      allocate(q((nerowm-1)*(ncol-1)))
      allocate(prex((nerowm-1)*(ncol-1)))
      allocate(expn((nerowm-1)*(ncol-1)))
      do i=1,ninc
      	read(2,108)ibcol,iecol,ibrow,ierow,coht,phit,denst,vmint,
     *	qt,prext,expnt
      	if(iecol.gt.ncol-1) then
     		print*,'### ERROR:'
     		print*,'    ending colm in changing mech props'
     		print*,'    is greater than ncol-1'
     		stop
     	endif	
     	if(ierow.gt.nerowm-1) then
     		print*,'### ERROR:'
     		print*,'    ending row in changing mech props'
     		print*,'    is greater than nerowm-1'
     		print*,'    ',ierow,nerowm-1
     		stop
     	endif
      	do icol=ibcol,iecol
      		do irow=ibrow,ierow
      			coh((icol-1)*(nerowm-1)+irow)=coht
      			phi((icol-1)*(nerowm-1)+irow)=phit
      			rhoc((icol-1)*(nerowm-1)+irow)=denst
      			vmin((icol-1)*(nerowm-1)+irow)=vmint
      			q((icol-1)*(nerowm-1)+irow)=qt
      			prex((icol-1)*(nerowm-1)+irow)=prext
      			expn((icol-1)*(nerowm-1)+irow)=expnt
      		end do
      	end do
      end do
c number of boundary layers
      do i=1,8
      	read(2,106)dummy
      end do	
      read(2,111)iblayt,iflgblt,tblayt
      read(2,106)dummy
      read(2,111)iblay,iflgbl,tblay
c variable thermal properties for mech domain       
      do i=1,5
      	read(2,106)dummy
      end do
      read(2,101)ntmchg
      allocate(tm_prop(ntmchg,9))
      read(2,106)dummy
      do i=1,ntmchg
      	read(2,108)itemp1,itemp2,itemp3,itemp4,tm_prop(i,5),
     * 	tm_prop(i,6),tm_prop(i,7),tm_prop(i,8),tm_prop(i,9)
     	if(itemp2.gt.ncol-1) then
     		print*,'### ERROR:'
     		print*,'    ending colm in changing thermal props'
     		print*,'    is greater than ncol-1'
     		stop
     	endif	
     	if(itemp4.gt.nerowm-1) then
     		print*,'### ERROR:'
     		print*,'    ending row in changing thermal props'
     		print*,'    is greater than nerowm-1'
     		print*,'    ',itemp4,nerowm-1
     		stop
     	endif
        tm_prop(i,1)=dble(itemp1)
        tm_prop(i,2)=dble(itemp2)
        tm_prop(i,3)=dble(itemp3)
        tm_prop(i,4)=dble(itemp4)
       end do 
c rigid viscosity
      read(2,106)dummy
      read(2,103)vrig
c compressibility
      read(2,106)dummy
      read(2,103)beta
      print*,beta
c flag for linear or non-linear eqns
      read(2,106)dummy
      read(2,101)linflag
c flag for just plastic def
      read(2,106)dummy
      read(2,101)iplasflg
c epsinv
      read(2,106)dummy
      read(2,105)epsinv
c tmax
      read(2,106)dummy
      read(2,103)tmax
c densities of fluid and mantle
      do i=1,4
      	read(2,106)dummy
      end do	
      read(2,103)rhof,rhom
c # of boundary nodes: vel, pressure, loaded sides, tan vel
      read(2,106)dummy
      read(2,106)dummy
      read(2,101)numvebn,numvetbn,numpbn,numsid,numvtbn
      if(numvtbn.ne.ncol) then
     	print*,'###################### ERROR ##########################'
     	print*,'## need to have the number of tan vel bcs match ncol ##'
     	print*,'###################### ERROR ##########################'
     	stop
      end if	
c # tsteps, out int, out int lagrangian, tstep length
      read(2,106)dummy
      read(2,109)ntst,intout,intoutl,delt
c min inter, max iter, num filtering passes, convergence tolerance
      read(2,106)dummy
      read(2,109)minint,maxint,npass,toler
c erosion parameters
      read(2,106)dummy
      read(2,103)erosl,erosr,peros,rpow
c sedimentation paramters 
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      read(2,112)ipkfill,ibasfill,isedl,isedr,sedmax
      if(isedl.gt.ncol.or.isedr.gt.ncol) then
      	print*,'######################################################'
      	print*,'## ERROR: bounds of sedimenation need to be within  ##'
      	print*,'##	the model bounds. ncol=',ncol
      	print*,'######################################################'
      	stop
      endif	
      if(ipkfill.eq.0.and.ibasfill.eq.1) then
      	print*,'######################################################'
      	print*,'## ERROR: cannot fill bounding basins w/o filling   ##'
      	print*,'##	between peaks                                   ##'
      	print*,'######################################################'
      	stop
      endif	
c basin tracking parameters
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      read(2,101)ibasflg,nbastary,nbastind,intmrkb
c maximum slope
      read(2,106)dummy
      read(2,106)dummy
      read(2,103)slpmax
c thermal runup
      read(2,106)dummy
      read(2,104)ntt2,deltt2
c variable thermal props
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      allocate(therm_prop(5,5))
      do i=1,5
      	read(2,103)therm_prop(i,1),therm_prop(i,2),
     *	therm_prop(i,3),therm_prop(i,4),therm_prop(i,5)
      end do
c BCs for thermal problem
c	ss bcs      
      do i=1,7
      	read(2,106)dummy
      end do	
      allocate(thermbcs(1,2))
      read(2,103)thermbcs(1,1),thermbcs(1,2)
c	cooling lithos bcs      
      do i=1,4
      	read(2,106)dummy
      end do
      read(2,104)iflgcl,agecl
c read in bcs for mech model
      read(2,106)dummy
      read(2,106)dummy
      read(2,106)dummy
      allocate(bc(5,numvebn+numvetbn+numpbn+numsid+numvtbn))
      index1=0
c	x vel at model edges      
      read(2,101)nsets
      idum=1
      do 300 is=1,nsets
      	read(2,110)num,istart,inc,val
      	ifin=istart+(num-1)*inc
      	do 320 i=istart,ifin,inc
c      		write(3,111) i,idum,val
      		index1=index1+1
      		bc(1,index1)=float(i)
      		bc(2,index1)=float(idum)
      		bc(3,index1)=val
  320 	continue
  300 continue
c	y vel at model edges
      read(2,101)nsets
      idum=2
      do 400 is=1,nsets
      	read(2,110)num,istart,inc,val
      	ifin=istart+(num-1)*inc
      	do 420 i=istart,ifin,inc
c	        write(3,111) i,idum,val
      		index1=index1+1
      		bc(1,index1)=float(i)
      		bc(2,index1)=float(idum)
      		bc(3,index1)=val
  420 	continue
  400 continue
c  pressure
      read(2,101)nsets
      do 500 is=1,nsets
      	read(2,110)num,istart,inc,val
      	ifin=istart+(num-1)*inc
      	do 520 i=istart,ifin,inc
c	        write(3,111) i,'0.0',val
      		index1=index1+1
      		bc(1,index1)=float(i)
      		bc(2,index1)=0.0
      		bc(3,index1)=val
  520 	continue
  500 continue
c tangential velocities at model edges
c	insetad of specifying x and y vel at the lhs and rhs model edges, 
c	specify a velocity tangential to the first two nodes at the base 
c	of the model
      read(2,101)nsets
      do 425 is=1,nsets
      	read(2,110)num,istart,inc,val
      	ifin=istart+(num-1)*inc
      	do 435 i=istart,ifin,inc
c	        write(3,111) i,val,'0.0'
      		index1=index1+1
      		bc(1,index1)=float(i)
      		bc(2,index1)=val
      		bc(3,index1)=0.0
  435 	continue
  425 continue
c	tangential velocity  
      read(2,101)nsets
      do 800 is=1,nsets
      	read(2,110)num,istart,inc,val,val2
      	ifin=istart+(num-1)*inc
      	do 820 i=istart,ifin,inc
c	        write(3,111) i,val,'0.0'
      		index1=index1+1
      		bc(1,index1)=float(i)
      		bc(2,index1)=val
      		bc(3,index1)=0.0
  820 	continue
  800 continue
c  loaded sides
      read(2,101)nsets
      if(nsets.eq.0)go to 601
      print*,'#######################################################'
      print*,'##  WARNING: mesh generator is not set up to process ##'
      print*,'##           loaded sides. also not in sub. output   ##'
      print*,'#######################################################'      
      stop
      do 600 is=1,nsets
      	read(2,110)num,istart,inc,val
      	ifin=istart+(num-1)*inc
      	do 620 i=istart,ifin,inc
      		ii=i+nrowc
      		iii=ii+nrowc
      		iiii=1
c	        write(3,111) i,ii,iii,iiii,val
      		index1=inde1+1
  620 	continue
  600 continue
  601 continue
c backstop nodes
      read(2,101)nsets
      if(nsets.eq.0)go to 701
      print*,'#######################################################'
      print*,'##  WARNING: mesh generator is not set up to process ##'
      print*,'##           backstop nodes. (also not in sub. output)#'
      print*,'#######################################################'      
      stop
      do 700 is=1,nsets
      	read(2,110)num,istart,inc
      	ifin=istart+(num-1)*inc
      	do 720 i=istart,ifin,inc
c      	  write(3,111) i
  720 	continue
  700 continue
  701 continue

c
c read in flags for output files
c

c numer of possible files
      do i=1,4
      	read(2,106)dummy
      end do	
      read(2,101)noutput
      allocate(output_flags(noutput))
      read(2,106)dummy
      icount=1
c coord
      read(2,113)dummy,iflag
      output_flags(icount)=iflag
c vel      
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c press      
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c stresses and stuff      
      read(2,106)dummy
      do i=1,7
      	read(2,113)dummy,iflag
      	icount=icount+1
      	output_flags(icount)=iflag
      end do	
c strain rates, dilitation      
      read(2,106)dummy
      do i=1,6
      	read(2,113)dummy,iflag
      	icount=icount+1
      	output_flags(icount)=iflag
      end do
c lmesh      
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c temp_mech      
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c visc_elem and visc_gp
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c erosion
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c temp_track      
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag      
c unvel      
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c exhum      
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c surf_prof      
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c duc_flag      
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c material props      
      read(2,106)dummy
      do i=1,7
      	read(2,113)dummy,iflag
      	icount=icount+1
      	output_flags(icount)=iflag
      end do
c basinfill
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c peakchop
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c basin_track
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c l_temp_all
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c coordt
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c velthermal_alt
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c velthermal
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c temp.dat (entire model)      
      read(2,106)dummy
      read(2,113)dummy,iflag
      icount=icount+1
      output_flags(icount)=iflag
c thermal props
      read(2,106)dummy
      do i=1,3
      	read(2,113)dummy,iflag
      	print*,dummy,iflag
      	icount=icount+1
      	output_flags(icount)=iflag
      end do
c update size of output flag array
      noutput=icount

      if(inflag.eq.1) close(9)
      close(2)
  101 format(9i5)
  103 format(6e10.2)
  104 format(i5,2e10.2)
  105 format(f7.2,i5)
  106 format(a10)
  107 format(i2,i4,e10.2)  
  108 format(4i4,9e10.2)
  109 format(3i5,e10.2)
  110 format(3i5,2d15.6)
  111 format(2i5,e10.2)
  112 format(4i5,e10.2)
  113 format(a15,i5)
  239 format(3e15.5)  

      end
      
c#######################################
c make x,y arrays for each plate
c#######################################
      subroutine mk_plates(ncol,npad,xsing,xadd,np1,nsing,nsing1,np2)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c PLATE 1
c calculate how many nodes to include for plate one to extend
c	xadd past the spoint
      if(xsing+xadd.gt.xbase(npad*2+ncol)) then
      	np1=ncol+npad*2
      	print*,'WARNING: Sub. extension goes past end of model'
      else 	
      	do i=1,ncol+npad*2
   			if(xbase(i).gt.xsing+xadd) then
  				if(abs(xsing+xadd-xbase(i)).
     *			ge.abs(xadd+xsing-xbase(i-1))) then
  					np1=i-1
  					exit
   				else
      				np1=i
	      			exit
	      		endif	
      		else if(xbase(i).eq.xadd+xsing) then
      			np1=i
      			exit
      		endif
      	end do	
      endif	
      	print*,'length of plate 1 with sub exten (desired,used):'
      	print*,'	',xsing+xadd,xbase(np1)
c make x array for plate 1
      allocate(xp1(np1))
      allocate(yp1(np1))
      allocate(dyinit1(np1))
      do i=1,np1
      	xp1(i)=xbase(np1)-xbase(np1-i+1)
      	yp1(i)=0.0
      	dyinit1(i)=0.0
      end do
      if(np1.gt.nsing+npad) then
      	nsing1=np1-(nsing+npad)+1
      else
      	nsing1=1
      endif	
c make initial mech. model thickness array for plate 1      
c	ensure thickness past spoint is zero
      do i=1,nsing+npad
      	dyinit1(np1-i+1)=dyinit(i)
      end do	
      do i=1,np1
      end do
c PLATE 2
c make x array
c 	note: plate 1 and 2 have the same begining node at s point
      np2=ncol+2*npad-(nsing+npad)+1 
      allocate(yp2(np2))
      allocate(xp2(np2))
      allocate(dyinit2(np2))
      do i=1,np2
      	xp2(i)=xbase(i+nsing+npad-1)-xbase(nsing+npad)
      	yp2(i)=0.0
      	dyinit2(i)=0.0
      end do	
c make initial mech. model thickness array for plate 2      
      do i=1,np2
      	dyinit2(i)=dyinit(nsing-1+i+npad)
      end do	
      end

c#########################################################
c main loop to calculate flexure of the two plates
c#########################################################
      subroutine calc_flex(nerowm,ncol,np1,np2,prigi,rrigi,rhom,
     *npad,nsing,nsing1,ctoler,sload,smomen,wdepth,wtoler,rhof,
     *ypmbase,yrmbase,wheight,inflag,dyc)
      use dyn_arrays
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
      g=9.8
      alpha1=(4.0*prigi/((rhom)*g))**0.25
      alpha2=(4.0*rrigi/((rhom)*g))**0.25
      plam1=1.0/alpha1
      plam2=1.0/alpha2
      fk=rhom*g

c########################################################
c caculate force on each node from mech. model thickness
c########################################################
      allocate(slen1(np1))
      allocate(fnode1(np1))
      allocate(slen2(np2))
      allocate(fnode2(np2))

c calculate average density of the colms at the LHS and RHS of model domain
c 	for use in calculating the force on the padded sections of the profiles.  
c	this is not really an average density since it does not take into account
c	differences in the size of elements, but it will work for this use
      rhosum=0.0
      do j=1,nerowm-1
      	rhosum=rhosum+rhoc(j)
      end do
      rhoavtl=rhosum/dble(nerowm-1)
      rhosum=0.0
      do j=1,nerowm-1
      	rhosum=rhosum+rhoc((nerowm-1)*(ncol-2)+j)
      end do
      rhoavtr=rhosum/dble(nerowm-1)
c plate 1
      call calc_force_p1(slen1,xp1,nerowm,rhoc,np1,fnode1,dyinit1,
     *npad,1,rhoavtl,nsing,g,nsing1,inflag,dencol)
c plate 2
      call calc_force_p2(slen2,xp2,nerowm,rhoc,np2,fnode2,dyinit2,
     *npad,2,rhoavtr,nsing,g,ncol,inflag,dencol)
c shift initial y positions for intial offsets defined above
      do i=1,np1
      	yp1(i)=-ypmbase
      end do
      do i=1,np2
      	yp2(i)=-yrmbase
      end do	
c calculate defection of plates from distributed load due to thickness of
c 	mech model, subduction end load/moment, shift in coupling point,
c	moment at coupling point and coupling constraint
      call deflect(np1,np2,xp1,xp2,yp1,yp2,fnode1,fnode2,plam1,plam2,
     *fk,ctoler,nsing1,sload,smomen,xbase,nsing,npad,dyc)
c calculate the deflection from the load of overlying water
      if(wdepth.gt.0.0) then
     	call deflectw(wdepth,xp1,xp2,slen1,slen2,wtoler,fnode1,
     *	fnode2,plam1,plam2,fk,rhof,np1,np2,npad,g,nsing1,ctoler,
     *	dyinit1,dyinit2,yp2,yp1,wheight)
      	dif=0.0
      endif

c###################################
c check that plates do not overlap
c###################################
      do i=1,nsing1-1
      	if(abs(yp1(i)).lt.abs(yp2(nsing1-i+1))) then
      		print*,'###############################'
      		print*,'## ERROR: plates overlap     ##'
      		print*,'###############################'
      		call profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
      		stop
      	endif
      end do	
      deallocate(fnode1)
      deallocate(fnode2)
      deallocate(slen1)
      deallocate(slen2)
      end

c###############################################################
c make arrays of flexure profile for model (mech and sub lithos)
c###############################################################
      subroutine mech_bndry(ncol,nsing1,nsing,np1,npad,np2,
     *plthick,athick,yshift)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
c make arrays mech model upper/lower boundary, sub. lithos
      allocate(ymbase(ncol),ymtop(ncol))
      allocate(yltemp(nsing1))
      allocate(xltemp(nsing1))
      allocate(xmbase(ncol))
c mech model boundary      
      do i=1,nsing
      	ymbase(i)=-yp1(np1-npad+1-i)
      	ymtop(i)=-yp1(np1-npad+1-i)+dyinit(i+npad)
      	xmbase(i)=xbase(i+npad)
      end do	
      do i=1,np2-npad-1
      	ymbase(i+nsing)=-yp2(i+1)
      	ymtop(i+nsing)=-yp2(i+1)+dyinit(nsing+i+npad)
      	xmbase(i+nsing)=xbase(i+nsing+npad)
      end do
c top of sub lithos past spoint      
      do i=1,nsing1
      	xltemp(i)=xmbase(nsing+i-1)
      	yltemp(i)=-yp1(nsing1-i+1)
      end do	
c shift all arrays so the y=0 is defined by model base on pro-side
      yshift=ymbase(1)-plthick-athick
      do i=1,ncol
      	ymbase(i)=ymbase(i)-yshift
      	ymtop(i)=ymtop(i)-yshift
      end do
      do i=1,nsing1
      	yltemp(i)=yltemp(i)-yshift
      end do	
      end

c#######################################################################
c find the closest node in the xdimension giving an x psoition
c#######################################################################
      subroutine find_node(node,ntype,xwant,ncol,npad,xuse)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      node=0
      if(ntype.eq.0) then
        node=int(xwant)
      else
        do i=1+npad,ncol+npad
            if(xbase(i).gt.xwant) then
                if(abs(xwant-xbase(i)).ge.abs(xwant-xbase(i-1))) then
                    node=i-1-npad
                    exit
                else
                    node=i-npad
                    exit
                endif
            else if(xbase(i).eq.xwant) then
                node=i-npad
                exit
            endif
        end do
      endif
      xuse=xbase(node+npad)
      end

c ####################################################################
c output mesh and parameters
c#####################################################################
      subroutine output(ncol,nerowm,ntrow,nsing,plvel,upvel,
     *iunflag,iunbeg,vrig,beta,epsinv,rhof,
     *rhom,iso,prigi,rrigi,sload,smomen,xadd,ctoler,wdepth,wtoler,
     *numvebn,numpbn,numsid,numvtbn,ntst,delt,intout,intoutl,minint,
     *maxint,npass,toler,erosl,erosr,peros,rpow,ntt2,deltt2,np1,np2,
     *nsing1,npad,nplbase,npltop,plscale,rlscale,blscale,dfact,slpmax,
     *tmax,nrowl,plthick,numvetbn,linflag,iplasflg,iblay,iblayt,isedl,
     *isedr,iexflg,ibasflg,nbastary,nbastind,intmrkb,ipkfill,ibasfill,
     *sedmax,ntbcs,noutput)
      use dyn_arrays
      implicit real*8(a-h,o-z)
      implicit integer(i-n)

      open(3,file='input/mesh',position='rewind')
      open(1,file='input/connections.dat',position='rewind')
c mech model (# nodes, # elements, l-mesh style)      
      nnodesm=ncol*nerowm
      nelem=(ncol-1)*(nerowm-1)
      write(3,102)nnodesm,nelem
c mech model (# rows/ colms of element verticies
      write(3,102)nerowm,ncol
c lmesh parameters
      write(3,103)plscale,rlscale,blscale,dfact
c thermal model (#nodes, # elements, # rows in lith)
      nnodest=ncol*ntrow
      nelet=(ncol-1)*(ntrow-1)*2
      write(3,102)nnodest,nelet
      write(3,102)ntrow,ncol,nrowl
c ref. thickness of lithosphere. used in thermal remeshing
      write(3,103)plthick
c spoint node
      write(3,102)nsing
c convergence and underplating velocity
      write(3,119)plvel,upvel
c underplating parameter
      if(iunflag.eq.2) iunflag=0
      write(3,107)iunflag,iunbeg
c rigid visc      
      write(3,103)vrig
c compressibility, epsinv, tmax
      write(3,103)beta,epsinv,tmax
c flag for using linear or non-linear eqns      
      write(3,102)linflag
c flag for allowing purely plastic def (no viscous)
      write(3,102)iplasflg
c overlying fluid and mantle density
      write(3,103)rhof,rhom
c flexural paramters
      write(3,102)iso
      write(3,103)prigi,rrigi,sload,smomen
      write(3,103)xadd,ctoler,wdepth,wtoler
c     if not outputting extension for plasti, change nsing1 and np1      
      if(iexflg.eq.1) then
      	write(3,102)np1-nsing1+1,np2,npad,1
      else
      	write(3,102)np1,np2,npad,nsing1
      endif	
c number of boundary conditions
      write(3,102)numvebn,numvetbn,numpbn,numsid,numvtbn
c number of timesteps, timestep length
      write(3,104)ntst,delt
c output interval for all, output interval for lmesh
      write(3,102)intout,intoutl
c min iter, max iter, # filtering passes, conv. toler.
      write(3,109)minint,maxint,npass,toler
c erosion parameters
      write(3,103)erosl,erosr,peros,rpow
c sedimentation parameters
      write(3,112)ipkfill,ibasfill,isedl,isedr,sedmax
c basin tracking parameters
      write(3,101)ibasflg,nbastary,nbastind,intmrkb
c maximum surface slope
      write(3,103)slpmax
c thermal runup parameters
      write(3,104)ntt2,deltt2
c number of bounadry layers
      write(3,102)iblay,iblayt
c output variable mat. props      
      do i=1,nelem
      	write(3,113)rhoc(i),phi(i),coh(i),vmin(i),q(i),prex(i),expn(i)
      end do	
c output node coordinates
      do i=1,nnodest
      	write(3,113)pos(i,1),pos(i,2)
      end do	
c output the slope on the lhs and rhs base
      dx=pos(ntrow-nerowm+1,1)-pos(ntrow*2-nerowm+1,1)
      dy=pos(ntrow-nerowm+1,2)-pos(ntrow*2-nerowm+1,2)
      slp1=atan(dy/dx)
      dx=pos(nnodest-nerowm+1,1)-pos(nnodest-ntrow-nerowm+1,1)
      dy=pos(nnodest-nerowm+1,2)-pos(nnodest-ntrow-nerowm+1,2)
      slp2=atan(dy/dx)
      print*,'#### Slope at base of model (deg)'
      print*,'####   lhs =',slp1*180.0/3.1416
      print*,'####   rhs =',slp2*180.0/3.1416
c output node connections for mech model
      do icol=1,ncol-1
      	do irow=1,nerowm-1
      		n1=(icol-1)*nerowm+irow
      		n2=n1+nerowm
      		n3=n2+1
      		n4=n1+1
      		write(3,102)n1,n2,n3,n4
      	end do
      end do	
c output boundary conditions
c	velocity boundary conditions on model edges
      index=0
      if(numvebn.gt.0) then
c		x and y vel      
      	do i=1,numvebn
      		index=index+1
      		write(3,114)int(bc(1,index)),int(bc(2,index)),bc(3,index)
      	end do
      end if	
c 	pressure boudary conditions
      if(numpbn.gt.0) then
      	do i=1,numpbn
      		index=index+1
      		write(3,115)int(bc(1,index)),bc(3,index)
      	end do
      end if
c edge tangential vel boundary conditions
      if(numvetbn.gt.0) then
      	do i=1,numvetbn
      		index=index+1
      		write(3,115)int(bc(1,index)),bc(2,index)
      	end do
      endif	
c basal tangential vel boundary conditions
      if(numvtbn.gt.0) then
      	do i=1,numvtbn
      		index=index+1
      		write(3,115)int(bc(1,index)),bc(2,index),bc(3,index)
      	end do
      endif	
c THERMAL Output
c number of nodes, number of elements, number of domains, number of temp BCs
      write(3,102)nnodest,nelet,5,ntbcs
c mesh connections
      write(1,102)nelet
      do icol=1,ncol-1
      	do irow=1,ntrow-1
      		n1=(icol-1)*(ntrow)+irow
      		n2=n1+ntrow
      		n3=n1+1
      		write(1,102)n1,n2,n3
      		n1=n1+1
      		n2=n2
      		n3=n2+1
      		write(1,102)n1,n2,n3
      	end do
      end do	
c domain definitions
      do i=1,nelet
      	write(3,102)ndomain(i)
      end do	
c Variable properties
      do i=1,nelet
      	write(3,113)therm_cond(i,1),therm_cond(i,2)
      end do
      do i=1,nelet
      	write(3,113)therm_rho(i)
      end do
      do i=1,nelet
      	write(3,113)spec_heat(i)
      end do
      do i=1,nelet
      	write(3,113)heat_prod(i)
      end do	
c BCs
c	const temp nodes
      do i=1,ntbcs
      	write(3,115)int(therm_bc(1,i)),therm_bc(2,i)
      end do

c output flexure arrays in x and y for use in isostacy calc in plasti 
c     if iexflg=1 do not ouput the extension for use in the plasti 
c     flexure calculation
      if(iexflg.eq.1) then
      	do i=nsing1,np1
      		write(3,113)xp1(i)-xp1(nsing1),yp1(i),dyinit1(i)
      	end do	
      else
      	do i=1,np1
      		write(3,113)xp1(i),yp1(i),dyinit1(i)
      	end do
      endif	
      do i=1,np2
      	write(3,113)xp2(i),yp2(i),dyinit1(i)
      end do	
c
c output file output flags
c
      open(4,file='output/output_flags',position='rewind')
      write(3,101)noutput
      write(4,101)noutput
      do i=1,noutput
      	write(3,101)output_flags(i)
      	write(4,101)output_flags(i)
      end do	
      close(3)      
      close(4)
      close(1)

c output the node numbers of the domain boundaries
      open(3,file='input/boundary_nodes.dat')
c	  top of mech model      
      write(3,102)ncol
      do i=1,ncol
      	write(3,102)mecht_nodes(i)
      end do
c     base of mech model
      write(3,102)ncol
      do i=1,ncol
      	write(3,102)mechb_nodes(i)
      end do
c     base of pro lith
      write(3,102)nplbase
      do i=1,nplbase
      	write(3,102)plithb_nodes(i)
      end do	
c     top of pro lith
      write(3,102)npltop
      do i=1,npltop
      	write(3,102)plitht_nodes(i)
      end do	
c     base of retro lith
      write(3,102)ncol-nsing
      do i=1,ncol-nsing
      	write(3,102)rlithb_nodes(i)
      end do	

  101 format(9i5)
  102 format(9i8)
  103 format(4e16.8)
  104 format(i5,2e16.8)
  107 format(i2,i4,e16.8)  
  109 format(3i5,e16.8)
  112 format(4i5,e10.2)  
  113 format(9e23.15)
  114 format(2i8,4e23.15)
  115 format(i8,4e23.15)
  119 format(2f8.1,e13.8)
      end

c######################################################
c make array of variable tehrm properties
c######################################################
      subroutine mk_therm_para(ntrow,ncol,ntmchg,nrowl,nrowa,
     *ntbcs,iflgcl,agecl,nplbase,npltop)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      real*4 arg

C THERMAL BCs
c	allocate space      
      if(iflgcl.eq.1) then
      	ntbcs=2*ncol+nplbase-npltop+ntrow-nrowa
      	allocate(therm_bc(2,ntbcs))
      else
      	ntbcs=2*ncol+nplbase-npltop
      	allocate(therm_bc(2,ntbcs))
      endif	
c	calculate the temp for 1-d,semi-infinite cooling lithos
      index=0
      if(iflgcl.eq.1) then
      	time=agecl*3.15578e13
      	rho=therm_prop(2,3)
      	cp=therm_prop(2,4)
      	tcond=therm_prop(2,2)
      	tsurf=thermbcs(1,1)
      	tmant=thermbcs(1,2)
      	ysurf=pos(ntrow,2)
      	do i=ntrow,nrowa,-1
      		index=index+1
      		depth=ysurf-pos(i,2)
      		arg=depth/(2.0*(tcond*time/(rho*cp))**(.5))
      		temp=erfc(arg)*(tsurf-tmant)+tmant
      		therm_bc(1,index)=dble(i)
      		therm_bc(2,index)=temp
      	end do	
      	print*,'## Temp in asthen:    ',tmant
      	print*,'## Temp at Lith base: ',therm_bc(2,index)
      endif
c	temp at surface
      do i=1,ncol
      	index=index+1
      	therm_bc(1,index)=dble(ntrow*i)
      	therm_bc(2,index)=thermbcs(1,1)
      end do
c	temp at base
      do i=1,nplbase
      	index=index+1
      	therm_bc(1,index)=dble(1+(i-1)*ntrow)
      	therm_bc(2,index)=thermbcs(1,2)
      end do
      do i=npltop+1,ncol
      	index=index+1
      	therm_bc(1,index)=dble(1+(i-1)*ntrow)
      	therm_bc(2,index)=thermbcs(1,2)
      end do

C THERMAL PROPS
      ntele=(ncol-1)*(ntrow-1)*2
      allocate(therm_cond(ntele,2))
      allocate(therm_rho(ntele))
      allocate(spec_heat(ntele))
      allocate(heat_prod(ntele))
c initial definitions from thermal domains
      do i=1,ntele
      	therm_cond(i,1)=therm_prop(ndomain(i),1)
      	therm_cond(i,2)=therm_prop(ndomain(i),2)
      	therm_rho(i)=therm_prop(ndomain(i),3)
      	spec_heat(i)=therm_prop(ndomain(i),4)
      	heat_prod(i)=therm_prop(ndomain(i),5)
      end do	
c change in therm props in the mech model
      do j=1,ntmchg
      	ibcol=int(tm_prop(j,1))
      	iecol=int(tm_prop(j,2))
      	ibrow=int(tm_prop(j,3))
      	ierow=int(tm_prop(j,4))
      	nele_row=(ntrow-1)*2
      	nele_al=(nrowa+nrowl)*2
      	do icol=ibcol,iecol
      		do irow=ibrow,ierow
      			iele1=nele_al+irow*2-1+(icol-1)*nele_row
      			iele2=nele_al+irow*2+(icol-1)*nele_row
      			therm_cond(iele1,1)=tm_prop(j,5)
      			therm_cond(iele2,1)=tm_prop(j,5)
      			therm_cond(iele1,2)=tm_prop(j,6)
      			therm_cond(iele2,2)=tm_prop(j,6)
      			therm_rho(iele1)=tm_prop(j,7)
      			therm_rho(iele2)=tm_prop(j,7)
      			spec_heat(iele1)=tm_prop(j,8)
      			spec_heat(iele2)=tm_prop(j,8)
      			heat_prod(iele1)=tm_prop(j,9)
      			heat_prod(iele2)=tm_prop(j,9)
      		end do
      	end do
      end do	
      end

c##############################################################
c output the domain boundaries for plotting
c##############################################################
      subroutine bndry_output(wheight,ncol,npad,wdepth,np1,np2,yshift
     *,nsing,nsing1,npltop,nrlbase,nplbase)
      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      open(2,file='profiles/sealevel',action='write')
      open(3,file='profiles/mech_top',action='write')
      open(4,file='profiles/plith_top',action='write')
      open(7,file='profiles/rlith_top',action='write')
      open(8,file='profiles/plith_base',action='write')
      open(9,file='profiles/rlith_base',action='write')
      open(10,file='profiles/mech_base_slope',action='write')
      open(11,file='profiles/mech_top_slope',action='write')
      open(12,file='profiles/mech_base_slope_rise_run',action='write')

c ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c output domains for plotting
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c output sealevel over extended model domain
      wref=wheight-yshift
      print*,'water height',wheight,yshift
      if(wdepth.gt.0.0) then
      	do i=1,ncol+2*npad
      		write(2,113)xbase(i)/1000.0,(wheight-yshift-wref)/1000.0
	    end do	
	  else
      	do i=1,ncol+2*npad
      		write(2,*)'0.0  0.0\n'
	    end do
	  endif   
c output top of mech model over plate 1 extended model domain	  
      index=0
      do i=1,nsing+npad
        ip=np1-i+1
        index=index+1
      	write(3,113)xbase(i)/1000.0,
     *	(-yp1(ip)+dyinit1(ip)-yshift-wref)/1000.0
      end do	
c ouput slope of top/bottom of model plate 1
      index2=npad
      do i=npad+1,nsing+npad
        index2=index2+1
        ip=np1-i+1
      	dx=xbase(i)-xbase(i-1)
      	dy=(-yp1(ip)+dyinit1(ip))-(-yp1(ip+1)+dyinit1(ip+1))
      	dy2=(-yp1(ip))-(-yp1(ip+1))
      	slp=atan(dy/dx)*180.0/3.141592
      	slp2=atan(dy2/dx)*180.0/3.141592
      	slp3=dy2/dx
      	write(11,113)(xbase(i)-dabs(dx-2.0))/1000.0,slp
      	write(10,113)(xbase(i)-dabs(dx-2.0))/1000.0,slp2
      	write(12,113)(xbase(i)-dabs(dx-2.0))/1000.0,slp3
      end do	
c output top of mech model over plate 2 extended model domain	  
      do i=1,np2
      	index=index+1
      	write(3,113)xbase(index-1)/1000.0,
     *	(-yp2(i)+dyinit2(i)-yshift-wref)/1000.0
      end do      
c output slope of top of model plate 2
      do i=1,np2-npad
      	index2=index2+1
      	dx=xbase(index2-1)-xbase(index2)
      	dy=(-yp2(i)+dyinit2(i))-(-yp2(i+1)+dyinit2(i+1))
      	dy2=(-yp2(i))-(-yp2(i+1))
      	slp=atan(dy/dx)*180.0/3.141592
      	slp2=atan(dy2/dx)*180.0/3.141592
      	slp3=dy2/dx
      	write(11,113)(xbase(index2-1)+dabs(dx/2.0))/1000.0,slp
      	write(10,113)(xbase(index2-1)+dabs(dx/2.0))/1000.0,slp2
      	write(12,113)(xbase(index2-1)+dabs(dx/2.0))/1000.0,slp3
      end do   
c output top of pro-lithosphere to end of plate 1
      do i=1,np1
      	ip=np1-i+1
	  	write(4,113)xbase(i)/1000.0,(-yp1(ip)-yshift-wref)/1000.0
	  end do	
c output top of pro-lith from end of plate 1 to model base	  
      do i=np1+1,npltop+npad
      	write(4,113)xbase(i)/1000.0,(ypltop(i-npad)-wref)/1000.0
      end do	
c ouput top of retro lith
      do i=1,np2
      	write(7,113)xbase(nsing-1+i+npad)/1000.0,
     *	(-yp2(i)-yshift-wref)/1000.0
      end do	
c output base of pro-lith
      do i=1,nplbase
      	write(8,113)xmbase(i)/1000.0,(yplbase(i)-wref)/1000.0
      end do	
c output base of retro-lith      
      do i=1,nrlbase
      	write(9,113)xrlbase(i)/1000.0,(yrlbase(i)-wref)/1000.0
      end do	

      print*,' Y pos of coupling point relative to sea level (km)'
     *,(-yp2(1)-yshift-wref)/1000.0

      close(9);close(2);close(3);close(4);close(7);close(8)
      close(10);close(11);close(12)
  113 format(4e16.8)
      end

c###########################################################
c make initial profiles from circular arcs defined by rigid
c###########################################################

      subroutine arc_prof(nerowm,ncol,np1,np2,prigp,rrigp,rhom,
     *npad,nsing,nsing1,wdepth,ypmbase,yrmbase,wheight,inflag,
     *dyc,itrench,sdip)

      use dyn_arrays
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

c slab dip in radians
      sdip=sdip*3.141592/180.0
c flexural parameters      
      g=9.8
      alpha1=(4.0*prigp/((rhom)*g))**0.25
      alpha2=(4.0*rrigp/((rhom)*g))**0.25
c define curvature radius as function of flexural parameter
      rad1=3.141592*alpha1/2.0
      rad2=3.141592*alpha2/2.0

C PLATE 1
c define initial depth from offsets in plates between model edge
c	and the trench (begining of arc)
      itrench1=np1-(itrench+npad)+1
      do i=itrench1,np1
      	yp1(i)=-ypmbase
      end do	
c set depth from trench landward from circular arc till prescribed
c 	dip is reached
      icatch=0
c     define center of arc(xnot,ynot)
      xnot=xp1(itrench1)
      ynot=-ypmbase-rad1
c     calculate arc      
      do i=itrench1-1,1,-1
      	yp1(i)=ynot+dsqrt(rad1**2-(xp1(i)-xnot)**2)
      	dip=datan((yp1(i+1)-yp1(i))/(xp1(i+1)-xp1(i)))
      	if(dip.gt.sdip) then
      		icatch=1
      		exit
      	endif
      end do
      if(icatch.eq.0) then
      	print*,'ERROR: circular arc on pro side did not reach sdip'
      	print*,'	increase the length of the extension'
      	call profdump(xbase,-yp1,-yp2,np1,np2,nsing,npad)
      	stop
      endif
c     set remaining with prescribed slope      
      do j=i,1,-1
      	yp1(j)=yp1(j+1)-(xp1(j+1)-xp1(j))*dtan(sdip)
      end do

C PLATE 2
c check that plate two is not below the level of the spoint, this
c	is not a valid solution
      if(yp1(nsing1).gt.-yrmbase) then
      	print*,'ERROR: the elevation of plate 2 is below the spoint'
      	print*,'	could increase slab dip, decrease plate offset,'
      	print*,'	do anything to lower the level of spoint.'
      	print*,'	y(s)=',yp1(nsing1),'y plate 2=',-yrmbase
      	yp2(1)=yp1(nsing1)
      	do i=1,np2
      		yp2(i)=-yrmbase
      	end do
      	call profdump(xbase,-yp1,-yp2,np1,np2,nsing,npad)
      	stop
      endif
c define center of arc for plate 2 (xnot,ynot)
      xnot=dsqrt(rad2**2-(-rad2-yp1(nsing1)-yrmbase)**2)
      ynot=-rad2-yrmbase
      yp2(1)=yp1(nsing1)
      icatch=0
c     calculate arc      
      do i=2,np2
      	yp2(i)=ynot+dsqrt(rad2**2-(xp2(i)-xnot)**2)
      	if(xp2(i).ge.xnot) then
      		icatch=1
      		exit
      	endif
      end do
      if(icatch.eq.0) then
      	print*,'ERROR: circular arc on retro side did not reach zero'
      	print*,'	slope.'
      	call profdump(xbase,yp1,yp2,np1,np2,nsing,npad)
      	stop
      endif
c     set remaining with initial position      
      do j=i,np2
      	yp2(j)=-yrmbase
      end do

c switch sign so that deflections down are positive      
      yp2=yp2*(-1.0)
      yp1=yp1*(-1.0)

c set water height for plotting
      wheight=-yp1(np1-npad)+dyinit1(np1-npad)+wdepth

      end


