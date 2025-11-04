c
c  heatflow/smoker
c  ----------------
c
c  heat, c1,c2
c  c1 can be age or suspended particles with retention (MK thesis), c2 does not contribute to density
c
c  black smoker version June/July 1993
c  1D vertical line elements added
c  read ilx1,ilx2,ily1,ily2,ilz1,ilz2,idim,iloc,aperture(2b)
c
c  2D planar fractures added June 1994
c  read ilx1,ilx2,ily1,ily2,ilz1,ilz2,idim,iloc,full aperture(2b)
c
c  This version: 
c      - debugged January 1995
c      - cosmetic changes Jan 2004, tecplot output, kprint removed
c      - porosity read by element, hardwired for fractures

c hardwires Oct 2004: v in fracture, thermal K of fracture
c Feb 2007: fixed flux bug for face 3 & 4 - surf
c March 2007 ... added first order decay for mass transport (kmass=1)
c                decay works for direct integration only
c                works for 3D porous blocks and fractures
c                decay validated against line2d and craflush
c
c July 2009 ... exact derivative for wu
c               latent heat reformulation
c
c August 2010 ... 
c   - gamma for mass transport density added - set gamma=true concentration for dense mass transport
c   - viscosity & density functions may be hardwired for persulphate
c   - check for decay ... 
c
c Feb 2021 ... use C/Co when using gamma=0.024
c
c  small bug fixed in inl() vs inline() for fractures - had not affected results.
c
c Feb2021 .... madiha version depc,vexp,smax
c Jan2022 .... depd (detachment coefficient)
c Feb2022 .... smax a vector
c
c  option for heat transport:                          - set kmass = 0
c  option to perform mass transport only (no heat)     - set kmass = 1
c  option to simulate age transport only:              - set kmass = 2
c  option to simulate heat & 1-comp. mass transport :  - set kmass = 3
c  option to simulate heat & 2-comp. mass transport :  - set kmass = 4
c
c 
c  Lichtner formulation added feb 2017 (read alh,alv,... and kdisp=1)
c  Second component c2 added Dec 2022; assumed dilute, does not contribute to density
c  to do: add c2 to all output files, separate pumping concentration c1 from c2
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c    h - e - a - t    f - l - o - w    m - o - d - e - l
c    ---------------------------------------------------
c
c     code developed by:  j.w. molson and e.o. frind 
c
c       waterloo centre for groundwater research
c               university of waterloo
c
c copyright 1993 e.o.frind / j.w.molson
c duplication of this program or any part thereof without
c the express written consent of the authors is prohibited.
c ----------------------------------------------------------
c
c three-dimensional finite element numerical model for simulating
c transient, density-dependent groundwater flow, and thermal energy
c transport in a porous medium.
c
c temperature-dependent fluid viscosity and relative density
c
c conjugate gradient solver for flow,
c conjugate gradient solver for transport.
c
c linear, isoparametric quadrilateral elements.
c choice of numerical or direct integration 
c
c subroutine object files required for linking:
c
c      -  heat            ... heatflow main program
c
c      -  flow            ... 3d flow model source
c      -  trans           ... 3d thermal transport model source
c
c            - prism      ... grid generation
c            - sp3lin     ... basis functions
c            - mindex     ... boundary arrays, c.code.
c            - surf       ... boundary surface nodal influence areas
c            - veloc1,2       ... velocity calculations
c            - deform         ... grid deformation for wt mounding
c            - moment         ... spatial moments
c            - volume         ... domain volume
c            - exint0,exint1,exint2,exint3  ... exact integration for flow and trans,1,2 and 3
c            - gquad0,1,2,3                 ... numerical integration
c
c      - precg            ... pre-cond. conjugate gradient solver
c                             (and other associated matrix solver subroutines)
c ____________________________________________________________________________
c
c
c this routine is the primary driver file, it sets up arrays,
c contains the time and non-linear loops, flow and transport.
c
c ######################################################################
c
c variables:
c ----------
c
c        {u0}/{t0} - head/temperature at old time step
c        {u1}/{t1} - head/temperature at last iteration
c        {u2}/{t2} - head/temperature at most recent solution
c
c        {vx,vy,vz} - elemental velocities
c        {cx,cy,cz} - elemental hydraulic conductivities (0 c)
c
c
c ######################################################################
c
c array sizes ...
c ---------------
c     maxne = number of elements
c     maxnn = number of nodes
c     maxn  = number of degrees of freedom
c     maxnb = non-zero bandwidth in case of conjugate gradient solver
c     maxnbb= bandwidth for gauss solver (transport)
c     maxna = total number of non-zero matrix entries in condensed matri
c           = 14*n (approximately)
c     nf    = maximum number of nodes on one face
c             (also used for fap(nf) for element fractures in 2d plane read from random2d.plt
c     nw    = bandwidth
c     laa   = 3*n + na
c     maxfrac= maximum number of fractures 
c     maxfx = max number of internal heat flux elements
c     maxbt = maximum number of wells at which breakthrough curves will
c              be generated (a "well" can either be a vertical strip of nodes,
c              where the breakthrough concentrations will be averaged over 
c              the well screen, or it can be a single node)
c     maxit = max time steps
c     maxair = max # of air temperature points read in
c     mxs = max surface zones for surfat
c     maxbdy = max # head data points for transient watertable, and for transient bz (ex. snow thickness)
c     maxbzz = max # nodal zones for top bz, bzflx conditions, for file bzin.data
c     ****************************************************************************************************
c
c 3d: 40x31x24 nodes ... grid size for Borden validation run
c     parameter(maxn=1500000,maxnn=maxn,maxne=1500000,nf=252000,nw=15,
c    + maxss=500,maxna=14*maxnn,laa=3*maxnn+maxna,np1=maxnn+1,
c    + mxgx=7,mxgy=7,mxgz=5,maxfrac=300000,maxfx=maxne,
c    + maxbt=50,maxit=2000000)
c
c     Base:
c      parameter(maxn=202202,maxnn=maxn,maxne=100000,nf=25100,nw=15,
c     + maxss=500,maxna=14*maxnn,laa=3*maxnn+maxna,np1=maxnn+1,
c     + mxgx=7,mxgy=7,mxgz=5,maxfrac=500000,maxfx=maxne,maxnez=500,
c     + maxbt=15,maxit=20000,maxair=10000,mxs=5,maxbdy=51001,maxbzz=8,
c     + maxnx=2000) 
c     Jonas maxit:
       parameter(maxn=200000,maxnn=maxn,maxne=200000,nf=25100,nw=15,
     + maxss=500,maxna=14*maxnn,laa=3*maxnn+maxna,np1=maxnn+1,
     + mxgx=7,mxgy=7,mxgz=5,maxfrac=500000,maxfx=maxne,maxnez=500,
     + maxbt=11,maxit=200000,maxair=1000,mxs=5,maxbdy=51001,maxbzz=8,
     + maxnx=2000) 
c
c   for Prague:
c     parameter(maxn=11030000,maxnn=maxn,maxne=10000000,nf=1003000,
c    + nw=15,
c    + maxss=500,maxna=14*maxnn,laa=3*maxnn+maxna,np1=maxnn+1,
c    + mxgx=7,mxgy=7,mxgz=5,maxfrac=10000000,maxfx=maxne,
c    + maxbt=1,maxit=5000)    
c      
c
c 3d: grid for ATES example in Appendix B of User Guide
c      parameter(maxn=115351,maxnn=maxn,maxne=108000,nf=9000,nw=15,
c    + maxss=10000,maxna=14*maxnn,laa=3*maxnn+maxna,np1=maxnn+1,
c    + mxgx=5,mxgy=5,mxgz=5,maxfrac=10001,maxfx=maxne)
c
c 3d: grid for test example in Section 5 of the User Guide
c      parameter(maxn=10000,maxnn=maxn,maxne=10000,nf=1500,nw=15,
c     + maxss=50,maxna=14*maxnn,laa=3*maxnn+maxna,np1=maxnn+1,
c     + mxgx=5,mxgy=5,mxgz=5,maxfrac=10001,maxfx=maxne)
c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
c
c     constant dimensions
c     -------------------
      character*80 title1,title2,title3,title4,title5,date
      character*80 ktitle1,ktitle2
      dimension xi(8),yi(8),ag(9),hag(9)
      dimension itim(3),pt(5)
      ! TODO Missing kb4 definition
      dimension se(8,8),pe(8,8),inl(8),kb1(7),kb2(7),kb3(7),kb4(7)
      dimension knox(5),knoy(5),knoz(5)
      logical lunsat,lplot,lflow,l3dplt,lmom,lmom2,lplot43,lplot49
      logical lwtc,leak1,leak2,leak3,lss,update_flow
      logical lneg,lheat,lmass,lage,lzero,l3dvplt,lwtfgt0
      logical lsubtr1,lsubtr2,lkr1top,lbzin
      logical lf2d(maxne,6)
      real*8 start,finish
c
c     variable grid ...
c     -----------------
      dimension xlim(mxgx),ylim(mxgy),zlim(mxgz)
      dimension nlx(mxgx),nly(mxgy),nlz(mxgz)
      integer map(maxnn),mpa(maxnn),link(nf)
      dimension dz(nf),sumdz(nf),dztot(nf),ztotsx(maxnx)
      dimension inb(nf,4),inbl(nf,4),inbr(nf,4)
c
c     direct element matrices
c     -----------------------
      dimension ixx(8,8),iyy(8,8),izz(8,8),ixy(8,8),ixz(8,8),iyz(8,8),
     + ivx(8,8),ivy(8,8),ivz(8,8),itc(8,8),itl(8,8),igz(8)
      dimension iwx(8),iwy(8),iwz(8),dwx(8),dwy(8),dwz(8)
      dimension hxx(8,8),hyy(8,8),hzz(8,8),hxy(8,8),hxz(8,8),hyz(8,8),
     + hvx(8,8),hvy(8,8),hvz(8,8),htmf(8,8),htmt(8,8),hgz(8),
     + rg(8,8),rg2df(4,4)
c
c     arrays of size nn
c     -----------------
      real*8 pq(maxnn),tq(maxnn),cq(maxnn)
      real*8 x(maxnn),y(maxnn),z(maxnn)
      real*8 fc(maxnn),fb(maxnn),fst(maxnn),fs(maxnn)
      dimension fc0(maxnn)
      integer   ic(maxnn),lc(maxnn),kssn(maxss)
c
c     arrays for internal heat source
c     --------------------------------
      dimension flux(maxfx),fluxv(maxfx),ifl(maxfx),nflux(maxfx,8)
      real*4 zdts(maxnez),tdts(maxnez,maxit)
c
c     unknown vectors: u0,u1,u2 - pressure, t0,t1,t2 - temperature
c                      c0,c1,c2 - concentration 1
c                      d0,d1,d2 - concentration 2
c ------------------------------------------------------------------------------------------
      real*8    u0(maxnn),u1(maxnn),t0(maxnn),t1(maxnn)
      dimension u2(maxnn),t2(maxnn),fx(maxnn)
      dimension c0(maxnn),c1(maxnn),c2(maxnn)
      dimension d0(maxnn),d1(maxnn),d2(maxnn)
c
c     arrays of size ne, (ne,8)
c     -------------------------
      real*8  cx(maxne),cy(maxne),cz(maxne)
      real*8  por(maxne),por0(maxne),sw(maxne),rhob(maxne),ps(maxne)
      real*8  vt(maxne),exl(maxne),eyl(maxne),ezl(maxne)
      real*8  vx(maxne),vy(maxne),vz(maxne)
      real*8  tclm(maxne),cot0(maxne),cot1(maxne)
      real*8  tkl(maxne),cpsm(maxne),dzl(maxne)
      real*8  depc(maxne),vexp(maxne),depd(maxne)
      real*8  decay(maxne),decayc2(maxne),spn1(maxne),spn0(maxne)
      real*8  spe1(maxne),spe0(maxne),sdepd(maxne),smax(maxne)
      integer icnt(maxnn)
      real*8 pp(maxne),qq(maxne),ppn(maxnn),qqn(maxnn)
c
c     1D or 2D fracture elements
c     --------------------------
      dimension vlin(maxfrac),vx2d(maxfrac),vy2d(maxfrac)
      dimension ckl(maxfrac),ifracl(maxfrac),inline(4),xarea(maxfrac)
      dimension lvert(maxfrac),ifdim(maxfrac)
      integer ihx(4,4),ihy(4,4),ihh(4,4),ihc(4,4),ivx2(4,4),ivy2(4,4) 
      integer ihl(4,4),igz2(4),iwell(maxfrac)
      dimension hx(4,4),hy(4,4),hh(4,4),ht(4,4),wx(4,4),wy(4,4),hgz2(4)
      integer ifrac2dz1(maxfrac)
      dimension fap(nf)    !fracture apertures from random2d.plt ...use nf as proxy for nexy
c
c     breakthrough points single precision
c     ---------------------------------------
      dimension ixw(maxbt),iyw(maxbt),izw1(maxbt),izw2(maxbt)
      real*4 ubt(maxbt,maxit),ubh(maxbt,maxit)
      real*4 btime(maxit),tpkwell(maxbt,maxit),tminwell(maxbt,maxit)
      real*4 cpkwell(maxbt,maxit,2),cminwell(maxbt,maxit)
      real*4 spnwell(maxbt,maxit)
      real*4 xts(maxit),zts1(maxit),zts2(maxit)
c
c     masses
c     ---------
      real*4 cmass(maxit),smass(maxit),fluxin(maxit),fluxout(maxit)
c
c     air temperatures, bdy heads, etc. 
c     -------------------------------------
      dimension timeair(maxair),tempair(maxair)
      dimension thbdy(maxbdy),hbdy(maxbdy)
      dimension tbzbdy(maxbdy),bzbdy(maxbdy,maxbzz),
     +          bzflx(maxbdy,maxbzz),bzq(maxbdy,maxbzz),bztdiff(nf)
      dimension bz(nf),bzf(nf),bztemp(maxnx),bzftemp(maxnx),bzqn(nf)
      dimension nbzz(maxbdy),ix1bz(maxbdy,maxbzz),ix2bz(maxbdy,maxbzz)
c     dimension ix1bz(mxs),ix2bz(mxs),surfminbz(mxs),ampbz(mxs),
c    +          phasebz(mxs),cutoffbz(mxs)
c     dimension ix1sat(mxs),ix2sat(mxs),surfminsat(mxs),ampsat(mxs),
c    +          phasesat(mxs),cutoffsat(mxs)
c
!     TODO missing in4 definition
      integer  in(maxne,8),in2(maxne,8),in4(maxne,8),ivt(maxnn)
c
c     array of size nf
c     ----------------
      dimension fbff(nf),fbft(nf),fbfc(nf),ara(nf)
c
c     arrays a(na),iaa(na),ind(n+1),ib(n,nw),
c            aa(laa) ... where laa=3n+na
c            for conjugate gradient solver only
c     ------------------------------------------------
      dimension a(maxna),aa(laa)
      integer   iaa(maxna),ind(np1),ib(maxn,nw)
      integer kblck(maxnn)
c
c     duplicate transport arrays for boundary conditions, indeces:
c     ------------------------------------------------------------
      integer ib2(maxn,nw),ind2(np1),iaa2(maxna),ic2(maxnn),lc2(maxnn)
      integer ib3(maxn,nw),ind3(np1),iaa3(maxna),ic3(maxnn),lc3(maxnn)       !component 1
      integer ib4(maxn,nw),ind4(np1),iaa4(maxna),ic4(maxnn),lc4(maxnn)       !component 2

      real*8 fb2(maxnn),fc2(maxnn)
      real*8 fb3(maxnn),fc3(maxnn)
      real*8 fb4(maxnn),fc4(maxnn)
c
      abs(fn)=dabs(fn)
c
c     definition of gauss points for numerical integration:
c     -----------------------------------------------------
      data ag/-.5773502692,.5773502692,-.7745966692,0.0,.7745966692,
     +        -.8611363116,-.3399810436,.3399810436,.8611363116/
      data hag/1.000000000,1.000000000,.5555555556,.8888888889,
     +         .5555555556,.3478548451,.6521451549,.6521451549,
     +         .3478548451/
c
c     input/output list:
c     -------------------
      open(unit=5,file='smoker.data',status='old')
      open(unit=6,file='smoker.lst',status='unknown')
      open(unit=7,file='therm_in.data',status='old',iostat=ierr1)
      open(unit=8,file='smoker_txz.plt',status='unknown')
      open(unit=9,file='smoker_hxz.plt',status='unknown')
      open(unit=10,file='hinput.data',status='old',iostat=ierr2)
      open(unit=11,file='tbck.out',status='unknown')
      open(unit=12,file='veloc2dxz.out',status='unknown')
      open(unit=13,file='zdat.out',status='unknown')
      open(unit=14,file='smoker_txy_cxy.plt',status='unknown')
      open(unit=15,file='break.plt',status='unknown')
      open(unit=16,file='smoker_tyz.plt',status='unknown')
      open(unit=17,file='moment.out',status='unknown')
      open(unit=18,file='peaktrace.plt',status='unknown')
      open(unit=19,file='htot_flux.plt',status='unknown')
      open(unit=20,file='t1dz.plt',status='unknown')
      open(unit=21,file='smoker_hyz.plt',status='unknown')
      open(unit=22,file='vfracxz.out',status='unknown')
      open(unit=23,file='var.out',status='unknown')
      open(unit=24,file='param.out',status='unknown')
      open(unit=25,file='tdif.out',status='unknown')
      open(unit=26,file='t3d1_smoker.out',status='unknown')
      open(unit=27,file='t3d2_smoker.out',status='unknown')
      open(unit=28,file='t3d3_smoker.out',status='unknown')
      open(unit=29,file='t3d4_smoker.out',status='unknown')
      open(unit=30,file='t3d5_smoker.out',status='unknown')
      open(unit=31,file='tc1dx.plt',status='unknown')
      open(unit=32,file='veloc3d.out',status='unknown')
      open(unit=33,file='vfracxy.out',status='unknown')
      open(unit=34,file='veloc2dxy.out',status='unknown')
      open(unit=35,file='t3d_smokermovie.out',status='unknown')
      open(unit=36,file='veloc2dyz.out',status='unknown')
      open(unit=40,file='smoker_mass_fracs_pm.plt',status='unknown')
      open(unit=41,file='smoker_tpeak.plt',status='unknown')
      open(unit=42,file='smoker_1dxfront.plt',status='unknown')
      open(unit=43,file='t1dz_tdum.plt',status='unknown')
      open(unit=44,file='heatflux_2Dsurface.plt',status='unknown')
      open(unit=45,file='tair.data',status='old',iostat=ierr8)
      open(unit=46,file='smoker_hxy.plt',status='unknown')
      open(unit=47,file='fgen92.asc',status='old',iostat=ierr3)    
      open(unit=48,file='hbdy.data',status='old',iostat=ierr4)    
      open(unit=49,file='t1dx_tdum.plt',status='unknown')
      open(unit=51,file='kxz.plt',status='unknown')            
      open(unit=52,file='vfracxyelemset1.plt',status='unknown') 
      open(unit=522,file='vfracxydiscrete.plt',status='unknown') 
      open(unit=53,file='veloc2d_in.data',status='old',iostat=ierr5) 
      open(unit=54,file='heatflux_bndy.plt',status='unknown')
      open(unit=55,file='heatercolumnavg_1Dz.plt',status='unknown')
      open(unit=56,file='heatercolumnavg_vs_t.plt',status='unknown')
      open(unit=57,file='bzin.data',status='old',iostat=ierr6)
      open(unit=58,file='smoker_1dzfront.plt',status='unknown')
      open(unit=59,file='smoker_Pe_Co_xz.plt',status='unknown')


      open(unit=60,file='t1dx_hflux_xl.plt',status='unknown')
      open(unit=61,file='t1dx_bzq.plt',status='unknown')      
      open(unit=62,file='smoker_cxz.plt',status='unknown')
      open(unit=63,file='smoker_mass_particles.plt',status='unknown')


      open(unit=82,file='grid3d.plt',status='old',iostat=ierr7)
      open(unit=88,file='smoker-debug_1dx.plt')      
c
c     open(99) and open(98) used below, do not open here      
c
c
c     define variable grid: (normal size- workstation use)
c     set z-elements up to ground surface (transport domain)
c     -------------------------------------------------------
c
c     base matrices for prismatic elements
c     ------------------------------------
      data ixx/ 4,-4,-2, 2, 2,-2,-1, 1,-4, 4, 2,-2,-2, 2, 1,-1,
     +         -2, 2, 4,-4,-1, 1, 2,-2, 2,-2,-4, 4, 1,-1,-2, 2,
     +          2,-2,-1, 1, 4,-4,-2, 2,-2, 2, 1,-1,-4, 4, 2,-2,
     +         -1, 1, 2,-2,-2, 2, 4,-4, 1,-1,-2, 2, 2,-2,-4, 4/
      data iyy/ 4, 2,-2,-4, 2, 1,-1,-2, 2, 4,-4,-2, 1, 2,-2,-1,
     +         -2,-4, 4, 2,-1,-2, 2, 1,-4,-2, 2, 4,-2,-1, 1, 2,
     +          2, 1,-1,-2, 4, 2,-2,-4, 1, 2,-2,-1, 2, 4,-4,-2,
     +         -1,-2, 2, 1,-2,-4, 4, 2,-2,-1, 1, 2,-4,-2, 2, 4/
      data izz/ 4, 2, 1, 2,-4,-2,-1,-2, 2, 4, 2, 1,-2,-4,-2,-1,
     +          1, 2, 4, 2,-1,-2,-4,-2, 2, 1, 2, 4,-2,-1,-2,-4,
     +         -4,-2,-1,-2, 4, 2, 1, 2,-2,-4,-2,-1, 2, 4, 2, 1,
     +         -1,-2,-4,-2, 1, 2, 4, 2,-2,-1,-2,-4, 2, 1, 2, 4/
      data ixy/ 2, 0,-2, 0, 1, 0,-1, 0, 0,-2, 0, 2, 0,-1, 0, 1,
     +         -2, 0, 2, 0,-1, 0, 1, 0, 0, 2, 0,-2, 0, 1, 0,-1,
     +          1, 0,-1, 0, 2, 0,-2, 0, 0,-1, 0, 1, 0,-2, 0, 2,
     +         -1, 0, 1, 0,-2, 0, 2, 0, 0, 1, 0,-1, 0, 2, 0,-2/
      data ixz/ 2, 0, 0, 1, 0,-2,-1, 0, 0,-2,-1, 0, 2, 0, 0, 1,
     +          0,-1,-2, 0, 1, 0, 0, 2, 1, 0, 0, 2, 0,-1,-2, 0,
     +          0, 2, 1, 0,-2, 0, 0,-1,-2, 0, 0,-1, 0, 2, 1, 0,
     +         -1, 0, 0,-2, 0, 1, 2, 0, 0, 1, 2, 0,-1, 0, 0,-2/
      data iyz/ 2, 1, 0, 0, 0, 0,-1,-2, 1, 2, 0, 0, 0, 0,-2,-1,
     +          0, 0,-2,-1, 1, 2, 0, 0, 0, 0,-1,-2, 2, 1, 0, 0,
     +          0, 0, 1, 2,-2,-1, 0, 0, 0, 0, 2, 1,-1,-2, 0, 0,
     +         -1,-2, 0, 0, 0, 0, 2, 1,-2,-1, 0, 0, 0, 0, 1, 2/
      data ivx/-4,-4,-2,-2,-2,-2,-1,-1, 4, 4, 2, 2, 2, 2, 1, 1,
     +          2, 2, 4, 4, 1, 1, 2, 2,-2,-2,-4,-4,-1,-1,-2,-2,
     +         -2,-2,-1,-1,-4,-4,-2,-2, 2, 2, 1, 1, 4, 4, 2, 2,
     +          1, 1, 2, 2, 2, 2, 4, 4,-1,-1,-2,-2,-2,-2,-4,-4/
      data ivy/-4,-2,-2,-4,-2,-1,-1,-2,-2,-4,-4,-2,-1,-2,-2,-1,
     +          2, 4, 4, 2, 1, 2, 2, 1, 4, 2, 2, 4, 2, 1, 1, 2,
     +         -2,-1,-1,-2,-4,-2,-2,-4,-1,-2,-2,-1,-2,-4,-4,-2,
     +          1, 2, 2, 1, 2, 4, 4, 2, 2, 1, 1, 2, 4, 2, 2, 4/
      data ivz/-4,-2,-1,-2,-4,-2,-1,-2,-2,-4,-2,-1,-2,-4,-2,-1,
     +         -1,-2,-4,-2,-1,-2,-4,-2,-2,-1,-2,-4,-2,-1,-2,-4,
     +          4, 2, 1, 2, 4, 2, 1, 2, 2, 4, 2, 1, 2, 4, 2, 1,
     +          1, 2, 4, 2, 1, 2, 4, 2, 2, 1, 2, 4, 2, 1, 2, 4/
c
c     consistent and lumped matrices
c     -------------------------------
      data itc/ 8, 4, 2, 4, 4, 2, 1, 2, 4, 8, 4, 2, 2, 4, 2, 1,
     +          2, 4, 8, 4, 1, 2, 4, 2, 4, 2, 4, 8, 2, 1, 2, 4,
     +          4, 2, 1, 2, 8, 4, 2, 4, 2, 4, 2, 1, 4, 8, 4, 2,
     +          1, 2, 4, 2, 2, 4, 8, 4, 2, 1, 2, 4, 4, 2, 4, 8/
      data itl/27, 0, 0, 0, 0, 0, 0, 0, 0,27, 0, 0, 0, 0, 0, 0,
     +          0, 0,27, 0, 0, 0, 0, 0, 0, 0, 0,27, 0, 0, 0, 0,
     +          0, 0, 0, 0,27, 0, 0, 0, 0, 0, 0, 0, 0,27, 0, 0,
     +          0, 0, 0, 0, 0, 0,27, 0, 0, 0, 0, 0, 0, 0, 0,27/
c
c     gravity term matrix
c     -------------------
      data igz/-1,-1,-1,-1, 1, 1, 1, 1/
c
c     2D plane fracture element matrices
c     -----------------------------------
      data ihx/2,-2,-1,1,-2,2,1,-1,-1,1,2,-2,1,-1,-2,2/
      data ihy/2,1,-1,-2,1,2,-2,-1,-1,-2,2,1,-2,-1,1,2/
      data ihh/1,0,-1,0,0,-1,0,1,-1,0,1,0,0,1,0,-1/
      data ihc/4,2,1,2,2,4,2,1,1,2,4,2,2,1,2,4/
      data ihl/9,0,0,0,0,9,0,0,0,0,9,0,0,0,0,9/
      data ivx2/-2,-2,-1,-1,2,2,1,1,1,1,2,2,-1,-1,-2,-2/
      data ivy2/-2,-1,-1,-2,-1,-2,-2,-1,1,2,2,1,2,1,1,2/
      data igz2/-1,-1,+1,+1/
c
      do 4 i=1,4
      hgz2(i)=igz2(i)
      do 4 j=1,4
      hx(i,j)=ihx(i,j)
      hy(i,j)=ihy(i,j)
      hh(i,j)=ihh(i,j)
      ht(i,j)=ihc(i,j)
      wx(i,j)=ivx2(i,j)
    4 wy(i,j)=ivy2(i,j)
c
c     derivatives for direct velocities ...
c     --------------------------------------
      data iwx/-1,+1,+1,-1,-1,+1,+1,-1/
      data iwy/-1,-1,+1,+1,-1,-1,+1,+1/
      data iwz/-1,-1,-1,-1,+1,+1,+1,+1/
c
c     transfer to real element arrays
c     -------------------------------
      do 5 i=1,8
      do 5 j=1,8
      hxx(i,j)=ixx(i,j)
      hyy(i,j)=iyy(i,j)
      hzz(i,j)=izz(i,j)
      hxy(i,j)=ixy(i,j)
      hxz(i,j)=ixz(i,j)
      hyz(i,j)=iyz(i,j)
      hvx(i,j)=ivx(i,j)
      hvy(i,j)=ivy(i,j)
      hvz(i,j)=ivz(i,j)
      htmf(i,j)=itl(i,j)
      htmt(i,j)=itc(i,j)
      rg(i,j)=0.d0       !initialize rg if kint=1
      if(i.eq.j) then
       hgz(i)=igz(i)
       dwx(i)=iwx(i)
       dwy(i)=iwy(i)
       dwz(i)=iwz(i)
      endif
    5 continue

c
c     initialize cpu times
c     --------------------
      fpcgt=0.
      tpcgt=0.
      fast=0.
      tast=0.
      vcpu=0.

      iplot3d=0
      izonexz = 0
      izonexy = 0
      izoneyz = 0
      izone3d = 0
c
c *********************************************************************
c
c     input data ...
c     --------------
      read (5,10) title1,title2,date
   10 format (a)
      write (6,9) title1,title2,date
      write (*,9) title1,title2,date
    9 format (10x,60(1h*)//10x,'heatflow/smoker',/10x,
     +                         'john.molson@ggl.ulaval.ca',//
     + 10x,'3-d density-dependent flow and thermal transport model'/
     + 10x,'linear isoparametric elements'//10x,60(1h*)//10x,a//10x,
     + a/,10x,a/10x,60(1h*)/)
c  __________________________________________________________________
c
      read (5,*,err=8454) 
     + kprt,kcntrl,kwt,kint,kintv,kgo,ksat,kmass,ktair,krk,kdz,ksubtr,
     + kr1top,kfreec,krfrac
c     -----------------------------------------------------------------
c
      write(6,389) 
     + kprt,kcntrl,kwt,kint,kintv,kgo,ksat,kmass,ktair,krk,kdz,ksubtr,
     + kr1top,kfreec,krfrac
  389 format(/10x,'run-time options ...',
     +       /10x,'kprt,kcntrl,kwt,kint,kintv,kgo,ksat,kmass,ktair,',
     +            'krk,kdz,ksubtr,kr1top,kfreec,krfrac'
     +       /10x,15i3)
      lkr1top=.false.
      if(kr1top.eq.1) lkr1top=.true.
c
      if(kmass.eq.0) then 
      write(6,834)
 834  format(/10x,'heat transport option in effect')
      lheat=.true.
      lmass=.false.
      lage=.false.
      endif
      if(kmass.eq.1) then 
      write(6,835)
 835  format(/10x,'mass transport option in effect')
      lheat=.false.
      lmass=.true.
      lage=.false.
      endif
      if(kmass.eq.2) then 
      write(6,836)
 836  format(/10x,'age transport option in effect')
      lheat=.false.
      lmass=.false.
      lage=.true.
      endif
      if(kmass.eq.3) then 
      write(6,837)
 837  format(/10x,'heat & 1-component mass transport option in effect')
      lheat=.true.
      lmass=.true.
      lage=.false.
      endif
      if(kmass.eq.4) then 
      write(6,838)
 838  format(/10x,'heat & 2-component mass transport option in effect')
      lheat=.true.
      lmass=.true.
      lage=.false.
      endif
c
c     ksol = 1   cholesky solver flow
c          = 2   conjugate gradient solver flow
c
c     ks   = 1   centred scheme transport
c          = 2   second order scheme transport
c
c     kprt = 1   coordinates, incidences, boundary arrays printed
c          = 0   otherwise
c
c     kcntrl: program control options 
c     ---------------------------------------------------------
c          = 0 ... fully coupled transient flow and transport, stop.
c          = 1 ... steady flow only, no temperature or mass solution
c          = 2 ... no flow, transport with uniform {v} field, stop.
c          = 3 ... transient flow only, no transport.
c
c     kwt  = 1 ... option to iterate for watertable mounding
c          = 0 ... watertable fixed by grid geometry
c
c    kint  = 1 ... full element numerical integration
c          = 0 ... direct integration
c    
c    kintv = 1 ... numerical derivatives for velocity calculation
c          = 0 ... direct derivatives
c
c    kgo options without reading grid:
c    kgo   = 1 ... read initial temp, heads from restart file 7 (therm_in.data), use heads only
c    kgo   = 2 ... read initial temp, heads from restart file 7, use temp only
c    kgo   = 3 ... read initial temp, heads from restart file 7, use heads, T and C
c
c    kgo options with grid:
c    kgo   = 4 ... read new grid x,y,z, from grid3d.plt but do not use heads, T or C
c
c          = 0 ... do not use therm_in.data
c
c    ksat   = 0 ... unsaturated zone included for transport
c           = 1 ... saturated domain only
c
c    with ksat = 0, an unsaturated zone is included, and is defined by the
c                   topmost grid sub-interval in the z-direction
c                   i.e. it is located between your last two 'zlim' values
c                   and the flow velocities here are interpolated between the
c                   watertable and the ground surface where vx=0.
c                   In this case, when entering hydraulic conductivity values,
c                   be sure to only define those elements within the saturated
c                   flow domain.
c
c    with ksat = 1, the entire domain is saturated
c
c    kmass = 0 ... heat transport
c           =1 ... mass transport
c           =2 ... age transport, steady flow only
c           =3 ... heat & 1-component mass transport
c           =4 ... heat & 2-component mass transport
c
c    ktair  = 0 ... read air temperature cosine curve parameters
c           = 1 ... read air temperatures from tair.data
c
c      krk  = 1 ... read K field from fgen92.asc file 
c           = 2 ... read porosity field from fgen92.asc file (located after K field, use dogfld=T in fgen.gen)
c           = 0 ... get K field from input file 
c
c      kdz  = 1 ... allow grid deformation dz due to freeze/thaw (phase change volume only) 
c           = 0 ... no deformation
c
c    ksubtr = 1 ... new subroutine trans (no kappa, no R, needed for latent heat formulation)
c    ksubtr = 2 ... original subroutine trans (kappa, R) /R
c
c    kr1top = 0 .... normal conditions
c           = 1 .... assign top element row to always have kr = 1 (unfrozen). Note Wu is not affected.
c
c   kfreec = 0 ... leave boundary conditions alone
c            1 ... switch from fixed concentration to zero-gradient if flow is outwards
c
c  --------------------------------------------------------------------------------------------
c   hardwire:
c   ---------
      ksol=2
      ks=2
c
      if (ksol.eq.2) write (6,8)
    8 format (//10x,'conjugate gradient solver for flow'/)
      if(kint.eq.0) write (6,88)
   88 format(/10x,'direct integration for element matrices')
      if(kint.eq.1) write (6,188)
  188 format(/10x,'numerical integration for element matrices')
      if(kintv.eq.0) write (6,89)
   89 format(/10x,'direct derivatives for velocities')
      if(kintv.eq.1) write (6,189)
  189 format(/10x,'numerical derivatives for velocities')
c
      if(kgo.gt.0) write(6,289) kgo 
  289 format(/10x,'restart option in effect ... kgo = ',i4)
      if(kgo.le.3.and.kgo.gt.0) write(6,290) 
  290 format(/10x,'restart using therm_in.data ... ',/)
      if(kgo.eq.1) write(6,2789)
 2789 format(10x,'(using heads only ...)')
      if(kgo.eq.2) write(6,2790)
 2790 format(10x,'(using temperatures only ...)')
      if(kgo.eq.3) write(6,2791)
 2791 format(10x,'(using both heads, T and C ...)')
      if(kgo.eq.4) write(6,2792)
 2792 format(10x,'(reading grid, no external initial conditions ...')
c
c     ----------------------------
      if (ks.eq.1) write (6,26)
      if (ks.eq.2) write (6,27)
   26 format (/10x,'transport: standard scheme centered in time'/)
   27 format (/10x,'transport: advective term lagged, augmentation',
     +        /21x,           'matrix at midpoint in time'/)
      if (ks.eq.1) wp=.5
      if (ks.eq.2) wp=1.
      wa=0.
      wb=0.
      if (ks.eq.2) wa=0.5
      wp1=1.-wp
      wa1=1.-wa
      wb1=1.-wb
      write(6,21) wp
   21 format (/10x,'weighting for physical dispersion',f8.3/)
      if (ks.ne.1) write (6,22) wa,wb
   22 format ( 10x,'weighting for augmentation term  ',f8.3,
     +        /10x,'weighting for cauchy boundary    ',f8.3/)
c --------------------------------------------------------------
c --------------------------------------------------------------
      lunsat = (ksat.eq.0)
      if(lunsat) write(6,911)
      if(.not.lunsat) write(6,916)
  911 format(/10x,'unsaturated transport zone included')
  916 format(/10x,'saturated domain only',/)

      lsubtr1 = .false.
      lsubtr2 = .false.
      lsubtr1 = (ksubtr.eq.1)
      lsubtr2 = (ksubtr.eq.2)
      if(lsubtr1) then 
      write(6,912)
  912 format(/10x,'using transport formulation v1 for latent heat')
      elseif(lsubtr2) then
      write(6,9112)
 9112 format(/10x,'using transport formulation v2 with kappa,R')
      else
      write(6,9113)
 9113 format(/10x,'error, ksubtr for thermal transport must be 1',
     +            ' (Co on rhs for L.Heat) or 2 (kappa,R)',
     +       /10x,'program stopping ...')
      stop
      
      endif
c
c     read variable grid spacing
c     --------------------------
      read(5,*) ngx,ngy,ngz
      write(6,726) ngx,ngy,ngz
 726  format(/10x,'grid definition:',/10x,16('-'),/10x,
     +            'ngx,ngy,ngz: ',/10x,3i6)
      if(ngx.gt.mxgx.or.ngy.gt.mxgy.or.ngz.gt.mxgz) then
        write(6,883)
 883    format(/10x,'dimension error in grid definition ',
     +     /10x,'you have too many grid sub-intervals',
     +     /10x,'recompile with larger mxgx,mxgy,mxgz',
     +     /10x,'or reduce your number of sub-intervals (ngx,ngy,ngz)',
     +     /10x,'... program ending ')
        stop
      endif
c
      if(lunsat.and.ngz.le.1) then
      write(6,915)
 915  format(/10x,'error in grid definition ',
     +       /10x,'if unsat. zone included, then ngz must be > 1 ',
     +       /10x,'with the topmost layer representing the unsat. zone')
      stop
      endif
c
      read(5,*) (xlim(i),i=1,ngx)
      read(5,*) (ylim(i),i=1,ngy)
      read(5,*) (zlim(i),i=1,ngz)
 167  format(5f10.0)
      read(5,*) (nlx(i),i=1,ngx)
      read(5,*) (nly(i),i=1,ngy)
      read(5,*) (nlz(i),i=1,ngz)
   11 format (10i5)
c
      write(6,285)
 285  format(/10x,'Global grid intervals for element spacing:',
     +       /10x,42('-'))
      write(6,286) (xlim(i),i=1,ngx)
      write(6,287) (ylim(i),i=1,ngy)
      write(6,288) (zlim(i),i=1,ngz)
 286  format(10x,'x-dimension grid intervals: ',/(10x,6f10.2))
 287  format(10x,'y-dimension grid intervals: ',/(10x,6f10.2))
 288  format(10x,'z-dimension grid intervals: ',/(10x,6f10.2))
      write(6,710)
  710 format(/10x,'element distribution (nlx,nly,nlz): ',/)
      write(6,711) (nlx(i),i=1,ngx)
      write(6,712) (nly(i),i=1,ngy)
      write(6,713) (nlz(i),i=1,ngz)
 711  format(10x,'x-dimension: ',5i7)
 712  format(10x,'y-dimension: ',5i7)
 713  format(10x,'z-dimension: ',5i7)
c     ------------------------
      nex=0
      ney=0
      nez=0
      do 937 i=1,ngx
 937  nex=nex+nlx(i)
      do 938 i=1,ngy
 938  ney=ney+nly(i)
      maxz = ngz
      if(lunsat) maxz = ngz-1       !if unsat zone incuded, nez=flow grid
      do 939 i=1,maxz
 939  nez=nez+nlz(i)
c
      xl=xlim(ngx)
      yl=ylim(ngy)
      if(lunsat) zl=zlim(ngz-1)
      if(.not.lunsat) zl=zlim(ngz)
c
      nb=nw
      nbb=nw
c     -------------------------------
      write (6,14) nex,ney,nez,nb,nbb
   14 format (//10x,'flow grid size, elements (nex,ney,nez):',3i5/
     1          10x,'estimated bandwidth (nb,nbb):     ',2i5)
      nx=nex+1
      ny=ney+1
      nz=nez+1
      nn=nx*ny*nz
      ne=nex*ney*nez
      nexy=nex*ney
      nxy=nx*ny
      nxz=nx*nz
      nyz=ny*nz
      write (6,6)  nx,ny,nz,nn,ne
    6 format (/10x,'flow grid size, nx,ny,nz:',3x,3i5/
     1         10x,'total nodes in flow grid   ',i10/
     2         10x,'total elements in flow grid',i10/)
c
c default before testing for unsat zone grid ...      
      nzt=nz
      nnt = nx*ny*nzt
      net = ne
c
c    check limits
c    -------------
      if(ne.gt.maxne .or. nn.gt.maxn) then
      write(6,8333)
 8333 format(/10x,'error, ne>maxne or nn>maxn, check dimensions',
     +       /10x,'program stopping ...')
      stop
      endif

      if(nxy.gt.nf) then
      write(6,8433)
 8433 format(/10x,'error, nxy>nf, check dimensions',
     +       /10x,'program stopping ...')
      stop
      endif
c
c     initialize dzl,sw,xts
c     ----------------------
      do i=1,ne
      dzl(i)=0.
      sw(i)=0.
      depc(i)=0.d0
      depd(i)=0.d0
      decay(i)=0.d0
      decayc2(i)=0.d0
      enddo
      do i=1,maxit
      xts(i)=0.
      zts1(i)=0.
      zts2(i)=0.
      enddo
c
c
c     read number of deformable rows (nwtl), grid datum (datum), 
c     and concentration/density coefficient (gamma). If you set your
c     maximum normalized concentration = 1, then use gamma to set the
c     concentration/density relationship (see User Guide)      
c     note: assumes all tl layers are uniform thickness
c     ----------------------------------------------------
      read (5,*) nwtl,datum,gamma
c  16 format (i5,5f10.0)
      write (6,17) nwtl,datum,gamma
   17 format (/10x,'number of deformable layers: ',i5,
     +        /10x,'grid datum:                ',f10.4,
     +        /10x,'conc./density coefficient: ',f10.4)
c      if(gamma.eq.0.) write(6,8829)
c8829  format(/10x,'this is a linear mass transport run')
      if(gamma.gt.0.) write(6,8828)
8828  format(/10x,'gamma > 0; this is a nonlinear mass transport run')
c      
      if(nwtl.gt.nez) then
      write(6,778)
 778  format(/10x,'error - the number of deformable watertable layers',
     +      ' (nwtl) is greater than your number of vertical elements',
     +      ' program stopping ...',/)
      stop
      endif
c
c     read index nodes for breakthrough data
c     ----------------------------------------
      iw=0
 992  iw=iw+1
      read(5,*,err=3030) ixw(iw),iyw(iw),izw1(iw),izw2(iw),more
      write(6,8825) ixw(iw),iyw(iw),izw1(iw),izw2(iw)
8825  format(10x,'well index data for breakthrough well:',4i5)
      if( (ixw(iw)*iyw(iw)*izw1(iw)*izw2(iw)) .eq. 0) iw=iw-1
      if(more.gt.0) goto 992
      nwells=iw
      write(6,8826) nwells
 8826 format(/10x,'number of monitor wells ... ',i6)
      if(nwells.gt.maxbt) then
       write(6,3748) maxbt
 3748  format(/10x,'max # monitor wells = ',i8,
     +        /10x,'Error - # monitor wells > maxbt',
     +        /10x,'program stopping .... ')
      stop
      endif 
c
c     prismatic flow grid generation
c     -------------------------------
c
      kcall=0
      call prism (x,y,z,in,xl,yl,zl,nx,ny,nz,nex,ney,nez,nn,ne,
     +        maxnn,maxne,xlim,ylim,zlim,nlx,nly,nlz,ngx,ngy,ngz,
     +        mxgx,mxgy,mxgz,exl,eyl,ezl,map,mpa,inb,inbl,inbr,
     +        kcall,nf,lunsat,kgo)
c
c     read external grid
c     --------------------
      if(kgo.eq.4) then
      if(ierr7.gt.0) then
       write(6,1116)
 1116  format(/10x,'!!! error detected in grid3d.plt file !!! ',
     + /10x,'you are trying to read an external grid (grid3d.plt)'
     + /10x,'but this file does not exist',
     + /10x,'use kgrid=0 or create grid3d.plt file',
     + /10x,'program stopping ...')
       stop
      endif
      
      write(6,1122)
 1122 format(/10x,'kgo=4, reading external grid grid3d.plt ...')
      read(82,1117,err=10002) title1,title2,nz2,ny2,nx2
 1117 format(a,/,a,/7x,i8,4x,i8,4x,i8)
      if(nx2.ne.nx.or.ny2.ne.ny.or.(nz2.ne.nz.and.ksat.eq.1)) then
       write(6,1118) nx2,ny2,nz2
 1118  format(/10x,'!!! error in grid3d.plt !!!',
     +    /10x,'grid dimensions are not the same as in smoker.data',
     +    /10x,'nx2,ny2,nz2 = ',3i8,
     +    /10x,'program stopping ...')
       stop
      endif
c
c     check unsat zone case
c     ---------------------
      if(ksat.eq.0) then
      write(6,3489) nx,ny,nz,nlz(ngz),nx2,ny2,nz2
3489  format(/10x,'unsat check with grid3d.plt ...',
     +       /10x,'nx,ny,nz,nlz(ngz),nx2,ny2,nz2: ',7i6)
      if(nz2.ne.(nz+nlz(ngz))) then
      write(6,3488)
3488  format(/10x,'error detected with unsat option and grid3d.plt ...')
      stop     
      endif 
      endif

      nzt=nz2
      nnt = nx*ny*nzt
      call flush(6)
c      
      do i=1,nnt
      read(82,*,end=1119) x(i),y(i),z(i),dum
      enddo     
      endif 
c      
c      if (kprt.eq.0) go to 20
c      write (6,15)
c   15 format (/10x,'node coordinates'/
c     1         8x,'node',8x,'x',11x,'y',11x,'z'/)
c      do 18 i=1,nn/10
c   18 write (6,19) i,x(i),y(i),z(i)
c   19 format(2(i10,3f7.2))
c   20 continue
c
      write(6,199) (i,x(i),y(i),z(i),i=1,nz)
      write(6,299) (i,x(i),y(i),z(i),i=nn-nz+1,nn)
 199  format(/10x,'vertical strip at (x,y)=(0,0) x,y,z coords:',/,
     + (2(i9,3e10.3)))
 299  format(/10x,'vertical strip at (x,y)=(xl,yl) x,y,z coords:',/,
     + (2(i9,3e10.3)))
c
c     1D line or 2D plane fractures
c     1D: read fracture element position, full aperture (cylinder assumed)
c     for a 2D plane fracture, read 2b = full aperture
c     vertical line elements are located at incidences 1 and 5 of element lf
c     horizontal line elements are located at incidences 1 and 2 
c     faces of 2D plane fractures defined 1-6
c     read idim (dimension of fracture = 1 or 2)
c     1D: xarea=pi*R^2
c     2D: xarea=2b
c     -----------------------------------------------------------------
c   hardwire for Cambridge: identify sealed well elements - make sure these are the same as low-K well elements
c   (for well elements from 46-47, set ixw2=48 to flag yz fractures on face 1)
c     ixw1=46;ixw2=48;iyw1=11;iyw2=12;nswe=0
      ixw1=1;ixw2=1;iyw1=1;iyw2=1;nswe=0
      do i=ixw1,ixw2
      do j=iyw1,iyw2
      do k=1,nez
      nswe=nswe+1
      l3d = (k-1)*nex*ney+(j-1)*nex+i
c     iwell(nswe) = l3d
      iwell(nswe) = 0                 !deactivate hardwire
      enddo
      enddo
      enddo


      open(99,file='smoker_fracgeo_xz.plt',status='unknown')
      open(98,file='smoker_fracgeo_xy.plt',status='unknown')

      do i=1,ne      
      do k=1,6
      lf2d(i,k)=.false.     !initialize all element surfaces (need to do same for 1D fractures ...)
      enddo
      enddo

      grav = 10. 
      tbk1 = 4.
      if(kmass.eq.1.or.kmass.eq.2.or.kmass.eq.3.or.kmass.eq.4) tbk1=0.
      tbk2 = 0.
      nfrac = 0
      nskip=0
      zk1 = den(tbk1,lheat,lmass,lage)*grav 
     +                    / ( 160.*rvisc(tbk2,lheat,lmass,lage)/86400.)

  727 read(5,*,err=6661) 
     +     ilx1,ilx2,ily1,ily2,ilz1,ilz2,idim,iloc,aperture,more

      if(ilx1*ilx2*ily1*ily2*ilz1*ilz2.eq.0) goto 8772
      
      if(idim.eq.1) xareax = 4.d0 * atan(1.d0) * (aperture/2.)**2
      if(idim.eq.2) xareax = aperture
      if(idim.eq.1) ckll = zk1*((aperture/2.d0)**2)/8.d0
      if(idim.eq.2) ckll = zk1*((aperture**2))/12.d0

      if(idim.eq.1) write(6,882) 
     + ilx1,ilx2,ily1,ily2,ilz1,ilz2,iloc,ckll,aperture,xareax
 882  format(/10x,'1D line fracture:',
     +       /10x,'element indeces:                ... ',6i5,
     +       /10x,'orientation (1,2,3 = z,x,y)     ... ',i12,
     +       /10x,'hyd. conductivity (m/s@ 0C) ... ',e12.4,
     +       /10x,'aperture (=pipe diameter)   ... ',e12.4,
     +       /10x,'x-sectional area            ... ',e12.4,/)
c
      if(idim.eq.2) write(6,881)
     + ilx1,ilx2,ily1,ily2,ilz1,ilz2,iloc,ckll,aperture,xareax
 881  format(/10x,'2D plane fracture:',
     +       /10x,'element indeces:                ... ',6i5,
     +       /10x,'face (1 = left face)            ... ',i12,
     +       /10x,'hyd. conductivity (m/s@ 0C) ... ',e12.4,
     +       /10x,'aperture (=2b)              ... ',e12.4,
     +       /10x,'x-sectional area (=2b)      ... ',e12.4,/)
c

      if(lunsat.and. (ilz1.gt.nez .or. ilz2.gt.nez))  then
      write(6,989)
 989  format(/10x,'fracture elements found within unsaturated zone',
     +       /10x,'Rerun with ksat = 1 or remove fracs from unsat zone',
     +       /10x,'program stopping ...')
      stop
      endif

      if( ilx1.gt.nex.or.ily1.gt.ney.or.ilz1.gt.nez
     + .or.ilx2.gt.nex.or.ily2.gt.ney.or.ilz2.gt.nez) then
              write(6,8375) ilx1,ilx2,ily1,ily2,ilz1,ilz2
 8375         format(/10x,'!! warning, fracture grid index truncated',
     +               /10x,'input indeces: ',6i6)
              endif

      if(ilx2.gt.nex) ilx2=nex
      if(ily2.gt.ney) ily2=ney
      if(ilz2.gt.nez) ilz2=nez
      if(ilx1.gt.nex) ilx1=nex
      if(ily1.gt.ney) ily1=ney
      if(ilz1.gt.nez) ilz1=nez
c
c     assign fracture element properties
c     ----------------------------------
      if(xareax.gt.0.) then
      do  i=ilx1,ilx2
      do  j=ily1,ily2
      do  k=ilz1,ilz2
      l3d=(k-1)*nex*ney + (j-1)*nex + i

c skip if already done  !warning: this is very slow with lots of fractures ....
c      do i2 = 1, nfrac
c      if(ifracl(i2).eq.l3d .and. ifdim(i2).eq.idim 
c     +    .and. lvert(i2).eq.iloc) then
c                 nskip=nskip+1 
c                 goto 3001
c                 endif
c      enddo

      if(lf2d(l3d,iloc)) then      !skip if already done (only checks 2d elements ...)
      nskip = nskip + 1
      goto 3001      
      endif

c hardwire: skip if fracture element lies within sealed well
      do iwe=1,nswe
      if(l3d.eq.iwell(iwe).and.idim.eq.2) then
         if(iloc.eq.1) goto 3001
         if(iloc.eq.5.and.i.ne.ixw2) goto 3001  !allow xy fracture on bottom of outside column
         endif
      enddo
      
      nfrac  = nfrac +1
      ckl(nfrac) = ckll 
      xarea(nfrac) = xareax
      lvert(nfrac) = iloc
      ifracl(nfrac) = l3d
      ifdim(nfrac) = idim
      lf2d(l3d,iloc) = .true.
 3001 continue
      enddo
      enddo
      enddo
      endif

       if(nfrac.gt.maxfrac) then
       write(6,8820) nfrac,nskip
 8820  format(/10x,'Error: nfrac > maxfrac, nfrac =',i10,'nskip = ',i10,
     +        /10x,'program stopping ...')
       stop
       endif

       nfracxy=0
       nfracyz=0
       nfracxz=0
       do i=1,nfrac
       if(ifdim(i).eq.2 .and. (lvert(i).eq.5 .or. lvert(i).eq.6) ) 
     +                                              nfracxy=nfracxy+1
       if(ifdim(i).eq.2 .and. (lvert(i).eq.1 .or. lvert(i).eq.2) ) 
     +                                              nfracyz=nfracyz+1
       if(ifdim(i).eq.2 .and. (lvert(i).eq.3 .or. lvert(i).eq.4) ) 
     +                                              nfracxz=nfracxz+1
       enddo
       write(6,8822) nfracxy,nfracyz,nfracxz
 8822  format(/10x,'number of 2D xy plane fractures: ',i7,
     +        /10x,'number of 2D yz plane fractures: ',i7,
     +        /10x,'number of 2D xz plane fractures: ',i7)
c
c     write fractures to geometry files for plotting in xz and xy sections
c ----------------------------------------------------------------------------
c     l1=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1      !x start at y1
c     l2=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx2      !x end   at y1
c     l3=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1    !y start
c     l4=(ilz1-1)*nex*ney + (ily2-1)*nex + ilx1    !y end
c     l5=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1    !z start
c     l6=(ilz2-1)*nex*ney + (ily1-1)*nex + ilx1    !z end
c
c    write xy plane fractures to geometry file in xz section
c    here, in() are still flow grid incidences
c    -----------------------------------------
      if(idim.eq.2 .and.iloc.eq.5) then
              l1=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1      !x start at y1
              l2=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx2      !x end   at y1
              x1=x(in(l1,1))
              z1=z(in(l1,1))
              x2=x(in(l2,2))
              z2=z(in(l2,2))
              write(99,341) x1,z1,x2,z2
c             write(99,*) 'ok1'
              endif

      if(idim.eq.2 .and.iloc.eq.6) then
              l1=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1      !x start at y1
              l2=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx2      !x end   at y1
              x1=x(in(l1,5))
              z1=z(in(l1,5))
              x2=x(in(l2,6))
              z2=z(in(l2,6))
              write(99,341) x1,z1,x2,z2
c             write(99,*) 'ok2'
              endif
c
c    write yz plane fractures to geometry file in xz section

      if(idim.eq.2.and.iloc.eq.1) then
              l1=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1      !z start at y1
              l2=(ilz2-1)*nex*ney + (ily1-1)*nex + ilx1      !z end   at y1
              x1=x(in(l1,1))
              z1=z(in(l1,1))
              x2=x(in(l2,5))
              z2=z(in(l2,5))
              write(99,341) x1,z1,x2,z2
c             write(99,5857) ilx1,ilx2,ily1,ily2,ilz1,ilz2 
c5857         format('ok1',6i9)             
              endif

      if(idim.eq.2.and.iloc.eq.2) then
              l1=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1      !z start at y1
              l2=(ilz2-1)*nex*ney + (ily1-1)*nex + ilx1      !z end   at y1              
              x1=x(in(l1,2))
              z1=z(in(l1,2))
              x2=x(in(l2,6))
              z2=z(in(l2,6))
              write(99,341) x1,z1,x2,z2
c             write(99,*) 'ok3'
              endif
c
c     write 1D vertical line fractures to geometry file for plotting in xz section
c     (need to add 1D horizontal ...)
      if(idim.eq.1.and.iloc.eq.1) then
              l1=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1      !z start at y1
              l2=(ilz2-1)*nex*ney + (ily1-1)*nex + ilx1      !z end   at y1              
              x1=x(in(l1,1))
              z1=z(in(l1,1))
              x2=x(in(l2,5))
              z2=z(in(l2,5))
              write(99,341) x1,z1,x2,z2
c             write(99,*) 'ok4' 
              endif
c
      if(idim.eq.1.and.iloc.eq.2) then
              l1=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1      !z start at y1
              l2=(ilz2-1)*nex*ney + (ily1-1)*nex + ilx1      !z end   at y1                            
              x1=x(in(l1,2))
              z1=z(in(l1,2))
              x2=x(in(l2,6))
              z2=z(in(l2,6))
              write(99,341) x1,z1,x2,z2 
              write(99,*) 'ok5'
              endif 
c
c     write(99,341) x1,z1,x2,z2   
 341  format('geometry x=0,y=0,t=line,c=black,cs=grid,lt=0.2',
     +       /,'1',/,'2',/,2f10.3,/,2f10.3)
c
c     write yz plane fractures to geometry file for plotting in xy section
      if(idim.eq.2.and.iloc.eq.1) then
              l1=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1      !y start at z1
              l2=(ilz1-1)*nex*ney + (ily2-1)*nex + ilx1      !y end   at z1          
              x1=x(in(l1,5))
              y1=y(in(l1,5))
              x2=x(in(l2,8))
              y2=y(in(l2,8))
              write(98,341) x1,y1,x2,y2   
              endif
c
c     write xz plane fractures to geometry file for plotting in xy section
      if(idim.eq.2.and.iloc.eq.3) then
              l1=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1      !x start at y1
              l2=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx2      !x end   at y1                        
              x1=x(in(l1,1))
              y1=y(in(l1,1))
              x2=x(in(l2,2))
              y2=y(in(l2,2))
              write(98,341) x1,y1,x2,y2   
              endif
c
c     write 1D line fractures to geometry file for plotting in xy section
      if(idim.eq.1.and.iloc.eq.1) then
              l1=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1      !x start at y1
              l2=(ilz1-1)*nex*ney + (ily1-1)*nex + ilx1      !x end   at y1                             
              x1=x(in(l1,1))
              y1=y(in(l2,1))
              write(98,341) x1,y1,x1,y1                 
              endif

8772  if(more.gt.0) goto 727                !get next fracture input line
c
      call flush(99)
      close(99)
      call flush(98)
      close(98)

      write(6,7382) nskip
 7382 format(/10x,'found ',i8,' duplicate fracture elements')


c --------------------------------
c
c     read 2D planar fracture distribution from random2d.plt
c     later: extend to 3D, add error checks
c     -------------------------------------

      if(krfrac.eq.1) then
      iplane=3                        !hardwire plane to 4rd plane (xy)
      izn = (nez/2) + 1                !hardwire place in midlle of z
      if(iplane.eq.3) neplane2d=nex*ney
      open(unit=64,file='random2d.plt',status='old',iostat=ierr9)
      if(ierr9.gt.0) then
      write(6,1126)
 1126  format(/10x,'!!! error detected in random2d.plt !!! ',
     + /10x,'you are trying to read an external aperture field'
     + /10x,'but this file does not exist',
     + /10x,'use krfrac=0 or create random2d.plt file',
     + /10x,'program stopping ...')
       stop
      endif

      read(64,10,err=1137) title1,title2,title3
      ! TODO Wrong pharentesis 
      read(64,*,err=1137) (dum1,dum2,fap(i),dum3,dum4,i=1,neplane2d)
      close(64)

      do i=1,neplane2d
      nfrac=nfrac+1
      ie3d = (izn-1)*neplane2d + i             !only for xy plane
      ifracl(nfrac) = ie3d
      ifdim(nfrac) = 2
      lvert(nfrac) = 5                         !on bottom of elements
      xarea(nfrac) = fap(i)
      zk1 = den(tbk1,lheat,lmass,lage)*grav 
     +                    / ( 160.*rvisc(tbk2,lheat,lmass,lage)/86400.)
      ckll = zk1*((fap(i)**2))/12.d0
      ckl(nfrac) = ckll
      enddo
      endif

      write(6,827) nfrac
 827  format(/10x,'number of fracture elements ... ',i7,/,
     +    /10x,'summary: listing of first 100 fracture elements:',
     +    /10x,'l       ck      xarea       lvert    ifracl    ifdim',
     +    /10x,'       (m/s)    (m^2)       (face)   (3D #)       '/)
      if(nfrac.gt.0) then 
      write(6,779) (l,ckl(l),xarea(l),lvert(l),ifracl(l),ifdim(l),
     +                                         l=1,min(100,nfrac))
 779  format(1x,i10,2e12.4,3i9)
      
      write(6,828) 
 828  format(/10x,'summary: listing of last 100 fracture elements:',
     +    /10x,'l       ck      xarea       lvert    ifracl    ifdim',
     +    /10x,'       (m/s)    (m^2)       (face)   (3D #)       '/)
      write(6,779) (l,ckl(l),xarea(l),lvert(l),ifracl(l),ifdim(l),
     +                                   l=max(1,nfrac-100),nfrac)
      endif
      call flush(6)
c
c ==============================================================
c
c     set up boundary arrays and condensation code for flow ...
c     ---------------------------------------------------------
c
      ktype = 0
      call mindex(maxnn,maxne,maxn,maxna,laa,np1,nn,ne,nf,n,
     + nx,ny,nz,in,ic,fb,fc,lc,fbff,ib,a,iaa,ind,nb,kprt,ktype,
     + nc,m,nb1,x,y,z,ara,datum,leak1,por,
     + wh,ampcos,ampsin,wpc,wps,slope,phasewt,map,kb1,lwtfgt0,ffc)
      call flush(6)
c
c     initialize and save initial fc arrays
c     --------------------------------------
      do i=1,nn
      fc0(i)=fc(i)
      enddo
c
c     write(6,9883) (fc(i),i=nz,nn,nz)
c
c ___________________________________________________________________
c
c     hydraulic conductivities (elemental at zero degrees celsius)
c     initial backgound conductivies must be set here,
c     add heterogeneities by layer in this read block, 
c     or using grid indeces (in next read statement) to define 'boxes' of K
c     local element axes follow principal directions
c     continue with 'more' = +1
c     end with      'more' = -1
c ----------------------------------------------------------------
c
      write(6,101)
  101 format(/10x,'specified hydraulic conductivity etc. @ 0 celsius:',
     +   /10x,45('-'),
     +   /10x,'element range        ckx          cky          ckz',
     + 7x,'tcx(sol)   cps(sol)    ps(sol)   porosity   saturation',
     + 1x,' depc_coef. depd_coef.  v_exp(n)    smax      rhob',
     + 1x,'      p     q')
   35 read (5,*,err=6662) j1,j2,ckx,cky,ckz,tcx,cps,ps0,theta,satw,
     +                    depc1,depd1,vexp1,smax0,pp0,qq0,more
c  33 format (2i9,6f10.0)
      rhob0 = (1.0-theta)*ps0                                           !!! calculate bulk density
c     rhop = ps0
      write(6,102) j1,j2,ckx,cky,ckz,tcx,cps,ps0,theta,satw,
     +              depc1,depd1,vexp1,smax0,rhob0,pp0,qq0
  102 format(7x,2i8,2x,3(e10.3,3x),10f11.3,2f6.2)
c
      if( j1.le.0.or.j2.le.0 .or. (j2.gt.ne.and. (.not.lunsat) )) then 
      write(6,924)
 924  format(/10x,'!! warning !! element range for hyd. conductivity',
     +     /10x,'is out of bounds; check data file and grid dimensions',
     +     /10x,'program stopping',/)
      stop
      endif
      if( j1.gt.ne.or.j2.gt.ne) then 
      write(6,9224)
9224  format(/10x,'!! warning !! K above wt. will not be used',
     +    /10x,'por,sw,tclm,cps*ps will be used if unsat zone exists')
      endif
      if(smax0.le.0.000) then 
      write(6,9235)
9235  format(/10x,'!! warning !! smax <=0 ',
     +    /10x,'program stopping')
      stop
      endif

c     if(rhob.le.0.000) then 
c     write(6,9236)
c9236  format(/10x,'!! warning !! bulk density rhob <=0 ',
c     +    /10x,'program stopping')
c      stop
c      endif
c
c    transfer to vectors.
c    smax is now a vector
c    ---------------------------
      do j=j1,j2
      por(j) = theta
      por0(j) = theta        !initial porosity
      sw(j)  = satw
      tclm(j) = tcx
      cpsm(j)  = cps*ps0
      depc(j) = depc1
      depd(j) = depd1
      vexp(j) = vexp1
      rhob(j) = rhob0 
      ps(j) = ps0
      pp(j) = pp0
      qq(j) = qq0 
      smax(j) = smax0 
      if(j.le.ne) then              !assign K's only for flow grid
      cx(j) = ckx
      cy(j) = cky
      cz(j) = ckz
      endif
      enddo
      if(more.gt.0) go to 35
   30 continue
c
c     override with indexed elements ... end with -1
c     -----------------------------------------------
      write(6,914)
 914  format(/10x,'indexed hydraulic k:',
     +       /10x,'--------------------')
 943  read(5,*,err=3041)
     +        i1,i2,j1,j2,k1,k2,ckx,cky,ckz,tcx,cps,ps0,theta,satw,
     + depc1,depd1,vexp1,smax0,pp0,qq0,more

      rhob0 = (1.0-theta)*ps0   !!! calculate bulk density
c     rhop = ps0

      write(6,8166) i1,i2,j1,j2,k1,k2
 8166 format(/10x,'element indices ix1-ix2,iy1-iy2,iz1-iz2 :',6i5)
      write(6,8167) ckx,cky,ckz,tcx,cps,ps0,theta,satw,
     +                      depc1,depd1,vexp1,smax0,rhob0,pp0,qq0
! TODO Wrong format 6i5 as second element
 8167 format(/10x,/10x,/25x,3(e10.3,3x),10f11.3,2f6.1)
c
      if(i1*i2*j1*j2*k1*k2.le.0) then
      write(6,8309)
 8309 format(/10x,'Indexed elements for K contain 0, data ignored ... ')
      goto 945
      endif
      if(i2.gt.nex.or.j2.gt.ney.or.
     +                    (k2.gt.nez.and.(.not.lunsat))) then
      write(6,8310)
 8310 format(/10x,'Indexed elements for K: bad entries, stopping ... ')
      stop
      endif
c
      do 944 k=k1,k2
      do 944 j=j1,j2
      do 944 i=i1,i2
      iel=(k-1)*nexy+(j-1)*nex+i
c      if(iel.lt.0 .or. iel.gt.(ne-nlz(ngz)*nexy)) goto 944   !bad, ne is flow grid
      if(iel.le.0 .or. (iel.gt.ne.and.(.not.lunsat))) goto 944
      por(iel)=theta
      por0(iel)=theta
      sw(iel) = satw
      tclm(iel)=tcx
      cpsm(iel)  = cps*ps0
      depc(iel) = depc1
      depd(iel) = depd1
      vexp(iel) = vexp1
      rhob(iel) = rhob0
      ps(iel) = ps0                !bug fix iel not j dec 2021
      pp(iel) = pp0
      qq(iel) = qq0 
      smax(iel) = smax0 
      if(iel.le.ne) then              !assign K's only for flow grid      
      cx(iel)=ckx
      cy(iel)=cky
      cz(iel)=ckz
      endif
 944  continue
      if(more.gt.0) goto 943

 945  continue
c
c     get nodal ppn(),qqn() arrays needed for printing
c     -------------------------------------------------
      do i=1,nnt
      icnt(i)=0
      ppn(i)=0.
      qqn(i)=0.
      enddo
      do l=1,net
      do j=1,8
      node = in(l,j)
      icnt(node) = icnt(node) + 1
      ppn(node) = ppn(node) + pp(l)
      qqn(node) = qqn(node) + qq(l)
      enddo
      enddo
      do i=1,nn
      ppn(i) = ppn(i)/float(icnt(i))
      qqn(i) = qqn(i)/float(icnt(i))
      enddo
      ppmax= -999.
      ppmin= +999.
      qqmax= -999.
      qqmin= +999.
      do i=1,nn
      if(ppn(i).gt.ppmax) ppmax=ppn(i)
      if(ppn(i).lt.ppmin) ppmin=ppn(i)
      if(qqn(i).gt.qqmax) qqmax=qqn(i)
      if(qqn(i).lt.qqmin) qqmin=qqn(i)
      enddo
      write(6,7788) ppmax,ppmin,qqmax,qqmin
7788  format(/10x,'unfrozen moisture curve parameters:',
     +       /10x,'pmax,pmin: ',2f6.1,
     +       /10x,'qmax,qmin: ',2f6.1)
c
c  
c ==================================================================
c
c     read random K & porosity field from fgen92.asc file
c     will over-ride K's and porosities read from smoker.data input file
c     --------------------------------------------------------------------
      if(krk.ge.1) then

      if(ierr3.gt.0) then
           write(6,8824)
 8824      format(/10x,'error opening fgen92.asc ...',
     +            /10x,'this file does not exist ...',
     +            /10x,'program stopping. ')
           stop
      endif
c
c     read K's (H field in fgen)
c     --------------------------
      write(6,8821) 
 8821 format(/10x,'reading K field from fgen92.asc ... ') 
      cmax=-999.
      cmin=999.
      read(47,10,err=3033) ktitle1
      read(47,*,err=3033) nxfg,nyfg,nzfg
      nefg=nxfg*nyfg*nzfg
      if(nefg.ne.ne) then
           write(6,8814) nefg
 8814      format(/10x,'error reading fgen92.asc (K) ...',
     +            /10x,'fgen elements (',i7,') <> ne in smoker ...',
     +            /10x,'program stopping ... ')
           stop
      endif

      do 7733 i=1,ne
      read(47,*,err=3033) cxx
      cx(i)=exp(cxx)
      cy(i)=exp(cxx)
      cz(i)=exp(cxx)
      if(cx(i).gt.cmax) cmax = cx(i)
      if(cy(i).gt.cmax) cmax = cy(i)
      if(cz(i).gt.cmax) cmax = cz(i)
      if(cx(i).lt.cmin) cmin = cx(i)
      if(cy(i).lt.cmin) cmin = cy(i)
      if(cz(i).lt.cmin) cmin = cz(i)
 7733  continue
c
      write(6,8823) cmax,cmin
 8823 format(/10x,'max, min conductivities from fgen input file ... ',
     +       /10x,'cmax m/s... ',e12.4,'  cmin m/s... ',e12.4)
c
c     read porosities (G field in fgen)
c     ----------------------------------
      if(krk.eq.2) then
      pormax=-999.
      pormin=999.
      read(47,10,end=3042) ktitle1
      read(47,*,err=3043,end=3043) nxfg,nyfg,nzfg
      nefg=nxfg*nyfg*nzfg
      if(nefg.ne.ne) then
           write(6,8815) nefg
 8815      format(/10x,'error reading fgen92.asc (porosities) ...',
     +            /10x,'fgen elements (',i7,') <> ne in smoker ...',
     +            /10x,'program stopping ... ')
           stop
      endif
      write(6,8811) 
 8811 format(/10x,'reading porosity field from fgen92.asc ... ')       
      do  i=1,ne
      read(47,*,err=3043,end=3043) porxx
      por(i)=exp(porxx)
      if(por(i).gt.pormax) pormax = por(i)
      if(por(i).lt.pormin) pormin = por(i)
      enddo
c
      write(6,8813) pormax,pormin
 8813 format(/10x,'max, min porosities from fgen input file ... ',
     +       /10x,'pormax... ',e12.4,'  pormin... ',e12.4)

      goto 9901
3042  continue
      write(6,5757) 
 5757 format(/10x,'end of file reading porosities: fgen92.asc file ...'
     +  /10x,'porosities do not exist (check dogfld: fgen92.gen ?) ...',
     +  /10x,'continuing with porosities from smoker.data ...')

 9901 continue
      endif     
      endif 
c
c     check conductivities and porosities
c     -------------------------------------
      ierrk=0
      ierrp=0
      do 955 i=1,ne
      if(cx(i).le.0..or.cy(i).le.0..or.cz(i).le.0.) ierrk=1
      if(por(i).le.0. .or. por(i).gt.1.) ierrp=1
  955 continue
      if(ierrk.eq.1) then
      write(6,825)
  825 format(/10x,'error, hydraulic conductivity for all elements',
     +       /10x,'must be > 0 ; check grid size, and element ranges',
     +       /10x,'for hydraulic conductivity distribution',
     +       /10x,'program stopping ... ')
      stop
      else
      write(6,843)
 843  format(/10x,'hydraulic conductivities are all positive ... ',
     +       /10x,'continuing ... ')
      endif

      if(ierrp.eq.1) then
      write(6,845)
  845 format(/10x,'error, porosity out of range, stopping ... ')
      stop
      else
      write(6,8465)
 8465 format(/10x,'porosities all positive, continuing ... ')
      endif


c
c     print conductivities
c     ---------------------
c      if(kprt.ne.0) then
c      write (6,36)
c   36 format (//10x,'element data'/10x,12(1h-)/10x,'element',15x,
c     1'incidence',25x,'conductivity'/)
c      write (6,37) (i,(in(i,j),j=1,8),cx(i),cy(i),cz(i),i=1,ne)
c   37 format (10x,i4,2x,8i5,4x,3f8.0)
c      endif
c

c     read storage term
c     ss=porosity for top layer
c     ---------------------------
      read (5,*,err=3034) ss
      write(6,125) ss
  125 format(/10x,'storage term ... ',e12.4)
c
c     initial conditions for flow ...
c     enter init =0 for uniform initial flow condition,
c                =1 to read the initial condition from data file #10.
c     -----------------------------------------------------------------
c
      read(5,*,err=3135) init,h0
  126 format(i5,f10.0)
      if(init.eq.1) goto 128
      write(6,127) h0
  127 format(/10x,'uniform initial condition for flow; h0= ',e12.4)
      do 120 i=1,nn
      u0(i) = h0
      u1(i) = h0
      u2(i) = h0
  120 continue
      goto 129
c
c     read head initial condition from unit 10 (hinput.dat)
c     hinput.data can be derived from a hxz.data file
c     just rename an old hxz.out file and use init=1
c     initial heads will be uniform in the y-direction
c     --------------------------------------------------------
  128 continue
      if(ierr2.gt.0) goto 1010
      read(10,10,end=1009) title3
      do 813 i=1,nx
      i1= (i-1)*nyz 
      read(10,809,end=1009) (xdum,zdum,u0(i1+k),k=1,nz)
c
c     extend in y-direction
c     ---------------------
      do 814 jj=2,ny
      do 814 kk=1,nz
      node=i1 + (jj-1)*nz + kk
      u0(node)=u0(i1+kk)
  814 continue
c
  813 continue
      write(6,231)
  231 format(/10x,'initial condition for flow has been read in',/)
      do 232 i=1,nn
      u1(i)= u0(i)
  232 u2(i)= u0(i)
      goto 129
c
 1009 write(6,1111)
 1111 format(/10x,'!!! end-of-file detected in hinput.data file !!! ',
     +       /10x,'the file is not compatible with your grid ',
     +       /10x,'check source of old hinput.data file, and',
     +       /10x,'check present grid dimensions',
     +       /10x,'program stopping ...')
      stop
 1010 write(6,1112)
 1112 format(/10x,'!!! error detected in hinput.data file !!! ',
     +       /10x,'you are trying to read an initial head condition'
     +       /10x,'but this file does not exist; ',
     +       /10x,'rename a hxz.out file to hinput.data and try again',
     +       /10x,'or use init=0 and use a uniform initial head'
     +       /10x,'program stopping ...')
      stop
c
  129 continue
c
c ___________________________________________________________________
c -------------------------------------------------------------------
c
c
c     enter transport parameters
c
c     redefine coords in transport:
c     x,y,z coords become defined on large transport grid
c     in2 is new incidence array for transport grid:
c     ----------------------------------------------
c
      kcall=1
c     nzt=nz
      nezt=nez
      net=ne
      nnt=nx*ny*nzt
c
      call prism (x,y,z,in2,xl,yl,zl,nx,ny,nzt,nex,ney,nezt,nnt,net,
     +       maxnn,maxne,xlim,ylim,zlim,nlx,nly,nlz,ngx,ngy,ngz,
     +       mxgx,mxgy,mxgz,exl,eyl,ezl,map,mpa,inb,inbl,inbr,
     +       kcall,nf,lunsat,kgo)
c
      nyzt=ny*nzt
      write(6,109) nzt,nezt,nnt,net
 109  format(//10x,'transport grid ...',/10x,18('-'),
     + /10x,'number of nodes in vertical (z) direction  (nzt) = ',i9,
     + /10x,'number of z-elements in transport grid:   (nezt) = ',i9,
     + /10x,'total number of nodes in transport grid:   (nnt) = ',i9,
     + /10x,'total number of elements in transport grid:(net) = ',i9)
c
c     set up boundary arrays and condensation codes for heat transport 
c     ----------------------------------------------------------------
      leak2=.false. 
      ffc2=0.
      if(kmass.eq.0.or.kmass.eq.3.or.kmass.eq.4) then
      ktype = 1
c
      call mindex(maxnn,maxne,maxn,maxna,laa,np1,nnt,net,nf,n2,
     + nx,ny,nzt,in2,ic2,fb2,fc2,lc2,fbft,ib2,a,iaa2,ind2,
     + nw,kprt,ktype,nc2,m2,nb2,x,y,z,ara,datum,leak2,por,
     + wh,ampcos,ampsin,wpc,wps,slope,phasewt,map,kb2,lwtfgt0,ffc2)
      endif
c
c     set up boundary arrays and condensation codes for mass transport: c1
c     ---------------------------------------------------------------------
      ktype = 2
c
      if(kmass.eq.1.or.kmass.eq.2.or.kmass.eq.3.or.kmass.eq.4) then
      call mindex(maxnn,maxne,maxn,maxna,laa,np1,nnt,net,nf,n3,
     + nx,ny,nzt,in2,ic3,fb3,fc3,lc3,fbfc,ib3,a,iaa3,ind3,
     + nw,kprt,ktype,nc3,m3,nb3,x,y,z,ara,datum,leak3,por,
     + wh,ampcos,ampsin,wpc,wps,slope,phasewt,map,kb3,lwtfgt0,ffc3)
c     print*,'done mindex c1'
      endif
c
c     set up boundary arrays and condensation codes for mass transport: C2
c     ---------------------------------------------------------------------
      ktype = 3
c
      if(kmass.eq.4) then
      call mindex(maxnn,maxne,maxn,maxna,laa,np1,nnt,net,nf,n4,
     + nx,ny,nzt,in2,ic4,fb4,fc4,lc4,fbfc,ib4,a,iaa4,ind4,
     + nw,kprt,ktype,nc4,m4,nb4,x,y,z,ara,datum,leak3,por,
     + wh,ampcos,ampsin,wpc,wps,slope,phasewt,map,kb4,lwtfgt0,ffc4)
c     print*,'done mindex c2'
      endif
c
c     determine elemental volumes:
c     ----------------------------
      call volume(x,y,z,maxnn,ag,hag,in2,vt,maxne,net,vol1,kint,
     +    exl,eyl,ezl)
c
c     read properties of leaky boundary:
c     -------------------------------------------------
      read(5,*,err=8484) tcs,porsurf,satsurf
      if(leak2) write(6,4833) 
 4833 format(/10x,'top surface thermally leaky')
      if(.not.leak2) write(6,4834) 
 4834 format(/10x,'top surface not thermally leaky')
      write(6,166) tcs,porsurf,satsurf
 166  format(/10x,'top surface thermal parameters',
     +   /10x,'thermal conductivity of surface zone:',e12.3,
     +   /10x,'porosity of surface zone:            ',e12.3,
     +   /10x,'saturation of unsaturated zone (deactivated):  ',e12.3)
c
c     read top surface air temperature parameters
c     --------------------------------------------
      read(5,*,err=8498) surfmin,amp,period,phase,cutoff,gradt,gradx
      write(6,177) surfmin,amp,period,phase,cutoff,gradt,gradx
 177  format(/10x,'top surface air temperature variation',
     +      /10x,'minimum air temperature (C)  ... ',f12.4,
     +      /10x,'cosine curve amplitude  (C) ...... ',f12.4,
     +      /10x,'cosine curve period (days)  ...... ',f12.4,
     +      /10x,'phase shift (days)     ........... ',f12.4,
     +      /10x,'temperature cutoff  .............. ',f12.4,
     +      /10x,'temperature gradient over time (degC/day) .. ',f12.5,
     +      /10x,'temperature gradient over x dist. (degC/m) .. ',f12.5)
      if(period.lt.1e-10) then
        write(6,9848) 
 9848   format(/10x,'error ... cosine period=0',
     +         /10x,'program stopping')
      stop
      endif
c
      if(surfmin.lt.-50..or.surfmin.gt.200.
     +                  .or.amp.gt.50..or.amp.lt.-0.1
     +                  .or.phase.gt.365..or.phase.lt.-365.
     +                  .or.gradt.gt.1..or.gradt.lt.-1
     +                  .or.gradx.gt.1..or.gradx.lt.-1) then
        write(6,9838) 
 9838   format(/10x,'error ... ',
     +         /10x,'non-physical sine/cosine temperature parameters',
     +         /10x,'program stopping')
      stop
      endif            
c
c     over-ride with air temperatures from external file tair.data
c     -------------------------------------------------------------
      if(ktair.gt.0) then
          if(ierr8.gt.0) then
          write(6,2041)
 2041   format(/10x,'error... you are trying to read tair.data',
     +         /10x,'(ktair>0), but this file does not exist',
     +         /10x,'program stopping ...')
          stop
          endif              
       read(45,10) title3
       read(45,*) nair
       if(nair.gt.maxair) then
        write(6,5533) nair,maxair
5533    format(/10x,'error, #Tair points (',i6,') > max (',i5,')',
     +         /10x,'program stopping ...')
        call flush(6)
        stop
        endif
       read(45,*,err=3035) (timeair(i),tempair(i),i=1,nair)
       write(6,3844) nair,(i,timeair(i),tempair(i),i=1,5),
     +                    (i,timeair(i),tempair(i),i=nair-5,nair)
 3844  format(/10x,'surface air temperatures read from tair.data',
     +        /10x,'# data points ... ',i8,
     +        /10x,'first 5 and last 5 time/temperatures ... ',/,
     +        (10x,i6,2e20.7))
      call flush(6)
      endif
c
c     read top surface conductance layer thickness variation parameters
c     2016: tgradx2 not yet active
c     ------------------------------------------------------------------
c     k=0
c2302 continue
c     read(5,*)    surfmin2,amp2,period2,phase2,cutoff2,more
      read(5,*,err=8498) 
     +     surfmin2,amp2,period2,phase2,cutoff2,gradt2,gradx2
      write(6,1177) surfmin2,amp2,period2,phase2,cutoff2,gradt2,gradx2
1177  format(/10x,'conductance layer thickness variation',
     +       /10x,'minimum thickness (m)   ........ ',f12.4,
     +       /10x,'cosine curve amplitude  (C) ...... ',f12.4,
     +       /10x,'cosine curve period (days)  ...... ',f12.4,
     +       /10x,'phase shift (days)     ......... ',f12.4,
     +       /10x,'thickness cutoff  .............. ',f12.4,
     +       /10x,'thickness gradient over time (m/day) .... ',f12.4,
     +       /10x,'thickness gradient over x dist. (degC/m) .. ',f12.4)
c
c     read top surface conductance layer saturation variation parameters
c     2016: tgradx3 not yet active
c     ------------------------------------------------------------------
c     k=0
c2303 continue
c     read(5,*)   ix1,ix2,surfmin3,amp3,period3,phase3,cutoff3,more
      read(5,*,err=8498)
     +              surfmin3,amp3,period3,phase3,cutoff3,gradt3,gradx3
      write(6,2177) surfmin3,amp3,period3,phase3,cutoff3,gradt3,gradx3
2177  format(/10x,'conductance layer saturation variation',
     +       /10x,'minimum saturation    .......... ',f12.4,
     +       /10x,'cosine curve amplitude  (C) ...... ',f12.4,
     +       /10x,'cosine curve period (days)  ...... ',f12.4,
     +       /10x,'phase shift (days)     ......... ',f12.4,
     +       /10x,'saturation cutoff  ............. ',f12.4,
     +       /10x,'thickness gradient over time (m/day) .... ',f12.4,
     +       /10x,'thickness gradient over x dist. (degC/m) .. ',f12.4)

c
c     read heads to apply at left boundary
c     --------------------------------------
      if(ierr4.gt.0) then
        write(6,2039)
 2039   format(/10x,'no hbdy.data file found ... continuing')
        else

      read(48,10) title3
      read(48,*) numhbdy
      if(numhbdy.gt.maxbdy) then
        write(6,5534) numhbdy,maxbdy
5534    format(/10x,'error, #hbdy points (',i8,') > max (',i8,')',
     +         /10x,'program stopping ...')
        stop
        else
      read(48,*,err=3036) (thbdy(i),hbdy(i),i=1,numhbdy)
      write(6,3845) numhbdy,(thbdy(i),hbdy(i),i=numhbdy-5,numhbdy)
 3845 format(/10x,'bdy heads read from hbdy.data',
     +        /10x,'# data points ... ',i8,
     +        /10x,'last 5 time/heads ... ',/,
     +        (10x,2e20.7))      
      endif
      endif
c      
c     read transient bz and heat flux(W/m2) if the file exists
c     ---------------------------------------------------------
      do i=1,nf
       bz(i)=0.
       bzf(i)=0.
       bztdiff(i)=0.       
      enddo
      do i=1,maxbdy
      do j=1,maxbzz
       bzflx(i,j)=0.
       bzq(i,j)=0. 
      enddo
      enddo
      if(ierr6.gt.0) then
        lbzin=.false.
        write(6,2049)
 2049   format(/10x,'no bzin.data file found ... ok, continuing')
        else
        lbzin=.true.
      read(57,10) title3
      read(57,*) numbzbdy
      if(numbzbdy.gt.maxbdy) then
        write(6,5535) numbzbdy,maxbdy
5535    format(/10x,'error, #bzbdy points (',i8,') > max (',i8,')',
     +         /10x,'program stopping ...')
        stop
        else
      read(57,*,err=3040) (tbzbdy(i),nbzz(i),
     +(ix1bz(i,j),ix2bz(i,j),bzbdy(i,j),bzflx(i,j),bzq(i,j),j=1,nbzz(i)
     +                                                               ),
     +                                        bztdiff(i),i=1,numbzbdy)
      write(6,3846) numbzbdy
 3846 format(/10x,'top bdy bz data read from bzin.data',
     +/10x,'# data points...',i8,
     +/10x,'last 5:  time,tdiff,nzones, (ix1-ix2,bz,bzflx,bzq)_1,2... ')
      do  i=numbzbdy-5,numbzbdy
      if(i.gt.0) then
      write(6,3847) tbzbdy(i),bztdiff(i),nbzz(i),
     + (ix1bz(i,j),ix2bz(i,j),bzbdy(i,j),bzflx(i,j),bzq(i,j),j=1,maxbzz)
      endif
      enddo
 3847 format(10x,2e15.4,2x,i7,/,(15x,'(',2i5,3e12.4,')'))
      endif

      endif      
      call flush(6)
c
c ==================================================================
c
c      uniform velocity field ...
c      ivel=0 ... to use {v} from flow solution, enter
c      ivel=1 ... to read uniform {v} ... vxl,vyl,vzl
c      ivel=2 ... to use uniform v only in fractures, and v=0 elsewhere
c      ivel=3 ... to read vx,vz from veloc2d_in.data (3 titles, x,z,vx,vz) from veloc2dxz.out; vy=0:  
c      ----------------------------------------------------------------------------------------------
c
      read (5,*) ivel,vxl,vyl,vzl,iflipv
c  12 format(i5,3f10.0)

      if(iflipv.lt.0) write(6,9022)
 9022 format(/10x,'iflipv=-1: reversing velocity field ...')
      if(iflipv.eq.0) then
          write(6,8585)
 8585     format(/10x,'error, iflipv=0, stopping ...')
          stop
      endif

      if(ivel.eq.0.and.kcntrl.eq.2) then
              write(6,9021)
 9021         format(/10x,'Error: kcntrl=2 and ivel=0',
     +               /10x,'program stopping ...')
              stop
      endif
      
      if(ivel.eq.1) then 
      write (6,38) vxl,vyl,vzl
   38 format (//10x,'fixed groundwater velocity (pre-flipped): ',
     +         /10x,' x-direction ',e12.5,
     +         /10x,' y-direction ',e12.5,
     +         /10x,' z-direction ',e12.5/)
      do l=1,net
      vx(l)=vxl*float(iflipv)
      vy(l)=vyl*float(iflipv)
      vz(l)=vzl*float(iflipv)
      enddo
      endif

      if(ivel.eq.2) then
      write(6,3362)
 3362 format(/10x,'assuming uniform fracture velocity & v=0 in matrix')
      write (6,38) vxl,vyl,vzl
      do l=1,net
      vx(l)=1.0e-20
      vy(l)=1.0e-20
      vz(l)=1.0e-20
      enddo
      endif

      if(ivel.eq.3) then
        if(ierr5.gt.0) then
        write(6,2040)
 2040   format(/10x,'error... you are trying to read veloc2d_in.data',
     +         /10x,'(ivel=3), but this file does not exist',
     +         /10x,'program stopping ...')
        stop
        endif
      write(6,39) 
   39 format(/10x,'ivel=3: reading vx,vz from veloc2d_in.data; vy=0. ')
      read(53,10,err=3037) title1,title2,title3
      do i=1,nex
      do j=1,ney
      do k=1,nezt
      iel = (k-1)*nex*ney + (j-1)*nex + i 
      read(53,*,err=3037) xdum,zdum,vx(iel),vz(iel)
      vy(iel)=0.d0
      enddo
      enddo
      enddo
      write(6,41) vx(net),vz(net)
   41 format(10x,'veloc2d_in read successful: vx(net),vz(net): ',2e15.5)
      endif
c
c 161 continue
c
c     longitudinal and transverse dispersivity, diffusion coefficient
c     also read retardation and decay for mass transport option ...
c     read agefx age growth term (=1 day/day or = 0 day/day), use with kmass=2
c     kdisp = 0: original Burnett-Frind dispersion tensor
c     kdisp = 1: Lichtner dispersion tensor
c     alh = longitudinal dispersivity for horizontal flow
c     alv = longitudinal dispersivity for vertical flow (Lichtner form)
c     if kdisp = 0, alv is read but ignored 
c     ------------------------------------------------------------------------
c
      read (5,*,err=1140) alh,alv,ath,atv,dd,
     +                    rtrans,rtransc2,decaylin,decaylin2,agefx,kdisp
c
      al=alh
c     al = conventional longitudinal dispersivity, B-F form  
c
      if(kdisp.eq.0) write (6,722)
722   format (/10x,'Original Burnett-Frind dispersion tensor')
      if(kdisp.eq.1) write (6,723)
723   format (/10x,'Lichtner dispersion tensor')
c      
      write (6,40) al,alv,ath,atv,dd,
     +                          rtrans,rtransc2,decaylin,decaylin2,agefx
   40 format (//10x,'longitudinal dispersivity (m):      ',f12.5/
     +   10x,'longitudinal transverse dispersivity (Lichtner): ',f12.5/
     +          10x,'transverse horizontal dispersivity: ',f12.5/
     +          10x,'transverse vertical   dispersivity: ',f12.5/
     +          10x,'molecular diffusion coef. (m^2/s):  ',e12.5/
     +          10x,'retardation for mass transport c1:     ',e12.5/
     +          10x,'retardation for mass transport c2:     ',e12.5/
     +          10x,'decay coef. for transport c1: (sec^-1) ',e12.5/
     +          10x,'decay coef. for transport c2: (sec^-1) ',e12.5/
     +          10x,'agefx age growth term: (-):         ',e12.5)
c      
c     detect errors
c     --------------
      if(kint.eq.1 .and. decaylin.gt.0.) then
      write(6,724)
 724  format(/10x,'decay only coded for direct integration option',
     +       /10x,'please contact authors for latest update',
     +       /10x,'in the meantime, use kint=1', 
     +       /10x,'program stopping')
      stop
      endif

      if(kmass.eq.0.and.(decaylin.gt.0..or.agefx.gt.0.))then
      write(6,725)
 725  format(/10x,'incompatible kmass-decay-agefx'
     +   /10x,'kmass=0 (heat) but either decay>0 or agefx>0 ',
     +   /10x,'   if agefx>0 (age transport), must set kmass=2',
     +   /10x,'program stopping')
      stop
      endif
      if(decaylin.gt.0. .and.agefx.gt.0.)  then
      write(6,735)
 735  format(/10x,'incompatible decay-agefx'
     +   /10x,'agefx>0 then decay must be =0',
     +   /10x,'program stopping')
      stop
      endif
c
c     apply decaylin to all elements
c     ------------------------------
      do i=1,net
      decay(i) = decaylin
      decayc2(i) = decaylin2
      enddo
c
c     thermal parameters:
c     set modelwu= 0 for no effect
c     ts = freezing front (start of freezing)
c     cm,pm now read by element ....
c     Jan 2016: Wu, kr separated in matrix / fractures
c     Aug 2021 read phi (bulk density conversion factor)
c              read alpha0 and alphap for relative k calculation from retained particles
c     Dec 2021: p,q removed here, read by element above
c     ----------------------------------------------------------------------------------
      read(5,*,err=3038) cf,pf,ci,pi,tclw,tcli,wlh,phi,alpha0,alphap
      read(5,*,err=3139) ts,omega,modelwu,modelkr
      read(5,*,err=3039) p2,q2,ts2,omega2,modelwu2,modelkr2                  !fractures
      write(6,50) cf,pf,ci,pi,tclw,tcli,wlh,phi,alpha0,alphap
   50 format(/10x,'specified thermal parameters:',
     +       /10x,30('-'),/,
     +        10x,'fluid specific heat (J/kg/C):     ',e10.3,/
     +        10x,'fluid density (kg/m3):            ',e10.3,/
     +        10x,'ice specific heat:                ',e10.3,/
     +        10x,'ice density                       ',e10.3,/
     +        10x,'water thermal conductivity (J/m/s/C): ',e10.3,/
     +        10x,'ice thermal conductivity:             ',e10.3,/
     +        10x,'latent heat of water (J/kg):          ',e10.3,/
     +        10x,'fraction phi for bulk density:        ',e10.3,/
     + 10x,'alpha0 specific surface area for clean bed (1/m):  ',e10.3,/
     + 10x,'alphap specific surf. area of particles (1/m)   :  ',e10.3)

      if(tclw.gt.100. .or. tcli.gt.100) then
              write(6,9033) tclw,tcli
9033          format(/10x,'error in thermal K of water or ice',
     +               /10x,'tclw,tcli ... ',2e15.5,
     +               /10x,'program stopping ...')
              stop
      endif     

       if(gamma.lt.1.e-10.and.phi.gt.998.) write(6,8829)
 8829  format(/10x,'gamma=0 and phi>999; linear mass transport run')

      write(6,511) ts,omega,modelwu,modelkr
  511 format(/10x,'parameters for functions Wu and kr in porous media:',
     +  /10x,'(unfrozen moisture curve parameters p,q read by element)',
     +       /10x,51('-'),/,
     +        10x,'freezing front temperature (Ts):         ',e10.3,/
     +        10x,'relative permeability exponent omega:    ',e10.3,/
     +        10x,'unfrozen moisture content model# wu:      ',i5,/
     +        10x,'relative permeability function model# kr: ',i5)

      write(6,522) p2,q2,ts2,omega2,modelwu2,modelkr2
  522 format(/10x,'parameters for functions Wu and kr in fractures:',
     +       /10x,48('-'),/,
     +        10x,'terminal unfrozen moisture content(p):   ',e10.3,/
     +        10x,'unfrozen moisture content curvature(q):  ',e10.3,/
     +        10x,'freezing front temperature (Ts):         ',e10.3,/
     +        10x,'relative permeability exponent omega:    ',e10.3,/
     +        10x,'unfrozen moisture content model# wu:      ',i5,/
     +        10x,'relative permeability function model# kr: ',i5)
      call flush(6)

      if(modelwu.lt.0 .or. modelwu.gt.4) then
           write(6,2288)
 2288      format(/10x,'Error in Wu: modelwu is not 0-4',
     +            /10x,'program stopping ...')
           stop
      endif
      if(modelkr.lt.0 .or. modelkr.gt.2) then
           write(6,2287)
 2287      format(/10x,'Error in kr: modelkr is not 0-4',
     +            /10x,'program stopping ...')
           stop
      endif

      do iel=1,net
      if( (modelwu.eq.2.or.modelwu.eq.3) .and. qq(iel) .lt.1e-20 )then
      write(6,2289) 
      stop
      endif 
      enddo

      if( (modelwu.eq.2.or.modelwu.eq.3) .and. q2.lt.1e-20 )then
           write(6,2289) 
 2289      format(/10x,'Error in Wu: q cannot=0 with Wu model 2 or 3',
     +            /10x,'Program stopping ...')
           stop
       endif
      if(modelwu.eq.0)then
           write(6,2291) 
 2291      format(/10x,'ModelWu = 0, so unfrozen moisture always = 1')
       endif       
      if(modelkr.eq.0)then
           write(6,2292) 
 2292      format(10x,'Modelkr = 0, so relative perm. kr always = 1')
       endif       
c
c     initial conditions: transport ...
c     ==================================
c added May22
      do i=1,nnt
      spn0(i)=0.0
      spn1(i)=0.0
      enddo
      do i=1,net
      spe0(i)=0.0
      spe1(i)=0.0
      enddo
c      
c     default initial condition is dirichlet boundary array:
c     ------------------------------------------------------
      do i=1,nnt
      if(fc2(i).gt.-900) then
      t0(i) = fc2(i)
      t1(i) = fc2(i)
      t2(i) = fc2(i)
      endif
      if(fc3(i).gt.-900) then
      c0(i) = fc3(i)
      c1(i) = fc3(i)
      c2(i) = fc3(i)
      endif
      if(fc4(i).gt.-900) then
      d0(i) = fc4(i)
      d1(i) = fc4(i)
      d2(i) = fc4(i)
      endif
      enddo
c
c     best to read T i.c. anyway - some transport params will depend on T
c     so need to set a uniform background T
c     if(kmass.eq.0.or.kmass.eq.3) 
c     ------------------------------------------
c     if(kmass.eq.0.or.kmass.eq.3.or.kmass.eq.4) 
      write (6,29)
   29 format (/10x,'initial temperature conditions'/,
     +         10x,'   x1   x2   y1   y2   z1   z2',4x,
     +         'temperature range')
c
c     over-ride default with specified initial conditions :
c     enter x-nodes i1 to i2, and bottom,top surface temperatures
c     temperatures are linearly interpolated, 
c     fc2 (Dirichlet temperature) is NOT overwritten unless fc2 = -999.
c     --------------------------------------------------------------
c
      tsurfmin=999.
      tsurfmax=-999.
      tdays=0.
   45 read (5,*,err=3032) i1,i2,j1,j2,k1,k2,uinb,uint,more
      write (6,34) i1,i2,j1,j2,k1,k2,uinb,uint,more
   34 format (10x,6i5,2e12.4,i5)
      if( i1.lt.1 .or. i2.gt.nx .or.
     +    j1.lt.1 .or. j2.gt.ny .or. k1.lt.1 .or. k2.gt.nzt
     +    .or.k2.lt.k1) then
      write(6,9385) nx,ny,nzt
 9385 format(/10x,'error in initial conditions for T transport ...',
     +       /10x,'nx,ny,nzt ...',3i8,
     +       /10x,'program stopping ...')
      stop
      endif
      call flush(6)
      do i=i1,i2
      do j=j1,j2
      nbot=(i-1)*nyzt + (j-1)*nzt + k1
      ntop=nbot+k2-k1
      do k=k1,k2
      node= nbot + k - k1
c      print*,'node= ',node,' ntop= ',ntop,' nbot= ',nbot
c      print*,'z(ntop)= ',z(ntop),'  z(nbot)= ',z(nbot)
c     tk=uinb
c     if( .not. leak2)
c    + tk = (z(node)-z(nbot))/(z(ntop)-z(nbot))*(uint-uinb) + uinb
c      ----------------
c     if(leak2 .and. k2.eq.nzt) then     !if top is leaky, use cosine curve for top temperature
c     uint =
c    +    surfat(surfmin,amp,period,phase,tdays,cutoff,gradt)
c    +                                         +x(ntop)*gradx
c     endif
      if(ntop.ne.nbot) then
            tk = (z(node)-z(nbot))/(z(ntop)-z(nbot))*(uint-uinb) + uinb
      else
            tk = uint
      endif
c      ------------------
c bug with abs() ? removed sept 2021
c     if(ic2(node).ne.1 .or. (abs(fc2(node)).lt.-900.)) then   !only for free nodes, or if Dirichlet T = -999
      if(ic2(node).ne.1 .or. (fc2(node)).lt.-900.) then   !only for free nodes, or if Dirichlet T = -999
c      print*,'ICs ... k,node= ',k,node,' tk= ',tk
        fc2(node)= tk
        t0(node) = tk
        t1(node) = tk
        t2(node) = tk
      endif
      enddo
      enddo
      enddo

      if(more.gt.0) go to 45

      do i=1,nnt
      tsurfmin = min(tsurfmin,t2(i))
      tsurfmax = max(tsurfmax,t2(i))           
      enddo
c
      write(6,3374) tsurfmin,tsurfmax
 3374 format(/10x,'initial min,max temps : ',2e15.5)
c
c     ===============================================================
c
c     over-ride default with specified initial conditions :  concentrations
c     enter x-nodes i1 to i2, and bottom,top surface concentrations
c     concentrations are linearly interpolated, 
c     fc3 (Dirichlet concentrations) is NOT overwritten unless fc2 = -999.
c     only read for mass transport
c     --------------------------------------------------------------
c
      if(kmass.eq.1.or.kmass.eq.2.or.kmass.eq.3.or.kmass.eq.4) then
              if(kmass.ne.2) write (6,209)
  209 format (/10x,'initial concentration conditions'/,
     +         10x,'   x1   x2   y1   y2   z1   z2',4x,
     +         'comp #   concentration range')
              if(kmass.eq.2) write (6,2099)
 2099  format (/10x,'initial age conditions'/,
     +         10x,'   x1   x2   y1   y2   z1   z2',4x,
     +         'comp #   age range')
c
      csurfmin=999.
      csurfmax=-999.
      csurfmin2=999.
      csurfmax2=-999.
      tdays=0.
 4045 read (5,*,err=3142) i1,i2,j1,j2,k1,k2,icomp,uinb,uint,more
      write (6,3004) i1,i2,j1,j2,k1,k2,icomp,uinb,uint,more
 3004 format (10x,6i5,i10,2e12.4,i5)
      if( i1.lt.1 .or. i2.gt.nx .or.
     +    j1.lt.1 .or. j2.gt.ny .or. k1.lt.1 .or. k2.gt.nzt) then
      write(6,9395) i1,i2,j1,j2,k1,k2
 9395 format(/10x,'error in initial c condition nodes for transport...',
     +       /10x,'i1,i2,j1,j2,k1,k2: ',6i5,
     +       /10x,'program stopping ...')
      stop
      endif
      call flush(6)
      do i=i1,i2
      do j=j1,j2
      nbot=(i-1)*nyzt + (j-1)*nzt + k1
      ntop=nbot+k2-k1
      do k=k1,k2
      node= nbot + k - k1
c      print*,'node= ',node,' ntop= ',ntop,' nbot= ',nbot
c      print*,'z(ntop)= ',z(ntop),'  z(nbot)= ',z(nbot)
      tk=uinb
      if(k2.gt.k1)
     + tk = (z(node)-z(nbot))/(z(ntop)-z(nbot))*(uint-uinb) + uinb
      csurfmin = min(csurfmin,tk)
      csurfmax = max(csurfmax,tk)           
c bug with abs() ? removed sept 2021
c     if(ic3(node).ne.1 .or. (abs(fc3(node)).lt.-900.)) then   !only for free nodes, or if Dirichlet T = -999
      if(icomp.eq.1.and.(ic3(node).ne.1 .or. (fc3(node).lt.-900.))) then   !only for free nodes, or if Dirichlet T = -999
        fc3(node)= tk
        c0(node) = tk
        c1(node) = tk
        c2(node) = tk
      endif
      if(icomp.eq.2.and.(ic4(node).ne.1 .or. (fc4(node).lt.-900.))) then   !only for free nodes, or if Dirichlet T = -999
        fc4(node)= tk
        d0(node) = tk
        d1(node) = tk
        d2(node) = tk
      endif
      enddo
      enddo
      enddo

      if(more.gt.0) go to 4045
c
      do i=1,nnt
      csurfmin = min(csurfmin,c2(i))
      csurfmax = max(csurfmax,c2(i))           
      csurfmin2 = min(csurfmin2,d2(i))
      csurfmax2 = max(csurfmax2,d2(i))           
      enddo
      write(6,3384) csurfmin,csurfmax
 3384 format(/10x,'initial min/max concs; component 1: ',2e15.5)
      if(kmass.eq.4) write(6,3385) csurfmin2,csurfmax2
 3385 format(/10x,'initial min/max concs; component 2: ',2e15.5)

      endif
c
c    if ffc2=-999, assign cosine temps for top Dirichlet T nodes
c    will override therm_in temps if kgo=1      
c    -----------------------------------------------------------
      if(ffc2.lt.-998.) then
      write(6,4747) 
 4747 format(/10x,'assigning cosine air temps as type-1 top nodes')   
      do ntop=nzt,nnt,nzt
      ttt=  surfat(surfmin,amp,period,phase,tdays,cutoff,gradt)
     +                                              +x(ntop)*gradx
c      write(6,*) x(ntop),surfmin,amp,ttt
        fc2(ntop)= ttt
        t0(ntop) = ttt
        t1(ntop) = ttt
        t2(ntop) = ttt
       enddo
c      write(6,4748) tdays,x(nzt),gradx,t0(nzt),t0(nnt) 
c 4748 format(10x,'Type-1 cosine temps: ',
c     + /10x,'tdays,x(nzt),gradx,t0(nzt),t0(nnt) : ',5e15.5)   
      endif
c
      tsurfmin=999.
      tsurfmax=-999.
      do i=nzt,nn,nzt
      tsurfmin = min(tsurfmin,t0(i))
      tsurfmax = max(tsurfmax,t0(i))    
      enddo
      write(6,3174) tsurfmin,tsurfmax
 3174 format(/10x,'min,max temps along surface : ',2e15.5)      
c
c    detect negative temperatures for non-linear iterations
c    -------------------------------------------------------
      lneg=.false.
      if(surfmin.lt.0.) lneg = .true.
      do i=1,nnt
      if(fc2(i).lt.0. .or. t2(i).lt.0.) lneg=.true.
      enddo
c
c    calculate initial temperature-dependent thermal conductivities
c    for transport grid
c    first assign porsurf and saturation to unsat zone
c    -------------------------------------------------------
c     if(ksat.eq.0) then
c       do l=ne+1,net
c       por(l)=porsurf
c       sw(l)=satsurf    !hardwire - use satsurf
c       tclm(l)=tclm(ne) !hardwire - assume same thermal K solids in unsat zone as at watertable
c       enddo
c     endif

      tclmin = 999.
      tclmax =-999.
      do l=1,net
      do 1101 i=1,8
 1101 inl(i)=in2(l,i)
      tav1=0.
      tav2=0.
      do i=1,8
        tav1=tav1+t0(inl(i))
        tav2=tav2+t1(inl(i))
      enddo
      tav1=tav1/8.0
      tav2=tav2/8.0
      tavg=(tav1+tav2)/2.0
      ppl=pp(l)
      qql=qq(l)
      ww=por(l)*sw(l)*wu(tavg,ppl,qql,modelwu,ts)
      wi=por(l)*sw(l)*(1.-wu(tavg,ppl,qql,modelwu,ts))
      wm=(1.-por(l))

      tkl(l) = (ww*sqrt(tclw)+wi*sqrt(tcli)+wm*sqrt(tclm(l)))**2        !more physically realistic
      if(l.eq.1) write(6,7746)
7746  format(/10x,'using square-root thermal K weighting')      

c     tkl(l) = ww*tclw + wi*tcli + wm*tclm(l)                       !hardwire Interfrost
c      if(l.eq.1) write(6,7746)
c7746  format(/10x,'using geometric (linear) thermal K weighting')      

      if(tkl(l).gt.tclmax) tclmax = tkl(l)
      if(tkl(l).lt.tclmin) tclmin = tkl(l)
      enddo

      write(6,9384) tclmin,tclmax
 9384 format(/10x,'Initial limits for bulk thermal K of porous medium:',
     +       /10x,'minimum thermal K: ',e12.3,
     +       /10x,'maximum thermal K: ',e12.3)
c
      if(tclmin.le.0.) then
       write(6,9855)
 9855  format(/10x,'error ... bulk thermal K <= 0',
     + /10x,'program stopping ...')
      stop
      endif      
c
c     sample derived parameters ... base values with pf=1000 kg/m**3
c     sat = saturation of 1st element
c     ----------------------------------------------------------
      tdays = 0.
      sat =
     +surfat(surfmin3,amp3,period3,phase3,tdays,cutoff3,gradt3)
c   
      xpor = por(1)
      cpf  = cf * pf
c     cpm  = cm * pm
      cpm = cpsm(1)       !reference element 1
      cpi  = ci * pi
      ppl = pp(1)         !reference element 1
      qql = qq(1) 
      ww   = xpor*sw(1)*wu(t0(1),ppl,qql,modelwu,ts)
      wi   = xpor*sw(1)*(1.-wu(t0(1),ppl,qql,modelwu,ts))
      cpe  = ww*cpf + wi*cpi + (1.-xpor)*cpm
      cpeu = porsurf*sat*wu(t0(nnt),ppl,qql,modelwu,ts)*cpf 
     +     + porsurf*sat*(1.-wu(t0(nnt),ppl,qql,modelwu,ts))*cpi 
     +     + (1.-porsurf)*cpm

      td   = tkl(1)/cpe
c     pd   = (1.-xpor)*pm    !dry density not needed

      clh = xpor * pi * wlh
      ret = -999.
      retu = -999.
      if(ww.gt.0.) then
          ret  = (cpe+clh*dwu(tavg,ppl,qql,modelwu,ts))/(ww*cpf)
          retu = cpeu/(ww*cpf)
      endif
      write(6,51) cpe,cpeu,td,ret,retu,clh,ww,wi
   51 format(/10x,'derived thermal parameters ',
     +            '(at tdays=0, ref temp for pf, element 1)',
     +       /10x,'(-999 = undefined)',/
     +     10x,65('-'),//,
     +     10x,'sat.   aquifer heat capacity (J/m3/C):  ',e10.3,/
     +     10x,'unsat. aquifer heat capacity:           ',e10.3,/
     +     10x,'aquifer thermal diffusivity (m2/s):     ',e10.3,/
     +     10x,'sat.   retardation coefficient:         ',e10.3,/
     +     10x,'unsat. retardation coefficient (top):   ',e10.3,/
     +     10x,'latent heat coefficient (L*p*ww):       ',e10.3,/
     +     10x,'ww: unfrozen water content:             ',e10.3,/,
     +     10x,'wi:   frozen water content:             ',e10.3,/)
c
c    ======================================================================
c
c     adjust all given fixed heads to 
c     equivalent freshwater heads for the density case:
c     H = Ho + gamma * c (Ztop - z)   see the User Guide.
c     (where Ho is original head)
c     -----------------------------------------------------------------
       do 772 i=1,nx
       do 772 j=1,ny
       nodetop = (i-1)*nyz + j*nz
       do 772 k=1,nz
       node=(i-1)*nyz + (j-1)*nz + k
       if(ic(node).eq.1) then
        depth = z(nodetop)-z(node)
        beta = gamma*c2(node)*depth           !only c2 (component 1) contributes to density, not d2
        fc(node) = fc(node) + beta
        u0(node) = u0(node) + beta
        u1(node) = u0(node)
        u2(node) = u0(node)
       endif
  772  continue
c
c    ======================================================================
c      
c     read restart file to start at specified time 
c     --------------------------------------------
      if(kgo.gt.0 .and. kgo.lt.4) then
      if(ierr1.gt.0) then
         write(6,1115)
 1115    format(/10x,'!!! error detected in therm_in.data file !!! ',
     +          /10x,'you are trying to restart a previous run but'
     +          /10x,'this file does not exist',
     +          /10x,'use kgo=0 or rename a t3d file to therm_in.data',
     +          /10x,'program stopping ...')
         stop
      endif
      xlast=x(nnt)
      read(7,10,end=2001) title3,title4,title5
      read(7, *,end=2001) (x(i),y(i),z(i),tq(i),cq(i),pq(i),i=1,nnt)    !use tq,cq,pq as temporary for t,c,h
      if(abs(xlast-x(nnt)).gt.1.e-5) write(6,2019)
c     read(7,10,end=2001) title3
c     read(7, *,end=2001) (dum1,dum2,dum3,pq(i),i=1,nn)    !use pq as temporary for heads
      goto 2002
 2001 write(6,2004)
 2004 format(/10x,'end of file detected in therm_in.data',
     +       /10x,'this file is either empty or not compatible with'
     +       /10x,'your existing data file - please check grid size',
     +       /10x,'and source of your existing therm_in.data file',
     +       /10x,'program stopping ...')
      stop
 2019 format(/10x,'!! x(nn) in restart file therm_in.data is not ',
     +       /10x,'compatible with the existing grid (x coords) !!',
     +       /10x,'... continuing anyway .. ')
 2002 continue
      if(kgo.eq.2.or.kgo.eq.3) then                        !temperatures
      do 434 i=1,nnt
      t0(i)=tq(i)
      t1(i)=t0(i)
      t2(i)=t0(i)
  434 continue
      endif
      if(kgo.eq.3) then                                    !concentrations
      do 435 i=1,nnt
      c0(i)=cq(i)
      c1(i)=c0(i)
      c2(i)=c0(i)
  435 continue
      endif
      if(kgo.eq.1.or.kgo.eq.3) then                        !heads
      do 436 i=1,nn
      u0(i)=pq(i)
      u1(i)=u0(i)
      u2(i)=u0(i)
      if(ffc.lt.-998.) fc(i)=u0(i)     !hardwire set heads from therm_in.data as fc Dirichlet nodes if ffc=-999.
  436 continue
      endif

      write(6,3002) x(nnt),y(nnt),z(nnt)
 3002 format(/10x,'restart read successfull:',
     +       /10x,'last x,y,z: ',3f14.3)
      call flush(6)
      endif
c
c     ---------------------------------------------------------------
c
c     check head i.c. = 0 ?
c     ------------------------
      do i=1,nn
      if(u0(i).gt.0.) lzero=.false.
      enddo
      if(lzero) then
      write(6,8322)
 8322 format(/10x,'head initial condition should not be all = 0',
     +       /10x,'program stopping ...')
      stop
      endif

c      
c ___________________________________________________________________
c
c      convergence criteria, iteration limits, weighting factors:
c      ----------------------------------------------------------------
      read(5,*,err=8001) ccp,cct,ccw,ccc,maxit1,maxit2
  150 format(3f10.0,2i5)
      write(6,151) ccp,cct,ccw,ccc,maxit1,maxit2
  151 format(/10x,'convergence criteria for pressure        : ',f10.5,
     +       /10x,'convergence criteria for temperature     : ',f10.5,
     +       /10x,'convergence criteria for watertable      : ',f10.5,
     +       /10x,'convergence criteria for concentrations  : ',f10.5,
     +       /10x,'iteration limit for non-linearity        : ',i10,
     +       /10x,'iteration limit for watertable mound     : ',i10,/)
c
      read(5,*) wfu,wft,tsa
      write(6,173) wfu,wft,tsa
  173 format(/10x,'extrapolation factor for heads:',f10.5,
     +       /10x,'extrapolation factor for temp.:',f10.5,
     +       /10x,'temperature search area for tpeak:',f10.5)
c
c     -------------------------------------------------------------------
c
c     enter node in y-direction for x/z contouring plane:
c     ---------------------------------------------------
      read(5,*) (knox(i),i=1,5)
      read(5,*) (knoy(i),i=1,5)
      read(5,*) (knoz(i),i=1,5)
      write(6,152) (knox(i),i=1,5)
  152 format(/10x,'x-nodes for y/z contouring plane   : ',5i5)
      write(6,252) (knoy(i),i=1,5)
  252 format(10x,'y-nodes for x/z contouring plane   : ',5i5)
      write(6,352) (knoz(i),i=1,5)
  352 format(10x,'z-nodes for x/y contouring plane   : ',5i5,/)
c
      do i=1,5
      if(knox(i).lt.0 .or. knox(i).gt.nx) then
              write(6,3744) knox(i)
 3744         format(/10x,'Error with knox ... ',i5,' stopping ... ')
              stop
      endif
            if(knoy(i).lt.0 .or. knoy(i).gt.ny) then
              write(6,3745) knoy(i)
 3745         format(/10x,'Error with knoy ... ',i5,' stopping ... ')
              stop
      endif
            if(knoz(i).lt.0 .or. knoz(i).gt.nzt) then
              write(6,3746) knoz(i)
 3746         format(/10x,'Error with knoz ... ',i5,' stopping ... ')
              stop
      endif
      enddo
c
c     read specific print times
c     -------------------------
      read(5,*) (pt(i),i=1,5)
      write(6,339) (pt(i),i=1,5)
 339  format(/10x,'3D print times (days) ...',/10x,5e12.3)
 229  format(5f10.0)
c
c     ==================================================================
c
c     write(17,*) 'moment     ... tdays,xbpk,zbpk,xbar,zbar' 
      write(18,*) 'variables="time(d)","peak temp(C)","minimum temp(C)"'
      write(18,*) 'zone i=9999, T="replace with real # on last line"'
      write(19,1989) 
 1989 format('Title="htotr=wrt_Tmpa;htota=wrt_Tref;hinr=pq;',
     +         'hts2=cum_leak2bdy;hqtop=cum_top_recharge,Vwat,Vice"')
      write(19,1990) 
 1990 format('variables="time","htotr(J)","htota(J)","hinr(J)"',
     +                       ',"hts2-leak(J)","hqtop(J)","Vwat","Vice"')
      write(19,*) 'zone i=9999, T="htot(J); replace 9999 with bottom #"'
      write(55,*) 'Title="Multi-Heater Column 12cm avg, z vs Tavg"'
      write(55,*) 'variables="z","12cmTavg"'
      write(56,*) 'Title="Multi-Heater Column 12cm avg, Tavg vs. time"'
      write(56,*) 'variables="tdays","12cmTavg"'
c     write(23,*) 'therm_var  ... xvar,zvar,xzvar'
c
c     initial heat content of aquifer:
c     --------------------------------
      sumin=0.
c     tref= -273.15
c     tref= 0.0
      tref= t2(nzt/2)      !for TH2, use T(mid-znode on left) as reference temperature for heatfluxr
      htota=0.
      htotr=0.
      hssa=0.
      hssr=0.
      hoc=0.
      hic=0.
      hts2=0.
      tzero = -273.15
      vice=0.0
      vwat=0.0
      do 61 l=1,net
      tavg=0.
      tavgr=0.
      hqtop = 0.

      do j=1,8
       tavg=tavg+ t2(in2(l,j))
c      tavgr=tavgr+ (t2(in2(l,j))-t2(mpa(in2(l,j))))
       tavgr=tavgr+ (t2(in2(l,j))+273.)
      enddo
      tavg=tavg/8.
      tavga=tavg-tref
      tavgr=tavgr/8.

      cpf=cf*den(tavg,lheat,lmass,lage)
      ppl = pp(l)
      qql = qq(l)
      ww   = por(l)*sw(l)*wu(tavg,ppl,qql,modelwu,ts)
      wi   = por(l)*sw(l)*(1.-wu(tavg,ppl,qql,modelwu,ts))
      cpe  = ww*cpf + wi*cpi + (1.-por(l))*cpsm(l)

      vice=vice+vt(l)*wi
      vwat=vwat+vt(l)*ww
c      if(l.gt.ne) cpe=por(l)*sat*cpf + (1.d0-por(l))*cpsm(l)

      htota=htota + tavga * vt(l)*cpe
      htotr=htotr + tavgr * vt(l)*cpe
   61 continue

      write(6,64) tref,htota,htotr,vwat,vice
   64 format(/10x,'Reference temperature for heat: tref   : ',e12.4,
     +       /10x,'initial absolute heat content (wrt Tref): ',e12.4,
     +       /10x,'initial relative heat content (wrt t2mpa: ',e12.4,
     +       /10x,'initial volume of water:                  ',e12.4,
     +       /10x,'initial volume of ice:                    ',e12.4)
      write(19,644) htotr,htota,vwat,vice
 644  format(13x,'0.0',2e16.6,3(13x,'0.0'),2e16.6)
c
c     ==========================
c
c     print initial head and temperature & c distribution
c     -----------------------------------------------
      write(8,8071)
      write(9,8072)
      write(25,8073)
      write(62,8074)
 8071 format('variables="x","z","Temp.","Wu"') 
 8072 format('variables="x","z","head"')
 8073 format('variables="x","z","T-Tbck"')
 8074 format('variables="x","z","c1","c2","Spn"')    !Spn: solid particle (retained) nodal
      tdays=0.
      do 1260 jny=1,5
      if(knoy(jny).eq.0) goto 1260
      izonexz=izonexz+1
c      write(8,808) nzt,nx,tdays,knoy(jny)
      write(8,2246) izonexz,(tdays/365.),nzt,nx,(tdays/365.),knoy(jny)
2246  format('text x=75., y=70.,f=helv,cs=frame,hu=point,h=20,zn=',i5,
     +   ',c=black,bx=filled,bxf=white,bxo=white,t="',f11.3,' yrs "',
     +   /,'zone i=',i5,', j=',i5,',T="t=',f14.3,' yrs,knoy=',i3,'"')
      if(lmass) 
     + write(62,2246) izonexz,(tdays/365.),nzt,nx,(tdays/365.),knoy(jny)
      write(25,808) nzt,nx,tdays,knoy(jny)
      write(9,810) nz,nx,tdays,knoy(jny)
      do 1805 i=1,nx
      i1= (i-1)*nyzt + (knoy(jny)-1)*nzt
      if1= (i-1)*nyz + (knoy(jny)-1)*nz
      write(8,8091) (x(i1+k),z(i1+k),t2(i1+k),
     +              wu(t2(i1+k),ppn(i1+k),qqn(i1+k),modelwu,ts),k=1,nzt)
      if(lmass) 
     +     write(62,8091) (x(i1+k),z(i1+k),c2(i1+k),d2(i1+k),spn0(i1+k),
     +                                                          k=1,nzt)
      write(25,809) (x(i1+k),z(i1+k),(t2(i1+k)-t2(k)),k=1,nzt)
      write(9,809) (x(map(if1+k)),z(map(if1+k)),u2(if1+k),k=1,nz)
 1805 continue
 1260 continue
c
c     yz head & temps, section data
c     -------------------------------
      do jnx=1,5
      if(knox(jnx).eq.0) goto 960
      izoneyz=izoneyz+1
      if(jnx.eq.1) then                           !first zone
      write(16,8011) izoneyz,tdays,nzt,ny,tdays,knox(jnx)
8011  format('title="yz temp"',
     +  /,'variables="y","z","Temp","Wu"',
     +  /,'text x=75., y=70.,f=helv,cs=frame,hu=point,h=20,zn=',i5,
     +    ',c=black,bx=filled,bxf=white,bxo=white,t="',f11.3,' days "',
     +  /,'zone i= ',i5,', j= ',i5,
     +    ',T= "t= ',f8.1,' days; yz section ',i3,'"')

      write(21,8409) nz,ny,tdays,knox(jnx)
8409  format('title="yz head"',
     +     /,'variables="y","z","head"',
     +     /,'zone i= ',i5,', j= ',i5,
     +       ',T= "t= ',f12.3,' days; yz section ',i3,'"')

c     -----------------------------
      do i=1,ny
      i1= (knox(jnx)-1)*nyzt + (i-1)*nzt
      write(16,8091) ((y(i1+k)),z(i1+k),t2(i1+k),
     +             wu(t2(i1+k),ppn(i1+k),qqn(i1+k),modelwu,ts),k=1,nzt)
      enddo
      do i=1,ny
      i1= (knox(jnx)-1)*nyz + (i-1)*nz
      write(21,809) ((y(map(i1+k))),z(map(i1+k)),u2(i1+k),k=1,nz)
      enddo

      else                      ! next zones
c                          =====================

      write(16,8012) nzt,ny,tdays,knox(jnx)
8012  format('zone i= ',i5,', j= ',i5,',D=(1,2)',
     +       ',T= "t= ',f12.3,' days; yz section ',i3,'"')

      write(21,8410) nz,ny,tdays,knox(jnx)
8410  format('zone i= ',i5,', j= ',i5,',D=(1,2)',
     +       ',T= "t= ',f12.3,' days; yz section ',i3,'"')
c
c     ------------------------------------
      do i=1,ny
      i1= (knox(jnx)-1)*nyzt + (i-1)*nzt
      write(16,8091) (t2(i1+k),wu(t2(i1+k),
     +                         ppn(i1+k),qqn(i1+k),modelwu,ts),k=1,nzt)
      enddo
      do i=1,ny
      i1= (knox(jnx)-1)*nyz + (i-1)*nz
      write(21,809) (u2(i1+k),k=1,nz)
      enddo

      endif                  
      enddo                !next knoz
  960 continue
c
c     xy initial head t=0 plot: 
c     -------------------------
c     
      do k=1,5
      if(knoz(k).lt.1.or.knoz(k).gt.nz) goto 8334
      iz = knoz(k)

      if(k.eq.1) then
      write(46,8449) ny,nx,tdays,iz
8449  format('title="xy head"',
     +     /,'variables="x","y","head"',
     +     /,'zone i= ',i5,', j= ',i5,
     +       ',T= "t= ',f8.1,' days; xy section ',i3,'"')
      do  i=1,nx
      i1=(i-1)*nyz 
      write(46,809)(x(map(i1+kz)),y(map(i1+kz)),u2(i1+kz),
     +               kz=iz,nyz,nz)
      enddo
      else

      write(46,8349) ny,nx,tdays,iz
8349  format('zone i= ',i5,', j= ',i5,', D=(1,2)',
     +       ',T= "t= ',f8.1,' days; xy section ',i3,'"')
      do  i=1,nx
      i1=(i-1)*nyz 
      write(46,809)(u2(i1+kz),kz=iz,nyz,nz)
      enddo
      endif
 8334 continue
c
c     xy initial temp t=0 plot:
c     ------------------------
      if(knoz(k).lt.1) goto 7334
      izonexy=izonexy+1
      if(k.eq.1) then
      write(14,8419) izonexy,tdays,ny,nx,tdays,knoz(k)
8419  format('title="xy temp"',
     +     /,'variables="x","y","Temp","Wu","Conc"',
     +     /,'text x=75., y=70.,f=helv,cs=frame,hu=point,h=20,zn=',i5,
     +   ',c=black,bx=filled,bxf=white,bxo=white,t="',f12.4,' days "',
     +     /,'zone i= ',i5,', j= ',i5,
     +       ',T= "t= ',f12.4,' days; xy section ',i3,'"')
      do  i=1,nx
      i1=(i-1)*nyzt 
      write(14,8091) (x(i1+kz),y(i1+kz),t2(i1+kz),
     +               wu(t2(i1+kz),ppn(i1+kz),qqn(i1+kz),modelwu,ts),
     +                                        c2(i1+kz),kz=iz,nyzt,nzt)
      enddo
      else
      write(14,8319) ny,nx,tdays,iz
8319  format('zone i= ',i5,', j= ',i5,', D=(1,2)',
     +       ',T= "t= ',f12.4,' days; xy section ',i3,'"')
      do  i=1,nx
      i1=(i-1)*nyzt 
      write(14,8091) 
     +  (t2(i1+kz),wu(t2(i1+kz),ppn(i1+kz),qqn(i1+kz),modelwu,ts),
     +                                        c2(i1+kz),kz=iz,nyzt,nzt)
      enddo
      endif

      enddo             !next knoz
 7334 continue
c
c     file headers for heatflux across surfaces
c     ------------------------------------------
      write(44,8171)
 8171 format('Title="2D conductive heatflux top & bot; + heat leaving"',
     +     /,'variables="x","y","fluxtop W/m^2","fluxbot W/m^2"')
      write(54,8181)
 8181 format('variables="time","fluxleft (W)","fluxright (W)",',
     +       '"Sum_fluxleft (J)","Sum_fluxright (J)"',
     +       '"hflxtop (W)","hflxbot (W)","Wu_vol (m3)",',
     +       '"Keqleft (m/s)","Keqright(m/s)"',
     +   /,'zone i= 9999',',T="Boundary Heatfluxes"')
c
c     Tpeak and Tmin (normally for probability and life expectancy results)
c    ------------------------------------------------------------------------
      do i=1,nxy
      dz(i)=0.
      dztot(i)=0.
      bzqn(i)=0.      
      enddo
      dzmin=0.
      dzmax=0.
      write(41,8421) ny,nx,tdays
8421  format('title="xy Tpeak and Tmin, and dz"',
     +     /,'variables="x","y","Tpk","Tmin","dz","dztot"',
     +     /,'zone i= ',i5,', j= ',i5,',T= "t= ',f10.1,' days"')
      do  i=1,nx
      do  j=1,ny
      tpeak = -1.e33
      tmin  =  1.e33
      do  k=1,nzt
      i1=(i-1)*nyzt+(j-1)*nzt+k
      t2i1 = t2(i1)
      if(t2i1.gt.tpeak) tpeak = t2i1 
      if(t2i1.lt.tmin)  tmin = t2i1
      enddo 
c     nodexy  = (j-1)*nx+i
      nodexy2 = (i-1)*ny+j
      write(41,8092) x(i1),y(i1),tpeak,tmin,dz(nodexy2),dztot(nodexy2)
 8092 format(6e15.6)
      enddo
      enddo
      call flush(41)
c
c     initial t=0 total liquid water volume in domain
c     ----------------------------------------------
      vwtot=0.
      do l=1,net
        tavg=0.
        do i=1,8
         tavg=tavg+t2(in2(l,i))
        enddo
      tavg=tavg/8.
      vwtot=vwtot+vt(l)*por(l)*wu(tavg,pp(l),qq(l),modelwu,ts)    !check por(l) defined for unsat elements ?  
      enddo
      heatfluxl=0.0
      heatfluxr=0.0
      heatfluxrcum=0.
      heatfluxlcum=0.
      hflxtop=0.      
      hflxbot=0.
      ckequivl = cx(1)                !for t=0
      ckequivr = cx(nex)              !for t=0
      tdays=0.
      write(54,5758) tdays,heatfluxl,heatfluxr,heatfluxlcum,
     +        heatfluxrcum,hflxtop,hflxbot,vwtot,ckequivl,ckequivr
 5758 format(10e15.5)
c
c     t1dx.plt:
c     (co and tkl only print for 1d systems)
c     ----------
      write(31,8373) 
 8373 format('variables="x","Temp","Wu","Conc."')
      write(60,8473) 
 8473 format('variables="x","hflux_xl(W/m2)"')
      write(88,8378) 
 8378 format('variables="x","Co","tkl"')
c
      kny = knoy(1)
      if(kny.eq.0) kny=1
      do jnz=1,5
      if(knoz(jnz).eq.0.or.kny.eq.0) goto 9776
      write(31,8374) nx,knoz(jnz)
      write(88,8374) nx-1,knoz(jnz),kny
 8374 format('zone i= ',i5,',T = "t = 0 days, knoz= ',i5,
     +                                     ', knoy= ',i5,'"')
c
c     1dx profile print
c     -----------------
      do  i=1,nx
      i1=(i-1)*nyzt+((kny-1)*nzt)+knoz(jnz)
      kxy=(i-1)*ny+kny
c
c     depth to Ts, search up
c     ------------------------
      nbot=(i-1)*nyzt + ((kny-1)*nzt)+1
      ntop=nbot+nzt-1
      ztots=0.
      do k=1,nzt
      node = nbot+k-1
      if(t2(node).ge.ts) then
                 ztots=z(ntop)-z(node)
                 goto 5999
         endif
      enddo
c
c     interpolate between this node and underlying node  
5999  if(k.gt.1) then
      tdiff1 = t2(node)-t2(node-1)
      if(tdiff1.gt.0.0) ztots = z(ntop) - 
     +      (z(node-1) + (ts-t2(node-1))*(z(node)-z(node-1))/tdiff1)
      endif
      ztotsx(i) = ztots      
c
      write(31,4809) x(i1),t2(i1),wu(t2(i1),ppn(i1),
     +                                     qqn(i1),modelwu,ts),c2(i1)
 4809 format(7e15.6)
      
      enddo
      ztotsx(i) = ztots      
c
c     debug print at t=0
c     ---------------------
      if(ney.eq.1.and.nezt.eq.1) then
      do k=1,nex
      node = (k-1)*ny*nzt+(kny-1)*nzt+knoz(jnz)
      ielem = (knoz(1)-1)*nex*ney + (kny-1)*nex + k
      write(88,809) x(node),cpe,tkl(ielem)             !debug - only for 1dx systems
      enddo
      endif

 9776 continue
      enddo            !jnz

c
c     time t0: write dz, bzq stuff to file 61 t1dx_bzq x-profile data at top node row at knoy(1)
c  ---------------------------------------------------------------------------------------------
      write(61,8573) 
 8573 format('variables="x","dztot(m)","depth_Ts(m)","bz(m)",',
     +           '"bzf(W/m2)","bzq(m/s)","qz(m/s)"')

      write(61,8384) nx,nzt
 8384 format('zone i= ',i5,',T=" t=0 days, znode= ',i5,'"')
      do  i=1,nx
      i1=(i-1)*nyzt+((knoy(1)-1)*nzt)+nzt     
      kxy=(i-1)*ny+knoy(1)
      l = net-nexy 
     +      + (min(max((knoy(1)-1),1),ney-1)*nex) + (max(i,nex))       !element
      tavg=0.0
      do j=1,8
      tavg=tavg+t2(in(l,j))
      enddo
      tavg=tavg/8.
      ppl=pp(l)
      qql=qq(l)
      qz = vz(l)*por(l)*sw(l)*wu(tavg,ppl,qql,modelwu,ts)
      write(61,4809)
     +     x(i1),dztot(kxy),ztotsx(i),bz(kxy),bzf(kxy),bzqn(kxy),qz
      enddo

      
c
c ================================================================================      
      write(20,998) 
 998  format('title="1d vertical head, T and C1,C2 profile"',/,
     + 'variables="z","h(m)","T(degC)","C1","C2",',
     + '"Spn","Wu","Wu*Sw","dzl(m)"')

      write(11,9963) 
9963  format('title="1d vertical T profile: background (i,j)=(1,1)"',/,
     +        'variables="z","T(degC)","Wu"')

      write(40,3340) 
3340  format('title="Mass in matrix, fractures and column"',/,
     + 'variables="tdays","Mass in matrix (kg)","Mass in fractures",
     + "Mass in column"',/,'zone i=999')

      write(43,9918) nzt 
9918  format('title="1D vertical temperature profile - time plot"',/,
     +        'variables="t(d)","z(m)","T(degC)","Wu"',/,
     +        'zone i=',i5,',j= 9999,','T="at knox1,knoy1"') 
      write(49,9919) nx  
9919  format('title="1D horizontal temperature profile - time plot"',/,
     +        'variables="t(d)","x(m)","T(degC)","Wu"',/,
     +        'zone i=',i5,',j= 9999,','T="at knoy1,knoz1"') 

c     
c     print initial K's in flow grid in plane knoy(1) (or in plane 1) to file 51
c     ---------------------------------------------------------------------------
      write(51,7808)
 7808 format('title="2d K section "',/,
     +  'variables="x","z","Ksat","K(T).kr","log(K(T).kr)","porosity"',
     +                                    ',"Thermal K_pm","Bulk den"')
      write(51,7809) nz,nx 
 7809 format('zone i= ',i5,', j= ',i5,',T="Initial K, por, TK, rhob"')


c---
      do i=1,nx
      j=max(knoy(1),1)
      do k=1,nz

      node = (i-1)*ny*nz + (j-1)*nz +k       !flow node
      l=(k-1)*nex*ney + (j-1)*nex + i        !element associated with node
      if(i.eq.nx.and.k.eq.nz) l=ne-ney+1
      if(i.eq.nx.and.k.lt.nz) l=k*nex*ney-ney+1
      if(i.lt.nx.and.k.eq.nz) l=(k-2)*nex*ney+(j-1)*nex + i  
      cxx= cx(l)
      porel= por(l)
      tklel= tkl(l)
      por0el= por0(l)

c---
c get real K from Wu
c     zkrw=1.d0
      swe=wu(t2(map(node)),pp(l),qq(l),modelwu,ts)
c     zkrw = ((swe - p)/(1.d0-p))**4
c     zkrw = max(zkrw,1.e-6)
c
      alpha = alpha0 + por(l)*alphap*(rhob(l)/ps(l))*spe1(l)
      zkrw = fnzkrw(swe,pp(l),omega,porel,modelkr,por0el,alpha0,alpha)             !relative permeability function
      if(lkr1top .and. l.ge.(ne-nexy+1)) zkrw=1.d0                    !hardwire keep kr=1 for top surface elements
c
c     update temp. dependent conductivity:
c     adjust for ice saturation dependent relative k
c     check: rden should give true K's (not needed in flow routine -see user guide equations)
c     ---------------------------------------------
c     cxx = cxx*den(t2(map(node)),lmass,lage)/1000.
c     cxx2=(cxx/rvisc(t2(map(node)),lmass,lage)) * zkrw
      cxx2=(cxx*(den(t2(map(node)),lheat,lmass,lage)/1000.)
     +                 /rvisc(t2(map(node)),lheat,lmass,lage)) * zkrw
c
c     write(6,4776) swe,zkrw,t2(map(node)),
c    +                   den(t2(map(node)),lheat,lmass,lage),
c    +                 rvisc(t2(map(node)),lheat,lmass,lage)
c4776 format(10x,'initial swe,zkrw,t2,den,rvisc: ',5e15.5)
c
      write(51,856) x(map(node)),z(map(node)),
     +              cxx,cxx2,log10(cxx2),porel,tklel,rhob(l)
c856  format(7e16.7)
      end do
      end do     
      call flush(51) 
c     ----------------------------------------------------
c
c     count 2D xy plane fractures attached to knoz(1)
c     ------------------------------------------------
      nfknoz1=0
      if(knoz(1).gt.0) then
      do lxy=1,nexy
      iel = (knoz(1)-1)*nexy + lxy
        do i=1,nfrac             !search all 2D fracs for those on top or bottom of this element
        if(ifracl(i).eq.iel 
     +     .and. ifdim(i).eq.2
     +     .and. (lvert(i).eq.5.or.lvert(i).eq.6)) then
               nfknoz1=nfknoz1+1
               ifrac2dz1(nfknoz1) = i
        endif
        enddo
      enddo
      endif
c
c     get connections for use in distributing dz frost heave/thaw settlement
c     ---------------------------------------------------------------------------
      do i=1,nxy
      link(i)=0
      enddo
      do i=1,nexy
      do k=1,4
      nodexy = (in2(i,k)+nzt-1)/nzt
      link(nodexy) = link(nodexy) + 1
      enddo
      enddo
c     write(6,*) (link(i),i=1,nxy)
c
c -------------------------------------------------------------------
c     plot co vs.  Temp and dCo vs. Temp
c -----------------------------------------------------------------------
      open(99,file='smoker_den_visc_T.plt')
      write(99,8221)
 8221 format('variables="Temp","pw(kg/m3)","rden(T,C=0)","rden(T,C=1)",
     + "rvisc_w(T)","Wu(T)","dWu(T)","tk(T)","Co(T)"',
     +  /,'zone i= 500, T="(at element 1)"')
      cavg0=0.0
      cavg1=1.0
      do i=1,500
      tt= -5. + float(i)*(6./100.)
      porsw = por(1)*sw(1)
      wu0 = wu(tt,pp(1),qq(1),modelwu,ts)
      ww0=porsw * wu0
      wi0=porsw * (1.d0-wu0)
      wm=(1.d0-por(1))  
      den0 = den(tt,lheat,lmass,lage)
      rden0 = rden(tt,cavg0,gamma,lheat,lmass,lage)
      rden1 = rden(tt,cavg1,gamma,lheat,lmass,lage)
      rvisc0 = rvisc(tt,lheat,lmass,lage)
      cpf0=cf*den0
      cpe0=ww0*cpf0 + wi0*cpi + wm*cpsm(1)
c
c      thermal conductivity tk trans1
c      -------------------------------
c      either linear implicit ...
c      tkl(l) = ww0*tclw + wi0*tcli + wm*tclm(l)            ! geometric mean    !hardwire Interfrost
c      linear - centre-weighted ....
c      tcl01 =  ww0*tclw + wi0*tcli 
c      tkl1 = tcl01 + wm*tclm(1)             !centre-weighting of TK
c
c      or sqrt-average implicit ...
c      ww=por(1)*sw(1)*wu(tt,p,q,modelwu,ts)
c      wi=por(1)*sw(1)*(1.-wu(tt,p,q,modelwu,ts))
c      wm=(1.-por(1))
      tk = (ww0*sqrt(tclw)+wi0*sqrt(tcli)+wm*sqrt(tclm(1)))**2     !more physically realistic
c     dlh=por(l)*den(tav1,lmass,lage)*wlh*dwu(tav1,p,q,modelwu,ts)
      dwu0 = dwu(tt,pp(1),qq(1),modelwu,ts)
      dlh0=por(1)*pi*wlh*dwu0            !pi (ice density) 
      cot00=(cpe0+dlh0)                !heat capacity 
c
      write(99,*)  tt,den0,rden0,rden1,rvisc0,wu0,dwu0,tk,cot00
      enddo
      close(99) 

      write(59,9849) 
 9849 format('variables="x","z","Pe_x","Pe_z","Co_x","Co_z"')
c
c     save initial fc3 array for applying source decay2
c     --------------------------------------------------
c     do i=1,nnt
c     fc30(i)=fc3(i)
c     enddo
c -------------------------------------------------------------------
      call flush(8)
      call flush(9)
      call flush(14)
      call flush(16)
      call flush(21)
      call flush(25)
      call flush(31)
      call flush(60)
      call flush(61)
      call flush(62)
c
c     begin new time loop with new dt ...
c     stop if dt = 0
c     ######################################
c
c     temperature solution is plotted  every 'kplot' time steps
c     flow solution is updated every 'kflow' time steps
c     3D temperature solution printed every 'k3d' time steps
c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
c
c     time step data are in entered in days
c     -------------------------------------------------------------
c
      igo=0
      igof=0
      igo3dv=0
      it= 0
      itot= 0
      momt = 0
      numpt43=0
      numpt49=0
      flxin=0.
      flxout=0.
      dflxin=0.
      dflxout=0.
 1000 continue
      read(5,*,err=5001,end=5001) 
     +     time0,time1,dt,kplot,kflow,k3d,kmom,more_time
  121 format(3f10.0,3i5)
c
      time0=time0*86400.
      time1=time1*86400.
      dt=dt*86400.
      write(6,122)
     + time0/86400.,time1/86400.,dt/86400.,kplot,kflow,k3d,kmom
  122 format(//10x,'time interval from',e15.7,' to ',e15.7,/10x,
     +             'time step dt =    ',e15.7,' (days)'/10x,60('-'),
     +  /10x,'2d plot frequency =      ',i8,
     +  /10x,'flow update frequency =  ',i8,
     +  /10x,'3d plot frequency =      ',i8,
     +  /10x,'moment/tpeak frequency = ',i8)
c
      if(time1.le.time0) then
      write(6,1555)
 1555 format(/10x,'error detected in time step data',
     +       /10x,'check and rerun ... program stopping ')
      stop
      endif
c
c     also specify increment  for time adjustment of {fc} 
c     (time-variable dirichlet potentials)
c     --------------------------------------------------------------
      hinc=0.
      rinc=0.
      read(5,*,err=3031) hinc,rinc,rtdiff,bzfinc
      write(6,124) hinc,rinc,rtdiff,bzfinc
  124 format(/10x,'head increment at boundary (+/-) :    ',f7.3,
     +       /10x,'recharge increment this interval (*): ',f7.3,
     +       /10x,'recharge temperature difference (wrt Tair) ',f7.3,
     +       /10x,'bzf multiplier for surface heat flux  ',f7.3)
c
c     adjust transient dirichlet boundary potentials ...
c     allows for uniform rise in watertable - (constant gradient)
c     allows for transient dirichlet air temperatures 
c     allows for transient recharge across top
c     note: increments to head and recharge are cumulative
c     -----------------------------------------------------------
      if(kb1(6).ne.6) then       !do not apply if top is cosine wt; hinc used as hinc*amp
      do 650 i=1,nn
      fc(i)=fc(i)+hinc
  650 continue
      endif

      if(.not.lbzin) then              !file bzin does not exist ...
      k=0
      do i=nz,nn,nz
      k=k+1
      fb(i)=fb(i)*rinc
      fbff(k) = fbff(k)*rinc
      enddo
      endif
      
c
c     read variable surface air temperature, defined by a surface patch, 
c     x node from i1 to i2, y node from j1 to j2, vtemp=surface temp
c     surface temp defined later by SURFAT() will not override this patch
c     vtemp must be > -999 to function, if vtemp=-999, uses function surfat
c     if kb(6)=1, this patch is Dirichlet; if kb(6)=4, this patch is leaky
c     -------------------------------------------------------------------
c     write(6,9303) (ic2(i),fc2(i),i=nzt,nnt,nzt)
c9303 format('ic2,fc2:',4(i5,e12.3))
      do 1166 i=1,nnt
 1166 ivt(i) = 0
 1165 read(5,*) i1,i2,j1,j2,vtemp,more
 1167 format(4i5,f10.0)

      write(6,1169) i1,i2,j1,j2,vtemp
 1169 format(/10x,'patch surface temp: ',4i6,f12.3)
      if(i1*i2*j1*j2 .eq.0) goto  7670      
      do 1168 i=i1,i2
      do 1168 j=j1,j2
      node= (i-1)*nyzt + j*nzt
      if(node.gt.0.and.node.le.nnt) then
c
c  check if node is located on 1st or 4th-type face ...
c  --------------------------------------------------------
         if(ic2(node).ne.1) write(6,9012) i,j
 9012    format(10x,'boundary patch node at i,j: ',2i7,' is leaky')
         if(ic2(node).eq.1) write(6,9013) i,j
 9013    format(10x,'boundary patch node at i,j: ',2i7,' is Dirichlet')
         call flush(6)
         ivt(node) = 1
         fc2(node) = vtemp
      else
      write(6,7669) i,j,node
 7669 format(/10x,'warning - surface temp patch i,j,node ',3i7,
     +       /10x,'is not valid',
     +       /10x,'program stopping ...')
      stop
      endif

 1168 continue
      if(more.gt.0) goto 1165

 7670 continue

c     read internal source elements, source flux, and decay term
c     July 1995  (decay units in 1/days, flux in J/m^2 s = W/m^2)
c     read element x-range (i1-i2) y-range (j1-j2) and z-range (k1-k2)
c     warning: for heater - should still be only 1- or 4- columns
c     otherwise, if kmass=0 can use these source elements to compute mass
c     fluxin1=W/m2 ... for thermal wells - flux applied to outer faces
c     fluxin2=W/m3 ... for internal heat ... flux applied to volume
c     ----------------------------------------------------------------
      numis=0
      numiscol=0
      decay2=0.
      sumallheaters=0.
      sumallheatersv=0.
      sumheatflux=0.
      wpermsum=0.
 489  read(5,*,err=9922) i1,i2,j1,j2,k1,k2,fluxin1,fluxin2,decay2,more
c      print*, 'test:', i1,i2,j1,j2,k1,k2,fluxin1,fluxin2,decay2,more
      call flush(6)

      if(i2.gt.nex .or. j2.gt.ney .or. k2 .gt.nezt   
     +   .or. (((j2-j1).gt.1 .or.(i2-i1).gt.1)
     +         .and.kmass.eq.0.and.fluxin1.gt.0.)) then
      write(6,688)
 688  format(/10x,'error defining internal heat source element range',
     +       /10x,'check for 1- or 4-column heater, i2.le.nex etc. ...',
     +       /10x,'program stopping ...')
      stop
      endif
      if(fluxin1.gt.0.d0 .and. fluxin2.gt.0.d0) then
      write(6,6688)
6688  format(/10x,'error... should not have both fluxin1 & fluxin2 > 0',
     +       /10x,'check fluxin1, fluxin2...',
     +       /10x,'program stopping ...')
      stop
      endif

      if(i1.lt.1.or.i2.lt.1.or.j1.lt.1.or.j2.lt.1.or.k1.lt.1.or.k2.lt.1)
     + then 
c      elimit = float((i1*i2)*(j1*j2)*(k1*k2))
c      print*,'elimit= ',elimit
c      if(elimit.lt.1.) then                   !detect 0 entries
      write(6,588)
 588  format(/10x,'0 entry for internal source elements, continuing...')
      else 
      numiscol=numiscol+1
      write(6,488) i1,i2,j1,j2,k1,k2,fluxin1,fluxin2,decay2
 488  format(/10x,'internal heat source this time interval ... ',
     +       /10x,'element range:            ',6i5,
     +       /10x,'heat flux (W/m^2):        ',e12.5,
     +       /10x,'heat flux (W/m^3):        ',e12.5,
     +       /10x,'exponential source decay: ',e12.5)
      call flush(6)
c
c     get source elements and nodes
c     hardwired for outer faces only of first and last columns
c     sarea and sumarea are the exposed (outer) areas of the column
c     so this will work only if source is one column or 4
c     fluxv should work for any number of elements
c     ---------------------------------------------------------
      sumarea = 0.d0
      sumvol = 0.d0
      zcol=0.d0
      numezcol=k2-k1+1
      zcolavg=0.d0
      do  i=i1,i2
      do  j=j1,j2
      zcol=0.
      do  k=k1,k2
      numis=numis+1
      lxyz = (k-1)*nex*ney + (j-1)*nex + i
      ifl(numis)=lxyz
      flux(numis) = fluxin1
      fluxv(numis) = fluxin2
      zcol=zcol+ezl(lxyz)
c 
c     only for fluxin1 ... W/m2
c     -------------------------
      if(fluxin1.gt.0.d0) then      
      if(i2.eq.i1 .and. j2.eq.j1) then        !for 1-column
           nflux(numis,1)=in2(ifl(numis),1)
           nflux(numis,2)=in2(ifl(numis),2)
           nflux(numis,3)=in2(ifl(numis),3)
           nflux(numis,4)=in2(ifl(numis),4)
           nflux(numis,5)=in2(ifl(numis),5)
           nflux(numis,6)=in2(ifl(numis),6)
           nflux(numis,7)=in2(ifl(numis),7)
           nflux(numis,8)=in2(ifl(numis),8)
      endif
      
      if(i2.gt.i1 .and. j2.gt.j1) then        !for 4-columns
      if(j.eq.j1) then
           nflux(numis,1)=in2(ifl(numis),1)
           nflux(numis,2)=in2(ifl(numis),2)
           nflux(numis,3)=in2(ifl(numis),5)
           nflux(numis,4)=in2(ifl(numis),6)
      endif
      if(j.eq.j2) then
           nflux(numis,1)=in2(ifl(numis),3)
           nflux(numis,2)=in2(ifl(numis),4)
           nflux(numis,3)=in2(ifl(numis),7)
           nflux(numis,4)=in2(ifl(numis),8)
      endif
      if(i.eq.i1) then
           nflux(numis,5)=in2(ifl(numis),4)
           nflux(numis,6)=in2(ifl(numis),1)
           nflux(numis,7)=in2(ifl(numis),5)
           nflux(numis,8)=in2(ifl(numis),8)
      endif
      if(i.eq.i2) then
           nflux(numis,5)=in2(ifl(numis),2)
           nflux(numis,6)=in2(ifl(numis),3)
           nflux(numis,7)=in2(ifl(numis),6)
           nflux(numis,8)=in2(ifl(numis),7)
      endif
      endif                                   !end 4 columns
      endif                               !end W/m2

      if(fluxin2.gt.0.d0) then            !for W/m3
        nflux(numis,1)=in2(ifl(numis),1)
        nflux(numis,2)=in2(ifl(numis),2)
        nflux(numis,3)=in2(ifl(numis),3)
        nflux(numis,4)=in2(ifl(numis),4)
        nflux(numis,5)=in2(ifl(numis),5)
        nflux(numis,6)=in2(ifl(numis),6)
        nflux(numis,7)=in2(ifl(numis),7)
        nflux(numis,8)=in2(ifl(numis),8)              
      endif

      l=ifl(numis)
      sarea = (exl(l)*ezl(l) + eyl(l)*ezl(l))
      if(nflux(numis,3) .eq. in2(l,3)) sarea = sarea*2.d0    !single column: two-face area*2
      sumarea = sumarea + sarea
      sumvol = sumvol+exl(l)*eyl(l)*ezl(l)
c
c     write(6,4959) numis,l,nflux(numis,3),in2(l,3),exl(l),eyl(l),ezl(l),
c    +              sarea
c4959 format(10x,'numis,l,nflux(),in2(),exl(l),eyl(l),ezl(l),sarea.. ',
c    +           4i6,4e12.3)
      
      enddo
      zcolavg=zcolavg+zcol
      enddo
      enddo
      ncolumns = (i2-i1+1)*(j2-j1+1)
      zcolavg=zcolavg/float(ncolumns)

      sumallheaters = sumallheaters + sumarea
      sumallheatersv = sumallheatersv + sumvol
      sumheatflux = sumheatflux + fluxin1*sumarea +fluxin2*sumvol
      wperm = fluxin1*sumarea/zcolavg                 !W/m for this heater group
      wpermsum = wpermsum + wperm

      write(6,5485) numezcol,zcol,sumarea,sumvol,wperm
 5485 format(10x,'number of element layers in this heater ... ',i5,
     +      /10x,'vertical height of heater column .......... ',e12.4,
     +      /10x,'total surface area of this heater column .. ',e12.4,
     +      /10x,'total volume of this heater column ........ ',e12.4,
     +      /10x,'average heat flux per m length (W/m) ...... ',e12.4)
      
      if(more.gt.0) goto 489
      endif

      write(6,5768) numiscol,numis,sumallheaters,sumallheatersv,
     +              sumheatflux,wpermsum
 5768 format(/10x,'heater column input data:',
     +      /10x,'total number of heater columns .... ',i5,
     +      /10x,'cumulative # of elements in all heater columns: ',i5,
     +      /10x,'total surface area of all heater columns ... ',e12.4,
     +      /10x,'total volume of all heater columns ......... ',e12.4,
     +      /10x,'total heat flux over all columns (W) ....... ',e12.4,
     +      /10x,'total flux/m length for all columns (W/m) .. ',e12.4)
c
c788  continue
c
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     enter source/sink parameters for this time interval:
c     ****************************************************
c     zero source/sink arrays: (numss=number of source/sink nodes)
c     read source/sink location (on transport grid): i (x grid position) ,j (y grid position),
c     and vertical interval k1,k2. pumping rates entered in m^3/s per node
c     so total pumping rate = Q = pqq*(k2-k1+1)
c     --------------------------------------------------------------------
      do 240 i=1,nnt
      pq(i)=0.
      tq(i)=0.
      cq(i)=0.
 240  continue      
      numss=0
c
  241 read(5,*) i,j,k1,k2,pqq,tqq,cqq,more
c     print*,'well::: i,j,k1,k2,pqq,tqq,cqq,more',i,j,k1,k2,pqq,tqq,more
c
c     determine node number range on transport grid
c     (but pq,tq must fall within flow grid)
c     ---------------------------------------------
      nq1 = (i-1)*nyzt + (j-1)*nzt + k1
      nq2 = nq1 + k2 - k1
c     print*,'nq1,nq2',nq1,nq2
      if(nq1.gt.0 .and. nq2.le.nn 
     + .and. k1.le.nz .and. k2.le.nz .and. j.le.ny .and. i.le.nx) then
      write(6,243) nq1,nq2,pqq,tqq,cqq
  243 format(/10x,'revised source/sink nodes ',i7,' to ',i7,
     +       /10x,'source/sink strength:       ',e12.4,
     +       /10x,'            temperature:    ',e12.4,/,
     +       /10x,'            concentration:  ',e12.4,/)
      do 244 i=nq1,nq2
         numss=numss+1
         kssn(numss)=i
         pq(i)=pqq
         tq(i)=tqq
         cq(i)=cqq
 244  continue         
      else 
      if(nq1.gt.nn .or. nq2.gt.nn) then
        write(6,419) i,j,k1,k2
 419    format(/10x,'error in source/sink nodal range ... ',
     +         /10x,4i7,
     +         /10x,'nodal range must be within the flow grid ... ',
     +         /10x,'program stopping ... ',/)
               stop
      endif
      endif
      if(more.gt.0) goto 241
      lss=numss.ne.0

      call flush(6)
c
c     set time ...
c     ------------
      time = time0
      dt1 = 1./dt
      dt2 = dt/2.
c
c     time loop with constant dt ...
c     ##############################
c
      itdt= 0
 1001 continue
      dumx=999.
      itdt= itdt + 1
      it = it + 1
      numit=it
      time = time + dt
      tdays= time/86400.
      write(6,123) it,tdays
c     write(*,123) it,tdays
  123 format(/10x,'time step #',i6,4x,'time = ',e15.7,' days',
     +       /10x,45('#'),/)
      print *,'         beginning solution for time= ',tdays,' days'

      if(it.gt.maxit) then
      write(6,4857) it,maxit
 4857 format(/10x,'error ... #time steps (',i7,')',' > maxit (',i7,')',
     +       /10x,'program stopping ...')
      stop
      endif
c
c     3D initial mass calc for kmass=1
c     --------------------------
c     mass in porous blocks:
c     -------------------------------------------------------
      if((kmass.eq.1 .or. kmass.eq.3.or.kmass.eq.4) 
     +         .and. kprt.eq.1 .and. it.eq.1) then
      xmass0=0.d0
      do l=1,net
      xe=0.d0
      do k=1,8
      xe = xe + c2(in2(l,k))
      enddo
      xmass0 = xmass0 + xe*vt(l)*por(l)/8.d0
      enddo
c
c     mass in fractures:
c     -----------------------
      fmass0 = 0.d0
      do kf=1,nfrac
      l = ifracl(kf)
      if(ifdim(kf).eq.2) then
      call frac_plane(maxfrac,maxne,lvert,inline,in,kf,l,fdimx,fdimy, 
     +             exl,eyl,ezl,rexel,sarea)
      fmass0 = fmass0 + (c2(map(inline(1)))+c2(map(inline(2)))
     +                +  c2(map(inline(3)))+c2(map(inline(4))))
     +                *sarea*xarea(kf)/4.d0
      endif

      if(ifdim(kf).eq.1) then
      call frac_line(maxfrac,maxne,lvert,inline,in,kf,l,fracdim,
     +             exl,eyl,ezl)
      fmass0 = fmass0 + (c2(map(inline(1)))+c2(map(inline(2))))
     +                            *fracdim*xarea(kf)/2.d0
      endif
      enddo

c     mass within internal source elements (even if fluxin=0)
c     just used for convenience - use isource to identify elements around injection well
c     ----------------------------------------------------------------------------------
      wmass0 = 0.d0
      do i=1,numis
      l=ifl(i)
      xe=0.d0
      do k=1,8
      xe = xe + c2(map(in(l,k)))
      enddo
      wmass0 = wmass0 + xe*vt(l)*por(l)/8.d0
      enddo
      endif
c     ================  end initial mass calc) ===================

c
c     get background surface air temperature for this interval
c     applies for surface area if given temp = -999. (either Dirichlet or Type-4)
c     for either Dirichlet or leaky, patch source air/ground temp was given as vtemp and passed to fc2
c     temp gradient in space gradx only applies for ktair = 0 (not if using tair file)
c     ------------------------------------------------------------------------------------------------
      if(ktair.eq.0) then
         ttemp =
     +   surfat(surfmin,amp,period,phase,tdays,cutoff,gradt)
       else
         if(tdays.gt.timeair(nair)) then
                 ttemp = tempair(nair)
                 goto 1279 
         endif
         do i=1,nair-1
         if(tdays.gt.timeair(i) .and. tdays.le.timeair(i+1)) then
                 dtair1 = tdays-timeair(i)
                 dtair2 = timeair(i+1) - timeair(i)
                 dtempair1 = tempair(i+1)-tempair(i)
                 ttemp = (dtair1*dtempair1/dtair2)+tempair(i)      !linear interpolation
                 goto 1279
         endif
         enddo
      endif                    !end ktair
 1279 continue     


c     if(leak2) write(6,356) ttemp,ttemp+rtdiff
c     if(.not.leak2) write(6,396) ttemp
c 356 format(/10x,'reference surface air temperature (x=0) =    ',f7.2,         ! at x=0 since no gradx*x correction to surfat
c    +       /10x,'reference recharge water temperature =       ',f7.2)
c 396 format(/10x,'reference surface ground temperature (x=0) = ',f7.2)
c     do 651 i=nzt,nnt,nzt
c     if(vtemp.lt.-998..and.ivt(i).eq.1)    fc2(i)=ttemp     !use surfat function for boundary temperature (without gradx)
c     if(leak2.and.ivt(i).eq.0)  fc2(i)=ttemp + x(i)*gradx     !update Leaky air temps outside of source (with gradx) 
c 651 continue    
c
c     get conductance layer thickness bz and bzflx
c     if bzin.data file exists, will use this transient bz and bzflx data
c     otherwise will use bz from input sine function and surfat function
c     also for bzq and bztdiff
c     --------------------------------------------------------------------

      if(lbzin) then
         if(tdays.gt.tbzbdy(numbzbdy)) then                !if t > tmax in bzin file ... just use last value
             do ibz=1,nbzz(numbzbdy)                        !# of spatial zones at last bz-time
             do i=ix1bz(numbzbdy,ibz),ix2bz(numbzbdy,ibz)
             do j=1,ny
               nodexy = (i-1)*ny + j
               bz(nodexy) = bzbdy(numbzbdy,ibz)
               bzf(nodexy) = bzflx(numbzbdy,ibz)
               bzqn(nodexy) = bzq(numbzbdy,ibz)
             enddo     !ny
             enddo     !ix1-ix2
             enddo     !nbzz zones     

c  all nodes for recharge ...
             k=0
             do i=nz,nn,nz
               k=k+1
               fb(i)=bzqn(k)*ara(k)
               fbff(k) = bzqn(k)
             enddo

        else                                                !apply linear interpolation between times

                do i=1,numbzbdy-1
         if(tdays.gt.tbzbdy(i) .and. tdays.le.tbzbdy(i+1)) then     !tdays check
                 dtbz1 = tdays       - tbzbdy(i)
                 dtbz2 = tbzbdy(i+1) - tbzbdy(i)
                 ftime = dtbz1/dtbz2                        !fractional time 0-1
                 do ibz=1,nbzz(i)                           !# space zones
                 do i2=ix1bz(i,ibz),ix2bz(i,ibz)    !x-index
                 do j=1,ny
                    nodexy = (i2-1)*ny + j

                 dbzbdy1 = bzbdy(i+1,ibz) - bzbdy(i,ibz)   
                 dbzflx1 = bzflx(i+1,ibz) - bzflx(i,ibz)  
                 dbzq1   = bzq(i+1,ibz) - bzq(i,ibz)
                 bz(nodexy)  = bzbdy(i,ibz) + (ftime*dbzbdy1)     !linear interpolation 
                 bzf(nodexy) = bzflx(i,ibz) + (ftime*dbzflx1)     !linear interpolation
                 bzqn(nodexy) = bzq(i,ibz) + (ftime*dbzq1)        !linear interpolation New sept 2017 for bzq variable in space

              enddo     !ny
             enddo     !ix1-ix2
            enddo     !nbzz zones in space at this time                 
c
c        bzq and bztdiff constant in space - interpolate in time
c        -------------------------------------------------------
         k=0
c          dbzq     = bzq(i+1)-bzq(i)
          dbztdiff = bztdiff(i+1)-bztdiff(i)
          rtdiff   = bztdiff(i) + ftime*dbztdiff
         do i2=nz,nn,nz
          k=k+1
          fbff(k) = bzqn(k)
          fb(i2)  = bzqn(k)*ara(k)     !bugfix () Aug 31,2017;   another old bug: bzq(i) .. what was i ? time. ok .
         enddo

         goto 1290                                 !  done, get out
        endif       !end tdays check

        enddo       !try next bztime interval
       endif        !end tdays > last value

1290  continue
      endif        !ierr6 if bzin.data exists

c
c     adjust bzf(nodexy) with bzfinc
c     -------------------------------
      if(lbzin) then
      do i=1,nx
      do j=1,ny
      nodexy = (i-1)*ny + j
      bzf(nodexy) = bzf(nodexy) * bzfinc
      enddo
      enddo
      endif
c
c     print
c     ---------
      if(leak2) write(6,356) ttemp,rtdiff,fbff(1)
  356 format(/10x,'reference surface air temperature (x=0) =    ',f7.2,         ! at x=0 since no gradx*x correction to surfat
     +       /10x,'reference temp. difference (Tq - Tair) =     ',f7.2,
     +       /10x,'recharge flux m/s =                          ',e10.2)
      if(.not.leak2) write(6,396) ttemp
  396 format(/10x,'reference surface ground temperature (x=0) = ',f7.2)
      do i=nzt,nnt,nzt
      if(vtemp.lt.-998..and.ivt(i).eq.1)    fc2(i)=ttemp     !use surfat function for boundary temperature (without gradx)
      if(leak2.and.ivt(i).eq.0)  fc2(i)=ttemp + x(i)*gradx     !update Leaky air temps outside of source (with gradx) 
      enddo   

c
c     get bz from surfat function
c     assume bzf = 0
c     -------------------------------
      if(.not.lbzin) then
      bzsurf = 
     +surfat(surfmin2,amp2,period2,phase2,tdays,cutoff2,gradt2)
      do i=1,nxy
      bz(i) = bzsurf
      bzf(i) = 0.
      enddo
      endif


c
c     get conductance layer saturation
c     ---------------------------------
      sat =
     +surfat(surfmin3,amp3,period3,phase3,tdays,cutoff3,gradt3)

      write(6,948) bz(1),bz(nxy/2),bz(nxy)
 948  format(10x,'surface layer thickness bz 1,mid,nx  = ',3f8.3)
      write(6,949) bzf(1),bzf(nxy/2),bzf(nxy)
 949  format(10x,'surface heat flux bzflx    1,mid,nx  = ',3e9.2)
c      write(6,948) (bz(i),i=1,nxy,ny)
c 948  format(10x,'surface layer thickness bz    = ',/(2x,15f8.2))
c      write(6,949) (bzf(i),i=1,nxy,ny)
c 949  format(10x,'surface layer heat flux bzflx = ',/(2x,15f8.4))
      write(6,947) sat
 947  format(10x,'surface layer saturation sat = ',f7.2)

      do i=1,nxy
      if(leak2.and.(bz(i).le.0.d0.or.sat.le.0.d0)) then
      write(6,4849) 
 4849 format(/10x,'error: bz or sat <= 0; check sin curve parameters',
     +       /10x,'or check bzin.data file for bz <= 0',
     +       /10x,'program stopping ...')
      stop
      endif
      enddo
c
c
      lplot  = ((itdt/kplot) * kplot.eq.itdt)
      lflow  = ((itdt/kflow) * kflow.eq.itdt)
      l3dplt = ((itdt/k3d)   * k3d.eq.itdt)
      lmom   = ((itdt/kmom)  * kmom.eq.itdt)
      lmom2  = ((itdt/(kmom*2))  * (kmom*2).eq.itdt)      
c
      vol=vol1
c
c  ===================================================================== 
c
c    get head value for left bdy at this time step - file hbdy.data
c    intended for cosine curve      
c    if data exist, they will be used, 
c     ... and if top is type 6, then'slope' from input file is applied
c             and right bdy heads are adjusted 
c    ---------------------------------------------------------------------
c      
      if(ierr4.le.0 .and. (.not.lwtfgt0)) then
         if(tdays.gt.thbdy(numhbdy)) then
                 hleft = hbdy(numhbdy)
                 goto 1289 
         endif
         do i=1,numhbdy-1
         if(tdays.gt.thbdy(i) .and. tdays.le.thbdy(i+1)) then
                 dth1 = tdays-thbdy(i)
                 dth2 = thbdy(i+1) - thbdy(i)
                 dhbdy1 = hbdy(i+1)- hbdy(i)
                 hleft = (dth1*dhbdy1/dth2)+hbdy(i)      !linear interpolation
                 goto 1289
         endif
         enddo
 1289 continue    
c
c     extend across watertable
c     ---------------------------
      pie = 4.d0 * datan(1.d0)
      plc = 2.d0 * wpc * pie/x(map(nn))
      pls = 2.d0 * wps * pie/x(map(nn))
c
c     regenerate wt array ... based on new hleft
c     Note: will not work with multiple wt functions along x (kb(6)=6)
c     (lwtfgt0)
c     ----------------------------------------------------------------
      do i3d=nz,nn,nz
      xx = x(map(i3d))
      cosx = dcos(plc*xx+phasewt)
      sinx = dsin(pls*xx+phasewt)
      fc(i3d) =  hleft + slope*xx + ampcos*cosx + ampsin*sinx
      enddo
c
c     adjust right type-1 bdy heads 
c     -------------------------------
      hright = fc(nn)
      do i=1,nyz
      node = (nx-1)*nyz + i
      if (ic(node).eq.1) fc(node)=hright
      enddo
      slope = (hright-hleft)/x(map(nn))
      write(6,1299) hleft,hright,slope
1299  format(10x,'left,right bdy head,slope: ',3f10.4)

      endif
c
c     regenerate fc using hinc as a multiplier for ampcos and ampsin for kb1(type6 on top)
c     (hardwire) uses wh so will not apply if hbdy.dat file is used to update hleft.
c     and will not work for multiplt wt functions in x (lwtfgt0)
c     ------------------------------------------------------------------------------------
      if(kb1(6).eq.6 .and. (.not.lwtfgt0)) then
      write(6,3987) 
 3987 format(/10x,'regenerating fc array ... ')      
      pie = 4.d0 * datan(1.d0)
      plc = 2.d0 * wpc * pie/x(map(nn))
      pls = 2.d0 * wps * pie/x(map(nn))
      do i3d=nz,nn,nz
      xx = x(map(i3d))
      cosx = dcos(plc*xx+phasewt)
      sinx = dsin(pls*xx+phasewt)
      fc(i3d) = wh + slope*xx + hinc*ampcos*cosx + hinc*ampsin*sinx

c  hardwire Domenico wt
c      BB1 = cosh(pie*z(nn)/x(nnt))
c      AA1 = x(nnt)/2. + BB1
c      fc(i3d) = AA1-BB1*cos(pie*xx/x(nnt))


      enddo
      write(6,1399) hinc
1399  format(10x,'wt amplitude updated by *hinc:',e10.2)
      endif 
c ____________________________________________________________
c ____________________________________________________________
c
c     begin non-linear iteration loop:
c     ################################
c _____________________________________________________________
c
      kit  = -1
 500  kit = kit + 1

      write(6,256) kit
  256 format(/10x,'non-linear iteration # (kit): ',i5,/10x,35('-'))
      update_flow=.false.

c     skip flow solution for kcntrl=2 (uses uniform vx,vy,vz)
c     --------------------------------------------------------
      if(kcntrl.eq.2 .or. ivel.ne.0) goto 502
c
c     skip flow if already converged - hardwire
c     ------------------------------
      if(dumx.le.ccp) goto 502

      update_flow=.false.
      if(itdt.lt.3)   update_flow=.true.   !always for first 3 time steps (of each new time interval)
c     if(itdt.lt.3)   update_flow=.true.
c !!! hardwire: after third time step of each time interval, don't do flow:
c     if(itdt.gt.3 .and. (ss .le. 1.e-9)) update_flow=.false.
c !!! hardwire: after first time step of entire solution, don't do flow:
c      if(it.gt.2) update_flow=.false.
c      if(itdt.lt.3) update_flow=.true.
      if(kcntrl.eq.1)  update_flow=.true.     !for flow solution only
      if(kcntrl.eq.3)  update_flow=.true.     !for transient flow solution only      
      if(lflow) update_flow=.true.
      if(gamma.gt.1.0e-5) update_flow=.true.
      if(phi.lt.999.) update_flow=.true. 

      if(.not. update_flow) goto 502

c
c     begin watertable deformation loop:
c     ##################################
      jit=-1
  501 jit=jit+1
c
c     write(6,257) jit
c     write(*,257) jit
c 257 format(10x,'watertable iteration # ',i5,/10x,28('-'))
c
      itot= itot+1
c
c     solve for flow ...
c     assembles with temperature {t1}, heads {u0}, solves for {u2}
c     ============================================================
c
c     tsurfmin=999.
c     tsurfmax=-999.
c     do i=1,nn
c     tsurfmin = min(tsurfmin,t2(i))
c     tsurfmax = max(tsurfmax,t2(i))           
c     enddo
c
c     write(6,3375) tsurfmin,tsurfmax
c3375 format(/10x,'call flow min,max temps : ',2e15.5)

      call flow(maxn,maxnn,maxne,nf,maxna,laa,np1,xi,yi,se,
     + pe,x,y,z,u2,u0,fc,fb,fx,ic,lc,cx,cy,cz,in,in2,
     + a,iaa,ind,ib,aa,nn,n,nb,ne,kprt,
     + dt1,fst,t0,t2,c0,c2,ss,pq,ag,hag,hxx,hyy,hzz,htmf,hgz,kint,
     + fast,fpcgt,tdays,map,exl,eyl,ezl,nexy,por,gamma,pi,pw,vt,
     + ckl,xarea,ifracl,lvert,nfrac,maxfrac,nex,ney,nez,ifdim,
     + hx,hy,ht,hgz2,kblck,pp,qq,modelwu,lheat,lmass,lage,omega,modelkr,
     + p2,q2,ts2,omega2,modelwu2,modelkr2,lkr1top,
     + por0,alpha0,alphap,rhob,ps,spn0,spn1,spe1)
c
      write(6,233) 
 233  format(10x,'flow solution complete ')
c     write(6,9883) (fc(i),i=nz,nn,nz)
c9883 format(10x,'fc array top nodes...',
c    +   /,(10e12.3))
c
c     check convergence of watertable:
c     --------------------------------
c     (lwtc=true) if converged, (.not.lwtc) if not, and deforms grid.
c     ---------------------------------------------------------------
      if(kwt.eq.0) goto 265
      call deform(maxnn,nn,nx,ny,nz,z,u2,nwtl,lwtc,ccw,datum,
     + map,ezl,nlz,ngz,mxgz,zlim,nzt,nnt,in2,net,nex,ney,maxne,lunsat)
c
c
      if((.not.lwtc).and.(jit.lt.maxit2)) goto 501
      if(lwtc) write(6,261) jit
      if(jit.ge.maxit2) write(6,262)
  261 format(10x,'watertable has converged in ',i5,' iterations')
  262 format(10x,'watertable iteration limit reached')
      call flush(6)

c
c     Colomac hardwire: for unsat'd zone option: 
c     adjust surface topography to be 2m above watertable
c     -----------------------------------------------------
c     if(lunsat) then
c     uthick = 2.
c     do i=1,nx
c     do j=1,ny
c     nodewt = (i-1)*nyzt + (j-1)*nzt + nz
c     do k=1,nlz(ngz)
c     dz = uthick/nlz(ngz)
c     node = nodewt + k
c     z(node) = z(nodewt) + k*dz
cc     if(i.eq.1.or.i.eq.nx) write(6,4756) i,j,k,nodewt,node,z(node)
c      enddo
c      enddo
c      enddo
c4756  format('Colomac wt: i,j,k,nodewt,node,z: ',5i6,f10.2)
c      endif
c
  265 continue
c
c      write(7,263)
c  263 format(/10x,'node #, surface z co-ordinates:',/)
c      write(7,264) (i,z(i),i=nzt,nn,nzt)
c  264 format(4(i5,e12.4))
c

      call flush(6)
c
c ================================================================
c
c ===================================================================
c ===================================================================
c
c     compute elemental velocities:
c     ==============================
c
c     call clock@(start)
c
c     numerically evaluated derivatives ...
c     --------------------------------------
      if(kintv.eq.1) 
     +call veloc1(maxnn,maxne,ne,nexy,x,y,z,in,u2,t2,c2,d2,vx,vy,vz,por,
     +  cx,cy,cz,gamma,map,pp,qq,modelwu,lheat,lmass,lage,ts,omega,
     +  modelkr,lkr1top,
     +  por0,alpha0,alphap,rhob,ps,spn0,spn1,spe1)
c
c     direct evaluation at centroid ...
c     flow ne, incidences passed - flow element velocities only
c     ---------------------------------------------------------
      if(kintv.eq.0)
     + call veloc2(maxnn,maxne,ne,nexy,x,y,z,in,u2,t2,c2,d2,vx,vy,vz,
     +   por,cx,cy,cz,dwx,dwy,dwz,exl,eyl,ezl,gamma,map,
     +   pp,qq,modelwu,lheat,lmass,lage,ts,omega,modelkr,lkr1top,
     +   por0,alpha0,alphap,rhob,ps,spn0,spn1,spe1)
c
c     get velocities for 1D line or 2D plane fracture elements
c     ----------------------------------------------------------
      if(nfrac.gt.0) then 
c      write(6,8888)
c8888 format(/10x,'calling vfracture ...')
      call vfracture(maxnn,maxne,ne,in,ezl,exl,eyl,u2,t2,c2,d2,vlin,ckl,
     +      map,xarea,ifracl,lvert,nfrac,nex,ney,nez,maxfrac,ifdim,
     +      vx2d,vy2d,gamma,lheat,lmass,lage,
     +      p2,q2,ts2,omega2,modelwu2,modelkr2,
     +      por0,alpha0,alphap,rhob,ps,spn0,spn1)
      call flush(6)

      endif
c
c     interpolate across capillary fringe only
c     hardwire: first layer above watertable = v at wt
c     could also interpolate to zero at ground surface:
c     ex: vx(l) = vx(lwt) - float((k-nez)/(nezt-nez))*vx(lwt)
c     --------------------------------------------------------
      if(ksat.eq.0) then
c new:
      do k=nez+1,nezt
      do lxy = 1,nexy
      l = (k-1)*nexy + lxy
c     lwt = (nez-1)*nexy + lxy      !element et wt
c     if(k.eq.nez+1) then
c     vx(l) = vx(lwt) 
c     vy(l) = vy(lwt)
c     vz(l) = vz(lwt)
c     else
      vx(l)=1.0e-20
      vy(l)=1.0e-20
      vz(l)=1.0e-20
c     endif
      enddo
      enddo

      endif
c     ---------------
c
c     call clock@(finish)
c     vcpu=vcpu+finish-start

c
c     hardwire !
      do l=1,ne
      if(dabs(vx(l)).lt.1.e-12) vx(l)=1.e-12
      if(dabs(vy(l)).lt.1.e-12) vy(l)=1.e-12
      if(dabs(vz(l)).lt.1.e-12) vz(l)=1.e-12
      enddo

c
c     Reverse velocity fields if needed:
c     -----------------------------------

      if(ivel.eq.0) then

      do i=1,ne
      vx(i) = vx(i)*float(iflipv)
      vy(i) = vy(i)*float(iflipv)
      vz(i) = vz(i)*float(iflipv)
      enddo

      do i=1,nfrac
      if(ifdim(i).eq.2) then
      vx2d(i) = vx2d(i)*float(iflipv)
      vy2d(i) = vy2d(i)*float(iflipv)
      endif
      if(ifdim(i).eq.1) vlin(i) = vlin(i)*float(iflipv)
      enddo

      write(6,234) tdays,vx(1),vx(ne),vz(1),vz(ne)
  234 format(10x,'velocity solution complete at time= ',e12.3,' days',
     +  /10x,'vx(1),vx(ne) = ',2e13.4,' vz(1),vz(ne) = ',2e14.4)

      if(nfrac.gt.0) then
      if(ifdim(1).eq.2) write(6,446) tdays,vx2d(1)
  446 format(10x,'fracture velocities complete t= ',f10.2,' days',
     +            ', vx2d(1): ',e15.5)
      if(ifdim(1).eq.1) write(6,437) tdays,vlin(1) 
  437 format(10x,'fracture velocities complete t= ',f10.2,' days',
     +            ', vlin(1): ',e15.5)
      endif
      endif    !end ivel=0

      call flush(6)
c
c ===================================================================
c     use kfreec ?
c     free up age nodes (assign ic2 = 0) on discharge areas
c     works only for saturated systems
c     by default, all top surface nodes should be defined as A=0 in input file
c     -------------------------------------------------------------------------
      if(kmass.eq.2 .and. ksat.eq.1) then
c     if(ksat.eq.1) then                               !do always, for heat, age or mass
      iflag=0
      do lxy=1,nexy 
      l3d = net -  nexy + lxy
      if(vz(l3d).gt.0.d0) then     !discharge zone Vz >0
        do i=5,8
        node = in2(l3d,i)
        if(ic2(node).gt.0) then
        iflag=iflag+1
        write(6,*) 'found vz>0 ',vz(l3d),' at x= ',x(node)
        ic2(node)=0
        endif
        enddo
      endif
      enddo
      if(iflag.gt.0) then 
              write(6,8550) iflag 
 8550 format(/10x,'converting',i5,
     +            ' t-1 age nodes on exit boundaries ...')
c              
c     re condensation code for age transport 
c     ------------------------------------------------------------
      ktype = 2
      call condense(maxnn,maxne,maxn,maxna,laa,np1,nnt,net,nf,n3,
     + nx,ny,nzt,in2,ic3,fb3,fc3,lc3,fbfc,ib3,a,iaa3,ind3,
     + nw,kprt,ktype,nc3,m3,nb3,x,y,z,ara,datum,leak2,
     + por)
      endif
c
      endif
c
c     skip transport solution if requested:
c     -------------------------------------
c
      if(kcntrl.eq.1 .or. kcntrl.eq.3) goto 302
c
c     short circuit soln - transport only
c     ------------------------------------
  502 continue
c
c
c     check for discharge along top surface and remove fixed T and c
c     later: check only for kmass transport ? 
c     ---------------------------------------------------------------
      if(kfreec.eq.1 .and. kit.lt.3 ) then          !hardwire, to help convergence, do not do type change after 2 iterations
c                                                   !but check ...
      iflag1 = 0
      do i=1,nexy
      do j=1,4
      node = inb(i,j)
      ie3d = ne-nexy + i
      if(vz(ie3d).gt.0.  .and. (ic3(node).eq.1.or.ic4(node).eq.1)) then     !check using vz and if c node is type-1
        ic2(node) = 0                                   !free-up T node
        if(ic3(node).eq.1) ic3(node) = 0                                   !free-up c node
        if(ic4(node).eq.1) ic4(node) = 0                                   !free-up c node
        iflag1 = iflag1+1 
      write(6,8857) i,j,node,ie3d,x(node),y(node),z(node),vz(ie3d)
 8857 format(10x,'type change +vzTOP; i,j,node,ie3d,x,y,z,vz:',
     +                                           4i7,3f12.2,e12.2)
      endif
      enddo
      enddo
c
c     left face: check for -vx    hardwire off
c     do i=1,ney*nez
c     do j=1,4
c     node = inbl(i,j)
c     ie3d = (i-1)*nex+1
c     if(vx(ie3d).lt.0. .and. ic2(node).eq.1) then
c       ic2(node) = 0
c       iflag = 1
c     write(6,8877) ie3d,j,node,ie3d,x(node),y(node),z(node),vx(ie3d)
c8877 format(10x,'type change typ1-0 -vxLEFT; i,j,node,ie3d,x,y,z,vx:',
c    +                                            4i7,3f9.2,e12.2)
c     endif
c     enddo
c     enddo
c
c     right face: check for +vx and free-up if needed
c     only checks ic3 (concentration) but applies to both T and C and D
c     ------------------------------------------------------------------
      iflag2 = 0
      iflag3 = 0
      do i=1,ney*nez
      do j=1,4               
      node = inbr(i,j)
      ie3d = i*nex
      if(vx(ie3d).gt.0.) then                 !velocity outward
        if(ic3(node).eq.1.and.j.ge.3) then    !check top 2 nodes of each element on rt face
        ic2(node) = 0                         !free-up nodes 
        ic3(node) = 0
        ic4(node) = 0                     !free ic4 anyway
        iflag2 = iflag2 + 1
      write(6,8867) i,j,node,ie3d,x(node),y(node),z(node),vx(ie3d)
 8867 format(10x,'type change typ1-0 +vxRIGHT;i,j,node,ie3d,x,y,z,vx:',
     +                                              4i7,3f9.2,e12.2)
      endif

      else                                    !velocity inward
        if(ic3(node).eq.0.and.j.le.2) then    !check bottom 2 nodes
        ic2(node) = 1                         !impose type-1
        ic3(node) = 1
        ic4(node) = 1                     !fix ic4 anyway
        iflag3 = iflag3 + 1
        t0(node) = fc2(node)                 !re-initialize with b.c.
        c0(node) = fc3(node)                 !re-initialize with b.c.
        d0(node) = fc4(node)                 !re-initialize with b.c.
      write(6,8877) i,j,node,ie3d,x(node),y(node),z(node),vx(ie3d)
 8877 format(10x,'type change typ1-0 +vxRIGHT;i,j,node,ie3d,x,y,z,vx:',
     +                                              4i7,3f9.2,e12.2)
      endif

      endif
      enddo
      enddo      
c
c     right face - flip back to type-1 if needed
c     check - applies for T and C
c     are T and C b.c. values still stored ? I think so ....
c     ------------------------------------------------------
c     iflag3 = 0
c     do i=1,ney*nez
c     do j=1,2               !check only bottom 2 nodes
c     node = inbr(i,j)
c     ie3d = i*nex
c     if(vx(ie3d).lt.0. .and. ic3(node).eq.0) then
c       ic2(node) = 1
c       ic3(node) = 1
c       iflag3 = iflag3 + 1
c       t0(node) = fc2(node)                 !re-initialize with b.c.
c       c0(node) = fc3(node)                 !re-initialize with b.c.
c     write(6,8877) i,j,node,ie3d,x(node),y(node),z(node),vx(ie3d)
c8877 format(10x,'type change typ0-1 -vxRIGHT;i,j,node,ie3d,x,y,z,vx:',
c    +                                              4i7,3f9.2,e12.2)
c     endif
c     enddo
c     enddo      
c
c     recondense transport if a fixed c node has changed type on an outflow boundary 
c -------------------------------------------------------------------------------------
c
      if(iflag1.gt.0.or.iflag2.gt.0.or.iflag3.gt.0) then 
              write(6,8551) iflag1,iflag2,iflag3 
 8551 format(/10x,'Converting',3i5,
     +            ' type-1 T or C nodes on exit boundaries ...')
c
c     re condensation code for heat transport 
c     ------------------------------------------------------------
      call condense(maxnn,maxne,maxn,maxna,laa,np1,nnt,net,nf,n2,
     + nx,ny,nzt,in2,ic2,fb2,fc2,lc2,fbft,ib2,a,iaa2,ind2,
     + nw,kprt,ktype,nc2,m2,nb2,x,y,z,ara,datum,leak2,
     + por)
c
c     re condensation code for mass transport 1
c     ------------------------------------------------------------
      call condense(maxnn,maxne,maxn,maxna,laa,np1,nnt,net,nf,n3,
     + nx,ny,nzt,in2,ic3,fb3,fc3,lc3,fbfc,ib3,a,iaa3,ind3,
     + nw,kprt,ktype,nc3,m3,nb3,x,y,z,ara,datum,leak2,
     + por)
c
c     re condensation code for mass transport 2 
c     ------------------------------------------------------------
      if(kmass.eq.4) 
     + call condense(maxnn,maxne,maxn,maxna,laa,np1,nnt,net,nf,n4,
     + nx,ny,nzt,in4,ic4,fb4,fc4,lc4,fbfc,ib4,a,iaa4,ind4,
     + nw,kprt,ktype,nc4,m4,nb4,x,y,z,ara,datum,leak2,
     + por)
c
      endif

      endif                          !end kfreec=1
c
c  ====
c     tsurfmin=999.
c     tsurfmax=-999.
c     fcmin=999.
c     fcmax=-999.
c     u0min=999.
c     u0max=-999.
c     do i=1,nn 
c     tsurfmin = min(tsurfmin,t2(i))
c     tsurfmax = max(tsurfmax,t2(i))           
c     fcmin = min(fcmin,fc2(i))
c     fcmax = max(fcmax,fc2(i))           
c     u0min = min(u0min,t0(i))
c     u0max = max(u0max,t0(i))           
c     enddo
c
c     write(6,3276) tsurfmin,tsurfmax,fcmin,fcmax,u0min,u0max
c3276 format(/10x,'before call trans2 min,max t2,fc,t0: ',6e15.5)

c  ====
c
c     uniform fracture v if ivel=2
c     assumes uniform v is in x-direction
c     ------------------------------------
      if(ivel.eq.2) then
      do i=1,nfrac
      vx2d(i) = vxl*float(iflipv)
      vy2d(i) = vyl*float(iflipv)
      vlin(i) = vxl*float(iflipv)
      enddo
      endif
c
c     ================================================================
c
c     solve for thermal transport ...
c     assembles with {t0}, and returns {t2} as most recent solution
c     {t1} is the intermediate solution at the last iteration
c     thermal parameters passed: tc,tcs,por,cf,cpsm,sat,ara
c     =============================================================
c
c     if(wlh.gt.0.0) then                                              !latent heat version: Co on rhs
      if(lsubtr1.and.lheat) then
      call trans1(maxn,maxnn,maxne,nf,maxna,laa,nw,
     + x,y,z,t0,t1,t2,fc2,fb2,fx,fs,ic,ic2,lc2,vx,vy,vz,in,in2,map,
     + a,aa,iaa2,ind2,ib2,alh,alv,ath,atv,dd,kdisp,n2,nbb,dt1,dt2,nzt,
     + net,nnt,wp,wa,wp1,wa1,pq,tq,cq,por,sw,ag,hag,ttemp,ara,leak2,
     + hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,htmt,kint,tast,tpcgt,tdays,
     +exl,eyl,ezl,wlh,pp,qq,modelwu,vt,tcs,bz,bzf,cpsm,cpi,cf,pi,sat,ne,
     + fbff,rinc,htmf,ifracl,lvert,nfrac,xarea,nex,ney,nez,vlin,maxfrac,
     + tclw,tcli,tclm,ifdim,hgz2,hx,hy,hh,ht,wx,wy,vx2d,vy2d,
     + cot0,cot1,tkl,numis,flux,fluxv,ifl,maxfx,rtrans,decay,decay2,
     + nflux,agefx,rg,rg2df,porsurf,lmass,lheat,lage,ts,rtdiff,
     + p2,q2,ts2,omega2,modelwu2,modelkr2,depc,vexp,smax,spn0,spn1,
     + ppn,qqn)

      write(6,237) tdays
  237 format(/10x,'heat transport(v1) solution complete at t= ',
     +                                                    e12.3,' days')
      endif
      if(lsubtr2.and.lheat) then
      call trans2(maxn,maxnn,maxne,nf,maxna,laa,nw,                             !Thermal R version /R
     + x,y,z,t0,t1,t2,fc2,fb2,fx,fs,ic,ic2,lc2,vx,vy,vz,in,in2,map,
     + a,aa,iaa2,ind2,ib2,alh,alv,ath,atv,dd,kdisp,n2,nbb,dt1,dt2,nzt,
     + net,nnt,wp,wa,wp1,wa1,pq,tq,cq,por,sw,ag,hag,ttemp,ara,leak2,
     + hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,htmt,kint,tast,tpcgt,tdays,
     +exl,eyl,ezl,wlh,pp,qq,modelwu,vt,tcs,bz,bzf,cpsm,cpi,cf,pi,sat,ne,
     + fbff,rinc,htmf,ifracl,lvert,nfrac,xarea,nex,ney,nez,vlin,maxfrac,
     + tclw,tcli,tclm,ifdim,hgz2,hx,hy,hh,ht,wx,wy,vx2d,vy2d,
     + cot0,cot1,tkl,numis,flux,fluxv,ifl,maxfx,rtrans,decay,decay2,
     + nflux,agefx,rg,rg2df,porsurf,lmass,lheat,lage,ts,rtdiff,
     + p2,q2,ts2,omega2,modelwu2,modelkr2,depc,vexp,smax,spn0,spn1,
     + ppn,qqn)

      write(6,238) tdays
  238 format(/10x,'heat transport(v2) solution complete at t= ',
     +                                                    e12.3,' days')
      endif
c
c     tsurfmin=999.
c     tsurfmax=-999.
c     do i=1,nn
c     tsurfmin = min(tsurfmin,t2(i))
c     tsurfmax = max(tsurfmax,t2(i))           
c     enddo
c
c     write(6,3376) tsurfmin,tsurfmax
c3376 format(/10x,'after trans2 min,max temps : ',2e15.5)
c      
c     optional free/thaw deformation
c     will distribute total net dz over top few nodes
c     and updates z and elz dimensions along each nodal column
c     --------------------------------------------------------
      if(kdz.gt.0) 
     + call settlement(maxnn,maxne,nf,nnt,nxy,nzt,nexy,nezt,net,mxgz,
     + nex,ney,t0,t2,z,dz,dztot,ezl,sumdz,por,sw,in2,p,q,modelwu,ts,
     + link,nlz,ngz,dzmin,dzmax,dzl,pi,lheat,lmass,lage,lunsat)
c
c
c     Mass Transport (component 1 'c0,c1,c2') or age transport
c     fbff ok (flow)  fbfc not used
c     ##############################################

      if((lmass.and.(kmass.eq.1.or.kmass.eq.3.or.kmass.eq.4))
     +                                                .or.lage) then 
      call trans3(maxn,maxnn,maxne,nf,maxna,laa,nw,                        !mass transport (was Thermal R version /R
     + x,y,z,c0,c1,c2,fc3,fb3,fx,fs,ic,ic3,lc3,vx,vy,vz,in,in2,map,
     + a,aa,iaa3,ind3,ib3,alh,alv,ath,atv,dd,kdisp,n3,nbb,dt1,dt2,nzt,
     + net,nnt,wp,wa,wp1,wa1,pq,tq,cq,por,sw,ag,hag,ttemp,ara,leak2,
     + hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,htmt,kint,tast,tpcgt,tdays,
     +exl,eyl,ezl,wlh,pp,qq,modelwu,vt,tcs,bz,bzf,cpsm,cpi,cf,pi,sat,ne,
     + fbff,rinc,htmf,ifracl,lvert,nfrac,xarea,nex,ney,nez,vlin,maxfrac,
     + tclw,tcli,tclm,ifdim,hgz2,hx,hy,hh,ht,wx,wy,vx2d,vy2d,
     + cot0,cot1,tkl,numis,flux,fluxv,ifl,maxfx,rtrans,decay,decay2,     !decay is linear decay for c1, decay2 is decay for heat source (or Dirichlet c1)
     + nflux,agefx,rg,rg2df,porsurf,lmass,lheat,lage,ts,rtdiff,
     + p2,q2,ts2,omega2,modelwu2,modelkr2,depc,depd,vexp,smax,spn0,spn1,
     + spe0,spe1,sdepd,ppn,qqn,rhob,t0,t1,icnt)

      write(6,239) tdays
  239 format(10x,'mass transport(v2 c1) solution complete at t= ',
     +                                                    e12.3,' days')

c  ----------------------------------------------------------------------------
c
c     Mass Transport (component 2 'd0,d1,d2') 
c     fbff ok (flow)  fbfc not used
c     ##############################################

      if(kmass.eq.4) then    
      call trans3(maxn,maxnn,maxne,nf,maxna,laa,nw,                        !mass transport (was Thermal R version /R
     + x,y,z,d0,d1,d2,fc4,fb4,fx,fs,ic,ic4,lc4,vx,vy,vz,in,in2,map,
     + a,aa,iaa4,ind4,ib4,alh,alv,ath,atv,dd,kdisp,n4,nbb,dt1,dt2,nzt,
     + net,nnt,wp,wa,wp1,wa1,pq,tq,cq,por,sw,ag,hag,ttemp,ara,leak2,
     + hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,htmt,kint,tast,tpcgt,tdays,
     +exl,eyl,ezl,wlh,pp,qq,modelwu,vt,tcs,bz,bzf,cpsm,cpi,cf,pi,sat,ne,
     + fbff,rinc,htmf,ifracl,lvert,nfrac,xarea,nex,ney,nez,vlin,maxfrac,
     + tclw,tcli,tclm,ifdim,hgz2,hx,hy,hh,ht,wx,wy,vx2d,vy2d,
     + cot0,cot1,tkl,numis,flux,fluxv,ifl,maxfx,rtransc2,decayc2,decay2,   !using rtransc2 and decayc2 for c2
     + nflux,agefx,rg,rg2df,porsurf,lmass,lheat,lage,ts,rtdiff,
     + p2,q2,ts2,omega2,modelwu2,modelkr2,depc,depd,vexp,smax,spn0,spn1,
     + spe0,spe1,sdepd,ppn,qqn,rhob,t0,t1,icnt)

      write(6,2399) tdays
 2399 format(10x,'mass transport(v2 c2) solution complete at t= ',
     +                                                    e12.3,' days')
      endif
c
c
c     write(6,3377) tsurfmin,tsurfmax
c3377 format(/10x,'after trans3 min,max temps : ',2e15.5)
c
c     adjust all given fixed heads to 
c     equivalent freshwater heads for the density case:
c     H = Ho + gamma * c (Ztop - z)   see the User Guide.
c     (where Ho is original head)
c     use gamma=0 for age transport, to turn this off...
c     -----------------------------------------------------------------
c     only update if kfreec=1
c     -----------------------------------------------------------------
      if(gamma.gt.0.d0.and.kcntrl.eq.0.and.kfreec.eq.1) then
       do 7772 i=1,nx
       do 7772 j=1,ny
       nodetop = (i-1)*nyz + j*nz
       do 7772 k=1,nz
       node=(i-1)*nyz + (j-1)*nz + k
       if(ic(node).eq.1) then
        depth = z(nodetop)-z(node)
        beta = gamma*c2(node)*depth     !only component 1 (c2) for density, not d2
        fc(node) = fc0(node) + beta
        u0(node) = fc(node)
        u1(node) = u0(node)
        u2(node) = u0(node)
       endif
 7772  continue
       endif

       endif                           !for transport
c     -------------------------------------------
c     calculate concentration of deposited particles 
c     spn: accumulate deposition of particles by node
c     later: if some get remobilized, need to keep track of this.
c     check - what happens if frozen ? what stops particles from being transferred (retained)? low kr hence low v ?
c --------------------------------------------------------------------------------------------------------

      do l=1,net
      tavg=0.d0
      cavg=0.d0
c     savg=0.d0
      do i=1,8
      node = in(l,i)
      tavg=tavg+(t2(node)+t0(node))/2.
      cavg=cavg+(c2(node)+c0(node))/2.
c     savg=savg+(spn1(node)+spn0(node))/2.
      enddo
      tavg=tavg/8.
      cavg=cavg/8.
c     savg=savg/8.
c     savg = (spe1(l)+spe0(l))/2.d0

c     tt2=0.d00
c     if(depd(l).gt.0.0d0)  tt2 = (depd(l)*rhob(l)*spe1(l)/por(l))
c    +                   * (sqrt(vx(l)**2+vy(l)**2+vz(l)**2))
c    +                   * (1.787d-3*rvisc(tavg,lheat,lmass,lage))        !multiply by viscosity
c    +                   * vt(l)                                          !  * volume
      spe1(l) = spe0(l) + (decay(l)*cavg*por(l)*sw(l)*dt/rhob(l))       !this is good
     +                  - ((sdepd(l)*por(l)*dt)/rhob(l))                !this is good
      enddo
c
c
c     do i=1,nn
c     if(spn1(i).gt.smax) spn1(i)=smax
c     enddo
      do i=1,net
      if(spe1(i).gt.smax(i)) spe1(i)=smax(i)
      enddo
c
c     update porosity from spn
c     check order - update spn1 first or por ?
c     -----------------------------------------
      do l=1,net
      por(l) = por0(l) - spe1(l)*rhob(l)/(phi*ps(l))

      if(por(l).le.0.) then
      write(6,4246) l,spe1(l),rhob(l),phi,ps(l)
 4246 format(/10x,'error, por<0',
     +        /10x,'l,spnavg,rhob,phi,ps... ',i6,4e12.4,
     +        /10x,'program stopping ......')
      stop
      endif

      enddo
c
c     check for convergence:
c     apply relaxation factors to heads and temperature:
c     update head and temperature array:
c     (to save solution at most recent iteration)
c     ---------------------------------------------------
c     ---------------------
      dumx=-999.
      dtmx=-999.
      dcmx=-999.
      ddmx=-999.
      idnode=-999
      idnodec=-999
      idnoded=-999
c
c     flow
c     ----
      do 300 i=1,nn
      adu = abs(u2(i)-u1(i))
      if(adu.gt.dumx) then
             dumx=adu
             dumx2=u2(i)-u1(i)
      endif
c     u2(i)= wfu*u2(i) + (1.-wfu)*u1(i)            !hardwire off
      u1(i)=u2(i)
  300 continue
c
c     heat transport
c     ---------------
      do 301 i=1,nnt
      adt = abs(t2(i)-t1(i))
      if(adt.gt.dtmx) then
            dtmx=adt
            dtmx2=t2(i)-t1(i)
            idnode=i
      endif
c     t2(i)= wft*t2(i) + (1.-wft)*t1(i)             !hardwire off
      t1(i)=t2(i)
  301 continue
c
c     mass transport
c     ---------------
      if(lmass) then
      do 3302 i=1,nnt
      adt = abs(c2(i)-c1(i))
      if(adt.gt.dcmx) then
            dcmx=adt
            dcmx2=c2(i)-c1(i)
            idnodec=i
      endif
c     c2(i)= wft*c2(i) + (1.-wft)*c1(i)             !hardwire off
      c1(i)=c2(i)
c
      if(kmass.eq.4) then               !component 2
      adt = abs(d2(i)-d1(i))
      if(adt.gt.ddmx) then
            ddmx=adt
            ddmx2=d2(i)-d1(i)
            idnoded=i
      endif
      d1(i)=d2(i)             !comp2
      endif
 3302 continue
      endif
c
      write(6,153) dumx2,dtmx2,dcmx2,ddmx2,
     +                         idnode,t2(idnode),z(idnode),
     +                        idnodec,c2(idnodec),z(idnodec),
     +                        idnoded,d2(idnoded),z(idnoded),
     +  dzmin,dzmax
 153  format(10x,'dumx= ',e12.4,4x,'dtmx= ',e12.4,5x,'dcmx(c1)= ',e12.4,
     +                                            5x,'dcmx(c2)= ',e12.4,
     + /10x,'node of max. temp. change: ',i7,' T = ',e13.5,' z=',f7.2,
     + /10x,'node of max. c1 change (c1): ',i7,' c = ',e13.5,' z=',f7.2,
     + /10x,'node of max. c2 change (c2): ',i7,' c = ',e13.5,' z=',f7.2,
     + /10x,'dzmin/dzmax deformation (m) (+dzmax = thaw) : ',2e10.2)

c     if(numss.eq.0) goto 157
c      write(6,154)
c 154  format(10x,'source/sink temps:')
c      write(6,155) (kssn(i),t2(kssn(i)),i=1,numss)
c 155  format(2(10x,i7,2x,f10.2))
c157  continue
c
c     iterate if non-linear solution has not converged
c     and iteration limit has not been exceeded:
c     ##################################################
c
c      if(kcntrl.eq.0 .and. ivel.eq.0) then
      if( (kcntrl .eq.0 .and. ivel.eq.0 .and. kmass.ne.2) .or.               !maxit doesn't apply for age 
     + ((wlh.gt.0. .or. lneg) .and. (kmass.eq.0.or.kmass.ge.3))) then        !only proceed if need to iterate ...

      if(((dumx.gt.ccp).or.(dtmx.gt.cct).or.(dcmx.gt.ccc))
     +                                   .and.(kit.lt.maxit1)) goto 500      !maxit not yet reached, keep iterating ...
      if(kit.ge.maxit1) write(6,162)
      if(kit.lt.maxit1) write(6,163) kit
c      if(kit.ge.maxit1) write(*,162)
c      if(kit.lt.maxit1) write(*,163) kit
 162  format(10x,'thermal/concentration iteration limit reached',
     +          ' at this time step')
 163  format(10x,'thermal solution has converged in',i5,' iterations')
      endif

c     after convergence, accumulate total thaw settlement dztot
c     +dz is thaw settlement, so -ve dztot is settlement
c     ----------------------------------------------------------
      if(kdz.gt.0) then
      dztotmax = -99999.
      dztotmin = +99999.
      do i=1,nxy
      dztot(i) = dztot(i) - dz(i)                  !cumulative dz
      if(dztot(i).gt.dztotmax) dztotmax = dztot(i) 
      if(dztot(i).lt.dztotmin) dztotmin = dztot(i) 
      enddo
      write(6,*) 'cumulative max dztot (+ve = heave     ) : ',dztotmax
      write(6,*) 'cumulative min dztot (-ve = settlement) : ',dztotmin
      endif

c
c     print Pe and Co values in xz plane - at lplot times at knoy(1)
c     to simplify: uses max retardation from thermal, c1 and c2
c     --------------------------------------------------------------
      if(lmom.or.itdt.eq.1) then 
      if(itdt.eq.1) write(59,9850) nezt,nex,tdays
 9850 format('zone i= ',i5,', j= ',i5,',T="',e15.6,' days"')
      if(itdt.gt.1) write(59,9851) nezt,nex,tdays
 9851 format('zone i= ',i5,', j= ',i5,', D=(1,2),T="',e15.6,' days"')
      rmax=max(ret,rtrans,rtransc2,1.d0)
      kny = knoy(1)
      if(knoy(1).lt.1) kny = 1      
      do i=1,nex
      do k=1,nezt
      iel = (k-1)*nex*ney+(kny-1)*nex + i
      xbar = (x(in2(iel,1))+x(in2(iel,2)))/2.d0
      zbar = (z(in2(iel,1))+z(in2(iel,5)))/2.d0
      pxx = (vx(iel)*exl(iel)/rmax)/((vx(iel)*al/rmax)+td)
      pzz = (vz(iel)*ezl(iel)/rmax)/((vz(iel)*al/rmax)+td)
      cxx= (vx(iel)/rmax)*dt/exl(iel)
      czz= (vz(iel)/rmax)*dt/ezl(iel)
      if(itdt.eq.1) write(59,9852) xbar,zbar,pxx,pzz,cxx,czz
      if(itdt.gt.1) write(59,9852) pxx,pzz,cxx,czz
 9852 format(6e15.7)      
      end do
      end do
      call flush(59)
      endif
c

c     -------------------------------------------      
c
c     skip to here for kcntrl=1 flow only
c     ------------------------------------
 302  continue
c
c ===================================================================
c     if(kcntrl.eq.1) then
c      write(12,235)
c  235 format(/10x,'elemental velocites:'/)
c     write(12,236) (vx(iel),vy(iel),vz(iel),iel=1,ne)
c 236 format(2(3e13.6))
c      write(12,337) 
c  337 format(/10x,'z- coordinates:'/)
c
c     write(13,338) (z(map(i)),i=1,nn)
c     write(13,338) (x(map(i)),y(map(i)),i=nz,nn,nz)
c 338 format(8f10.4)
c     endif
c     endif
c
c ===================================================================
c print velocity vectors if required:
c if kprt=1 then at each kplot, but also at each 3dprint time
c ===================================================================
      tdif1=abs(tdays-pt(1))
      tdif2=abs(tdays-pt(2))
      tdif3=abs(tdays-pt(3))
      tdif4=abs(tdays-pt(4))
      tdif5=abs(tdays-pt(5))
      l3dvplt = ((tdif1.lt.1.0e-5).or.
     +           (tdif2.lt.1.0e-5).or.
     +           (tdif3.lt.1.0e-5).or.
     +           (tdif4.lt.1.0e-5).or.
     +           (tdif5.lt.1.0e-5))
c     if((kprt.eq.1.and.lplot) .or. l3dvplt) then
      if(lplot .or. l3dvplt) then

      if(lplot.and.igo.eq.0) then 
      write(6,3990)
 3990 format(/10x,'writing veloc2d files ...')
      write(12,5820) 
 5820 format('title="2d xz vectors"',
     +       /,'variables="x","z","vx","vz"')
      write(34,5823) 
 5823 format('title="2d xy vectors"',
     +       /,'variables="x","y","vx","vy"')
      write(36,5221) 
 5221 format('title="2d yz vectors"',
     +       /,'variables="y","z","vy","vz"')
      igo=1
      endif
      if(lplot.and.igo.eq.1) then
      write(12,5821) nezt,nex,tdays
 5821 format('zone i= ',i5,', j= ',i5,', T="',e12.5,' days"')
      write(34,5821) ney,nex,tdays
      write(36,5821) nezt,ney,tdays
      endif
      
      if(l3dvplt.and.igo3dv.eq.0) then
      write(32,5920) 
 5920 format('title="3d vectors"',
     +       /,'variables="x","y","z","vx","vy","vz"')
      igo3dv=1      
      endif
      if(l3dvplt.and.igo3dv.eq.1) then
      write(32,5921) nezt,ney,nex,tdays
 5921 format('zone i= ',i5,', j= ',i5,', k= ',i5,', T="',e12.5,'"')
      endif
c
c   prints only the first knoy xz-plane of velocities
c   ---------------------------------------------------
      if(lplot) then
      kny = knoy(1)
      if(knoy(1).lt.1) kny = 1
      do 8891 i=1,nex
      do 8891 k=1,nezt
      iel = (k-1)*nex*ney+(kny-1)*nex + i
      xbar = (x(in2(iel,1))+x(in2(iel,2)))/2.d0
      zbar = (z(in2(iel,1))+z(in2(iel,5)))/2.d0
      write(12,236) xbar,zbar,vx(iel),vz(iel)
  236 format(6e14.6)
 8891 continue
c
c   prints only the first knox yz-plane of velocities
c   ---------------------------------------------------
      if(knox(1).gt.0) then
      iel_yl = (nezt-1)*nex*ney+ knox(1) + (ney-1)*nex
      yll=y(in2(iel_yl,4))
      do 8991 i=1,ney
      do 8991 j=1,nezt
      iel = (j-1)*nex*ney+ knox(1) + (i-1)*nex
      ybar = (y(in2(iel,1))+y(in2(iel,4)))/2.d0
      zbar = (z(in2(iel,1))+z(in2(iel,5)))/2.d0
      write(36,236) yll-ybar,zbar,vy(iel),vz(iel)
 8991 continue
      endif

c  xy at knoz(1)
      iez = max(min(knoz(1),nezt),1)
c      iez = max(iez,1)
      do 8992 i=1,nex
      do 8992 j=1,ney
      iel = (iez-1)*nex*ney+(j-1)*nex+ i 
      xbar = (x(in2(iel,1))+x(in2(iel,2)))/2.d0
      ybar = (y(in2(iel,1))+y(in2(iel,4)))/2.d0
      write(34,236) xbar,ybar,vx(iel),vy(iel)
 8992 continue

c
c     fracture velocities in xy plane for krfrac=1
c     for a full xy plane of fractures ( krfrac=1):
c
c     for discrete random fractures in xy plane: (but cannot do particle tracks)
c     write(52,1677) tdays,nfknoz1  then:   /,'zone i= ',i10) 
c     -------------------------------------------------
      if(krfrac.eq.1) then
      if(igof.eq.0) write(52,1477) tdays,ney,nex,tdays/365.
 1477 format(' title="xy fracture velocities on knoz1 element surface, '
     +         'tdays= ',f10.3,'"',
     +       /,'variables="x","y","vx","vy"',
     +       /,'zone i= ',i10,', j= ',i10,', T="',f12.4,' yrs"')
      if(igof.gt.0) write(52,1478) ney,nex,tdays/365.
 1478 format('zone i= ',i10,', j= ',i10,', T="',f12.4,' yrs"')

      vzero=0.
      do 8725 kf=1,nfrac
      l = ifracl(kf)
      call frac_plane(maxfrac,maxne,lvert,inline,in,kf,l,fdimx,fdimy,
     +             exl,eyl,ezl,rexel,sarea)
      xbar= 0.
      ybar= 0.
c     zbar= 0.
      do 8724 j=1,4
      xbar=xbar+x(map(inline(j)))
      ybar=ybar+y(map(inline(j)))
c     zbar=zbar+z(map(inline(j)))
 8724 continue
      xbar=xbar/4.
      ybar=ybar/4.
c     zbar=zbar/4.

      do kk=1,nfknoz1
      if(kf.eq.ifrac2dz1(kk)) write(52,8623) xbar,ybar,vx2d(kf),vy2d(kf)    !2D horizontal fractures shown in xy plane at knoz1
      enddo

 8725 continue
      endif                    !end krfrac=1
      endif                    !end lplot
      igof=1
      call flush(52)
c
c
      endif                         !end lplot and l3dv
c
c     3d velocity output
c     -------------------
      if(l3dvplt) then 
      do i=1,nex
      do j=1,ney
      do k=1,nezt
      iel = (k-1)*nex*ney+(j-1)*nex + i 
      xbar = (x(in2(iel,1))+x(in2(iel,2)))/2.
      ybar = (y(in2(iel,1))+y(in2(iel,4)))/2.
      zbar = (z(in2(iel,1))+z(in2(iel,5)))/2.
      write(32,236) xbar,ybar,zbar,vx(iel),vy(iel),vz(iel)
      enddo
      enddo
      enddo
      endif

      call flush(6)
      call flush(12)
      call flush(32)
      call flush(34)
      call flush(36)
c
c     update {u0}, {t0} and {c0} and spn vectors for next time step:
c     =============================================================== 
      do 303 i=1,nnt
      if(i.le.nn) u0(i)=u2(i)
      t0(i)=t2(i)
      c0(i)=c2(i)
      d0(i)=d2(i)
      spn0(i)=spn1(i)
      spn1(i) = 0.0
 303  continue

      do l=1,net
      spe0(l)=spe1(l)
      enddo     
c
c     back to nodes
c     accumulate element spe to nodes
c     --------------------------------
      do l=1,net
      do i=1,8
      node=in(l,i)
      spn1(node) = spn1(node)+ spe1(l) 
      enddo
      enddo

      do i=1,nnt
      spn1(i) = spn1(i)/float(icnt(i))
      enddo
c
c     heat gain/loss at source/sinks:
c     -------------------------------
      hina=0.
      houta=0.
      hinr=0.
      houtr=0.
      hin2=0.
      hout2=0.

      if(lss) then
      do 145 i=1,nnt
      if(pq(i).gt.0.) then
           hina = hina+pq(i)*(tq(i)-tref)
           hinr = hinr+pq(i)*(tq(i)-t2(mpa(i)))    !injection_T - aquifer_T
           hin2= hin2+pq(i)
      else
           houta = houta+pq(i)*(t2(i)-tref)
           houtr = houtr+pq(i)*(t2(i)-t2(mpa(i)))
           hout2= hout2+pq(i)
      endif
 145  continue
      endif

c
c     skip for flow only:
c     --------------------
      if(kcntrl.eq.1 .or. kcntrl.eq.3) goto 304

      
      cpfdt = cpf*dt
      if(kmass.eq.1.or.kmass.eq.2) cpfdt = dt
      hina=hina*cpfdt
      houta=houta*cpfdt
      hssa=hssa+hina+houta
c
      hinr=hinr*cpfdt
      houtr=houtr*cpfdt
      hssr=hssr+hinr+houtr
c      write(6,146) hina,houta,hssa,hinr,houtr,hssr
c 146  format(/10x,'absolute heat gain at source/sinks: (this step)',
c     +       /10x,'gain at sources  ..................  ',e12.3,
c     +       /10x,'loss at sinks    ..................  ',e12.3,
c     +       /10x,'net cumulative gain/loss ..........  ',e12.3,
c     +       /10x,'relative heat gain at source/sinks: (this step)',
c     +       /10x,'gain at sources  ..................  ',e12.3,
c     +       /10x,'loss at sinks    ..................  ',e12.3,
c     +       /10x,'net cumulative gain/loss ..........  ',e12.3)
      hic=hic+hinr
      hoc=hoc+houtr
c
c     heat gain/loss through top surface:
c     -----------------------------------
      hts=0.
      if(leak2) then
      do 147 i=1,nxy
 147  hts=hts+ara(i)*(ttemp-tref)/bz(i)
      hts=hts*tcs*dt
c      write(6,148) hts
c 148  format(10x, 'heat transfer across top surface ..  ',e12.3)
      endif
      hts2=hts2+hts
      heat1=htota+(hina+houta+hts)

c
c     heat gain/loss from recharge across top surface:
c     ------------------------------------------------
      do k=1,nxy
      taq = t2(k*nzt)
      tempq = max(0.,ttemp+rtdiff-taq)
      hqtop = hqtop + ara(k)*fbff(k)*tempq*cpfdt
c      write(6,4748) k,taq,ara(k),fbff(k),tempq,cpfdt
c4748  format('hqtop:k,taq,ara(k),fbff(k),tempq,cpfdt:',i4,5e12.4)
      enddo
c
c     update elemental volumes, and calculate heat content of aquifer
c     ----------------------------------------------------------------
      if(update_flow) call volume(x,y,z,maxnn,ag,hag,in2,vt,maxne,net,
     +  vol1,kint,exl,eyl,ezl)
c
c     total heat and ice content:
c     (only every kplot or kmom intervals)
c     -------------------------------------
      if((lplot.or.lmom).and.kmass.eq.0.or.kmass.ge.3) then
      htota=0.
      htotr=0.
      vice=0.
      vwat=0.
      do 361 l=1,net
      tavg=0.
      tavgr=0.
      do 362 j=1,8
      tavg=tavg+ t2(in2(l,j))
c     tavgr=tavgr+ (t2(in2(l,j))-t2(mpa(in2(l,j))))
      tavgr=tavgr+ (t2(in2(l,j))+273.)
  362 continue
      tavg=tavg/8.d0
      tavga=tavg-tref
      tavgr=tavgr/8.d0
      ppl=pp(l)
      qql=qq(l)
      cpf = cf*den(tavg,lheat,lmass,lage)
      ww  = por(l)*sw(l)*wu(tavg,ppl,qql,modelwu,ts)
      wi  = por(l)*sw(l)*(1.-wu(tavg,ppl,qql,modelwu,ts))
      cpe = ww*cpf + wi*cpi + (1.-por(l))*cpsm(l)     !also ok for unsat zone with sw(l) ? or use sat? 

c     if(l.gt.ne) cpe = por(l)*sat*wu(tavg,p,q,modelwu,ts)*cpf 
c    +              +   por(l)*sat*(1.-wu(tavg,p,q,modelwu,ts))*cpi 
c    +              +   (1.d0-por(l))*cpsm(l)    !sat is from surfat (variable surface layer sat)

      htota=htota + tavga * vt(l)*cpe
      htotr=htotr + tavgr * vt(l)*cpe
      vice=vice+vt(l)*wi
      vwat=vwat+vt(l)*ww
  361 continue

c      write(6,364) htotr,htota
c  364 format(10x,'relative heat content of aquifer      .....  ',e12.3,
c     +      /10x,'absolute heat content of aquifer      .....  ',e12.3)
c
c     sum in at bottom boundary ...
c     hardwire for Jianwen test single fracture
c     ------------------------------------------
c     sumin = sumin + vy2d(1)
c
c     heat balance:
c     -------------
      hbal= (htota-heat1)/htota * 100.
      write(6,149) hbal
 149  format(10x, 'heat balance (%) ..................  ',f12.3)
      endif

c
c     mass balance aqueous phase and retained particle phase
c     ------------------------------------------------------
      cmass0=0.0
      smass0=0.0
      do l=1,net 
      cavg=0.0
c     spnavg=0.0
      do i=1,8                      !get avg elemental concentrations
      cavg = cavg+c2(in(l,i))       !
c     spnavg = spnavg+spn1(in(l,i)) !
      enddo
      cavg=cavg/8.
c     savg=spnavg/8.                ! cumulative particles over time
      cmass0 = cmass0 + cavg*vt(l)*por(l)*sw(l)*1.d0   !hardwire *1  assuming Co=1 kg_p/m3_w, later: read Co
c     smass0 = smass0 + spe1(l)*vt(l)*rhob(l)*sw(l)*1.d0    !  
      smass0 = smass0 + spe1(l)*vt(l)*rhob(l)*1.d0          !  
      enddo
      cmass(it)=cmass0
      smass(it)=smass0
c
c     advective fluxes at inflow and outflow for 1D case
c     assumes flow is downward from top source
c     advective mass flux J = q*Co*Ao*dt, assuming vz is negative down
c     dispersive mass flux J = por*D*dc/dz*Ao*dt, assuming vz is negative down
c     remember vz is negative if downward
c     -------------------------------------------------------------------------
      areaxy = exl(1)*eyl(1)         !assuming vertical columns so area top = area bottom
      cavg=c2(nzt)
c     cavg=1.0
      gradc1 = (c2(nzt)-c2(nzt-1))/ezl(nezt)
      gradc2 = (c2(2)-c2(1))/ezl(1)
      flxin  = flxin  - vz(nezt)*por(nezt)*sw(nezt)*cavg*areaxy*dt      !top node
      flxout = flxout - vz(1)*por(1)*sw(1)*c2(1)*areaxy*dt                 !bottom node
      dflxin = dflxin - por(nezt)*alh*vz(nezt)*gradc1*areaxy*dt      !dispersive fluxmon top
      dflxout = dflxout - por(1)*alh*vz(1)*gradc2*areaxy*dt      !dispersive fluxmon top
      fluxin(it) = flxin + dflxin 
      fluxout(it) = flxout + dflxout 
c     write(63,3772) tdays,cmass(it),smass(it)
c3772 format(3e12.4)

c----
c     continue flow only
c     -------------------
 304  continue

c
c     fluid mass balance:
c     -------------------
      if(lss) then
      fminj=hin2*dt
      fmpur=hout2*dt
      write(6,247) fminj,fmpur
  247 format(/10x,'fluid volume injected  this step: ',e12.3,
     +       /10x,'fluid volume withdrawn this step: ',e12.3)
      endif
c
c     calculate first and second moments of 'above background' plume
c     also - peak temp trace
c     *************************************************************** 
c
      if(lmom) then
        call moment(maxnn,maxne,in2,nx,ny,nzt,nyzt,x,z,t2,knoy,net,mpa,
     +            xbar,zbar,tpeak,tmin,xvar,zvar,xzvar,xbpk,zbpk,tsa)
c
c     print trace of first moment ... (file 17)
c     print trace of peak temp    ... (file 18)
c     print trace of thermal loss across the unsat. zone (file 19)
c     --------------------------------------------------------------
c     write(17,186) tdays,xbpk,zbpk,xbar,zbar
c      write(*,*) 'tdays,xbpk,zbpk ... ',tdays,xbpk,zbpk
c      write(*,*) '      xbar,zbar ... ',xbar,zbar
      write(18,186) tdays,tpeak,tmin
  186 format(8e16.8)    
      momt=momt+1 
      endif 
c
      write(19,186) tdays,htotr,htota,hinr,hts2,hqtop,vwat,vice
c     write(23,186) tdays,xvar,zvar,xzvar,hbal
c     write(*,*) 'x,z,xzvar ... ',xvar,zvar,xzvar
c
c     write data to plotting files ...
c     ################################
      if(.not.lplot) goto 6101
      if(kcntrl.eq.1.or.kcntrl.eq.3) goto 305      !flow only
c
      write(6,9020)
 9020 format(10x,'plot time ...')      
c
c     1-d vertical strips (to file #11) ...
c     prints 1d vertical profile @ (x,y)=(0,0) (background)
c     ##################################################################
      write(11,608) nzt,tdays
  608 format('zone i=',i5,', T= " vertical strip ',f12.3,' days "')
      write(11,804) 
     +      (z(k),t2(k),wu(t2(k),ppn(k),qqn(k),modelwu,ts),k=1,nzt)
  804 format(3e15.6)
c
c     2-d x-z section for contouring:
c     'knoy' is the y-node section location of the x-z section
c     ########################################################
c
c     x/z temperature data to file #8 ...
c     x/z C and Spn   data to file #62 ...
c     -----------------------------------
      do 260 jny=1,5
      if(knoy(jny).eq.0) goto 260
      izonexz=izonexz+1
      if(kwt.eq.0.and.kmass.ne.0)
     + write(8,2247)izonexz,(tdays/365.),nzt,nx,(tdays/365.),knoy(jny)
2247  format('text x=75., y=70.,f=helv,cs=frame,hu=point,h=20,zn=',i5,
     +   ',c=black,bx=filled,bxf=white,bxo=white,t="',f12.4,' yrs "',
     +   /,'zone i=',i5,' j=',i5,', d=(1,2)',
     +     ',T="t=',f14.3,' yrs,knoy=',i3,'"')
      if(kwt.eq.0.and.lmass)
     + write(62,2247) izonexz,(tdays/365.),nzt,nx,(tdays/365.),knoy(jny)
      if(kwt.eq.1)
     +   write(8,2248)izonexz,(tdays/365.),nzt,nx,(tdays/365.),knoy(jny)
      if(kwt.eq.1.and.lmass)
     +  write(62,2248)izonexz,(tdays/365.),nzt,nx,(tdays/365.),knoy(jny)
2248  format('text x=75., y=70.,f=helv,cs=frame,hu=point,h=20,zn=',i5,
     +   ',c=black,bx=filled,bxf=white,bxo=white,t="',f12.4,' yrs "',
     +   /,'zone i=',i5,', j=',i5,',T="t=',f14.3,' yrs,knoy=',i3,'"')
 808  format('zone i= ',i5,', j= ',i5,',T="',f14.3,' days; ',i3,'"')
c
c     xz temps 
c     --------
      do 805 i=1,nx
      i1= (i-1)*nyzt + (knoy(jny)-1)*nzt
      if(kwt.eq.0) 
     +      write(8,8091) (t2(i1+k),wu(t2(i1+k),
     +                    ppn(i1+k),qqn(i1+k),modelwu,ts),k=1,nzt)
      if(kwt.eq.1) write(8,8091) (x(i1+k),z(i1+k),t2(i1+k),
     +              wu(t2(i1+k),ppn(i1+k),qqn(i1+k),modelwu,ts),k=1,nzt)

c     xz concs
c     --------
      if(lmass) then
      if(kwt.eq.0) 
     +      write(62,8091) (c2(i1+k),d2(i1+k),spn1(i1+k),k=1,nzt)
      if(kwt.eq.1) 
     +      write(62,8091) (x(i1+k),z(i1+k),c2(i1+k),d2(i1+k),spn1(i1+k)
     +                                                         ,k=1,nzt)
      endif
c      write(25,809) (x(i1+k),z(i1+k),(t2(i1+k)-t2(k)),k=1,nzt)
 809  format(3e15.6)
 8091 format(5e15.6)
  805 continue
  260 continue
c
c     y/z temperature data to file #16 ...
c     (transverse section looking downgradient)
c     -----------------------------------------
      do 961 jnx=1,5
      if(knox(jnx).eq.0) goto 961
      izoneyz=izoneyz+1
      write(16,908) izoneyz,tdays,nzt,ny,tdays,knox(jnx)
 908  format('text x=75., y=70.,f=helv,cs=frame,hu=point,h=20,zn=',i5,
     +   ',c=black,bx=filled,bxf=white,bxo=white,t="',f11.3,' days "',
     +   /,'zone i= ',i5,', j= ',i5,
     +     ',T= "t= ',f9.3,' days; yz section ',i3,'"')
c
c     output flow watertable location ...
c     -----------------------------------
c     loop1=knox(jnx)*nyz
c     loop2=loop1-nyz+nz
c     write(16,804) (yl-y(map(i)),z(map(i)),i=loop1,loop2,-nz)
c
c     yz temps
c     --------
      do 905 i=1,ny
      i1= (knox(jnx)-1)*nyzt + (i-1)*nzt
      write(16,8091) ((y(i1+k)),z(i1+k),t2(i1+k),
     +              wu(t2(i1+k),ppn(i1+k),qqn(i1+k),modelwu,ts),k=1,nzt)
  905 continue
  961 continue
c
c     1d vertical temp and c profiles ...
c     (within xz section @ knoy(1))
c     Note hardwire: heads (u2) are also printed on the same nodes - 
c     so only works with saturated option (flow grid=transport grid)
c     heads are not printed to t1dz file if unsat zone is included (replaced by -999)
c     ===============================================================
      kny=knoy(1)
      if(kny.eq.0) kny=1
      do ii=1,5
      k1d = knox(ii)
      if(k1d.eq.0) goto 829
      nd1=(k1d-1)*nyzt + (kny-1)*nzt + 1
      nd2=nd1+nzt-1
      nd1f=(k1d-1)*nyz + (kny-1)*nz + 1    !flow nodes
      nd2f=nd1f+nz-1
      write(20,897) nzt,k1d,tdays
      do k=nd1,nd2
      wuu = wu(t2(k),ppn(k),qqn(k),modelwu,ts)

      knxx = knox(ii)
      knyy = knoy(1)
      if(knox(ii).eq.nx) knxx = nx-1
      if(knoy(1).eq.ny)  knyy = ny-1
      iel = (k-nd1)*nexy + (knyy-1)*nex + knxx                          !element #
      if(ksat.eq.0) then
         if( (k-nd1+1).le.nz) u22=u2(nd1f + k-nd1)
         if( (k-nd1+1).gt.nz) u22=u2(nd2f)                  !extend head at wt into unsat zone
      write(20,898)
     +      z(k),u22,t2(k),c1(k),c2(k),spn1(k),wuu,wuu*sw(iel),dzl(iel)
      endif
      if(ksat.eq.1) write(20,898) z(k),u2(k),t2(k),c1(k),c2(k),spn1(k),
     +                            wuu,wuu*sw(iel),dzl(iel)
      enddo
 897  format('zone i=',i5,',t="k1dx=',i4,1x,f15.4,' days"')
 898  format(9e15.5)
      enddo
 829  continue
c
c     1d horizontal profiles at knoy(1),knoz(i)
c     -----------------------------------------
c
c     if(lplot) ... 1d horizontal profiles at knoy(1),knoz(i)
c     --------------------------------------------------------
      do 8906 i=1,5
      if(knoz(i).gt.nzt) goto 8906
      if(knoz(i).lt.1) goto 8906
      nd1= (knoy(1)-1)*nzt + knoz(i)
      nd2=nd1+(nx-1)*nyzt
      write(31,8897) nx,tdays,knoz(i)

c     if(i.eq.1) write(61,8385) nx,tdays
c8385 format('zone i=',i6,',T="',f12.3,', days, q(m/s) top"')

      ix=0
      do k=nd1,nd2,nyzt
      ix=ix+1
c     kxy = (k-knoz(1) + nzt)/nzt 
      write(31,8898) x(k),t2(k),wu(t2(k),ppn(k),qqn(k),modelwu,ts),c2(k)
      enddo
8897  format('zone i=',i6,',T="',f12.3,', days, kz=',i5,'"')
8898  format(7e15.5)
8906  continue
c     =========================================
c     ============================= 
c
c     1d horizontal profiles at knoy(1),nzt
c     -----------------------------------------      
      write(61,8385) nx,tdays,nzt
 8385 format('zone i=',i6,',T="t=',f12.3,', days, kz=',i5,'"')

      nd1= knoy(1)*nzt 
      nd2=nd1+(nx-1)*nyzt
      write(60,8887) nx,tdays,nzt
8887  format('zone i=',i6,',T="',f12.3,
     +       ', days, heatflux across transfer layer kz=',i5,'"')      

      ix=0
      do k=nd1,nd2,nyzt     !loop in x along top nodes at knoy(1)
      ix=ix+1
      kxy = k/nzt 
c
c     calculate Darcy flux across top
      l = net-nexy 
     +      + ((min(max((knoy(1)-1),1),ney-1)-1)*nex) + (min(ix,nex))       !element
      tavg=0.0
      do j=1,8
      tavg=tavg+t2(in(l,j))
      enddo
      tavg=tavg/8.
      ppl=pp(l)
      qql=qq(l)
      qz = vz(l)*por(l)*sw(l)*wu(tavg,ppl,qql,modelwu,ts)
c
c     depth to Ts, search up
c     ------------------------
      nbot=(ix-1)*nyzt + ((knoy(1)-1)*nzt)+1
      ntop=nbot+nzt-1
      ztots=0.
      do kz=2,nzt
      node = nbot+kz-1
      if(t2(node).ge.0.0) then
c
c     interpolate between this node and underlying node 
c     -------------------------------------------------
      tdiff1 = t2(node)-t2(node-1)
      if(tdiff1.gt.0.01)  then
             ztots = z(ntop) - 
     +       (z(node-1) + (ts-t2(node-1))*(z(node)-z(node-1))/tdiff1)
      else
             ztots = z(ntop)-(z(node)+z(node-1))/2.
      endif
      goto 5998
      endif
      enddo
c     
5998  continue

      ztotsx(ix) = ztots
c
      write(61,4809) 
     +    x(k),dztot(kxy),ztotsx(ix),bz(kxy),bzf(kxy),bzqn(kxy),qz

c      enddo
c8897  format('zone i=',i6,',T="',f12.3,', days, kz=',i5,'"')
c8898  format(7e15.5)
c8906  continue
c      
c     conductive heat flux across transfer layer J = -tkl()*gradt + (-bzf)
c     - need to add advective flux

      if(leak2) then 
      gradtt = (ttemp-t2(ntop))/bz(kxy)
c     hflux1dx =  -tcs * gradtt + (-bzf(kxy))  ! check sign of bzf; hflux1dx<0 heat flux down; bzf>0 heat flux down
      hflux1dx =  -tcs * gradtt   ! 
      write(60,8898) x(k),hflux1dx
      endif

      enddo


      write(88,8897) nx-1,tdays,knoz(1)         !debug - only for 1dx systems
      if(knoz(1).gt.0.and.knoy(1).gt.0) then 
      do k=1,nex
      node = (k-1)*ny*nzt+(knoy(1)-1)*nzt+knoz(1)
      ielem = (knoz(1)-1)*nex*ney + (knoy(1)-1)*nex + k 
      write(88,809) x(node),cot1(ielem),tkl(ielem)             
      enddo
      endif
      call flush(31)
      call flush(88)
c
c     head data ( to file #9) hxz.out ...
c     ------------------------------------------

 305  continue

      do 360 jny=1,2
      if(knoy(jny).eq.0) goto 360
      write(9,810) nz,nx,tdays,knoy(jny)
 810  format('zone i= ',i5,', j= ',i5,',t="',f12.3,' days; ',i3,'"')
      do 811 i=1,nx
      i1= (i-1)*nyz + (knoy(jny)-1)*nz
      write(9,809) (x(map(i1+k)),z(map(i1+k)),u2(i1+k),k=1,nz)
  811 continue
  360 continue
c
c     y/z head data to file #21 ...
c     (transverse section looking downgradient)
c     -----------------------------------------
      do 9960 jnx=1,5
      if(knox(jnx).eq.0) goto 9960
      write(21,9908) nz,ny,tdays,knox(jnx)
 9908 format('zone i= ',i5,', j= ',i5,
     +       ',T= "t= ',f9.3,' days; yz section ',i3,'"')
c
c     yz head
c     --------
      do 9905 i=1,ny
      i1= (knox(jnx)-1)*nyz + (i-1)*nz
      write(21,809) ((y(map(i1+k))),z(map(i1+k)),u2(i1+k),k=1,nz)
 9905 continue
 9960 continue
c
c     print xy plane solution ...
c     ---------------------------------------------------
c
c     heads....
c     ----------
      do ik=1,5
      if(knoz(ik).lt.1) goto 8335
      izonexy=izonexy+1
      iz = knoz(ik)
      write(46,910) izonexy,tdays,ny,nx,tdays,iz
      do  i=1,nx
      i1=(i-1)*nyz 
      write(46,809) (u2(i1+k),k=iz,nyz,nz)
      enddo

c     temperatures....
c     -----------------
      write(14,910) izonexy,tdays,ny,nx,tdays,knoz(ik)
 910  format('text x=75., y=70.,f=helv,cs=frame,hu=point,h=20,zn=',i5,
     +   ',c=black,bx=filled,bxf=white,bxo=white,t="',f12.4,' days "', 
     +    /,'zone i= ',i5,', j= ',i5,',D=(1,2)',
     +       ',T= "t= ',f12.4,' days; xy section ',i3,'"')
      do 812 i=1,nx
      i1=(i-1)*nyzt 
      write(14,8091) (t2(i1+k),
     +       wu(t2(i1+k),ppn(i1+k),qqn(i1+k),modelwu,ts),c2(i1+k),
     +                                              k=knoz(ik),nyzt,nzt)
 812  continue
      enddo          !for knoz(1-5)

 8335 continue       !for knoz(i)=0
c
c
c     Tpeak and Tmin (normally for probability and life expectancy results)
c    ------------------------------------------------------------------------
      write(41,8422) ny,nx,tdays
8422  format('zone i= ',i5,', j= ',i5,',D=(1,2), T="t= ',f10.1,' days"')
      do  i=1,nx
      do  j=1,ny
      tpeak = -1.e33
      tmin  =  1.e33
      do  k=1,nzt
      i1=(i-1)*nyzt+(j-1)*nzt+k
      t2i1 = t2(i1)
      if(t2i1.gt.tpeak) tpeak = t2i1 
      if(t2i1.lt.tmin)  tmin = t2i1
      enddo 
c     nodexy = (j-1)*nx + i   ????
      nodexy2 = (i-1)*ny+j  
      write(41,8093) tpeak,tmin,dz(nodexy2),dztot(nodexy2)
8093  format(4e15.6)      
      enddo
      enddo
      call flush(41)

 6101 continue              ! end kplot prints      


c     if lmom ... print K's in flow grid in plane knoy(1) (or in plane 1) to file 51
c     ------------------------------------------------------------------
      if(lmom) then
      write(51,7810) nz,nx,tdays 
c7810 format('zone i= ',i5,', j= ',i5,',D=(1,2),T="t= ',f12.4,' days"')
 7810 format('zone i= ',i5,', j= ',i5,',T="t= ',f12.4,' days"')

      do i=1,nx
      j=max(knoy(1),1)
      do k=1,nz

      node = (i-1)*ny*nz + (j-1)*nz +k       !flow node
      l=(k-1)*nex*ney + (j-1)*nex + i        !element associated with node
      if(i.eq.nx.and.k.eq.nz) l=ne-ney+1
      if(i.eq.nx.and.k.lt.nz) l=k*nex*ney-ney+1
      if(i.lt.nx.and.k.eq.nz) l=(k-2)*nex*ney+(j-1)*nex + i  
      cxx= cx(l)
      porel= por(l)
      tklel= tkl(l)
      por0el= por0(l)
c
c get real K from Wu
c     zkrw=1.d0
      swe=wu(t2(map(node)),pp(l),qq(l),modelwu,ts)
c      zkrw = ((swe - p)/(1.d0-p))**4
c      zkrw = max(zkrw,1.e-6)
      alpha = alpha0 + por(l)*alphap*(rhob(l)/ps(l))*spe1(l)
      zkrw = fnzkrw(swe,pp(l),omega,porel,modelkr,por0el,alpha0,alpha)             !relative permeability function
      if(lkr1top .and. l.ge.(ne-nexy+1)) zkrw=1.d0                    !hardwire keep kr=1 for top surface elements

c     if(i.eq.nx/2 .and. k.eq.1) then
c      write(6,9229) 
c    + por(l),por0(l),rhob(l),ps(l),spn1(node),alphap,alpha0,alpha,zkrw
c9229 format(2x,'3:por,por0,rhob,ps,spn1,alphap,alpha0,alpha,zkrw:',
c    +       9e10.3)
c     endif 
c
c     update temp. dependent conductivity:
c     adjust for ice saturation dependent relative k
c     check: rden should give true K's (not needed in flow routine -see user guide equations)
c     ---------------------------------------------
c     cxx = cxx*den(t2(map(node)),lmass,lage)/1000.
      cxx2=(cxx*(den(t2(map(node)),lheat,lmass,lage)/1000.)
     +                 /rvisc(t2(map(node)),lheat,lmass,lage)) * zkrw
c
      write(51,856) x(map(node)),z(map(node)),
c     write(51,856)                                                                ! removed x,z (D=1,2
     +              cxx,cxx2,log10(cxx2),porel,tklel,rhob(l)
 856  format(8e16.7)

c     write(6,4844) node,map(node),den(t2(map(node)),lheat,lmass,lage),
c    +              rvisc(t2(map(node)),lheat,lmass,lage),
c    +              zkrw,t2(map(node))
c4844  format(/10x,'node,map(node),den,rvisc,zkrw,t2:',2i8,4e12.3)

      end do
      end do     
      call flush(51) 
      endif


c
c     heatflow across top and bottom surfaces - by element
c     J = - K* gradT    W/m^2   conduction only
c     (W/m^2 * area = J/s)
c     -----------------------------------------------------
      if(lmom) then
      hflxtopw=0.
      hflxbotw=0.
      do i=1,nex
      do j=1,ney
       lt = (net-nexy)+(j-1)*nex+i
       lb = (j-1)*nex+i
       in2lt1 = in2(lt,1)
       in2lt2 = in2(lt,2)
       in2lt3 = in2(lt,3)
       in2lt4 = in2(lt,4)
       in2lt5 = in2(lt,5)
       in2lt6 = in2(lt,6)
       in2lt7 = in2(lt,7)
       in2lt8 = in2(lt,8)
       in2lb1 = in2(lb,1)
       in2lb2 = in2(lb,2)
       in2lb3 = in2(lb,3)
       in2lb4 = in2(lb,4)
       in2lb5 = in2(lb,5)
       in2lb6 = in2(lb,6)
       in2lb7 = in2(lb,7)
       in2lb8 = in2(lb,8)
       tavgt1 = (t2(in2lt5)+t2(in2lt6)+t2(in2lt7)+t2(in2lt8))/4.   !surface t1 is above t2
       tavgt2 = (t2(in2lt1)+t2(in2lt2)+t2(in2lt3)+t2(in2lt4))/4.
       tavgb1 = (t2(in2lb5)+t2(in2lb6)+t2(in2lb7)+t2(in2lb8))/4.   !surface b1 is above b2
       tavgb2 = (t2(in2lb1)+t2(in2lb2)+t2(in2lb3)+t2(in2lb4))/4.   

       xavg =   (x(in2lt5)+x(in2lt6)+x(in2lt7)+x(in2lt8))/4.
       yavg =   (y(in2lt5)+y(in2lt6)+y(in2lt7)+y(in2lt8))/4.
       zavgt1 = (z(in2lt5)+z(in2lt6)+z(in2lt7)+z(in2lt8))/4.   !surface t1 is above t2
       zavgt2 = (z(in2lt1)+z(in2lt2)+z(in2lt3)+z(in2lt4))/4.
       zavgb1 = (z(in2lb5)+z(in2lb6)+z(in2lb7)+z(in2lb8))/4.   !surface b1 is above b2
       zavgb2 = (z(in2lb1)+z(in2lb2)+z(in2lb3)+z(in2lb4))/4.
      tcont = tkl(lt)
      tconb = tkl(lb)
      atop = exl(lt)*eyl(lt)
      abot = exl(lb)*eyl(lb)
      heatfluxtop = -tcont*((tavgt1-tavgt2)/(zavgt1-zavgt2))     !+flux = heat leaving top J/s/m2=W/m2
      heatfluxbot = -tconb*((tavgb2-tavgb1)/(zavgb1-zavgb2))     !+flux = heat leaving bot J/s/m2=W/m2
      hflxtopw = hflxtopw + heatfluxtop * atop           !surface flux top (W)
      hflxbotw = hflxbotw + heatfluxbot * abot           !surface flux bot (W)
      if(lplot) then
      write(44,8173) ney,nex,tdays
 8173 format('zone i= ',i6,',j= ',i6,', T= "t= ',f7.1,' days"')
      write(44,8172) xavg,yavg,heatfluxtop,heatfluxbot       !flux J/m2/s = W/m2
 8172 format(4e15.6)
      endif
      enddo
      enddo
      call flush(44)
c
c     heatflux across right and left bdy, every kmom time step (W/m^2 * area = J/s) and cumulative (J) (kmom should = 1)
c     conductive flux -ve for leftward conductive flux on right bdy
c     conductive flux +ve for leftward conductive flux on left bdy
c     so assumes vx is +ve to the right, and across both left & rt boundaries, a +ve gradT gives a negative conductive flux
c -------------------------------------------------------------------------------------------------------------------------
c     average heads left and right    
      havgleft=0.
      havgright=0.
      do ileft=1,nyz
      iright = nn-nyz+ileft
      havgleft = havgleft + u2(ileft)
      havgright = havgright + u2(iright)
      enddo
      havgleft = havgleft/nyz
      havgright= havgright/nyz
      hgrad = (havgleft-havgright)/(x(nnt)-x(1))      !use average bdy heads for TH2,TH3 benchmarks ok

      do k=1,2
      heatflux=0.
      qflux = 0.
      arealeft=0.
      arearight=0.

      do lr = nex,net,nex
      if(k.eq.1) l = lr-nex+1        !k=1: left bdy element
      if(k.eq.2) l = lr              !k=2: right bdy element
c     sign = float((k*2)-3)          !k=1 = -    k=2 = +
      areal = eyl(l)*ezl(l)
      if(k.eq.1) arealeft = arealeft + areal
      if(k.eq.2) arearight = arearight + areal
      tavg1 = (t2(in2(l,1))+t2(in2(l,4))+t2(in2(l,5))+t2(in2(l,8))) /4.
      tavg2 = (t2(in2(l,2))+t2(in2(l,3))+t2(in2(l,6))+t2(in2(l,7))) /4.
      tavgl = (tavg1+tavg2)/2.
      xavg1 = (x(in2(l,1))+x(in2(l,4))+x(in2(l,5))+x(in2(l,8)))/4. 
      xavg2 = (x(in2(l,2))+x(in2(l,3))+x(in2(l,6))+x(in2(l,7)))/4.

      ww   = por(l)*wu(tavgl,pp(l),qq(l),modelwu,ts)   *sw(l)                          !elemental *sw(lr) ??? 
      heatflux = heatflux + (-tkl(l)*((tavg2-tavg1)/(xavg2-xavg1))      !conductive flux 
     +  + vx(l)*ww*(tavgl-tref)*cf*den(tavgl,lheat,lmass,lage)) * areal    !convective flux (+ for T>tref for vx+)
      qflux = qflux + vx(l) * areal * ww        !(Q=v*por*A) 
      enddo

      if(k.eq.1) then
        ckequivl = 1.0e-20              
        heatfluxl = heatflux
        darcyfluxl = qflux/arealeft                               !q = Q/A 
        heatfluxlcum = heatfluxlcum + heatflux * dt     
        if(hgrad.gt.1.0e-20) ckequivl = darcyfluxl / hgrad      !equivalent K  
c        write(6,8668) 
c     +  havgleft,havgright,arealeft,arearight,qflux,hgrad,ckequivl
c8668   format( 'k=1: hleft,hright,aleft,aright,qflux,hgrad,ckequivl:',
c     +       /, 7e12.4) 
      endif          
      if(k.eq.2) then
        ckequivr = 1.0e-20
        heatfluxr = heatflux
        darcyfluxr = qflux/arearight
        heatfluxrcum = heatfluxrcum + heatflux * dt
        if(hgrad.gt.1.0e-20)  ckequivr = darcyfluxr / hgrad     !equivalent K 
c       write(6,8778) 
c     +  havgleft,havgright,arealeft,arearight,qflux,hgrad,ckequivr
c8778   format( 'k=2: hleft,hright,aleft,aright,qflux,hgrad,ckequivr:',
c     +       /, 7e12.4) 
      endif

      enddo
c
c     total liquid water volume in domain
c     -----------------------------------
      vwtot=0.
      do l=1,net
        tavg=0.
        do i=1,8
         tavg=tavg+t2(in2(l,i))
        enddo
      tavg=tavg/8.
      vwtot=vwtot+vt(l)*por(l)*wu(tavg,pp(l),qq(l),modelwu,ts)  
      enddo

      write(54,5758) tdays,heatfluxl,heatfluxr,heatfluxlcum,
     + heatfluxrcum,hflxtopw,hflxbotw,vwtot,ckequivl,ckequivr

      endif                !kmom
c
c     1d vertical profiles in time ...
c     (within xz section @ knoy(1) at knox(1))
c     hardwire iskip - every iskip1 time steps
c     =========================================
      ii=1
      iskip1 = 2
      lplot43  = ((itdt/iskip1) * iskip1.eq.itdt)
      k1dx = knox(ii)
      k1dy = knoy(ii)
      if(k1dx.eq.0 .or. k1dy.eq.0 .or. .not.lplot43) goto 8290
      numpt43=numpt43+1
      nd1=(k1dx-1)*nyzt + (k1dy-1)*nzt + 1
      nd2=nd1+nzt-1
      write(43,803)(tdays,z(k),t2(k),wu(t2(k),
     +                            ppn(k),qqn(k),modelwu,ts),k=nd1,nd2)
 803  format(4e15.5)
8290  continue     
c
c     1d horizontal profiles in time ...
c     (within xz section @ knoy(1) at knoz(1))
c     hardwire iskip - every iskip2 time steps
c     =========================================
      ii=1
      iskip2 = 100
      lplot49  = ((itdt/iskip2) * iskip2.eq.itdt)
      k1dz = knoz(ii)
      k1dy = knoy(ii)
      if(k1dz.eq.0 .or. k1dy.eq.0 .or. .not.lplot49) goto 8291
      numpt49=numpt49+1
      nd1= (k1dy-1)*nzt + k1dz
      nd2=nd1 + (nx-1)*nyzt
      write(49,803)(tdays,x(k),t2(k),wu(t2(k),ppn(k),qqn(k),modelwu,ts),
     +              k=nd1,nd2,nyzt)
8291  continue            
c
c     assemble breakthrough data - loop over wells
c     --------------------------------------------
      do 850 iw=1,nwells
      tavw = 0.
      havw = 0.
      tpk  = -999.
      tmin = +999.
      cpkc1  = -999.
      cpkc2  = -999.
      cmin = +999.
      spk = -999.
c
c     loop over vertical screen of this well
c     cavw:  avg. conc in well
c     tpkwell = peak temp.
c     cpkwell = peak conc.
c     --------------------------------------
      havg=0.
      do 851 k=izw1(iw),izw2(iw)
      node =  (ixw(iw)-1)*nyzt + (iyw(iw)-1)*nzt + k
      nodef=  (ixw(iw)-1)*nyz  + (iyw(iw)-1)*nz  + k       !flow node
      nodefb= (ixw(iw)-1)*nyz  + (iyw(iw)-1)*nz  + 1       !flow node at bottom of grid
      if(node.lt.1 .or. node.gt.nnt) then
              write(6,8439) node
 8439 format(/10x,'error in breakthrough curve nodes ixw() etc.',
     +       /10x,'btc node = ',i6,
     +       /10x,'stopping ...')
              stop
      endif
      tavw  = tavw + t2(node)
      if(ksat.eq.1.or.k.le.nz)  havw  = havw + u2(nodef)     !!! note: will only work for ksat=1 (saturated) 
      if(ksat.eq.0.) then
         if(k.le.nz) havg = havg + u2(nodef)
         if(k.gt.nz) havg = havg + u2(nodefb+nz-1)     !extend wt head into unsat zone
      endif
      if(t2(node) .gt. tpk ) tpk = t2(node)
      if(t2(node) .lt. tmin) tmin = t2(node)
      if(c2(node) .gt. cpkc1 ) cpkc1 = c2(node)
      if(d2(node) .gt. cpkc2 ) cpkc2 = d2(node)
      if(c2(node) .lt. cmin) cmin = c2(node)
      if(spn1(node) .gt. spk) spk = spn1(node)
 851  continue
      tavw=tavw/(izw2(iw)-izw1(iw)+1)
      havw=havw/(izw2(iw)-izw1(iw)+1)
      ubt(iw,it)  = tavw
      ubh(iw,it)  = havw
      tpkwell(iw,it) = tpk
      tminwell(iw,it) = tmin
      cpkwell(iw,it,1) = cpkc1
      cpkwell(iw,it,2) = cpkc2
      spnwell(iw,it) = spk
      cminwell(iw,it) = cmin
 850  continue
      btime(it) = tdays
c     -------------------------------------------------------------
c
c     3d temperature and head output
c     any of these files can be used later for a restart - use kgo=1
c     (rename to therm_in.data)
c     ==========================================================
      tdif1=abs(tdays-pt(1))
      tdif2=abs(tdays-pt(2))
      tdif3=abs(tdays-pt(3))
      tdif4=abs(tdays-pt(4))
      tdif5=abs(tdays-pt(5))
      if(tdif1.le.1.0e-7) write (26,131) tdays,nzt,ny,nx
      if(tdif2.le.1.0e-7) write (27,131) tdays,nzt,ny,nx
      if(tdif3.le.1.0e-7) write (28,131) tdays,nzt,ny,nx
      if(tdif4.le.1.0e-7) write (29,131) tdays,nzt,ny,nx
      if(tdif5.le.1.0e-7) write (30,131) tdays,nzt,ny,nx
  131 format (10x,'Title="T,C,H, at t= ',e12.4,' days"',
     +       /10x,'variables="x","y","z","T","C","H"',
     +       /10x,'zone i= ',i7,', j= ',i7,', k= ',i7,', T="T,C,H"')
c
      if(tdif1.le.1.0e-7) 
     + write(26,132) (x(i),y(i),z(i),t2(i),c2(i),u2(i),i=1,nnt)
      if(tdif2.le.1.0e-7) 
     + write(27,132) (x(i),y(i),z(i),t2(i),c2(i),u2(i),i=1,nnt)
      if(tdif3.le.1.0e-7) 
     + write(28,132) (x(i),y(i),z(i),t2(i),c2(i),u2(i),i=1,nnt)
      if(tdif4.le.1.0e-7) 
     + write(29,132) (x(i),y(i),z(i),t2(i),c2(i),u2(i),i=1,nnt)
      if(tdif5.le.1.0e-7) 
     + write(30,132) (x(i),y(i),z(i),t2(i),c2(i),u2(i),i=1,nnt)
  132 format (6e15.6)
c
c     if(tdif1.le.1.0e-7) write (26,133)
c    +            nz,ny,nx,(x(map(i)),y(map(i)),z(map(i)),u2(i),i=1,nn)
c     if(tdif2.le.1.0e-7)  write (27,133)
c    +            nz,ny,nx,(x(map(i)),y(map(i)),z(map(i)),u2(i),i=1,nn)
c     if(tdif3.le.1.0e-7)  write (28,133)
c    +            nz,ny,nx,(x(map(i)),y(map(i)),z(map(i)),u2(i),i=1,nn)
c     if(tdif4.le.1.0e-7)  write (29,133)
c    +            nz,ny,nx,(x(map(i)),y(map(i)),z(map(i)),u2(i),i=1,nn)
c     if(tdif5.le.1.0e-7)  write (30,133)
c    +            nz,ny,nx,(x(map(i)),y(map(i)),z(map(i)),u2(i),i=1,nn)

c 133 format ('zone i=',i7,', j=',i7,', k=',i7,', T="head"',
c    +        /,(4e15.6))

c
c     3d movie file
c     -------------
      if(l3dplt) then
      izone3d = izone3d +1
      if(iplot3d.eq.0) write (35,1031) nzt,ny,nx,tdays
 1031 format (10x,'Title="temperature 3D"',
     +   /10x,'variables="x","y","z","Temperature"',
     +   /10x,'zone i= ',i7,', j= ',i7,', k= ',i7,',T= "',f6.1,' days"')
      if(iplot3d.ne.0) write (35,1131) nzt,ny,nx,tdays
 1131 format (10x,'zone i= ',i7,', j= ',i7,', k= ',i7,
     +            ',D=(1,2,3)',',T= "',f6.1,' days"')
c
      if(iplot3d.eq.0) write (35,132) (x(i),y(i),z(i),t2(i),i=1,nnt)
      if(iplot3d.ne.0) write (35,132) (t2(i),i=1,nnt)
      write (35,1231)  tdays,izone3d
 1231 format('text x=80., y=90.,f=helv,cs=frame,hu=point,',
     + 'h=20,c=black,bx=filled,bxf=white,bxo=white,t="',f11.3,' days "',
     + ',zn=',i5)
      iplot3d=1
      call flush(35)
      endif
c
c     3D mass output for kmass=1 at this time step
c     ----------------------------------------------
c     mass in porous blocks:
c     -------------------------------------------------------
      if((kmass.eq.1.or.kmass.eq.3) .and. kprt.eq.1 .and. lmom) then
      xmass=0.d0
      do l=1,net
      xe=0.d0
      do k=1,8
      xe = xe + t2(in(l,k))
      enddo
      xmass = xmass + xe*vt(l)*por(l)/8.d0
      enddo

c     mass in fractures:
c     -----------------------
      fmass = 0.d0
      do kf=1,nfrac
      l = ifracl(kf)
      if(ifdim(kf).eq.2) then
      call frac_plane(maxfrac,maxne,lvert,inline,in,kf,l,fdimx,fdimy,
     +             exl,eyl,ezl,rexel,sarea)
      fmass = fmass + (t2(map(inline(1)))+t2(map(inline(2)))
     +                +t2(map(inline(3)))+t2(map(inline(4))))
     +                *sarea*xarea(kf)/4.d0
      endif

      if(ifdim(kf).eq.1) then
      call frac_line(maxfrac,maxne,lvert,inline,in,kf,l,fracdim,
     +             exl,eyl,ezl)
      fmass = fmass + (t2(map(inline(1)))+t2(map(inline(2))))
     +                                *fracdim*xarea(kf)/2.d0
      endif
      enddo

c     mass within internal source elements (even if fluxin=0)
c     just used for convenience - use isource to identify elements around injection well
c     ----------------------------------------------------------------------------------
      wmass = 0.d0
      do i=1,numis
      l=ifl(i)
      xe=0.d0
      do k=1,8
      xe = xe + t2(in(l,k))
      enddo
      wmass = wmass + xe*vt(l)*por(l)/8.d0
      enddo
      write(40,2240) tdays,xmass-xmass0,fmass-fmass0,wmass-wmass0
 2240 format(4e15.6)

      endif
c
c     find freezing front in knoy(1) plane at knoz(1)
c     searches from left to right for first T > Ts
c     only for 1D x systems
c     -----------------------------------------------
      if(ny.eq.2 .and. nzt.eq.2 .and. lheat ) then
      xts(it) = 0.
      kny=knoy(1)
      knz=knoz(1)
      if(kny.eq.0) kny=1
      if(knz.eq.0) knz=1
      do i=1,nx
      node2 = (i-1)*ny*nzt+(kny-1)*nzt + knz
      if(i.gt.1 .and. t2(node2).gt.ts) then
         node1 = node2 - ny*nzt
         xts(it) = ((x(node2)-x(node1))/(t2(node2)-t2(node1)))
     +                      *(ts-t2(node1)) + x(node1)
         goto 1445
      endif
      enddo
 1445 continue
      endif
c
c     find thaw and freezing front in knoy(1) plane at knox(1)
c     searches from top to bottom for first T > Ts
c     zts1: z_T>0    zts2:  z_T<0
c     -----------------------------------------------
      if(ny.eq.2 .and. nx.eq.2 .and. lheat ) then
      kny=knoy(1)
      knx=knox(1)
      if(kny.eq.0) kny=1
      if(knx.eq.0) knx=1
c   zts1 ...      
      do i=1,nzt
      node2 = (knx-1)*ny*nzt+(kny-1)*nzt + nzt-i+1
      zts1(it) = z(node2)
      if(i.gt.1 .and. t2(node2).gt.ts) then
         node1 = node2 + 1
         if(dabs(t2(node2)-t2(node1)).gt.1.0e-5)
     +    zts1(it) = z(nzt)- (z(node2)-
     +              ((t2(node2)-ts)/(t2(node2)-t2(node1)) 
     +             * (z(node2)-z(node1))))
c        write(6,5009) node1,node2,ts,z(node1),t2(node1),
c    +                 z(node2),t2(node2),zts1(it)
c5009  format(10x,'freeze front debug: n1,n2,ts,(z1,T1),(z2,T2),zts:',
c    +            2i5,6f10.3)
         goto 1446
      endif
      enddo
 1446 continue    
c   zts2 ...   
      do i=1,nzt
      node2 = (knx-1)*ny*nzt+(kny-1)*nzt + nzt-i+1
      zts2(it) = z(node2)              
      if(i.gt.1 .and. t2(node2).lt.ts) then
         node1 = node2 + 1
         if(dabs(t2(node2)-t2(node1)).gt.1.0e-5)
     +    zts2(it) = z(nzt)- (z(node2)-
     +              ((t2(node2)-ts)/(t2(node2)-t2(node1)) 
     +             * (z(node2)-z(node1))))
         goto 1447
      endif
      enddo
 1447 continue
     
      endif 
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     moving average for heater columns
c     hardwire only for >0 columns (to target Carlos' case)
c     assumes columns extend from 1-nez
c     ------------------------------------------------------
      if(numiscol.gt.0) then
      kkk=0
      do k=1,nezt
      tavg=0.
      do i=1,numiscol      !loop horizontally
      lis = (i-1)*nezt+k   !local element index for 1-numis
      l = ifl(lis)         !l is true 3D element number
      if(z(in2(l,1)).lt.0.12) goto 4746
      if((zl-z(in2(l,8))).lt.0.12) goto 4746
      lbasegrid=mod(l,nexy)
      if(lbasegrid.eq.0) lbasegrid=nexy
      do kk=l,lbasegrid,-nexy             !find bottom element
      if( (z(in2(l,8))-z(in2(kk,1))).ge.0.06) goto 7444 
      enddo
7444  lbot=kk
      ltopgrid = lbasegrid+(nezt-1)*nexy
      do kk=l,ltopgrid,nexy               !find top element
      if( (z(in2(kk,8))-z(in2(l,1))).ge.0.06) goto 7445
      enddo
7445  ltop=kk
c      write(6,*) 'k,i,l,lbasegrid,lbot,ltopgrid,ltop: ',
c     +            k,i,l,lbasegrid,lbot,ltopgrid,ltop 
      do kk=lbot,ltop,nexy
      do j=1,8
      tavg=tavg+t2(in2(kk,j))
      enddo
      enddo
      enddo                        !end loop horizontal
      kkk=kkk+1
      zdts(kkk) = (z(in2(l,1))+z(in2(l,8)))/2.
      tdts(kkk,it)=tavg/float((numiscol*(((ltop-lbot)/nexy)+1))*8)
4746  continue
      enddo                        !end loop nez
      numavgpts=kkk
      write(55,2829) kkk,tdays,(zdts(i),tdts(i,it),i=1,numavgpts)
2829  format('zone i=',i5,',T=" ',e12.4,' days"',/, (2e15.5))
      call flush(55)
      endif
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c

      call flush(6)
      call flush(8)
      call flush(9)
      call flush(14)
      call flush(16)
      call flush(18)
      call flush(19)
      call flush(20)
      call flush(21)
      call flush(25)
      call flush(26)
      call flush(27)
      call flush(28)
      call flush(29)
      call flush(30)
      call flush(40)
      call flush(46)
      call flush(54)
      call flush(60)
      call flush(62)

c     check for end of time loop with constant dt ...
c     -----------------------------------------------
c
      if((time1-time).gt.(time1*1.0e-6)) goto 1001
      if(more_time.gt.0) goto 1000
c
c     end of program ...
c     ------------------
c
      write(6,4040) numit
 4040 format(/10x,'# time steps numit ... ',i6)
c
c    'breakthrough' of column average data
c     hardwire ... every iskip3 time steps      
c     ------------------------------------
      if(numiscol.gt.0) then
      iskip3=1
      do i=1,numavgpts
      write(56,9225) numit/iskip3,zdts(i)
 9225 format('zone i= ',i6,', T="z= ',f12.4,'"') 
      write(56,9226) (btime(k),tdts(i,k),k=1,numit/iskip3)
 9226 format(2e15.6)
      enddo
      endif
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     breakthrough at selected nodes to file #15
c     print heads if system is saturated (ksat=1)
c     ---------------------------------------------
      write(15,377) nwells
 377  format('title="breakthrough data set at ',i5,' wells"',
     +     /,'variables="time","Havg","Tavg","Tpeak","Tmin",'
     +     /,'"cpkwellc1","cpkwellc2","spnwell"')
c477  format('title="breakthrough data set at ',i5,' wells"',
c    +      /,'variables="time","Tavg","Tpeak"')
      do i=1,nwells
      write(15,378) numit,i
 378  format('zone i= ',i9,', T="monitor ',i5,'"')
      write(15,6205)
     +(btime(k),ubh(i,k),ubt(i,k),tpkwell(i,k),tminwell(i,k),
     +          cpkwell(i,k,1),cpkwell(i,k,2),spnwell(i,k),k=1,numit)
c    +write(15,6206) (btime(k),ubt(i,k),tpkwell(i,k),k=1,numit)
 6205 format(8e15.6)
      enddo
      call flush(15)
c
c     print freezing front location: x
c     ---------------------------------
      write(42,2234) numit,ts
 2234 format('title="freezing front location: 1dx at knoy(1),knoz(1)"',
     +  /,'variables="time(days)","x(m)"',
     +  /,'zone i= ',i7, ', T="ts=',f10.2,'"')
      write(42,2235) (btime(k),xts(k),k=1,numit)
 2235 format(2e15.5)     
c
c     print freezing/thawing front location: z
c     -----------------------------------------
      write(58,21234) numit,ts
21234 format(
     +   'title="freeze/thaw front location: 1dz at knox(1),knoy(1)"',
     +  /,'variables="time(days)"," z_T>0","z_T<0"', 
     +  /,'zone i= ',i7,', T="ts=',f10.2,'"')
      write(58,21235) (btime(k),zts1(k),zts2(k),k=1,numit)
21235 format(3e15.5)     
c
c     terminate files with line totals
c     --------------------------------
      write(18,3737) momt
      write(19,3737) numit
      write(43,3737) numpt43
      write(49,3737) numpt49
      write(54,3737) numit
      close(43)
      close(49)
      open(43,file='t1dz_tdum.plt')
      open(49,file='t1dx_tdum.plt')
      open(143,file='t1dz_time.plt')
      open(149,file='t1dx_time.plt')
c
c     fix file 43 (t1dz_time.plt)
c     -------------------------
      read(43,10) title1,title2,title3
      write(143,10) title1,title2
      write(143,3743) nzt,numpt43
 3743 format('zone i=',i5,',j= ',i8,',T="at knox1,knoy1"')       

      do i=1,1000000
      read(43,*,err=99143) dum1,dum2,dum3,dum4
      write(143,803)  dum1,dum2,dum3,dum4
      enddo
99143 continue
      close(143)
c
c     fix file 49 (t1dx_time.plt)
c     --------------------------
      read(49,10) title1,title2,title3
      write(149,10) title1,title2
      write(149,3749) nx,numpt49
 3749 format('zone i=',i5,',j= ',i8,',T="at knoy1,knoz1"')       

      do i=1,1000000
      read(49,*,err=99149) dum1,dum2,dum3,dum4
      write(149,803)  dum1,dum2,dum3,dum4
      enddo
99149 continue
      close(149)

 3737 format('number of print times: ',i6,' (replace 9999 above)')

c
c     print 2D section velocities at final time
c     ------------------------------------------
c      if(ney.eq.1) then
c      write(12,678) tdays
c  678 format(' 2d velocity field tdays= ',f10.2)
c      do 7722 l=1,net
c      xbar=0.
c      zbar=0.
c      do 7723 j=1,8
c      xbar=xbar+x(in(l,j))
c      zbar=zbar+z(in(l,j))
c 7723 continue
c      xbar=xbar/8.
c      zbar=zbar/8.
c      write(12,6623) xbar,zbar,vx(l),vz(l)
c 6623 format(4e17.6)
c 7722 continue
c      endif
c
c     print fracture velocities at final time
c     hardwire - print all velocities
c     later: separate into slices, check for orientation
c     ----------------------------------------------------
      write(22,1678) tdays,nfrac
 1678 format(' title="xz fracture velocity field tdays= ',f10.3,'"',
     +       /,'variables="x","z","vx","vz"',
     +       /,'zone i= ',i10)
      write(33,1679) tdays,nfrac
 1679 format(' title="xy fracture velocity field tdays= ',f10.3,'"',
     +       /,'variables="x","y","vx","vy"',
     +       /,'zone i= ',i10)

c
c     for discrete random fractures in xy plane:
c     (but cannot do particle tracks)
c     write(52,1677) tdays,nfknoz1
c     then:   /,'zone i= ',i10) 
c     -----------------------------
      write(522,1677) tdays,nfknoz1
 1677 format('title="xy discrete fracture velocities on knoz1 surface, '
     +         'tdays= ',f10.3,'"',
     +       /,'variables="x","y","vx","vy"',
     +       /,'zone i= ',i10)

      vzero=0.
      do 8722 kf=1,nfrac
      l = ifracl(kf)
      if(ifdim(kf).eq.2) then      
      call frac_plane(maxfrac,maxne,lvert,inline,in,kf,l,fdimx,fdimy,
     +             exl,eyl,ezl,rexel,sarea)
      xbar= 0.
      ybar= 0.
      zbar= 0.
      do 8723 j=1,4
      xbar=xbar+x(map(inline(j)))
      ybar=ybar+y(map(inline(j)))
      zbar=zbar+z(map(inline(j)))
 8723 continue
      xbar=xbar/4.
      ybar=ybar/4.
      zbar=zbar/4.
      if(lvert(kf).eq.1.or.lvert(kf).eq.2)      
     +write(22,8623) xbar,zbar,vzero,vy2d(kf)      !2D vertical yz fractures shown in xz plane 
      if(lvert(kf).eq.5.or.lvert(kf).eq.6)      
     +write(22,8623) xbar,zbar,vx2d(kf),vzero      !2D horizontal fracture shown in xz plane

 8623 format(4e15.6)
      if(lvert(kf).eq.1.or.lvert(kf).eq.2)      
     +write(33,8623) xbar,ybar,vzero,vx2d(kf)       !2D vertical yz fractures shown in xy plane            
      if(lvert(kf).eq.3.or.lvert(kf).eq.4)      
     +write(33,8623) xbar,ybar,vx2d(kf),vzero       !2D vertical xz fractures shown in xy plane

      do kk=1,nfknoz1
      if(kf.eq.ifrac2dz1(kk))write(522,8623) xbar,ybar,vx2d(kf),vy2d(kf)    !2D horizontal fractures shown in xy plane at knoz1
      enddo

      else
      call frac_line(maxfrac,maxne,lvert,inline,in,kf,l,fracdim,
     +             exl,eyl,ezl)
      xbar= 0.
      ybar= 0.
      zbar= 0.
      do j=1,2
      xbar=xbar+x(map(inline(j)))
      ybar=ybar+y(map(inline(j)))
      zbar=zbar+z(map(inline(j)))
      enddo
      xbar=xbar/2.
      ybar=ybar/2.
      zbar=zbar/2.
      if(lvert(kf).eq.1) write(22,8623) xbar,zbar,vzero,vlin(kf)
      if(lvert(kf).eq.2) write(22,8623) xbar,zbar,vlin(kf),vzero
      endif      
 8722 continue
c
c     thermal or mass recovery:
c     -------------------------
      if(kcntrl.eq.0.and.lss) then
      if(abs(hic).lt.1.e-20) then 
        write(6,272)
  272   format(/10x,'thermal energy recovery undefined',/)
      else 
      if(abs(hic).ge.1.e-20) tr=abs(hoc/hic) * 100.
      if(kmass.eq.0.or.kmass.eq.3) write(6,172) hic,hoc,tr
  172 format(/10x,'total thermal energy injected  above ambient:',e12.3,
     +       /10x,'total thermal energy recovered above ambient:',e12.3,
     +       /10x,'thermal recovery efficiency (%) ... ',e12.3)
      if(kmass.eq.1.or.kmass.eq.3) write(6,372) hic,hoc,tr
  372 format(/10x,'total mass injected  above background :',e12.3,
     +       /10x,'total mass recovered above background :',e12.3,
     +       /10x,'mass recovery efficiency (out/in %)   :',e12.3)
      endif
      endif
c     
c     mass balance over time
c     ----------------------
      write(63,*) 'Title="Mass balance, mobile and retained particles"'
      write(63,3399) numit
 3399 format('variables="time(days)","kg suspended particles",',
     + '"kg retained particles","Total mass (kg)","fluxin(nzt)kg",',
     + '"fluxout(1)kg","M_balance(kg)"',
     + /,'zone i= ',i7)
      do k=1,numit
      sum=cmass(k)+smass(k)
      balance = fluxin(k) - fluxout(k) - sum
      write(63,3772) btime(k),cmass(k),smass(k),sum,
     +                fluxin(k),fluxout(k),balance 
 3772 format(7e12.4)
      enddo
c     ---------------------------     
c
      write (6,170) itot
  170 format(/10x,'total iterations performed: ',i5)
c
c     write(6,270) fast,fpcgt,tast,tpcgt,vcpu
c270  format(/10x,'cpu run times (secs) ...',
c    +       /10x,'flow matrix assembly ... ',e12.3,
c    +       /10x,'flow pcg solution   .... ',e12.3,
c    +       /10x,'transport matrix assembly ... ',e12.3,
c    +       /10x,'transport pcg solution   .... ',e12.3,
c    +       /10x,'velocity subroutine ......... ',e12.3)
      write (6,171)
      write (*,171)
  171 format (//10x,'**** normal exit ****'//)
c
c     closing ...
c     -----------
      close(5)
      close(6)
      close(7)
      close(8)
      close(9)
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)
      close(20)
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)
      close(30)
      close(31)
      close(33)
      close(34)
      close(35)

      close(43)
      close(44)
      close(45)
      close(49)
      close(51)
      close(52)
      close(522)
      close(54)
      close(57)
c

      close(60)
      close(61)
      close(62)
      close(63)
c
      stop
5001  write(6,5002) 
5002  format(/10x,'error reading time step data ...',
     +       /10x,'program stopping ...')
      stop

8484  write(6,8485)
8485  format(/10x,'error reading tcon,porsurf,satsurf',
     +       /10x,'program stopping ...')
      stop     
8454  write(6,8455)
8455  format(/10x,'error reading kprt,kcntrl,kwt,kint,kintv,kgo,ksat',
     +            ',kmass,ktair,krk,kdz,ksubtr,kr1top,kfreec,krfrac',
     +       /10x,'program stopping ...')
      stop            
3030  write(6,8486)
8486  format(/10x,'error reading breakthrough points',
     +       /10x,'program stopping ...')
      stop     
3031  write(6,8487)
8487  format(/10x,'error reading hinc,rinc,rtdiff',
     +       /10x,'program stopping ...')
      stop     
3032  write(6,8488)
8488  format(/10x,'error reading temperature initial conditions',
     +       /10x,'program stopping ...')
      stop    
3142  write(6,8588)
8588  format(/10x,'error reading concentration initial conditions',
     +       /10x,'program stopping ...')
      stop           
3033  write(6,8489)
8489  format(/10x,'error reading fgen92.asc K field',
     +       /10x,'program stopping ...')
      stop     
3034  write(6,8490)
8490  format(/10x,'error reading storage term ss',
     +       /10x,'program stopping ...')
      stop    
3135  write(6,8590)
8590  format(/10x,'error reading init,ho (initial head options)',
     +       /10x,'program stopping ...')
      stop    
3035  write(6,8491)
8491  format(/10x,'error reading air temperatures in file tair.data'
     +       /10x,'program stopping ...')
      stop    
3036  write(6,8492)
8492  format(/10x,'error reading bdy heads in file hbdy.data'
     +       /10x,'program stopping ...')
      stop   
3037  write(6,8493)
8493  format(/10x,'error reading veloc2d_in.data; ivel=3'
     +       /10x,'program stopping ...')
      stop    
3038  write(6,8494)
8494  format(/10x,'error reading cf,pf ... phi, alpha0, alphap '
     +       /10x,'program stopping ...')
      stop           
3139  write(6,8595)
8595  format(/10x,'error reading (no p,q) ts,omega,modelwu,modelkr'
     +       /10x,'program stopping ...')
      stop           
3039  write(6,8495)
8495  format(/10x,
     +          'error reading fracture p2,q2,ts2,omega,modelwu,modelkr'
     +       /10x,'program stopping ...')
      stop           
3040  write(6,8496)
8496  format(/10x,'error reading bz data in file bzin.data'
     +       /10x,'program stopping ...')
      stop
3041  write(6,8497)
8497  format(/10x,'error reading indexed K data in smoker.data'
     +       /10x,'program stopping ...')
      stop
3043  write(6,8598)
8598  format(/10x,'error reading fgen92.asc porosity field',
     +       /10x,'program stopping ...')
      stop     

8498  write(6,8499)
8499  format(/10x,'error reading surfat variables (tair,bz or sat) ...',
     +       /10x,'surfmin,amp,period,phase,cutoff,gradt,gradx',      
     +       /10x,'program stopping ...')
      stop

 1119 write(6,1120)
 1120 format(/10x,'end of file detected in grid3d.plt',
     +       /10x,'this file is either empty or not compatible with',
     +       /10x,'your existing data file - please check grid size',
     +       /10x,'and source of your existing data file',
     +       /10x,'program stopping ...')
      stop

10002 write(6,8577)
 8577 format(/10x,'error in reading grid3d.plt',
     +     /10x,'check spacing for data read:',
     +     /10x,'zone i=    nz, j=    ny, k=    nx (7x,i8,4x,i8,4x,i8)'
     +     /10x,'program stopping')
      stop   
c
9922  write(6,8578)
 8578 format(/10x,'error in reading internal heat flux',
     +     /10x,'check i1,i2,j1,j2,k1,k2,fluxin1,fluxin2,decay2,more',
     +     /10x,'program stopping')
      stop        
 6661 write(6,8579)  
 8579 format(/10x,'error in reading fracture data',
     +/10x,'check ilx1,ilx2,ily1,ily2,ilz1,ilz2,idim,iloc,aperture,more'
     +  ,/10x,'program stopping')
      stop       
c
c   35 read (5,*,err=6662) j1,j2,ckx,cky,ckz,tcx,cps,ps0,theta,satw,
c     +                    depc1,depd1,vexp1,smax0,pp0,qq0,more
 
 6662 write(6,8580)  
 8580 format(/10x,'error reading line with elemental Kx,Ky,Kz,'
     +       /10x,'check j1,j2,ckx,cky,ckz,tcx,cps,ps,theta,satw,'
     +            'depc,depd,vexp,smax,pp,qq,more',
     +       /10x,'program stopping')
      stop        
c
 8001 write(6,8581)  
 8581 format(/10x,'error in reading convergence criteria,',
     +       /10x,'check ccp,cct,ccw,ccc,maxit1,maxit2',
     +       /10x,'program stopping')
      stop

 1137 write(6,1138) neplane2d
 1138 format(/10x,'!!! error detected in random2d.plt !!! ',
     + /10x,'you are trying to read an external aperture field'
     + /10x,'but this file is not compatible',
     + /10x,'read format: 3 titles, then:',
     + /10x,'             (dum1,dum2,fap(i),dum3,dum4),i=1,neplane2d)',
     + /10x,'             neplane2d = ',i8,
     + /10x,'use krfrac=0 or check random2d.plt file',
     + /10x,'program stopping ...')
       stop
c
 1140 write(6,8582)  
 8582 format(/10x,'error in reading dispersivities,rtrans etc,',
     + /10x,'(need rtrans c1,c2 and decay c1,c2)'
     + /10x,'check alh,alv,ath,atv,dd,rtrans_c1,rtrans_c2,',
     +'decay_c1,decay_c2,agefx,kdisp',
     +       /10x,'program stopping')
      stop
c     
      end
c
c *******************************************************************
c
c *******************************************************************
c
c
c  flow.for
c
c  ********  f l o w    ******
c
c  three-dimensional transient density-dependent flow
c  linear isoparametric quadrilateral elements
c  conjugate gradient solver
c
c     e.o. frind, j.w. molson
c     centre for groundwater research
c     university of waterloo
c     january 1990
c
c   ******************************************************************
c
c       copyright 1990:  emil o. frind, j.w. molson
c
c       duplication of this program or any part thereof without
c       the express written consent of the authors is prohibited.
c
c   ******************************************************************
c
c   this package contains ...
c
c           flow.for    ...    main flow subroutine
c           gquad1      ...    gauss quadrature integration subroutine
c   __________________________________________________________________
c
c     maxne = number of elements
c     maxnn = number of nodes
c     maxn  = number of degrees of freedom
c     maxnb = non-zero bandwidth in case of conjugate gradient solver
c     maxna = total number of non-zero matrix entries in condensed matri
c     nf    = maximum number of nodes on one face
c     laa   = 3*n + na
c
c   ********************************************************************
c
      subroutine flow(maxn,maxnn,maxne,nf,maxna,laa,np1,xi,yi,se,
     + pe,x,y,z,u,u0,fc,fb,f,ic,lc,cx,cy,cz,in,in2,
     + a,iaa,ind,ib,aa,nn,n,nb,ne,kprt,
     + dt1,fst,temp0,temp,c0,c2,ss,pq,ag,hag,hxx,hyy,hzz,htmf,hgz,
     + kint,fast,fpcgt,tdays,map,exl,eyl,ezl,nexy,por,gamma,pi,pw,vt,
     + ckl,xarea,ifracl,lvert,nfrac,maxfrac,nex,ney,nez,ifdim,
     + hx,hy,ht,hgz2,kblck,pp,qq,modelwu,lheat,lmass,lage,omega,modelkr,
     + p2,q2,ts2,omega2,modelwu2,modelkr2,lkr1top,
     + por0,alpha0,alphap,rhob,ps,spn0,spn1,spe1)

      implicit real*8(a-h,o-z)
      external wu,dwu
c
c -------------------------------------------------------------------
c
c     constant dimensions
c     -------------------
      real*8 start,finish
      dimension xi(8),yi(8),se(8,8),pe(8,8),inl(8),ge(8)
      dimension ag(9),hag(9)
      dimension hxx(8,8),hyy(8,8),hzz(8,8),htmf(8,8),hgz(8)
c
c
c     arrays of size nn
c     -----------------
      real*8 pq(maxnn)
      real*8 x(maxnn),y(maxnn),z(maxnn)
      real*8 u0(maxnn),fc(maxnn),fb(maxnn),fst(maxnn)
      real*8 exl(maxne),eyl(maxne),ezl(maxne)
      integer   ic(maxnn),lc(maxnn),map(maxnn)
      dimension u(maxnn),f(maxnn),temp(maxnn),temp0(maxnn),vt(maxne)
      dimension c0(maxnn),c2(maxnn)
      dimension spn0(maxnn),spn1(maxnn),spe1(maxne)
c
c     1D line or 2D plane fracture elements
c     --------------------------------------
      dimension rel(4,4),inline(4),sel(4,4),gel(4)
      dimension ckl(maxfrac),xarea(maxfrac)
      dimension ifracl(maxfrac),ifdim(maxfrac)
      dimension lvert(maxfrac)
      dimension hx(4,4),hy(4,4),ht(4,4),hgz2(4)
c
c     arrays of size ne, (ne,8)
c     -------------------------
      real*8  cx(maxne),cy(maxne),cz(maxne),por(maxne)
      integer   in(maxne,8),in2(maxne,8),inl2(8)
      real*8 rhob(maxne),por0(maxne),ps(maxne)
      real*8 pp(maxne),qq(maxne)
c
c     array of size nf
c     ----------------
c     dimension fbf(nf)
c
c     arrays a(na),iaa(na),ind(n+1),ib(n,nb),aa(laa)
c     ... where laa=3n+na, for conjugate gradient solver only
c     -------------------------------------------------------
      dimension a(maxna),aa(laa)
      integer   iaa(maxna),ind(np1),ib(maxn,nb)
      integer kblck(maxnn)

      logical lheat,lmass,lage,lkr1top
c
c ====================================================
c
c function defining relative viscosity and relative density
c  of the fluid as a function of temperature (c).
c viscosity is defined relative to 0 celsius, and density
c  is relative to 4 celsius (temp. of maximum density).
c
c *********************************************************************
c

      lblo=maxnn

c     initialize all arrays ...
c     -------------------------
      do 45 i=1,nn
      kblck(i)=0
      f(i) = 0.
      fst(i) = 0.
      if (ic(i).eq.1) go to 45
      ii=i-lc(i)
      f(ii)=fb(i)
   45 continue
c
c     conjugate gradient matrices
c     ---------------------------
  150 continue
c      m=(nb-1)*n
      do 47 i=1,maxna
      a(i)=0.
   47 continue
c
c  =====================
c     tsurfmin=999.
c     tsurfmax=-999.
c     do i=1,nn
c     tsurfmin = min(tsurfmin,temp(i))
c     tsurfmax = max(tsurfmax,temp(i))           
c     enddo
c
c     write(6,3374) tsurfmin,tsurfmax
c3374 format(/10x,'in flow min,max temps : ',2e15.5)
c  =====================

c
c     loop over elements
c     ------------------
c     call clock@(start)
c
      do 100 l=1,ne
c
c     elemental matrices
c     -------------------
c     spe=0.
c     do i=1,8
c     inl(i)=in(l,i)
c     inl2(i) =map(inl(i))
c     spe = spe + (spn0(inl(i))+spn1(inl(i)))/2.
c     enddo
c     spe = spe/8. 
      spe=spe1(l)
c
c     update ice saturation dependent conductivity (from dnapl):
c     p is residual unfrozen water content
c     ---------------------------------------------------------
c     zkrw=1.d0
      swe = 0.d0
       do i=1,8
       inl(i)=in(l,i)
       inl2(i) =map(inl(i))
       swe=swe+wu(temp(inl2(i)),pp(l),qq(l),modelwu,ts)
       enddo
      swe=swe/8.d0
c      zkrw = ((swe - p)/(1.d0-p))**4
c      zkrw = max(zkrw,1.e-6)

c     print*, l,alpha0,alphap
c     print*, por(l)
c     print*, rhob(l)
c     print*, ps(l)
c     print*, spe


      alpha = alpha0 + por(l)*alphap*(rhob(l)/ps(l))*spe

c     if(l.eq.ne/2) write(6,9229) 
c    +           por(l),por0(l),rhob(l),ps(l),spe,alphap,alpha0,alpha
c9229 format(2x,'1:por,por0,rhob,ps,spe,alphap,alpha0,alpha:',8e10.3)
 
      zkrw = fnzkrw(swe,pp(l),omega,por(l),modelkr,por0(l),alpha0,alpha)             !relative permeability function
      if(lkr1top .and. l.ge.(ne-nexy+1)) zkrw=1.d0                    !hardwire keep kr=1 for top surface elements
      ss1 = ss*swe

c     if(l.eq.ne/2) write(6,9239) 
c    + por(l),por0(l),rhob(l),ps(l),spe,alphap,alpha0,alpha,zkrw
c9239 format(2x,'2:por,por0,rhob,ps,spe,alphap,alpha0,alpha,zkrw:',
c    +      9e10.3)

c      
c     set storage = porosity for top layer
c     ------------------------------------
c     if(l.gt.(ne-nexy).and.ss.gt.0.0) ss1=por(ne)        !deactivated 
c
c     update temp. dependent conductivity:
c     adjust for ice saturation dependent relative k
c     ---------------------------------------------
      tavg = 0.d0
      cavg = 0.d0
      do  i=1,8
      tavg = tavg + temp(inl2(i))
      cavg = cavg + c2(inl2(i))
      enddo
      tavg = tavg * 0.125d0
      cavg = cavg * 0.125d0
      tx=(cx(l)/rvisc(tavg,lheat,lmass,lage)) * zkrw
      ty=(cy(l)/rvisc(tavg,lheat,lmass,lage)) * zkrw
      tz=(cz(l)/rvisc(tavg,lheat,lmass,lage)) * zkrw

c     if(l.eq.ne/2) then
c       write(6,3833) l,modelkr,swe,p,omega,por(l),por0(l),alpha0,
c    +          alpha,zkrw,cx(l),tavg,cavg,rvisc(tavg,lheat,lmass,lage)
c3833 format(10x,'l,modelkr,swe,p,omega,por(l),por0(l),alpha0,alpha,
c    +zkrw,cx(l),tavg,cavg,rvisc(tavg,lheat,lmass,lage):'
c    +,2i5,12e12.3)
c     endif
c
c     numerical (kint=1) or exact (kint=0) integration ... flow
c     ----------------------------------------------------------
       if(kint.eq.1)
     + call gquad0(tx,ty,tz,x,y,z,inl2,se,pe,l,maxnn,ss1,dt1,ge,         ! inl2 ok because x,y,z are needed
     +  ag,hag)
c
       if(kint.eq.0)
     + call exint0(tx,ty,tz,x,y,z,inl,se,pe,l,maxnn,ss1,dt1,ge,
     + hxx,hyy,hzz,htmf,hgz,exl,eyl,ezl,maxne)
c
c      if((l.eq.8218.or.l.eq.10330).and.(tdays.eq.6..or.tdays.eq.20.)) 
c     + then
c      write(6,281) tdays
c  281 format(/10x,' flow matrix check at tdays= ',f10.3)
c      write(6,181) l,(i,(se(i,j),j=1,8),i=1,8)
c      write(6,181) l,(i,(pe(i,j),j=1,8),i=1,8)
c      write(6,181) l,(i,(ge(i),i=1,8))
c 181  format(/,' flow matrix check- se,pe,ge ...',i6,/(i2,8(e11.3,1x)))
c      endif
c
c     assembly for conjugate gradient solver
c     ---------------------------------------
      do 81 i=1,8
      ki = inl(i)
      if (ic(ki).eq.1) go to 81
      ii = ki-lc(ki)
      do 82 j=1,8
      kj = inl(j)
      if (ic(kj).eq.1) go to 85
      jj=kj-lc(kj)
      if (ii.gt.jj) go to 82
      do 83 k=1,nb
      if (ii.ne.ib(jj,k)) go to 83
      kk=k
      go to 84
   83 continue
   84 continue
      iv=ind(jj)+kk-1
c     write (6,738) l,i,j,ii,jj,ib(jj,1),ib(jj,kk),ind(jj),kk,iv,
c    1a(iv),se(i,j)
c 738 format (10i5,2f14.8)
c     -----------------------------------------------------------------
c     lhs assembly {a}, rhs storage+density term {fst}, and
c     linking term [se]*{fc}:
c     -------------------------------
      a(iv)= a(iv) + se(i,j) + pe(i,j)
      if(i.eq.j)fst(ii)=fst(ii)+pe(i,j)*u0(kj)-ge(i)*tz
     +                         * rden(tavg,cavg,gamma,lheat,lmass,lage)
      go to 82
   85 continue
      f(ii)=f(ii)-se(i,j)*fc(kj)
   82 continue
   81 continue
c
   90 continue
c
c     end of element loop
c     -------------------
  100 continue
c
c     fracture element assembly - flow
c     --------------------------------
      porfrac=1.d0
      do 882 kf=1,nfrac
c
c     1D line elements:
c     ------------------------------------------------------------------
      if(ifdim(kf).eq.1) then
      l = ifracl(kf)
      call frac_line(maxfrac,maxne,lvert,inline,in,kf,l,fracdim,
     +             exl,eyl,ezl)
c
      tbar = 0.50 * (temp(map(inline(1))) + temp(map(inline(2))))
      cbar = 0.50 * (c2(map(inline(1))) + c2(map(inline(2))))
      spe = 0.50 * (spn0(map(inline(1))) + spn1(map(inline(2))))
c     gamma=tbar/1000.  !outxx
c     viscosity = 160. * rvisc(tbar,lmass,lage)/ 86400.
      gel(1)=0.
      if(lvert(kf).eq.1 .or. lvert(kf).eq.4) gel(1) = - xarea(kf)
      gel(2) =  - gel(1)

c     hardwire 1D fracture Kr: set p=0.01, q=5, porosity =1.0 (Interfrost:p,q=0.05,0.5)
      swe = wu(tbar,p2,q2,modelwu2,ts2)
      alpha = alpha0 + por(l)*alphap*(rhob(l)/ps(l))*spe
      zkrw = fnzkrw(swe,p2,omega2,porfrac,modelkr,por0(l),alpha0,alpha)      !relative permeability function- fracture

      dfel = xarea(kf) 
     +       * (ckl(kf)/(rvisc(tbar,lheat,lmass,lage)*fracdim))*zkrw
      ssel = ss * fracdim * (dt1/2.d0)  * xarea(kf)
c
c      write(6,333) l,dfel,ssel
c 333  format(/10x,'working on line element: element# ',i5,
c     +       /10x,' dfel: ',e12.3,' ssel: ',e12.3) 
c
      rel(1,1) = dfel 
      rel(1,2) = -dfel
      rel(2,1) = -dfel
      rel(2,2) = dfel
      sel(1,1) = ssel
      sel(1,2) = 0.d0
      sel(2,1) = 0.d0
      sel(2,2) = ssel
c
      do 2012 i=1,2
      ki = inline(i) 
      if(ic(ki).eq.1) goto 2012
      ii= ki-lc(ki)
        do 2013 j=1,2
        kj = inline(j)
        if(ic(kj).eq.1) goto 2014
        jj = kj - lc(kj)
        if(ii.gt.jj) goto 2013  
c
      do 2083 k=1,nb
      if (ii.ne.ib(jj,k)) go to 2083
      kk=k
      go to 2084
 2083 continue
 2084 continue
c
      iv=ind(jj)+kk-1
      a(iv)= a(iv) + rel(i,j) + sel(i,j)
      if(i.eq.j) 
     +fst(ii)=fst(ii) + sel(i,j)*u0(kj)-gel(i)*ckl(kf) * zkrw
     +   * rden(tbar,cbar,gamma,lheat,lmass,lage)
     +                                    /rvisc(tbar,lheat,lmass,lage)
      goto 2013
 2014 continue
      f(ii) = f(ii) - rel(i,j)*fc(kj)
c
 2013 continue
 2012 continue
 2011 continue
c
c     end 1D line elements
c     ---------------------
      endif
c
c     2D plane elements - flow
c     ------------------------
      if(ifdim(kf).eq.2) then
      l = ifracl(kf)
      call frac_plane(maxfrac,maxne,lvert,inline,in,kf,l,fdimx,fdimy,
     +             exl,eyl,ezl,rexel,sarea)
c     write(6,2200) l,fdimx,fdimy,rexel,sarea
c2200 format(/10x,'assembly of 2D fracture elements - flow, #',i5,
c    +       /10x,'fdimx= ',e10.2,' fdimy= ',e10.2,
c    +       /10x,'rexel= ',e10.2,' sarea= ',e10.2)
c
      tbar = 0.25 * (temp(map(inline(1))) + temp(map(inline(2)))
     +              + (temp(map(inline(3))) + temp(map(inline(4)))))
      cbar = 0.25 * (c2(map(inline(1))) + c2(map(inline(2)))
     +              + (c2(map(inline(3))) + c2(map(inline(4)))))
      spe = 0.25 * (spn1(map(inline(1))) + spn1(map(inline(2)))
     +              + (spn1(map(inline(3))) + spn1(map(inline(4)))))
c changes March 1997: sarea removed from denomin dfel, /36 in ssel

c     hardwire 2D fracture Kr: set p=0.01, q=5, porosity=1.0 (Interfrost:p,q=0.05,0.5)
      swe = wu(tbar,p2,q2,modelwu2,ts2)
      alpha = alpha0 + por(l)*alphap*(rhob(l)/ps(l))*spe
      zkrw = fnzkrw(swe,p2,omega2,porfrac,modelkr2,por0(l),alpha0,alpha)      !relative permeability function- fracture

      dfel = xarea(kf) * (ckl(kf)/rvisc(tbar,lheat,lmass,lage)) * zkrw
      ssel = ss * dt1  * xarea(kf) * sarea/36.d0
      rex = dfel*rexel/6.
      rey = dfel/(6.*rexel)
c
c     write(6,333) l,dfel,ssel,rex,rey,tbar
c 333 format(/10x,'working on 2D plane element# ',i5,
c    +       /10x,' dfel: ',e12.3,' ssel: ',e12.3,
c    +       /10x,' rex:  ',e12.3,' rey:  ',e12.3,' tbar: ',e12.3)
c
c     assemble elemental matrices - no gravity term for face 5 or 6
c     -------------------------------------------------------------
      do 445 i=1,4
      gel(i) = xarea(kf) * hgz2(i) * fdimx/2.
      if(lvert(kf).eq.5 .or. lvert(kf).eq.6 ) gel(i) = 0.
      do 445 j=1,4
      rel(i,j) = rex*hx(i,j) + rey*hy(i,j)
      sel(i,j) = ssel * ht(i,j)
 445  continue
c
      do 3012 i=1,4
      ki = inline(i)
      if(ic(ki).eq.1) goto 3012
      ii= ki-lc(ki)
        do 3013 j=1,4
        kj = inline(j)
        if(ic(kj).eq.1) goto 3014
        jj = kj - lc(kj)
        if(ii.gt.jj) goto 3013
c
      do 3083 k=1,nb
      if (ii.ne.ib(jj,k)) go to 3083
      kk=k
      go to 3084
 3083 continue
 3084 continue
c
      iv=ind(jj)+kk-1
      a(iv)= a(iv) + rel(i,j) + sel(i,j)
      if(i.eq.j)
     +fst(ii)=fst(ii) + sel(i,j)*u0(kj)-gel(i)*ckl(kf) * zkrw
     +   * rden(tbar,cbar,gamma,lheat,lmass,lage)
     +                                 /rvisc(tbar,lheat,lmass,lage)
      goto 3013
 3014 continue
      f(ii) = f(ii) - rel(i,j)*fc(kj)
c
 3013 continue
 3012 continue
 3011 continue
c
c     end 2D plane elements
c     ----------------------
      endif
c
c     continue to next fracture element
c     ---------------------------------
 882  continue

c      write(6,8798) zkrmin,zkrmax
c8798  format(10x,'flow: zkrmin,zkrmax ...',2e12.4)      
c
c     assemble final flux vector
c     --------------------------
      do 51 i=1,nn
      if(ic(i).eq.1) goto 51
      ii=i-lc(i)
      f(ii) = f(ii) + fst(ii) + pq(map(i)) 
   51 continue
c
c     water saturation gradient term ... test lumped
c     ----------------------------------------------
c     dswmax = -999.
c     do i=1,nn
c     if(ic(i).eq.1) goto 511
c     pw=den(temp(map(i)),lmass,lage)
c     denr = (pi-pw)/pw 
c     dsw =  wu(temp(map(i)),p,q,modelwu,ts)
c    +                          - wu(temp0(map(i)),p,q,modelwu,ts)
c     ii=i-lc(i)
c     f(ii) = f(ii) + por(1)*denr*dsw*dt1*vt(1)  !will need nodal porosity,volume  +/- ?
c     dswmax=max(dsw,dswmax)
c511  continue
c     enddo
c     write(6,4848) por(1),dt1,denr,dswmax
c4848 format(10x,'por(1),dt1,denr,dswmax: ',4e15.4)   
c
c     do 733 i=1,maxna
c 733 write (6,734) i,iaa(i),a(i)
c 734 format (2i5,f14.8)
c
c     conjugate gradient solution - flow
c     -----------------------------------
      neu=1
c      call precg (a,n,iaa,ind,f,u,aa,laa,neu)
      call precgn(a,n,iaa,ind,f,u,aa,laa,kblck,lblo,neu,maxnn)
c
  130 continue
c
c     expand solution vector
c     ----------------------
      do 121 ii=1,nn
      i = nn-ii+1
      if (ic(i).eq.1) go to 122
      k = i-lc(i)
      u(i) = u(k)
      go to 121
  122 u(i) = fc(i)
  121 continue
c
c      write (7,132) (i,u(i),i=1,nn)
c  132 format (3(i10,e11.4))
c
      return
      end
c ==================================================================
c
      subroutine gquad0 (tx,ty,tz,x,y,z,inl2,se,pe,l,maxnn,ss,dt1,ge,
     +                   ag,hag)
c
c     generation of isoparametric element matrices
c     full gauss quadrature numerical integration
c     --------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension detj(27),f(8),ge(8),dgx(8),dgy(8),dgz(8),ag(9),hag(9),
     +       se(8,8),pe(8,8),ff(8,27),dx(8,27),dy(8,27),dz(8,27)
      real*8   x(maxnn),y(maxnn),z(maxnn)
      dimension inl2(8)
c
      m=8
c
c     number of gauss points
c     ----------------------
c      np=2
c      npx=np
c      npy=np
c      npz=np
c      np2=npx*npy*npz
c      ns=0
c      if (np.eq.3) ns=2
c
c     placement of gauss points
c     -------------------------
      do 310 k=1,2
      zi=ag(k)
      hz=hag(k)
      do 310 j=1,2
      yi=ag(j)
      hy=hag(j)
      do 310 i=1,2
      xi=ag(i)
      hx=hag(i)
      kl=(k-1)*4+(j-1)*2+i
c     if (l.eq.1) write (6,741) ii,ij,ik,kl,xi,yi,zi,hx,hy,hz
c 741 format (/' gauss point',4i5,6f14.7)
c
c     basis functions and derivatives
c     -------------------------------
      call sp3lin (x,y,z,inl2,m,xi,yi,zi,det,l,dgx,dgy,dgz,f,maxnn)
c
      do 300 jj=1,m
      ff(jj,kl)=f(jj)
      dx(jj,kl)=dgx(jj)
      dy(jj,kl)=dgy(jj)
  300 dz(jj,kl)=dgz(jj)
      detj(kl)=det*hz*hy*hx
  310 continue
c     if (l.ne.1) go to 311
c     write (6,314) l
c 314 format (/' basis functions element',i5)
c     write (6,313) (k,(dx(i,k),i=1,8),k=1,8)
c     write (6,313) (k,(dy(i,k),i=1,8),k=1,8)
c     write (6,313) (k,(dz(i,k),i=1,8),k=1,8)
c 311 continue
c 313 format (i5,8f14.8)
c
c     numerical integration of elemental matrices
c     --------------------------------------------
      do 350 i=1,m
      do 351 j=i,m
      s=0.d0
      p=0.d0
      g=0.d0
      do 325 k=1,8
c
c     coefficient matrix term:
c     ------------------------
      s = s +(tx*dx(i,k)*dx(j,k)+ty*dy(i,k)*dy(j,k)+tz*dz(i,k)
     +          *dz(j,k))*detj(k)
c
c     storage matrix term: lumped used here
c     p= p + ff(i,k)*ss*ff(j,k)*detj(k)  (use this for consistent)
c     -------------------------------------
      if(i.eq.j) p= p + ff(i,k)*detj(k)
c
c     density term:
c     -------------
      if(i.eq.j) g = g + dz(i,k)*detj(k)
  325 continue
      se(i,j)= s
      pe(i,j)= p
      if(i.eq.j) pe(i,j)= pe(i,j)*dt1*ss
      if(i.eq.j) ge(i) = g
  351 continue
  350 continue
c
c       if(l.le.5) then
c       write(6,887) l
c  887  format(1x,'ge array ... element',i5,/)
c       write(6,888) (ge(i),i=1,8)
c       write(6,889) l
c  889  format(1x,'se array ... element',i5/)
c       write(6,888) ((se(i,j),j=1,8),i=1,8)
c       write(6,889) l
c  889  format(1x,'pe array ... element',i5/)
c       write(6,888) ((pe(i,j),j=1,8),i=1,8)
c  888  format(4(/1x,2e10.3))
c       endif
c       if(l.eq.5) stop
c
c     fill lower half of se and pe arrays
c     -----------------------------------
       do 380 i=2,m
       i1=i-1
       do 380 j=1,i1
c      pe(i,j)=pe(j,i)     !hardwire-not sure why this was commented out. pe(6,3) becomes undefined
  380  se(i,j)=se(j,i)
c
c     write (6,21) l,(inl(i),i=1,8)
c  21 format (/' element',i5,5x,8i5)
c     write (6,22) (i,(se(i,j),j=1,8),i=1,8)
c  22 format (i10,8f12.6)
c
      return
      end
c
c ==================================================================
c
      subroutine exint0(tx,ty,tz,x,y,z,inl,se,pe,l,maxnn,ss,dt1,ge,
     +    hxx,hyy,hzz,htmf,hgz,exl,eyl,ezl,maxne)
c
      implicit real*8(a-h,o-z)
      real*8    x(maxnn),y(maxnn),z(maxnn)
      real*8 exl(maxne),eyl(maxne),ezl(maxne)
      dimension se(8,8),pe(8,8),ge(8)
      dimension inl(8)
      dimension hxx(8,8),hyy(8,8),hzz(8,8),htmf(8,8),hgz(8)
c
c     assembly of element matrices - exact integration
c     rectangular or orthogonal quadrilateral elements
c     lumped mass matrix for flow assumed
c     --------------------------------------------------
c
      exll=exl(l)
      eyll=eyl(l)
      ezll=ezl(l)
      txx=tx*eyll*ezll/exll/36.d0
      tyy=ty*exll*ezll/eyll/36.d0
      tzz=tz*exll*eyll/ezll/36.d0
      rtt=dt1*ss*exll*eyll*ezll/216.d0
      rgg=exll*eyll/4.
c
      do 156 i=1,8
      do 156 j=1,8
      pe(i,j)=0.d0
      se(i,j)=txx*hxx(i,j)+tyy*hyy(i,j)+tzz*hzz(i,j)
      if(i.eq.j) pe(i,j)=htmf(i,j)*rtt
      if(i.eq.j)   ge(i)=hgz(i)*rgg
  156 continue
c
      return     
      end
c **************************************************************
c
c *******************************************************************      New version 2014 for latent heat
c
c trans.for
c
c    ******** subroutine  transport  **********
c
c  three-dimensional advective-dispersive transport
c  linear isoparametric quadrilateral elements
c
c     e.o. frind, j.w. molson
c     university of waterloo
c     (c) 1993
c
c   *******************************************************************
c
c       copyright 1987 emil o. frind, 1988 j.w. molson
c
c       duplication of this program or any part thereof without
c       the express written consent of the author is prohibited.
c
c   *******************************************************************
c
      subroutine trans1(maxn,maxnn,maxne,nf,maxna,laa,nw,                     !heat transport ksubtr=1
     + x,y,z,u0,u1,u2,fc,fb,f,fs,icflow,ic,lc,vx,vy,vz,inflow,in,map,
     + a,aa,iaa,ind,ib,alh,alv,ath,atv,dd,kdisp,n,nbb,dt1,dt2,nz,
     + ne,nn,wp,wa,wp1,wa1,pq,tq,cq,por,sw,ag,hag,ttemp,ara,leak2,
     + hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,htmt,kint,tast,tpcgt,tdays,
     + exl,eyl,ezl,wlh,pp,qq,modelwu,vt,tcs,bz,bzf,cpsm,cpi,cf,pi,sat,
     + nef,
     + fbf,rinc,htmf,ifracl,lvert,nfrac,xarea,nex,ney,nez,vlin,maxfrac,
     + tclw,tcli,tclm,ifdim,hgz2,hx,hy,hh,ht,wx,wy,vx2d,vy2d,
     + cot0,cot1,tkl,numis,flux,fluxv,ifl,maxfx,rtrans,decay,decay2,
     + nflux,agefx,rg,rg2df,porsurf,lmass,lheat,lage,ts,rtdiff,
     + p2,q2,ts2,omega2,modelwu2,modelkr2,depc,vexp,smax,spn0,spn1,
     + ppn,qqn)
c      
      implicit real*8(a-h,o-z)
      external wu,dwu
      real*8 start,finish
c
c     constant dimensions
c     direct element matrices:
c     -------------------------
      dimension hxx(8,8),hyy(8,8),hzz(8,8),hxy(8,8),hxz(8,8),hyz(8,8),
     + hvx(8,8),hvy(8,8),hvz(8,8),htmt(8,8),htmf(8,8),
     + rg(8,8),rg2df(4,4)
c
      dimension re(8,8),ro(8,8),inl(8),kb(7),ag(9),hag(9)
      logical leak2,lmass,lheat,lage
c
c     ne = number of elements
c     nn = number of nodes
c     n  = number of degrees of freedom
c     na = total number of non-zero matrix entries in condensed matrix
c        = 14*n (approximately)
c     nw = width of matrix ib ( = 15 for prismatic grid )
c     nf = maximum number of nodes on one face
c     laa = 3*n + na
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     variable dimensions:
c     arrays of size nn
c     -----------------
      real*8 tq(maxnn),pq(maxnn),cq(maxnn)
      real*8  x(maxnn),y(maxnn),z(maxnn),por(maxne),sw(maxne)
      real*8 u0(maxnn),u1(maxnn),fc(maxnn),fb(maxnn),fs(maxnn)
      dimension  f(maxnn),u2(maxnn),fbf(nf)
      integer   ic(maxnn),icflow(maxnn),lc(maxnn),map(maxnn)
c
c     1D line or 2D plane elements - transport 
c     -----------------------------------------
      dimension rel(4,4),rol(4,4),ifracl(maxfrac),inline(4)
      dimension lvert(maxfrac),ifdim(maxfrac)
      real*8 vlin(maxfrac),vx2d(maxfrac),vy2d(maxfrac),xarea(maxfrac)
      dimension hx(4,4),hy(4,4),hh(4,4),ht(4,4),wx(4,4),wy(4,4),hgz2(4)
c
c     arrays of size ne, (ne,8)
c     ------------------------
      real*8   vx(maxne),vy(maxne),vz(maxne)
      real*8   exl(maxne),eyl(maxne),ezl(maxne)
      real*8   vt(maxne),tclm(maxne),cot0(maxne),cot1(maxne)
      real*8   tkl(maxne),cpsm(maxne)
      real*8    depc(maxne),vexp(maxne),decay(maxne),smax(maxne)
      integer  inflow(maxne,8),in(maxne,8)
      real*8   pp(maxne),qq(maxne),ppn(maxnn),qqn(maxnn)
      real*8   spn0(maxnn),spn1(maxnn)
c
c     arrays for internal heat source
c     --------------------------------
      dimension flux(maxfx),fluxv(maxfx),ifl(maxfx),nflux(maxfx,8)
c
c     array of size nf
c     ----------------
      dimension ara(nf)
      dimension bz(nf),bzf(nf)
c
c     arrays a(na),iaa(na),ind(n+1),ib(n,15),aa(3n+na)
c            for conjugate gradient solver only
c     ------------------------------------------------
      dimension a(maxna),aa(laa)
      integer iaa(maxna),ind(maxnn),ib(maxnn,nw)
c
c     function defining -absolute- water density as a function of t
c     for use in latent heat term, thermal parameters
c     --------------------------------------------------------------
c
      nx=nex+1
      ny=ney+1
      nexy = nex*ney

c     initialize arrays,
c     condense non-zero flux array:
c     -------------------------------
      do 45 i=1,nn
      fs(i)=0.d0
      f(i)=0.d0
      if(ic(i).eq.1) goto 45
      ii=i-lc(i)
      f(ii)=fb(i)                !for type2 transport bdy, fb contains heat flux
  45  continue
c
      do 247 i=1,maxna
  247 a(i)=0.d0
c
c     _________________________________________
c
c     matrix generation loop
c     loop over elements:
c     =========================================
c     =========================================
c
c     call clock@(start)
c
c     write(6,903)
c903  format(56x,'l     tavg       por        sw           ww1'
c    +'      wi1          wm         wu1       dwu1      dlh1'
c    +'      Co       tk')

      eps = 0.000001         !increment dT for derivative of co

      do 100 l=1,ne
      ppl=pp(l)
      qql=qq(l)
c
      do 101 i=1,8
  101 inl(i)=in(l,i)
c
c     average elemental temperature:
c     tav0: avg. at start of time step / tav1: avg. at latest iteration
c     -----------------------------------------------------------------
      tav0=0.d0
      tav1=0.d0
      do 103 i=1,8
      tav0=tav0+u0(inl(i))
      tav1=tav1+u1(inl(i))
 103  continue
      tav0=tav0/8.d0
      tav1=tav1/8.d0
      tavg=(tav0+tav1)/2.0
c     tavg=tav2
c
c     fixed retardation = rtrans for mass transport or age only
c     ----------------------------------------------------------
c     if(lmass.or.lage) then

c     cot0(l)=rtrans        !for mass transport, cot0 and cot1 represent retardation, and go on rhs 
c     cot1(l)=rtrans
c     cotl0=rtrans
c     cotl1=rtrans
c     tccx=0.d0
c     tccy=0.d0
c     tccz=0.d0

c     wwavg=1.d0            !check for other terms below for heat that need to be defined for transport
c     cpfavg=1.d0
c     cotlavg = (cotl0+cotl1)/2.d0
c    
c     else 
c
c     identify temp.-dependent thermal parameters for sat/unsat zone:
c     and elemental latent heat contribution term
c     thermal K tccx adjusted for water,ice,solids fractions
c     need to adjust cpi for temperature-dependence
c     cpi=c_i * rho_i = 2000.*916.7
c     cpm=c_m * rho_m = ...
c     sqrt-mean thermal K: see Roy et al. 1981 (in Mottaghy & Rath GJI 2006)
c     Non-linear thermal R ... 
c     -----------------------------------------------------------------------
      porsw = por(l)*sw(l)
      wu0 = wu(tav0,ppl,qql,modelwu,ts)
      wu1 = wu(tav1,ppl,qql,modelwu,ts)
      ww0=porsw * wu0
      ww1=porsw * wu1
      wi0=porsw * (1.d0-wu0)
      wi1=porsw * (1.d0-wu1)
      wm=(1.d0-por(l))  

      den0 = den(tav0,lheat,lmass,lage)
      den1 = den(tav1,lheat,lmass,lage)
      cpf0=cf*den0
      cpf1=cf*den1
      cpe0=ww0*cpf0 + wi0*cpi + wm*cpsm(l)
      cpe1=ww1*cpf1 + wi1*cpi + wm*cpsm(l)

c      thermal conductivity tk trans1
c      -------------------------------
c      either linear implicit ...
c      tkl(l) = ww1*tclw + wi1*tcli + wm*tclm(l)            ! geometric mean    !hardwire Interfrost

c      linear - centre-weighted ....
c      tcl01 =  ww0*tclw + wi0*tcli 
c      tcl02 =  ww1*tclw + wi1*tcli 
c      tkl(l) = ((tcl01 + tcl02)/2.)   + wm*tclm(l)             !centre-weighting of TK
c
c
c      or sqrt-average implicit ...
      tkl(l) = (ww1*sqrt(tclw)+wi1*sqrt(tcli)+wm*sqrt(tclm(l)))**2     !more physically realistic

      tk = tkl(l)

c     dlh=por(l)*den(tav1,lmass,lage)*wlh*dwu(tav1,p,q,modelwu,ts)
      dwu0 = dwu(tav0,ppl,qql,modelwu,ts)
      dwu1 = dwu(tav1,ppl,qql,modelwu,ts)
      dlh0=por(l)*pi*wlh*dwu0            !pi (ice density) 
      dlh1=por(l)*pi*wlh*dwu1
      cot0(l)=(cpe0+dlh0)                !heat capacity 
      cot1(l)=(cpe1+dlh1)                
      cotl0=cot0(l)
      cotl1=cot1(l)
      cotlavg = (cotl0+cotl1)/2.d0
c
c     derivative term T*dCo/dT
c    ------- ------------------
c      co1 = (porsw*wu(tavg,p,q,modelwu,ts)*cf*den(tavg,lmass,lage)) +
c     +       (porsw*(1.-wu(tav0,p,q,modelwu,ts))*cpi) +
c     +       (wm*cpsm(l)) + 
c     +       (por(l)*pi*wlh*dwu(tavg,p,q,modelwu,ts))
c      tavg2 = tavg+eps 
c      co2 = (porsw*wu(tavg2,p,q,modelwu,ts)*cf*den(tavg2,lmass,lage)) +
c     +       (porsw*(1.-wu(tavg2,p,q,modelwu,ts))*cpi) +
c     +       (wm*cpsm(l)) + 
c     +       (por(l)*pi*wlh*dwu(tavg2,p,q,modelwu,ts))
c      tdcodt =  (tavg)*(co2-co1)/eps                      !
c     if(tavg.lt.-0..and.tavg.gt.-3.) then
c            write(6,4898) tavg,co1,co2,tdcodt,cotlavg
c4898         format(5x,'tavg,co1,co2,tdcodt,cotlavg ',5e15.4)
c      endif           
c       cotlavg = cotlavg + tdcodt

c     wwavg = time-centred porosity*Sw*Wu
c     cpfavg = time-centred cfpf
c     tccx   = thermal conductivity lambda
c     November 2013: keep tk on left hand side, keep Co on right hand side
c     uses average ww and cpf for advective term
c     ------------------------------------------------------------------------------------------------------
      wwavg = (ww0+ww1)/2.d0
      cpfavg = (cpf0+cpf1)/2.d0
      tccx = tk                 
      tccy = tccx
      tccz = tccx
c
c hardwire for Lunardini tests - fixed parameters in 3 zones
c ------------------------------------------------------------
c  Lunardini parameters
c     tk = 5.0182e-6                                        !ice
c     if(tavg.ge.-4. .and. tavg .le.0.) tk = 1.6512e-7      !mushy
c     if(tavg.ge.0.) tk=3.5030e-6                           !unfrozen     

c  Sutra parameters
c     tk = 1.3178e-6                                        !ice
c     if(tavg.gt.-4.d0 .and. tavg .lt.0.d0) tk = 4.6515e-7      !mushy
c     if(tavg.ge.0.d0) tk=1.2188e-6                           !unfrozen     

c Lunardini L=0 3 zones
c      tk = 5.0182e-6                                      !ice
c      if(tavg.ge.-4. .and. tavg .le.0.) tk = 4.26e-6      !mushy : lamda/C = 2.94135/690360
c      if(tavg.ge.0.) tk=3.5030e-6                         !unfrozen      

c Lunardini1 (3-zone, no convection): k and Co split; with latent heat included - use this for test
c      tk = 3.464                                 !ice
c      cotl0 = 690360.
c      if(tavg.ge.-4. .and. tavg .le.0.) then     !mushy
c         tk = 2.941                 
c         cotl0 = 1.78e7 
c      endif
c      if(tavg.ge.0.) tk=2.418                    !unfrozen
c
c Lunardini2 (2-zone, + convection): k and Co split; with latent heat included - use this for test
c      tk = 2.57                                 !ice
c      cotl0 = 2.06e6
c      if(tavg.ge.0) then
c       tk=1.839                                  !unfrozen
c       cotl0 = 3.201e6
c      endif
c     --------------------------------------------------------
c     final values used in equation for Lunardini 3-zone hardwire:
c     if setting wwavg=1 and cpfavg=1, must make sure that all dispersion (and molecular diffusion) terms =0
c     activate for all Lunardini solutions
c     ------------------------------------------------------------------------------------------------------
c      cotl1=cotl0
c      cotlavg = (cotl0+cotl1)/2.d0
c      wwavg = 1.d0
c      cpfavg = 1.d0
c      tccx = tk                 
c      tccy = tccx
c      tccz = tccx
c
c     debug print
c     ------------
c     if(l.eq.1 .or. l.eq.5 .or. l.eq.10.or.l.eq.ne)           
c    +write(6,756) l,tavg,por(l),sw(l),ww1,wi1,wm,wu1,dwu1,dlh1,cotl1,tk
c756  format(5x,'l,tavg,por,sw,ww1,wi1,wm,wu1,dwu1,dlh1,Co,tk:',
c    +            i7,11e11.3)

c     if(l.le.3.or. (ne-l).lt.3) then
c     if(l.eq.1) then
c     write(6,9333)
c    + l,tav1,cpe,cpm,dlh,ww,wi,por(l),sw(l),wu(tav1,p,q,modelwu,ts),
c    +                             dwu(tav1,p,q,modelwu,ts),col1,tk,tccx
c9333 format(2x,'l='i6,',tav1=',e10.3,',cpe=',e10.3,',cpm=',e10.3,
c    + ', dlh=',e10.3,', ww=',e10.3,', wi=',e10.3,
c    + ', por=',e10.3,', sw=',e10.3,', wu=',e10.3,', dwu=',e10.3,
c    + ', Co1=',e10.3,', tk=',e10.3', tccx=',e10.3)
c     endif

c     end hardwire
c --------------------------------------------------------------

c     endif                     !for heat transport

c     dispersion tensor
c     (retardation term is entered during elemental matrix assembly)
c     Dij * wwavg*cpfavg (not divided out now)
c     --------------------------------------------------------------
      vxl=vx(l)
      vyl=vy(l)
      vzl=vz(l)
      vx2=vxl*vxl
      vy2=vyl*vyl
      vz2=vzl*vzl
      vxy=vxl*vyl
      vxz=vxl*vzl
      vyz=vyl*vzl
      v= sqrt(vx2+vy2+vz2)
c
c     BF or Lichtner dispersion formulation:
c     ---------------------------------------
      if (kdisp.eq.0) then 
      al=alh
c  original B-F version
      dxx=((al*vx2+ath*vy2+atv*vz2)/v + dd) + tccx 
      dyy=((ath*vx2+al*vy2+atv*vz2)/v + dd) + tccy
      dzz=((atv*vx2+atv*vy2+al*vz2)/v + dd) + tccz
      dxy=(al-ath)*vxl*vyl/v
      dxz=(al-atv)*vxl*vzl/v
      dyz=(al-atv)*vyl*vzl/v
      
      else
c
c  Lichtner version (Lichtner et al, 2002)
c  added March 12, 2013, EOF
c
      costh=vzl/v
      cos2th=costh*costh
      vxy2=vx2+vy2
      all=alh+cos2th*(alv-alh)
      att=atv+cos2th*(ath-atv)

      dxx=((all*vx2+ath*vy2*(1.0+vz2/vxy2)+att*vz2*vx2/vxy2)/v+dd)+tccx
      dyy=((ath*vx2*(1.0+vz2/vxy2)+all*vy2+att*vz2*vy2/vxy2)/v+dd)+tccy
      dzz=((att*vxy2+all*vz2)/v+dd)+tccz
      dxy=(all-ath*(1.0+vz2/vxy2)+att*vz2/vxy2)*vxl*vyl/v
      dxz=(all-att)*vxl*vzl/v
      dyz=(all-att)*vyl*vzl/v

      endif

      wucpf = wwavg*cpfavg

c     dxx=((al*vx2+ath*vy2+atv*vz2)/v + dd)  + tccx
c     dyy=((ath*vx2+al*vy2+atv*vz2)/v + dd)  + tccy
c     dzz=((atv*vx2+atv*vy2+al*vz2)/v + dd)  + tccz
c     dxy=(al-ath)*vxl*vyl/v  
c     dxz=(al-atv)*vxl*vzl/v 
c     dyz=(al-atv)*vyl*vzl/v  

c     hardwire
c      if(dxy.lt.1.0e-5) write(6,8837) dxy
c      if(dxz.lt.1.0e-5) write(6,8838) dxy
c      if(dyz.lt.1.0e-5) write(6,8839) dxy
c 8837 format(10x,'dxy: ',e15.4)
c 8838 format(10x,'dxz: ',e15.4)
c 8839 format(10x,'dyz: ',e15.4)
c      call flush(6)
c
c     elemental matrix, scheme 2
c     --------------------------
      rex=wp*dxx+wa*vx2*dt2  * wucpf          !wucpf ok - tested with craflush and oneD and lunardini July10
      rey=wp*dyy+wa*vy2*dt2  * wucpf            !check - not needed in transv2 ???
      rez=wp*dzz+wa*vz2*dt2  * wucpf
      rxy=wp*dxy+wa*vxy*dt2  * wucpf
      rxz=wp*dxz+wa*vxz*dt2  * wucpf
      ryz=wp*dyz+wa*vyz*dt2  * wucpf
      pex=wp1*dxx+wa1*vx2*dt2  * wucpf
      pey=wp1*dyy+wa1*vy2*dt2  * wucpf
      pez=wp1*dzz+wa1*vz2*dt2  * wucpf
      pxy=wp1*dxy+wa1*vxy*dt2  * wucpf
      pxz=wp1*dxz+wa1*vxz*dt2  * wucpf
      pyz=wp1*dyz+wa1*vyz*dt2  * wucpf
      pvx=vxl        * wucpf                    !new
      pvy=vyl        * wucpf                    !new
      pvz=vzl        * wucpf                    !new
  155 continue
c
c
c     numerical(kint=1) or direct(kint=0) integration ... heat transport1
c     returns {re},{ro}
c     ===================================================================
c     if (l.eq.1) write (6,712) rex,rey,rez,rxy,rxz,ryz,
c    1pex,pey,pez,pxy,pxz,pyz,pvx,pvy,pvz
c 712 format (9f14.7)
c
      if(kint.eq.1)
     + call gquad1(rex,rey,rez,rxy,rxz,ryz,pex,pey,pez,
     +      pxy,pxz,pyz,pvx,pvy,pvz,dt1,x,y,z,inl,re,ro,l,
     +      maxnn,cotl0,cotl1,cotlavg,ag,hag,tk)
c
       if(kint.eq.0)
     + call exint1(rex,rey,rez,rxy,rxz,ryz,pex,pey,pez,
     +   pxy,pxz,pyz,pvx,pvy,pvz,hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,
     +   dt1,x,y,z,inl,re,ro,maxnn,htmt,htmf,cotl0,cotl1,cotlavg,
     +   agefx,rg,exl,eyl,ezl,maxne,l,tk,vx,vy,vz)
c
c      if((l.eq.8218.or.l.eq.10330).and.(tdays.eq.6..or.tdays.eq.20.))
c     + then
c      write(6,281) tdays
c  281 format(/10x,'transport matrix check at tdays = ',f10.3)
c      write(6,181) l,(i,(re(i,j),j=1,8),i=1,8)
c      write(6,181) l,(i,(ro(i,j),j=1,8),i=1,8)
c 181  format(/,' transport matrix check - re,ro: ',i6,/(i2,8(e11.3,1x)))
c      endif
c
c     assembly for conjugate gradient solver
c     ======================================
      do 81 i=1,8
      ki = inl(i)
      if (ic(ki).eq.1) go to 81
      ii = ki-lc(ki)
      do 82 j=1,8
      kj = inl(j)
      if (ic(kj).eq.1) go to 85
      jj=kj-lc(kj)
      if (ii.gt.jj) go to 86
      do 83 k=1,nw
      if (ii.ne.ib(jj,k)) go to 83
      kk=k
      go to 84
   83 continue
   84 continue
      iv=ind(jj)+kk-1
c
c     write (6,738) l,i,j,ii,jj,ib(jj,1),ib(jj,kk),ind(jj),kk,iv,
c    1a(iv),se(i,j)
c 738 format (10i5,2f14.8)
c     -------------------------
c     left side coefficient vector
c     ----------------------------
      a(iv)=a(iv)+re(i,j)
c
c     right side contribution
c     -----------------------
   86 continue
      fs(ii)=fs(ii) - ro(i,j)*u0(kj) + rg(i,j)
      go to 82
c
c     linking term:
c     -------------
   85 continue
      f(ii)=f(ii)-(re(i,j)+ro(i,j))*u0(kj)
   82 continue
   81 continue
c
c     end of element loop
c     ===================
  100 continue
c
c     1D line and 2D plane element assembly - transport
c     Note: fracture thermal conductivity assumed =  0.5 J/msC (water)
c     hardwire porfrac
c     ----------------------------------------------------------------
c     porfrac=1.
c     vxmax=-999.
c     vymax=-999.
c     vxmin=+999.
c     vymin=+999.
c     do kf=1,nfrac
c     if(ifdim(kf).eq.2) then
c     vxmax=max(vx2d(kf),vxmax)
c     vymax=max(vy2d(kf),vymax)
c     vxmin=min(vx2d(kf),vxmin)
c     vymin=min(vy2d(kf),vymin)
c     endif
c     enddo
c     write(6,7822) vxmax,vymax,vxmin,vymin
c7822 format(/10x,'fractures: vxmax,vymax,vxmin,vymin: ',4e15.4)

c     all fracture elements:
c     -----------------------
      do 555 kf=1,nfrac
c       
c     1D fractures - transport
c     ------------------------
      if(ifdim(kf).eq.1) then
      vl = vlin(kf)            
      l =  ifracl(kf)
c     cotl0 = cot0(l)
c     cotl1 = cot1(l)
c
      call frac_line(maxfrac,maxne,lvert,inline,inflow,kf,l,fracdim,
     +              exl,eyl,ezl)
c
c     update tcc for frozen state -
c     --------------------------------------
      tb = 0.25d0
     +       * (u1(map(inline(1)))+u1(map(inline(2)))
     +         +u0(map(inline(1)))+u0(map(inline(2))))
c      tcc = tkl(l)/(cf*den(tb,lmass,lage))
c     tcc = 0.6d0/(cf*den(tb,lmass,lage))
      wufrac = wu(tb,p2,q2,modelwu2,ts2)
      cpfavg = (cf*den(tb,lheat,lmass,lage))
      dwufrac = dwu(tb,p2,q2,modelwu2,ts2)
c     dlhfrac = por(1)*pi*wlh*dwufrac            !pi (ice density)       
      dlhfrac =        pi*wlh*dwufrac            !pi (ice density)    por = 1. in fracs !
      cotl = wufrac*cpfavg + (1.d0-wufrac)*cpi +dlhfrac            !Co+L
      tcc = (wufrac*sqrt(0.6d0) + (1.d0-wufrac)*sqrt(tcli))**2     !tk
c     if(lmass.or.lage) then            !use else-if to skip above for mass transport
c         tcc=0.d0
c         cotl = rtrans
c         wufrac=1.d0
c         cpfavg=1.d0
c     endif
      dxx = (alh*vl + dd) * (wufrac*cpfavg) + tcc              !tcc on lhs
      sqa = xarea(kf)/por(1)              !check por
      rezl = sqa  * ( wp*dxx  +  wa*vl*vl*dt2 ) / fracdim  
      pezl = sqa  * ( wp1*dxx + wa1*vl*vl*dt2 ) / fracdim 
      pvzl = sqa  *  (vl/2.d0)  * (wufrac*cpfavg) 
c
c     lumped -     use htl = xarea * fracture_length*dt1/2.
c     consistent - use htl = xarea * fracture_length*dt1/6.
c     check decay rdec (March2007): /fracdim or *fracdim
c  -----------------------------------------------------------------------------
      htl =  cotl* (sqa  *  fracdim * dt1/2.d0)   !cotl (Co) on rhs, do not divide by R from blocks
c     rdec =  sqa  *  fracdim * decay(1) 
      rdec =  0.0d0                          !nodecay for heat transport 
c
      rel(1,1) = rezl + htl + rdec
      rel(1,2) = -rezl + htl + rdec
      rel(2,1) = rel(1,2)
      rel(2,2) = rel(1,1)
      rol(1,1) = (+pezl - pvzl) - htl
      rol(1,2) = (-pezl + pvzl) 
      rol(2,1) = (-pezl - pvzl)
      rol(2,2) = (+pezl + pvzl) - htl
c
      do 2012 i=1,2
      ki = map(inline(i))
      if(ic(ki).eq.1) goto 2012
      ii= ki-lc(ki)
        do 2013 j=1,2
        kj = map(inline(j))
        if(ic(kj).eq.1) goto 2014
        jj = kj - lc(kj)
        if(ii.gt.jj) goto 2015  
c
      do 2083 k=1,nw
      if (ii.ne.ib(jj,k)) go to 2083
      kk=k
      go to 2084
 2083 continue
 2084 continue
c
      iv=ind(jj)+kk-1
      a(iv)= a(iv) + rel(i,j) 
 2015 continue
      fs(ii) = fs(ii) - rol(i,j)*u0(kj) 
      goto 2013
 2014 continue
      f(ii) = f(ii) - (rel(i,j)+rol(i,j))*u0(kj)
c
 2013 continue
 2012 continue
 2011 continue
      endif
c       
c     2D fractures - transport
c     ------------------------
c     ----------------------------------------------------------------
      if(ifdim(kf).eq.2) then
c
c     fracture velocities: 
c     use 2 velocities: vxl,vyl where x,y could also be x,z or y,z
c     --------------------------------------------------------------

c hardwire v
c     vx2d(kf)=1.16e-3
c     vy2d(kf)=1.0e-15


      vx2 = vx2d(kf)**2
      vy2 = vy2d(kf)**2
      vxy = vx2d(kf)*vy2d(kf)
      vv = sqrt(vx2+vy2)
      l =  ifracl(kf)
c
c      if(kf.eq.1.or.kf.eq.nfrac)
c    + write(6,9922) l,vx2d(kf),vy2d(kf),vv,xarea(kf),dd
c9922 format(/10x,'2D fracture transport l= ',i5,
c    +       /10x,'vx2 = ',e12.4,' vy2 = ',e12.4,' vv= ',e12.4,
c    +       /10x,'xarea ... ',e12.5,' dd ... ',e12.5)
c
      call frac_plane(maxfrac,maxne,lvert,inline,inflow,kf,l,
     +                fdimx,fdimy,exl,eyl,ezl,rexel,sarea)
c
      mapinl1 = map(inline(1))
      mapinl2 = map(inline(2))
      mapinl3 = map(inline(3))
      mapinl4 = map(inline(4))
      tb = 0.125 * ( u1(mapinl1)+u1(mapinl2)+u1(mapinl3)+u1(mapinl4) 
     +              +u0(mapinl1)+u0(mapinl2)+u0(mapinl3)+u0(mapinl4) ) 

c     hardwire thermal diffusivity of fracture
c     divide by por (from porous matrix equation since we use v not q)
c     check: must assume uniform porosity, correct for water content
c     ----------------------------------------------------------------
c     tcc = 0.6d0/(cf*den(tb,lmass,lage))
      wufrac = wu(tb,p2,q2,modelwu2,ts2)      !por ?
      cpfavg = cf*den(tb,lheat,lmass,lage)     
      dwufrac = dwu(tb,p2,q2,modelwu,ts2)      !bug fixed p2,q2,ts2
c     dlhfrac = por(1)*pi*wlh*dwufrac            !pi (ice density) check por 
      dlhfrac =        pi*wlh*dwufrac            !pi (ice density) check por Jan 2016 por removed (por=1.0) 
      cotl = wufrac*cpfavg + (1.d0-wufrac)*cpi + dlhfrac        !Co+latent heat, no solids in fracture   hardwire - por

      tcc = (wufrac*sqrt(0.6d0) + (1.d0-wufrac)*sqrt(tcli))**2     !tk
c     cotl =cotl/por(1)

c     if(lmass.or.lage) then            !use else-if to skip above for mass transport
c         tcc=0.d0
c         cotl = rtrans
c         wufrac=1.d0
c         cpfavg=1.d0
c     endif      

c hardwire for sfrac-h no conduction in fracture
      tcc=0.
c
      sqa = xarea(kf)/por(1)       !por ok
c     if(kf.eq.1.or.kf.eq.nfrac) 
c    +   write(6,*) 'kf,wufrac,cpfavg,por(1): ',kf,wufrac,cpfavg,por(1)
      wucpf = wufrac*cpfavg
c
c     in 2D fractures, always use al=alh and alt=ath
c     if(lvert(kf).eq.3 .or. lvert(kf).eq.4) alphat = atv
      alphat = ath
c
c     check later ... use alh for horizontal frcs, alv for vertical ?
      dxx = ((alh*vx2 + alphat*vy2)/vv + dd)  + tcc           !check later: should not have tcc*wucp
      dyy = ((alphat*vx2 + alh*vy2)/vv + dd)  + tcc
      dxy = ((alh-alphat)*vx2d(kf)*vy2d(kf)/vv)  
c
c     if(kf.eq.1.or.kf.eq.nfrac)
c    +  write(6,7722) kf,xarea(kf),rexel,fdimx,fdimy,dt1,dt2
c7722 format(/10x,'fracture plane # ... ',i7,
c    +       /10x,'xarea: ',e12.5,' rexel: ',e12.5,
c    +       /10x,'fdimx: ',e12.5,' fdimy: ',e12.5,
c    +       /10x,'dt1,dt2: ',2e15.5)
      rex = sqa * ( wp*dxx  +  wa*vx2*dt2 ) * (rexel/6.d0)  * wucpf
      rey = sqa * ( wp*dyy  +  wa*vy2*dt2 ) / (rexel*6.d0)  * wucpf
      rxy = sqa * ((wp*dxy  +  wa*vxy*dt2 ) / 2.d0)         * wucpf
      pex = sqa * ( wp1*dxx + wa1*vx2*dt2 ) * (rexel/6.d0)  * wucpf
      pey = sqa * (( wp1*dyy + wa1*vy2*dt2 ) / (rexel*6.d0))* wucpf
      pxy = sqa * (( wp1*dxy + wa1*vxy*dt2) / 2.d0)         * wucpf
      pvx = sqa * (vx2d(kf)*fdimy/12.d0)                    * wucpf
      pvy = sqa * (vy2d(kf)*fdimx/12.d0)                    * wucpf
c
c     lumped - 
c     consistent - 
c     --------------------------------------------
      ret0 = cotl*(sqa * sarea * dt1/36.d0)            !/por for cotl already in sqa
      ret1 = ret0
c     rdec = sqa * sarea * decay(1) /36.d0
      rdec =  0.0d0                                     !no decay for heat
      rage = sqa * agefx * sarea / 36.d0
c
c     write(6,8888) l,dxx,dyy,dxy,rex,rey,rxy,ret1
c8888 format(/10x,'2D fracture - transport, l= ',i5,
c    +       /10x,' dxx,dyy,dxy: ',3e12.3,
c    +       /10x,' rex,rey,rxy: ',3e12.3,' ret1: ',e12.3)
c
      do 56 i=1,4
      do 56 j=1,4
      rd = rdec*ht(i,j)
      rel(i,j) = (rex*hx(i,j) + rey*hy(i,j) + rxy*hh(i,j))
     +                                      + ret1*ht(i,j) + rd
      rol(i,j) = (pex*hx(i,j) + pey*hy(i,j) + pxy*hh(i,j)
     +           + pvx*wx(i,j) + pvy*wy(i,j))  - ret0*ht(i,j)
      rg2df(i,j) = rage * ht(i,j)
  56  continue
c
      do 3012 i=1,4
      ki = map(inline(i))
      if(ic(ki).eq.1) goto 3012
      ii= ki-lc(ki)
        do 3013 j=1,4
        kj = map(inline(j))
        if(ic(kj).eq.1) goto 3014
        jj = kj - lc(kj)
        if(ii.gt.jj) goto 3015  
c
      do 3083 k=1,nw
      if (ii.ne.ib(jj,k)) go to 3083
      kk=k
      go to 3084
 3083 continue
 3084 continue
c
      iv=ind(jj)+kk-1
      a(iv)= a(iv) + rel(i,j) 
 3015 continue
      fs(ii) = fs(ii) - rol(i,j)*u0(kj) + rg2df(i,j)
      goto 3013
 3014 continue
      f(ii) = f(ii) - (rel(i,j)+rol(i,j))*u0(kj)
c
 3013 continue
 3012 continue
 3011 continue
c
      endif
c
c     get next 2D transport fracture
c     ------------------------------
 555  continue
c
c     contribution from thermal leakage boundary:
c     new: use fc2 (fc here) to store variable surface temp vtemp
c     check with new formulation Nov 2013
c     -------------------------------------------------------
      if(leak2) then
      k=0

      do 102 i=nz,nn,nz
      tcsbz=tcs/bz(i/nz)

c     skip if 1st type temp node
      if(ic(i).eq.1) goto 102

c     xxx = x(i)
      i2=i/nz
      i3 = mod(i2,ny)
      if(i3.eq.0) i3=ny-1
      i4 = i2/ny +1
      if(mod(i2,ny).eq.0) i4=i4-1
      if(i4.eq.nx) i4=i4-1
      ielem = (nez-1)*nexy+(i3-1)*nex+i4     !need element # for porosity

      if(ielem.gt.ne.or.ielem.lt.1) then
      write(6,3984)
 3984 format(/10x,'error in assigning element# in trans ...',
     +       /10x,'stopping ... ')
      stop
      endif

      tavg=(u1(i)+u0(i))/2.
      wuleak = wu(tavg,p,q,modelwu,ts)
      ww=porsurf*sat*wuleak
      wi=porsurf*sat*(1.0-wuleak)
      wm=(1.0-porsurf) 
      dwuleak = dwu(tavg,ppn(i),qqn(i),modelwu,ts)
      dlhleak=por(ielem)*pi*wlh*dwuleak            !pi (ice density), don't use porsurf ?

      cpf=cf*den(tavg,lheat,lmass,lage)
c     cpe=ww*cpf + wi*cpi + wm*cpm
      cpeu = ww*cpf + wi*cpi + (1.-porsurf)*cpsm(ielem) + dlhleak 
c     p1cpm=(1.-porsurf)*cpsm(ielem)

c     rtemp=fc(i)
      rtemp = max(0.,ttemp+rtdiff)          !recharge water temperature rtemp, min=0
      cfd=cf*den(rtemp,lheat,lmass,lage)          !
      k=k+1
      ii=i-lc(i)

c     bdyflux = fbf(k)
c     if(icflow(i).eq.1) bdyflux = -vz(ielem)*por(ielem)     !for type1 flow node, use v*por as flux
      bdyflux = -vz(ielem)*por(ielem)     !always use v*por as flux
      if(bdyflux.lt.0.) bdyflux=0.

c     if(icflow(i).eq.1) then
c             write(6,4848) i,ielem,fbf(k),bdyflux
c4848          format(10x,'i,ielem,fbf(k),bdyflux: ',2i6,2e12.3)
c       endif


c     write(6,5676) k,ielem,fbf(k),ara(k)
c5676 format(i5,'leak2: processing element # ',i7,',fbf,ara= ',2e10.3)
c     cpeu=por(ielem)*sat*cf*den(tavg,lmass) + p1cpm
c     cpe=ww*cpf + wi*cpi + wm*cpm
c=====
c     check with new formulation Co on rhs, ww*cpf on lhs?
c     check ... lhs units are m^2*(J/m/s/C)/m * C = J/s = W
c     --------------------------------------------------------------
       arak=ara(k)
       a(ind(ii)) = a(ind(ii)) + (arak*(tcsbz + bdyflux*cfd))      !new for recharge, but then removed from lhs, ok.
       f(ii)=f(ii) + arak*(fc(i)*tcsbz + rtemp*bdyflux*cfd + bzf(i/nz))  ! added ara(k)*bzf from bzin.data

c      write(6,4949) 
c    +  k,ara(k),fc(i),tcsbz,fc(i)*tcsbz,rtemp*bdyflux*cfd,bzf(i/nz),
c    +  ii,f(ii)
c4949  format('k,arak,fc(i),tcsbz,fc(i)*tcsbz,rtemp*bdyflux*cfd,',
c    +        'bzf(i/nz),ii,f(ii):' ,i5,6e12.4,i8,e15.5)

c  ok to march 2015, arak not factored, too slow ...
c      a(ind(ii)) = a(ind(ii)) + (arak*tcsbz) 
c     +                        + (arak*bdyflux*cfd)      !new for recharge, but then removed from lhs, ok.
c      f(ii)=f(ii) + (arak*fc(i))*tcsbz + (arak*rtemp*bdyflux*cfd)
c     +            + (arak*bzf(i/nz))
c original:
c      a(ind(ii)) = a(ind(ii)) + ara(k)*tcs/(bz*cpeu)
c      f(ii)=f(ii)+ (ara(k)*rtemp/cpeu)*(tcsbz+fbf(k)*cfd)      
  102 continue
c=====
      endif

c      write(6,4747) (k,ara(k),bzf(k),k=1,nn/nz)
c4747  format('  k,ara(k),bzf(k) ',i5,2e15.5)      
c
c      write (6,695) (i,fs(i),fc(i),i=nz,nn,nz)
c 695  format(/10x,'i,fs(i),fc(i)...',/(2(i6,2e12.3)))
c
c     source term contribution
c     need to multiply pqq by wwavg*cpf as in 3D blocks
c     (pqq*tq) units are (m^3/s)*(J/m^3.C)*C = J/s
c     --------------------------------------------------
      do 556 i=1,nn
      ii=i-lc(i)
      if(pq(i).le.0.) goto 557
      tavg=(u1(i)+u0(i))/2.      
      wwavg = por(1)*sw(1)*wu(tavg,ppn(i),qqn(i),modelwu,ts)     !check - need element porosity, Sw
c     write(6,34857) 
c    +    por(1),sw(1),wu(tavg,p,q,modelwu,ts),cf,den(tavg,lmass,lage)
34857 format('por(1),sw(1),wu,cf,den: ',5e15.5)

      cpfavg = (cf*den(tavg,lheat,lmass,lage))            
c     cpfavg = cf            
c     if(lmass.or.lage) then
c              wwavg = 1.d0
c              cpfavg=1.d0
c     endif
      pqq= pq(i)*wwavg*cpfavg/(por(1))       !check - need wwavg*cpfavg  - 
c     pqq= pq(i)*wwavg*cpfavg       !check - need wwavg*cpfavg  - 
      a(ind(ii)) = a(ind(ii)) + pqq
      f(ii)=f(ii)+pqq*tq(i)
  557 continue
  556 continue
c
c     internal heat flux elements:
c     Feb. 2005: assumes unit element column, or 4 columns 
c     hardwire for flux across outside faces only. 
c     - distribute over 4 nodes on each of the two outside faces, = (total area/8 nodes)
c     use elemental cpe calculated from stored cot0(l) and cot1(l) arrays
c     needs to be updated for freezing
c     fluxv (W/m3) applied to element volume (1/8 to each node)
c     -----------------------------------------------------------------------------------
      if(numis.gt.0) then
      do 667 i=1,numis
      l = ifl(i)
c     cpe = (cot0(l) + cot0(l))/2.
      
c     wufrac = wu(tb,p,q,modelwu,ts)
c     cpfavg = (cf*den(tb,lmass,lage))
c     cotl = wufrac*cpfavg + (1.d0-wufrac)*cpi                     !Co

c     original
c     sarea = 0.25*(exl(l)*ezl(l) + exl(l)*eyl(l) + eyl(l)*ezl(l))
c     sarea = (2./8.)*(exl(l)*ezl(l) + eyl(l)*ezl(l))
      sarea = (exl(l)*ezl(l) + eyl(l)*ezl(l))
      if(nflux(i,3).eq.in(l,3)) sarea = sarea*2.    !single column: area*2

      sarea = sarea/8.d0                                !distribute over 8 nodes
c     cpe=por(l)*cpf + (1.d0-por(l))*cpm    !need to adjust for freezing
c     if(l.gt.nef) cpe=por(l)*sat*cpf + (1.d0-por(l))*cpm
      do 668 j=1,8
      node = nflux(i,j)
      ii = node - lc(node)
c      f(ii) = f(ii) + (flux(i)*sarea/cpe) * exp(-decay2*tdays)  !hardwire decay off
c     f(ii) = f(ii) + (flux(i)*sarea/cpe)
      f(ii) = f(ii) + flux(i)*sarea + fluxv(i)*exl(l)*eyl(l)*ezl(l)/8.d0       !W/m2*area + W/m3*volume

c     write(6,8475) l,ii,flux(i),sarea,cpe
c8475 format(10x,'debug: l,ii,flux,sarea,cpe: ',2i10,3e15.5) 

 668  continue
 667  continue

      endif
c
c     assemble flux vector
c     --------------------
      do 121 i=1,n
      f(i)= f(i) + fs(i) 
  121 continue
c
c     if (it.gt.1) go to 735
c     do 733 i=1,maxna
c 733 write (6,734) i,iaa(i),a(i)
c 734 format (2i5,f14.8)
c 735 continue
c
c     conjugate gradient solution - transport1
c     ========================================
      neu=1
c     call clock@(start)
      call precg (a,n,iaa,ind,f,u2,aa,laa,neu)
c      write(*,123)
c  123 format(10x,'transport pcg solution complete ...')
c     call clock@(finish)
c     tpcgt=tpcgt+finish-start
c
   90 continue
c
c     expand solution vector:
c     -----------------------
      do 141 ii=1,nn
      i = nn-ii+1
      if(ic(i).eq.1) goto 142
      k=i-lc(i)
c      if(lmass .and. u2(k).lt.1.d-50) u2(k)=0.d0  !hardwire check for transport
      u2(i)=u2(k)
      goto 141
  142 u2(i)=fc(i)
  141 continue
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      return
      end
c
c
c ******************************************************************
c ******************************************************************
c
      subroutine gquad1(rex,rey,rez,rxy,rxz,ryz,pex,pey,
     +    pez,pxy,pxz,pyz,pvx,pvy,pvz,dt1,x,y,z,inl,re,ro,l,
     +    maxnn,rtd0,rtd1,cotlavg,ag,hag,tk)
c
c     rtd = thermal retardation term for mass transport, or heat capacity for heat transport
c     generation of isoparametric element matrices
c     full gauss quadrature numerical integration
c     Nov2013: not coded yet for non-linear R (i.e. freeze-thaw)
c     ----------------------------------------------------------
c
      implicit real*8(a-h,o-z)
c
      dimension detj(27),f(8),dgx(8),dgy(8),dgz(8)
      real*8    x(maxnn),y(maxnn),z(maxnn)
      dimension ag(9),hag(9),re(8,8),ro(8,8),inl(8)
      dimension ff(8,27),dx(8,27),dy(8,27),dz(8,27)
c
      m=8
c     number of gauss points
c     ----------------------
c      np=2
c      npx=np
c      npy=np
c      npz=np
c      np2=npx*npy*npz
c      ns=0
c      if (np.eq.3) ns=2
c
c     placement of gauss points
c     -------------------------
      do 310  k=1,2
      zi=ag(k)
      hz=hag(k)
      do 310 j=1,2
      yi=ag(j)
      hy=hag(j)
      do 310 i=1,2
      xi=ag(i)
      hx=hag(i)
      kl=(k-1)*4+(j-1)*2+i
c     if (l.eq.1) write (6,741) i,j,k,kl,xi,yi,zi,hx,hy,hz
c 741 format (/' gauss point',4i5,6f14.7)
c
c     basis functions and derivatives
c     ===============================
      call sp3lin (x,y,z,inl,m,xi,yi,zi,det,l,dgx,dgy,dgz,f,maxnn)
      do 300 jj=1,m
      ff(jj,kl)=f(jj)
      dx(jj,kl)=dgx(jj)
      dy(jj,kl)=dgy(jj)
  300 dz(jj,kl)=dgz(jj)
      detj(kl)=det*hz*hy*hx
  310 continue
c     if (l.ne.1) go to 311
c     write (6,314) l
c 314 format (/' basis functions element',i5)
c     write (6,313) (k,(dx(i,k),i=1,8),k=1,8)
c     write (6,313) (k,(dy(i,k),i=1,8),k=1,8)
c     write (6,313) (k,(dz(i,k),i=1,8),k=1,8)
c 311 continue
c 313 format (i5,8f14.8)
c
c     elemental matrices, scheme 2
c     ----------------------------
      do 51 i=1,m
      do 51 j=1,m
      res=0.
      ros=0.
      do 52 k=1,8
      rd=rex*dx(i,k)*dx(j,k)+rey*dy(i,k)*dy(j,k)+rez*dz(i,k)*dz(j,k)
     1  +rxy*(dx(i,k)*dy(j,k)+dy(i,k)*dx(j,k))
     2  +rxz*(dx(i,k)*dz(j,k)+dz(i,k)*dx(j,k))
     3  +ryz*(dy(i,k)*dz(j,k)+dz(i,k)*dy(j,k))
      rp=pex*dx(i,k)*dx(j,k)+pey*dy(i,k)*dy(j,k)+pez*dz(i,k)*dz(j,k)
     1  +pxy*(dx(i,k)*dy(j,k)+dy(i,k)*dx(j,k))
     2  +pxz*(dx(i,k)*dz(j,k)+dz(i,k)*dx(j,k))
     3  +pyz*(dy(i,k)*dz(j,k)+dz(i,k)*dy(j,k))
      pv=(pvx*dx(j,k)+pvy*dy(j,k)+pvz*dz(j,k))*ff(i,k)
c
c     consistent formulation used here
c     Nov2013: cotlavg is on rhs, as R for mass transport or Co for heat
c     -------------------------------------------------------------------
      rt=cotlavg*dt1*ff(i,k)*ff(j,k)
      res=res+( rd + rt)*detj(k)
      ros=ros+( (rp+pv) - rt)*detj(k)

c      pre Nov2013
c      rt=dt1*ff(i,k)*ff(j,k)
c      res=res+( rd/rtd1 + rt)*detj(k)
c      ros=ros+( (rp+pv)/rtd1 - rt)*detj(k)
   52 continue
c
      re(i,j)=res
      ro(i,j)=ros
   51 continue
  210 continue
c
      return
      end
c
c *********************************************************************
c
      subroutine exint1(rex,rey,rez,rxy,rxz,ryz,pex,pey,pez,
     +     pxy,pxz,pyz,pvx,pvy,pvz,hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,
     +     dt1,x,y,z,inl,re,ro,maxnn,htmt,htmf,rtd0,rtd1,cotlavg,
     +     agefx,rg,exl,eyl,ezl,maxne,l,tk,vx,vy,vz)
c
c     rtd = thermal retardation term
c     assembly of element matrices - exact integration
c     rectangular or orthogonal quadrilateral elements
c     --------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension hxx(8,8),hyy(8,8),hzz(8,8),hxy(8,8),hxz(8,8),hyz(8,8),
     + hvx(8,8),hvy(8,8),hvz(8,8),htmt(8,8),htmf(8,8)
      real*8   x(maxnn),y(maxnn),z(maxnn)
      real*8  exl(maxne),eyl(maxne),ezl(maxne)
      real*8    vx(maxne),vy(maxne),vz(maxne)
      dimension ro(8,8),re(8,8),rg(8,8)
      dimension inl(8)
c
      exll=exl(l)
      eyll=eyl(l)
      ezll=ezl(l)
      rex=rex*eyll*ezll/exll/36.
      rey=rey*exll*ezll/eyll/36.
      rez=rez*exll*eyll/ezll/36.
      rxy=rxy*ezll/12.
      rxz=rxz*eyll/12.
      ryz=ryz*exll/12.
      pex=pex*eyll*ezll/exll/36.
      pey=pey*exll*ezll/eyll/36.
      pez=pez*exll*eyll/ezll/36.
      pxy=pxy*ezll/12.
      pxz=pxz*eyll/12.
      pyz=pyz*exll/12.
      pvx=pvx*eyll*ezll/72.
      pvy=pvy*exll*ezll/72.

      pvz=pvz*exll*eyll/72.
c     rdec=decay(l)*exll*eyll*ezll/216.
      rtt=cotlavg* dt1*exll*eyll*ezll/216.             !split: R or Co (cotlavg) on rhs 
c     rttl=rtt+ rtd1* dt1*exll*eyll*ezll/216.      !new ; doesn't work ...
c     rttr=rtt+ rtd0* dt1*exll*eyll*ezll/216.      !new
c     rage=0.0d0
c
c     lumped: use htmf
c     consistent: use htmt
c     ------------------------
      do 156 i=1,8
      do 156 j=1,8
      re(i,j)=( rex*hxx(i,j)+rey*hyy(i,j)+rez*hzz(i,j)
     1       +rxy*hxy(i,j)+rxz*hxz(i,j)+ryz*hyz(i,j))
     2       +rtt*htmt(i,j)                                  !rdec removed 
      ro(i,j)=( pex*hxx(i,j)+pey*hyy(i,j)+pez*hzz(i,j)
     1       +pxy*hxy(i,j)+pxz*hxz(i,j)+pyz*hyz(i,j)
     2       +pvx*hvx(i,j)+pvy*hvy(i,j)+pvz*hvz(i,j))
     3       -rtt*htmt(i,j)
      rg(i,j)=0.0d0                                       !rage(agefx)removed
  156 continue
c
      return     
      end
c **************************************************************
c =======================================
c
c subpak.for
c
c this subroutine package contains the following routines
c for access from 'heat.for' ...
c -------------------------------------------------------
c
c  1  - prism   ... finite element grid generator
c  2  - sp3lin  ... generation of basis functions and derivatives
c  3  - mindex  ... generation of boundary arrays and condensation codes
c  4  - surf    ... nodal influence areas
c  5  - veloc   ... calculation of elemental velocities:
c                   veloc1 - numerical derivatives
c                   veloc2 - direct derivatives
c  6  - deform  ... grid deformation for watertable mounding
c  7  - volume  ... element volume
c  8  - cpu     ... elapsed cpu execution time
c  9  - second  ... time conversion to h:m:s
c 10  - moment  ... calculates first and second moments
c 11  - wu      ... external function defining unfrozen moisture content
c                            as a function of temperature.
c =====================================================================
c
      subroutine prism (x,y,z,in,xl,yl,zl,nx,ny,nz,nex,ney,nez,nn,ne,
     +       maxnn,maxne,xlim,ylim,zlim,nlx,nly,nlz,ngx,ngy,ngz,
     +       mxgx,mxgy,mxgz,exl,eyl,ezl,map,mpa,inb,inbl,inbr,
     +       kcall,nf,lunsat,kgo)
c
      implicit real*8(a-h,o-z)
c
      real*8  x(maxnn),y(maxnn),z(maxnn)
      real*8  exl(maxne),eyl(maxne),ezl(maxne)
      dimension xlim(mxgx),ylim(mxgy),zlim(mxgz)
      dimension nlx(mxgx),nly(mxgy),nlz(mxgz)
      integer   in(maxne,8),map(maxnn),mpa(maxnn)
      dimension inb(nf,4),inbl(nf,4),inbr(nf,4)
      
      logical lunsat
c     float(i)=dflotj(i)
c
      nyz=ny*nz
      nez=nz-1
c
c     transport arrays ... 
c     add unsaturated zone to transport grid if required
c     --------------------------------------------------


      if(kcall.eq.1.and.lunsat) then
      nz= nz + nlz(ngz)
      nez= nz-1
      ne= ne + nlz(ngz)*nex*ney
      nn= nn + nlz(ngz)*nx*ny
      nyz=ny*nz
      endif

      if(kgo.eq.4) goto 500            !skip if grid was read from grid3d.plt      

c
c     x-regions
c     ---------
      do 501 i=1,nyz
 501  x(i)=0.
      isumx=0
      do 100 irx=1,ngx
      isumx=isumx+nlx(irx)
      dx=xlim(irx)/float(nlx(irx))
      if(irx.gt.1)dx=(xlim(irx)-xlim(irx-1))/float(nlx(irx))
      do 200 i=(isumx-nlx(irx)+2),(isumx+1)
      do 200 j=1,ny
      do 200 k=1,nz
      node=(i-1)*nyz+(j-1)*nz+k
  200 x(node)=x(node-nyz)+dx
  100 continue
c
c     y-regions
c     ----------
      do 502 i=1,nx
      do 502 k=1,nz
      node=(i-1)*nyz+k
 502  y(node)=0.
      isumy=0
      do 101 iry=1,ngy
      isumy=isumy+nly(iry)
      dy=ylim(iry)/float(nly(iry))
      if(iry.gt.1)dy=(ylim(iry)-ylim(iry-1))/float(nly(iry))
      do 201 j=(isumy-nly(iry)+2),(isumy+1)
      do 201 i=1,nx
      do 201 k=1,nz
      node=(i-1)*nyz+(j-1)*nz+k
  201 y(node)=y(node-nz)+dy
  101 continue
c
c     z-regions
c     top zone for unsaturated transport grid only - 
c     ------------------------------------------------
c     print*,'nn,nz ...',nn,nz
      do 503 k=1,nn-nz+1,nz
 503  z(k)=0.
      isumz=0
      loop=ngz
      if(kcall.eq.0.and.lunsat) loop=ngz-1
      do 102 irz=1,loop
      isumz=isumz+nlz(irz)
      dz=zlim(irz)/float(nlz(irz))
      if(irz.gt.1)dz=(zlim(irz)-zlim(irz-1))/float(nlz(irz))
      do  k=(isumz-nlz(irz)+2),(isumz+1)
      do  j=1,ny
      do  i=1,nx
      node=(i-1)*nyz+(j-1)*nz+k
      z(node)=z(node-1)+dz
      enddo
      enddo
      enddo
  102 continue
c
c     incidences
c     ----------
  500 continue
      nexy=nex*ney
      do 22 k=1,nez
      do 22 j=1,ney
      do 22 i=1,nex
      l=(k-1)*nexy+(j-1)*nex+i
      in1=(i-1)*nyz+(j-1)*nz+k
      in2=in1+nyz
      in3=in2+nz
      in4=in1+nz
      in(l,1)=in1
      in(l,2)=in2
      in(l,3)=in3
      in(l,4)=in4
      in(l,5)=in1+1
      in(l,6)=in2+1
      in(l,7)=in3+1
      in(l,8)=in4+1
c
c
c     flow and transport grid element length arrays ...
c     -------------------------------------------------
      i1=in(l,1)
      i2=in(l,2)
      i3=in(l,3)
      i4=in(l,4)
      i5=in(l,5)
      i6=in(l,6)
      i7=in(l,7)
      i8=in(l,8)
      exl(l)=(x(i2)+x(i3)+x(i7)+x(i6)-x(i1)-x(i4)-x(i8)-x(i5))/4.
      eyl(l)=(y(i4)+y(i3)+y(i7)+y(i8)-y(i1)-y(i2)-y(i6)-y(i5))/4.
      ezl(l)=(z(i5)+z(i6)+z(i7)+z(i8)-z(i1)-z(i2)-z(i3)-z(i4))/4.
c
   22 continue
c
c
c     top surface incidence array {inb}: (uses a LOCAL element index)
c     inb array numbered left->right, then in y
c     -------------------------------------------
      do 23 j=1,ney
      do 23 i=1,nex
      l=(j-1)*nex+i
      in1=(i-1)*nyz + j*nz
      in2=in1+nyz
      in3=in2+nz
      in4=in1+nz
      inb(l,1)=in1
      inb(l,2)=in2
      inb(l,3)=in3
      inb(l,4)=in4
  23  continue
      
c     left surface incidence array {inbl}:
c     inbl array numbered in y then z
c     ---------------------------------------
      do 24 k=1,nez
      do 24 j=1,ney
c     l=(k-1)*nexy+(j-1)*nex + 1   !(global)
      l=(k-1)*ney+j
      in1=(j-1)*nz+k
      in2=in1+nz
      in3=in1+1
      in4=in2+1
      inbl(l,1)=in1
      inbl(l,2)=in2
      inbl(l,3)=in3
      inbl(l,4)=in4
  24  continue     

c     right surface incidence array {inbl}:
c     inbl array numbered in y then z
c     ---------------------------------------
      do 25 k=1,nez
      do 25 j=1,ney
c     l=(k-1)*nexy+(j-1)*nex + nex
      l=(k-1)*ney+j
      in1=(nx-1)*nyz+(j-1)*nz+k
      in2=in1+nz
      in3=in1+1
      in4=in2+1
      inbr(l,1)=in1
      inbr(l,2)=in2
      inbr(l,3)=in3
      inbr(l,4)=in4
  25  continue            
c
c     generate mapping array flow <- transport
c     for a flow domain node i, map(i) gives the transport node #
c     ------------------------------------------------------------
      if(kcall.eq.0) then
      do 343 i=1,nn
      if(lunsat) map(i)=i+((i-1)/nz)*nlz(ngz)
      if(.not.lunsat) map(i)=i
  343 continue
      endif
c
c     second mapping array... for reference vertical strip
c     for a transport node i, mpa(i) gives the equivalent # for 1d strip
c     ------------------------------------------------------------------
      if(kcall.eq.1) then
      do 344 i=1,nn
      mpa(i)= mod(i,nz)
      if(mpa(i).eq.0) mpa(i)=nz
  344 continue
      endif
c
      return
      end
c
c ===========================================================
c =======================================================
c
      subroutine sp3lin (x,y,z,inl,m,xi,yi,zi,det,l,dwx,dwy,dwz,
     +  f,maxnn)
c
c     generation of linear basis functions and derivatives
c
c     ----------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension dwx(8),dwy(8),dwz(8),inl(8),dgx(8),dgy(8),dgz(8),f(8),
     *          sum(9),sinv(9)
      real*8  x(maxnn),y(maxnn),z(maxnn)
c
c     linear integration point value
c     ------------------------------
      xi1=1.-xi
      xi2=1.+xi
      yi1=1.-yi
      yi2=1.+yi
      zi1=1.-zi
      zi2=1.+zi
c
c     corner node shape functions, basic part
c     ---------------------------------------
       f(1)=.125*xi1*yi1*zi1
       f(2)=.125*xi2*yi1*zi1
       f(3)=.125*xi2*yi2*zi1
       f(4)=.125*xi1*yi2*zi1
       f(5)=.125*xi1*yi1*zi2
       f(6)=.125*xi2*yi1*zi2
       f(7)=.125*xi2*yi2*zi2
       f(8)=.125*xi1*yi2*zi2
c
c     shape function derivatives
c     --------------------------
c         dni/dx
      dgx(1)=-.125*yi1*zi1
      dgx(2)=-dgx(1)
      dgx(3)= .125*yi2*zi1
      dgx(4)=-dgx(3)
      dgx(5)=-.125*yi1*zi2
      dgx(6)=-dgx(5)
      dgx(7)= .125*yi2*zi2
      dgx(8)=-dgx(7)
c         dni/dy
      dgy(1)=-.125*xi1*zi1
      dgy(2)=-.125*xi2*zi1
      dgy(3)=-dgy(2)
      dgy(4)=-dgy(1)
      dgy(5)=-.125*xi1*zi2
      dgy(6)=-.125*xi2*zi2
      dgy(7)=-dgy(6)
      dgy(8)=-dgy(5)
c         dni/dz
      dgz(1)=-.125*xi1*yi1
      dgz(2)=-.125*xi2*yi1
      dgz(3)=-.125*xi2*yi2
      dgz(4)=-.125*xi1*yi2
      dgz(5)=-dgz(1)
      dgz(6)=-dgz(2)
      dgz(7)=-dgz(3)
      dgz(8)=-dgz(4)
c     if (l.ne.1) go to 72
c     write (6,71) (dgx(i),i=1,8)
c     write (6,71) (dgy(i),i=1,8)
c     write (6,71) (dgz(i),i=1,8)
c  71 format (9f14.8)
c  72 continue
c
c     jacobian
c     --------
      do 245 lz=1,9
  245 sum(lz)=0.
      do 260 i=1,m
      ki=inl(i)
      sum(1)=sum(1)+dgx(i)*x(ki)
      sum(2)=sum(2)+dgx(i)*y(ki)
      sum(3)=sum(3)+dgx(i)*z(ki)
      sum(4)=sum(4)+dgy(i)*x(ki)
      sum(5)=sum(5)+dgy(i)*y(ki)
      sum(6)=sum(6)+dgy(i)*z(ki)
      sum(7)=sum(7)+dgz(i)*x(ki)
      sum(8)=sum(8)+dgz(i)*y(ki)
      sum(9)=sum(9)+dgz(i)*z(ki)
  260 continue
c     if (l.eq.1) write (6,71) (sum(i),i=1,9)
      det = sum(1)*(sum(5)*sum(9)-sum(6)*sum(8))
     *     -sum(4)*(sum(2)*sum(9)-sum(3)*sum(8))
     *     +sum(7)*(sum(2)*sum(6)-sum(3)*sum(5))
      det1=1./det
c     write (6,71) det,det1
      sinv(1)= det1*(sum(5)*sum(9)-sum(6)*sum(8))
      sinv(2)=-det1*(sum(2)*sum(9)-sum(3)*sum(8))
      sinv(3)= det1*(sum(2)*sum(6)-sum(3)*sum(5))
      sinv(4)=-det1*(sum(4)*sum(9)-sum(6)*sum(7))
      sinv(5)= det1*(sum(1)*sum(9)-sum(3)*sum(7))
      sinv(6)=-det1*(sum(1)*sum(6)-sum(3)*sum(4))
      sinv(7)= det1*(sum(4)*sum(8)-sum(5)*sum(7))
      sinv(8)=-det1*(sum(1)*sum(8)-sum(2)*sum(7))
      sinv(9)= det1*(sum(1)*sum(5)-sum(2)*sum(4))
c     if (l.eq.1) write (6,71) (sinv(i),i=1,9)
c
c     shape function derivatives - global
c     -----------------------------------
      do 270 j=1,m
      dwx(j)=sinv(1)*dgx(j)+sinv(2)*dgy(j)+sinv(3)*dgz(j)
      dwy(j)=sinv(4)*dgx(j)+sinv(5)*dgy(j)+sinv(6)*dgz(j)
      dwz(j)=sinv(7)*dgx(j)+sinv(8)*dgy(j)+sinv(9)*dgz(j)
  270 continue
c     if (l.ne.1) go to 73
c     write (6,71) (dwx(i),i=1,8)
c     write (6,71) (dwy(i),i=1,8)
c     write (6,71) (dwz(i),i=1,8)
c  73 continue
      return
      end
c
c
c
c
c**********************************************************************
c
c subroutine mindex:
c ------------------
c
c routine to assemble boundary condition arrays {ic},{fb},{fc},{lc}
c from inputted boundary conditions (flow and transport).
c
c note for flow:       pass nb --> nb; ktype = 0
c      for heat transport:  pass nw --> nb2; ktype = 1
c      for mass transport c1:  pass nw --> nb3; ktype = 2
c      for mass transport c2:  pass nw --> nb4; ktype = 3
c
c     call mindex(maxnn,maxne,maxn,maxna,laa,np1,nnt,net,nf,n3,
c     + nx,ny,nzt,in2,ic3,fb3,fc3,lc3,fbfc,ib3,a,iaa3,ind3,
c     + nw,kprt,ktype,nc3,m3,nb3,x,y,z,ara,datum,leak3,por,
c     + wh,ampcos,ampsin,wpc,wps,slope,phasewt,map,kb3,lwtfgt0,ffc3)

c ======================================================================
c
      subroutine mindex(maxnn,maxne,maxn,maxna,laa,np1,nn,ne,nf,n,
     + nx,ny,nz,in,ic,fb,fc,lc,fbf,ib,a,iaa,ind,nb,kprt,ktype,
     + nc,m,nb1,x,y,z,ara,datum,leak,por,
     + wh,ampcos,ampsin,wpc,wps,slope,phasewt,map,kb,lwtfgt0,ffc)
c
      implicit real*8(a-h,o-z)
c
      dimension in(maxne,8),ic(maxnn),lc(maxnn),map(maxnn)
      dimension a(maxna),iaa(maxna),ind(np1),ib(maxn,nb)
      dimension fbf(nf)
      dimension kb(7),ara(nf)
      real*8     x(maxnn),y(maxnn),z(maxnn)
      real*8    fb(maxnn),fc(maxnn),por(maxne)
      logical leak,lwtfgt0
c
c
      if(ktype.eq.0) write(6,200)
      if(ktype.eq.1) write(6,201)
      if(ktype.eq.2) write(6,203)
      if(ktype.eq.3) write(6,204)
 200  format(//10x,'mindex routine ... flow boundary arrays ...',/)
 201  format(//10x,'mindex ... heat transport boundary arrays ...',/)
 203  format(//10x,'mindex ... mass transport c1 boundary arrays ...',/)
 204  format(//10x,'mindex ... mass transport c2 boundary arrays ...',/)
c
      nex=nx-1
      ney=ny-1
      nez=nz-1
      nexy=nex*ney
      nyz=ny*nz
c
c     boundary conditions
c     -------------------
c     face 1  x = 0
c          2  x = xl
c          3  y = 0
c          4  y = yl
c          5  z = 0
c          6  z = zl
c          7  internal nodes
c
c     boundary condition type 0  impermeable flow, zero-gradient transport
c                        type 1  dirichlet, fixed head or temperature
c                        type 2  neumann non-zero flux or temp gradient
c                        type 3  cauchy (transport, NOT supported)
c                        type 4  linear leakage boundary on top
c     (type 3 and type 4 presently only for transport, top boundary)
c                        type 5  internal fixed node
c                        type 6  cosine wt flow only
c     ----------------------------------------------------------------
c
      read(5,*) (kb(i),i=1,7)
   11 format(10i5)
      write (6,13) (kb(i),i=1,7)
   13 format (10x,'boundary codes',3x,7i5)
c
c     zero boundary arrays ...
c     ------------------------
      do 913 i=1,nf
      fbf(i)= 0.
 913  continue
c
      do 41 i=1,nn
      ic(i)=0.
      fc(i)=0.
   41 fb(i)=0.
      nxy=nx*ny
c     iwtxlast = 1
c
c     loop over each face and assemble arrays:
c     ----------------------------------------
      do 50 l=1,7
      if((kb(l).eq.3.or.kb(l).eq.4).and.l.ne.6) then
      write(6,122)
 122  format(//10x,'incompatible boundary conditions',
     +        /10x,'only top surface boundary can be type 3 or 4',
     +        /10x,'program stopping')
      stop
      endif

      if(kb(l).eq.6.and.(l.ne.6.or.ktype.ne.0)) then
      write(6,123)
 123  format(//10x,'incompatible boundary conditions',
     +        /10x,'can only have type 6 on top bdy for flow',
     +        /10x,'program stopping')
      stop
      endif
c
c     define number of nodes on this face
c     -----------------------------------
      if(l.eq.1) nnf=ny*nz
      if(l.eq.2) nnf=ny*nz
      if(l.eq.3) nnf=nx*nz
      if(l.eq.4) nnf=nx*nz
      if(l.eq.5) nnf=nx*ny
      if(l.eq.6) nnf=nx*ny
c
c     print*,'ktype:',ktype,' ok 1'

      if (kb(l).eq.0.or.kb(l).eq.3) go to 50
      if (kb(l).eq.2) go to 95
      if (kb(l).eq.4) go to 96
      if (kb(l).eq.5) go to 98
      if (kb(l).eq.6) go to 198
c
c     dirichlet condition:
c     set ffc = -999 to use initial conditions as bounday head/conc
c     or      =-999 to use air temps for Diriclet surface temps.
c     ------------------------------------------------------------
c
 2233 read(5,*) isx1,isx2,isy1,isy2,ffc,more 
 1111 format(4i5,f10.0)
      if(ktype.eq.0) write(6,776) l,isx1,isx2,isy1,isy2,ffc
      if(ktype.eq.1) write(6,777) l,isx1,isx2,isy1,isy2,ffc
      if(ktype.eq.2) write(6,778) l,isx1,isx2,isy1,isy2,ffc
      if(ktype.eq.3) write(6,779) l,isx1,isx2,isy1,isy2,ffc
  776 format(/10x,'Flow: Dirichlet patch source on surface:',i4,
     +    /10x,'indeces: x1,x2, y/z 1, y/z 2 ... ',4i4,'  head: ',f9.4)
  777 format(/10x,'Heat: Dirichlet patch source on surface:',i4,
     +    /10x,'indeces: x1,x2, y/z 1, y/z 2 ... ',4i4,'  temp: ',f9.4)
  778 format(/10x,'Mass: Dirichlet patch source on surface:',i4,
     +    /10x,'indeces: x1,x2, y/z 1, y/z 2 ... ',4i4,'  c1: ',f9.4)
  779 format(/10x,'Mass: Dirichlet patch source on surface:',i4,
     +    /10x,'indeces: x1,x2, y/z 1, y/z 2 ... ',4i4,'  c2: ',f9.4)
      call flush(6)
c
      if( ((l.eq.1.or.l.eq.2).and.isx2.gt.ny) .or.
     +    ((l.eq.1.or.l.eq.2).and.isy2.gt.nz) .or.
     +    ((l.eq.3.or.l.eq.4).and.isx2.gt.nx) .or. 
     +    ((l.eq.3.or.l.eq.4).and.isy2.gt.nz)) then
      write(6,8475) ktype,l
 8475 format(/10x,'Error in Dirichlet boundary indices; ktype = ',i4,
     +       /10x,'face ',i5,/10x,'Program stopping ...')
      stop 
      endif
c      
      do 2231 i=isx1,isx2
      do 2231 j=isy1,isy2
      if(l.eq.1) iglobal_node= (i-1)*nz  + j
      if(l.eq.2) iglobal_node= (nx-1)*nyz + (i-1)*nz  + j
      if(l.eq.3) iglobal_node= (i-1)*nyz  + j
      if(l.eq.4) iglobal_node= (i)*nyz  + j - nz 
      if(l.eq.5) iglobal_node= (i-1)*nyz  + (j-1)*nz + 1
      if(l.eq.6) iglobal_node= (i-1)*nyz  + j*nz
c      print*,'iglobal_node: ,fixed head ',iglobal_node,ffc
      if(iglobal_node.eq.0.or.iglobal_node.gt.nn) then
      write(6,7778) iglobal_node,nn
7778  format(/10x,'Error in Dirichlet b.c.: iglobal_node <0 or >nn)'
     +       /10x,'iglobal_node,nn = ',2i8,
     +       /10x,'Program stopping ...')
      stop
      endif
      fc(iglobal_node)= ffc
c     if(ffc.lt.-998.)  fc(iglobal_node)=z(iglobal_node)    !hardwire -999 = use topo for heads
      ic(iglobal_node)= 1
 2231 continue

c      write(6,9883) (fc(i),i=nz,nn,nz)
c9883  format(' fc array in mindex...',
c     +       (10e12.3))
      call flush(6)
      if(more.gt.0) goto 2233
      go to 96

c     internal fixed nodes
c     ---------------------
 98   continue
      read(5,*) icx1,icx2,icy1,icy2,icz1,icz2,fixed,more

      if(ktype.eq.0) write(6,8352) icx1,icx2,icy1,icy2,icz1,icz2,fixed
 8352 format(/10x,'internal fixed flow nodes: ',
     +       /10x,'icx1,icx2,icy,icy2,icz1,icz2,head: ',6i4,e15.5)
      if(ktype.eq.1) write(6,8353) icx1,icx2,icy1,icy2,icz1,icz2,fixed
 8353 format(/10x,'internal fixed heat transport nodes: ',
     +       /10x,'icx1,icx2,icy,icy2,icz1,icz2,temp: ',6i4,e15.5)
      if(ktype.eq.2) write(6,8354) icx1,icx2,icy1,icy2,icz1,icz2,fixed
 8354 format(/10x,'internal fixed mass transport c1 nodes: ',
     +       /10x,'icx1,icx2,icy,icy2,icz1,icz2,conc: ',6i4,e15.5)
      if(ktype.eq.3) write(6,8355) icx1,icx2,icy1,icy2,icz1,icz2,fixed
 8355 format(/10x,'internal fixed mass transport c2 nodes: ',
     +       /10x,'icx1,icx2,icy,icy2,icz1,icz2,conc: ',6i4,e15.5)
      call flush(6)

      if(icx1.gt.nx.or.icx2.gt.nx.or.icy1.gt.ny.or.icy2.gt.ny
     +             .or.icz1.gt.nz.or.icz2.gt.nz) then
      write(6,7779)
7779  format(/10x,'Error in internal fixed node indices'
     +       /10x,'Program stopping ...')
      stop
      endif

      do i=icx1,icx2
      do j=icy1,icy2
      do k=icz1,icz2
      node = (i-1)*ny*nz + (j-1)*nz + k
      fc(node) = fixed
      ic(node) = 1
      enddo
      enddo
      enddo
      if(more.gt.0) goto 98
      goto 96
c
c     kb(6) = 6  type 1 cosine watertable, for flow 
c     allow multiple wt functions along x
c     if reading hleft in time, can only have 1 wt function
c     -------------------------------------------------------------------
 198  continue
      write(6,9874) 
 9874 format(/10x,'top kb(6)=6, Dirichlet cosine head condition')      
      lwtfgt0 = .false.

c     if(ktype.eq.0) then            !flow only
 199  continue
      read (5,*,err=9900) iwtx1,iwtx2,wh,ampcos,ampsin,
     +                                       wpc,wps,slope,phasewt,more
      if(more.gt.0) lwtfgt0=.true.
c      
c     iwtx1,iwtx2 = x-node indices for wt function
c     wh = watertable rise = cosine curve amplitude
c     amp = sine curve amplitude factor
c     wpc = fraction of wave length for cosine curve
c     wps = fraction of wave length for sine curve
c     slope = linear slope
c
c     print watertable data ...
c     --------------------------
      write(6,317) iwtx1,iwtx2,wh,ampcos,ampsin,wpc,wps,slope,phasewt
  317 format (/10x,'watertable function data:',/
     +         10x,'x-node indices ix1-ix2:     ',2i6/
     +         10x,'watertable rise             ',f16.4/
     +         10x,'cosine curve amplitude      ',f16.4/
     +         10x,'sine curve amplitude        ',f16.4/
     +         10x,'wave fraction, cosine curve ',f16.4/
     +         10x,'wave fraction, sine curve   ',f16.4/
     +         10x,'linear slope                ',f16.4/
     +         10x,'phase shift                 ',f16.4)
      call flush(6)
c
c     watertable curve - wt array
c
      nfirst = iwtx1*ny*nz
      nlast = iwtx2*ny*nz
      pie = 4.d0 * datan(1.d0)
      plc = 2.d0 * wpc * pie/(x(nlast)-x(nfirst))    !x(nn) is last flow node, nn not map(nn)
      pls = 2.d0 * wps * pie/(x(nlast)-x(nfirst))    !x(nn) is last flow node
c
c     generate wt array ...
c     ----------------------------------------------
      do ix=iwtx1,iwtx2
      i3d1 = (ix-1)*ny*nz + nz
      i3d2 = (ix)*ny*nz
      do i3d=i3d1,i3d2,nz
      xx = x(i3d)         !x(i3d) is flow node (transport x not defined yet)
      cosx = dcos(plc*(xx-x(nfirst))+phasewt)
      sinx = dsin(pls*(xx-x(nfirst))+phasewt)
c
c     general curve
c     ---------------
      fc(i3d) =  wh + slope*(xx-x(nfirst)) + ampcos*cosx + ampsin*sinx 
c
c  hardwire Domenico wt
c      BB = cosh(pie*z(nn)/x(nn))
c      AA = x(nn)/2. + BB
c      fc(i3d) = AA-BB*cos(pie*xx/x(nn))

      ic(i3d) = 1
      enddo 
      enddo 

c     iwtxlast = iwtx2
      if(more.gt.0) goto 199
c
c  print wt array ...
c --------------------------------
      write(6,131)
  131 format (/10x,'default watertable heads x,h'/)
      write(6,132) (i,x(i),fc(i),i=nz,nn,nz)
  132 format (2(i7,2e15.5,2x))

c     endif

c     if(ktype.eq.1) then          !transport
c     do i=nzt,nn,nzt
c     fc(i)=-999.
c     ic(i)=1
c     enddo
c     endif


      goto 96

9900  write(6,8479) 
 8479 format(/10x,'error reading cosine watertable parameters',
     +       /10x,'program stopping ...')
      stop
c
c #####################################################################
c
   95 continue
c
c     read recharge or thermal gradient:
c     ----------------------------------
c
 4433 read(5,*) isx1,isx2,isy1,isy2,fff,more 
      if(ktype.eq.0) write(6,886) l,isx1,isx2,isy1,isy2,fff
      if(ktype.eq.1) write(6,887) l,isx1,isx2,isy1,isy2,fff
      if(ktype.eq.2) write(6,888) l,isx1,isx2,isy1,isy2,fff
      if(ktype.eq.3) write(6,889) l,isx1,isx2,isy1,isy2,fff
  886 format(/10x,'Neumann recharge on surface:',i4,
     + /10x,'grid indeces: x1,x2,y1,y2 ... ',4i4,'  Darcy flux: ',e9.2)
  887 format(/10x,'Neumann non-zero gradient on surface:',i4,
     + /10x,'grid indeces: x1,x2,y1,y2 ... ',4i4,' temp gradient:',e9.2)
  888 format(/10x,'Neumann non-zero gradient on surface:',i4,
     + /10x,'grid indeces: x1,x2,y1,y2 ... ',4i4,' c1 gradient:',e9.2)
  889 format(/10x,'Neumann non-zero gradient on surface:',i4,
     + /10x,'grid indeces: x1,x2,y1,y2 ... ',4i4,' c2 gradient:',e9.2)
      do 4431 i=isx1,isx2
      do 4431 j=isy1,isy2
      if(l.eq.1.or.l.eq.2) nj=nz
      if(l.eq.3.or.l.eq.4) nj=nz 
      if(l.eq.5.or.l.eq.6) nj=ny
      index = (i-1)*nj + j
      fbf(index)= fff
 4431 continue
      if(more.gt.0) goto 4433
 4432 continue
c
c #####################################################################
c
   96 continue
c
c     generate index limits for boundary faces ...
c     and generate {ara} array for nodal influence areas
c     --------------------------------------------------
      do 202 i=1,nnf
  202 ara(i)=1.          !hardwire .... needed ?
c
      go to (51,52,53,54,55,56), l
   51 i1=1
      i2=(ny-1)*nz+1
      i3=nz
      j2=nz-1
      j3=1
      call surf(maxnn,nf,nz,1,ara,y,z,l,i1,i2,i3,j2,j3,nnf,kb)
      go to 57
   52 i1=(nx-1)*ny*nz+1
      i2=nn-nz+1
      i3=nz
      j2=nz-1
      j3=1
      call surf(maxnn,nf,nz,1,ara,y,z,l,i1,i2,i3,j2,j3,nnf,kb)
      go to 57
   53 i1=1
      i2=(nx-1)*ny*nz+1
      i3=ny*nz
      j2=nz-1
      j3=1
      call surf(maxnn,nf,ny,1,ara,x,z,l,i1,i2,i3,j2,j3,nnf,kb) 
      go to 57
   54 i1=(ny-1)*nz+1
      i2=nn-nz+1
      i3=ny*nz
      j2=nz-1
      j3=1
      call surf(maxnn,nf,ny,1,ara,x,z,l,i1,i2,i3,j2,j3,nnf,kb)  
      go to 57
   55 i1=1
      i2=(nx-1)*ny*nz+1
      i3=ny*nz
      j2=(ny-1)*nz
      j3=nz
      call surf(maxnn,nf,ny,nz,ara,x,y,l,i1,i2,i3,j2,j3,nnf,kb)
      go to 57
   56 i1=nz
      i2=(nx-1)*ny*nz+nz
      i3=ny*nz
      j2=(ny-1)*nz
      j3=nz
      call surf(maxnn,nf,ny,nz,ara,x,y,l,i1,i2,i3,j2,j3,nnf,kb)
      if(kb(6).eq.4) goto 50
c
c     write (6,777) l,i1,i2,i3,j2,j3
c 777 format (' boundary',6i5)
c
c     assemble {ic},{fb},{fc} for all nodes on this face:
c     ---------------------------------------------------
  57  continue
      if(kb(l).eq.2) then
      k=0
      do 58 i=i1,i2,i3
      jj=j2+i
      do 58 j=i,jj,j3
      k=k+1
      fb(j)=fb(j)+fbf(k)*ara(k)
   58 continue
      endif
c
c     next face:
c     ----------
   50 continue

c     print*,'ktype:',ktype,' ok 1'
c
c     set flag to identify top leaky layer:
c     -------------------------------------
      leak = kb(6).eq.4
c
c      hardwire for smoker cell with impermeable walls
c      if all flow boundaries are impermeable, set the top left corner
c      node to fixed head of 1 m.
c     --------------------------------------------------
c       if(ktype.eq.0) then
c        iflag=0
c        do 778 i=1,nn
c        if(ic(i).eq.1) iflag=1
c 778    continue 
c        if(iflag.eq.0) then
c        fc(nz) = 1.
c        ic(nz) = 1
c       endif
c       endif
c
c      if(kb(6).eq.2.and.ktype.eq.0) write (6,93) (i,fb(i),i=nz,nn,nz)
c   93 format (/10x,'top face neumann boundary,  non-zero',
c     +        ' darcy flux recharge:'/,(4(i7,e12.3)))
c
c ##################################################################
c #####################################################################
c
c     condensation code
c     -----------------
      lc(1)=ic(1)
      do 39 i=2,nn
   39 lc(i)=lc(i-1)+ic(i)
      nc = lc(nn)
      n = nn-nc
      write (6,42) nc,n
   42 format (/10x,'number of dirichlet nodes: ',i5/10x,
     +              'number of degrees of freedom: ',i9)
      call flush(6)
c      if (kprt.eq.0) go to 40
c      write (6,38)
c   38 format (/10x,'boundary arrays'/)
c      do 43 i=1,nn/10
c   43 write (6,44) i,ic(i),lc(i),fc(i),fb(i)
c   44 format (i10,5x,2i5,5x,2f12.2)
c   40 continue
c
c     conjugate gradient matrices
c     ---------------------------
      do 46 i=1,n
      do 46 j=1,nb
   46 ib(i,j)=0
      m=(nb-1)*n
      do 47 i=1,m
      a(i)=0.
   47 iaa(i)=0
      do 48 i=1,n
      ib(i,1)=i
   48 ind(i)=0
c
c     column indices
c     --------------
      nb1=0
      do 60 l=1,ne
      do 61 i=1,8
      ki=in(l,i)
      if (ic(ki).eq.1) go to 61
      ii=ki-lc(ki)
      do 62 j=1,8
      kj=in(l,j)
      if (ic(kj).eq.1) go to 62
      jj=kj-lc(kj)
c     write (6,728) l,i,j,ii,jj,(ib(jj,kk),kk=1,nb)
c 728 format (5i5,5x,10i5)
      if (ii.ge.jj) go to 64
      do 63 k=1,nb
      if (ii.eq.ib(jj,k)) go to 64
      if (ib(jj,k).ne.0) go to 65
      ib(jj,k)=ii
      if (k.gt.nb1) nb1=k
      go to 64
   65 if (k.lt.nb) go to 63
      write (6,166) jj,k,nb
  166 format (10x,'bandwidth insufficient: row',i5,', bw = ',
     +                                         i4,',nb= ',i3)
      nb1=nb+1
   63 continue
   64 continue
c     write (6,729) (ib(jj,kk),kk=1,nb)
c 729 format (30x,10i5)
   62 continue
   61 continue
   60 continue
      if (nb1.le.nb) go to 267
      write (6,268) nb1
  268 format (/10x,'bandwidth required',i5/)
      stop
  267 continue
      write (6,181) nb1
  181 format(10x,'final bandwidth ',i5)
c     do 722 i=1,n
c 722 write (6,721) i,(ib(i,j),j=1,nb)
c 721 format (i5,20i4)
c
c     indices ordering
c     ----------------
      do 66 i=1,n
      do 67 j=2,nb
      if (ib(i,j).eq.0) go to 69
      im=0
      do 68 k=j,nb
      if (ib(i,k).le.im) go to 68
      im=ib(i,k)
      kk=k
   68 continue
      ib(i,kk)=ib(i,j)
      ib(i,j)=im
   67 continue
   69 continue
c     write (6,721) i,(ib(i,j),j=1,nb)
   66 continue
c
c     diagonal position
c     -----------------
      ind(1)=1
      do 71 i=2,n
      do 72 j=1,nb
      if (ib(i-1,j).ne.0) go to 72
      ibm=j-1
      go to 73
   72 continue
   73 ind(i)=ind(i-1)+ibm
   71 continue
c     write (6,731) (i,ind(i),i=1,n)
c 731 format (2i5)
c
c     row indices
c     -----------
      k=0
      do 76 i=1,n
      do 77 j=1,nb
      if (ib(i,j).eq.0) go to 76
      k=k+1
      iaa(k)=ib(i,j)
   77 continue
   76 continue
      m=k
      ind(n+1)=m+1
      write (6,78) m
   78 format (10x,'vector length',i10)
c     write (6,731) (i,iaa(i),i=1,m)
      call flush(6)
c
c     return
      end
c
c ==================================================
c ==================================================
c
c**********************************************************************
c
c subroutine condense:
c ------------------
c
c routine to assemble boundary condition arrays {ic},{fb},{fc},{lc}
c from inputted boundary conditions (flow and transport).
c
c ======================================================================
c
      subroutine condense(maxnn,maxne,maxn,maxna,laa,np1,nn,ne,nf,n,
     + nx,ny,nz,in,ic,fb,fc,lc,fbf,ib,a,iaa,ind,nb,kprt,ktype,
     + nc,m,nb1,x,y,z,ara,datum,leak,por)
c
      implicit real*8(a-h,o-z)
c
      dimension in(maxne,8),ic(maxnn),lc(maxnn)
      dimension a(maxna),iaa(maxna),ind(np1),ib(maxn,nb)
      dimension fbf(nf)
      dimension kb(7),ara(nf)
      real*8     x(maxnn),y(maxnn),z(maxnn)
      real*8    fb(maxnn),fc(maxnn),por(maxne)
      logical leak
c
c     write(6,301)
c301  format(10x,'mindex routine ... transport recondense ...')
c
      nex=nx-1
      ney=ny-1
      nez=nz-1
      nexy=nex*ney
      nyz=ny*nz
c
c ##################################################################
c #####################################################################
c
c     condensation code
c     -----------------
      lc(1)=ic(1)
      do 39 i=2,nn
   39 lc(i)=lc(i-1)+ic(i)
      nc = lc(nn)
      n = nn-nc
      write (6,42) nc,n
   42 format (10x,'recondensing: new # dirichlet nodes: ',i5,
     +                         '; new # degrees of freedom: ',i9)
c
c     conjugate gradient matrices
c     ---------------------------
      do 46 i=1,n
      do 46 j=1,nb
   46 ib(i,j)=0
      m=(nb-1)*n
      do 47 i=1,m
      a(i)=0.
   47 iaa(i)=0
      do 48 i=1,n
      ib(i,1)=i
   48 ind(i)=0
c
c     column indices
c     --------------
      nb1=0
      do 60 l=1,ne
      do 61 i=1,8
      ki=in(l,i)
      if (ic(ki).eq.1) go to 61
      ii=ki-lc(ki)
      do 62 j=1,8
      kj=in(l,j)
      if (ic(kj).eq.1) go to 62
      jj=kj-lc(kj)
c     write (6,728) l,i,j,ii,jj,(ib(jj,kk),kk=1,nb)
c 728 format (5i5,5x,10i5)
      if (ii.ge.jj) go to 64
      do 63 k=1,nb
      if (ii.eq.ib(jj,k)) go to 64
      if (ib(jj,k).ne.0) go to 65
      ib(jj,k)=ii
      if (k.gt.nb1) nb1=k
      go to 64
   65 if (k.lt.nb) go to 63
      write (6,166) jj,k,nb
  166 format (10x,'bandwidth insufficient: row',i5,', bw = ',
     +                                         i4,',nb= ',i3)
      nb1=nb+1
   63 continue
   64 continue
c     write (6,729) (ib(jj,kk),kk=1,nb)
c 729 format (30x,10i5)
   62 continue
   61 continue
   60 continue
      if (nb1.le.nb) go to 267
      write (6,268) nb1
  268 format (/10x,'bandwidth required',i5/)
      stop
  267 continue
c     write (6,181) nb1
c 181 format(10x,'final bandwidth ',i5)
c     do 722 i=1,n
c 722 write (6,721) i,(ib(i,j),j=1,nb)
c 721 format (i5,20i4)
c
c     indices ordering
c     ----------------
      do 66 i=1,n
      do 67 j=2,nb
      if (ib(i,j).eq.0) go to 69
      im=0
      do 68 k=j,nb
      if (ib(i,k).le.im) go to 68
      im=ib(i,k)
      kk=k
   68 continue
      ib(i,kk)=ib(i,j)
      ib(i,j)=im
   67 continue
   69 continue
c     write (6,721) i,(ib(i,j),j=1,nb)
   66 continue
c
c     diagonal position
c     -----------------
      ind(1)=1
      do 71 i=2,n
      do 72 j=1,nb
      if (ib(i-1,j).ne.0) go to 72
      ibm=j-1
      go to 73
   72 continue
   73 ind(i)=ind(i-1)+ibm
   71 continue
c     write (6,731) (i,ind(i),i=1,n)
c 731 format (2i5)
c
c     row indices
c     -----------
      k=0
      do 76 i=1,n
      do 77 j=1,nb
      if (ib(i,j).eq.0) go to 76
      k=k+1
      iaa(k)=ib(i,j)
   77 continue
   76 continue
      m=k
      ind(n+1)=m+1
c     write (6,78) m
c  78 format (10x,'new vector length',i10)
c     write (6,731) (i,iaa(i),i=1,m)
c
c     return
      end
c
c ==================================================
c ==================================================
c
c subroutine surf: generates nodal influence area array for
c neumann flow boundary condition:
c
c true influence area: uses sum of 2 triangles
c c1,c2 are coordinates (either x,y or z arrays)
c
      subroutine surf(maxnn,nf,n1,n2,ara,c1,c2,l,i1,i2,i3,j2,j3,nnf,kb)
      implicit real*8(a-h,o-z)
      real*8 c1(maxnn),c2(maxnn)
      dimension ara(nf),kb(7)
c
      if(kb(l).eq.1) return
      k=0
c      write(6,756)
c  756 format(10x,'subroutine surf')
      do 58 i=i1,i2,i3
      jj=j2+i
      do 58 node=i,jj,j3
      k=k+1
      a1=c1(node)
      a2=a1
      a3=a1
      a4=a1
      b1=c2(node)
      b2=b1
      b3=b1
      b4=b1
      if(i.ne.i1) then
              a1 = c1(node-i3)
              b1 = c2(node-i3)
              endif
      if(i.ne.i2) then
              a2 = c1(node+i3)
              b2 = c2(node+i3)
              endif
      if(node.ne.i)  then
              a3 = c1(node-n2)
              b3 = c2(node-n2)
              endif
      if(node.ne.jj) then
              a4 = c1(node+n2)
              b4 = c2(node+n2)   
              endif
c
c      simplified for rectangular elements:
c      ara(k)=sqrt((a2-a1)**2+(b2-b1)**2)*sqrt((a4-a3)**2+(b4-b3)**2)/4.
c
c     true influence area (based on 1/4 sum of 2 triangles)
c     -----------------------------------------------------
      ara(k)=0.25*( (a2*b4-a4*b2) - (a1*b4-a4*b1)
     +          +  (a3*b2-a2*b3) + (a1*b3-a3*b1) )
c
c
  58  continue
c      write(6,50) (ara(k),k=1,nnf)
c  50  format(/1x,'{ara} ...',/(6e11.4))
      return
      end
c
c ==================================================
c ==================================================
c
c velocity calculation subroutine (called by 'heat')
c veloc1: numerically evaluated derivatives ...
c ice saturation relative k not implemented here yet
c ==================================================
c
      subroutine veloc1(maxnn,maxne,ne,nexy,x,y,z,in,u,t,c2,d2,vx,vy,vz,
     +   por,cx,cy,cz,gamma,map,pp,qq,modelwu,lheat,lmass,lage,ts,omega,
     +   modelkr,lkr1top,
     +   por0,alpha0,alphap,rhob,ps,spn0,spn1,spe1)
c
      implicit real*8(a-h,o-z)
      external wu,dwu
c
      logical lheat,lmass,lage,lkr1top
      real*8   x(maxnn),y(maxnn),z(maxnn)
      real*8  vx(maxne),vy(maxne),vz(maxne)
      real*8 cx(maxne),cy(maxne),cz(maxne),por(maxne)
      real*8 rhob(maxne),por0(maxne),ps(maxne)
      dimension u(maxnn),t(maxnn),map(maxnn),c2(maxnn),d2(maxnn)
      dimension dgx(8),dgy(8),dgz(8),inl(8),f(8)
      dimension in(maxne,8),ag(2),inl2(8)
      dimension pp(maxne),qq(maxne)
      dimension spn0(maxnn),spn1(maxnn),spe1(maxne)
c
      data ag/-.5773502692,.5773502692/
      m=8
c
c     loop over each flow element ...
c     --------------------------------
      do 100 l=1,ne
      do 110 i=1,8
      inl(i) = in(l,i)
 110  inl2(i) = map(inl(i))
c
c     velocities are defined at the gauss points,
c     then averaged for elemental velocities.
c     -------------------------------------------
      sumx=0.d0
      sumy=0.d0
      sumz=0.d0
      tsum=0.d0
      csum=0.d0
c     spnsum=0.d0
c
      do 101 ik=1,2
      zi=ag(ik)
      do 101 ij=1,2
      yi=ag(ij)
      do 101 ii=1,2
      xi=ag(ii)
c
      kl=(ik-1)*4 + (ij-1)*2 + ii
c
      call sp3lin(x,y,z,inl2,m,xi,yi,zi,det,l,dgx,dgy,dgz,f,maxnn)
c
c     summation loop:
c     ---------------
      do 200 i=1,8
      sumx=sumx+u(inl(i))*dgx(i)
      sumy=sumy+u(inl(i))*dgy(i)
      sumz=sumz+u(inl(i))*dgz(i)
 200  continue
c
      tsum=tsum+t(inl2(kl))
      csum=csum+c2(inl2(kl))
c     spnsum=spnsum+spn1(inl2(kl))
c
 101  continue
c
      tavg= tsum * 0.125d0
      cavg= csum * 0.125d0
c     spnavg= spnsum * 0.125d0
c
c     ice saturation relative k
c     -------------------------
      swe = wu(tavg,pp(l),qq(l),modelwu,ts) 
c      zkrw = ((swe-p)/(1.-p))**4
c      zkrw = max(zkrw,1.d-6)      
      alpha = alpha0 + por(l)*alphap*(rhob(l)/ps(l))*spe1(l)
      zkrw = fnzkrw(swe,pp(l),omega,por(l),modelkr,por0(l),alpha0,alpha)             !relative permeability function
      if(lkr1top .and. l.ge.(ne-nexy+1)) zkrw=1.d0                    !hardwire keep kr=1 for top surface elements
c
c     elemental velocities (avg. of gauss point values):
c     --------------------------------------------------
      porvs = por(l)*swe*rvisc(tavg,lheat,lmass,lage)
      vx(l) = -0.125d0 *(cx(l)*zkrw/porvs) * sumx
      vy(l) = -0.125d0 *(cy(l)*zkrw/porvs) * sumy
      vz(l) = -0.125d0 *(cz(l)*zkrw/porvs) *
     +              (sumz + 8.d0*rden(tavg,cavg,gamma,lheat,lmass,lage))
c
c     close element loop
c     ------------------
 100  continue
c
      return
      end
c
c ################################################### 
c
c velocity calculation subroutine (called by 'heat')
c veloc2: directly evaluated derivatives ... kint=0
c ===================================================
c
      subroutine veloc2(maxnn,maxne,ne,nexy,x,y,z,in,u,t,c2,d2,vx,vy,vz,
     + por,cx,cy,cz,dwx,dwy,dwz,exl,eyl,ezl,gamma,map,pp,qq,modelwu,
     + lheat,lmass,lage,ts,omega,modelkr,lkr1top,
     + por0,alpha0,alphap,rhob,ps,spn0,spn1,spe1)
c
      implicit real*8(a-h,o-z)
      external wu,dwu
c
      logical lheat,lmass,lage,lkr1top
      real*8   x(maxnn),  y(maxnn),  z(maxnn)
      real*8  vx(maxne), vy(maxne), vz(maxne)
      real*8  cx(maxne), cy(maxne), cz(maxne),por(maxne)
      real*8 rhob(maxne),por0(maxne),ps(maxne)
      real*8 exl(maxne),eyl(maxne),ezl(maxne)
      dimension u(maxnn),t(maxnn),c2(maxnn),d2(maxnn)
      dimension inl(8),inl2(8),in(maxne,8)
      dimension dwx(8),dwy(8),dwz(8)
      integer map(maxnn)
      dimension pp(maxne),qq(maxne)
      dimension spn0(maxnn),spn1(maxnn),spe1(maxne)
c
c     loop over each element ...
c     --------------------------
      do 100 l=1,ne
c
      do 110 i=1,8
      inl(i) = in(l,i)
 110  inl2(i) = map(inl(i))
c
c     velocities are defined at the centroids
c     -------------------------------------------
      sumx=0.d0
      sumy=0.d0
      sumz=0.d0
      tsum=0.d0
      csum=0.d0
c     spnsum=0.d0
c
c     summation loop:
c     ---------------
      do 200 i=1,8
      sumx=sumx+u(inl(i))*dwx(i)
      sumy=sumy+u(inl(i))*dwy(i)
      sumz=sumz+u(inl(i))*dwz(i)
      tsum=tsum+t(inl2(i))
      csum=csum+c2(inl2(i))
c     spnsum=spnsum+spn1(inl2(i))
 200  continue
c
      tavg= tsum * 0.125d0
      cavg= csum * 0.125d0
c     spnavg= spnsum * 0.125d0
c
c     ice saturation relative k
c     -------------------------
      swe = wu(tavg,pp(l),qq(l),modelwu,ts) 
c     zkrw = ((swe-p)/(1.-p))**4
c     zkrw = max(zkrw,1.d-6)
      alpha = alpha0 + por(l)*alphap*(rhob(l)/ps(l))*spe1(l)
      zkrw = fnzkrw(swe,pp(l),omega,por(l),modelkr,por0(l),alpha0,alpha)             !relative permeability function
      if(lkr1top .and.l.ge.(ne-nexy+1)) zkrw=1.d0                    !hardwire keep kr=1 for top surface elements
c
c     elemental velocities 
c     ---------------------
      porvs = por(l)*swe*rvisc(tavg,lheat,lmass,lage)
      vx(l) = - (cx(l)*zkrw/porvs) * (sumx/(4.d0*exl(l)))
      vy(l) = - (cy(l)*zkrw/porvs) * (sumy/(4.d0*eyl(l)))
      vz(l) = - (cz(l)*zkrw/porvs) * (sumz/(4.d0*ezl(l))
     +                      + rden(tavg,cavg,gamma,lheat,lmass,lage))
      
c     hardwire v - minimum v: (v (only in fracs ?) becomes 0 or very small if Ss > 1e-5) pcg error
      if(dabs(vx(l)).lt.1.d-40) vx(l) = 1.d-40
      if(dabs(vy(l)).lt.1.d-40) vy(l) = 1.d-40
      if(dabs(vz(l)).lt.1.d-40) vz(l) = 1.d-40
c
c     close element loop
c     ------------------
 100  continue
c
      return
      end
c
c ##################################################################
c
c     this routine determines the elemental velocities for the 
c     1D line or 2D plane fracture elements - assume line porosity = 1
c 
c     flux equation from Bear 1972
c     =============================
      subroutine vfracture(maxnn,maxne,ne,in,ezl,exl,eyl,u,t,c2,d2,vlin,
     +      ckl,map,xarea,ifracl,lvert,nfrac,nex,ney,nez,maxfrac,ifdim,
     +      vx2d,vy2d,gamma,lheat,lmass,lage,
     +      p2,q2,ts2,omega2,modelwu2,modelkr2,
     +      por0,alpha0,alphap,rhob,ps,spn0,spn1)
c
      implicit real*8(a-h,o-z)
      logical lheat,lmass,lage
      real*8 vlin(maxfrac),vx2d(maxfrac),vy2d(maxfrac),xarea(maxfrac)
      real*8 ezl(maxne),exl(maxne),eyl(maxne)
      real*8 rhob(maxne),por0(maxne),ps(maxne)
      dimension map(maxnn),u(maxnn),t(maxnn),c2(maxnn),d2(maxnn)
      dimension spn0(maxnn),spn1(maxnn)
      integer in(maxne,8),inline(4)
c
c     fracture element velocities
c     ----------------------------
      dimension ckl(maxfrac),ifracl(maxfrac)
      dimension lvert(maxfrac),ifdim(maxfrac)
c
      do 2011 kf=1,nfrac
c
c     1D fractures:
c     -------------
      if(ifdim(kf).eq.1) then
      l =  ifracl(kf) 
      call frac_line(maxfrac,maxne,lvert,inline,in,kf,l,dim,
     +                                       exl,eyl,ezl)
      tb = 0.50d0 * (t(map(inline(1))) + t(map(inline(2))))      !bug fix ok:   map(inline()) everywhere
      cb = 0.50d0 * (c2(map(inline(1))) + c2(map(inline(2))))      !bug fix ok:   map(inline()) everywhere
      spnavg = 0.50d0 * (spn1(map(inline(1))) + spn1(map(inline(2))))      !
c
c     ice saturation relative k
c     -------------------------
      porfrac = 1.0
      swe = wu(tb,p2,q2,modelwu2,ts2) 
      alpha = alpha0 + porfrac*alphap*(rhob(l)/ps(l))*spnavg
      zkrw = fnzkrw(swe,p2,omega2,porfrac,modelkr2,porfrac,alpha0,alpha)             !relative permeability function

      rdd = 0.d0
      if(lvert(kf).eq.1.or.lvert(kf).eq.4) 
     +                          rdd = rden(tb,cb,gamma,lheat,lmass,lage)
c     viscosity = 160.d0 * rvisc(tb,lmass,lage)/ 86400.d0
      vlin(kf)= - ((ckl(kf)*zkrw)/rvisc(tb,lheat,lmass,lage)) * 
     +                     ((u(inline(2))-u(inline(1)))/dim+rdd)   !flow nodes inline so don't need map()
c -----------------------------------------------------------------
c       if(kf.lt. 5) 
c     + write(6,992) l,u(inline(2)),u(inline(1)),tb,rden(tb),rvisc(tb,lmass,lage),
c     +                 dim,vlin(kf)
c  992 format(10x,' element # ',i5,' u2,u1 ... ',2e12.4,
c     +       /10x,' avg temp',e10.3,' rden: ',e10.3,' rvisc: ',e10.3,
c     +            'dim: ',e10.3,/10x,'vlin(kf) ... ',e12.3,//)
c      write(6,222) nfrac,(vlin(i),i=1,nfrac/100)
c 222  format(/10x,'1D line velocities: nfrac= ',i5,/,(10x,5e12.3))
c -----------------------------------------------------------------
      endif
c
c      2D plane fractures - velocities
c      --------------------------------
      if(ifdim(kf).eq.2) then
      l =  ifracl(kf) 
c
      call frac_plane(maxfrac,maxne,lvert,inline,in,kf,l,fdimx,fdimy,
     +             exl,eyl,ezl,rexel,sarea)
c
      tb = 0.25d0*(t(map(inline(1)))+t(map(inline(2)))                  !bugfix map()
     +            +t(map(inline(3)))+t(map(inline(4))))
      cb = 0.25d0*(c2(map(inline(1)))+c2(map(inline(2)))                  !bugfix map()
     +            +c2(map(inline(3)))+c2(map(inline(4))))
      spnavg = 0.25d0*(spn1(map(inline(1)))+spn1(map(inline(2)))                  !bugfix map()
     +            +spn1(map(inline(3)))+spn1(map(inline(4))))
c
c     viscosity = 160.d0 * rvisc(tb,lmass,lage)/ 86400.d0
c
c     ice saturation relative k
c     Feb 2022
c     -------------------------
      porfrac = 1.0
      swe = wu(tb,p2,q2,modelwu2,ts2) 
      alpha = alpha0 + porfrac*alphap*(rhob(l)/ps(l))*spnavg
      zkrw = fnzkrw(swe,p2,omega2,porfrac,modelkr2,porfrac,alpha0,alpha)             !relative permeability function

      rdd = 0.d0
      if(lvert(kf).ne.5 .and. lvert(kf).ne.6) 
     +                            rdd=rden(tb,cb,gamma,lheat,lmass,lage)
      dhdlx = (u(inline(2))-u(inline(1))+u(inline(3))-u(inline(4)))/
     +        (2.d0*fdimx)
      dhdly = (u(inline(4))-u(inline(1))+u(inline(3))-u(inline(2)))/
     +        (2.d0*fdimy)
c
      vx2d(kf)= -((ckl(kf)*zkrw)/rvisc(tb,lheat,lmass,lage)) * dhdlx 
      vy2d(kf)= -((ckl(kf)*zkrw)/rvisc(tb,lheat,lmass,lage))*(dhdly+rdd)

c     hardwire v - minimum fracture v: (v (only in fracs ?) becomes 0 or very small if Ss > 1e-5) pcg error
      if(dabs(vx2d(kf)).lt.1.d-40) vx2d(kf) = 1.d-40
      if(dabs(vy2d(kf)).lt.1.d-40) vy2d(kf) = 1.d-40
c
c     if(kf.le.2) write(6,8888) l,tb,viscosity,fdimx,fdimy,dhdlx,dhdly,
c    +              rdd,vx2d(kf),vy2d(kf)
c8888 format(/10x,'in vfracture ... l= ',i5,'  tbar= ',e12.3,
c    +       /10x,'viscosity= ',e12.3,' fdimx= ',e12.3,' fdimy= ',e12.3,
c    +       /10x,'dhdlx= ',e12.4,' dhdly= ',e12.4,' rdd= ',e12.4,
c    +       /10x,' vx,vy: ',2e12.3)
c     if(kf.le.2) write(6,8887) u(inline(1)),u(inline(2)),u(inline(3)),u(inline(4))
c8887 format(/10x,'u(1,2,3,4) ... ',/10x,4e12.4)
c
      endif
c
c      next fracture element:
c      ----------------------
 2011  continue
c
      return
      end
c
c ==================================================================
c
c routine to check convergence of the watertable, and to
c   deform the finite element grid vertically to account
c   for watertable mounding.
c __________________________________________________________________
c
      subroutine deform(maxnn,nn,nx,ny,nz,z,u,nwtl,lwtc,ccw,
     + datum,map,ezl,nlz,ngz,mxgz,zlim,nzt,nnt,in2,
     + net,nex,ney,maxne,lunsat)
c
c     nwtl is the number of vertical layers which may be deformed
c     j is node number of top surface node
c     -----------------------------------------------------------
c
      implicit real*8(a-h,o-z)
c
      real*8 z(maxnn),ezl(maxne)
      integer map(maxnn),nlz(mxgz),in2(maxne,8)
      dimension u(maxnn),zlim(mxgz)
      logical lwtc,lunsat
c
c     float(i)=dflotj(i)
c     ---------------------------
c
c     loop over top surface nodes to find maximum difference:
c     -------------------------------------------------------
c
      dwmx=-999.
c
      do 101 j=nz,nn,nz
      adw1= u(j)-z(map(j))+datum
      adw2= abs(adw1)
      if(adw2.gt.dwmx) then
         dwmx=adw2
         dwmx2=adw1
         inode=j
      endif
  101 continue
c
      lwtc= (dwmx.lt.ccw)
      write(6,200) dwmx2,inode
c      write(*,200) dwmx2,inode
 200  format(10x,'max. difference at watertable: ',e12.4,' node: ',i7)
      call flush(6)

      if(lwtc) return
c
c     watertable has not converged,
c     loop over vertical nodes to deform grid ...
c     also find minimum sat'd thickness
c     -------------------------------------------
c
      smin = 999.
      do 201 j=nz,nn,nz
c
      sat= u(j) + datum - z(map(j-nz+1))
      if(sat.lt.smin) smin=sat
      dz = (u(j)+datum-z(map(j-nwtl))) / (float(nwtl))
      if(dz.le.0.) then
        write(6,1010) j
 1010   format(/10x,'negative dz detected ... surface node# ',i8,
     +  /10x,'your flow grid is trying to deform down too much',
     +  /10x,'check flow boundary conditions, number of deformable',
     +  /10x,'layers (nwtl), datum, grid coordinates and dimensions',
     +  //10x,'program stopping ...')
        stop
      endif
      do 202 k=1,nwtl
      node= j-nwtl+k
      z(map(node))= z(map(node-1)) + dz
  202 continue
c
  201 continue
c     write(6,350) smin
c350  format(10x,'min. saturated thickness: ',f10.3)
c
c     adjust remaining z-coords in conductive layer...
c     'i' is the transport grid node
c     ------------------------------------------------
      if(lunsat) then
      do 444 i=map(nz),map(nn),nzt
      dz= (zlim(ngz)-z(i))/nlz(ngz)
      do 443 k=(i+1),(i+nlz(ngz))
      z(k)=z(k-1)+dz
  443 continue
  444 continue
      endif
c
c     adjust z-element lengths from nwtl up
c     -------------------------------------
      lstrt=net-nex*ney*nwtl + 1
      if(lunsat) lstrt=net-nex*ney*(nwtl+nlz(ngz)) + 1
      do 445 l=lstrt,net
      i1=in2(l,1)
      i2=in2(l,2)
      i3=in2(l,3)
      i4=in2(l,4)
      i5=in2(l,5)
      i6=in2(l,6)
      i7=in2(l,7)
      i8=in2(l,8)
      ezl(l)=(z(i5)+z(i6)+z(i7)+z(i8)-z(i1)-z(i2)-z(i3)-z(i4))/4.
  445 continue
c      
      return
      end
c
c ======================================================================
c
c     optional frost heave or thaw settlement due to phase change only
c     applied to selected elements at top only ... adjust later
c     -----------------------------------------------------------

      subroutine settlement(maxnn,maxne,nf,nnt,nxy,nzt,nexy,nezt,net,
     + mxgz,nex,ney,t0,t2,z,dz,dztot,ezl,sumdz,por,sw,in2,p,q,modelwu,
     + ts,link,nlz,ngz,dzmin,dzmax,dzl,pi,lheat,lmass,lage,lunsat)
c
      implicit real*8(a-h,o-z)
      dimension t0(maxnn),t2(maxnn),z(maxnn)
      dimension por(maxne),sw(maxne),ezl(maxne),dzl(maxne)
      dimension dz(nf),sumdz(nf),dztot(nf)
      integer link(nf),in2(maxne,8),nlz(mxgz)
      logical lmass,lage,lunsat,lheat
c     
c     initialize dz=0
      dzmin = 0.
      dzmax = 0.
      do i=1,nxy
      dz(i)=0.
      enddo
c      
c     loop over xy elements
      do i=1,nexy
      sumdz(i)=0.
c     loop over vertical elements
      do j=1,nezt
      l = (j-1)*nexy + i
      tavg0=0.
      tavg2=0.
      do k=1,8
      tavg0 = tavg0 + t0(in2(l,k))      !start of time step
      tavg2 = tavg2 + t2(in2(l,k))      !latest iteration
      enddo
      tavg0=tavg0/8.
      tavg2=tavg2/8.
      tavg = (tavg0+tavg2)/2.
      pf = den(tavg,lheat,lmass,lage)
      deltawu = wu(tavg2,p,q,modelwu,ts)-wu(tavg0,p,q,modelwu,ts)     ! + means thaw = settlement
      dzl(l) = por(l)*sw(l)*deltawu*ezl(l)*(pf-pi)/pi        ! dz this element: + means thaw = settlement
      sumdz(i)=sumdz(i) + dzl(l)                             ! cumulative from bottom to top, i = local nexy #
c     write(6,2828) i,j,l,pf,pi,ezl(l),por(l),sw(l),tavg0,tavg2,
c    +              deltawu,dzl(j),sumdz(i)
c2828 format(10x,'settlement: i,j,l,pf,pi,ezl,por,sw,t0,t2,deltawu,',
c    +           'dzl,sumdz(i): ',3i5,10e12.4)
      enddo   !j=1,nezt
c     apply to 4 xy nodes in this element column to allow contour plotting in 2D surface

      do k=1,4
      nodexy = (in2(net-nexy+i,k+4))/nzt            ! nodexy = 1,2,3...nxy
      dz(nodexy) = dz(nodexy) + sumdz(i)            ! +dz is settlement
c     write(6,4098) nodexy,dz(nodexy),sumdz(i)
c4098 format(10x,'settlement ... nodexy,dz,sumdz ... ',i5,2e15.5)      
      enddo

      enddo       !i=1,nexy
c      
c     get average dz for each nodal column based on # element connections link(i)
      do i=1,nxy
      dz(i) = dz(i)/float(link(i))
      if(dz(i).gt.dzmax) dzmax=dz(i)     ! max thaw settlement
      if(dz(i).lt.dzmin) dzmin=dz(i)     ! min thaw settlement
      enddo
c     apply total dz to minimum of top 120 nodes or nzt (adjust later)
      ndz = min(120,nzt)
      do i=nzt,nnt,nzt
      icolxy = i/nzt
      do k=i-ndz+1,i
      z(k) = z(k) - dz(icolxy)/ndz        !thaw is -dz so decrease z
      enddo
      enddo
c
c     update z-element lengths
c     -------------------------------------
      lstrt=net-nex*ney*(ndz-1) + 1
      if(lunsat) lstrt=net-nex*ney*(ndz-1+nlz(ngz)) + 1
      do l=lstrt,net
      i1=in2(l,1)
      i2=in2(l,2)
      i3=in2(l,3)
      i4=in2(l,4)
      i5=in2(l,5)
      i6=in2(l,6)
      i7=in2(l,7)
      i8=in2(l,8)
      ezl(l)=(z(i5)+z(i6)+z(i7)+z(i8)-z(i1)-z(i2)-z(i3)-z(i4))/4.
      enddo

      return
      end
c
c ********************************************************************
c ***********************************************************
c ----------------------------------------------------------------------
c
c subroutine volume - determines elemental volumes of
c hexahedral elements. (uses routine sp3lin)
c
      subroutine volume(x,y,z,maxnn,ag,hag,in,vt,maxne,ne,vol,kint,
     +               exl,eyl,ezl)
c
      implicit real*8 (a-h,o-z)
      dimension detj(27),f(8),dgx(8),dgy(8),dgz(8)
      real*8    x(maxnn),y(maxnn),z(maxnn)
      real*8   vt(maxne),exl(maxne),eyl(maxne),ezl(maxne)
      dimension ag(9),hag(9),inl(8),ff(8,27),in(maxne,8)
c
      vol=0.
      do 22 l=1,ne
c
      do 23 i=1,8
 23   inl(i)=in(l,i)
c
c     option if direct integration used:
c     ##################################
      if(kint.eq.0) then
      vt(l)=exl(l)*eyl(l)*ezl(l)
      vol=vol+vt(l)
      else
c
c     calculate exact element volume determined numerically
c     #####################################################
c     number of gauss points m=8
c     --------------------------
      m=8
      np=2
      npx=np
      npy=np
      npz=np
      np2=npx*npy*npz
      ns=0
      if (np.eq.3) ns=2
c
c     placement of gauss points
c     -------------------------
      do 310  ik=1,npz
      k=ik+ns
      zi=ag(k)
      hz=hag(k)
      do 310 ij=1,npy
      j=ij+ns
      yi=ag(j)
      hy=hag(j)
      do 310 ii=1,npx
      i=ii+ns
      xi=ag(i)
      hx=hag(i)
      kl=(ik-1)*npx*npy+(ij-1)*npx+ii
c
c     basis functions and derivatives
c     -------------------------------
      call sp3lin (x,y,z,inl,m,xi,yi,zi,det,l,dgx,dgy,dgz,f,maxnn)
      do 300 jj=1,m
  300 ff(jj,kl)=f(jj)
      detj(kl)=det*hz*hy*hx
  310 continue
c
c     elemental volumes
c     -----------------
      vt(l)=0.
      do 51 i=1,m
      do 51 k=1,np2
      vt(l) = vt(l) + ff(i,k) * detj(k)
   51 continue
c
      vol=vol+vt(l)
c
      endif
c
c      if(l.le.5) write(6,77) l,vt(l)
c 77   format(1x,'vt(l) for element ',i6,5x,e10.4)
c      if(l.eq.3) then
c      print*,'stop after 3 elements ...'
c      stop
c      endif
c
   22 continue
c
       write(6,78) vol
  78   format(10x,'total domain volume ...              ',e12.3)
c
      return
      end
c ***********************************************************
c ==============================================================
c **************************************************************
c
c     this routine determines first- and second- moments
c     of the heat anomaly (above background)
c     --------------------------------------
c
      subroutine moment(maxnn,maxne,in2,nx,ny,nzt,nyzt,x,z,t2,knoy,net,
     + mpa,xbar,zbar,tpeak,tmin,xvar,zvar,xzvar,xbpk,zbpk,tsa)
c
      implicit real*8(a-h,o-z)
      real*8 x(maxnn),z(maxnn)
      dimension t2(maxnn)
      integer in2(maxne,8),knoy(5),mpa(maxnn)
c
      nex=nx-1
      nez=nzt-1
      ney=ny-1
c
c     find centre of mass in xz section at knoy(1):
c     (first moment) - output to file #17
c     ==============================================
c     peak & minimum temp:
c     --------------------
      tpeak= -999.
      tmin = +999.
      if(knoy(1).eq.0) return
      do 184 i=1,nx
      do 184 k=1,nzt
      node=(i-1)*nyzt+(knoy(1)-1)*nzt+k
      if(t2(node).gt.tpeak) tpeak=t2(node)
      if(t2(node).lt.tmin)  tmin =t2(node)
  184 continue
c
c     moment calculation: 
c     ====================
      sumxt=0.
      sumzt=0.
      sumxx=0.
      sumzz=0.
      sumxz=0.
      sumt=0.
      do 185 l=1,net
      xb=0.
      zb=0.
      tb=0.
      tbt=0.
c
c     get elemental x,z centres, and elemental temperature
c     ----------------------------------------------------
      do 187 i=1,8
      xb=xb+x(in2(l,i))
      zb=zb+z(in2(l,i))
      tb=tb+ (t2(in2(l,i))-t2(mpa(in2(l,i))))
      tbt=tbt+t2(in2(l,i))
 187  continue
c
      xb=xb/8.
      zb=zb/8.
      tb=tb/8.
      sumxt=sumxt+xb*tb
      sumzt=sumzt+zb*tb
      sumxx=sumxx+xb*xb*tb
      sumzz=sumzz+zb*zb*tb
      sumxz=sumxz+xb*zb*tb
      sumt=sumt+tb
  185 continue
c
      xbar= sumxt/sumt
      zbar= sumzt/sumt
      xvar= (sumxx/sumt-xbar*xbar)
      zvar= (sumzz/sumt-zbar*zbar)
      xzvar=(sumxz/sumt-xbar*zbar)
c
c     position of plume peak: centre 2 xz element planes only, 
c     within tsa*tpeak
c     --------------------------------------------------------------
      sumxb=0.
      sumzb=0.
      sumtb=0.
      tw=tsa*tpeak
      do 201 kl=1,nez
      do 202 il=1,nex
      do 203 jl=(knoy(1)-1),knoy(1)
c      kny=knoy(1)
c      if(kny.eq.ny) kny=ny-1
      l=(kl-1)*nex*ney + (jl-1)*nex + il
      if(l.lt.1) goto 203
      xb=0.
      zb=0.
      tb=0.
      ik=0
      do 204 jj=1,8
      if(t2(in2(l,jj)).ge.tw) then
      xb=xb+x(in2(l,jj))
      zb=zb+z(in2(l,jj))
      tb=tb+t2(in2(l,jj))
      ik=ik+1
      endif
  204 continue
      if(ik.ne.0) then
      sumxb=sumxb+xb*tb/((ik)**2)
      sumzb=sumzb+zb*tb/((ik)**2)
      sumtb=sumtb+tb/ik
      endif
 203  continue
 202  continue
 201  continue
c
      xbpk=sumxb/sumtb
      zbpk=sumzb/sumtb
c
      return
      end
c
c ***************************************
c ============================================================
c
c     external function defining unfrozen moisture content
c     as a function of temperature.
c     currently, shift ts is only used for model 3
c     ------------------------------------------------
      function wu(tavg,p,q,modelwu,ts)
      implicit real*8(a-h,o-z)

      wu=1.0              !default for no effect - set modelwu=0
c     return          !hardwire
      
c function 1 (Molson et al WRR 1992)
      if(modelwu.eq.1) then
      if(tavg.lt.0.) wu = p + (1.0-p)*exp(q*tavg)
      return
      endif

c function 2 (SHEMAT)
      if(modelwu.eq.2) then
      if(tavg.lt.0.) wu = exp(-(tavg/q)**2)
      return
      endif

c function 3 (linear, p=residual Wu, slope = q ;  
c for neumann solution with only 2 zones, use (p-1)/q = Ts)
      if(modelwu.eq.3) then
      if(tavg.lt.ts) then
      Bt = ts + ((p-1.)/q)         !temperature for residual Wu (p)
      wu = p + q * (tavg-Bt)
      endif
      if(tavg.lt.Bt) wu=p
      return
      endif
      
c function 4 (Interfrost 2014, p=residual Wu, q=W ;  
c Interfrost 2014 benchmark has p=0.05 and q=W=0.5
c -------------------------------------------------
      if(modelwu.eq.4) then
      if(tavg.lt.0.0) wu = (1.0-p) * exp(-((tavg/q)**2)) + p
      return
      endif      
c
      return
      end
c ***************************************
c ============================================================
c
c     external function defining change in unfrozen moisture content
c     as a function of temperature.
c     ts only used in model 3
c     ------------------------------
      function dwu(tavg,p,q,modelwu,ts)
      implicit real*8(a-h,o-z)

      dwu=0.d0
c     return          !hardwire      
c
c function 1 (WRR)
c     exact derivative: 
      if(modelwu.eq.1) then
      if(tavg.lt.0.0) dwu = (1.d0-p)*q*exp(q*tavg)
      return
      endif

c function 2 (SHEMAT)
c bug corrected Dec 2014 in dwu
      if(modelwu.eq.2) then
      if(tavg.lt.0.0) dwu= -(2.*(tavg/q)*exp(-(tavg/q)**2))
      return
      endif

c function 3 (linear)
      if(modelwu.eq.3) then
      if(tavg.lt.ts .and. tavg.gt.(ts+(p-1.)/q)) dwu= q
      return
      endif

c function 4 (Interfrost 2014)
      if(modelwu.eq.4) then
      if(tavg.lt.0.0) dwu= (1.-p) * (-2.*(tavg/q**2)*exp(-(tavg/q)**2))
      return
      endif

c      write(6,4848) z1,z2,dwu,dwu2
c 4848 format(' z1,z2,dwu,dwu2: ',4e12.3)

      return
      end
c
c     function for relative permeability
c     Interfrost: kr_min = 1e-6
c     -------------------------------------
      function fnzkrw(swe,p,omega,por,modelkr,por0,alpha0,alpha)
      implicit real*8(a-h,o-z)

      fnzkrw=1.d0                   !default for no effect - set modelkr=0
c     return          !hardwire

c     Original (from dnapl)
      if(modelkr.eq.1) then
      fnzkrw = ((swe - p)/(1.d0-p))**4  
      fnzkrw = fnzkrw * ((por/por0)**3 * (alpha0/alpha)**2)
      fnzkrw = max(fnzkrw,1.e-6)
      return
      endif

c     Interfrost 2014:
      if(modelkr.eq.2) then
      fnzkrw = 10**(-por*omega*(1.-swe))
      fnzkrw = fnzkrw * ((por/por0)**3 * (alpha0/alpha)**2)
      fnzkrw = max(fnzkrw,1.e-6)
      return
      endif

      return
      end

c ==============================================================
c *******************************************************************
c
c     function defining variation of surface air temp with time
c     also used for bz and sat
c
c     gradt added March 2016 ... be careful with restarts ... uses total days
c  -----------------------------------------------------------------------------
c
c     parameters:
c                    surfmin   ..... minimum surface air temperature
c                    amp    ........ sine curve half-amplitude
c                    period   ...... sine curve period (days)
c                    phase  ........ phase shift (days)
c                    day  .......... time (days)
c
      function surfat(surfmin,amp,period,phase,day,cutoff,gradt)
      implicit real*8(a-h,o-z)
c
c     hardwire sin or cosine
c      surfat = (amp+surfmin) + amp*sin(2.*3.1415926/period* (day+phase))
      surfat = ( (amp+surfmin) 
     +       + amp*cos(2.*3.1415926/period* (day+phase))) + (gradt*day)
      surfat = max(surfat,cutoff)
c 
      return
      end
c
c *******************************************************************
c
c *******************************************************************
c function defining relative viscosity and relative density
c    of the fluid as a function of temperature (c).
c    viscosity is defined relative to 0 celsius, and density
c    is relative to 4 celsius (temp. of maximum density).
c    for no density effect, set rvisc = 1. and rden = 0.
c    see Fontaine et al. Earth & Planetary Sc. Letters 184 (2001) 407-425 for best 0-400 range
c
c    tmp needs to be temperature
c ---------------------------------------------------------------------------------------------------
c
c     relative viscosity change:
c     ===========================
c
      function rvisc(tmp,lheat,lmass,lage) 
      implicit real*8(a-h,o-z)
      logical lheat,lmass,lage
c   
      rvisc=1.d0
c     return           !hardwire test

      if(lheat) 
c    + rvisc=1.d0           !hardwire for Lunardini
     + rvisc =  dexp((-.03073d0 + 1.303d-04 * tmp) * tmp)     !best for 0-100 degrees
c     rvisc =  dexp((-.0254d0  + 0.443d-04 * tmp) * tmp)     !original for 0-400 degrees: bad for >200
c    +rvisc =  (2.414e-5*(10.**(247.8/(tmp+133.)))) / 1.d-3  !best for 0-400 Fontaine et al.

c     if(lage) rvisc   = 1.d0                                         !uniform viscosity
c
c     if(lmass) 
c    +rvisc=(0.89+4.4e-3*sqrt(tmp)-1.2e-3*tmp+1.5e-4*tmp**1.5)/0.89  !Na persulphate from data sheet at 25C
c    +rvisc=1.d0  !uniform viscosity

      return
      end
c     ********************************************************************
c
c     relative density change:
c     ===========================
c
      function rden(tmp,conc,gamma,lheat,lmass,lage) 
      implicit real*8(a-h,o-z)
      logical lheat,lmass,lage
c
      rden=0.d0
      if(lheat) 
c    + rden=0.0d0         !hardwire for Lunardini
     + rden = (-.6562d-5 + .2166d-7*(tmp-4.d0))*(tmp-4.d0)**2           !best for 0-100 degrees:
c    +rden  = (-.4350d-5 + .8380d-8*(tmp-4.d0))*(tmp-4.d0)**2           !original  for 0-400 degrees
c    +rden   = (1035. - 1415.e-4*(tmp) - 2384.e-6*(tmp**2))/1000. - 1.  !best for 0-400C Fontaine et al.
c
c     if(lage) rden = 0.0d0                                                          !hardwire for a uniform density
c
c     if(lmass) rden = 0.0d0                ! no density 
c
c     if(lmass) rden = tmp/1000.d0     ! when using gamma for mass transport + density, set max source conc = true concentration

      if(lmass) rden = rden + gamma*conc      ! new gamma read in, check: sum rden ???
c     if(lmass) rden = gamma*conc      ! new gamma read in, check: sum rden ???
c
c     if(lmass) rden = 0.000671*tmp - 1.5e-6*(tmp)**1.5     !!Na persulphate from data sheet at 25C

      return
      end
c     ********************************************************************
c
c     absolute density:
c     -------------------
      function den(tmp,lheat,lmass,lage)
      implicit real*8(a-h,o-z)
      logical lheat,lmass,lage
c
      den=1000.d0
      if(lheat) 
c    + den=1000.d0      !hardwire for Lunardini
     +    den=1000.d0*(1.+(-.6562d-5+.2166d-7*(tmp-4.d0))*(tmp-4.d0)**2)     !original for 0-100 degrees
c     den  = 1000.d0 * 
c    +         (1.d0+ (-.4350d-5 + .8380d-8*(tmp-4.d0))*(tmp-4.d0)**2)  !original for 0-400 degrees, bad for >200
c    +den   = 1035. - 1415.e-4*(tmp) - 2384.e-6*(tmp**2)                !best for 0-400C Fontaine et al.
c
c     if(lage) den = 1000.d0                                            !c     hardwire constant fluid density
c     if(lmass) den = 1000.d0                                            !c     hardwire constant fluid density
c      
c     if(lmass) den = 1000. + 0.671*tmp - 1.5e-3*(tmp)**1.5   !Na persulphate from data sheeet at 25C
c     if(lmass) den = 1000. + tmp              !generic for cl (e.g. for seawater with cmax = 24 kg/m^3, fluid density = 1000*24 = 24,000 kg/m^3)

      return
      end
c
c *******************************************************************
c  #########################################################
c
c     This routine determines the incidences of the 1D line
c     elements  (incidences on flow grid) 
c     -------------------------------------------------------------
c
      subroutine frac_line(maxfrac,maxne,lvert,inline,in,kf,l,fracdim,
     +               exl,eyl,ezl)
c
      implicit real*8(a-h,o-z)
      real*8 exl(maxne),eyl(maxne),ezl(maxne)
      dimension lvert(maxfrac),inline(4),in(maxne,8)
c
      if(lvert(kf).eq.1) inline(1) = in(l,1)
      if(lvert(kf).eq.2) inline(1) = in(l,1)
      if(lvert(kf).eq.3) inline(1) = in(l,1)
      if(lvert(kf).eq.4) inline(1) = in(l,4)
      if(lvert(kf).eq.5) inline(1) = in(l,4)
      if(lvert(kf).eq.6) inline(1) = in(l,2)
c
      if(lvert(kf).eq.1) inline(2) = in(l,5)
      if(lvert(kf).eq.2) inline(2) = in(l,2)
      if(lvert(kf).eq.3) inline(2) = in(l,4)
      if(lvert(kf).eq.4) inline(2) = in(l,8)
      if(lvert(kf).eq.5) inline(2) = in(l,3)
      if(lvert(kf).eq.6) inline(2) = in(l,3)
c
      if(lvert(kf).eq.1) fracdim = ezl(l)
      if(lvert(kf).eq.2) fracdim = exl(l)
      if(lvert(kf).eq.3) fracdim = eyl(l)
      if(lvert(kf).eq.4) fracdim = ezl(l)
      if(lvert(kf).eq.5) fracdim = exl(l)
      if(lvert(kf).eq.6) fracdim = eyl(l)
c
      return
      end
c  #########################################################
c
c     This routine determines the incidences of the 2D line
c     elements  (incidences on flow grid) 
c     -------------------------------------------------------------
c
      subroutine frac_plane(maxfrac,maxne,lvert,inline,in,kf,l,
     +          fdimx,fdimy,exl,eyl,ezl,rexel,sarea)
c
      implicit real*8(a-h,o-z)
      real*8 exl(maxne),eyl(maxne),ezl(maxne)
      dimension lvert(maxfrac),inline(4),in(maxne,8)
c
      if(lvert(kf).eq.1) then
      inline(1) = in(l,1)
      inline(2) = in(l,4)
      inline(3) = in(l,8)
      inline(4) = in(l,5)
      fdimx = eyl(l)
      fdimy = ezl(l)
      sarea = fdimx*fdimy
      rexel = ezl(l)/eyl(l)
      endif
      if(lvert(kf).eq.2) then
      inline(1) = in(l,2)
      inline(2) = in(l,3)
      inline(3) = in(l,7)
      inline(4) = in(l,6)
      fdimx = eyl(l)
      fdimy = ezl(l)
      sarea = fdimx*fdimy
      rexel = ezl(l)/eyl(l)
      endif
      if(lvert(kf).eq.3) then
      inline(1) = in(l,1)
      inline(2) = in(l,2)
      inline(3) = in(l,6)
      inline(4) = in(l,5)
      fdimx = exl(l)
      fdimy = ezl(l)
      sarea = fdimx*fdimy
      rexel = ezl(l)/exl(l)
      endif
      if(lvert(kf).eq.4) then
      inline(1) = in(l,4)
      inline(2) = in(l,3)
      inline(3) = in(l,7)
      inline(4) = in(l,8)
      fdimx = exl(l)
      fdimy = ezl(l)
      sarea = fdimx*fdimy
      rexel = ezl(l)/exl(l)
      endif
      if(lvert(kf).eq.5) then
      inline(1) = in(l,1)
      inline(2) = in(l,2)
      inline(3) = in(l,3)
      inline(4) = in(l,4)
      fdimx = exl(l)
      fdimy = eyl(l)
      sarea = fdimx*fdimy
      rexel = eyl(l)/exl(l)
      endif
      if(lvert(kf).eq.6) then
      inline(1) = in(l,5)
      inline(2) = in(l,6)
      inline(3) = in(l,7)
      inline(4) = in(l,8)
      fdimx = exl(l)
      fdimy = eyl(l)
      sarea = fdimx*fdimy
      rexel = eyl(l)/exl(l)
      endif
c
      return
      end
c
c *************************************************
c
c *************************************************
c**********************************************************************
c    gleichungsloeser fuer  a*x = b  mit konjugierten gradienten
c    vorkonditionierung durch unvollstaendige lu-zerlegung
c***********************************************************************
c
      subroutine precg (a,n, iaa,ind, b,x, aa,laa, neu)
      implicit real*8(a-h,o-z)
c
c                                            d.braess     august 1988
c   schneller loeser von a*x = b             ruhr-universitaet bochum
c   a wird in kondensierter form erwartet und nicht zerstoert
c   ind enthaelt positionen der diagonalelemente von a
c   aa ist hilfsspeicher der laenge laa.ge. 3*n+dim(a)
c   neu = 0:  neue matrix, startwert x=0, b und x duerfen identisch sein
c       = 1:  neue matrix, startwert fuer x vorgegeben
c       = 2:  alte matrix, startwert fuer x vorgegeben
c
c     aufteilung des hilfsfeldes: aa(kmax), g(n), h(n), d(n)
c********************************************************************
c
c                   hinweise zum programm cgicc
c                 ===============================
c
c   das programm steht den universitaeten und staatl. forschungsinstitut
c   kostenlos zur verfuegung. andere benutzer muessen eine lizenz
c   erwerben, wenn sie das programm ueber einen test hinaus verwenden.
c
c       warnungen
c       ---------
c
c   1) pcg erzielt in ... schritten
c      nur eine relative genauigkeit von ...
c
c      --  wenn diese meldung allein erscheint, reicht die vorkondi-
c          tionierung nicht aus, um die gewuenschte genauigkeit in der
c          angegebenen schrittzahl zu erzielen. man setze die angegebene
c          schrittzahl 'iter' hoch oder waehle ein anderes verfahren.
c
c   2) genauigkeit in der metrik der vorkonditionierung schon gut.
c
c      --  die genauigkeit ist um den faktor 'genau/100' reduziert,
c          wenn man sie in der metrik der vorkonditionierung misst.
c          es wird abgebrochen, um das iterative verbessern im bereich
c          der rundungsfehler zu vermeiden.
c
c   3) werte von beta bis ...
c      weisen auf eine instabile vorkonditionierung hin.
c
c      --  bei der berechnung der konjugierten richtungen treten kleine
c          nenner auf. die vorkonditionierung ist knapp ausreichend.
c
c   4) hinweis: bei der ilu-zerlegung wird eine reduktion
c      der diagonalelemente um den faktor ... verlangt.
c
c      --  durch die (unvollstaendige) choleski-zerlegung werden die
c          diagonalelemente reduziert. wenn der reduktionsfaktor nahe
c          unter 1 liegt, weist das auf eine grosse konditionszahl hin.
c          zur stabilisierung wird die reduktion gedaempft. evtl. werden
c          viele iterationsschritte benoetigt.
c
c      fehlermeldungen
c      ---------------
c   5) pcg mit icc erwartet ein hilfsfeld der laenge ...
c      hilfsspeicher ist zu klein!
c
c      --  das programm erfordert ein hilfsfeld der laenge 3n + anzahl
c          der effektiven matrixelemente.
c
c   6) pcg wird eine matrix mit ...-tem diagonalelement < 0 vorgelegt!
c
c      --  der algorithmus ist fuer positiv definite matrizen vorgesehen
c          und die diagonalelemente muessen positiv sein.
c
c   7) die gegebene matrix erweist sich im schritt ... als nicht definit
c
c      --  die auswertung mit dem korrekturvektor liefert einen
c          nicht positiven wert von d'ad. die matrix ist nicht positiv
c          definit, oder die auswertung ist anfaellig gegenueber
c          rundungsfehlern.
c
c   8) schrittweite zu gross!
c
c      --  die auswertung von d'ad fuehrt auf kleine nenner. die rechnun
c          ist anfaellig gegenueber rundungsfehlern.
c
c   9) nach halber maximaler schrittzahl noch keine echte reduktion!
c
c      --  es ist nicht zu erwarten, dass in der zugelassenen schrittzah
c          die gewuenschte genauigkeit erreicht wird.
c          die vorkonditionierung reicht nicht aus.
c
c   weiter auskuenfte erhalten sie vom aufsteller (bitnet: p150204@dboru
c********************************************************************
      integer n, neu
      integer ind(1),iaa(1), laa,laamin,kmax, kg,kh,kd
      real*8    a(1),aa(1), b(n),x(n)
c
      kmax   = ind(n+1) - 1
      laamin = 3*n + kmax
      if (laa .lt. laamin) go to 80
      kg   = kmax + 1
      kh   = kg + n
      kd   = kh + n
      call cgicc (a,n,iaa,ind, b,x, aa,aa(kg),aa(kh),aa(kd), neu)
      go to 90
   80 write (6,86) laamin
   86 format (' pcg mit icc erwartet ein hilfsfeld der laenge',
     2        i6 / ' hilfsspeicher ist zu klein!')
   90 return
      end
c------------------------------------------------------------
      subroutine cgicc (a,n,iaa,ind, b,x, aa, g,h,d, neu)
c------------------------------------------------------------
c   allgemeines cg-verfahren mit vorkonditioierung
c   hier speziell fuer icc-zerlegung, soonst ersetze icc/reicc
c
c   g0 = hx0 - b,  h0 = c**-1g0,  d0 = - h0
c   tau = gh / dhd
c   x' = x + tau * d
c   g' = g + tau * hd
c   h' = c**-1 * g'
c   beta = g'h' / gh
c   d' = - h' + beta * d
c
      implicit real*8(a-h,o-z)
c
      integer n,i,k,iter,      iaa(1), ind(1)
      real*8    a(1),b(n),x(n),  aa(1), g(n),h(n),d(n)
      real*8    gh,ghalt,dhd,gg, gg1,gh1,sqgg, tau,beta, genau,genau2
c     dimension df(200)
c
c     amax1(x,y)=max(x,y)
c     amin1(x,y)=min(x,y)
c     sqrt(x)=   dsqrt(x)
c
c     parameter
      genau  = 1.e-7
      genau2 = genau**2
      iter = 400
c     startwerte
      bb     = 0.
      if (neu .ne. 0) call mavmuk (a,n,iaa,ind, x,g)
      do  5 i=1,n
      bb   =   bb + b(i)*b(i)
      if (neu .ne. 0) go to 3
      g(i) = - b(i)
      x(i) =   0.
      go to 4
    3 g(i) =   g(i) - b(i)
    4 h(i) =   g(i)
      d(i) =   0.
      k = ind(i)
      if (a(k) .gt. 0.) go to 5
      write (6,81) i
   81 format(' pcg wird eine matrix mit',i6,
     2       '-tem diagonalelement < 0 vorgelegt!')
      stop
c      go to 99
    5 continue
c
c     kontrolle und aufruf der unvollstaendigen zerlegung
      if (neu .gt. 1) go to 7
      k = ind(n+1) - 1
      do 6 i=1,k
      aa(i) = a(i)
    6 continue
      call  icc (aa,n, iaa,ind)
c
c     iteration
    7 gh   =   0.
      beta =   0.
      betamx = 0.
      do 60 k=1,iter
      kk = k
c     jetzt wird h = c**-1 * g
      call reicc (aa,n, iaa,ind, h)
      ghalt = gh
      gh    = 0.
      gg    = 0.
      do 14 i=1,n
      gh = gh + g(i)*h(i)
      gg = gg + g(i)*g(i)
   14 continue
      if (k .eq. 1)   gg1  = gg
      if (k .eq. 1)   gh1  = gh
      if (k .ne. 1)   beta = gh / ghalt
      betamx =  dmax1(beta,betamx)
      do 16 i=1,n
      d(i) = - h(i) + beta * d(i)
   16 continue
c     jetzt wird h = a*d
      call mavmuk (a,n, iaa,ind, d,h)
      dhd = 0.
      do 20 i=1,n
      dhd = dhd + d(i)*h(i)
   20 continue
c     +++++++++++++  die naechsten 6 statements bei kontrolle aktivieren
c        sqgg = sqrt (gg)
c        df(k)   = gh * gh / dhd
c        if (k .eq. 1) write (6,22) n, ind(n+1), iter, genau
c        write (6,21)  k, gh, df(k), sqgg
c  21    format (' ',i3,1p3e12.3)
c  22    format (' n:',i5,' groesse:',i6,' max. schrittzahl:',i4,
c    2           ' genauigkeit:',1p1e9.1 /
c    3           '          gh         df       gg**1/2')
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (dhd .ge. 1.e-6*gh) go to 25
      write (6,24)
   24 format (' schrittweite zu gross!')
      go to 90
   25 if (dhd .gt. 0.) go to 27
      write (6,26) k
   26 format (' die gegebene matrix erweist sich in schritt', i4,
     2        ' als nicht definit!')
      go to 90
   27 tau = gh / dhd
      do 30 i=1,n
      x(i) = x(i) + tau * d(i)
      g(i) = g(i) + tau * h(i)
      h(i) = g(i)
   30 continue
      if (gg .le. dmin1(.01*gg1, genau2*bb)) go to 99
      if (gg .le. .0001*genau2*bb)           go to 99
      if (gh .gt. 0.0001*genau2 * gh1)       go to 55
      write (6,53)
   53 format (' genauigkeit in der metrik der vorkonditionierung',
     2        ' schon gut')
      go to 90
   55 if (2*k.lt.iter .or. gg.lt.0.25*gg1)   go to 60
      write (6,58)
   58 format
     2 (' nach halber max. schrittzahl noch keine echte reduktion!')
      go to 90
   60 continue
   90 if (gg .le. genau2 * bb)               go to 99
      sqgg = sqrt (gg/bb)
      write (6,87) kk, sqgg
   87 format (' pcg erzielt in', i4,' schritten'/
     2        ' nur eine relative genauigkeit im defekt von',1p1e11.2)
      if (betamx .gt. 20.) write (6,88) betamx
   88 format (' werte von beta bis',1pe11.3 /
     2        ' weisen auf eine instabile vorkonditionierung hin')
   99 continue
c     do 91 i=1,kk-1
c     k = kk-i
c  91 df(k) = df(k) + df(k+1)
c     do 92 k=1,kk
c  92 df(k) = sqrt (df(k)-df(kk))
c     write (6,93) (k, df(k), k=1,kk)
c  93 format (5(i4,1p1e11.3))
      return
      end
c----------------------------------------------------------------------
      subroutine icc (aa,n, iaa,indaa)
      implicit real*8(a-h,o-z)
      integer n, iaa(1),indaa(1), i,j,k,kugr,kogr,kj,k1,m,mm,l,ll,ii
      real*8    aa(1), diared,dimare, s,hv
c     amax1(x,y)=max(x,y)
c     amin1(x,y)=min(x,y)
c----------------------------------------------------------------------
c   unvollstaendige lu-zerlegung
c
      diared = 0.
      dimare = 0.9
      do 60 j=1,n
      kj = indaa(j)
      kugr = kj+1
      kogr = indaa(j+1)-1
      s = 0.
      if (kugr .gt. kogr) go to 40
      k = kogr + 1
      do 20 k1=kugr,kogr
      k = k-1
      i = iaa(k)
      s = 0.
c
c     aufsummation der skalarprodukte fuer nicht diagonale
      m  = kogr
      mm = indaa(i+1) - 1
   12 l  = iaa(m)
      if (l .ge. i) go to 18
   14 ll = iaa(mm)
      if (ll.ge. i) go to 18
   15 if (l .ne.ll) go to 16
c
      s = s + aa(m) * aa(mm)
c
      m  = m - 1
      mm = mm- 1
      go to 12
   16 if (l .ge. ll) go to 17
      m  = m - 1
      l  = iaa(m)
      if (l .ge. i)  go to 18
      go to 15
   17 mm = mm - 1
      go to 14
   18 continue
      aa(k) = aa(k) - s
   20 continue
c
c     skalarprodukt fuer diagonalelemente
c     renormierung oberhalb der diagonale
      s = 0
      do 30 k=kugr,kogr
      i  = iaa(k)
      ii = indaa(i)
      hv = aa(k)
      aa(k) = hv / aa(ii)
   30 s  = s + aa(k) * hv
   40 continue
c
c     neues diagonalelement
      diared = dmax1 (diared, s / aa(kj))
      s      = dmin1 (0.98*s, 0.5 * (s + dimare * aa(kj)))
      aa(kj) = aa(kj) - s
   60 continue
      if (diared .ge. dimare)  write (6,95) diared
   95 format (' hinweis: bei der ilu zerlegung wird eine reduktion' /
     2        ' der diagonalelemente um den faktor', f7.2,' verlangt.')
   99 return
      end
c
c
      subroutine reicc (aa,n,iaa,indaa, z)
      implicit real*8(a-h,o-z)
c
c   multiplikation mit inverser aus icc-zerlegung
c
      integer n, iaa(1),indaa(1), i,j,k,kugr,kogr
      real*8    aa(1),z(n), s
c     multiplikation mit l**-1
      do 30 i=1,n
      kugr = indaa(i)   + 1
      kogr = indaa(i+1) - 1
      if (kugr .gt. kogr) go to 30
      s = 0.
      do 20 k=kugr,kogr
      j = iaa(k)
   20 s = s + aa(k) * z(j)
      z(i) = z(i) - s
   30 continue
c     multiplikation mit d*--1
      do 40 j=1,n
      k = indaa(j)
   40 z(j) = z(j) / aa(k)
c     multiplikation mit r**-1
      do 60 jj=1,n
      j = n + 1 - jj
      kugr = indaa(j)    + 1
      kogr = indaa(j+1)  - 1
      if (kugr .gt. kogr) go to 60
      do 50 k=kugr,kogr
      i = iaa(k)
   50 z(i) = z(i) - aa(k) * z(j)
   60 continue
      return
      end
c
c
      subroutine mavmuk (a,n, iaa,inda, v,z)
      implicit real*8(a-h,o-z)
c
c   multiplikation von kondensierter matrix mit vektor
c
      integer n, iaa(1),inda(1), i,j,k,kugr,kogr
      integer m,mj
      real*8    a(1),v(n),z(n)
c
      do 10 i=1,n
      z(i) = 0.
   10 continue
      do 30 j=1,n
      mj   = inda(j)
      z(j) = z(j) + a(mj)*v(j)
      kugr = mj + 1
      kogr = inda(j+1)- 1
      if (kugr .gt. kogr) go to 30
      do 20 k=kugr,kogr
      i = iaa(k)
      z(j) = z(j) + a(k) * v(i)
      z(i) = z(i) + a(k) * v(j)
   20 continue
   30 continue
      return
      end
c ------------------------------------------------------------------
c
c  NEW SOLVER added June 2005
c
c ------------------------------------------------------------------
c
      subroutine precgn (a,n, iaa,ind, b,x, aa,laa, kblock,lblo, neu,
     +                   maxnn)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c    preconditioned conjugate gradient solver for a*x = b - conjugate
c    gradient preconditioning through incomplete lu-decomposition
c***********************************************************************
c
c
c                                            d.braess      sept. 1993
c   fast solver for [a]*{x} = {b}            ruhr-universitaet bochum
c   [a]   expected in condensed form, will not be destroyed
c   a     wird in kondensierter form erwartet und nicht zerstoert
c   ind: contains positions of diagonal elements of [a]
c   ind  enthaelt positionen der diagonalelemente von a
c   [aa]   additional storage array of length laa.ge. 6*n
c   aa     ist hilfsspeicher der laenge laa.ge. 6*n
c   kblock ist hilfsspeicher der laenge lblo
c          empfehlung  lblo = n
c
c   new = 0:  new matrix, start value x=0, b and x may be identical
c       = 1:  new matrix, defined start value for x 
c       = 2:  old matrix, defined start value for x
c   neu = 0:  startwert x=0, b und x duerfen identisch sein
c       = 1:  startwert fuer x vorgegeben
c
c     sharing of additional storage array: aa(kmax), g(n), h(n), d(n)
c********************************************************************
c
c                    program cgicc - notes
c                 ===============================
c   the program is available for universities and government research-
c   institutions for free.  other users have to purchase a license, if
c   they use the program exceeding tests. 
c
c       warnings
c       --------
c
c   1) pcg used ... steps
c      for a relative accuracy of only ...
c
c      --  if this warning appears exclusively, preconditioning is not
c          sufficient to reach the desired accuracy using the defined 
c          number of (iteration) steps. 
c          remedy: increase the defined number of iteration steps
c          ('iter') or choose a different solver.
c
c   2) accuracy in ?metrik? of preconditioning already good
c
c      --  the accuracy is reduced by the factor 'genau/100',
c          if measured in ?metrik? of preconditioning.
c          process is stopped, to improve the iterative procedure
c          and to avoid truncation errors.
c
c   3) values from beta to ...
c      indicate instable preconditioning.
c
c      --  during the computation of conjugated directions small
c          denominators arise.  preconditioning is barely sufficient.
c
c   4) note: during ilu-decomposition, a reduction
c      of diagonal elements by the factor ... is required.
c
c      --  due to the (incomplete) choleski-decomposition, diagonal elements
c          get reduced.  a factor of reduction close to 1 indicates
c          a high number of conditioning processes.
c          for stabilization, reduction is dampened.  this might lead to a
c          higher number of iteration steps.
c
c      errors
c      ------
c   5) pcg with icc expects additional storage array of length ...
c      additional storage array too small!
c
c      --  the program requires an additional storage array of length
c          3n + number of effective matrix elements.
c
c   6) pcg has been assigned a matrix with ...-th diagonal element < 0!
c                                                                     
c      --  the algorithm is designated for positive definite matrices.
c          the diagonal elements have to be positive.
c
c   7) step ... determined the matrix as indefinite.
c
c      --  analysis using the correction vector results in a value
c          for d'ad, which is not positive.  
c          the matrix is not positive definite or the analysis is
c          sensitive with respect to truncation errors.
c
c   8) step - spacing too large!
c
c      --   analysis of d'ad results in small denominators, calculation
c           is sensitive with respect to truncation errors.
c
c   9) 50% of maximum number of steps completed - reduction not satisfactory!
c
c      --  it can't be expected, that required accuracy will be reached
c          using the defined number of steps.  preconditioning not
c          sufficient.
c
c   for additional information contact: (bitnet: p150204@dboru)
c**********************************************************************
c    gleichungsloeser fuer  a*x = b  mit konjugierten gradienten
c    vorkonditionierung durch blockweise ssor-relaxation
c    tridiagonale bloecke werden bestimmt
c***********************************************************************
c
c
c     aufteilung des hilfsfeldes: g(n), h(n), d(n), ax(n)
      implicit none
      integer*4  maxnn,lblo
      integer*4  n, neu, ind(*),iaa(*), laa, kblock(lblo)
      integer*4  laamin,kmax, kg,kh,kd,kax,kdia,ksup
      integer*4  nz,irsl
      real*8   a(1),aa(*), b(n),x(n)
c
      kmax   = ind(n+1) - 1
      laamin = 6*n
      if (laa .lt. laamin) go to 80
      kg   = 1
      kh   = kg + n
      kd   = kh + n
      kax  = kd + n
      kdia = kax  + n
      ksup = kdia + n
      call cgssortri (a,n,iaa,ind, b,x, aa(kg),aa(kh),aa(kd),
     & aa(kax),aa(kdia),aa(ksup),kblock,lblo,neu,maxnn)
      go to 90
   80 write (6,86) laamin
   86 format (' pcg with ssor expects additonal storage array of length
     2 ',     i6 / ' additional storage array too small! ')
c   86 format (' pcg mit ssor erwartet ein hilfsfeld der laenge',
c     2        i6 / ' hilfsspeicher ist zu klein')
   90 return
      end
c============================================================
      subroutine cgssortri (a,n,iaa,ind, b,x, g,h,d,ax,
     *          dia,sup, kblock,lblo, neu,maxnn)
c============================================================
c   allgemeines cg-verfahren mit vorkonditionierung
c   hier speziell fuer block-ssor-relaxation
c
c   g0 = hx0 - b,  h0 = c**-1g0,  d0 = - h0
c   tau = gh / dhd
c   x' = x + tau * d
c   g' = g + tau * hd
c   h' = c**-1 * g'
c   beta = g'h' / gh
c   d' = - h' + beta * d
c
      implicit none
      integer*4   itmax
      real*8    genau,genau2
c     parameter zur steuerung des endes der iteration
      parameter (itmax=500, genau=2.e-7, genau2=genau**2)
c
      integer*4  lblo,maxnn
      integer*4 n,iaa(*),ind(1),neu,kblock(lblo)
      real*8    a(*),b(n),x(n),  g(n),h(n),d(n),ax(n),dia(n),sup(n)
      real*8    gh,ghalt,dhd,gg, gg0,gg1,gh1,sqgg, tau,beta,betamx, bb
      integer*4  nz,irsl,iremainder,iadd
      integer*4   i,k,iter
c++++++++++++++++++++++
c     real*8    df(itmax)
c++++++++++++++++++++++
c
c     startwerte und kontrolle der gegebenen matrix
      bb    = 0.
      if (neu .ne. 0) call mavmukn (a,n,iaa,ind, x,g)
      do  5 i=1,n
      bb   =   bb + b(i)*b(i)
      if (neu .eq. 0) then
         g(i) = - b(i)
         x(i) =   0.
      else
         g(i) =   g(i) - b(i)
      end if
      h(i) =   g(i)
      d(i) =   0.
      if (a(ind(i)) .le. 0.) then
         write (*,81) i
         write (6,81) i
   81 format(' pcg has been assigned a matrix with',i6,
     2       '-th diagonal entry < 0 !')
c   81    format(' pcg wird eine matrix mit',i6,
c     2          '-tem diagonalelement < 0 vorgelegt!')
          stop
c         go to 99
      end if
    5 continue
c     bei trivialer rechter seite wird nicht gerechnet
      if (dabs(bb).lt.1.d-20) then
         write (*,83)
         write (6,83)
83       format('precg has a right-hand side entry = 0 : no solution')
c     * ' precg wird eine gleichung mit rechter seite = 0 vorgelegt'
         return
      end if
c     gg0 = gg
c
      call zerltri (a,n, iaa,ind, dia,sup, kblock,lblo)
c----------------
c     iteration
c----------------
      gh   =   0.
      beta =   0.
      betamx = 0.
      do 60 k=1,itmax
      iter = k
c     jetzt wird h = c**-1 * g
      call ssortri (a,n, iaa,ind, g,h, ax, dia,sup, kblock,lblo)
      ghalt = gh
      gh    = 0.
      gg    = 0.
      do 14 i=1,n
      gh = gh + g(i)*h(i)
      gg = gg + g(i)*g(i)
   14 continue
      if (k .eq. 1) then
          gg1  = gg
          gh1  = gh
      else
          beta   = gh / ghalt
          betamx = max (beta,betamx)
      end if
      if (k .eq. 2)  gh1  = gh
      do 16 i=1,n
   16 d(i) = - h(i) + beta * d(i)
c     jetzt wird h = a*d
      call mavmukn (a,n, iaa,ind, d,h)
      dhd = 0.
      do 20 i=1,n
   20 dhd = dhd + d(i)*h(i)
c     +++++++++++++  die naechsten 6 statements bei kontrolle aktivieren ++++
c        sqgg = sqrt (gg)
c        df(k)   = gh * gh / dhd
c        if (k .eq. 1) write (*,22) n, ind(n+1), itmax,genau, sqrt(gg0)
c        write (*,21)  k, gh, df(k), sqgg
c  21    format (' ',i3,1p,3e12.3)
c  22    format (' vorkonditionierung durch block-ssor'/
c    1           ' n:',i6,' groesse:',i6,' max. schrittzahl:',i4,
c    2           ' genauigkeit:',1p,1e9.1 /
c    3           '          gh         df        gg**1/2' / 28x,e12.3)
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (dhd .le. 0.) then
          write (*,26) k
          write (6,26) k
   26 format (' step ', i4,' determined the matrix to be indefinite!')
c   26     format (' pcg: die gegebene matrix erweist sich in schritt',
c     2              i4,' als nicht definit!')
          stop
c          go to 90
      end if
      if (dhd .lt. 1.e-6*gh) then
          write (*,*) ' pcg: schrittweite zu gross!'
          go to 90
      end if
c
      tau = gh / dhd
      do 30 i=1,n
      x(i) = x(i) + tau * d(i)
      g(i) = g(i) + tau * h(i)
      h(i) = g(i)
   30 continue
      if (gg .le. min (.01*gg1, genau2*bb)
     2   .or. gg .le. .0001*genau2*bb)           go to 99
      if (gh .lt. 0.00001*genau2 * gh1) then
         if (gh .lt. -genau2 * gh1) then
           write(*,*) ' preconditioning is not good - problem to tough'
c             write (*,*) ' pcg: vorkonditionierung ist nicht definit'
         else
           write (*,*) ' pcg - warning:',
     *   ' matrix already preconditioned - input data is probably wrong'
c             write (*,*) ' pcg - warnung:',
c     *   '  genauigkeit in der metrik der vorkonditionierung schon gut'
         end if
         go to 90
      end if
      if (2*iter.gt.itmax .and. gg.gt.0.25*gg1) then
         write (*,*)
     2 ' pcg: nach halber max. schrittzahl noch keine echte reduktion!'
         go to 90
      end if
   60 continue
c-------------
   90 if (gg .ge. genau2 * bb) then
         sqgg = sqrt (gg/bb)
         write (*,87) iter, sqgg
         write (6,87) iter, sqgg
   87 format (10x,'pcg used ', i4,' steps'/10x,
     2        'for a relative accuracy of only ',1p1e11.2)
c  87    format (' pcg erzielt in', i4,' schritten'/
c     2        ' nur eine relative genauigkeit im defekt von',1p,1e11.2)
c        if (betamx .gt. 50.) write (*,88) betamx
c  88 format (' values from beta to',1pe11.3 /
c    2        ' indicate unstable preconditioning')
c   88    format (' werte von beta bis',1p,e11.3 /
c     2        ' weisen auf eine instabile vorkonditionierung hin')
      end if
   99 continue
c      write(6,100)iter
c100   format(10x,'pcg solver used ',i5,' iterations')
c     do 91 i=1,kk-1
c     k = kk-i
c  91 df(k) = df(k) + df(k+1)
c     do 92 k=1,kk
c  92 df(k) = sqrt (df(k)-df(kk))
c     write (*,93) (k, df(k), k=1,kk)
c  93 format (5(i4,1p,1e11.3))
      return
      end
c
c=============================================
      subroutine mavmukn (a,n, iaa,inda, v,z)
c=============================================
c
c   multiplikation von kondensierter matrix mit vektor
c
      implicit none
      integer*4 n
      real*8   a(*),v(n),z(n)
      integer*4  iaa(*),inda(1)
      integer*4  k,kugr,kogr, i,j
      real*8   vj,zj
c
c     wegen i<j ist keine initialisierung noetig, kopiere n
      kogr = inda(1) - 1
      do 30 j=1,n
      kugr = kogr + 2
      kogr = inda(j+1)- 1
c     diagonale
      vj   = v(j)
      zj   = a(kugr-1) * vj
      do 20 k=kugr,kogr
      i = iaa(k)
      zj   = zj   + a(k) * v(i)
      z(i) = z(i) + a(k) * vj
   20 continue
      z(j) = zj
   30 continue
      return
      end
c=============================================
      subroutine ssortri (a,n, iaa,ind, b,x, ax,dia,sup, kblock,lblo)
c     blockssor mit tridiagonalen bloecken
c     individuelle blockgroesse wird bestimmt
c=============================================
c
c   1 schritt der symmetrischen overrelaxation zu a * x = b
c   wobei bloecke simultan behandelt werden,
c   solange dies ueber tridiagonale matrizen geht
c   b und x duerfen nicht identisch sein
c   ax ist hilfsvektor der laenge n -- hinterher: dx/omega
c
      implicit none
      real*8    omega
      parameter (omega=1.3d0)
c
      integer*4   n, iaa(*), kogr, jogr
      integer*4   ind(1), mi,k
      integer*4   j,i, nschicht, jj, lblo,kblock(lblo), kkb,ib
      real*8    a(*), b(n),x(n),ax(n), dia(n),sup(n),  s,xj
      save      kkb
c
c     vorwaertsiteration
      kogr = ind(1)-1
      do 40 ib=1,kkb
      jj       = kblock(ib)
      jogr     = kblock(ib+1)-1
      nschicht = kblock(ib+1)-jj
      do 32 j=jj,jogr
      mi   = kogr + 1
      kogr = ind(j+1) - 1
      s    = b(j)
      do 30 k=mi+min(2,j-jj+1), kogr
   30 s = s - a(k)*x(iaa(k))
      ax(j) = s
   32 continue
c     im block vorwaerts
      x(jj) = ax(jj)
      do 34 j=jj+1,jogr
      x(j)  = ax(j) - sup(j) * x(j-1)
   34 continue
c     im block rueckwaerts
      x(jogr) = x(jogr) * dia(jogr)
      do 36 j=jogr,jj+1,-1
   36 x(j-1) = x(j-1) * dia(j-1) - sup(j) * x(j)
   40 continue
c
c     rueckwaertsiteration
      mi = ind(n+1)
      do 50 ib=kkb,1,-1
      jj       = kblock(ib)
      jogr     = kblock(ib+1)-1
      nschicht = kblock(ib+1)-jj
c     im block vorwaerts
      x(jj)   = ax(jj)
      do 44 j=jj+1,jogr
   44 x(j)    = ax(j) - sup(j) * x(j-1)
c     im block rueckwaerts
      x(jogr) = x(jogr) * dia(jogr)
      do 46 j=jogr,jj+1,-1
   46 x(j-1) = x(j-1) * dia(j-1) - sup(j) * x(j)
      do 48 j=jogr,jj,-1
      kogr = mi - 1
      mi   = ind(j)
      xj   = x(j)
      do 48 k=mi+min(2,j-jj+1), kogr
      i = iaa(k)
      ax(i) = ax(i) - a(k)*xj
   48 continue
   50 continue
      return
c     -----------------------------------------
      entry zerltri (a,n, iaa,ind, dia,sup, kblock,lblo)
c     cholesky-zerlegung der diagonalen bloecke
c     -----------------------------------------
c
c     aufspaltung in tridiagonale diagonalbloecke
      jj  = 0
      kkb = 0
      do 14 j=1,n
      mi = ind(j)
      if (ind(j+1).eq.mi+1 .or. iaa(mi+1).ne.j-1 .or.
     &   (ind(j+1).gt.mi+2 .and.iaa(mi+2).ge.jj)) then
c        neuer block
         kkb = kkb + 1
         jj  = j
         kblock(kkb) = j
         if (kkb .ge. lblo) then
            write (*,*)
     *    ' cgtriadapt: dimension fuer zahl der bloecke zu klein!'
            stop
         end if
      end if
      dia(j) = a(mi)
      sup(j) = a(mi+1)
   14 continue
      kblock(kkb+1) = n+1
c
      do 24 ib=1,kkb
      jj      = kblock(ib)
      jogr    = kblock(ib+1)-1
      dia(jj) = 1. / dia(jj)
c     sup(jj) = 0.
      do 24 j=jj+1, jogr
      s      = sup(j)
      sup(j) = sup(j) * dia(j-1)
      dia(j) = 1. / (dia(j) - sup(j)*s)
   24 continue
c     overrelaxtion
      do 26 j=1,n
   26 dia(j) = omega * dia(j)
      return
      end
c__________________________________________________________________________
c
c *******************************************************************
c
c trans.for
c
c    ******** subroutine  transport 2  **********   Original thermal R  version
c
c  three-dimensional advective-dispersive transport
c  linear isoparametric quadrilateral elements
c
c     e.o. frind, j.w. molson
c     university of waterloo
c     (c) 1993
c
c   *******************************************************************
c
c       copyright 1987 emil o. frind, 1988 j.w. molson
c
c       duplication of this program or any part thereof without
c       the express written consent of the author is prohibited.
c
c   *******************************************************************
c
      subroutine trans2(maxn,maxnn,maxne,nf,maxna,laa,nw,
     + x,y,z,u0,u1,u2,fc,fb,f,fs,icflow,ic,lc,vx,vy,vz,inflow,in,map,
     + a,aa,iaa,ind,ib,alh,alv,ath,atv,dd,kdisp,n,nbb,dt1,dt2,nz,
     + ne,nn,wp,wa,wp1,wa1,pq,tq,cq,por,sw,ag,hag,ttemp,ara,leak2,
     + hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,htmt,kint,tast,tpcgt,tdays,
     + exl,eyl,ezl,wlh,pp,qq,modelwu,vt,tcs,bz,bzf,cpsm,cpi,cf,pi,sat,
     + nef,
     + fbf,rinc,htmf,ifracl,lvert,nfrac,xarea,nex,ney,nez,vlin,maxfrac,
     + tclw,tcli,tclm,ifdim,hgz2,hx,hy,hh,ht,wx,wy,vx2d,vy2d,
     + cot0,cot1,tkl,numis,flux,fluxv,ifl,maxfx,rtrans,decay,decay2,
     + nflux,agefx,rg,rg2df,porsurf,lmass,lheat,lage,ts,rtdiff,
     + p2,q2,ts2,omega2,modelwu2,modelkr2,depc,vexp,smax,spn0,spn1,
     + ppn,qqn)
c
      implicit real*8(a-h,o-z)
      external wu,dwu
      real*8 start,finish
c
c     constant dimensions
c     direct element matrices:
c     -------------------------
      dimension hxx(8,8),hyy(8,8),hzz(8,8),hxy(8,8),hxz(8,8),hyz(8,8),
     + hvx(8,8),hvy(8,8),hvz(8,8),htmt(8,8),htmf(8,8),
     + rg(8,8),rg2df(4,4)
c
      dimension re(8,8),ro(8,8),inl(8),kb(7),ag(9),hag(9)
      logical leak2,lmass,lheat,lage
c
c     ne = number of elements
c     nn = number of nodes
c     n  = number of degrees of freedom
c     na = total number of non-zero matrix entries in condensed matrix
c        = 14*n (approximately)
c     nw = width of matrix ib ( = 15 for prismatic grid )
c     nf = maximum number of nodes on one face
c     laa = 3*n + na
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     variable dimensions:
c     arrays of size nn
c     -----------------
      real*8 tq(maxnn),pq(maxnn),cq(maxnn)
      real*8  x(maxnn),y(maxnn),z(maxnn),por(maxne),sw(maxne)
      real*8 u0(maxnn),u1(maxnn),fc(maxnn),fb(maxnn),fs(maxnn)
      dimension  f(maxnn),u2(maxnn),fbf(nf)
      integer   ic(maxnn),icflow(maxnn),lc(maxnn),map(maxnn)
c
c     1D line or 2D plane elements - transport 
c     -----------------------------------------
      dimension rel(4,4),rol(4,4),ifracl(maxfrac),inline(4)
      dimension lvert(maxfrac),ifdim(maxfrac)
      real*8 vlin(maxfrac),vx2d(maxfrac),vy2d(maxfrac),xarea(maxfrac)
      dimension hx(4,4),hy(4,4),hh(4,4),ht(4,4),wx(4,4),wy(4,4),hgz2(4)
c
c     arrays of size ne, (ne,8)
c     ------------------------
      real*8   vx(maxne),vy(maxne),vz(maxne)
      real*8   exl(maxne),eyl(maxne),ezl(maxne)
      real*8   vt(maxne),tclm(maxne),cot0(maxne),cot1(maxne)
      real*8   tkl(maxne),cpsm(maxne)
      real*8    depc(maxne),vexp(maxne),decay(maxne),smax(maxne)
      integer  inflow(maxne,8),in(maxne,8)
      real*8   pp(maxne),qq(maxne),ppn(maxnn),qqn(maxnn)
      real*8   spn0(maxnn),spn1(maxnn)
c
c     arrays for internal heat source
c     --------------------------------
      dimension flux(maxfx),fluxv(maxfx),ifl(maxfx),nflux(maxfx,8)
c
c     array of size nf
c     ----------------
      dimension ara(nf)
      dimension bz(nf),bzf(nf)
c
c     arrays a(na),iaa(na),ind(n+1),ib(n,15),aa(3n+na)
c            for conjugate gradient solver only
c     ------------------------------------------------
      dimension a(maxna),aa(laa)
      integer iaa(maxna),ind(maxnn),ib(maxnn,nw)
c
c     function defining -absolute- water density as a function of t
c     for use in latent heat term, thermal parameters
c     --------------------------------------------------------------
c
      nx=nex+1
      ny=ney+1
      nexy = nex*ney

c     initialize arrays,
c     condense non-zero flux array:
c     -------------------------------
      do 45 i=1,nn
      fs(i)=0.d0
      f(i)=0.d0
      if(ic(i).eq.1) goto 45
      ii=i-lc(i)
      f(ii)=fb(i) 
  45  continue
c
      do 247 i=1,maxna
  247 a(i)=0.d0
c
c     tsurfmin=999.
c     tsurfmax=-999.
c     fcmin=999.
c     fcmax=-999.
c     u0min=999.
c     u0max=-999.
c     do i=1,nn 
c     tsurfmin = min(tsurfmin,u2(i))
c     tsurfmax = max(tsurfmax,u2(i))           
c     fcmin = min(fcmin,fc(i))
c     fcmax = max(fcmax,fc(i))           
c     u0min = min(u0min,u0(i))
c     u0max = max(u0max,u0(i))           
c     enddo
c
c     write(6,3376) tsurfmin,tsurfmax,fcmin,fcmax,u0min,u0max
c3376 format(/10x,'in trans2 min,max temps,fc,u0: ',6e15.5)
c     _________________________________________
c
c     matrix generation loop
c     loop over elements:
c     =========================================
c     =========================================
c
c     call clock@(start)
c
c     write(6,903)
c903  format(56x,'l     tavg       por        sw           ww1'
c    +'      wi1          wm         wu1       dwu1      dlh1'
c    +'      Co       tk')

      do 100 l=1,ne
c
      do 101 i=1,8
  101 inl(i)=in(l,i)
c
c     average elemental temperature:
c     tav0: avg. at start of time step / tav1: avg. at latest iteration
c     -----------------------------------------------------------------
      tav0=0.d0
      tav1=0.d0
      do 103 i=1,8
      tav0=tav0+u0(inl(i))
      tav1=tav1+u1(inl(i))
 103  continue
      tav0=tav0/8.d0
      tav1=tav1/8.d0
      tavg=(tav0+tav1)/2.0
c     tavg=tav2
c
c     fixed retardation = rtrans for mass transport or age only
c     ----------------------------------------------------------
c     if(lmass.or.lage) then

c     cot0(l)=rtrans        !for mass transport, cot0 and cot1 represent retardation, and go on rhs 
c     cot1(l)=rtrans
c     cotl0=rtrans
c     cotl1=rtrans
c     tccx=0.d0
c     tccy=0.d0
c     tccz=0.d0

c     wwavg=1.d0            !check for other terms below for heat that need to be defined for transport
c     cpfavg=1.d0
c     cotlavg = (cotl0+cotl1)/2.d0
c     rett=cotlavg 
c     else 
c
c     identify temp.-dependent thermal parameters for sat/unsat zone:
c     and elemental latent heat contribution term
c     thermal K tccx adjusted for water,ice,solids fractions
c     need to adjust cpi for temperature-dependence
c     cpi=c_i * rho_i = 2000.*916.7
c     cpm=c_m * rho_m = ...
c     sqrt-mean thermal K: see Roy et al. 1981 (in Mottaghy & Rath GJI 2006)
c     Non-linear thermal R ... 
c     -----------------------------------------------------------------------
      porsw = por(l)*sw(l)
      ppl=pp(l)
      qql=qq(l)
      wu0 = wu(tav0,ppl,qql,modelwu,ts)
      wu1 = wu(tav1,ppl,qql,modelwu,ts)
      ww0=porsw * wu0
      ww1=porsw * wu1
      wi0=porsw * (1.d0-wu0)
      wi1=porsw * (1.d0-wu1)
      wm=(1.d0-por(l))  

      den0 = den(tav0,lheat,lmass,lage)
      den1 = den(tav1,lheat,lmass,lage)
      cpf0=cf*den0
      cpf1=cf*den1
      cpe0=ww0*cpf0 + wi0*cpi + wm*cpsm(l)
      cpe1=ww1*cpf1 + wi1*cpi + wm*cpsm(l)

c      thermal conductivity tk trans2
c      ------------------------------
       tkl(l) = (ww1*sqrt(tclw)+wi1*sqrt(tcli)+wm*sqrt(tclm(l)))**2     !more physically realistic
c      tkl(l) = ww1*tclw + wi1*tcli + wm*tclm(l)                       !hardwire Interfrost ! geometric mean
       tk = tkl(l)           

c     dlh=por(l)*den(tav1,lmass,lage)*wlh*dwu(tav1,p,q,modelwu,ts)
      dwu0 = dwu(tav0,ppl,qql,modelwu,ts)
      dwu1 = dwu(tav1,ppl,qql,modelwu,ts)
      dlh0=por(l)*pi*wlh*dwu0            !pi (ice density) 
      dlh1=por(l)*pi*wlh*dwu1
      cot0(l)=(cpe0+dlh0)                !heat capacity 
      cot1(l)=(cpe1+dlh1)                
      cotl0=cot0(l)
      cotl1=cot1(l)
c      cotlavg = (cotl0+cotl1)/2.d0

c     wwavg = time-centred porosity*Sw*Wu
c     cpfavg = time-centred cfpf
c     tccx   = thermal conductivity lambda
c     November 2013: keep tk on left hand side, keep Co on right hand side
c     uses average ww and cpf for advective term
c     ------------------------------------------------------------------------------------------------------
      wwavg = (ww0+ww1)/2.d0
      cpfavg = (cpf0+cpf1)/2.d0
      cotlavg=((cpe0+cpe1)/2.)/(wwavg*cpfavg)
      rett = cotlavg
      tccx = tk/(wwavg*cpfavg)                 
      tccy = tccx
      tccz = tccx
c
c hardwire for Lunardini tests - fixed parameters in 3 zones
c ------------------------------------------------------------
c  Lunardini parameters
c     tk = 5.0182e-6                                        !ice
c     if(tavg.ge.-4. .and. tavg .le.0.) tk = 1.6512e-7      !mushy
c     if(tavg.ge.0.) tk=3.5030e-6                           !unfrozen     

c  Sutra parameters
c     tk = 1.3178e-6                                        !ice
c     if(tavg.gt.-4.d0 .and. tavg .lt.0.d0) tk = 4.6515e-7      !mushy
c     if(tavg.ge.0.d0) tk=1.2188e-6                           !unfrozen     

c Lunardini L=0 3 zones
c      tk = 5.0182e-6                                      !ice
c      if(tavg.ge.-4. .and. tavg .le.0.) tk = 4.26e-6      !mushy : lamda/C = 2.94135/690360
c      if(tavg.ge.0.) tk=3.5030e-6                         !unfrozen      

c Lunardini k and Co split; with latent heat included - use this for test
c      tk = 3.464                                 !ice
c      cotl0 = 690360.
c      if(tavg.ge.-4. .and. tavg .le.0.) then     !mushy
c         tk = 2.941                 
c         cotl0 = 1.78e7 
c      endif
c      if(tavg.ge.0.) tk=2.418                    !unfrozen
c
c     final values used in equation for Lunardini 3-zone hardwire:
c     if setting wwavg=1 and cpfavg=1, must make sure that all dispersion (and molecular diffusion) terms =0
c     ------------------------------------------------------------------------------------------------------
c      cotl1=cotl0
c      cotlavg = (cotl0+cotl1)/2.d0
c      wwavg = 1.d0
c      cpfavg = 1.d0
c      tccx = tk                 
c      tccy = tccx
c      tccz = tccx
c
c     debug print
c     ------------
c     if(l.eq.1 .or. l.eq.5 .or. l.eq.10.or.l.eq.ne)           
c    +write(6,756) l,tavg,por(l),sw(l),ww1,wi1,wm,wu1,dwu1,dlh1,cotl1,tk
c756  format(5x,'l,tavg,por,sw,ww1,wi1,wm,wu1,dwu1,dlh1,Co,tk:',
c    +            i7,11e11.3)

c     if(l.le.3.or. (ne-l).lt.3) then
c     if(l.eq.1) then
c     write(6,9333)
c    + l,tav1,cpe,cpm,dlh,ww,wi,por(l),sw(l),wu(tav1,p,q,modelwu,ts),
c    +                             dwu(tav1,p,q,modelwu,ts),col1,tk,tccx
c9333 format(2x,'l='i6,',tav1=',e10.3,',cpe=',e10.3,',cpm=',e10.3,
c    + ', dlh=',e10.3,', ww=',e10.3,', wi=',e10.3,
c    + ', por=',e10.3,', sw=',e10.3,', wu=',e10.3,', dwu=',e10.3,
c    + ', Co1=',e10.3,', tk=',e10.3', tccx=',e10.3)
c     endif

c     end hardwire
c --------------------------------------------------------------

c     endif                     !for heat transport

c     dispersion tensor
c     (retardation term is entered during elemental matrix assembly)
c     Dij * wwavg*cpfavg (not divided out now)
c     --------------------------------------------------------------
      vxl=vx(l)
      vyl=vy(l)
      vzl=vz(l)
      vx2=vxl*vxl
      vy2=vyl*vyl
      vz2=vzl*vzl
      vxy=vxl*vyl
      vxz=vxl*vzl
      vyz=vyl*vzl
      v= sqrt(vx2+vy2+vz2)

c
c     BF or Lichtner dispersion formulation:
c     ---------------------------------------
      if (kdisp.eq.0) then 
      al=alh
c  original B-F version
      dxx=((al*vx2+ath*vy2+atv*vz2)/v + dd) + tccx 
      dyy=((ath*vx2+al*vy2+atv*vz2)/v + dd) + tccy
      dzz=((atv*vx2+atv*vy2+al*vz2)/v + dd) + tccz
      dxy=(al-ath)*vxl*vyl/v
      dxz=(al-atv)*vxl*vzl/v
      dyz=(al-atv)*vyl*vzl/v
      
      else
c
c  Lichtner version (Lichtner et al, 2002)
c  added March 12, 2013, EOF
c
      costh=vzl/v
      cos2th=costh*costh
      vxy2=vx2+vy2
      all=alh+cos2th*(alv-alh)
      att=atv+cos2th*(ath-atv)
      dxx=((all*vx2+ath*vy2*(1.0+vz2/vxy2)+att*vz2*vx2/vxy2)/v+dd)+tccx
      dyy=((ath*vx2*(1.0+vz2/vxy2)+all*vy2+att*vz2*vy2/vxy2)/v+dd)+tccy
      dzz=((att*vxy2+all*vz2)/v+dd)+tccz
      dxy=(all-ath*(1.0+vz2/vxy2)+att*vz2/vxy2)*vxl*vyl/v
      dxz=(all-att)*vxl*vzl/v
      dyz=(all-att)*vyl*vzl/v

      endif

c      wucpf = wwavg*cpfavg
c      dxx=((al*vx2+ath*vy2+atv*vz2)/v + dd)  + tccx
c      dyy=((ath*vx2+al*vy2+atv*vz2)/v + dd)  + tccy
c      dzz=((atv*vx2+atv*vy2+al*vz2)/v + dd)  + tccz
c      dxy=(al-ath)*vxl*vyl/v  
c      dxz=(al-atv)*vxl*vzl/v 
c      dyz=(al-atv)*vyl*vzl/v  

c     hardwire
c      if(dxy.lt.1.0e-5) write(6,8837) dxy
c      if(dxz.lt.1.0e-5) write(6,8838) dxy
c      if(dyz.lt.1.0e-5) write(6,8839) dxy
c 8837 format(10x,'dxy: ',e15.4)
c 8838 format(10x,'dxz: ',e15.4)
c 8839 format(10x,'dyz: ',e15.4)
c      call flush(6)
c
c     elemental matrix, scheme 2
c     --------------------------
      rex=wp*dxx+wa*vx2*dt2            !wucpf ok - tested with craflush and oneD and lunardini July10
      rey=wp*dyy+wa*vy2*dt2   
      rez=wp*dzz+wa*vz2*dt2   
      rxy=wp*dxy+wa*vxy*dt2  
      rxz=wp*dxz+wa*vxz*dt2   
      ryz=wp*dyz+wa*vyz*dt2   
      pex=wp1*dxx+wa1*vx2*dt2   
      pey=wp1*dyy+wa1*vy2*dt2  
      pez=wp1*dzz+wa1*vz2*dt2   
      pxy=wp1*dxy+wa1*vxy*dt2   
      pxz=wp1*dxz+wa1*vxz*dt2  
      pyz=wp1*dyz+wa1*vyz*dt2  
      pvx=vxl                           
      pvy=vyl                          
      pvz=vzl                           
  155 continue
c
c
c     numerical(kint=1) or direct(kint=0) integration heat transport2
c     returns {re},{ro}
c     now calls exint2 
c     ===========================================================================
c     if (l.eq.1) write (6,712) rex,rey,rez,rxy,rxz,ryz,
c    1pex,pey,pez,pxy,pxz,pyz,pvx,pvy,pvz
c 712 format (9f14.7)
c
      if(kint.eq.1)
     + call gquad2(rex,rey,rez,rxy,rxz,ryz,pex,pey,pez,
     +      pxy,pxz,pyz,pvx,pvy,pvz,dt1,x,y,z,inl,re,ro,l,
     +      maxnn,cotl0,cotl1,cotlavg,ag,hag,tk)                 !cotlavg=R
c
       if(kint.eq.0)
     + call exint2(rex,rey,rez,rxy,rxz,ryz,pex,pey,pez,
     +   pxy,pxz,pyz,pvx,pvy,pvz,hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,
     +   dt1,x,y,z,inl,re,ro,maxnn,htmt,htmf,cotl0,cotl1,cotlavg,
     +   agefx,rg,exl,eyl,ezl,maxne,l,tk,vx,vy,vz)
c
c      if((l.eq.8218.or.l.eq.10330).and.(tdays.eq.6..or.tdays.eq.20.))
c     + then
c      write(6,281) tdays
c  281 format(/10x,'transport matrix check at tdays = ',f10.3)
c      write(6,181) l,(i,(re(i,j),j=1,8),i=1,8)
c      write(6,181) l,(i,(ro(i,j),j=1,8),i=1,8)
c 181  format(/,' transport matrix check - re,ro: ',i6,/(i2,8(e11.3,1x)))
c      endif
c
c     assembly for conjugate gradient solver
c     ======================================
      do 81 i=1,8
      ki = inl(i)
      if (ic(ki).eq.1) go to 81
      ii = ki-lc(ki)
      do 82 j=1,8
      kj = inl(j)
      if (ic(kj).eq.1) go to 85
      jj=kj-lc(kj)
      if (ii.gt.jj) go to 86
      do 83 k=1,nw
      if (ii.ne.ib(jj,k)) go to 83
      kk=k
      go to 84
   83 continue
   84 continue
      iv=ind(jj)+kk-1
c
c     write (6,738) l,i,j,ii,jj,ib(jj,1),ib(jj,kk),ind(jj),kk,iv,
c    1a(iv),se(i,j)
c 738 format (10i5,2f14.8)
c     -------------------------
c     left side coefficient vector
c     ----------------------------
      a(iv)=a(iv)+re(i,j)
c
c     right side contribution
c     -----------------------
   86 continue
      fs(ii)=fs(ii) - ro(i,j)*u0(kj) + rg(i,j)
      go to 82
c
c     linking term:
c     -------------
   85 continue
      f(ii)=f(ii)-(re(i,j)+ro(i,j))*u0(kj)
   82 continue
   81 continue
c
c     end of element loop
c     ===================
  100 continue
c
c     1D line and 2D plane element assembly - transport
c     Note: fracture thermal conductivity assumed =  0.5 J/msC (water)
c     hardwire porfrac
c     ----------------------------------------------------------------
c     porfrac=1.
c     vxmax=-999.
c     vymax=-999.
c     vxmin=+999.
c     vymin=+999.
c     do kf=1,nfrac
c     if(ifdim(kf).eq.2) then
c     vxmax=max(vx2d(kf),vxmax)
c     vymax=max(vy2d(kf),vymax)
c     vxmin=min(vx2d(kf),vxmin)
c     vymin=min(vy2d(kf),vymin)
c     endif
c     enddo
c     write(6,7822) vxmax,vymax,vxmin,vymin
c7822 format(/10x,'fractures: vxmax,vymax,vxmin,vymin: ',4e15.4)

c     all fracture elements:
c     -----------------------
      do 555 kf=1,nfrac
c       
c     1D fractures - transport
c     ------------------------
      if(ifdim(kf).eq.1) then
      vl = vlin(kf)            
      l =  ifracl(kf)
c     cotl0 = cot0(l)
c     cotl1 = cot1(l)
c
      call frac_line(maxfrac,maxne,lvert,inline,inflow,kf,l,fracdim,
     +              exl,eyl,ezl)
c
c     update tcc for frozen state -
c     --------------------------------------
      tb = 0.25d0
     +       * (u1(map(inline(1)))+u1(map(inline(2)))
     +         +u0(map(inline(1)))+u0(map(inline(2))))
c      tcc = tkl(l)/(cf*den(tb,lmass,lage))
c      tcc = 0.6d0/(cf*den(tb,lmass,lage))              
      wufrac = wu(tb,p2,q2,modelwu2,ts2)
      cpfavg = (cf*den(tb,lheat,lmass,lage))
      dwufrac = dwu(tb,p2,q2,modelwu2,ts2)
c     dlhfrac = por(1)*pi*wlh*dwufrac            !pi (ice density)       
      dlhfrac =        pi*wlh*dwufrac            !pi (ice density)       jan 2016 por(l) removed, por(1) = 1 !
      cotl = wufrac*cpfavg + (1.d0-wufrac)*cpi +dlhfrac            !Co+L
c      tcc = (wufrac*sqrt(0.6d0) + (1.d0-wufrac)*sqrt(tcli))**2     !tk
      tcc = (wufrac*sqrt(tclw) + (1.d0-wufrac)*sqrt(tcli))**2     !tk
c     if(lmass.or.lage) then            !use else-if to skip above for mass transport
c         tcc=0.d0
c         cotl = rtrans
c         wufrac=1.d0
c         cpfavg=1.d0
c     endif

c     use highest dispersivities for fractures
c     al = dmax1(alv,alh)
c     always use al=alh

      dxx = (alh*vl + dd)  + tcc           !tcc on lhs - check normally tcc=0 for fracture
      sqa = xarea(kf)/por(1)              !check por
      rezl = sqa  * ( wp*dxx  +  wa*vl*vl*dt2 ) / fracdim  
      pezl = sqa  * ( wp1*dxx + wa1*vl*vl*dt2 ) / fracdim 
       pvzl = sqa  *  (vl/2.d0)  
c      pvzl = sqa  *  (vl/2.d0)  * (wufrac*cpfavg) 
c
c     lumped -     use htl = xarea * fracture_length*dt1/2.
c     consistent - use htl = xarea * fracture_length*dt1/6.
c     check decay rdec (March2007): /fracdim or *fracdim
c  -----------------------------------------------------------------------------
c      htl =  cotl* (sqa  *  fracdim * dt1/2.d0)   !cotl (Co) on rhs, do not divide by R from blocks
      htl =   (sqa  *  fracdim * dt1/2.d0)/rett    !cotl (Co) on rhs,        divide by R from blocks (R assumed uniform)
c     rdec =  sqa  *  fracdim * decay(1) 
      rdec =  0.0d0                            !no decay in heat
c
      rel(1,1) = rezl/rett + htl + rdec
      rel(1,2) = -rezl/rett + htl + rdec
      rel(2,1) = rel(1,2)
      rel(2,2) = rel(1,1)
      rol(1,1) = (+pezl - pvzl)/rett - htl
      rol(1,2) = (-pezl + pvzl)/rett 
      rol(2,1) = (-pezl - pvzl)/rett
      rol(2,2) = (+pezl + pvzl)/rett - htl
c
      do 2012 i=1,2
      ki = map(inline(i))
      if(ic(ki).eq.1) goto 2012
      ii= ki-lc(ki)
        do 2013 j=1,2
        kj = map(inline(j))
        if(ic(kj).eq.1) goto 2014
        jj = kj - lc(kj)
        if(ii.gt.jj) goto 2015  
c
      do 2083 k=1,nw
      if (ii.ne.ib(jj,k)) go to 2083
      kk=k
      go to 2084
 2083 continue
 2084 continue
c
      iv=ind(jj)+kk-1
      a(iv)= a(iv) + rel(i,j) 
 2015 continue
      fs(ii) = fs(ii) - rol(i,j)*u0(kj) 
      goto 2013
 2014 continue
      f(ii) = f(ii) - (rel(i,j)+rol(i,j))*u0(kj)
c
 2013 continue
 2012 continue
 2011 continue
      endif
c       
c     2D fractures - transport
c     ------------------------
c     ----------------------------------------------------------------
      if(ifdim(kf).eq.2) then
c
c     fracture velocities: 
c     use 2 velocities: vxl,vyl where x,y could also be x,z or y,z
c     --------------------------------------------------------------

c hardwire v
c     vx2d(kf)=1.16e-3
c     vy2d(kf)=1.0e-15


      vx2 = vx2d(kf)**2
      vy2 = vy2d(kf)**2
      vxy = vx2d(kf)*vy2d(kf)
      vv = sqrt(vx2+vy2)
      l =  ifracl(kf)
c
c      if(kf.eq.1.or.kf.eq.nfrac)
c    + write(6,9922) l,vx2d(kf),vy2d(kf),vv,xarea(kf),dd
c9922 format(/10x,'2D fracture transport l= ',i5,
c    +       /10x,'vx2 = ',e12.4,' vy2 = ',e12.4,' vv= ',e12.4,
c    +       /10x,'xarea ... ',e12.5,' dd ... ',e12.5)
c
      call frac_plane(maxfrac,maxne,lvert,inline,inflow,kf,l,
     +                fdimx,fdimy,exl,eyl,ezl,rexel,sarea)
c
      mapinl1 = map(inline(1))
      mapinl2 = map(inline(2))
      mapinl3 = map(inline(3))
      mapinl4 = map(inline(4))
      tb = 0.125 * ( u1(mapinl1)+u1(mapinl2)+u1(mapinl3)+u1(mapinl4) 
     +              +u0(mapinl1)+u0(mapinl2)+u0(mapinl3)+u0(mapinl4) ) 

c     hardwire thermal diffusivity of fracture
c     divide by por (from porous matrix equation since we use v not q)
c     check: must assume uniform porosity, correct for water content
c     ----------------------------------------------------------------
c     tcc = 0.6d0/(cf*den(tb,lmass,lage))
      wufrac = wu(tb,p2,q2,modelwu2,ts2)      !por ?
      cpfavg = cf*den(tb,lheat,lmass,lage)     
      dwufrac = dwu(tb,p2,q2,modelwu2,ts2)
c     dlhfrac = por(1)*pi*wlh*dwufrac            !pi (ice density) check por 
      dlhfrac =        pi*wlh*dwufrac            !pi (ice density) check por  jan 2016 por(l) removed, por=1.
      cotl = wufrac*cpfavg + (1.d0-wufrac)*cpi + dlhfrac        !Co+latent heat, no solids in fracture   hardwire - por

c     tcc = (wufrac*sqrt(0.6d0) + (1.d0-wufrac)*sqrt(tcli))**2     !tk
      tcc = (wufrac*sqrt(tclw) + (1.d0-wufrac)*sqrt(tcli))**2     !tk
c     cotl =cotl/por(1)

c     if(lmass.or.lage) then            !use else-if to skip above for mass transport
c         tcc=0.d0
c         cotl = rtrans
c         wufrac=1.d0
c         cpfavg=1.d0
c     endif      

c hardwire for sfrac-h no conduction in fracture
      tcc=0.
c
      sqa = xarea(kf)/por(1)       !por ok
c      if(kf.eq.1.or.kf.eq.nfrac) 
c     +   write(6,*) 'kf,wufrac,cpfavg,por(1): ',kf,wufrac,cpfavg,por(1)
      wucpf = wufrac*cpfavg

c     use highest dispersivities for fractures
c      al = dmax1(alv,alh)
c      alphat = dmax1(atv,ath)
      alphat=ath
             
      dxx = ((alh*vx2 + alphat*vy2)/vv + dd)  + tcc           !check later: should not have tcc*wucp
      dyy = ((alphat*vx2 + alh*vy2)/vv + dd)  + tcc
      dxy = ((alh-alphat)*vx2d(kf)*vy2d(kf)/vv)  
c
c     if(kf.eq.1.or.kf.eq.nfrac)
c    +  write(6,7722) kf,xarea(kf),rexel,fdimx,fdimy,dt1,dt2
c7722 format(/10x,'fracture plane # ... ',i7,
c    +       /10x,'xarea: ',e12.5,' rexel: ',e12.5,
c    +       /10x,'fdimx: ',e12.5,' fdimy: ',e12.5,
c    +       /10x,'dt1,dt2: ',2e15.5)
      rex = sqa * ( wp*dxx  +  wa*vx2*dt2 ) * (rexel/6.d0)   
      rey = sqa * ( wp*dyy  +  wa*vy2*dt2 ) / (rexel*6.d0)   
      rxy = sqa * ((wp*dxy  +  wa*vxy*dt2 ) / 2.d0)         
      pex = sqa * ( wp1*dxx + wa1*vx2*dt2 ) * (rexel/6.d0)   
      pey = sqa * (( wp1*dyy + wa1*vy2*dt2 ) / (rexel*6.d0)) 
      pxy = sqa * (( wp1*dxy + wa1*vxy*dt2) / 2.d0)         
      pvx = sqa * (vx2d(kf)*fdimy/12.d0)                    
      pvy = sqa * (vy2d(kf)*fdimx/12.d0)                     
c
c     lumped - 
c     consistent - 
c     --------------------------------------------
c      ret0 = cotl*(sqa * sarea * dt1/36.d0)            !/por for cotl already in sqa
       ret0 =      (sqa * sarea * dt1/36.d0)/rett            !/por for cotl already in sqa
      ret1 = ret0
      rdec = sqa * sarea * decay(1) /36.d0
      rage = sqa * agefx * sarea / 36.d0
c
c     write(6,8888) l,dxx,dyy,dxy,rex,rey,rxy,ret1
c8888 format(/10x,'2D fracture - transport, l= ',i5,
c    +       /10x,' dxx,dyy,dxy: ',3e12.3,
c    +       /10x,' rex,rey,rxy: ',3e12.3,' ret1: ',e12.3)
c
      do 56 i=1,4
      do 56 j=1,4
      rd = rdec*ht(i,j)
      rel(i,j) = (rex*hx(i,j) + rey*hy(i,j) + rxy*hh(i,j))/rett
     +                                      + ret1*ht(i,j) + rd
      rol(i,j) = (pex*hx(i,j) + pey*hy(i,j) + pxy*hh(i,j)
     +           + pvx*wx(i,j) + pvy*wy(i,j))/rett  - ret0*ht(i,j)
      rg2df(i,j) = rage * ht(i,j)
  56  continue
c
      do 3012 i=1,4
      ki = map(inline(i))
      if(ic(ki).eq.1) goto 3012
      ii= ki-lc(ki)
        do 3013 j=1,4
        kj = map(inline(j))
        if(ic(kj).eq.1) goto 3014
        jj = kj - lc(kj)
        if(ii.gt.jj) goto 3015  
c
      do 3083 k=1,nw
      if (ii.ne.ib(jj,k)) go to 3083
      kk=k
      go to 3084
 3083 continue
 3084 continue
c
      iv=ind(jj)+kk-1
      a(iv)= a(iv) + rel(i,j) 
 3015 continue
      fs(ii) = fs(ii) - rol(i,j)*u0(kj) + rg2df(i,j)
      goto 3013
 3014 continue
      f(ii) = f(ii) - (rel(i,j)+rol(i,j))*u0(kj)
c
 3013 continue
 3012 continue
 3011 continue
c
      endif
c
c     get next 2D transport fracture
c     ------------------------------
 555  continue
c
c     contribution from thermal leakage boundary:
c     new: use fc2 (fc here) to store variable surface temp vtemp
c     check with new formulation Nov 2013
c     -------------------------------------------------------
      if(leak2) then
      k=0

      do 102 i=nz,nn,nz
      tcsbz=tcs/bz(i/nz)

c     skip if 1st type temp node
      if(ic(i).eq.1) goto 102

      i2=i/nz
      i3 = mod(i2,ny)
      if(i3.eq.0) i3=ny-1
      i4 = i2/ny +1
      if(mod(i2,ny).eq.0) i4=i4-1
      if(i4.eq.nx) i4=i4-1
      ielem = (nez-1)*nexy+(i3-1)*nex+i4     !need element # for porosity

      if(ielem.gt.ne.or.ielem.lt.1) then
      write(6,3984)
 3984 format(/10x,'error in assigning element# in trans ...',
     +       /10x,'stopping ... ')
      stop
      endif

      tavg=(u1(i)+u0(i))/2.
      wuleak = wu(tavg,ppn(i),qqn(i),modelwu,ts)
      ww=porsurf*sat*wuleak
      wi=porsurf*sat*(1.0-wuleak)
      wm=(1.0-porsurf) 
      dwuleak = dwu(tavg,ppn(i),qqn(i),modelwu,ts)
      dlhleak=por(ielem)*pi*wlh*dwuleak            !pi (ice density), don't use porsurf ?

      cpf=cf*den(tavg,lheat,lmass,lage)
c     cpe=ww*cpf + wi*cpi + wm*cpm
      cpeu = ww*cpf + wi*cpi + (1.-porsurf)*cpsm(ielem) + dlhleak 
c     p1cpm=(1.-porsurf)*cpsm(ielem)

c     rtemp=fc(i)
      rtemp = max(0.,ttemp+rtdiff)          !recharge water temperature rtemp, min=0
      cfd=cf*den(rtemp,lheat,lmass,lage)          !
      k=k+1
      ii=i-lc(i)
      
      bdyflux = fbf(k)
      if(icflow(i).eq.1) bdyflux = -vz(ielem)*por(ielem)     !for type1 flow node, use v*por as flux
      if(bdyflux.lt.0.) bdyflux=0.

c     write(6,5676) k,ielem,fbf(k),ara(k)
c5676 format(i5,'leak2: processing element # ',i7,',fbf,ara= ',2e10.3)
c     cpeu=por(ielem)*sat*cf*den(tavg,lmass) + p1cpm
c     cpe=ww*cpf + wi*cpi + wm*cpm
c
c     check with new formulation Co on rhs, ww*cpf on lhs?
c     check ... lhs units are m^2*(J/m/s/C)/m = J/s
c     --------------------------------------------------------------
       arak=ara(k)
       a(ind(ii)) = a(ind(ii)) + (arak*(tcsbz + bdyflux*cfd))      !new for recharge, but then removed from lhs, ok.
       f(ii)=f(ii) + arak*(fc(i)*tcsbz + rtemp*bdyflux*cfd + bzf(i/nz))  ! added ara(k)*bzf from bzin.data

c      write(6,4949) k,ara(k),ii,f(ii),i,nz,bzf(i/nz)
c4949  format(' check: k,ara(k),ii,f(ii),i,nz,bzf(i/nz):',
c    +                 i5,e15.5,i8,e15.5,2i8,e15.5)


c  ok to march 2015, arak not factored, too slow ...
c      a(ind(ii)) = a(ind(ii)) + (arak*tcsbz) 
c     +                        + (arak*bdyflux*cfd)      !new for recharge, but then removed from lhs, ok.
c      f(ii)=f(ii) + (arak*fc(i))*tcsbz + (arak*rtemp*bdyflux*cfd)
c     +            + (arak*bzf(i/nz))
c original:
c      a(ind(ii)) = a(ind(ii)) + ara(k)*tcs/(bz*cpeu)
c      f(ii)=f(ii)+ (ara(k)*rtemp/cpeu)*(tcsbz+fbf(k)*cfd)      
  102 continue
      endif

c      write(6,4747) (k,ara(k),bzf(k),k=1,nn/nz)
c4747  format(i5,2e15.5)      
c
c      write (6,695) (i,fs(i),fc(i),i=nz,nn,nz)        !fs undefined 2017
c 695  format(/10x,'i,fs(i),fc(i)...',/(2(i6,2e12.3)))
c
c     source term contribution
c     need to multiply pqq by wwavg*cpf as in 3D blocks
c     (pqq*tq) units are (m^3/s)*(J/m^3.C)*C = J/s
c     --------------------------------------------------
      do 556 i=1,nn
      ii=i-lc(i)
      if(pq(i).le.0.) goto 557
      tavg=(u1(i)+u0(i))/2.      
      wwavg = por(1)*sw(1)*wu(tavg,ppn(i),qqn(i),modelwu,ts)     !check - need element porosity, Sw
c     write(6,34857) 
c    +    por(1),sw(1),wu(tavg,p,q,modelwu,ts),cf,den(tavg,lmass,lage)
34857 format('por(1),sw(1),wu,cf,den: ',5e15.5)

      cpfavg = (cf*den(tavg,lheat,lmass,lage))            
c     cpfavg = cf            
c     if(lmass.or.lage) then
c              wwavg = 1.d0
c              cpfavg=1.d0
c     endif
      pqq= pq(i)*wwavg*cpfavg/(por(1))       !check - need wwavg*cpfavg  - 
c     pqq= pq(i)*wwavg*cpfavg       !check - need wwavg*cpfavg  - 
      a(ind(ii)) = a(ind(ii)) + pqq
      f(ii)=f(ii)+pqq*tq(i)
  557 continue
  556 continue
c
c     internal heat flux elements:
c     Feb. 2005: assumes unit element column, or 4 columns 
c     hardwire for flux across outside faces only. 
c     - distribute over 4 nodes on each of the two outside faces, = (total area/8 nodes)
c     use elemental cpe calculated from stored cot0(l) and cot1(l) arrays
c     needs to be updated for freezing
c     -----------------------------------------------------------------------------------
      if(numis.gt.0) then
      do 667 i=1,numis
      l = ifl(i)
c     cpe = (cot0(l) + cot0(l))/2.
      
c     wufrac = wu(tb,p,q,modelwu,ts)
c     cpfavg = (cf*den(tb,lmass,lage))
c     cotl = wufrac*cpfavg + (1.d0-wufrac)*cpi                     !Co

c     original
c     sarea = 0.25*(exl(l)*ezl(l) + exl(l)*eyl(l) + eyl(l)*ezl(l))
c     sarea = (2./8.)*(exl(l)*ezl(l) + eyl(l)*ezl(l))
      sarea = (exl(l)*ezl(l) + eyl(l)*ezl(l))
      if(nflux(i,3).eq.in(l,3)) sarea = sarea*2.    !single column: area*2

      sarea = sarea/8.d0                                !distribute over 8 nodes
c     cpe=por(l)*cpf + (1.d0-por(l))*cpm    !need to adjust for freezing
c     if(l.gt.nef) cpe=por(l)*sat*cpf + (1.d0-por(l))*cpm
      do 668 j=1,8
      node = nflux(i,j)
      ii = node - lc(node)
c      f(ii) = f(ii) + (flux(i)*sarea/cpe) * exp(-decay2*tdays)  !hardwire decay off
c     f(ii) = f(ii) + (flux(i)*sarea/cpe)
      f(ii) = f(ii) + flux(i)*sarea

c     write(6,8475) l,ii,flux(i),sarea,cpe
c8475 format(10x,'debug: l,ii,flux,sarea,cpe: ',2i10,3e15.5) 

 668  continue
 667  continue

      endif
c
c     assemble flux vector
c     --------------------
      do 121 i=1,n
      f(i)= f(i) + fs(i) 
  121 continue
c
c     if (it.gt.1) go to 735
c     do 733 i=1,maxna
c 733 write (6,734) i,iaa(i),a(i)
c 734 format (2i5,f14.8)
c 735 continue
c     tsurfmin=999.
c     tsurfmax=-999.
c     fcmin=999.
c     fcmax=-999.
c     u0min=999.
c     u0max=-999.
c     do i=1,nn 
c     tsurfmin = min(tsurfmin,u2(i))
c     tsurfmax = max(tsurfmax,u2(i))           
c     fcmin = min(fcmin,fc(i))
c     fcmax = max(fcmax,fc(i))           
c     u0min = min(u0min,u0(i))
c     u0max = max(u0max,u0(i))           
c     enddo
c
c     write(6,3476) tsurfmin,tsurfmax,fcmin,fcmax,u0min,u0max
c3476 format(/10x,'before pcg in trans2 min,max temps,fc,u0: ',6e15.5)

c
c     conjugate gradient solution - transport2
c     =========================================
      neu=1
c     call clock@(start)
      call precg (a,n,iaa,ind,f,u2,aa,laa,neu)
c      write(*,123)
c  123 format(10x,'transport pcg solution complete ...')
c     call clock@(finish)
c     tpcgt=tpcgt+finish-start
c
   90 continue
c     tsurfmin=999.
c     tsurfmax=-999.
c     fcmin=999.
c     fcmax=-999.
c     u0min=999.
c     u0max=-999.
c     do i=1,nn 
c     tsurfmin = min(tsurfmin,u2(i))
c     tsurfmax = max(tsurfmax,u2(i))           
c     fcmin = min(fcmin,fc(i))
c     fcmax = max(fcmax,fc(i))           
c     u0min = min(u0min,u0(i))
c     u0max = max(u0max,u0(i))           
c     enddo
c
c     write(6,3556) tsurfmin,tsurfmax,fcmin,fcmax,u0min,u0max
c3556 format(/10x,'after pcg in trans2 min,max temps,fc,u0: ',6e15.5)

c
c     expand solution vector:
c     -----------------------
      do 141 ii=1,nn
      i = nn-ii+1
      if(ic(i).eq.1) goto 142
      k=i-lc(i)
c     if(lmass .and. u2(k).lt.1.d-50) u2(k)=0.d0  !hardwire check for transport
      u2(i)=u2(k)
      goto 141
  142 u2(i)=fc(i)
  141 continue
c
c     tsurfmin=999.
c     tsurfmax=-999.
c     fcmin=999.
c     fcmax=-999.
c     u0min=999.
c     u0max=-999.
c     do i=1,nn 
c     tsurfmin = min(tsurfmin,u2(i))
c     tsurfmax = max(tsurfmax,u2(i))           
c     fcmin = min(fcmin,fc(i))
c     fcmax = max(fcmax,fc(i))           
c     u0min = min(u0min,u0(i))
c     u0max = max(u0max,u0(i))           
c     enddo
c
c     write(6,3377) tsurfmin,tsurfmax,fcmin,fcmax,u0min,u0max
c3377 format(/10x,'after trans2 min,max temps,fc,u0: ',6e15.5)
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      return
      end
c
c
c ******************************************************************
c ******************************************************************
c *******************************************************************
c
c trans.for
c
c    ******** subroutine  transport 3  **********   Mass transport based on Original thermal R  version
c
c  three-dimensional advective-dispersive transport
c  linear isoparametric quadrilateral elements
c
c     e.o. frind, j.w. molson
c     university of waterloo
c     (c) 1993
c
c   *******************************************************************
c
c       copyright 1987 emil o. frind, 1988 j.w. molson
c
c       duplication of this program or any part thereof without
c       the express written consent of the author is prohibited.
c
c   *******************************************************************
c
      subroutine trans3(maxn,maxnn,maxne,nf,maxna,laa,nw,
     + x,y,z,u0,u1,u2,fc,fb,f,fs,icflow,ic,lc,vx,vy,vz,inflow,in,map,
     + a,aa,iaa,ind,ib,alh,alv,ath,atv,dd,kdisp,n,nbb,dt1,dt2,nz,
     + ne,nn,wp,wa,wp1,wa1,pq,tq,cq,por,sw,ag,hag,ttemp,ara,leak2,
     + hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,htmt,kint,tast,tpcgt,tdays,
     + exl,eyl,ezl,wlh,pp,qq,modelwu,vt,tcs,bz,bzf,cpsm,cpi,cf,pi,sat,
     + nef,
     + fbf,rinc,htmf,ifracl,lvert,nfrac,xarea,nex,ney,nez,vlin,maxfrac,
     + tclw,tcli,tclm,ifdim,hgz2,hx,hy,hh,ht,wx,wy,vx2d,vy2d,
     + cot0,cot1,tkl,numis,flux,fluxv,ifl,maxfx,rtrans,decay,decay2,
     + nflux,agefx,rg,rg2df,porsurf,lmass,lheat,lage,ts,rtdiff,
     + p2,q2,ts2,omega2,modelwu2,modelkr2,depc,depd,vexp,smax,spn0,spn1,
     + spe0,spe1,sdepd,ppn,qqn,rhob,t0,t1,icnt)
c
      implicit real*8(a-h,o-z)
      external wu,dwu
      real*8 start,finish
c
c     constant dimensions
c     direct element matrices:
c     -------------------------
      dimension hxx(8,8),hyy(8,8),hzz(8,8),hxy(8,8),hxz(8,8),hyz(8,8),
     + hvx(8,8),hvy(8,8),hvz(8,8),htmt(8,8),htmf(8,8),
     + rg(8,8),rg2df(4,4)
c
      dimension re(8,8),ro(8,8),inl(8),kb(7),ag(9),hag(9)
      logical leak2,lmass,lheat,lage
c
c     ne = number of elements
c     nn = number of nodes
c     n  = number of degrees of freedom
c     na = total number of non-zero matrix entries in condensed matrix
c        = 14*n (approximately)
c     nw = width of matrix ib ( = 15 for prismatic grid )
c     nf = maximum number of nodes on one face
c     laa = 3*n + na
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     variable dimensions:
c     arrays of size nn
c     -----------------
      real*8 tq(maxnn),pq(maxnn),cq(maxnn)
      real*8  x(maxnn),y(maxnn),z(maxnn),por(maxne),sw(maxne)
      real*8 u0(maxnn),u1(maxnn),fc(maxnn),fb(maxnn),fs(maxnn)
      real*8 t0(maxnn),t1(maxnn)
      dimension  f(maxnn),u2(maxnn),fbf(nf)
      integer   ic(maxnn),icflow(maxnn),lc(maxnn),map(maxnn)
      integer icnt(maxnn)
c
c     1D line or 2D plane elements - transport 
c     -----------------------------------------
      dimension rel(4,4),rol(4,4),ifracl(maxfrac),inline(4)
      dimension lvert(maxfrac),ifdim(maxfrac)
      real*8 vlin(maxfrac),vx2d(maxfrac),vy2d(maxfrac),xarea(maxfrac)
      dimension hx(4,4),hy(4,4),hh(4,4),ht(4,4),wx(4,4),wy(4,4),hgz2(4)
c
c     arrays of size ne, (ne,8)
c     ------------------------
      real*8   vx(maxne),vy(maxne),vz(maxne)
      real*8   exl(maxne),eyl(maxne),ezl(maxne)
      real*8   vt(maxne),tclm(maxne),cot0(maxne),cot1(maxne)
      real*8   tkl(maxne),cpsm(maxne)
      integer  inflow(maxne,8),in(maxne,8)
      real*8   depc(maxne),depd(maxne),vexp(maxne),decay(maxne)
      real*8   spn0(maxnn),spn1(maxnn),rhob(maxne)
      real*8   spe0(maxne),spe1(maxne),sdepd(maxne),smax(maxne)
      real*8   pp(maxne),qq(maxne),ppn(maxnn),qqn(maxnn)
c
c     arrays for internal heat source
c     --------------------------------
      dimension flux(maxfx),fluxv(maxfx),ifl(maxfx),nflux(maxfx,8)
c
c     array of size nf
c     ----------------
      dimension ara(nf)
      dimension bz(nf),bzf(nf)
c
c     arrays a(na),iaa(na),ind(n+1),ib(n,15),aa(3n+na)
c            for conjugate gradient solver only
c     ------------------------------------------------
      dimension a(maxna),aa(laa)
      integer iaa(maxna),ind(maxnn),ib(maxnn,nw)
c
c     function defining -absolute- water density as a function of t
c     for use in latent heat term, thermal parameters
c     --------------------------------------------------------------
c
      nx=nex+1
      ny=ney+1
      nexy = nex*ney

c     initialize arrays,
c     condense non-zero flux array:
c     -------------------------------
      do 45 i=1,nn
      fs(i)=0.d0
      f(i)=0.d0
      if(ic(i).eq.1) goto 45
      ii=i-lc(i)
      f(ii)=fb(i) 
  45  continue
c
      do 247 i=1,maxna
  247 a(i)=0.d0

c
c     decay of source Dirichlet concentrations
c     u0 is used as linking term, fc here is fc3 (=C0)
c     later: apply to Cauchy ?
c     --------------------------------------------------
      do i=1,nn
c     fc(i)=fc(i)*exp(-decay2*tdays*86400)
      if(ic(i).eq.1)  u0(i)=fc(i)*exp(-decay2*tdays*86400)
      enddo
c
c     _________________________________________
c
c     matrix generation loop
c     loop over elements:
c     =========================================
c     =========================================
c
c     call clock@(start)
c
c     write(6,903)
c903  format(56x,'l     tavg       por        sw           ww1'
c    +'      wi1          wm         wu1       dwu1      dlh1'
c    +'      Co       tk')

      do 100 l=1,ne
c
c     element average T (for viscosity in exint3)
c     and conc
c     --------------------------------------------
      tavg=0.d0
      cavg=0.d0
      spnavg=0.d0
      do     i=1,8
      inl(i)=in(l,i)
      tavg=tavg+ (t0(inl(i))+t1(inl(i)))
      cavg=cavg+ (u0(inl(i))+u1(inl(i)))
      spnavg=spnavg+ (spn0(inl(i))+spn1(inl(i)))
      enddo
      tavg=tavg/16.d0    !divide by 2 for centre-weighting,then by 8 nodes
      cavg=cavg/16.d0    !divide by 2 for centre-weighting,then by 8 nodes
      spnavg=spnavg/16.d0    !divide by 2 for centre-weighting,then by 8 nodes
c
c     average elemental temperatures
c     tav0: avg. at start of time step / tav1: avg. at latest iteration
c     -----------------------------------------------------------------
c     tav0=0.d0
c     tav1=0.d0
c     do 103 i=1,8
c     tav0=tav0+u0(inl(i))
c     tav1=tav1+u1(inl(i))
c103  continue
c     tav0=tav0/8.d0
c     tav1=tav1/8.d0
c     tavg=(tav0+tav1)/2.0
c     tavg=tav2
c
c     fixed retardation = rtrans for mass transport or age only 
c     which is always the case for this saubroutine 
c     ----------------------------------------------------------
c     if(lmass.or.lage) then

c     cot0(l)=rtrans        !for mass transport, cot0 and cot1 represent retardation, and go on rhs 
c     cot1(l)=rtrans
      cotl0=rtrans
      cotl1=rtrans
c      if(l.eq.1 .or.l.eq.ne) write(6,*) 'l,rtrans=', l,rtrans
      tccx=0.d0             !thermal k = 0 for mass transport
      tccy=0.d0
      tccz=0.d0

c     wwavg=1.d0            !check for other terms below for heat that need to be defined for transport
c     cpfavg=1.d0
c     cotlavg = (cotl0+cotl1)/2.d0
c     rett=cotlavg
      cotlavg = rtrans
      rett = rtrans

c     dispersion tensor
c     (retardation term is entered during elemental matrix assembly)
c     Dij * wwavg*cpfavg (not divided out now)
c     --------------------------------------------------------------
      vxl=vx(l)
      vyl=vy(l)
      vzl=vz(l)
      vx2=vxl*vxl
      vy2=vyl*vyl
      vz2=vzl*vzl
      vxy=vxl*vyl
      vxz=vxl*vzl
      vyz=vyl*vzl
      v= sqrt(vx2+vy2+vz2)
c
c     BF or Lichtner dispersion formulation:
c     ---------------------------------------
      if (kdisp.eq.0) then 
      al=alh
c     original B-F version
c     wucpf = wwavg*cpfavg
      dxx=((al*vx2+ath*vy2+atv*vz2)/v + dd)  + tccx
      dyy=((ath*vx2+al*vy2+atv*vz2)/v + dd)  + tccy
      dzz=((atv*vx2+atv*vy2+al*vz2)/v + dd)  + tccz
      dxy=(al-ath)*vxl*vyl/v  
      dxz=(al-atv)*vxl*vzl/v 
      dyz=(al-atv)*vyl*vzl/v  

      else
c
c  Lichtner version (Lichtner et al, 2002)
c  added Feb 2017
      costh=vzl/v
      cos2th=costh*costh
      vxy2=vx2+vy2
      all=alh+cos2th*(alv-alh)
      att=atv+cos2th*(ath-atv)
c
      dxx=(all*vx2+ath*vy2*(1.0+vz2/vxy2)+att*vz2*vx2/vxy2)/v+dd
      dyy=(ath*vx2*(1.0+vz2/vxy2)+all*vy2+att*vz2*vy2/vxy2)/v+dd
      dzz=(att*vxy2+all*vz2)/v+dd
      dxy=(all-ath*(1.0+vz2/vxy2)+att*vz2/vxy2)*vxl*vyl/v
      dxz=(all-att)*vxl*vzl/v
      dyz=(all-att)*vyl*vzl/v

      endif
c      
c     hardwire
c      if(dxy.lt.1.0e-5) write(6,8837) dxy
c      if(dxz.lt.1.0e-5) write(6,8838) dxy
c      if(dyz.lt.1.0e-5) write(6,8839) dxy
c 8837 format(10x,'dxy: ',e15.4)
c 8838 format(10x,'dxz: ',e15.4)
c 8839 format(10x,'dyz: ',e15.4)
c      call flush(6)
c
c     elemental matrix, scheme 2
c     --------------------------
      rex=wp*dxx+wa*vx2*dt2            !wucpf ok - tested with craflush and oneD and lunardini July10
      rey=wp*dyy+wa*vy2*dt2   
      rez=wp*dzz+wa*vz2*dt2   
      rxy=wp*dxy+wa*vxy*dt2  
      rxz=wp*dxz+wa*vxz*dt2   
      ryz=wp*dyz+wa*vyz*dt2   
      pex=wp1*dxx+wa1*vx2*dt2   
      pey=wp1*dyy+wa1*vy2*dt2  
      pez=wp1*dzz+wa1*vz2*dt2   
      pxy=wp1*dxy+wa1*vxy*dt2   
      pxz=wp1*dxz+wa1*vxz*dt2  
      pyz=wp1*dyz+wa1*vyz*dt2  
      pvx=vxl                           
      pvy=vyl                          
      pvz=vzl                           
  155 continue
c
c
c     numerical(kint=1) or direct(kint=0) integration: trans3
c     returns {re},{ro}
c     ========================================================
c     if (l.eq.1) write (6,712) rex,rey,rez,rxy,rxz,ryz,
c    1pex,pey,pez,pxy,pxz,pyz,pvx,pvy,pvz
c 712 format (9f14.7)
c
      if(kint.eq.1)
     + call gquad3(rex,rey,rez,rxy,rxz,ryz,pex,pey,pez,
     +      pxy,pxz,pyz,pvx,pvy,pvz,dt1,x,y,z,inl,re,ro,l,
     +      maxnn,cotl0,cotl1,cotlavg,ag,hag,tk)                 !cotlavg=R
c
       if(kint.eq.0)
     + call exint3(rex,rey,rez,rxy,rxz,ryz,pex,pey,pez,
     + pxy,pxz,pyz,pvx,pvy,pvz,hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,
     + dt1,x,y,z,inl,re,ro,maxnn,htmt,htmf,cotl0,cotl1,cotlavg,decay,
     + agefx,rg,exl,eyl,ezl,maxne,l,tk,depc,depd,vexp,smax,spnavg,
     + vx,vy,vz,rhob,tavg,cavg,por,sdepd,lheat,lmass,lage)
c
c----
c      if((l.eq.8218.or.l.eq.10330).and.(tdays.eq.6..or.tdays.eq.20.))
c     + then
c      write(6,281) tdays
c  281 format(/10x,'transport matrix check at tdays = ',f10.3)
c      write(6,181) l,(i,(re(i,j),j=1,8),i=1,8)
c      write(6,181) l,(i,(ro(i,j),j=1,8),i=1,8)
c 181  format(/,' transport matrix check - re,ro: ',i6,/(i2,8(e11.3,1x)))
c      endif
c
c     assembly for conjugate gradient solver
c     ======================================
      do 81 i=1,8
      ki = inl(i)
      if (ic(ki).eq.1) go to 81
      ii = ki-lc(ki)
      do 82 j=1,8
      kj = inl(j)
      if (ic(kj).eq.1) go to 85
      jj=kj-lc(kj)
      if (ii.gt.jj) go to 86
      do 83 k=1,nw
      if (ii.ne.ib(jj,k)) go to 83
      kk=k
      go to 84
   83 continue
   84 continue
      iv=ind(jj)+kk-1
c
c     write (6,738) l,i,j,ii,jj,ib(jj,1),ib(jj,kk),ind(jj),kk,iv,
c    1a(iv),se(i,j)
c 738 format (10i5,2f14.8)
c     -------------------------
c     left side coefficient vector
c     ----------------------------
      a(iv)=a(iv)+re(i,j)
c
c     right side contribution
c     -----------------------
   86 continue
      fs(ii)=fs(ii) - ro(i,j)*u0(kj) + rg(i,j)
      go to 82
c
c     linking term:
c     -------------
   85 continue
      f(ii)=f(ii)-(re(i,j)+ro(i,j))*u0(kj)
   82 continue
   81 continue
c
c     end of element loop
c     ===================
  100 continue
c
c     1D line and 2D plane element assembly - transport
c     Note: fracture thermal conductivity assumed =  0.5 J/msC (water)
c     hardwire porfrac
c     ----------------------------------------------------------------
c     porfrac=1.
c     vxmax=-999.
c     vymax=-999.
c     vxmin=+999.
c     vymin=+999.
c     do kf=1,nfrac
c     if(ifdim(kf).eq.2) then
c     vxmax=max(vx2d(kf),vxmax)
c     vymax=max(vy2d(kf),vymax)
c     vxmin=min(vx2d(kf),vxmin)
c     vymin=min(vy2d(kf),vymin)
c     endif
c     enddo
c     write(6,7822) vxmax,vymax,vxmin,vymin
c7822 format(/10x,'fractures: vxmax,vymax,vxmin,vymin: ',4e15.4)

c     all fracture elements:
c     -----------------------
      do 555 kf=1,nfrac
c       
c     1D fractures - transport
c     ------------------------
      if(ifdim(kf).eq.1) then
      vl = vlin(kf)            
      l =  ifracl(kf)
c     cotl0 = cot0(l)
c     cotl1 = cot1(l)
c
      call frac_line(maxfrac,maxne,lvert,inline,inflow,kf,l,fracdim,
     +              exl,eyl,ezl)
c
c     update tcc for frozen state -
c     --------------------------------------
      tb = 0.25d0
     +       * (u1(map(inline(1)))+u1(map(inline(2)))
     +         +u0(map(inline(1)))+u0(map(inline(2))))
c     if(lmass.or.lage) then            !use else-if to skip above for mass transport
          tcc=0.d0
          cotl = rtrans
          wufrac=1.d0
          cpfavg=1.d0
c     endif

c     al = dmax1(alv,alh)
      dxx = (alh*vl + dd)  + tcc           !tcc on lhs - check normally tcc=0 for fracture
      sqa = xarea(kf)/por(1)              !check por
      rezl = sqa  * ( wp*dxx  +  wa*vl*vl*dt2 ) / fracdim  
      pezl = sqa  * ( wp1*dxx + wa1*vl*vl*dt2 ) / fracdim 
       pvzl = sqa  *  (vl/2.d0)  
c      pvzl = sqa  *  (vl/2.d0)  * (wufrac*cpfavg) 
c
c     lumped -     use htl = xarea * fracture_length*dt1/2.
c     consistent - use htl = xarea * fracture_length*dt1/6.
c     check decay rdec (March2007): /fracdim or *fracdim
c  -----------------------------------------------------------------------------
c      htl =  cotl* (sqa  *  fracdim * dt1/2.d0)   !cotl (Co) on rhs, do not divide by R from blocks
      htl =   (sqa  *  fracdim * dt1/2.d0)/rett    !cotl (Co) on rhs,        divide by R from blocks (R assumed uniform)
      rdec =  sqa  *  fracdim * decay(1)           !hardwire decay for fractures (ok if linear & homogeneous)
c
      rel(1,1) = rezl/rett + htl + rdec
      rel(1,2) = -rezl/rett + htl + rdec
      rel(2,1) = rel(1,2)
      rel(2,2) = rel(1,1)
      rol(1,1) = (+pezl - pvzl)/rett - htl
      rol(1,2) = (-pezl + pvzl)/rett 
      rol(2,1) = (-pezl - pvzl)/rett
      rol(2,2) = (+pezl + pvzl)/rett - htl
c
      do 2012 i=1,2
      ki = map(inline(i))
      if(ic(ki).eq.1) goto 2012
      ii= ki-lc(ki)
        do 2013 j=1,2
        kj = map(inline(j))
        if(ic(kj).eq.1) goto 2014
        jj = kj - lc(kj)
        if(ii.gt.jj) goto 2015  
c
      do 2083 k=1,nw
      if (ii.ne.ib(jj,k)) go to 2083
      kk=k
      go to 2084
 2083 continue
 2084 continue
c
      iv=ind(jj)+kk-1
      a(iv)= a(iv) + rel(i,j) 
 2015 continue
      fs(ii) = fs(ii) - rol(i,j)*u0(kj) 
      goto 2013
 2014 continue
      f(ii) = f(ii) - (rel(i,j)+rol(i,j))*u0(kj)
c
 2013 continue
 2012 continue
 2011 continue
      endif
c       
c     2D fractures - transport
c     ------------------------
c     ----------------------------------------------------------------
      if(ifdim(kf).eq.2) then
c
c     fracture velocities: 
c     use 2 velocities: vxl,vyl where x,y could also be x,z or y,z
c     --------------------------------------------------------------

c hardwire v
c     vx2d(kf)=1.16e-3
c     vy2d(kf)=1.0e-15


      vx2 = vx2d(kf)**2
      vy2 = vy2d(kf)**2
      vxy = vx2d(kf)*vy2d(kf)
      vv = sqrt(vx2+vy2)
      l =  ifracl(kf)
c
c      if(kf.eq.1.or.kf.eq.nfrac)
c    + write(6,9922) l,vx2d(kf),vy2d(kf),vv,xarea(kf),dd
c9922 format(/10x,'2D fracture transport l= ',i5,
c    +       /10x,'vx2 = ',e12.4,' vy2 = ',e12.4,' vv= ',e12.4,
c    +       /10x,'xarea ... ',e12.5,' dd ... ',e12.5)
c
      call frac_plane(maxfrac,maxne,lvert,inline,inflow,kf,l,
     +                fdimx,fdimy,exl,eyl,ezl,rexel,sarea)
c
      mapinl1 = map(inline(1))
      mapinl2 = map(inline(2))
      mapinl3 = map(inline(3))
      mapinl4 = map(inline(4))
      tb = 0.125 * ( u1(mapinl1)+u1(mapinl2)+u1(mapinl3)+u1(mapinl4) 
     +              +u0(mapinl1)+u0(mapinl2)+u0(mapinl3)+u0(mapinl4) ) 

c     hardwire thermal diffusivity of fracture
c     divide by por (from porous matrix equation since we use v not q)
c     check: must assume uniform porosity, correct for water content
c     ----------------------------------------------------------------

c     if(lmass.or.lage) then            !use else-if to skip above for mass transport
          tcc=0.d0
          cotl = rtrans
          wufrac=1.d0
          cpfavg=1.d0
c     endif      

c hardwire for sfrac-h no conduction in fracture
c     tcc=0.
c
      sqa = xarea(kf)/por(1)       !por ok
c      if(kf.eq.1.or.kf.eq.nfrac) 
c     +   write(6,*) 'kf,wufrac,cpfavg,por(1): ',kf,wufrac,cpfavg,por(1)
      wucpf = wufrac*cpfavg

c     use highest dispersivities for fractures
c     al = dmax1(alv,alh)
c     alphat = dmax1(atv,ath)
c     alphat = ath
c     if(lvert(kf).eq.3 .or. lvert(kf).eq.4) alphat = atv
      alphat=ath

      dxx = ((alh*vx2 + alphat*vy2)/vv + dd)  + tcc           !check later: should not have tcc*wucp
      dyy = ((alphat*vx2 + alh*vy2)/vv + dd)  + tcc
      dxy = ((alh-alphat)*vx2d(kf)*vy2d(kf)/vv)  
c
c     if(kf.eq.1.or.kf.eq.nfrac)
c    +  write(6,7722) kf,xarea(kf),rexel,fdimx,fdimy,dt1,dt2
c7722 format(/10x,'fracture plane # ... ',i7,
c    +       /10x,'xarea: ',e12.5,' rexel: ',e12.5,
c    +       /10x,'fdimx: ',e12.5,' fdimy: ',e12.5,
c    +       /10x,'dt1,dt2: ',2e15.5)
      rex = sqa * ( wp*dxx  +  wa*vx2*dt2 ) * (rexel/6.d0)   
      rey = sqa * ( wp*dyy  +  wa*vy2*dt2 ) / (rexel*6.d0)   
      rxy = sqa * ((wp*dxy  +  wa*vxy*dt2 ) / 2.d0)         
      pex = sqa * ( wp1*dxx + wa1*vx2*dt2 ) * (rexel/6.d0)   
      pey = sqa * (( wp1*dyy + wa1*vy2*dt2 ) / (rexel*6.d0)) 
      pxy = sqa * (( wp1*dxy + wa1*vxy*dt2) / 2.d0)         
      pvx = sqa * (vx2d(kf)*fdimy/12.d0)                    
      pvy = sqa * (vy2d(kf)*fdimx/12.d0)                     
c
c     lumped - 
c     consistent - 
c     --------------------------------------------
c      ret0 = cotl*(sqa * sarea * dt1/36.d0)            !/por for cotl already in sqa
       ret0 =      (sqa * sarea * dt1/36.d0)/rett            !/por for cotl already in sqa
      ret1 = ret0
      rdec = sqa * sarea * decay(1) /36.d0          !hardwire decay for fractures (ok if linear & homogeneous)
      rage = sqa * agefx * sarea / 36.d0
c
c     write(6,8888) l,dxx,dyy,dxy,rex,rey,rxy,ret1
c8888 format(/10x,'2D fracture - transport, l= ',i5,
c    +       /10x,' dxx,dyy,dxy: ',3e12.3,
c    +       /10x,' rex,rey,rxy: ',3e12.3,' ret1: ',e12.3)
c
      do 56 i=1,4
      do 56 j=1,4
      rd = rdec*ht(i,j)
      rel(i,j) = (rex*hx(i,j) + rey*hy(i,j) + rxy*hh(i,j))/rett
     +                                      + ret1*ht(i,j) + rd
      rol(i,j) = (pex*hx(i,j) + pey*hy(i,j) + pxy*hh(i,j)
     +           + pvx*wx(i,j) + pvy*wy(i,j))/rett  - ret0*ht(i,j)
      rg2df(i,j) = rage * ht(i,j)
  56  continue
c
      do 3012 i=1,4
      ki = map(inline(i))
      if(ic(ki).eq.1) goto 3012
      ii= ki-lc(ki)
        do 3013 j=1,4
        kj = map(inline(j))
        if(ic(kj).eq.1) goto 3014
        jj = kj - lc(kj)
        if(ii.gt.jj) goto 3015  
c
      do 3083 k=1,nw
      if (ii.ne.ib(jj,k)) go to 3083
      kk=k
      go to 3084
 3083 continue
 3084 continue
c
      iv=ind(jj)+kk-1
      a(iv)= a(iv) + rel(i,j) 
 3015 continue
      fs(ii) = fs(ii) - rol(i,j)*u0(kj) + rg2df(i,j)
      goto 3013
 3014 continue
      f(ii) = f(ii) - (rel(i,j)+rol(i,j))*u0(kj)
c
 3013 continue
 3012 continue
 3011 continue
c
      endif
c
c     get next 2D transport fracture
c     ------------------------------
 555  continue
c
c      write(6,4747) (k,ara(k),bzf(k),k=1,nn/nz)
c4747  format(i5,2e15.5)      
c
c      write (6,695) (i,fs(i),fc(i),i=nz,nn,nz)        !fs undefined 2017
c 695  format(/10x,'i,fs(i),fc(i)...',/(2(i6,2e12.3)))
c
c     source term contribution
c     need to multiply pqq by wwavg*cpf as in 3D blocks
c     (pqq*cq) units are (m^3/s)*(J/m^3.C)*C = J/s
c     --------------------------------------------------
      do 556 i=1,nn
      ii=i-lc(i)
      if(pq(i).le.0.) goto 557
c     tavg=(u1(i)+u0(i))/2.      

c     wwavg = por(1)*sw(1)*wu(tavg,p,q,modelwu,ts)     !check - need element porosity, Sw
c     write(6,34857) 
c    +    por(1),sw(1),wu(tavg,p,q,modelwu,ts),cf,den(tavg,lmass,lage)
34857 format('por(1),sw(1),wu,cf,den: ',5e15.5)
c     cpfavg = (cf*den(tavg,lheat,lmass,lage))            
c     cpfavg = cf            
c     if(lmass.or.lage) then
c              wwavg = 1.d0
c              cpfavg=1.d0
c     endif
      pqq= pq(i)/(por(1))           !check hardwire - need real nodal porosity here 
c     pqq= pq(i)*wwavg*cpfavg       !check - need wwavg*cpfavg  - 
      a(ind(ii)) = a(ind(ii)) + pqq
      f(ii)=f(ii)+pqq*cq(i)
  557 continue
  556 continue
c
c     mobilization source term
c     time-weight spn 
c     -------------------------
c     do l=1,ne
c     if(depd(l).gt.0.d0) then
c     tavg = 0.d0
c     do j=1,8
c      i=in(l,j)
c      tavg = tavg + (t0(i)+t1(i))/2.d0
c     enddo
c     tavg=tavg/8.d0
c     speavg = (spe1(l)+spe0(l))/2.d0
c     tt2 = speavg * (depd(l)*rhob(l)
c    +           * (sqrt(vx(l)**2+vy(l)**2+vz(l)**2))
c    +           * (1.787d-3*rvisc(tavg,lheat,lmass,lage)))
c    +           * vt(l)    ! * volume
c
c     back to nodes
c     --------------
c     do j=1,8
c      i=in(l,j)
c      ii = i-lc(i)
c      f(ii) = f(ii) - tt2/float(icnt(i))
c     enddo
c     
c     endif
c     enddo
c
c     assemble flux vector
c     --------------------
      do 121 i=1,n
      f(i)= f(i) + fs(i) 
  121 continue
c
c     if (it.gt.1) go to 735
c     do 733 i=1,maxna
c 733 write (6,734) i,iaa(i),a(i)
c 734 format (2i5,f14.8)
c 735 continue
c
c     conjugate gradient solution - transport3 mass or age
c     =====================================================
      neu=1
c     call clock@(start)
      call precg (a,n,iaa,ind,f,u2,aa,laa,neu)
c      write(*,123)
c  123 format(10x,'transport pcg solution complete ...')
c     call clock@(finish)
c     tpcgt=tpcgt+finish-start
c
   90 continue
c
c     expand solution vector:
c     -----------------------
      do 141 ii=1,nn
      i = nn-ii+1
      if(ic(i).eq.1) goto 142
      k=i-lc(i)
      if(u2(k).lt.1.d-50) u2(k)=0.d0  !hardwire check for transport
      u2(i)=u2(k)
      goto 141
  142 continue
c     u2(i)=fc(i)
      u2(i)=u0(i)      !Nov2021 using u0 as linking term, u0 decays with decay2
  141 continue
c
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      return
      end
c
c
c ******************************************************************
c ******************************************************************
c ******************************************************************
c
      subroutine gquad2(rex,rey,rez,rxy,rxz,ryz,pex,pey,
     +    pez,pxy,pxz,pyz,pvx,pvy,pvz,dt1,x,y,z,inl,re,ro,l,
     +    maxnn,rtd0,rtd1,cotlavg,ag,hag,tk)
c
c     rtd = thermal retardation term for mass transport, or heat capacity for heat transport
c     generation of isoparametric element matrices
c     full gauss quadrature numerical integration
c     Nov2013: not coded yet for non-linear R (i.e. freeze-thaw)
c     ----------------------------------------------------------
c
      implicit real*8(a-h,o-z)
c
      dimension detj(27),f(8),dgx(8),dgy(8),dgz(8)
      real*8    x(maxnn),y(maxnn),z(maxnn)
      dimension ag(9),hag(9),re(8,8),ro(8,8),inl(8)
      dimension ff(8,27),dx(8,27),dy(8,27),dz(8,27)
c
      m=8
c     number of gauss points
c     ----------------------
c      np=2
c      npx=np
c      npy=np
c      npz=np
c      np2=npx*npy*npz
c      ns=0
c      if (np.eq.3) ns=2
c
c     placement of gauss points
c     -------------------------
      do 310  k=1,2
      zi=ag(k)
      hz=hag(k)
      do 310 j=1,2
      yi=ag(j)
      hy=hag(j)
      do 310 i=1,2
      xi=ag(i)
      hx=hag(i)
      kl=(k-1)*4+(j-1)*2+i
c     if (l.eq.1) write (6,741) i,j,k,kl,xi,yi,zi,hx,hy,hz
c 741 format (/' gauss point',4i5,6f14.7)
c
c     basis functions and derivatives
c     ===============================
      call sp3lin (x,y,z,inl,m,xi,yi,zi,det,l,dgx,dgy,dgz,f,maxnn)
      do 300 jj=1,m
      ff(jj,kl)=f(jj)
      dx(jj,kl)=dgx(jj)
      dy(jj,kl)=dgy(jj)
  300 dz(jj,kl)=dgz(jj)
      detj(kl)=det*hz*hy*hx
  310 continue
c     if (l.ne.1) go to 311
c     write (6,314) l
c 314 format (/' basis functions element',i5)
c     write (6,313) (k,(dx(i,k),i=1,8),k=1,8)
c     write (6,313) (k,(dy(i,k),i=1,8),k=1,8)
c     write (6,313) (k,(dz(i,k),i=1,8),k=1,8)
c 311 continue
c 313 format (i5,8f14.8)
c
c     elemental matrices, scheme 2
c     ----------------------------
      do 51 i=1,m
      do 51 j=1,m
      res=0.
      ros=0.
      do 52 k=1,8
      rd=rex*dx(i,k)*dx(j,k)+rey*dy(i,k)*dy(j,k)+rez*dz(i,k)*dz(j,k)
     1  +rxy*(dx(i,k)*dy(j,k)+dy(i,k)*dx(j,k))
     2  +rxz*(dx(i,k)*dz(j,k)+dz(i,k)*dx(j,k))
     3  +ryz*(dy(i,k)*dz(j,k)+dz(i,k)*dy(j,k))
      rp=pex*dx(i,k)*dx(j,k)+pey*dy(i,k)*dy(j,k)+pez*dz(i,k)*dz(j,k)
     1  +pxy*(dx(i,k)*dy(j,k)+dy(i,k)*dx(j,k))
     2  +pxz*(dx(i,k)*dz(j,k)+dz(i,k)*dx(j,k))
     3  +pyz*(dy(i,k)*dz(j,k)+dz(i,k)*dy(j,k))
      pv=(pvx*dx(j,k)+pvy*dy(j,k)+pvz*dz(j,k))*ff(i,k)
c
c     consistent formulation used here
c     Nov2013: cotlavg is on rhs, as R for mass transport or Co for heat
c     -------------------------------------------------------------------
c     rt=cotlavg*dt1*ff(i,k)*ff(j,k)
c     res=res+( rd + rt)*detj(k)
c     ros=ros+( (rp+pv) - rt)*detj(k)
      rt=dt1*ff(i,k)*ff(j,k)
      res=res+( rd/cotlavg + rt)*detj(k)
      ros=ros+( (rp+pv)/cotlavg - rt)*detj(k)

c      pre Nov2013
c      rt=dt1*ff(i,k)*ff(j,k)
c      res=res+( rd/rtd1 + rt)*detj(k)
c      ros=ros+( (rp+pv)/rtd1 - rt)*detj(k)
   52 continue
c
      re(i,j)=res
      ro(i,j)=ros
   51 continue
  210 continue
c
      return
      end
c
c *********************************************************************
c ******************************************************************
c ******************************************************************
c
      subroutine gquad3(rex,rey,rez,rxy,rxz,ryz,pex,pey,
     +    pez,pxy,pxz,pyz,pvx,pvy,pvz,dt1,x,y,z,inl,re,ro,l,
     +    maxnn,rtd0,rtd1,cotlavg,ag,hag,tk)
c
c     rtd = thermal retardation term for mass transport, or heat capacity for heat transport
c     generation of isoparametric element matrices
c     full gauss quadrature numerical integration
c     Nov2013: not coded yet for non-linear R (i.e. freeze-thaw)
c     ----------------------------------------------------------
c
      implicit real*8(a-h,o-z)
c
      dimension detj(27),f(8),dgx(8),dgy(8),dgz(8)
      real*8    x(maxnn),y(maxnn),z(maxnn)
      dimension ag(9),hag(9),re(8,8),ro(8,8),inl(8)
      dimension ff(8,27),dx(8,27),dy(8,27),dz(8,27)
c
      m=8
c     number of gauss points
c     ----------------------
c      np=2
c      npx=np
c      npy=np
c      npz=np
c      np2=npx*npy*npz
c      ns=0
c      if (np.eq.3) ns=2
c
c     placement of gauss points
c     -------------------------
      do 310  k=1,2
      zi=ag(k)
      hz=hag(k)
      do 310 j=1,2
      yi=ag(j)
      hy=hag(j)
      do 310 i=1,2
      xi=ag(i)
      hx=hag(i)
      kl=(k-1)*4+(j-1)*2+i
c     if (l.eq.1) write (6,741) i,j,k,kl,xi,yi,zi,hx,hy,hz
c 741 format (/' gauss point',4i5,6f14.7)
c
c     basis functions and derivatives
c     ===============================
      call sp3lin (x,y,z,inl,m,xi,yi,zi,det,l,dgx,dgy,dgz,f,maxnn)
      do 300 jj=1,m
      ff(jj,kl)=f(jj)
      dx(jj,kl)=dgx(jj)
      dy(jj,kl)=dgy(jj)
  300 dz(jj,kl)=dgz(jj)
      detj(kl)=det*hz*hy*hx
  310 continue
c     if (l.ne.1) go to 311
c     write (6,314) l
c 314 format (/' basis functions element',i5)
c     write (6,313) (k,(dx(i,k),i=1,8),k=1,8)
c     write (6,313) (k,(dy(i,k),i=1,8),k=1,8)
c     write (6,313) (k,(dz(i,k),i=1,8),k=1,8)
c 311 continue
c 313 format (i5,8f14.8)
c
c     elemental matrices, scheme 2
c     ----------------------------
      do 51 i=1,m
      do 51 j=1,m
      res=0.
      ros=0.
      do 52 k=1,8
      rd=rex*dx(i,k)*dx(j,k)+rey*dy(i,k)*dy(j,k)+rez*dz(i,k)*dz(j,k)
     1  +rxy*(dx(i,k)*dy(j,k)+dy(i,k)*dx(j,k))
     2  +rxz*(dx(i,k)*dz(j,k)+dz(i,k)*dx(j,k))
     3  +ryz*(dy(i,k)*dz(j,k)+dz(i,k)*dy(j,k))
      rp=pex*dx(i,k)*dx(j,k)+pey*dy(i,k)*dy(j,k)+pez*dz(i,k)*dz(j,k)
     1  +pxy*(dx(i,k)*dy(j,k)+dy(i,k)*dx(j,k))
     2  +pxz*(dx(i,k)*dz(j,k)+dz(i,k)*dx(j,k))
     3  +pyz*(dy(i,k)*dz(j,k)+dz(i,k)*dy(j,k))
      pv=(pvx*dx(j,k)+pvy*dy(j,k)+pvz*dz(j,k))*ff(i,k)
c
c     consistent formulation used here
c     Nov2013: cotlavg is on rhs, as R for mass transport or Co for heat
c     -------------------------------------------------------------------
c     rt=cotlavg*dt1*ff(i,k)*ff(j,k)
c     res=res+( rd + rt)*detj(k)
c     ros=ros+( (rp+pv) - rt)*detj(k)
      rt=dt1*ff(i,k)*ff(j,k)
      res=res+( rd/cotlavg + rt)*detj(k)
      ros=ros+( (rp+pv)/cotlavg - rt)*detj(k)

c      pre Nov2013
c      rt=dt1*ff(i,k)*ff(j,k)
c      res=res+( rd/rtd1 + rt)*detj(k)
c      ros=ros+( (rp+pv)/rtd1 - rt)*detj(k)
   52 continue
c
      re(i,j)=res
      ro(i,j)=ros
   51 continue
  210 continue
c
      return
      end
c
c *********************************************************************

c *********************************************************************
c      
c exint2c for thermal transport call to exint2c only ... no 'decay' term for particles
c same as exint2b but no depc,vexp,smax,spe
c later: do same for gquad
c *********************************************************************
c
      subroutine exint2(rex,rey,rez,rxy,rxz,ryz,pex,pey,pez,
     +     pxy,pxz,pyz,pvx,pvy,pvz,hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,
     +     dt1,x,y,z,inl,re,ro,maxnn,htmt,htmf,rtd0,rtd1,cotlavg,
     +     agefx,rg,exl,eyl,ezl,maxne,l,tk,vx,vy,vz)
c
c---
c     rtd = thermal retardation term
c     assembly of element matrices - exact integration
c     rectangular or orthogonal quadrilateral elements
c     --------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension hxx(8,8),hyy(8,8),hzz(8,8),hxy(8,8),hxz(8,8),hyz(8,8),
     + hvx(8,8),hvy(8,8),hvz(8,8),htmt(8,8),htmf(8,8)
      real*8   x(maxnn),y(maxnn),z(maxnn)
      real*8  exl(maxne),eyl(maxne),ezl(maxne)
      real*8  vx(maxne),vy(maxne),vz(maxne)
      dimension ro(8,8),re(8,8),rg(8,8)
      dimension inl(8)
c
      exll=exl(l)
      eyll=eyl(l)
      ezll=ezl(l)
      rex=rex*eyll*ezll/exll/36.
      rey=rey*exll*ezll/eyll/36.
      rez=rez*exll*eyll/ezll/36.
      rxy=rxy*ezll/12.
      rxz=rxz*eyll/12.
      ryz=ryz*exll/12.
      pex=pex*eyll*ezll/exll/36.
      pey=pey*exll*ezll/eyll/36.
      pez=pez*exll*eyll/ezll/36.
      pxy=pxy*ezll/12.
      pxz=pxz*eyll/12.
      pyz=pyz*exll/12.
      pvx=pvx*eyll*ezll/72.
      pvy=pvy*exll*ezll/72.

      pvz=pvz*exll*eyll/72.
c     rdec=decay(l)*exll*eyll*ezll/216.       !removed 2021 ... no decay for heat transport
c      rtt=cotlavg* dt1*exll*eyll*ezll/216.             !split: R or Co (cotlavg) on rhs
       rtt= dt1*exll*eyll*ezll/216.             !split: R or Co (cotlavg) on rhs
      rage=agefx*exll*eyll*ezll/216.
c
c     lumped: use htmf
c     consistent: use htmt
c     ------------------------
      do 156 i=1,8
      do 156 j=1,8
      re(i,j)=( rex*hxx(i,j)+rey*hyy(i,j)+rez*hzz(i,j)
     1       +rxy*hxy(i,j)+rxz*hxz(i,j)+ryz*hyz(i,j))/cotlavg           !cotlavg = thermal or mass R
     2       +rtt*htmt(i,j) 
c    2       +rtt*htmt(i,j) + rdec*htmt(i,j)                            !removed 2021 - no rdec for heat transport
      ro(i,j)=( pex*hxx(i,j)+pey*hyy(i,j)+pez*hzz(i,j)
     1       +pxy*hxy(i,j)+pxz*hxz(i,j)+pyz*hyz(i,j)
     2       +pvx*hvx(i,j)+pvy*hvy(i,j)+pvz*hvz(i,j))/cotlavg
     3       -rtt*htmt(i,j)
c     rg(i,j)=rage*htmf(i,j)
      rg(i,j)=0.d0
  156 continue
c
      return     
      end
c **************************************************************
c *********************************************************************
c
      subroutine exint3(rex,rey,rez,rxy,rxz,ryz,pex,pey,pez,
     + pxy,pxz,pyz,pvx,pvy,pvz,hxx,hyy,hzz,hxy,hxz,hyz,hvx,hvy,hvz,
     + dt1,x,y,z,inl,re,ro,maxnn,htmt,htmf,rtd0,rtd1,cotlavg,decay,
     + agefx,rg,exl,eyl,ezl,maxne,l,tk,depc,depd,vexp,smax,spnavg,
     + vx,vy,vz,rhob,tavg,cavg,por,sdepd,lheat,lmass,lage)
c
c---
c     rtd = thermal retardation term
c     assembly of element matrices - exact integration
c     rectangular or orthogonal quadrilateral elements
c     --------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension hxx(8,8),hyy(8,8),hzz(8,8),hxy(8,8),hxz(8,8),hyz(8,8),
     + hvx(8,8),hvy(8,8),hvz(8,8),htmt(8,8),htmf(8,8)
      real*8   x(maxnn),y(maxnn),z(maxnn)
      real*8  exl(maxne),eyl(maxne),ezl(maxne)
      real*8  decay(maxne),depc(maxne),depd(maxne),vexp(maxne)
      real*8  rhob(maxne),por(maxne),sdepd(maxne),smax(maxne)
      real*8  vx(maxne),vy(maxne),vz(maxne)
      dimension ro(8,8),re(8,8),rg(8,8)
      dimension inl(8)
      logical lheat,lmass,lage
c
c
c     for Madiha: apply deposition coefficient
c     will over-ride linear decay term "decay" if depc > 0
c     need to * (1.-s/smax)
c     -------------------------------------------------------
c     tt1=0.d0
      tt1=decay(l)
      if(depc(l).gt.0.0d0) then 
      tt1 = depc(l) * (1.d0 - (spnavg/smax(l))) 
     +                   * (sqrt(vx(l)**2+vy(l)**2+vz(l)**2))**vexp(l)
      endif
      decay(l) = tt1 

      tt2=0.d0
      if(depd(l).gt.0.0d0) then
      tt2 = spnavg * (depd(l)*rhob(l)
     +           * (sqrt(vx(l)**2+vy(l)**2+vz(l)**2))
     +           * (1.787d-3*rvisc(tavg,lheat,lmass,lage)))/por(l)        ! /porosity
      endif
      sdepd(l) = tt2

c     if(l.eq.1.or.l.eq.20.or.l.eq.3131) 
c    +write(6,8288) l,depc(l),rhob(l),spnavg,smax(l),vz(l),vexp(l),
c    +tavg,rvisc(tavg,lheat,lmass,lage),decay(l),depd(l),sdepd(l)
c8288 format('l,depc(l),rhob(l),spnavg,smax,vz(l),vexp(l),',
c    +       'tavg,rvisc(l),decay(l),depd(l),sdepd(l) : ',i5,11e11.3)
      exll=exl(l)
      eyll=eyl(l)
      ezll=ezl(l)
      rex=rex*eyll*ezll/exll/36.
      rey=rey*exll*ezll/eyll/36.
      rez=rez*exll*eyll/ezll/36.
      rxy=rxy*ezll/12.
      rxz=rxz*eyll/12.
      ryz=ryz*exll/12.
      pex=pex*eyll*ezll/exll/36.
      pey=pey*exll*ezll/eyll/36.
      pez=pez*exll*eyll/ezll/36.
      pxy=pxy*ezll/12.
      pxz=pxz*eyll/12.
      pyz=pyz*exll/12.
      pvx=pvx*eyll*ezll/72.
      pvy=pvy*exll*ezll/72.

      pvz=pvz*exll*eyll/72.
c
c     decay and mobilization:
      rdec=decay(l)*exll*eyll*ezll/216.
      rdepd=sdepd(l)*exll*eyll*ezll/216.

c      rtt=cotlavg* dt1*exll*eyll*ezll/216.             !split: R or Co (cotlavg) on rhs
       rtt= dt1*exll*eyll*ezll/216.             !split: R or Co (cotlavg) on rhs
      rage=agefx*exll*eyll*ezll/216.

c
c     lumped: use htmf
c     consistent: use htmt
c     ------------------------
      do 156 i=1,8
      do 156 j=1,8
      re(i,j)=( rex*hxx(i,j)+rey*hyy(i,j)+rez*hzz(i,j)
     1       +rxy*hxy(i,j)+rxz*hxz(i,j)+ryz*hyz(i,j))/cotlavg            !cotlavg = thermal or mass R
     2       +rtt*htmt(i,j) + rdec*htmt(i,j)
      ro(i,j)=( pex*hxx(i,j)+pey*hyy(i,j)+pez*hzz(i,j)
     1       +pxy*hxy(i,j)+pxz*hxz(i,j)+pyz*hyz(i,j)
     2       +pvx*hvx(i,j)+pvy*hvy(i,j)+pvz*hvz(i,j))/cotlavg
     3       -rtt*htmt(i,j)
      rg(i,j)= (rage+rdepd)*htmf(i,j)                                   !age term + source term from mobilization
  156 continue
c
      return     
      end
c **************************************************************

