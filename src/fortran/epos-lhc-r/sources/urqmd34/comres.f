c $Id: comres.f,v 1.15 2003/06/29 14:26:36 weber Exp $
c
cdes This file contains definitions for the collision term
c
      integer maxbar,maxbra,minbar
      integer offmeson,maxmeson,pimeson,maxbrm,minnuc,mindel
      integer maxbrs1,maxbrs2
      integer numnuc,numdel,nucleon,maxnuc,maxdel
      integer minmes,maxmes


      parameter (minnuc=1) ! lowest baryon particle ID 
      parameter (minmes=100) ! lowest meson particle ID
corig      parameter (maxmes=132) ! hightest meson particle ID
c add 5 charmed mesons to array
      parameter (maxmes=139) ! hightest meson particle ID

c number of resonances of a kind
      parameter (numnuc=16) ! number of nucleon resonances
      parameter (numdel=10) ! number of delta resonances
c indices of minimal and maximal itype of a kind (redundant but nice)
      parameter (maxnuc=minnuc+numnuc-1) ! highest nucleon ID
      parameter (mindel=minnuc+maxnuc)   ! lowest delta ID
      parameter (maxdel=mindel+numdel-1) ! highest delta ID

c minres & maxres define the range of nonstable & nonstrange baryons
      integer minres,maxres
      parameter (minres=minnuc+1) ! lowest baryon resonance ID
      parameter (maxres=maxdel)   ! highest (nonstrange) baryon 
                                  ! resonance ID

c strangenes.ne.0 baryon resonances
      integer minlam,minsig,mincas,minome
      integer numlam,numsig,numcas,numome
      integer maxlam,maxsig,maxcas,maxome
      parameter (numlam=13) ! number of lambda states
      parameter (numsig=9)  ! number of sigma states
      parameter (numcas=6)  ! number of cascade states
      parameter (numome=1)  ! number of omega states
      parameter (minlam=mindel+numdel)   ! ID of lowest lambda state
      parameter (maxlam=minlam+numlam-1) ! ID of highest lambda state
      parameter (minsig=minlam+numlam)   ! ID of lowest sigma state
      parameter (maxsig=minsig+numsig-1) ! ID of highest sigma state
      parameter (mincas=minsig+numsig)   ! ID of lowest cascade state
      parameter (maxcas=mincas+numcas-1) ! ID of highest cascade state
      parameter (minome=mincas+numcas)   ! ID of lowest omega state
      parameter (maxome=minome+numome-1) ! ID of highest omega state

c minbar & maxbar define the range of all baryons
      parameter (minbar=minnuc) ! ID of lowest baryon state
      parameter (maxbar=maxome) ! ID of highest baryon state

      parameter (offmeson=minmes) ! offset between zero and lowest 
                                  ! meson state
      parameter (maxmeson=maxmes) ! ID of highest meson state
c... these variables are in principal obsolete and should be exchanged 
c were referenced 

c... avoid hard coded itypes
      integer itrho,itome,iteta,itkaon,itphi,itetapr
      parameter (itkaon=106)   ! ID of kaon
      parameter (itrho=104)    ! ID of rho meson 
      parameter (itome=103)    ! ID of omega meson
      parameter (iteta=102)    ! ID of eta
      parameter (itphi=109)    ! ID of phi
      parameter (itetapr=107)  ! ID of eta'
      parameter (pimeson=101)  ! ID of $\pi$
      parameter (nucleon=minnuc) ! ID of nucleon

      integer itmin,itmax
      parameter (itmin=minnuc)  ! lowest defined ID
      parameter (itmax=maxmes)  ! highest defined ID
c
      parameter (maxbra=11)  ! decay channels for $s=0$ baryon resonances
corig      parameter (maxbrm=25) ! decay channels for meson resonances
c add 2 decay channels for D*
      parameter (maxbrm=27) ! decay channels for meson resonances
      parameter (maxbrs1=10)! decay channels for $s=1$ baryon resonances
      parameter (maxbrs2=3) ! decay channels for $s=2$ baryon resonances

c 
       integer mlt2it(maxmes-minmes) ! meson IDs sorted by multipletts


      real*8 massoff,mresmin,mresmax
      parameter (massoff=1d-4)      ! offset for mass generation
      parameter (mresmin=1.0765d0)  ! minimum baryon resonance mass
      parameter (mresmax=3.5d0)       ! maximum baryon resonance mass

      character*45 versiontag
      common /versioning/ versiontag

      real*8 massres(minbar:maxbar),widres(minbar:maxbar)
      real*8 branmes(0:maxbrm,minmes+1:maxmes)
      real*8 branres(0:maxbra,minnuc+1:maxdel)
      real*8 branbs1(0:maxbrs1,minlam+1:maxsig)
      real*8 branbs2(0:maxbrs2,mincas+1:maxcas)
      integer Jres(minbar:maxbar)
      integer Jmes(minmes:maxmes)
      integer pares(minbar:maxbar),pames(minmes:maxmes)
      integer Isores(minbar:maxbar), Isomes(minmes:maxmes)
      integer brtype(4,0:maxbra),bmtype(4,0:maxbrm)
      integer bs1type(4,0:maxbrs1),bs2type(4,0:maxbrs2)
      real*8 massmes(minmes:maxmes)
      real*8 mmesmn(minmes:maxmes)
      real*8 widmes(minmes:maxmes)
      integer strres(minbar:maxbar),strmes(minmes:maxmes)
      integer chrmres(minbar:maxbar),chrmmes(minmes:maxmes)

      integer lbr(0:maxbra,minnuc+1:maxdel)
      integer lbs1(0:maxbrs1,minlam+1:maxsig)
      integer lbs2(0:maxbrs2,mincas+1:maxcas)
      integer lbm(0:maxbrm,minmes+1:maxmes)

      common /resonances/ massres,widres,massmes,widmes,mmesmn,
     ,                    branres,branmes,branbs1,branbs2,
     ,                    bs1type,bs2type,lbs1,lbs2,lbm,
     ,                    jres,jmes,lbr,brtype,pares,pames,
     ,                    bmtype,
     ,                    Isores,Isomes,strres,strmes,mlt2it,
     ,                    chrmres,chrmmes

c     massres   : baryon mass table
c     widres    : baryon decay width table
c     massmes   : meson mass table
c     widmes    : meson decay width table
c     mmesmn    : table of minimum masses for meson resonances
c     branres   : branching ratios for $s=0$ baryon resonances
c     branmes   : branching ratios for meson resonances
c     branbs1   : branching ratios for $s=1$ baryon resonances
c     branbs2   : branching ratios for $s=2$ baryon resonances
c     brtype    : definitions of the decay branches for $s=0$ baryon resonances
c     bmtype    : definitions of the decay branches for meson resonances
c     bs1type   : definitions of the decay branches for $s=1$ baryon resonances
c     bs2type   : definitions of the decay branches for $s=2$ baryon resonances
c     lbr       : decay angular momenta for $s=0$ baryon decays
c     lbm       : decay angular momenta for meson decays
c     lbs1      : decay angular momenta for $s=1$ baryon decays
c     lbs2      : decay angular momenta for $s=2$ baryon decays
c     jres      : spin table for baryons
c     jmes      : spin table for mesons
c     pares     : parity table for baryons
c     pames     : parity table for mesons
c     isores    : isospin table for baryons
c     isomes    : isospin table for mesons
c     strres    : strangeness table for baryons
c     strmes    : strangeness table for mesons 
c     
c
ccccccccccccccccccccc sigtab-declarations cccccccccccccccccccccccccccccccccc

      integer itblsz,nsigs,maxreac,maxpsig,sigver
c     VERSION NUMBER of SIGTAB
      parameter (sigver = 1000) ! version number for collision 
                                ! term and tables
ccccccccccccccccccccccccccccccccccccccc
c

      parameter (maxreac = 15) ! maximum number of collision classes
      parameter (maxpsig = 21) ! maximum number of cross 
                               ! sections per class
      parameter (nsigs = 10)   ! number of tabulated cross sections

      PARAMETER (ITBLSZ= 100)  ! table size of cross section tables


c
      integer sigmaLN(maxpsig,2,maxreac)
      integer sigmainf(nsigs,20)
      real*8 sigmascal(nsigs,5), sigmas(nsigs,itblsz)
c
      common/sigtabi/sigmaln,sigmainf
      common/sigtabr/sigmas,sigmascal

c     sigmaln   : pointer array connecting collision classes with cross sections
c     sigmainf  : information related to tabulated cross sections
c     sigmascal : information related to tabulated cross sections
c     sigmas    : tabulated cross sections

