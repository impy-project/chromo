c $Id: comres.f,v 1.15 2003/06/29 14:26:36 weber Exp $
c
cdes This file contains definitions for the collision term
c
      integer maxbar,maxbra,minbar
      integer offmeson,maxmeson,pimeson,maxbrm,minnuc,mindel
      integer maxbrs1,maxbrs2
      integer numnuc,numdel,nucleon,maxnuc,maxdel
      integer minmes,maxmes


      parameter (minnuc=1)
      parameter (minmes=100)
corig      parameter (maxmes=132)
c add 5 charmed mesons to array
      parameter (maxmes=139)

c number of resonances of a kind
      parameter (numnuc=16)
      parameter (numdel=10)
c indices of minimal and maximal itype of a kind (redundant but nice)
      parameter (maxnuc=minnuc+numnuc-1)
      parameter (mindel=minnuc+maxnuc)
      parameter (maxdel=mindel+numdel-1)

c minres & maxres define the range of nonstable & nonstrange baryons
      integer minres,maxres
      parameter (minres=minnuc+1)
      parameter (maxres=maxdel)
                                  ! resonance ID

c strangenes.ne.0 baryon resonances
      integer minlam,minsig,mincas,minome
      integer numlam,numsig,numcas,numome
      integer maxlam,maxsig,maxcas,maxome
      parameter (numlam=13)
      parameter (numsig=9)
      parameter (numcas=6)
      parameter (numome=1)
      parameter (minlam=mindel+numdel)
      parameter (maxlam=minlam+numlam-1)
      parameter (minsig=minlam+numlam)
      parameter (maxsig=minsig+numsig-1)
      parameter (mincas=minsig+numsig)
      parameter (maxcas=mincas+numcas-1)
      parameter (minome=mincas+numcas)
      parameter (maxome=minome+numome-1)

c minbar & maxbar define the range of all baryons
      parameter (minbar=minnuc)
      parameter (maxbar=maxome)

      parameter (offmeson=minmes)
                                  ! meson state
      parameter (maxmeson=maxmes)
c... these variables are in principal obsolete and should be exchanged 
c were referenced 

c... avoid hard coded itypes
      integer itrho,itome,iteta,itkaon,itphi,itetapr
      parameter (itkaon=106)
      parameter (itrho=104)
      parameter (itome=103)
      parameter (iteta=102)
      parameter (itphi=109)
      parameter (itetapr=107)
      parameter (pimeson=101)
      parameter (nucleon=minnuc)

      integer itmin,itmax
      parameter (itmin=minnuc)
      parameter (itmax=maxmes)
c
      parameter (maxbra=11)
corig      parameter (maxbrm=25)
c add 2 decay channels for D*
      parameter (maxbrm=27)
      parameter (maxbrs1=10)
      parameter (maxbrs2=3)

c 
       integer mlt2it(maxmes-minmes)


      real*8 massoff,mresmin,mresmax
      parameter (massoff=1d-4)
      parameter (mresmin=1.0765d0)
      parameter (mresmax=3.5d0)

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
      parameter (sigver = 1000)
                                ! term and tables
ccccccccccccccccccccccccccccccccccccccc
c

      parameter (maxreac = 15)
      parameter (maxpsig = 21)
                               ! sections per class
      parameter (nsigs = 10)

      PARAMETER (ITBLSZ= 100)


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

