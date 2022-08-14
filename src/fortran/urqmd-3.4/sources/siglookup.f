c $Id: siglookup.f,v 1.4 1999/01/18 09:57:15 ernst Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      real*8 function SIGLOOKUP(iline,sqrts)
c
c     Revision : 1.0
c
cinput iline : line\# for information in {\tt sigmainf} array
cinput sqrts : $\sqrt{s}$ of collision
c
c output : returns tabulated cross section value
c
c     This function returns the cross section stored in line ISIGLINE of
c     the SIGMAS array at the respective value of sqrt(s) (SQRTS)
c     Optional scaling according to SIGMASCAL is performed.
C
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none

      real*8 sqrts,xpt,x1,y1in,y2in,slin,tlow,tblstp,thigh
      integer isigline,index,iline
c
c later, sigtab.f should only be included into main.f and all relevant tables
c should be accessed via common blocks
      include 'comres.f'
c
c line for sigmas(line,*) array
      isigline=sigmainf(iline,1)
c
c log of lowest sqrt(s) entry in the table
      tlow=log(sigmascal(iline,2))

c
c exp(tblstp) is Delta(sqrts) for table
      tblstp=sigmascal(iline,3)
c
c maximum log of sqrt(s) enty in the table
      thigh=itblsz*tblstp+tlow
c
c table-lookup sequence (modified IQMD/RQMD)
c first generate x-coordinate from sqrt(s) for the table
      XPT=LOG(sqrts) 
c check if sqrt(s) is below lower boundary
      IF(XPT.LT.TLOW) THEN
         siglookup=0.d0
         return
      END IF
c now get index (which is the column in sigtab)
      INDEX=INT((XPT-TLOW)/TBLSTP) + 1
c check if sqrt(s) is above upper boundary
      IF(INDEX.GE.ITBLSZ) THEN
c also this solution is not clean...
         INDEX=ITBLSZ-1
         xpt=thigh
      END IF
C FIND SLOPES AND CROSSECTIONS
c now make a straighforward interpolation
c in the sigtab array
c bracket for sqrt(s) value (between X1 and XPT)
      X1=(INDEX-1)*TBLSTP + TLOW 
c sigmas(isigline,INDEX) is the c.s. array, the wanted c.s. is stored
c in line isigline
      Y1IN=sigmas(isigline,INDEX)
      Y2IN=sigmas(isigline,INDEX+1)
c get the slope
      SLIN=(Y2IN-Y1IN)/TBLSTP
c get the cross section and store it in SIGLOOKUP 
      SIGLOOKUP=SLIN*(XPT-X1) + Y1IN
c scale the cross section
      siglookup=sigmascal(iline,1)*siglookup
      return
      end
