C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C


c#############################################################################
c                   Managing particle overflow (i > mxptl)
c#############################################################################
c particles  are  dumped to  and restored from  cccptl  (C++ object)
c dump partile i to cccptl[i-1-mxptl] 
c restore partile i from cccptl[i-1-mxptl]  
c#############################################################################

c-----------------------------------------------------------------------
c Before calling dumpcccptl,restorecccptl  
c one needs  
c                            call createcccptl(n)  via    call checkcccptl
c and at the end
c                            call destroycccptl(n)
c 
c-----------------------------------------------------------------------

      subroutine dumpcccptl(j,i)
      idum1=i
      idum2=j
      end
      
      subroutine restorecccptl(i,j)
      idum1=i
      idum2=j
      end
      
      subroutine checkcccptl(n)            
      include "epos.inc"
      if(n.gt.mxptl)call utstop("checkccptl: n larger than mxptl !&")
      end

      subroutine closecccptl     
      end     
      
      !-----------------------------------------------------
      !     set
      !-----------------------------------------------------

      subroutine setnptl(ival)
      include "epos.inc"
      nptl=ival
      end
     
      subroutine setistptl(i,ival)
      include "epos.inc"
        istptl(i)=ival
      end

      subroutine setidptl(i,ival)            
      include "epos.inc"
        idptl(i)=ival
      end

      subroutine setityptl(i,ival)            
      include "epos.inc"
        ityptl(i)=ival
      end

      subroutine setiorptl(i,ival)            
      include "epos.inc"
        iorptl(i)=ival
      end

      subroutine setjorptl(i,ival)            
      include "epos.inc"
        jorptl(i)=ival
      end

      subroutine setifrptl(i,i1,i2)            
      include "epos.inc"
        ifrptl(1,i)=i1
        ifrptl(2,i)=i2
      end

      subroutine setibptl(i,i1,i2,i3,i4)            
      include "epos.inc"
        ibptl(1,i)=i1
        ibptl(2,i)=i2
        ibptl(3,i)=i3
        ibptl(4,i)=i4
      end
      subroutine set2ibptl(i,j,ival)            
      include "epos.inc"
        ibptl(j,i)=ival
      end

      subroutine setdesptl(i,val)            
      include "epos.inc"
        desptl(i)=val
      end

      subroutine setradptl(i,val)            
      include "epos.inc"
        radptl(i)=val
      end
      
      subroutine setrinptl(i,val)            
      include "epos.inc"
        rinptl(i)=val
      end
      
      subroutine setpptl(i,p1,p2,p3,p4,p5)
      include "epos.inc"
         pptl(1,i)=p1
         pptl(2,i)=p2
         pptl(3,i)=p3
         pptl(4,i)=p4
         pptl(5,i)=p5
      end

      subroutine setxorptl(i,x1,x2,x3,x4)
      include "epos.inc"
         xorptl(1,i)=x1
         xorptl(2,i)=x2
         xorptl(3,i)=x3
         xorptl(4,i)=x4
      end

      subroutine settivptl(i,t1,t2)
      include "epos.inc"
         tivptl(1,i)=t1
         tivptl(2,i)=t2
      end
      
      subroutine setzpaptl(i,t1,t2)
      include "epos.inc"
         zpaptl(1,i)=t1
         zpaptl(2,i)=t2
      end

      !-----------------------------------------------------
      !     get
      !-----------------------------------------------------

      subroutine getnptl(ival)
      include "epos.inc"
      ival=nptl
      end

      subroutine getistptl(i,ival)
      include "epos.inc"
        ival=istptl(i)
      end

      subroutine getidptl(i,ival)
      include "epos.inc"
         ival=idptl(i)
      end

      subroutine getityptl(i,ival)            
      include "epos.inc"
         ival=ityptl(i)
      end

      subroutine getiorptl(i,ival)            
      include "epos.inc"
         ival=iorptl(i)
      end

      subroutine getjorptl(i,ival)            
      include "epos.inc"
         ival=jorptl(i)
      end

      subroutine getifrptl(i,i1,i2)            
      include "epos.inc"
         i1=ifrptl(1,i)
         i2=ifrptl(2,i)
      end

      subroutine getibptl(i,i1,i2,i3,i4)            
      include "epos.inc"
         i1=ibptl(1,i)
         i2=ibptl(2,i)
         i3=ibptl(3,i)
         i4=ibptl(4,i)
      end
      subroutine get2ibptl(i,j,ival)            
      include "epos.inc"
         ival=ibptl(j,i)
      end

      subroutine getdesptl(i,val)              
      include "epos.inc"
         val=desptl(i)
      end

      subroutine getradptl(i,val)              
      include "epos.inc"
         val=radptl(i)
      end

      subroutine getrinptl(i,val)
      include "epos.inc"
        val=rinptl(i)
      end

      subroutine getpptl(i,p1,p2,p3,p4,p5)
      include "epos.inc"
         p1=pptl(1,i)
         p2=pptl(2,i)
         p3=pptl(3,i)
         p4=pptl(4,i)
         p5=pptl(5,i)
      end

      subroutine getxorptl(i,x1,x2,x3,x4)
      include "epos.inc"
         x1=xorptl(1,i)
         x2=xorptl(2,i)
         x3=xorptl(3,i)
         x4=xorptl(4,i)
      end

      subroutine gettivptl(i,t1,t2)
      include "epos.inc"
         t1=tivptl(1,i)
         t2=tivptl(2,i)
      end

      subroutine getzpaptl(i,t1,t2)
      include "epos.inc"
         t1=zpaptl(1,i)
         t2=zpaptl(2,i)
      end

