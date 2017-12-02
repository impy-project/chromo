      PROGRAM MAIN
      IMPLICIT NONE
      double precision eps10,one,s_rndm,xrndm
      integer nevt,ii
      data eps10,one,nevt /1.d-10,1.d0,100000000/
      call rnd_ini
      print *,'testing random number generator'
      print *,'nevt:',nevt
      DO II=1,NEVT
c         print*, ii
         xrndm = s_rndm(ii)
         if(xrndm.lt.eps10) write(6,*)'ii,s_rndm:',ii,xrndm
         if(xrndm.gt.one-eps10) write(6,*)'ii,s_rndm:',ii,xrndm
      ENDDO
      END
