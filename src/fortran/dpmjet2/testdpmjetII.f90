      Program testdpmjetII
!**********************************************************************

!     example program calling PHOJET without reading input cards

!**********************************************************************
      Implicit Double Precision (A-H, O-Z)
      Character *132 datdir
      Save

      datdir = '/lustre/fs17/group/that/af/m2m/iamdata/'
      Call dpmjin(12345, datdir)
!  number of events
      neve = 200000
      itry = 0

      Do k = 1, neve
        If (mod(k, 1000) == 0) Write(6,*) 'event', k
        Call dpmjet2_event(158.D0,1,1,1,1,1,2)
      End Do
      End Program testdpmjetII