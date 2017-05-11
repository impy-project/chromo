
      DOUBLE PRECISION FUNCTION sigela_pn(plab)
C-----------------------------------------------------------------------
C
C     low-energy pn/np elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 02/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  pn elastic cross section
      DATA (PTPP(K),K=    1,   18) /
     &  -1.0128E+00,-8.8365E-01,-7.8000E-01,-6.8973E-01,-5.7462E-01,
     &  -4.2138E-01,-2.9384E-01,-1.1581E-01, 1.1309E-01, 5.3273E-01,
     &   9.6497E-01, 1.4860E+00, 2.0449E+00, 2.6798E+00, 3.5939E+00,
     &   4.9903E+00, 6.2215E+00, 6.8942E+00/
      DATA (STPP(K),K=    1,   18) /
     &1.0001E+02,8.2414E+01,6.5819E+01,5.4660E+01,4.7794E+01,4.0500E+01,
     &3.5781E+01,3.3208E+01,2.9921E+01,2.3919E+01,1.8633E+01,1.4206E+01,
     &1.1068E+01,9.0752E+00,7.5167E+00,6.6817E+00,6.8455E+00,6.8568E+00/


C  initialize cross section tables

      if(init) then
        N = 18
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGELA_PN: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_pn = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGELA_PN: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_pn = FV(1) 

      END


      DOUBLE PRECISION FUNCTION sigela_pp(plab)
C-----------------------------------------------------------------------
C
C     low-energy pp elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 02/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  pp elastic cross section
      DATA (PTPP(K),K=    1,   20) /
     &  -1.0548E+00,-9.9070E-01,-8.2516E-01,-6.8608E-01,-4.7199E-01,
     &  -2.7085E-01,-1.0784E-01,-7.6152E-03, 1.6806E-01, 3.3154E-01,
     &   5.4551E-01, 8.2275E-01, 1.3768E+00, 2.0058E+00, 2.9862E+00,
     &   3.7151E+00, 4.3182E+00, 5.1348E+00, 5.6750E+00, 6.2152E+00/
      DATA (STPP(K),K=    1,   20) /
     &4.2555E+01,3.7310E+01,2.8426E+01,2.4873E+01,2.2758E+01,2.2166E+01,
     &2.3350E+01,2.4450E+01,2.5212E+01,2.4535E+01,2.2927E+01,1.9459E+01,
     &1.4213E+01,1.0745E+01,8.4602E+00,7.3604E+00,6.8528E+00,6.6836E+00,
     &6.6836E+00,6.6836E+00/


C  initialize cross section tables

      if(init) then
        N = 20
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGELA_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_pp = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGELA_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_pp = FV(1) 

      END


      DOUBLE PRECISION FUNCTION sigtot_pn(plab)
C-----------------------------------------------------------------------
C
C     low-energy pn and np total cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 02/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  pn total cross section
      DATA (PTPP(K),K=    1,   17) /
     &  -1.0129E+00,-8.4520E-01,-7.4136E-01,-5.3626E-01,-3.3210E-01,
     &  -1.2859E-01, 8.7237E-02, 3.1519E-01, 6.7022E-01, 1.0889E+00,
     &   1.5714E+00, 2.0792E+00, 2.6760E+00, 3.9453E+00, 4.9226E+00,
     &   5.6207E+00, 6.7629E+00/
      DATA (STPP(K),K=    1,   17) /
     &1.0000E+02,7.9053E+01,6.0976E+01,4.5194E+01,3.6729E+01,3.3429E+01,
     &3.3142E+01,3.7303E+01,4.0316E+01,4.1607E+01,4.0746E+01,3.9885E+01,
     &3.8594E+01,3.8307E+01,3.8881E+01,3.9168E+01,4.1320E+01/


C  initialize cross section tables

      if(init) then
        N = 17
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGTOT_PN: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_pn = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGTOT_PN: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_pn = FV(1) 

      END


      DOUBLE PRECISION FUNCTION sigtot_pp(plab)
C-----------------------------------------------------------------------
C
C     low-energy pp 
C     (based on spline interpolations)
C
C                                              (R.Engel 02/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  pp total cross section
      DATA (PTPP(K),K=    1,   23) /
     &  -1.4202E+00,-1.2583E+00,-1.0464E+00,-8.3253E-01,-6.0471E-01,
     &  -3.6376E-01,-8.4289E-02, 6.8739E-02, 1.9666E-01, 3.2471E-01,
     &   4.2673E-01, 5.5375E-01, 7.5675E-01, 1.0737E+00, 1.5176E+00,
     &   2.1393E+00, 2.7230E+00, 3.5353E+00, 4.3223E+00, 5.1728E+00,
     &   5.7949E+00, 6.2392E+00, 6.9122E+00/
      DATA (STPP(K),K=    1,   23) /
     &9.2081E+01,7.0000E+01,4.2437E+01,2.8579E+01,2.3858E+01,2.2335E+01,
     &2.3858E+01,2.8883E+01,3.5888E+01,4.3807E+01,4.7157E+01,4.7766E+01,
     &4.7157E+01,4.4569E+01,4.1523E+01,3.9695E+01,3.8782E+01,3.8173E+01,
     &3.8173E+01,3.8477E+01,3.9391E+01,4.0000E+01,4.1523E+01/


C  initialize cross section tables

      if(init) then
        N = 23
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGTOT_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_pp = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGTOT_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_pp = FV(1) 

      END


      DOUBLE PRECISION FUNCTION sigela_pipp(plab)
C-----------------------------------------------------------------------
C
C     low-energy pi+p elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  pi+p elastic cross section
      DATA (PTPP(K),K=    1,   24) /
     &  -9.1117E-01,-8.4887E-01,-7.8656E-01,-6.6196E-01,-5.3736E-01,
     &  -4.4390E-01,-3.6083E-01,-2.6738E-01,-1.8431E-01,-5.9706E-02,
     &   5.4515E-02, 1.3758E-01, 2.4142E-01, 3.5564E-01, 4.0756E-01,
     &   5.1140E-01, 6.9830E-01, 1.0410E+00, 1.6225E+00, 2.2455E+00,
     &   2.9620E+00, 3.7407E+00, 4.6026E+00, 5.5163E+00/
      DATA (STPP(K),K=    1,   24) /
     &7.3812E+01,5.8453E+01,4.5967E+01,3.1602E+01,2.2652E+01,1.6133E+01,
     &1.2044E+01,9.2818E+00,8.3978E+00,9.9448E+00,1.2818E+01,1.4144E+01,
     &1.6354E+01,1.8011E+01,1.7238E+01,1.2928E+01,1.0055E+01,7.1823E+00,
     &5.5249E+00,4.6409E+00,3.6464E+00,2.9834E+00,3.2044E+00,3.0939E+00/


C  initialize cross section tables

      if(init) then
        N = 24
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGTOT_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_pipp = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGTOT_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_pipp = FV(1) 

      END


      DOUBLE PRECISION FUNCTION sigtot_pipp(plab)
C-----------------------------------------------------------------------
C
C     low-energy pi+p total cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  pi+p total cross section
      DATA (PTPP(K),K=    1,   37) /
     &  -9.2155E-01,-8.6963E-01,-8.0733E-01,-7.2426E-01,-5.4774E-01,
     &  -4.7505E-01,-4.1275E-01,-3.6083E-01,-3.0891E-01,-2.2585E-01,
     &  -1.7393E-01,-8.0473E-02, 2.3363E-02, 1.5835E-01, 2.3104E-01,
     &   2.9334E-01, 3.1411E-01, 3.5564E-01, 4.1794E-01, 4.2833E-01,
     &   4.9063E-01, 5.7370E-01, 6.7754E-01, 7.2945E-01, 8.1252E-01,
     &   8.8521E-01, 9.9943E-01, 1.1033E+00, 1.4044E+00, 1.7782E+00,
     &   2.1313E+00, 2.6712E+00, 3.2942E+00, 3.8342E+00, 4.6441E+00,
     &   5.4748E+00, 5.8382E+00/
      DATA (STPP(K),K=    1,   37) /
     &7.3812E+01,6.4420E+01,5.0939E+01,3.7790E+01,2.3867E+01,1.8674E+01,
     &1.6022E+01,1.5138E+01,1.4365E+01,1.5138E+01,1.7127E+01,2.0773E+01,
     &2.4420E+01,2.7845E+01,3.3591E+01,3.9116E+01,4.0773E+01,4.1215E+01,
     &4.0000E+01,3.8232E+01,3.3370E+01,3.0608E+01,2.9061E+01,2.8619E+01,
     &2.9834E+01,3.0829E+01,3.0497E+01,2.9061E+01,2.7514E+01,2.5746E+01,
     &2.4862E+01,2.3646E+01,2.3094E+01,2.2873E+01,2.3204E+01,2.3978E+01,
     &2.4420E+01/


C  initialize cross section tables

      if(init) then
        N = 37
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGTOT_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_pipp = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGTOT_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_pipp = FV(1) 

      END




      DOUBLE PRECISION FUNCTION sigela_pimp(plab)
C-----------------------------------------------------------------------
C
C     low-energy pi-p elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  pi-p elastic cross section
      DATA (PTPP(K),K=    1,   56) /
     &  -1.8980E+00,-1.5458E+00,-1.4323E+00,-1.3602E+00,-1.2880E+00,
     &  -1.2571E+00,-1.1845E+00,-1.1531E+00,-1.1112E+00,-1.0691E+00,
     &  -1.0063E+00,-9.1252E-01,-8.2935E-01,-7.0477E-01,-6.0118E-01,
     &  -4.6652E-01,-4.1489E-01,-3.9435E-01,-3.6334E-01,-3.4267E-01,
     &  -3.0100E-01,-2.6966E-01,-2.4866E-01,-2.1741E-01,-1.6542E-01,
     &  -1.1357E-01,-9.2992E-02,-8.2923E-02,-4.1875E-02,-1.1054E-02,
     &   3.0281E-02, 7.2145E-02, 8.2958E-02, 1.1458E-01, 1.5645E-01,
     &   2.6051E-01, 3.4368E-01, 3.8539E-01, 4.7900E-01, 5.3080E-01,
     &   6.3455E-01, 7.4898E-01, 9.1527E-01, 1.1023E+00, 1.3412E+00,
     &   1.5594E+00, 1.9541E+00, 2.4007E+00, 2.7122E+00, 3.0653E+00,
     &   3.4392E+00, 3.8130E+00, 4.2387E+00, 5.0175E+00, 5.3602E+00,
     &   5.8897E+00/
      DATA (STPP(K),K=    1,   56) /
     &2.9793E+00,9.7103E+00,1.5007E+01,1.9862E+01,2.3393E+01,2.5269E+01,
     &2.6041E+01,2.4276E+01,2.1076E+01,1.6772E+01,1.3021E+01,1.0372E+01,
     &9.6000E+00,9.8207E+00,1.1697E+01,1.4234E+01,1.6441E+01,1.8207E+01,
     &1.9310E+01,2.0083E+01,1.8979E+01,1.7545E+01,1.5779E+01,1.5007E+01,
     &1.4455E+01,1.5007E+01,1.6441E+01,1.8869E+01,2.2621E+01,2.5159E+01,
     &2.6703E+01,2.4166E+01,2.0855E+01,1.7214E+01,1.4676E+01,1.2910E+01,
     &1.2138E+01,1.0814E+01,9.6000E+00,1.0483E+01,1.1145E+01,9.6000E+00,
     &8.3862E+00,7.5034E+00,6.6207E+00,6.0690E+00,4.9655E+00,4.4138E+00,
     &4.4138E+00,3.7517E+00,3.3103E+00,3.2000E+00,3.3103E+00,3.3103E+00,
     &3.3103E+00,3.5310E+00/


C  initialize cross section tables

      if(init) then
        N = 56
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGTOT_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_pimp = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGTOT_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_pimp = FV(1) 

      END


      DOUBLE PRECISION FUNCTION sigtot_pimp(plab)
C-----------------------------------------------------------------------
C
C     low-energy pi-p total cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  pi-p total cross section
      DATA (PTPP(K),K=    1,   53) /
     &  -1.9302E+00,-1.8269E+00,-1.6617E+00,-1.5490E+00,-1.4577E+00,
     &  -1.3146E+00,-1.2630E+00,-1.2211E+00,-1.1686E+00,-1.1364E+00,
     &  -1.0937E+00,-1.0305E+00,-9.4645E-01,-8.5245E-01,-7.6915E-01,
     &  -6.7584E-01,-5.2057E-01,-4.3813E-01,-4.0781E-01,-3.6669E-01,
     &  -3.1507E-01,-2.8372E-01,-2.6240E-01,-2.0995E-01,-1.7861E-01,
     &  -1.1661E-01,-9.6329E-02,-7.6149E-02,-3.5817E-02,-5.0811E-03,
     &   1.5958E-02, 5.8095E-02, 1.1175E-01, 1.7444E-01, 1.9540E-01,
     &   2.8868E-01, 3.7173E-01, 4.5500E-01, 5.4845E-01, 6.4176E-01,
     &   7.1436E-01, 8.3919E-01, 9.6397E-01, 1.3069E+00, 1.7018E+00,
     &   2.0447E+00, 2.5952E+00, 3.1249E+00, 3.6130E+00, 4.1426E+00,
     &   4.8175E+00, 5.3159E+00, 5.9284E+00/
      DATA (STPP(K),K=    1,   53) /
     &1.1145E+01,1.5007E+01,2.2179E+01,3.4428E+01,5.0428E+01,6.7862E+01,
     &7.0952E+01,6.7972E+01,6.3007E+01,5.5393E+01,4.6566E+01,3.9614E+01,
     &3.1779E+01,2.7586E+01,2.5821E+01,2.6924E+01,3.0676E+01,3.5531E+01,
     &4.1931E+01,4.5131E+01,4.7448E+01,4.5903E+01,4.1600E+01,3.7517E+01,
     &3.6083E+01,3.8400E+01,4.2152E+01,4.6676E+01,5.5945E+01,5.9145E+01,
     &5.7048E+01,5.2414E+01,3.9062E+01,3.6083E+01,3.4538E+01,3.5862E+01,
     &3.6083E+01,3.4538E+01,3.4538E+01,3.5641E+01,3.6303E+01,3.4538E+01,
     &3.3214E+01,3.1117E+01,2.8690E+01,2.7145E+01,2.5600E+01,2.4717E+01,
     &2.4166E+01,2.4166E+01,2.3945E+01,2.4055E+01,2.5159E+01/


C  initialize cross section tables

      if(init) then
        N = 53
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGTOT_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_pimp = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGTOT_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_pimp = FV(1) 

      END


      DOUBLE PRECISION FUNCTION sigela_kpp(plab)
C-----------------------------------------------------------------------
C
C     low-energy K+p elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  K+p elastic cross section
      DATA (PTPP(K),K=    1,   22) /
     &  -1.1500E+00,-8.0733E-01,-5.4774E-01,-4.1275E-01,-2.5700E-01,
     &  -8.0474E-02, 7.5281E-02, 2.5180E-01, 3.7641E-01, 5.3216E-01,
     &   6.8792E-01, 8.4368E-01, 1.0929E+00, 1.5913E+00, 1.9340E+00,
     &   2.3182E+00, 2.8166E+00, 3.2215E+00, 3.4708E+00, 3.9276E+00,
     &   4.6233E+00, 5.5475E+00/
      DATA (STPP(K),K=    1,   22) /
     &1.2227E+01,1.2570E+01,1.2499E+01,1.2498E+01,1.2428E+01,1.2012E+01,
     &1.1183E+01,1.0284E+01,9.4544E+00,8.2796E+00,6.8977E+00,5.9300E+00,
     &4.6854E+00,3.6461E+00,3.2293E+00,3.0193E+00,2.6704E+00,2.4602E+00,
     &2.3203E+00,2.0407E+00,2.2426E+00,2.5809E+00/


C  initialize cross section tables

      if(init) then
        N = 22
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGTOT_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_kpp = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGTOT_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_kpp = FV(1) 

      END


      DOUBLE PRECISION FUNCTION sigtot_kpp(plab)
C-----------------------------------------------------------------------
C
C     low-energy K+p total cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  K+p total cross section
      DATA (PTPP(K),K=    1,   20) /
     &  -1.0981E+00,-7.1388E-01,-4.7505E-01,-3.1930E-01,-1.7393E-01,
     &  -8.0474E-02, 2.3363E-02, 9.6049E-02, 1.9989E-01, 3.2449E-01,
     &   4.6986E-01, 6.2562E-01, 8.3329E-01, 1.0825E+00, 1.4355E+00,
     &   2.1001E+00, 2.6920E+00, 3.5434E+00, 4.6337E+00, 5.7448E+00/
      DATA (STPP(K),K=    1,   20) /
     &1.2158E+01,1.2362E+01,1.2429E+01,1.2428E+01,1.3187E+01,1.4429E+01,
     &1.5809E+01,1.7327E+01,1.8224E+01,1.8430E+01,1.7945E+01,1.7806E+01,
     &1.7459E+01,1.7250E+01,1.7041E+01,1.7381E+01,1.7446E+01,1.7853E+01,
     &1.8881E+01,2.0529E+01/


C  initialize cross section tables

      if(init) then
        N = 20
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGTOT_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_kpp = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGTOT_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_kpp = FV(1) 

      END




      DOUBLE PRECISION FUNCTION sigela_kmp(plab)
C-----------------------------------------------------------------------
C
C     low-energy K-p elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  K-p elastic cross section
      DATA (PTPP(K),K=    1,   36) /
     &  -1.7871E+00,-1.4709E+00,-1.2813E+00,-1.1867E+00,-1.0179E+00,
     &  -8.8055E-01,-8.0666E-01,-7.9648E-01,-7.7560E-01,-6.5951E-01,
     &  -5.6450E-01,-4.7995E-01,-3.9539E-01,-3.4256E-01,-2.7894E-01,
     &  -2.4691E-01,-2.0439E-01,-1.1952E-01,-1.3598E-02, 6.0479E-02,
     &   1.1311E-01, 1.4462E-01, 2.0784E-01, 2.6053E-01, 3.2387E-01,
     &   4.4022E-01, 5.5672E-01, 6.9424E-01, 8.6348E-01, 1.2127E+00,
     &   1.6678E+00, 2.3770E+00, 3.2133E+00, 3.9226E+00, 4.6425E+00,
     &   5.1612E+00/
      DATA (STPP(K),K=    1,   36) /
     &6.8962E+01,5.6135E+01,4.7307E+01,4.0271E+01,3.5582E+01,3.2549E+01,
     &3.0480E+01,2.6617E+01,2.3858E+01,2.0410E+01,1.7927E+01,1.6549E+01,
     &1.5308E+01,1.4343E+01,1.5310E+01,1.7794E+01,1.9451E+01,2.1108E+01,
     &2.1661E+01,2.1386E+01,1.8490E+01,1.6144E+01,1.3386E+01,1.1041E+01,
     &9.3860E+00,8.4219E+00,8.8376E+00,7.8738E+00,6.4965E+00,4.7080E+00,
     &3.8869E+00,3.3456E+00,2.6682E+00,2.5409E+00,2.6896E+00,2.6974E+00/


C  initialize cross section tables

      if(init) then
        N = 36
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGTOT_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_kmp = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGTOT_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_kmp = FV(1) 

      END


      DOUBLE PRECISION FUNCTION sigtot_kmp(plab)
C-----------------------------------------------------------------------
C
C     low-energy K-p total cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      implicit double precision (A-H,O-Z)
      save

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      data init /.true./

C  K-p total cross section
      DATA (PTPP(K),K=    1,   43) /
     &  -1.3500E+00,-1.2345E+00,-9.8216E-01,-8.2491E-01,-7.4143E-01,
     &  -6.1508E-01,-4.5679E-01,-3.7223E-01,-2.9802E-01,-2.6595E-01,
     &  -1.7037E-01,-1.0660E-01,-2.1599E-02,-2.5037E-04, 6.3445E-02,
     &   8.4428E-02, 1.3703E-01, 1.5769E-01, 1.8898E-01, 2.4156E-01,
     &   3.3667E-01, 3.5796E-01, 4.1106E-01, 5.1700E-01, 5.9099E-01,
     &   6.5431E-01, 6.9651E-01, 7.7067E-01, 8.5538E-01, 9.6104E-01,
     &   1.1303E+00, 1.3209E+00, 1.4266E+00, 1.5853E+00, 1.8075E+00,
     &   1.9769E+00, 2.4743E+00, 3.0353E+00, 3.5222E+00, 4.0515E+00,
     &   4.6550E+00, 5.1949E+00, 5.7455E+00/
      DATA (STPP(K),K=    1,   43) /
     &9.7669E+01,8.8840E+01,7.2700E+01,5.8076E+01,4.6625E+01,4.0142E+01,
     &3.5315E+01,3.4074E+01,3.5041E+01,3.7939E+01,4.0838E+01,4.3185E+01,
     &4.6084E+01,4.7740E+01,4.9397E+01,4.7603E+01,4.4430E+01,3.9601E+01,
     &3.5186E+01,3.1876E+01,3.0221E+01,3.1325E+01,3.2982E+01,3.3674E+01,
     &3.2571E+01,3.0640E+01,2.9261E+01,2.9814E+01,2.9953E+01,2.8023E+01,
     &2.6922E+01,2.6924E+01,2.5684E+01,2.4859E+01,2.4034E+01,2.3761E+01,
     &2.2112E+01,2.1155E+01,2.0472E+01,2.0480E+01,2.0627E+01,2.0773E+01,
     &2.1472E+01/


C  initialize cross section tables

      if(init) then
        N = 43
        M = 0
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      'SIGTOT_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_kmp = 0.
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        call SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      'SIGTOT_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_kmp = FV(1) 

      END




      SUBROUTINE SPLIN3(X,Y,DERIV,N,NC,Z,FVALUE,FDERIV,M,MC,IOP)
C***********************************************************************
C
C     CERN LIBRARY PROGRAM NO E-209.
C
C     REVISED VERSION JULY 1973.
C
C     CHANGED BY R.ENGEL (10/10/93) TO CONFORM WITH F77 STANDARD
C
C     PURPOSE = TO COMPUTE A NATURAL SPLINE APPROXIMATION OF THIRD ORDER
C               FOR A FUNCTION Y(X) GIVEN IN THE N POINTS (X(I),Y(I)) ,
C               I=1(1)N.
C
C     PARAMETERS (IN LIST).
C
C     X       = AN ARRAY STORING THE INPUT ARGUMENTS.DIMENSION X(N).
C     Y       = AN ARRAY STORING THE INPUT FUNCTION VALUES.THE ELEMENT
C               Y(I) REPRESENT THE FUNCTION VALUE Y(X) FOR X=X(I).
C     DERIV   = AN ARRAY USED FOR STORING THE COMPUTED DERIVATIVES OF
C               THE FUNCTION Y(X).IN DERIV(I,1) AND DERIV(I,2) ARE STOR-
C               ED THE FIRST-AND SECOND ORDER DERIVATIVES OF Y(X) FOR
C               X=X(I) RESPECTIVELY.
C     N       = NUMBER OF INPUT FUNCTION VALUES.
C     NC      = ARRAY DERIV IS DIMENSIONED DERIV(NC,2) IN CALLING
C               PROGRAM.
C     Z       = AN ARRAY STORING THE ARGUMENTS FOR THE INTERPOLATED
C               VALUES TO BE COMPUTED.
C     FVALUE  = AN ARRAY STORING THE COMPUTED INTERPOLATED VALUES.
C               FVALUE(J) REPRESENT THE FUNCTION VALUE FVALUE(Z) FOR
C               Z=Z(J).
C     FDERIV    = AN ARRAY USED FOR STORING THE DERIVATIVES OF THE COM-
C               PUTED INTERPOLATED VALUES.EXPLANATION AS FOR DERIV.
C     M       = NUMBER OF INTERPOLATED VALUES TO BE COMPUTED.
C     MC      = ARRAY FDERIV IS DIMENSIONED FDERIV(MC,2) IN CALLING
C               PROGRAM.
C     IOP     = OPTION PARAMETER.FOR IOP.LE.0 THE DERIVATIVES FOR EACH
C               SUB-INTERVAL IN THE SPLINE APPROXIMATION ARE COMPUTED.
C                                  IOP=-1, THE SECOND ORDER END-POINT
C                                          DERIVATIVES ARE COMPUTED BY
C                                          LINEAR EXTRAPOLATION.
C                                  IOP=0 , THE SECOND ORDER END-POINT
C                                          DERIVATIVES ASSUMED TO BE GI-
C                                          VEN (SEE COMMON /SPAPPR/).
C                                  IOP=1 , COMPUTE SPLINE APPROXIMATIONS
C                                          FOR THE ARGUMENTS GIVEN IN
C                                          THE ARRAY Z,THE DERIVATIVES
C                                          BEEING ASSUMED TO HAVE BEEN
C                                          CALCULATED IN A PREVIOUS CALL
C                                          ON THE ROUTINE.
C
C     PARAMETERS (IN COMMON BLOCK / SPAPPR /).
C
C     SECD1   = VALUE OF THE SECOND DERIVATIVE D2Y(X)/DX2 FOR THE INPUT
C               ARGUMENT X=X(1).
C     SECDN   = VALUE OF THE SECOND DERIVATIVE D2Y(X)/DX2 FOR THE INPUT
C               ARGUMENT X=X(N).
C               NB. VALUES HAVE TO BE ASSIGNED TO SECD1 AND SECDN IN THE
C               CALLING PROGRAM.IF A NATURAL SPLINE FIT IS WANTED PUT
C               SECD1=SECDN=0.
C     VOFINT  = COMPUTED APPROXIMATION FOR THE INTEGRAL OF Y(X) TAKEN
C               FROM X(1) TO X(N).
C     IERR    = ERROR PARAMETER.IERR=0,NO ERRORS OCCURED.
C                               IERR=1,THE NUMBER OF POINTS TOO SMALL
C                                      I.E.N LESS THAN 4.
C                               IERR=2,THE ARGUMENTS X(I) NOT IN INCREA-
C                                      SING ORDER.
C                               IERR=3,ARGUMENT TO BE USED IN INTERPOLA-
C                                      TION ABOVE RANGE.
C                               IERR=4,ARGUMENT TO BE USED IN INTERPOLA-
C                                      TION BELOW RANGE.
C     NXY     = N (SEE ABOVE),HAS TO BE STORED FOR ENTRIES CORRESPONDING
C               TO IOP=1.
C
C**********************************************************************
      implicit double precision (A-H,O-Z)
      save
C
      DIMENSION X(NC) , Y(NC) , DERIV(NC,2) , Z(MC) , FVALUE(MC) ,
     1          FDERIV(MC,2)
C
      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY
      DATA ZERO,HALF,ONE,THREE/0.,.5,1.,3./
      DATA THIRD , SIXTH / .333333333333333 , .166666666666667 /
C
 1000 IF (IOP.GT.0) GO TO 1110
C
      IERR=0
C
C     CHECK IF ENOUGH DATA-POINTS ARE AVAILABLEI.E. IF N LESS THAN 4 NO
C     THIRD ORDER SPLINE APPROXIMATION IS POSSIBLE.
C
      IF (N.GE.4) GO TO 1010
C
      IERR=1
      GO TO 2000
C
C     START CALCULATION OF COEFFICIENTS TO BE USED IN THE SYSTEM OF EQU-
C     ATIONS FOR THE SECOND ORDER DERIVATIVES OF Y(X).
C
 1010 IF (IOP.NE.-1) GO TO 1015
      SECD1=ZERO
      SECDN = ZERO
      BET1=ONE/(ONE+HALF*(X(2)-X(1))/(X(3)-X(2)))
      ALF1=BET1*(ONE- ((X(2)-X(1))/(X(3)-X(2)))**2)
      BETN=ONE/(ONE+HALF*(X(N)-X(N-1))/(X(N-1)-X(N-2)))
      ALFN=BETN*(ONE- ((X(N)-X(N-1))/(X(N-1)-X(N-2)))**2)
C
 1015 DERIV(1,2)=SECD1
      DERIV(N,2)=SECDN
      DERIV(1,1)=ZERO
      DXPLUS=X(2)-X(1)
C
C     CHECK IF ARGUMENTS ARE IN INCREASING ORDER.IF NOT PRINT ERROR
C     MESSAGE AND STOP.
C
      IF ( DXPLUS.GT.ZERO) GO TO 1020
      IN=1
      IERR=2
      GO TO 2000
C
 1020 DYPLUS=(Y(2)-Y(1))/DXPLUS
      IU=N-1
      DO 1040 I=2,IU
      DXMIN =DXPLUS
      DYMIN =DYPLUS
      DXPLUS=X(I+1)-X(I)
C
C     CHECK IF ARGUMENTS ARE IN INCREASING ORDER.IF NOT PRINT ERROR
C     MESSAGE AND STOP.
C
      IF (DXPLUS.GT.ZERO) GO TO 1030
C
      IN=I
      IERR=2
      GO TO 2000
C
 1030 DXINV =ONE/(DXPLUS+DXMIN)
      DYPLUS=(Y(I+1)-Y(I))/DXPLUS
      DIVDIF=DXINV*(DYPLUS-DYMIN)
      ALF   =HALF*DXINV*DXMIN
      BET   =HALF-ALF
C
      IF (I.EQ.2)  DIVDIF=DIVDIF-THIRD*ALF*DERIV(1,2)
      IF (I.EQ.IU) DIVDIF=DIVDIF-THIRD*BET*DERIV(N,2)
      IF (I.EQ.2) ALF=ZERO
C
      IF (IOP.NE.-1) GO TO 1035
      IF (I.NE.2) GO TO 1032
      BET=BET*ALF1
      DIVDIF=DIVDIF*BET1
      GO TO 1035
 1032 IF (I.NE.IU) GO TO 1035
      ALF=ALF*ALFN
      DIVDIF=DIVDIF*BETN
C
 1035 DXINV =ONE/(ONE+ALF*DERIV(I-1,1))
      DERIV(I,1)=-DXINV*BET
      DERIV(I,2)= DXINV*(THREE*DIVDIF-ALF*DERIV(I-1,2))
 1040 CONTINUE
C
C     COMPUTE THE SECOND DERIVATIVES BY BACKWARDS RECURRENCE RELATION.
C     THE SECOND ORDER DERIVATIVES FOR X=X(N-1) ALREADY COMPUTED.
C
 1050 DO 1060 I=2,IU
      J=N-I
      DERIV(J,2)=DERIV(J,1)*DERIV(J+1,2)+DERIV(J,2)
 1060 CONTINUE
C
      IF (IOP.NE.-1) GO TO 1070
      DERIV(1,2)=((X(3)-X(1))/(X(3)-X(2)))*DERIV(2,2)-((X(2)-X(1))/(X(3)
     1-X(2)))*DERIV(3,2)
      DERIV(N,2)=-((X(N)-X(N-1))/(X(N-1)-X(N-2)))*DERIV(N-2,2)+((X(N)-X(
     1N-2))/(X(N-1)-X(N-2)))*DERIV(N-1,2)
C
C     CALCULATION OF THE SECOND ORDER DERIVATIVES FINISHED.START CAL-
C     CULATION OF THE FIRST ORDER DERIVATIVES AND OF THE INTEGRAL.
C
 1070 VOFINT=ZERO
      DO 1080 I=1,IU
      DXPLUS=X(I+1)-X(I)
      DYPLUS=Y(I+1)-Y(I)
      DIVDIF=DYPLUS/DXPLUS
      DERIV(I,1)=DIVDIF-DXPLUS*(THIRD*DERIV(I,2)+SIXTH*DERIV(I+1,2))
      DXPLUS=HALF*DXPLUS
      VOFINT=VOFINT+DXPLUS*(Y(I+1)+Y(I)-THIRD*(DERIV(I+1,2)+DERIV(I,2))*
     1DXPLUS**2)
 1080 CONTINUE
C
C     COMPUTE THE LAST FIRST ORDER DERIVATIVE.
C
      DXPLUS=X(N)-X(N-1)
      DYPLUS=Y(N)-Y(N-1)
      DIVDIF=DYPLUS/DXPLUS
      DERIV(N,1)=DIVDIF+DXPLUS*(SIXTH*DERIV(N-1,2)+THIRD*DERIV(N,2))
C
C     CALCULATION OF FIRST ORDER DERIVATIVES AND INTEGRAL FINISHED.
C
C     SET VALUE OF N IN COMMON BLOCK / SPAPPR /.
C
      NXY=N
C
C     COMPUTE INTERPOLATED VALUES IF ANY.
C
 1110 IF (M.LT.1) RETURN
C
      XL=X(1)
      XU=X(2)
      IP=3
      IL=0
C
 1120 DO 1160 J=1,M
      ARG=Z(J)
      IF (ARG.GT.XU) GO TO 1170
      IF (ARG.LT.XL) GO TO 1190
C
C     ARGUMENT IN CORRECT RANGE.CHECK IF POLYNOMIAL COEFFICIENTS HAVE
C     TO BE CALCULATED.
C
 1130 IF (IL.GT.0) GO TO 1150
C
C     COMPUTE POLYNOMIAL COEFFICIENTS.
C
 1140 II=IP-2
      A0=Y(II)
      A1=DERIV(II,1)
      A4=DERIV(II,2)
      A6=(DERIV(II+1,2)-A4)/(XU-XL)
      A2=HALF*A4
      A3=SIXTH*A6
      A5=HALF*A6
      IL=1
C
C     CALCULATION OF POLYNOMIAL COEFFICIENTS FINISHED.COMPUTE VALUES.
C
 1150 ARG=ARG-XL
      FVALUE(J)=((A3*ARG+A2)*ARG+A1)*ARG+A0
      FDERIV(J,1)=(A5*ARG+A4)*ARG+A1
      FDERIV(J,2)=A6*ARG+A4
C
 1155 CONTINUE
      GOTO 1160
C
C     RANGE MOVING
C
C
C     ARGUMENT ABOVE PRESENT RANGE.SHIFT RANGE UPWARDS.
C
 1170 IF(IP.GT.NXY) GO TO 1185
      IPP=IP
      DO 1180 I=IPP,NXY
      IF (ARG.GT.X(I)) GO TO 1180
      XL=X(I-1)
      XU=X(I)
      IP=I+1
      IL=0
      GO TO 1140
C
 1180 CONTINUE
C
C     ARGUMENT  OUT OF RANGE,I.E. ARG GREATER THAN X(N).
C
 1185 IERR=3
      IP=NXY+1
      GO TO 2010
C
C     ARGUMENT BELOW PRESENT RANGE.SHIFT DOWNWARDS.
C
 1190 IPP=IP
      DO 1200 I=1,IPP
      II=IP-I-2
      IF (II.EQ.0) GO TO 1210
      IF (ARG.LT.X(II)) GO TO 1200
      XL=X(II)
      XU=X(II+1)
      IP=II+2
      IL=0
      GO TO 1140
C
 1200 CONTINUE
C
C     ARGUMENT OUT OF RANGE,I.E. ARG LESS THAN X(1).
C
 1210 IERR=4
      IP=3
      GO TO 2010
C
 2010 WRITE(6,3000)  IERR , ARG
C
      FVALUE(J)=ZERO
      FDERIV(J,1)=ZERO
      FDERIV(J,2)=ZERO
C
      II=IP-2
      XL=X(II)
      XU=X(II+1)
      IL=0
      GO TO 1155
C
C
C     END OF INTERPOLATION LOOP
C
 1160 CONTINUE
C
C     CALCULATION OF INTERPOLATED VALUES FINISHED.
C
      RETURN
C
C     PRINT ERROR MESSAGES.
C
 2000 IF (IERR.EQ.1) WRITE(6,3000)  IERR
      IF (IERR.EQ.2) WRITE(6,3000)  IERR , X(IN) , X(IN+1)
      RETURN
C
 3000 FORMAT(//5X,'*** SUBROUTINE SPLIN3 ERROR NO ',I2,' ***',
     1       2(4X,E21.14))
C
      END
