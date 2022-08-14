    Double Precision Function xshn(ecm, na, nb, ijproj, ntarg)
      ! **********************************************************************
      !*
      ! Inelastic cross sections according to *
      !
      ! Glauber's approach.                                                  *
      !
      ! NA / NB     mass numbers of proj./target nuclei                     *
      ! IJPROJ      bamjet-index of projectile (=1 in case of proj.nucleus) *
      ! ECMI kinematical variables   E_cm                      *
      ! IE       indices of energy
      ! NTARG       index of target nucleus set o NTARG=1 here             *
      ! This version dated 17.3.98  is written by S. Roesler mod by J.R.    *
      ! **********************************************************************
      !*
      Implicit Double Precision (A-H, O-Z)
      Save
      Parameter (lout=6, llook=9)
      Complex *16 czero, cone, ctwo
      Parameter (zero=0.0D0, one=1.0D0, two=2.0D0, three=3.0D0, &
        onethi=one/three, tiny25=1.0D-25)
      Parameter (twopi=6.283185307179586454D+00, pi=twopi/two, &
        gev2mb=0.38938D0, gev2fm=0.1972D0, alphem=one/137.0D0, & ! proton mass
        amp=0.938D0, amp2=amp**2, & ! approx. nucleon radius
        rnucle=1.12D0, &           ! number of bins in b-space
        ksiteb=200)
      Character *8 aname
      Common /dpar/aname(210), aam(210), gam(210), tau(210), iich(210), &
        iibar(210), ka1(210), ka2(210)
      Parameter (ncompx=1, neb=50)
      Common /dshmm/rash, rbsh(ncompx), bmax(ncompx), bstep(ncompx), sigsh, &
        rosh, gsh, bsite(0:neb, ncompx, ksiteb), nstatb, nsiteb
      Common /glaber/ecmnn(neb), ecmnow, xstot(neb), xsela(neb), xsqep(neb), &
        xsqet(neb), xsqe2(neb), xspro(neb), xetot(neb), xeela(neb), &
        xeqep(neb), xeqet(neb), xeqe2(neb), xepro(neb), bslope, elabb(neb)
      Common /vdmpar/rl2, epspol, intrge(2), idpdf, modega, ishad(3)
      Common /glapar/jstatb
      Complex *16 c, ca, ci
      Common /damp/ca, ci, ga
      Common /xsecnu/ecmuu, ecmoo, ngritt, nevtt
      Common /kglaub/jglaub
      Parameter (maxncl=210)
      Complex *16 pp11, pp12, pp21, pp22, ompp11, ompp12, ompp21, ompp22
      Dimension coop1(3, maxncl), coot1(3, maxncl), coop2(3, maxncl), &
        coot2(3, maxncl), bprod(ksiteb), sigshh(neb), sigto(neb), sigel(neb), &
        sigin(neb), sigsd(neb), sigdif(neb)

      jglaub = 1
      czero = dcmplx(zero, zero)
      cone = dcmplx(one, zero)
      ctwo = dcmplx(two, zero)
      ! re-define kinematics
      ie = 1
      ec111 = ecm
      s = ec111**2
      ecmnn(ie) = ec111
      Write (6, *) 'IE,EC111,S', ie, ec111, s
      ! parameters determining statistics in evaluating Glauber-xsection
      jstatb = nevtt
      nstatb = jstatb
      nsiteb = ksiteb
      ! set up interaction geometry (common /DSHM/)
      ! projectile/target radii
      rash = rnucle*dble(na)**onethi
      rbsh(ntarg) = rnucle*dble(nb)**onethi
      If (jglaub==1) Then
        If (na==9) rash = 2.52D0
        If (na==10) rash = 2.45D0
        If (na==11) rash = 2.37D0
        If (na==12) rash = 2.45D0
        If (na==13) rash = 2.44D0
        If (na==14) rash = 2.55D0
        If (na==15) rash = 2.58D0
        If (na==16) rash = 2.71D0
        If (na==17) rash = 2.66D0
        If (na==18) rash = 2.71D0
        If (nb==9) rbsh(ntarg) = 2.52D0
        If (nb==10) rbsh(ntarg) = 2.45D0
        If (nb==11) rbsh(ntarg) = 2.37D0
        If (nb==12) rbsh(ntarg) = 2.45D0
        If (nb==13) rbsh(ntarg) = 2.44D0
        If (nb==14) rbsh(ntarg) = 2.55D0
        If (nb==15) rbsh(ntarg) = 2.58D0
        If (nb==16) rbsh(ntarg) = 2.71D0
        If (nb==17) rbsh(ntarg) = 2.66D0
        If (nb==18) rbsh(ntarg) = 2.71D0
      End If
      ! maximum impact-parameter
      bmax(ntarg) = 4.0D0*(rash+rbsh(ntarg))
      bstep(ntarg) = bmax(ntarg)/dble(nsiteb-1)
      ! slope, rho ( Re(f(0))/Im(f(0)) )
      If (ijproj<=12) Then
        bslope = 8.5D0*(1.0D0+0.065D0*log(s))
        If (ecmnn(ie)<=3.0D0) Then
          rosh = -0.43D0
        Else If ((ecmnn(ie)>3.0D0) .And. (ecmnn(ie)<=50.D0)) Then
          rosh = -0.63D0 + 0.175D0*log(ecmnn(ie))
        Else If (ecmnn(ie)>50.0D0) Then
          rosh = 0.1D0
        End If
      Else
        bslope = 6.0D0*(1.0D0+0.065D0*log(s))
        rosh = 0.01D0
      End If
      ! projectile-nucleon xsection (in fm)
      elab = (s-aam(ijproj)**2-amp2)/(two*amp)
      elabb(ie) = elab/1000.
      plab = sqrt((elab-aam(ijproj))*(elab+aam(ijproj)))
      ! SIGSH = SHNTOT(IJPROJ,1,ZERO,PLAB)/10.0D0
      sigsh = dshpto(ijproj, plab)/10.D0
      sigshh(ie) = sigsh*10.D0
      Write (6, *) ' NSTATB,NSITEB,RASH,RBSH(NTARG),BMAX(NTARG), &
        &          BSLOPE,ROSH,SIGSH,ECM ELAB', nstatb, nsiteb, rash, &
        rbsh(ntarg), bmax(ntarg), bslope, rosh, sigsh, ec111, elab
      ! initializations
      Do i = 1, nsiteb
        bsite(0, ntarg, i) = zero
        bsite(ie, ntarg, i) = zero
        bprod(i) = zero
      End Do
      stot = zero
      stot2 = zero
      sela = zero
      sela2 = zero
      sqep = zero
      sqep2 = zero
      sqet = zero
      sqet2 = zero
      sqe2 = zero
      sqe22 = zero
      spro = zero
      spro2 = zero
      facn = one/dble(nstatb)
      ipnt = 0
      rpnt = zero
      ! ------------------------------------------------------
      ! cross sections averaged over NSTATB nucleon configurations
      Do is = 1, nstatb
        stotn = zero
        selan = zero
        sqepn = zero
        sqetn = zero
        sqe2n = zero
        spron = zero
        Call conuclx(coop1, na, rash, 0)
        Call conuclx(coot1, nb, rbsh(ntarg), 1)
        Call conuclx(coop2, na, rash, 0)
        Call conuclx(coot2, nb, rbsh(ntarg), 1)
        ! integration over impact parameter B
        Do ib = 1, nsiteb - 1
          stotb = zero
          selab = zero
          sqepb = zero
          sqetb = zero
          sqe2b = zero
          sprob = zero
          sdir = zero
          b = dble(ib)*bstep(ntarg)
          facb = 10.0D0*twopi*b*bstep(ntarg)
          ! integration over M_V^2 for photon-proj.
          ! DO 14 IM=1,JPOINT
          pp11 = cone
          pp12 = cone
          pp21 = cone
          pp22 = cone
          shi = zero
          facm = one
          dcoh = 1.0D10
          ! ------------------------------------------------------------
          gsh = 10.0D0/(two*bslope*gev2mb)
          ! common /DAMP/
          ga = gsh
          rca = ga*sigsh/twopi
          fca = -rosh*rca
          ca = dcmplx(rca, fca)
          ci = cone
          Do ina = 1, na
            kk1 = 1
            kk2 = 1
            Do inb = 1, nb
              x11 = b + coot1(1, inb) - coop1(1, ina)
              y11 = coot1(2, inb) - coop1(2, ina)
              xy11 = ga*(x11*x11+y11*y11)
              x12 = b + coot2(1, inb) - coop1(1, ina)
              y12 = coot2(2, inb) - coop1(2, ina)
              xy12 = ga*(x12*x12+y12*y12)
              x21 = b + coot1(1, inb) - coop2(1, ina)
              y21 = coot1(2, inb) - coop2(2, ina)
              xy21 = ga*(x21*x21+y21*y21)
              x22 = b + coot2(1, inb) - coop2(1, ina)
              y22 = coot2(2, inb) - coop2(2, ina)
              xy22 = ga*(x22*x22+y22*y22)
              If (xy11<=15.0D0) Then
                c = cone - ca*exp(-xy11)
                ar = dble(pp11)
                ai = dimag(pp11)
                If (abs(ar)<tiny25) ar = zero
                If (abs(ai)<tiny25) ai = zero
                pp11 = dcmplx(ar, ai)
                pp11 = pp11*c
                ar = dble(c)
                ai = dimag(c)
                shi = shi + log(ar*ar+ai*ai)
              End If
              If (xy12<=15.0D0) Then
                c = cone - ca*exp(-xy12)
                ar = dble(pp12)
                ai = dimag(pp12)
                If (abs(ar)<tiny25) ar = zero
                If (abs(ai)<tiny25) ai = zero
                pp12 = dcmplx(ar, ai)
                pp12 = pp12*c
              End If
              If (xy21<=15.0D0) Then
                c = cone - ca*exp(-xy21)
                ar = dble(pp21)
                ai = dimag(pp21)
                If (abs(ar)<tiny25) ar = zero
                If (abs(ai)<tiny25) ai = zero
                pp21 = dcmplx(ar, ai)
                pp21 = pp21*c
              End If
              If (xy22<=15.0D0) Then
                c = cone - ca*exp(-xy22)
                ar = dble(pp22)
                ai = dimag(pp22)
                If (abs(ar)<tiny25) ar = zero
                If (abs(ai)<tiny25) ai = zero
                pp22 = dcmplx(ar, ai)
                pp22 = pp22*c
              End If
            End Do
          End Do
          ompp11 = czero
          ompp21 = czero
          ompp11 = ompp11 + (cone-pp11)
          ompp21 = ompp21 + (cone-pp21)
          ompp12 = czero
          ompp22 = czero
          ompp12 = ompp12 + (cone-pp12)
          ompp22 = ompp22 + (cone-pp22)
          stotm = dble(ompp11+ompp22)
          selam = dble(ompp11*dconjg(ompp22))
          sprom = one - exp(shi)
          sqepm = dble(ompp11*dconjg(ompp21)) - selam
          sqetm = dble(ompp11*dconjg(ompp12)) - selam
          sqe2m = dble(ompp11*dconjg(ompp11)) - selam - sqepm - sqetm
          stotb = stotb + facm*stotm
          selab = selab + facm*selam
          If (nb>1) sqepb = sqepb + facm*sqepm
          If (na>1) sqetb = sqetb + facm*sqetm
          If ((na>1) .And. (nb>1)) sqe2b = sqe2b + facm*sqe2m
          sprob = sprob + facm*sprom
          ! 14       CONTINUE
          stotn = stotn + facb*stotb
          selan = selan + facb*selab
          sqepn = sqepn + facb*sqepb
          sqetn = sqetn + facb*sqetb
          sqe2n = sqe2n + facb*sqe2b
          spron = spron + facb*sprob
          bprod(ib+1) = bprod(ib+1) + facn*facb*sprob
        End Do
        stot = stot + facn*stotn
        stot2 = stot2 + facn*stotn**2
        sela = sela + facn*selan
        sela2 = sela2 + facn*selan**2
        sqep = sqep + facn*sqepn
        sqep2 = sqep2 + facn*sqepn**2
        sqet = sqet + facn*sqetn
        sqet2 = sqet2 + facn*sqetn**2
        sqe2 = sqe2 + facn*sqe2n
        sqe22 = sqe22 + facn*sqe2n**2
        spro = spro + facn*spron
        spro2 = spro2 + facn*spron**2
      End Do
      ! final cross sections
      ! 1) total
      xstot(ie) = stot
      ! 2) elastic
      xsela(ie) = sela
      ! 3) quasi-el.: A+B-->A+X (excluding 2)
      xsqep(ie) = sqep
      ! 4) quasi-el.: A+B-->X+B (excluding 2)
      xsqet(ie) = sqet
      ! 5) quasi-el.: A+B-->X (excluding 2-4)
      xsqe2(ie) = sqe2
      ! 6) production (= STOT-SELA-SQEP-SQET-SQE2!)
      xspro(ie) = spro
      xshn = spro
      Write(6,*) 'stot,spro,sela', stot,spro,sela

    End Function xshn

    Subroutine dpmjet2_event(epn, mp, mpch, mt, mtch, ijprojbam, ievframe)

    ! ----------------------------------------------------------------------
    ! dpmjet2_event

    ! Generates one DPMJET II.5 hadron-Air event
    ! Arguments:
    !     epn = energy per nucleon [GeV]
    !     
    !     ijprojbam = BAMJet code of hadron
    !     
    ! ----------------------------------------------------------------------

    Implicit Double Precision (A-H, O-Z)


    Parameter (intmd=252)
    Parameter (intmx=2488)
    Parameter (maxpro=8)
    Parameter (mxafbk=16)
    Parameter (mxffbk=6)
    Parameter (mxnfbk=10)
    Parameter (mxpsfb=41000)
    Parameter (mxpsst=300)
    Parameter (mxzfbk=9)
    Parameter (nmxhkk=49998)
    Parameter (nxafbk=mxafbk+1)
    Parameter (nxnfbk=12)
    Parameter (nxzfbk=11)
    Dimension dsig1(0:maxpro)
    Common /cmhico/cmhis
    Common /collap/s3, ijproj1, ijtar1, ptthr1, ptthr3, iophrd1, ijprlu1, &
      ijtalu1
    Common /collis/ss, ijprox, ijtar, ptthr, ptthr2, iophrd, ijprlu, ijtalu
    Common /diqi/ipvq(248), ippv1(248), ippv2(248), itvq(248), ittv1(248), &
      ittv2(248), ipsq(intmx), ipsq2(intmx), ipsaq(intmx), ipsaq2(intmx), &
      itsq(intmx), itsq2(intmx), itsaq(intmx), itsaq2(intmx), kkproj(248), &
      kktarg(248)
    Character *8 aname
    Common /dpar/aname(210), aam(210), ga(210), tau(210), iich(210), &
      iibar(210), k1(210), k2(210)
    Common /dprin/ipri, ipev, ippa, ipco, init, iphkk, itopd, ipaupr
    Common /evappp/ievap
    Common /extevt/idres(nmxhkk), idxres(nmxhkk), nobam(nmxhkk), &
      idbam(nmxhkk), idch(nmxhkk), npoint(10)
    Common /final/ifinal
    Logical lfrmbk, lncmss
    Common /frbkcm/amufbk, eexfbk(mxpsst), amfrbk(mxpsst), exfrbk(mxpsfb), &
      sdmfbk(mxpsfb), coufbk(mxpsfb), exmxfb, r0frbk, r0cfbk, c1cfbk, &
      c2cfbk, ifrbkn(mxpsst), ifrbkz(mxpsst), ifbksp(mxpsst), &
      ifbkpr(mxpsst), ifbkst(mxpsst), ipsind(0:mxnfbk, 0:mxzfbk, 2), &
      jpsind(0:mxafbk), ifbind(0:nxnfbk, 0:nxzfbk, 2), jfbind(0:nxafbk), &
      ifbcha(5, mxpsfb), iposst, iposfb, ifbstf, ifbfrb, nbufbk, lfrmbk, &
      lncmss
    Common /hkkevt/nhkk, nevhkk, isthkk(nmxhkk), idhkk(nmxhkk), &
      jmohkk(2, nmxhkk), jdahkk(2, nmxhkk), phkk(5, nmxhkk), &
      vhkk(4, nmxhkk), whkk(4, nmxhkk)
    Common /ifroto/ifrovp(248), itovp(248), ifrosp(intmx), ifrovt(248), &
      itovt(248), ifrost(intmx), jsshs(intmx), jtshs(intmx), jhkknp(248), &
      jhkknt(248), jhkkpv(intmx), jhkkps(intmx), jhkktv(intmx), &
      jhkkts(intmx), mhkkvv(intmx), mhkkss(intmx), mhkkvs(intmx), &
      mhkksv(intmx), mhkkhh(intmx), mhkkdv(248), mhkkvd(248), mhkkds(intmd), &
      mhkksd(intmd)
    Common /inpflg/iang, ifiss, ib0, igeom, istrag, keydk
    Common /nncms/gamcm, bgcm, umo, pcm, eproj, pproj
    Common /nucc/it, itz, ip, ipz, ijproj, ibproj, ijtarg, ibtarg
    Common /nucimp/prmom(5, 248), tamom(5, 248), prmfep, prmfen, tamfep, &
      tamfen, prefep, prefen, taefep, taefen, prepot(210), taepot(210), &
      prebin, taebin, fermod, etacou
    Logical ldiffr, linctv, levprt, lheavy, ldeexg, lgdhpr, lpreex, lhlfix, &
      lprfix, lparwv, lpower, lsngch, llvmod, lschdf
    Common /parevt/dpower, fsprd0, fshpfn, rn1gsc, rn2gsc, ldiffr(39), &
      lpower, linctv, levprt, lheavy, ldeexg, lgdhpr, lpreex, lhlfix, &
      lprfix, lparwv, ilvmod, jlvmod, llvmod, lsngch, lschdf
    Common /pomtyp/ipim, icon, isig, lmax, mmax, nmax, difel, difnu
    Common /pydat1/mstu(200), paru(200), mstj(200), parj(200)
    Common /rptshm/rproj, rtarg, bimpac
    Common /seasu3/seasq
    Common /shmakl/jssh(intmx), jtsh(intmx), inter1(intmx), inter2(intmx)
    Common /sigma/sigsof, bs, zsof, sighar, fill(7)
    Common /strufu/istrum, istrut
    Common /taufo/taufor, ktauge, itauve, incmod
    Character *80 titled
    Character *8 projty, targty
    Common /user1/titled, projty, targty
    Common /user2/cmener, sdfrac, ptlar, istruf, isingx, idubld
    Common /xseadi/xseacu, unon, unom, unosea, cvq, cdq, csea, ssmima, &
      ssmimq, vvmthr
    Common /xsecpt/ptcut, sigs, dsigh

    ! Common to store frame information
    Integer iframe
    Common /dpm2frame/ iframe

    ! DOUBLE PRECISION PPNPN(4)
    Logical debug
    Parameter(debug = .False.)
    Save
    ! ----------------------------------------------------------------------
    If (ievframe < 1 .Or. ievframe > 2) Then
      Write(6,*) "dpmjet2_event(): Wrong iframe parameter", ievframe
      Stop
    End If
    iframe = ievframe

    ipri = 0
    ipev = 0
    ippa = 0
    ipco = -2
    init = 0
    iphkk = 0
    mstu(22) = 0
    mstu(26) = 0

    it = mt
    itz = mtch

    !Simplify and assume air is pure nitrogen.
    If (it==0) Then
      it= 14
      itz = 7
    End If

    kkmat = 1
    amtar = aam(1)
    ijproj = ijprojbam
    If (ijproj==12 .Or. ijproj==19) Then
      ! TRANSFORM K(0)L/S  TO  K(0) OR ANTI-K(0)
      If (RNDM()>=0.5D0) Then
        ijproj = 24
      Else
        ijproj = 25
      End If
    End If
    ijprox = ijproj
    ibproj = iibar(ijproj)
    ip = mp
    ipz = mpch
    isingd = 1
    isingx = 1
    idubld = 0
    sdfrac = 1.D0

    ifinal = 1
    ievap = 0
    taufor = 105.
    ktauge = 0
    levprt = .False.
    ilvmod = 1
    ldeexg = .False.
    lheavy = .False.
    lfrmbk = .False.
    ifiss = 0

    amproj = aam(ijproj)
    ppn = sqrt((epn-amproj)*(epn+amproj))
    pproj = ppn
    eproj = epn
    umo = sqrt(amproj**2+amtar**2+2.D0*amtar*eproj)
    cmener = umo
    ss = umo**2

    If (istrut==1) Then
      ptthr = 2.1D0 + 0.15D0*(log10(cmener/50.))**3
      ptthr2 = ptthr
    Else If (istrut==2) Then
      ptthr = 2.5D0 + 0.12D0*(log10(cmener/50.))**3
      ptthr2 = ptthr
    End If
    gamcm = (eproj+amtar)/umo
    bgcm = pproj/umo
    pcm = gamcm*pproj - bgcm*eproj
    If (debug) Write (mdebug, *) 'DPMJLK: AMPROJ,PPN=', amproj, ppn
    alfa = 1.076D0
    a = 37.8D0
    sigsof = a*ss**(alfa-1.D0)

    isingd = 1
    isingx = 1

    ptcut = ptthr
    If (debug) Write (mdebug, *) 'DPMJLK: BEFORE CSJ1MI'
    Call csj1mi(ptcut, dsig1)
    If (debug) Write (mdebug, *) 'DPMJLK: AFTER CSJ1MI'
    sig1 = dsig1(0)
    dsigh = sig1
    ! INITIALIZE TRANSVERSE MOMENTA FOR SOFT SCATTERING
    Call samppt(0, pt)
    If (debug) Write (mdebug, *) 'DPMJLK: AFTER SAMPPT'

    iit = it
    iitz = itz
    iip = ip
    iipz = ipz
    iiproj = ijproj
    e000 = epn

    ! INITIALIZE HARD SCATTERING
    s3 = ss
    ijproj1 = ijprox
    ijtar1 = ijtar
    ptthr1 = ptthr
    iophrd1 = iophrd
    ijprlu1 = ijprlu
    ijtalu1 = ijtalu
    ptthr3 = ptthr2

    ! ----------------------------------------------------------------------
    If (ipev>=6) Then
      Write (6, '(/A,I5/5X,A)') ' DPMJLK: IP=', ip, &
        ' J,IPVQ(J),IPPV1(J),IPPV2(J),ISTHKK,KKPROJ,PRMOM'
      Do j = 1, ip
        Write (6, '(I4,5I3,5(1P,E11.3))') j, isthkk(j), kkproj(j), ipvq(j), &
          ippv1(j), ippv2(j), (prmom(jj,j), jj=1, 5)
      End Do

      Write (6, '(/A,I5/5X,A)') ' KKEVT : IT=', it, &
        ' J,ITVQ(J),ITTV1(J),ITTV2(J),ISTHKK,KKTARG,TAMOM'
      ihkk = ip
      Do j = 1, it
        ihkk = ihkk + 1
        Write (6, '(I4,5I3,5(1P,E11.3))') j, isthkk(ihkk), kktarg(j), &
          itvq(j), ittv1(j), ittv2(j), (tamom(jj,j), jj=1, 5)
      End Do
    End If

    If (debug) Write (mdebug, *) 'DPMJLK: NOW KKINC IS CALLED WITH IJPROJ=', &
      ijproj, ' IT=', it

    ! NOW THE REAL WORK STARTS
    ! 200  CONTINUE
    If (iip==1) Then
      ! FOR HADRON-NUCLEUS COLLISIONS
      elabt = epn
      iiipro = iiproj
      iiip = iip
      iiipz = iipz
      iiit = iit
      iiitz = iitz
      Call kkinc(epn,iit,iitz,1,iipz,iiproj,kkmat,1, nhkkh1,irej)
    Else
      ! FOR NUCLEUS-NUCLEUS COLLISIONS WITH PROJECTILE FRAGMENTATION
      elabt = epn*0.001D0
      iiipro = iiproj
      iiip = iip
      iiipz = iipz
      iiit = iit
      iiitz = iitz
      Call dpmevt(elabt, iiipro, iiip, iiipz, iiit, iiitz, kkmat, nhkkh1)
    End If

    !Replace K0 and K0BAR by with K0S
    Do j = 0, nhkk
      If (isthkk(j) /= 1) Cycle
      If (idhkk(j)==311 .Or. idhkk(j)==-311) Then
        If (RNDM()>=0.5D0) Then
          idhkk(j) = 310
        Else
          idhkk(j) = 130
        End If
      End If
      If (idbam(j) <= 210) Then
        idch(j) = iich(idbam(j))
      End If
    End Do

    ! STORE SECONDARY PARTICLES TO STACK
!       Call dpmjst

    Return
    End Subroutine dpmjet2_event

    Subroutine dpmjet2_lt2lab

    !***********************************************************************
    ! Lorentz-transformation to lab-system. This subroutine scans DTEVT1   *
    ! for final state particles/fragments defined in nucleon-nucleon-cms   *
    ! and transforms them to the lab.                                      *
    ! This version dated 07.01.96 is written by S. Roesler                 *
    !***********************************************************************

      implicit double precision (a-h,o-z)
      save

! event history
      Parameter (nmxhkk=49998)
      Common /hkkevt/nhkk, nevhkk, isthkk(nmxhkk), idhkk(nmxhkk), &
      jmohkk(2, nmxhkk), jdahkk(2, nmxhkk), phkk(5, nmxhkk), &
      vhkk(4, nmxhkk), whkk(4, nmxhkk)

      do i=0,nhkk
         if ((abs(isthkk(i)).eq.1).or.(isthkk(i).eq.1000).or. &
                                  (isthkk(i).eq.1001)) then
            Write(6,*) "transforming", phkk(3,i),phkk(4,i)
            call ltnuc(phkk(3,i),phkk(4,i),pz,pe,-3)
            phkk(3,i) = pz
            phkk(4,i) = pe
            Write(6,*) "to", phkk(3,i),phkk(4,i)
         endif
      End Do

      End Subroutine dpmjet2_lt2lab

    ! -- Author :    D. HECK IK FZK KARLSRUHE       11/11/1996
    ! =======================================================================

    Subroutine dpjsig(plab, icz, paircs)

      ! ----------------------------------------------------------------------
      ! D(UAL) P(ARTON) J(ET MODEL) SIG(MA)

      ! CALCULATES INELASTIC CROSS-SECTION.
      ! THIS SUBROUTINE IS CALLED FROM BOX2.
      ! ARGUMENTS:
      ! PLAB   =  LABORATORY MOMENTUM (GEV)
      ! ICZ    =  HADRON TYPE: 1 = NUCLEON
      ! 2 = PION
      ! 3 = KAON
      ! >99 = NUCLEUS
      ! ----------------------------------------------------------------------

      Implicit None

      Common /crparpar/curpar, secpar, prmpar, outpar, c, & 
        e00, e00pn, ptot0, ptot0n, thickh, itype, levl

      Double Precision curpar(0:16), secpar(0:16), prmpar(0:16), outpar(0:16), &

        c(50), e00, e00pn, ptot0, ptot0n, thickh
      Integer itype, levl


      Common /crrunpar/fixhei, thick0, hiloecm, hiloelb, sig1i, targ1i, &
        stepfc, & 
        sigmaq, & 
        nrrun, nshow, mpatap, moniin, moniou, mdebug, nucnuc, mtabout, &
        mlongout, iseed1i, & 
        propmod, lstck, & 
        lstck1, lstck2, & 
        ishowno, ishw, nopart, nrecs, nblks, maxprt, ndebdl, n1sttr, mdbase, &

        debdel, debug, fdecay, fegs, firsti, fixinc, fixtar, fix1i, fmuadd, &
        fnkg, fprint, fdbase, fparout, ftabout, flongout, gheish, ghesig, &
        gheisdb, uselow, tmargin & 
        , foutfile, ifinam & 
        , fflatout
      Common /crrunpac/datdir, dsn, dsntab, dsnlong, host, user & 
        , lstdsn, filout
      Double Precision fixhei, thick0, hiloecm, hiloelb, sig1i, targ1i, stepfc

      Double Precision sigmaq(4)

      Integer nrrun, nshow, mpatap, moniin, moniou, mdebug, nucnuc, ishowno, &
        ishw, nopart, nrecs, nblks, maxprt, ndebdl, n1sttr, mdbase, mtabout, &
        mlongout, iseed1i(3)

      Integer propmod
      Integer lstck & 
        , lstck1, lstck2
      Character *132 filout

      Character *255 dsn, dsntab, dsnlong
      Character *132 datdir
      Character *20 host, user

      Character *9 lstdsn
      Logical debdel, debug, fdecay, fegs, firsti, fixinc, fixtar, fix1i, &
        fmuadd, fnkg, fprint, fdbase, fparout, ftabout, flongout, gheish, &
        ghesig, gheisdb, uselow, tmargin & 
        , fflatout

      Logical foutfile
      Integer ifinam


      Common /crsigm/sigma, sigann, sigair, fractn, frctno, sigairs
      Double Precision sigma, sigann, sigair, fractn, frctno, sigairs


      Double Precision plab, plablg, sgfrcn(4), sgfrno(4), sgpicn(4), &
        sgpino(4), sgkcn(4), sgkno(4), sgpair(4), sgpiair(4), sgkair(4), &
        sigpp(5)
      ! *                 ,SIGPIP(5),SIGKP(5)
      Integer i, icz

      Double Precision paircs
!f2py intent(out) paircs
      Save

      ! FITTED VALUES FOR CROSS-SECTION FUNCTION BY D. HECK, DEC. 4, 2000
      ! FOR PROTON NITROGEN
      Data sgfrcn/87.972D0, 35.978D0, -0.61458D0, 101.55D0/
      ! FOR PROTON NITROGEN+OXYGEN
      Data sgfrno/116.04D0, 46.477D0, -0.80115D0, 130.69D0/
      ! FOR PROTON AIR
      Data sgpair/117.73D0, 46.875D0, -0.81046D0, 131.734D0/
      ! FITTED VALUES FOR CROSS-SECTION FUNCTION BY D. HECK, DEC. 4, 2000
      ! FOR PION NITROGEN
      Data sgpicn/44.493D0, 34.177D0, -0.75053D0, 87.036D0/
      ! FOR PION NITROGEN+OXYGEN
      Data sgpino/59.134D0, 44.289D0, -0.97918D0, 112.29D0/
      ! FOR PION AIR
      Data sgpiair/60.117D0, 44.719D0, -0.99175D0, 113.29D0/
      ! FITTED VALUES FOR CROSS-SECTION FUNCTION BY D. HECK, DEC. 4, 2000
      ! FOR KAON NITROGEN
      Data sgkcn/32.754D0, 38.879D0, -1.0726D0, 60.358D0/
      ! FOR KAON NITROGEN+OXYGEN
      Data sgkno/43.931D0, 50.388D0, -1.3971D0, 77.539D0/
      ! FOR KAON AIR
      Data sgkair/44.772D0, 50.878D0, -1.4137D0, 78.189D0/

      ! FOR PROTON-PROTON (04.12.2000)
      Data sigpp/ -257.36D0, 27.545D0, -0.53336D0, 2334.7D0, 7.9204D0/
      ! FOR PION PROTON (04.12.2000)
      ! DATA SIGPIP /-159.05D0, 18.743D0, -0.3839D0, 1282.7D0, 6.9655D0/
      ! FOR KAON PROTON (04.12.2000)
      ! DATA SIGKP /-102.06D0, 15.253D0, -0.29995D0, 726.89D0, 6.1179D0/
      ! ----------------------------------------------------------------------

      If (debug) Write (mdebug, *) 'DPJSIG: PLAB=', sngl(plab), ' ICZ=', icz
      ! DECADIC LOGARITH OF LABORATORY MOMENTUM
      plablg = log10(plab)

      If (icz==1) Then
        ! FOR BARYON PROJECTILES
        sigair = sgpair(1) + sgpair(2)*plablg + sgpair(3)*plablg**2 + &
          sgpair(4)/plablg
        fractn = sgfrcn(1) + sgfrcn(2)*plablg + sgfrcn(3)*plablg**2 + &
          sgfrcn(4)/plablg
        frctno = sgfrno(1) + sgfrno(2)*plablg + sgfrno(3)*plablg**2 + &
          sgfrno(4)/plablg
        sigma = 0.D0

      Else If (icz==2) Then
        ! FOR PION PROJECTILES
        sigair = sgpiair(1) + sgpiair(2)*plablg + sgpiair(3)*plablg**2 + &
          sgpiair(4)/plablg
        fractn = sgpicn(1) + sgpicn(2)*plablg + sgpicn(3)*plablg**2 + &
          sgpicn(4)/plablg
        frctno = sgpino(1) + sgpino(2)*plablg + sgpino(3)*plablg**2 + &
          sgpino(4)/plablg
        sigma = 0.D0

      Else If (icz==3) Then
        ! FOR KAON PROJECTILES
        sigair = sgkair(1) + sgkair(2)*plablg + sgkair(3)*plablg**2 + &
          sgkair(4)/plablg
        fractn = sgkcn(1) + sgkcn(2)*plablg + sgkcn(3)*plablg**2 + &
          sgkcn(4)/plablg
        frctno = sgkno(1) + sgkno(2)*plablg + sgkno(3)*plablg**2 + &
          sgkno(4)/plablg
        sigma = 0.D0

      Else If (icz>=200) Then
        ! FOR NUCLEUS PROJECTILES DETERMINE ONLY NN CROSS-SECTION
        sigma = sigpp(1) + sigpp(2)*plablg + sigpp(3)*plablg**2 + &
          sigpp(4)/(plablg+sigpp(5))
        sigair = 0.D0
        fractn = 0.D0
        frctno = 0.D0

      Else

        Write (moniou, 100)(curpar(i), i=0, 9)
100     Format (' DPJSIG: CURPAR=', 1P, 10E11.3)

        Write (moniou, *) 'DPJSIG: ILLEGAL PROJECTILE TYP =', icz
        Stop
      End If

      If (debug) Write (mdebug, *) 'DPJSIG: SIGMA=', sngl(sigma), ' SIGAIR=', &
        sngl(sigair)
      paircs = sigair
      Return
    End Subroutine dpjsig


    ! -- Author :    D. HECK IK FZK KARLSRUHE       14/02/1996
    ! =======================================================================

    Block Data dpmdat

      ! ----------------------------------------------------------------------
      ! DPM(JET) DAT(A)

      ! SETS PARTICLE CODE TABLES FOR DPMJET-II LINK
      ! ----------------------------------------------------------------------

      Implicit None

      Common /crdpmlin/ictabl
      Integer ictabl(200)

      Common /kainit/umoda, inich1, inich2, inich3, inich4, inich5, inich6, &
        inich7, inich8, inich9, iretur, kais, inipri, nusept, musept, istart
      Double Precision umoda
      Integer inich1, inich2, inich3, inich4, inich5, inich6, inich7, inich8, &
        inich9, iretur, kais, inipri, nusept, musept
      Logical istart

      ! ICTABL CONVERTS CORSIKA PARTICLES INTO DPMJET PARTICLES
      Data ictabl/7, 4, 3, 0, 10, 11, 23, 13, 14, 12, & ! CHARMED PARTICLES
        15, 16, 8, 1, 2, 19, 31, 17, 21, 22, & 
        20, 97, 98, 109, 9, 18, 99, 100, 101, 102, & ! RESET ALL INDICATORS
                                                     ! FOR FIRST
                                                     ! INITIALIZATION OF
                                                     ! COMMON /KAINIT/
        103, 115, 68*0, &          ! 10
        10*0, &                    ! 20
        0, 0, 0, 0, 0, 116, 117, 118, 119, 120, & ! 30
        121, 122, 123, 124, 125, 126, 127, 128, 0, 130, & ! 110
        131, 132, 133, 134, 0, 0, 137, 138, 139, 140, & ! 120
        141, 142, 143, 144, 145, 0, 0, 0, 149, 150, & ! 130
        151, 152, 153, 154, 155, 156, 157, 0, 0, 0, & ! 140
        161, 162, 163, 0, 0, 0, 0, 0, 0, 0, & ! 150
        171, 172, 173, 0, 0, 25*0/
      ! 160
      ! 170
      Data umoda/0.D0/, inich1/0/, inich2/0/, inich3/0/, inich4/0/
      Data inich5/0/, inich6/0/, inich7/0/, inich8/0/, inich9/0/
      Data iretur/0/, kais/0/, inipri/0/, nusept/0/, musept/0/
      Data istart/.True./

    End Block Data dpmdat

    ! -- Author :    D. HECK IK FZK KARLSRUHE       07/03/1996
    ! =======================================================================

    Subroutine dpmjin(iseedin, datadir)

      ! ----------------------------------------------------------------------
      ! DPMJ(ET) IN(ITIALIZE)

      ! INITIALIZES DPMJET-II.
      ! THIS SUBROUTINE IS CALLED FROM START.
      ! ----------------------------------------------------------------------

      Implicit Double Precision (A-H, O-Z)


      Common /crrunpar/fixhei, thick0, hiloecm, hiloelb, sig1i, targ1i, &
        stepfc, & 
        sigmaq, & 
        nrrun, nshow, mpatap, moniin, moniou, mdebug, nucnuc, mtabout, &
        mlongout, iseed1i, & 
        propmod, lstck, & 
        lstck1, lstck2, & 
        ishowno, ishw, nopart, nrecs, nblks, maxprt, ndebdl, n1sttr, mdbase, &
        debdel, debug, fdecay, fegs, firsti, fixinc, fixtar, fix1i, fmuadd, &
        fnkg, fprint, fdbase, fparout, ftabout, flongout, gheish, ghesig, &
        gheisdb, uselow, tmargin & 
        , foutfile, ifinam & 
        , fflatout
      Common /crrunpac/datdir, dsn, dsntab, dsnlong, host, user & 
        , lstdsn, filout
      Double Precision fixhei, thick0, hiloecm, hiloelb, sig1i, targ1i, stepfc

      Double Precision sigmaq(4)

      Integer nrrun, nshow, mpatap, moniin, moniou, mdebug, nucnuc, ishowno, &
        ishw, nopart, nrecs, nblks, maxprt, ndebdl, n1sttr, mdbase, mtabout, &
        mlongout, iseed1i(3)

      Integer propmod
      Integer lstck & 
        , lstck1, lstck2
      Character *132 filout

      Character *255 dsn, dsntab, dsnlong
      Character *132 datdir
      Character *20 host, user

      Character *9 lstdsn
      Logical debdel, debug, fdecay, fegs, firsti, fixinc, fixtar, fix1i, &
        fmuadd, fnkg, fprint, fdbase, fparout, ftabout, flongout, gheish, &
        ghesig, gheisdb, uselow, tmargin & 
        , fflatout

      Logical foutfile
      Integer ifinam


      Common /crdpmjet/levldb, idpmver, fdpmjt, fdpjsg
      Integer levldb, idpmver
      Logical fdpmjt, fdpjsg

      Common /kainit/umoda, inich1, inich2, inich3, inich4, inich5, inich6, &
        inich7, inich8, inich9, iretur, kais, inipri, nusept, musept, istart
      Double Precision umoda
      Integer inich1, inich2, inich3, inich4, inich5, inich6, inich7, inich8, &
        inich9, iretur, kais, inipri, nusept, musept
      Logical istart


      Parameter (intmx=2488)
      Parameter (mxafbk=16)
      Parameter (mxffbk=6)
      Parameter (mxnfbk=10)
      Parameter (mxzfbk=9)
      Parameter (mxpsst=300)
      Parameter (mxpsfb=41000)
      Parameter (nxafbk=mxafbk+1)
      Parameter (nxnfbk=12)
      Parameter (nxzfbk=11)
      Common /casadi/casaxx, icasad
      Common /cmhico/cmhis
      Common /cronin/cronco, mkcron
      Common /colle/nevhad, nvers, ihadrz, nfile
      Common /collis/s, ijprox, ijtar, ptthr, ptthr2, iophrd, ijprlu, ijtalu
      Common /coulo/icoul
      Common /diffra/isingd, idiftp, ioudif, iflagd
      Common /diqsum/ndvuu, ndvus, ndvss, nvduu, nvdus, nvdss, ndsuu, ndsus, &
        ndsss, nsduu, nsdus, nsdss, ndzuu, ndzus, ndzss, nzduu, nzdus, nzdss, &
        nadvuu, nadvus, nadvss, navduu, navdus, navdss, nadsuu, nadsus, &
        nadsss, nasduu, nasdus, nasdss, nadzuu, nadzus, nadzss, nazduu, &
        nazdus, nazdss
      Common /diquax/amedd, idiqua, idiquu
      Common /diqrej/idiqre(7), idvre(3), ivdre(3), idsre(3), isdre(3), &
        idzre(3), izdre(3), idiqrz(7)
      Common /dprin/ipri, ipev, ippa, ipco, init, iphkk, itopd, ipaupr
      Common /dropjj/dropjt, dropva
      Logical intpt, fermp, ihadss, ihadsv, ihadvs, ihadvv, ihada, ipadis, &
        ishmal, lpauli
      Common /droppt/intpt, fermp, ihadss, ihadsv, ihadvs, ihadvv, ihada, &
        ipadis, ishmal, lpauli
      Common /edens/ieden
      Common /evappp/ievap
      Common /final/ifinal
      Common /fluctu/ifluct
      Common /secint/isecin
      Logical lfrmbk, lncmss
      Common /frbkcm/amufbk, eexfbk(mxpsst), amfrbk(mxpsst), exfrbk(mxpsfb), &
        sdmfbk(mxpsfb), coufbk(mxpsfb), exmxfb, r0frbk, r0cfbk, c1cfbk, &
        c2cfbk, ifrbkn(mxpsst), ifrbkz(mxpsst), ifbksp(mxpsst), &
        ifbkpr(mxpsst), ifbkst(mxpsst), ipsind(0:mxnfbk, 0:mxzfbk, 2), &
        jpsind(0:mxafbk), ifbind(0:nxnfbk, 0:nxzfbk, 2), jfbind(0:nxafbk), &
        ifbcha(5, mxpsfb), iposst, iposfb, ifbstf, ifbfrb, nbufbk, lfrmbk, &
        lncmss
      Common /gluspl/nugluu, nsgluu
      Common /hadthr/ehadth, inthad
      Common /hdjase/nhse1, nhse2, nhse3, nhase1, nhase2, nhase3
      Common /hettp/nhstp, nbertp, iosub, insrs
      Common /ifragm/ifrag
      Common /infore/ifrej
      Common /inpflg/iang, ifiss, ib0, igeom, istrag, keydk
      Common /kglaub/jglaub
      Common /nstari/nstart
      Common /ncshxx/ncouxh, ncouxt
      Common /nncms/gamcm, bgcm, umo, pcm, eproj, pproj
      Common /nucc/it, itz, ip, ipz, ijproj, ibproj, ijtarg, ibtarg
      Common /nuccc/jt, jtz, jp, jpz, jjproj, jbproj, jjtarg, jbtarg
      Common /nucimp/prmom(5, 248), tamom(5, 248), prmfep, prmfen, tamfep, &
        tamfen, prefep, prefen, taefep, taefen, prepot(210), taepot(210), &
        prebin, taebin, fermod, etacou
      Common /nuclea/pfermp(2), pfermn(2), fermdd, ebindp(2), ebindn(2), &
        epot(2, 210), etacoo(2), icoull
      Logical ldiffr, linctv, levprt, lheavy, ldeexg, lgdhpr, lpreex, lhlfix, &
        lprfix, lparwv, lpower, lsngch, llvmod, lschdf
      Common /parevt/dpower, fsprd0, fshpfn, rn1gsc, rn2gsc, ldiffr(39), &
        lpower, linctv, levprt, lheavy, ldeexg, lgdhpr, lpreex, lhlfix, &
        lprfix, lparwv, ilvmod, jlvmod, llvmod, lsngch, lschdf
      Common /pomtab/ipomta
      Common /pomtyp/ipim, icon, isig, lmax, mmax, nmax, difel, difnu
      Common /popcor/pdb, ajsdef
      Common /popcck/pdbck, pdbse, pdbseu, ijpock, irejck, ick4, ihad4, ick6, &
        ihad6, irejse, ise4, ise6, irejs3, ise43, ise63, irejs0, ihada4, &
        ihada6, irejsa, isea4, isea6, ireja3, isea43, isea63, irejao
      Common /inxdpm/intdpm
      Common /projk/iprojk
      Common /promu/ipromu
      Common /pshow/ipshow
      Common /ptlarg/xsmax
      Common /ptsamp/isampt
      Common /pydat1/mstu(200), paru(200), mstj(200), parj(200)
      Common /recom/irecom
      Logical lseadi
      Common /seadiq/lseadi
      Common /seaqxx/seaqx, seaqxn
      Common /seasu3/seasq
      Common /sincha/isicha
      Common /stars/istar2, istar3
      Common /strufu/istrum, istrut
      Common /taufo/taufor, ktauge, itauve, incmod
      Character *80 titled
      Character *8 projty, targty
      Common /user1/titled, projty, targty
      Common /user2/cmener, sdfrac, ptlar, istruf, isingx, idubld
      Common /vxsvd/vxsp(50), vxst(50), vxsap(50), vxsat(50), vxvp(50), &
        vxvt(50), vxdp(50), vxdt(50), nxsp, nxst, nxsap, nxsat, nxvp, nxvt, &
        nxdp, nxdt
      Common /xseadi/xseacu, unon, unom, unosea, cvq, cdq, csea, ssmima, &
        ssmimq, vvmthr
      Common /zentra/icentr

      integer          iseedin
      Character *132   datadir 
      Save
      Data ncount/0/
      ! ----------------------------------------------------------------------

      Call INIT_RMMARD(iseedin)
      debug = .False.
      mdebug = 6

      datdir = datadir

      If (debug) Then
        Write (mdebug, *) 'DPMJIN:'
        ! SET PRINT PARAMETERS (DEFAULT SETTING IN BLOCK DATA BLKD41)
        ipri = levldb
        ipev = levldb
        ippa = levldb
        ipco = levldb
        init = levldb
        iphkk = levldb
        mstu(22) = 10
        mstu(26) = 10
      Else
        ! OUTLEVEL
        ipri = 0
        ipev = 0
        ippa = 0
        ipco = -2
        init = 0
        iphkk = 0
        mstu(22) = 0
        mstu(26) = 0
      End If

      idiqre(1) = 0
      idiqre(2) = 0
      idiqre(3) = 0
      idiqre(4) = 0
      idiqre(5) = 0
      idiqre(6) = 0
      idiqre(7) = 0
      idiqrz(1) = 0
      idiqrz(2) = 0
      idiqrz(3) = 0
      idiqrz(4) = 0
      idiqrz(5) = 0
      idiqrz(6) = 0
      idiqrz(7) = 0
      idvre(1) = 0
      idvre(2) = 0
      idvre(3) = 0
      ivdre(1) = 0
      ivdre(2) = 0
      ivdre(3) = 0
      idsre(1) = 0
      idsre(2) = 0
      idsre(3) = 0
      isdre(1) = 0
      isdre(2) = 0
      isdre(3) = 0
      idzre(1) = 0
      idzre(2) = 0
      idzre(3) = 0
      izdre(1) = 0
      izdre(2) = 0
      izdre(3) = 0
      ndvuu = 0
      ndvus = 0
      ndvss = 0
      nvduu = 0
      nvdus = 0
      nvdss = 0
      ndsuu = 0
      ndsus = 0
      ndsss = 0
      nsduu = 0
      nsdus = 0
      nsdss = 0
      ndzuu = 0
      ndzus = 0
      ndzss = 0
      nzduu = 0
      nzdus = 0
      nzdss = 0
      nadvuu = 0
      nadvus = 0
      nadvss = 0
      navduu = 0
      navdus = 0
      navdss = 0
      nadsuu = 0
      nadsus = 0
      nadsss = 0
      nasduu = 0
      nasdus = 0
      nasdss = 0
      nadzuu = 0
      nadzus = 0
      nadzss = 0
      nazduu = 0
      nazdus = 0
      nazdss = 0
      nhse1 = 0
      nhse2 = 0
      nhse3 = 0
      nhase1 = 0
      nhase2 = 0
      nhase3 = 0

      ! 1000 CONTINUE
      ncount = ncount + 1

      ! Parton pt distribution plot initializing
      Call parpt(1, pt1, pt2, ipt, nevt)

      ! Initialization of BAMJET, DECAY, and HADRIN
      Call ddatar
      Call dhadde
      Call dchant
      Call dchanh

      ! Print the title
      Write (moniou, 100)
100   Format ('   ***************************************************', /, &
        '   *         DPMJET VERSION 2.55 (Nov. 2001)         *', /, &
        '   * DUAL PARTON MODEL FOR HADRON NUCLEUS COLLISIONS *', /, &
        '   *         AND NUCLEUS NUCLEUS COLLISIONS          *', /, &
        '   * INCLUDING A FORMATION ZONE INTRANUCLEAR CASCADE *', /, &
        '   *   MINIJETS AND DTUJET LIKE MULTIPLE SOFT JETS   *', /, &
        '   *                                                 *', /, &
        '   * AUTHOR: J. RANFT, email: Johannes.Ranft@cern.ch *', /, &
        '   ***************************************************')
      idpmver = 25
      Call defaul(epn, ppn)
      Call defaux(epn, ppn)

      ! ********************************************************************
      icoul = 1
      icoull = 1
      ! EDENSITY
      ieden = 0
      ! TOPDRAW (option removed)
      itopd = 0
      ! TAUFOR
      ! ---formation zone intranuclear cascade
      ! TAUFOR = 105.   ! for pp interactions
      ! KTAUGE = 0      ! for pp interactions
      taufor = 5.D0
      ktauge = 25
      itauve = 1
      incmod = 1
      ! SEADISTR
      ! ---definition of soft quark distributions
      xseaco = 1.D0
      xseacu = 1.05D0 - xseaco
      unon = 3.5D0
      unom = 1.11D0
      unosea = 5.D0
      ! FERMI
      fermp = .True.
      fermod = 0.6D0
      fermdd = 0.6D0
      iferfo = 1
      ! PAULI
      ipaupr = 0
      lpauli = .True.
      ! XCUTS
      ! ---cutoff parameters for x-sampling
      cvq = 1.8D0
      cdq = 2.D0
      csea = 0.5D0
      ssmima = 1.201D0
      ssmimq = ssmima**2
      vvmthr = 0.D0
      ! NOFINALE
      ! IFINAL = 1          !  NO FINALE CALL
      ! ==
      ifinal = &                   ! RECOMBIN
        0
      ! FINALE CALL
      irecom = 0
      lseadi = .True.
      ! SEASU3
      seasq = 0.5D0
      ! CRONINPT
      ! MKCRON = 1
      ! CRONCO = 0.64
      mkcron = 0
      cronco = 0.D0
      ! ALLPART
      ihada = .True.
      ! INTERDPM
      intdpm = 0
      iroeh = 0
      ! POPCORCK
      pdbck = 0.D0
      ijpock = 0
      ! CASADIQU
      icasad = 1
      ! CASAXX = 0.5D0
      casaxx = &                   ! POPCORSE
        0.05D0
      ! corrected Nov. 2001
      pdbse = &                    ! with baryon stopping
        0.45D0
      pdbseu = &                   ! PDBSE  = 0.D0               ! without
                                   ! baryon stopping
        0.45D0
      ! PDBSEU = 0.D0               ! without baryon stopping

      ! with baryon stopping
      irejck = 0
      irejse = 0
      irejs3 = 0
      irejs0 = 0
      ick4 = 0
      ise4 = 0
      ise43 = 0
      ihad4 = 0
      ick6 = 0
      ise6 = 0
      ihad6 = 0
      irejsa = 0
      ireja3 = 0
      ireja0 = 0
      isea4 = 0
      isea43 = 0
      ihada4 = 0
      isea6 = 0
      isea63 = 0
      ihada6 = 0
      ! POPCORN
      pdb = 0.1D0
      ajsdef = 0.D0
      ! FLUCTUAT
      ifluct = 0
      ! INTPT
      intpt = .True.

      ! HADRONIZ
      ihadrz = 2
      ifrag = 1

      ! charmed particles do not decay at the interaction vertex
      ipromu = 0
      If (ihadrz>=2) Then
        ifrag = ihadrz - 1
        ! LUNDIN must be initialized at each call of DPMJLNK
        ! to let charm decay in specialized decay routine
        ! CALL LUNDIN
        ! CALL TESLUN
      End If
      ! DIQUARKS
      idiqua = 1
      idiquu = 1
      amedd = 0.9D0
      ! SINGLECH
      isicha = 0
      ! EVAPORAT
      ! IEVAP = 0
      ievap = 1
      ! SEAQUARK
      ! sea quarks in multiple chains
      seaqx = 0.5D0
      ! sea quarks in Glauber events
      seaqxn = 0.5D0
      ! GLAUBERI
      ! GLAUBERA
      jglaub = 2
      ! HADRINTH
      ehadth = 5.D0
      ! POMTABLE
      ! IPOMTA = 0
      ipomta = 1

      ! CMHISTO
      cmhis = &                    ! CENTRAL
        1
      ! Lab System
      icentr = 0
      ! STRUCFUN
      istruf = 222
      istrum = 0
      istrut = istruf/100
      istruf = istruf - istrut*100
      istrum = istruf
      ! SINGDIFF
      isingd = 1
      isingx = 1
      idubld = 0
      sdfrac = 1.
      ! START
      ! NEVNTS passed to DTMAI as argument, NEVHAD in COMMON
      nfile = 0
      nstart = 1
      istar2 = 0
      istar3 = 0
      ptlar = 2.D0
      iglaub = 0

      ! the following initialization involve input cards for internal use
      ! (in alphabetical order)

      ! GLUSPLIT
      nugluu = 1
      nsgluu = 0
      ! PARTEV
      ! ITEST  = 0    (see SIGMAPOM)
      npev = 30
      nvers = 1
      ! SAMPT
      isampt = 4
      ! SELHARD

      dropjt = 0.D0
      ! ITEST  = 0    (see SIGMAPOM)
      iophrd = 2
      ptthr = 3.D0
      ptthr2 = ptthr
      cmener = 100.D0
      ! WRITE(6,*)' CMENER ',CMENER
      If (istrut==1) Then
        ptthr = 2.1D0 + 0.15D0*(log10(cmener/50.))**3
        ptthr2 = ptthr
      Else If (istrut==2) Then
        ptthr = 2.5D0 + 0.12D0*(log10(cmener/50.))**3
        ptthr2 = ptthr
      End If
      ! SIGMAPOM

      itest = 0
      ipim = 2
      icon = 48
      isig = 10
      ! default changed in relation to DPT =4
      lmax = 30
      mmax = 100
      nmax = 2
      difel = 0.D0
      difnu = 1.D0
      ! some options use special routines for MMAX=0 so far NMAX=0
      ! PSHOWER
      ipshow = 1
      ! SECINTER
      isecin = 0
      ! EVAPORATE
      ! set default if EVAP requested without "what-values"
      levprt = .True.
      ilvmod = 1
      ldeexg = .True.
      lheavy = .True.
      ! LFRMBK = .FALSE.
      lfrmbk = .True.
      ifiss = 0
      ! INITIALIZATION OF EVAPORATION MODULE
      nbertp = 14
      lunber = 14
      If (debug) Write (mdebug, *) 'DPMJIN: before NUCLEAR.BIN opening'
      ! INITIALIZATION OF EVAPORATION MODULE
      ! dh   04/05/2012
      ! File NUCLEAR.BIN is a binary file of length depending on system
      ! 198472  in 32-bit mode
      ! 198768  in 64-bit mode
      ! To use 32-bit mode files on 64-bit machines, use in compilation
      ! the flag:
      ! -frecord-marker=4  (Default is -frecord-marker=8)
!       datdir = '/Users/afedynitch/Documents/KIT/artifacts/matrix_method/iamdata/'

      Open (Unit=lunber, File=datdir(1:index(datdir, &
        ' ')-1)//'NUCLEAR.BIN', Status='OLD', Form='UNFORMATTED')
      If (debug) Write (mdebug, *) 'DPMJIN: before BERTTP NUCLEAR.BIN', &
        ' opened LUNBER= ', lunber
      Call berttp
      If (debug) Write (mdebug, *) 'DPMJIN: before INCINI'
      Call incini
      If (debug) Write (mdebug, *) 'DPMJIN: after INCINI'
      Close (Unit=nbertp)
      If (debug) Write (mdebug, *) 'DPMJIN: NUCLEAR.BIN closed'
      ! starting parameters not read in

      ! ISTART = .FALSE.
      xsmax = 0.8D0
      itopd = 0

      If (debug) Then
        ! Printout of important Parameters (defaults and input cards)
        Write (mdebug, *) ' Printout of important Parameters before DPMJET.', &
          ' Please note for DPMJET input all numbers are floating point!'
        Write (mdebug, *) 'PROJPAR  ', ip, ipz
        Write (mdebug, *) 'TARPAR   ', it, itz
        Write (mdebug, *) 'MOMENTUM ', ppn
        Write (mdebug, *) 'ENERGY   ', epn
        Write (mdebug, *) 'CMENERGY ', umo
        Write (mdebug, *) 'NOFINALE ', ifinal
        Write (mdebug, *) 'EVAPORAT ', ievap
        Write (mdebug, *) 'OUTLEVEL ', ipri, ipev, ippa, ipco, init, iphkk
        auauau = rd2out(iseed1, iseed2)
        Write (mdebug, *) 'RANDOMIZ ', iseed1, iseed2, &
          ' Initial RNDM(RM48) seeds'
        Write (mdebug, *) 'STRUCFUN ', istruf + 100*istrut
        Write (mdebug, *) 'SAMPT    ', isampt
        Write (mdebug, *) 'SELHARD  ', 0, iophrd, 0, dropjt, ptthr, ptthr2
        Write (mdebug, *) 'SIGMAPOM ', 0, isig, ipim + 10*icon, imax, mmax, &
          nmax
        Write (mdebug, *) 'PSHOWER  ', ipshow
        Write (mdebug, *) 'CENTRAL  ', icentr
        Write (mdebug, *) 'CMHISTO  ', cmhis
        Write (mdebug, *) 'SEASU3   ', seasq
        Write (mdebug, *) 'RECOMBIN ', irecom
        Write (mdebug, *) 'SINGDIFF ', isingd
        Write (mdebug, *) 'TAUFOR   ', taufor, ktauge, itauve
        Write (mdebug, *) 'POPCORN  ', pdb
        Write (mdebug, *) 'POPCORCK ', ijpock, pdbck
        Write (mdebug, *) 'POPCORSE ', pdbse, pdbseu
        Write (mdebug, *) 'CASADIQU ', icasad, casaxx
        Write (mdebug, *) 'DIQUARKS ', idiqua, idiquu, amedd
        Write (mdebug, *) 'HADRONIZ ', ihadrz
        Write (mdebug, *) 'INTPT    ', intpt
        Write (mdebug, *) 'PAULI    ', lpauli
        Write (mdebug, *) 'FERMI    ', fermp, fermod
        Write (mdebug, *) 'CRONINPT ', mkcron, cronco
        Write (mdebug, *) 'SEADISTR ', xseacu + 0.95D0, unon, unom, unosea
        Write (mdebug, *) 'SEAQUARK ', seaqx, seaqxn
        Write (mdebug, *) 'SECINTER ', isecin
        Write (mdebug, *) 'XCUTS    ', cvq, cdq, csea, ssmima
        Write (mdebug, *) 'START    ', ncases
        Write (mdebug, *) ' Printout of important Parameters before DPMJET', &
          ' Please note for DPMJET input all numbers are floating point!'
        ! Printout of important Parameters (defaults and input cards)
      End If
      ! LAB SYSTEM
      Call distr(1, ijproj, ppn, idummy)

      ! UNIT 47 CONTAINS THE MATERIAL DEPENDEND VALUES OF 'BSITEN'
      Open (Unit=47, File=datdir(1:index(datdir,' ')-1)//'GLAUBTAR.DAT', &
        Status='UNKNOWN')
      ! Read Glauber Data from file GLAUBTAR for p on Nitrogen
      ip = 1
      ipz = 1
      it = 14
      itz = 7
      ! KKMAT = 1
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for p on Oxygen
      it = 16
      itz = 8
      ! KKMAT = 2
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for p on Argon
      it = 40
      itz = 18
      ! KKMAT = 3
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for 4He on Nitrogen
      ip = 4
      ipz = 2
      it = 14
      itz = 7
      ! KKMAT = 4
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for 4He on Oxygen
      it = 16
      itz = 8
      ! KKMAT = 5
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for 4He on Argon
      it = 40
      itz = 18
      ! KKMAT = 6
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for 14N on Nitrogen
      ip = 14
      ipz = 7
      it = 14
      itz = 7
      ! KKMAT = 7
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for 14N on Oxygen
      it = 16
      itz = 8
      ! KKMAT = 8
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for 14N on Argon
      it = 40
      itz = 18
      ! KKMAT = 9
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for Fe on Nitrogen
      ip = 56
      ipz = 26
      it = 14
      itz = 7
      ! KKMAT = 10
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for Fe on Oxygen
      it = 16
      itz = 8
      ! KKMAT = 11
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for Fe on Argon
      it = 40
      itz = 18
      ! KKMAT = 12
      Call shmakf(ip, ipz, it, itz)
      ! Read Glauber Data from file GLAUBTAR for p on Proton
      ip = 1
      ipz = 1
      it = 1
      itz = 1
      ! KKMAT = 13
      Call shmakf(ip, ipz, it, itz)
      Close (Unit=47)
      If (debug) Write (mdebug, *) 'DPMJIN: GLAUBTAR closed'

      cmener = 1000.D0
      If (ipim==2) Call prblm2(cmener)
      ! initialize hard scattering
      Call jtdtu(0)
      ! initialize transverse momenta for soft scattering
      Call samppt(0, pt)

      Return
    End Subroutine dpmjin

!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      Integer Function mcihad(mcind)
      Implicit Double Precision (a-h,o-z)
!*** calculation of the particle index according to the pdg proposal.
!    mcihad particle number as in decay, bamjet, hadevt, fluka etc.
!    mcind  pdg particle number
      integer hamcin
      common /hamcin/ hamcin(410)
      
      ih=0
      mcihad=0
      
      If((mcind == 0).or.(mcind > 70000)) Return
      
      Do i=1,410
        ih=i
        if (hamcin(i) == mcind) Exit
      End Do

      mcihad = ih
      
      Return
      End Function mcihad

!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      Integer Function mpdgha(mcind)
      Implicit Double Precision (a-h,o-z)
!*** calculation of the particle index according to the pdg proposal.
!    mcind  particle number as in decay, bamjet, hadevt, fluka etc.
!    mpdgha pdg particle number
      Integer hamcin
      Common /hamcin/ hamcin(410)
      
      mpdgha = hamcin(mcind)
      
      Return
      End Function mpdgha

      Double precision function siinel(kproj,ktarg,umo)
      implicit double precision (a-h,o-z)
!     inelastic cross section
      siinel=dshnto(kproj,ktarg,umo)-dshnel(kproj,ktarg,umo)
      return
      end Function siinel
    