      SUBROUTINE FORCE_VECTORS(XRATE,N1,N2)
      IMPLICIT NONE
      SAVE
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

c     external types
      double precision xrate
      integer n1,n2

c     internal types
      integer ipi2vec,lcon,lreschex,ll,la,la_new,i,j,kba
      DIMENSION IPI2VEC(99)
      double precision pts,pz2,xmts,xf,xfs,s_rndm
      
      DIMENSION LCON(6:43),LRESCHEX(6:39)
c     charge exchange map, i.e. pip -> pi0 ...
      DATA LCON /7,6,6,22,21,9,9,14,13,4*0,20,19,10,9,23,24,27,27,25,
     &     31,30,29,28,32,33,35,34,35,38,37,39,41,42,41,42/
c     charge and spin exchange map, i.e. pip -> rho0
c     approximate, proton and neutron should go to N(1520) not Delta
      DATA LRESCHEX /26,27,27,31,30,9,9,42,41,19*0,45,44,45,48,47,39/ 
      INTEGER IFIRST
      DATA IFIRST /0/

      if(ifirst.eq.0)then
         print *,'initializing..'
         do j=1,99
            IPI2VEC(J) = J
         enddo
         IPI2VEC(6) = 27
         IPI2VEC(7) = 25
         IPI2VEC(8) = 26
         ifirst = 1
      endif

      KBA = IABS(KB)
      
      select case(ipar(45))

c     trivial exchange model      
      case (1)
         do i=N1,N2
c     replace pions with vector mesons
            ll = mod(llist(I),10000)
            la = abs(ll)
            IF(s_rndm(la).lt.xrate)then
c     put back on mass shell
               la_new = IPI2VEC(la)
               xmts = p(i,1)**2 + p(i,2)**2 + am2(la_new)
               pz2 = p(i,4)**2 - xmts
               if(pz2.gt.EPS8)then
                  p(i,3) = sign(sqrt(pz2),p(i,3))
                  p(i,5) = am(la_new)
                  LLIST(I) = ISIGN(la_new,ll)
               endif
            endif
         enddo
         
c     large xf only, neutral pions only
      case (2)
         do i=N1,N2
            ll = mod(llist(I),10000)
            la = abs(ll)
            IF(la.eq.6)then
               xf = TWO*p(i,3)/SQS
               IF(s_rndm(la).lt.xrate*xf)then
c     exhcange and put back on mass shell
                  la_new = IPI2VEC(la)
                  xmts = p(i,1)**2 + p(i,2)**2 + am2(la_new)
                  pz2 = p(i,4)**2 - xmts
                  if(pz2.gt.EPS8)then
                     p(i,3) = sign(sqrt(pz2),p(i,3))
                     p(i,5) = am(la_new)
                     LLIST(I) = ISIGN(la_new,ll)
                  endif
               endif
            endif
         enddo

c     large xf only, charge and spin exchange
      case (3)
         do i=N1,N2
            ll = mod(llist(I),10000)
            la = abs(ll)
            IF(ll.eq.LCON(KBA))then
               xf = TWO*p(i,3)/sqs
               IF(s_rndm(la).lt.xrate*xf)then
c     replace charge exchange product of beam with
c     charge and spin exchange product, i.e.
c     pip-beam -> rho0 instead of pip-beam -> pi0
c     so replace pi0 with rho0 in final state
                  la_new = LRESCHEX(KBA)
c     put back on mass shell
                  xmts = p(i,1)**2 + p(i,2)**2 + am2(la_new)
                  pz2 = p(i,4)**2 - xmts
                  if(pz2.gt.EPS8)then
                     p(i,3) = sign(dsqrt(pz2),p(i,3))
                     p(i,5) = am(la_new)
                     LLIST(I) = ISIGN(la_new,ll)
                  endif
               endif
            endif
         enddo

c     large xf only, charge and spin exchange
      case (4)
         do i=N1,N2
            ll = mod(llist(I),10000)
            la = abs(ll)
            IF(ll.eq.LCON(KBA))then
               xf = TWO*p(i,3)/sqs
               xfs = xf ** 2
               IF(s_rndm(la).lt.xrate*xfs)then
c     replace charge exchange product of beam with
c     charge and spin exchange product, i.e.
c     pip-beam -> rho0 instead of pip-beam -> pi0
c     so replace pi0 with rho0 in final state
                  la_new = LRESCHEX(KBA)
c     put back on mass shell
                  xmts = p(i,1)**2 + p(i,2)**2 + am2(la_new)
                  pz2 = p(i,4)**2 - xmts
                  if(pz2.gt.EPS8)then
                     p(i,3) = sign(dsqrt(pz2),p(i,3))
                     p(i,5) = am(la_new)
                     LLIST(I) = ISIGN(la_new,ll)
                  endif
               endif
            endif
         enddo

      end select
      if(ndebug.ge.5) call sib_list(6)
      END
