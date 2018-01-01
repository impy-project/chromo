   MODULE Terms

     IMPLICIT NONE
     REAL(KIND(0D0)) :: mass,Tem,mui,dsttilde,dsxtilde,dsztilde,gamma,vtilde
     INTEGER :: pm

   END MODULE Terms


   SUBROUTINE Peaks(fmax,fpxmax,fpymax,fpzmax,m,a,T,mu,  &
                    vx,vy,vz,dst,dsx,dsy,dsz)

     IMPLICIT NONE
     REAL(KIND(0D0)) :: fmax,fpxmax,fpymax,fpzmax
     REAL(KIND(0D0)) :: m, T, mu
     REAL(KIND(0D0)) :: vx,vy,vz,dst,dsx,dsy,dsz
     INTEGER :: a, info

     CALL CFmax(fmax,fpxmax,fpymax,fpzmax,m,a,T,mu,                    &
                vx,vy,vz,dst,dsx,dsy,dsz,1,info)
     IF ((fmax .lt. 0D0).or.(info .eq. 0).or.(info .gt. 3)) THEN
        CALL CFmax(fmax,fpxmax,fpymax,fpzmax,m,a,T,mu,                 &
                   vx,vy,vz,dst,dsx,dsy,dsz,2,info)
        IF ((fmax .lt. 0D0).or.(info .eq. 0).or.(info .gt. 3)) THEN
           CALL BackUp(fmax,fpxmax,fpymax,fpzmax,m,a,T,mu,             &
                       vx,vy,vz,dst,dsx,dsy,dsz)
        END IF
     END IF

   END SUBROUTINE Peaks


   SUBROUTINE CFmax(fmax,fpxmax,fpymax,fpzmax,m,a,             &
                    T,mu,vx,vy,vz,dst,dsx,dsy,dsz,round,info)

     USE Terms
     IMPLICIT NONE
     REAL(KIND(0D0)) :: fmax,fpxmax,fpymax,fpzmax
     REAL(KIND(0D0)) :: m, T, mu
     REAL(KIND(0D0)) :: vx,vy,vz,dst,dsx,dsy,dsz
     INTEGER         :: a,round,info

     REAL(KIND(0D0)) :: vT
     REAL(KIND(0D0)) :: cosPhi,sinPhi,cosTheta,sinTheta,cosPhi2,sinPhi2
     REAL(KIND(0D0)) :: dsxprime,dsyprime,dszprime
     REAL(KIND(0D0)) :: dsx2prime,dsy2prime,dsz2prime
     REAL(KIND(0D0)) :: dsytilde
     REAL(KIND(0D0)) :: pxtilde,pytilde,pztilde,px2,py2,pz2,px1,py1,pz1
     REAL(KIND(0D0)) :: fp, p0, alpha

     REAL(KIND(0D0)),DIMENSION(2) :: P,F
     REAL(KIND(0D0)),DIMENSION(10) :: Work
     REAL(KIND(0D0)),PARAMETER :: Ftol = 1D-9, Xtol = 1D-9
     REAL(KIND(0D0)) :: Ft, Xt, dsig
     INTEGER  :: i
     EXTERNAL :: Derivatives

     mass = m
     pm   = a
     Tem  = T
     mui  = mu

! Rotate the system to make it easier to handle
     vtilde = vx**2+vy**2+vz**2
     gamma = 1/Sqrt(1-vtilde)
     vtilde = Sqrt(vtilde)
     vT = Sqrt(vx**2+vy**2)

     IF (ABS(vT) .gt. Ftol) THEN
        cosPhi = vx/vT
        sinPhi = vy/vT
     ELSE
        cosPhi = 1D0
        sinPhi = 0D0
     END IF

     dsig = dst + Sqrt(dsx**2+dsy**2+dsz**2)

     dsttilde =  dst/dsig
     dsxprime = ( cosPhi*dsx + sinPhi*dsy)/dsig
     dsyprime = (-sinPhi*dsx + cosPhi*dsy)/dsig
     dszprime =  dsz/dsig

     IF (ABS(vtilde) .gt. Ftol) THEN
        cosTheta = vz/vtilde
        sinTheta = vT/vtilde
     ELSE
        cosTheta = 1D0
        sinTheta = 0D0
     END IF

     dsx2prime = cosTheta*dsxprime - sinTheta*dszprime
     dsy2prime = dsyprime
     dsz2prime = sinTheta*dsxprime + cosTheta*dszprime

     dsxtilde = Sqrt(dsx2prime**2+dsy2prime**2)
     dsytilde = 0D0
     dsztilde = dsz2prime
     IF (ABS(dsxtilde) .gt. Ftol) THEN
        cosPhi2 = dsx2prime/dsxtilde
        sinPhi2 = dsy2prime/dsxtilde
     ELSE
        cosPhi2 = 1D0
        sinPhi2 = 0D0
     END IF

! Find the location

     IF (round .eq. 1) THEN
        P(1) = 0D0
        P(2) = m*vtilde*gamma
        F = 0D0
     ELSE
        P(1) = dsxtilde/Sqrt(dsxtilde**2+dsztilde**2)
        P(2) = dsztilde/Sqrt(dsxtilde**2+dsztilde**2)
     END IF

     CALL DSNLEQ90(2,P,F,Ftol,Xtol,250,0,Info,Derivatives,Work)

     pxtilde = P(1)
     pytilde = 0D0
     pztilde = P(2)

! rotate back to the orginal frame

     px2 = cosPhi2*pxtilde - sinPhi2*pytilde
     py2 = sinPhi2*pxtilde + cosPhi2*pytilde
     pz2 = pztilde

     px1 =  cosTheta*px2 + sinTheta*pz2
     py1 =  py2
     pz1 = -sinTheta*px2 + cosTheta*pz2

     fpxmax = cosPhi*px1 - sinPhi*py1
     fpymax = sinPhi*px1 + cosPhi*py1
     fpzmax = pz1

! evaluate the maximum

     p0 = Sqrt(m**2+fpxmax**2+fpymax**2+fpzmax**2)
     fmax = fp(p0,fpxmax,fpymax,fpzmax,vx,vy,vz,gamma,dst,dsx,dsy,dsz,a,T,mu)

   END SUBROUTINE CFmax


   SUBROUTINE Derivatives(i,P,F,K)

     USE Terms
     IMPLICIT NONE
     INTEGER :: i,K
     REAL(KIND(0D0)) :: P(*),F(*)
     REAL(KIND(0D0)) :: px,pz,E,dsigmupmu,pdotu,fprime

     px = P(1)
     pz = P(2)

     E = Sqrt(mass**2+px**2+pz**2)
     dsigmupmu = dsttilde*E + px*dsxtilde + pz*dsztilde
     pdotu = gamma*(E-pz*vtilde)
     fprime = gamma/Tem*E/(1+pm*Exp(-(pdotu-mui)/Tem))

     SELECT CASE (K)
     CASE (1)
        F(1) = (E**2-px**2)*dsxtilde - pz*px*dsztilde    &
              -dsigmupmu*px*fprime
     CASE (2) 
        F(2) = (E**2-pz**2)*dsztilde - pz*px*dsxtilde    &
              -dsigmupmu*(pz-vtilde*E)*fprime
     END SELECT

   END SUBROUTINE Derivatives


   SUBROUTINE BackUp(fmax,fpxmax,fpymax,fpzmax,m,a,T,mu,  &
                     vx,vy,vz,dst,dsx,dsy,dsz)

     implicit none
     real*8 fmax,fpxmax,fpymax,fpzmax
     real*8 vx,vy,vz,g,dst,dsx,dsy,dsz,T,mu,m
     integer a
     real*8 intpxmin,intpymin,intpzmin,dp,fp
     real*8 dist,p0,px,py,pz
     integer ix,iy,iz

     fmax=0d0    
     intpxmin=-2.d0
     intpymin=-2.d0
     intpzmin=-8.d0
     dp=0.1d0
     g = 1/Sqrt(1-vx**2-vy**2-vz**2)

     do ix=1,40
        do iy=1,40
           do iz=1,160  
              px=intpxmin+ix*dp
              py=intpymin+iy*dp
              pz=intpzmin+iz*dp
              p0=sqrt(m**2+px**2+py**2+pz**2)
      
              dist=fp(p0,px,py,pz,vx,vy,vz,g,dst,dsx,dsy,dsz,a,T,mu)
              if(dist.gt.fmax)then
                 fmax=dist
                 fpxmax=px
                 fpymax=py
                 fpzmax=pz
              endif
           end do
        end do
     end do

!.. ensure maximum is large enough 
     fmax=fmax*1.2d0        

   END SUBROUTINE BackUp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine DSNLEQ of CERNLIB. Written in the free format of Fortran90
! by P. Huovinen, Frankfurt, July 2012
!
! $Id: snleq64.F,v 1.1.1.1 1996/04/01 15:01:52 mclareni Exp $
!
! $Log: snleq64.F,v $
! Revision 1.1.1.1  1996/04/01 15:01:52  mclareni
! Mathlib gen
!
!
    SUBROUTINE DSNLEQ90(N,X,F,Ftol,Xtol,MaxF,Iprt,Info,SUB,W)
!
! $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $

! $Log: imp64.inc,v $
! Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
! Mathlib gen
!
!
! imp64.inc
!
      IMPLICIT NONE
!     Based on   J.J. More  and  M.Y. Cosnard
!
!       ALGORITHM 554 BRENTM, A Fortran Subroutine for the
!       Numerical Solution of Systems of Nonlinear Equations [C5]
!
!     ACM Trans. Math. Software 6 (1980) 240-251.
      INTEGER :: N
      REAL(KIND(0D0)),DIMENSION(N) :: X, F
      REAL(KIND(0D0)) :: Xtol, Ftol
      INTEGER :: MaxF, Iprt, Info
      REAL(KIND(0D0)),DIMENSION(N,*) :: W

      LOGICAL :: LCV
      REAL(KIND(0D0)),PARAMETER :: Z1 = 1D0, Scale = 1D1
      REAL(KIND(0D0)),PARAMETER :: P05 = 5*Z1/1D2
!**** EPS = SQRT(SMALLEST FP.NUMBER)
!     EPS = 1 / SQRT( 16D0**13 )
      REAL(KIND(0D0)),PARAMETER :: Eps =  0.149011611938476600D-7
!      INTEGER,DIMENSION(288),PARAMETER :: MPT =                             &
!           (/ 1,2,(3,i=1,3),(4,i=1,3),(5,i=1,4)
!           (/ 1* 1,1* 2,3* 3,3* 4,4* 5,4* 6,4* 7,4* 8,5* 9,5*10,5*11,5*12,  &
!              5*13,5*14,6*15,6*16,5*17,6*18,6*19,6*20,7*21,6*22,6*23,7*24,  &
!              6*25,7*26,6*27,7*28,7*29,7*30,7*31,7*32,7*33,7*34,7*35,7*36,  &
!              8*37,7*38,7*39,8*40,7*41,8*42,7*43,8*44,8*45,7*46,8*47,8*48 /)
      INTEGER :: i, j, k, m
      INTEGER :: Mopt, Iflag, Numf, Nfcall, Nier6, Nier7, Nier8, Nsing 
      REAL(KIND(0D0)) :: H, Temp, Fnorm, Fnorm1, Xnorm, Difit, Difit1, Delta
      REAL(KIND(0D0)) :: Fky, Fkz, Eta, Sknorm
      INTEGER,DIMENSION(288) :: MPT
      DATA (MPT(i),i=1,288)                                          &
      /1* 1,1* 2,3* 3,3* 4,4* 5,4* 6,4* 7,4* 8,5* 9,5*10,5*11,5*12,  &
       5*13,5*14,6*15,6*16,5*17,6*18,6*19,6*20,7*21,6*22,6*23,7*24,  &
       6*25,7*26,6*27,7*28,7*29,7*30,7*31,7*32,7*33,7*34,7*35,7*36,  &
       8*37,7*38,7*39,8*40,7*41,8*42,7*43,8*44,8*45,7*46,8*47,8*48/

      Info = 0
      IF ((N .le. 0).or.(Ftol .le. 0).or.(Xtol .le. 0)) RETURN
!
!     Find optimal Mopt for iterative refinement
!
      IF(N .le. 288) THEN
         Mopt = MPT(N)
      ELSE
         H = 0D0
         DO i = 49,N
            Temp = LOG(i+Z1)/(N+2*i+1)
            IF(Temp .lt. H) THEN
               Mopt = i-1
               EXIT
            ENDIF
            H = Temp
         END DO
      ENDIF
      Iflag = 0
      Numf = 0
      Nfcall = 0
      Nier6 = -1
      Nier7 = -1
      Nier8 = 0
      Fnorm = 0
      Difit = 0
      Xnorm = MAXVAL(ABS(X))
      Delta = Scale*Xnorm
      IF(Xnorm .eq. 0) Delta = Scale
      Iterate: DO
20       CONTINUE
         IF (Iprt .ne. 0) WRITE(6,'(1X,I5,D25.14)') (i,X(i),i=1,N)
         Nsing = N
         Fnorm1 = Fnorm
         Difit1 = Difit
         Fnorm = 0D0
!
!     Compute step H for the divided difference which approximates
!     the K-th row of the Jacobian matrix
!
         H = Eps*Xnorm
         IF (H .eq. 0D0) H = Eps
         DO j = 1,N
            W(:,j+3) = 0D0
            W(j,j+3) = H
            W(j,2) = X(j)
         END DO
!
!     Enter a subiteration
!
         DO k = 1,N
            Iflag = k
            CALL SUB(N,W(1,2),F,Iflag)
            Fky = F(k)
            Nfcall = Nfcall+1
            Numf = Nfcall/N
            IF(Iflag .lt. 0) GO TO 230
            Fnorm = MAX(Fnorm,ABS(Fky))

!     Compute the K-th row of the Jacobian matrix
!
            DO j = k,N
               W(:,3) = W(:,2) + W(:,j+3)
               CALL SUB(N,W(1,3),F,Iflag)
               Fkz = F(k)
               Nfcall = Nfcall + 1
               Numf = Nfcall/N
               IF (Iflag .lt. 0) GO TO 230
               W(j,1) = Fkz - Fky
            END DO
            F(k) = Fky
!
!     Compute the Householder transformation to reduce the K-th row
!     of the Jacobian matrix to a multiple of the K-th unit vector
!
            Eta = MAXVAL(ABS(W(k:N,1)))
            IF (Eta .ne. 0D0) THEN
               Nsing = Nsing - 1
               W(k:N,1) = W(k:N,1)/Eta
               Sknorm = SQRT(SUM(W(k:N,1)**2))
               IF(W(k,1) .lt. 0D0) Sknorm = -Sknorm
               W(k,1) = W(k,1) + Sknorm
!
!     Apply the transformation
!
               W(:,3) = 0D0
               DO j = k,N
                  W(:,3) = W(:,3) + W(j,1)*W(:,j+3)
               END DO
               DO j = k,N
                  Temp = W(j,1)/(Sknorm*W(k,1))
                  W(:,j+3) = W(:,j+3) - Temp*W(:,3)
               END DO
!
!     Compute the subiterate
!
               W(k,1) = Sknorm*Eta
               Temp = Fky/W(k,1)
               IF(H*ABS(Temp) .gt. Delta) Temp = SIGN(Delta/H,Temp)
               W(:,2) = W(:,2) + Temp*W(:,k+3)
            END IF
         END DO
!
!     Compute the norms of the iterate and correction vector
!
         Xnorm = MAXVAL(ABS(W(:,2)))
         Difit = MAXVAL(ABS(X-W(:,2)))
         X = W(:,2)
!
!     Update the bound on the correction vector
!
         Delta = MAX(Delta,Scale*Xnorm)
!
!     Determine the progress of the iteration
!
         LCV = ((Fnorm .lt. Fnorm1).and.(Difit .lt. Difit1).and.(Nsing .eq. 0))
         Nier6 = Nier6 + 1
         Nier7 = Nier7 + 1
         Nier8 = Nier8 + 1
         IF (LCV) Nier6 = 0
         IF ((Fnorm .lt. Fnorm1).or.(Difit .lt. Difit1)) Nier7 = 0
         IF (Difit .gt. Eps*Xnorm) Nier8 = 0
!
!     Tests for convergence
!
         IF (Fnorm .le. Ftol) Info = 1
         IF ((Difit .le. Xtol*Xnorm).and.(LCV)) Info = 2
         IF ((Fnorm .le. Ftol).and.(Info .eq. 2)) Info = 3
         IF (Info .ne. 0) EXIT Iterate
!
!     Tests for termination
!
         IF (Numf .ge. MaxF) Info = 4
         IF (Nsing .eq. N) Info = 5
         IF (Nier6 .eq. 5) Info = 6
         IF (Nier7 .eq. 3) Info = 7
         IF (Nier8 .eq. 4) Info = 8
         IF (Info .ne. 0) EXIT Iterate
         IF ((LCV).and.(Difit .le. P05*Xnorm)) THEN
!           IF ((.NOT.LCV).or.(Difit .gt. P05*Xnorm)) GO TO 20
!
!     Iterative refinement  (if the iteration is converging)
!
            DO m = 2,Mopt
               Fnorm1 = Fnorm
               Fnorm = 0D0
               DO k = 1,N
                  Iflag = k
                  CALL SUB(N,W(1,2),F,Iflag)
                  Fky = F(k)
                  Nfcall = Nfcall+1
                  Numf = Nfcall/N
                  IF (Iflag .lt. 0) GO TO 230
                  Fnorm = MAX(Fnorm,ABS(Fky))
!
!     Iterative refinement is terminated if it does not give a
!     reduction on residuals
!
                  IF(Fnorm .ge. Fnorm1) THEN
                     Fnorm = Fnorm1
                     GO TO 20
                  ENDIF
                  Temp = Fky/W(k,1)
                  W(:,2) = W(:,2) + Temp*W(:,k+3)
               END DO
!
!     Compute the norms of the iterate and correction vector
!
               Xnorm = MAXVAL(ABS(W(:,2)))
               Difit = MAXVAL(ABS(X-W(:,2)))
               X = W(:,2)
!
!     Stopping criteria for iterative refinement
!
               IF (Fnorm .le. Ftol) Info = 1
               IF (Difit .le. Xtol*Xnorm) Info = 2
               IF ((Fnorm .le. Ftol).and.(Info .eq. 2)) Info = 3
               IF ((Numf .ge. MaxF).and.(Info .eq. 0)) Info = 4
               IF (Info .ne. 0) EXIT Iterate
            END DO
         END IF
         IF (Info .ne. 0) EXIT Iterate
      END DO Iterate

230   IF (Iflag .lt. 0) Info = Iflag

    END SUBROUTINE DSNLEQ90
