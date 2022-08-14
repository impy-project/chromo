   SUBROUTINE Cornelius(E0,HyperCube,Ngp,dSigma,Nsurf,Vmid,dt,dx,dy,dz,  &
                        Nambi,Ndisc)
!
!  version 1.3 zeta
!
! This routine searches for a 3-dimensional isosurface of constant X in
! a 4-dimensional space when the values of X are known at the vertices
! (=corners) of the hypercube and X is interpolated linearly between the
! vertices. I.e. the usual problem of finding the freeze-out surface.
!
!    Variables: E0:      the value of X on the surface (INPUT)
!               HyperCube:   4D cube (2x2x2x2) to store the values of X (INPUT)
!               Ngp:     number of hypercube corners, i.e. gridpoints
!                        within the freeze-out surface, i.e. X >= E0 (INPUT)
!               dSigma:  4-vector table to store the normal vector(s) of
!                        the surface element(s), |dSigma| is the volume
!                        (=hyperarea) of the surface element (OUTPUT)
!               Nsurf:   Number of surface elements within the hypercube (OUT)
!               Vmid:    coordinates of the approximate centroid(s) of the
!                        element(s) if the origin is at (0,0,0,0) corner of
!                        the hypercube (OUTPUT)
!               dt, dx, dy, dz: the lengths of the edges of the hypercube (IN)
!               Nambi:   number of ambiguous faces on surfaces (INOUT)
!                        (do not set this to zero between successive calls)
!               Ndisc:   number of disconnected surface-elements so far (INOUT)
!                        (do not set this to zero between successive calls)
!
!        -- P. Huovinen, Jyvaskyla-Frankfurt, March 2005-Jan 2009 --
!
! Changes in version 1.1:
! - new name
! - rule for treating ambiguous surfaces
! - ability to recognize more than one surface within the hypercube
! - renamed variables, rewrote the code using more subroutines
!        -- P. Huovinen, Seattle-Ames-Frankfurt, June-July 2010 --
!
! Changes in version 1.2:
! - checks that the value Ngp makes sense
! - variables' intent (input or output) specified
! - allows 1-15 corners to be exactly at the freeze-out temperature by assuming
!   that such a corner is within the surface and setting the surface a distance
!   1D-9*dr from the corner towards the lower value
!        -- PH, Frankfurt-Ryanair-VR-Varkaus, Oct 2010 --
!
! Changes in version 1.3:
! - improves the handling of case where a corner is exactly at the freeze-out
!   temperature and the face is ambiguous. In previous version the behaviour
!   dependend on which corner was at the FO temperature, now the treatment is
!   consistent depending only on the value in the center of the face
!        -- PH, Frankfurt, July 2011 --
!
! The ordering of values in HyperCube(t,i,j,k):
! first index time, second x, third y, fourth z
!  HyperCube(0,0,0,0) <-> t=0,x=0,y=0,z=0
!  HyperCube(0,0,0,1) <-> t=0,x=0,y=0,z=1
!  HyperCube(0,0,1,0) <-> t=0,x=0,y=1,z=0
!  HyperCube(0,0,1,1) <-> t=0,x=0,y=1,z=1
!  HyperCube(0,1,0,0) <-> t=0,x=1,y=0,z=0
!  HyperCube(0,1,0,1) <-> t=0,x=1,y=0,z=1
!  HyperCube(0,1,1,0) <-> t=0,x=1,y=1,z=0
!  HyperCube(0,1,1,1) <-> t=0,x=1,y=1,z=1
!  HyperCube(1,0,0,0) <-> t=1,x=0,y=0,z=0
!  HyperCube(1,0,0,1) <-> t=1,x=0,y=0,z=1
!  HyperCube(1,0,1,0) <-> t=1,x=0,y=1,z=0
!  HyperCube(1,0,1,1) <-> t=1,x=0,y=1,z=1
!  HyperCube(1,1,0,0) <-> t=1,x=1,y=0,z=0
!  HyperCube(1,1,0,1) <-> t=1,x=1,y=0,z=1
!  HyperCube(1,1,1,0) <-> t=1,x=1,y=1,z=0
!  HyperCube(1,1,1,1) <-> t=1,x=1,y=1,z=1
!
     IMPLICIT NONE

     REAL(KIND(0D0)),INTENT(IN) :: E0
     REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1,0:1),INTENT(IN) :: HyperCube
     INTEGER,INTENT(IN)         :: Ngp
     REAL(KIND(0D0)),DIMENSION(0:3,8),INTENT(OUT) :: dSigma
     INTEGER,INTENT(OUT)        :: Nsurf
     REAL(KIND(0D0)),DIMENSION(0:3,8),INTENT(OUT) :: Vmid
     REAL(KIND(0D0)),INTENT(IN) :: dt,dx,dy,dz
     INTEGER,INTENT(INOUT)      :: Nambi, Ndisc

     REAL(KIND(0D0)),DIMENSION(0:3,3,96) :: Edge    ! Table for ends of edges  
                                              ! i.e. corners of the polyhedra
     INTEGER :: Nedge                               ! # of edges
     LOGICAL :: Ambig                               ! ambiguous structure
     REAL(KIND(0D0)),DIMENSION(0:3,96) :: Out       ! outside direction
     LOGICAL :: Pathological                        ! Too difficult structure
     LOGICAL :: Suspicious                 ! Necessary to check connectedness ?
     INTEGER,DIMENSION(9) :: EdgeSet ! which edges belong to different surfaces 
     INTEGER :: i



     Nedge = 0
     dSigma = 0D0
     Ambig = .false.
     Pathological = .false.
     IF ((Ngp .lt. 1).or.(Ngp.gt.15)) CALL StrangeN(Ngp)

     CALL Cubes(E0,HyperCube,Edge,Nedge,Out,dt,dx,dy,dz,Ambig,Pathological)

     IF (Pathological) CALL DeadEnd(HyperCube,E0)

     IF (Ambig) THEN
        Nambi = Nambi + 1
!        Write(*,*) 'Ambiguous structure number',Nambi
     END IF

     IF ((Nedge .gt. 96).or.(Nedge .lt. 12)) THEN
        WRITE(*,*) 'Too many (>96) or too few (<12) tetrahedra'
        WRITE(*,*) 'Ntetra =',Nedge
        CALL DeadEnd(HyperCube,E0)
     END IF

     IF (.not.(Ambig)) Ambig = Suspicious(Nedge,Ngp)

     IF (Ambig) THEN
        CALL Disconnected3D(Nedge,Edge,Out,Nsurf,EdgeSet)
     ELSE
        Nsurf = 1
        EdgeSet(1) = 1
        EdgeSet(2) = Nedge+1
     END IF
     IF (Nsurf .gt. 1) THEN
        Ndisc = Ndisc+1
!        Write(*,*) 'Disconnected surface element number',Ndisc
     END IF

     DO i = 1,Nsurf
        Nedge = EdgeSet(i+1)-EdgeSet(i)
        IF (Nedge .lt. 12) THEN
           WRITE(*,*) Nsurf,' surfaces in hypercube.'
           WRITE(*,*) i,'. of them has',Nedge,' edges. Weird.'
           CALL DeadEnd(HyperCube,E0)
        END IF
        CALL NormalVector(i,EdgeSet(i),EdgeSet(i+1)-1,Nedge,Edge,Out,  &
                          dSigma,Vmid)
     END DO

!     CALL PutPut(Edge,EdgeSet,Vmid,Nsurf,dt,dx,dz,dz)


   END SUBROUTINE Cornelius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE Cubes(E0,HyperCube,Edge,Nedge,Out,dt,dx,dy,dz,Ambig,Pathological)

     IMPLICIT NONE
     REAL(KIND(0D0)) :: E0,dt,dx,dy,dz
     REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1,0:1) :: HyperCube
     REAL(KIND(0D0)),DIMENSION(0:3,3,96) :: Edge
     REAL(KIND(0D0)),DIMENSION(0:3,96) :: Out
     INTEGER         :: Nedge
     LOGICAL         :: Ambig,Pathological

     REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1) :: Cube
     REAL(KIND(0D0)),DIMENSION(0:2,3,12)    :: CubEdge
     REAL(KIND(0D0)),DIMENSION(0:2,12)      :: Ut
     INTEGER :: NcEdge, i

     DO i = 0,1 ! t = 0 or 1
        Cube = HyperCube(i,:,:,:)
        IF (ANY(Cube .ge. E0) .and. ANY(Cube .le. E0)) THEN
           CALL CubeCut(E0,Cube,CubEdge,NcEdge,Ut,dx,dy,dz,Ambig,Pathological)
           IF (Pathological) RETURN
           CALL StoreCube(Edge,Nedge,Out,CubEdge,NcEdge,Ut,i*dt,0,(/1,2,3/))
        END IF
     END DO

     DO i = 0,1 ! x = 0 or 1
        Cube = HyperCube(:,i,:,:)
        IF (ANY(Cube .ge. E0) .and. ANY(Cube .le. E0)) THEN
           CALL CubeCut(E0,Cube,CubEdge,NcEdge,Ut,dt,dy,dz,Ambig,Pathological)
           IF (Pathological) RETURN
           CALL StoreCube(Edge,Nedge,Out,CubEdge,NcEdge,Ut,i*dx,1,(/0,2,3/))
        END IF
     END DO

     DO i = 0,1 ! y = 0 or 1
        Cube = HyperCube(:,:,i,:)
        IF (ANY(Cube .ge. E0) .and. ANY(Cube .le. E0)) THEN
           CALL CubeCut(E0,Cube,CubEdge,NcEdge,Ut,dt,dx,dz,Ambig,Pathological)
           IF (Pathological) RETURN
           CALL StoreCube(Edge,Nedge,Out,CubEdge,NcEdge,Ut,i*dy,2,(/0,1,3/))
        END IF
     END DO

     DO i = 0,1 ! z = 0 or 1
        Cube = HyperCube(:,:,:,i)
        IF (ANY(Cube .ge. E0) .and. ANY(Cube .le. E0)) THEN
           CALL CubeCut(E0,Cube,CubEdge,NcEdge,Ut,dt,dx,dy,Ambig,Pathological)
           IF (Pathological) RETURN
           CALL StoreCube(Edge,Nedge,Out,CubEdge,NcEdge,Ut,i*dz,3,(/0,1,2/))
        END IF
     END DO

   END SUBROUTINE Cubes


   SUBROUTINE CubeCut(E0,Cube,Edge,Nedge,Ut,dt,dx,dy,Ambig,Patho)

     IMPLICIT NONE
     REAL(KIND(0D0)) :: E0,dx,dy,dt
     REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1) :: Cube
     REAL(KIND(0D0)),DIMENSION(0:2,3,12)    :: Edge
     REAL(KIND(0D0)),DIMENSION(0:2,12)      :: Ut
     INTEGER :: Nedge
     LOGICAL :: Ambig, Patho

     REAL(KIND(0D0)),DIMENSION(0:2) :: Vmid
     INTEGER,DIMENSION(5) :: EdgeSet
     INTEGER :: i,j, Nsurfs
     LOGICAL :: Outo

     Outo = .false.
     Nedge = 0
     Ut = 0D0

     CALL Edges(E0,Cube,Edge,Ut,dt,dx,dy,Nedge,Outo,Patho)
     IF (Patho) RETURN

     IF ((.not.(Outo)).and.(Nedge .eq. 6))                           &
          Outo = (((Cube(0,0,0)-E0)*(Cube(1,1,1)-E0).gt.0D0).and.    &
                  ((Cube(1,0,0)-E0)*(Cube(0,1,1)-E0).gt.0D0).and.    &
                  ((Cube(1,1,0)-E0)*(Cube(0,0,1)-E0).gt.0D0).and.    &
                  ((Cube(0,1,0)-E0)*(Cube(1,0,1)-E0).gt.0D0))
     IF (Outo) THEN
        CALL Disconnected2D(Nedge,Edge,Ut,Nsurfs,EdgeSet)
        Ambig = ((Ambig).or.(Nsurfs .gt. 1))
     ELSE
        Nsurfs = 1
        EdgeSet(1) = 1
        EdgeSet(2) = Nedge+1
     END IF

! Calculate the average of intersection points

     DO j = 1,Nsurfs
        Vmid = 0D0
        DO i = Edgeset(j),Edgeset(j+1)-1
           Vmid = Vmid + Edge(:,1,i) + Edge(:,2,i)
        END DO
        Vmid = Vmid/(2*(Edgeset(j+1)-Edgeset(j)))
        DO i = Edgeset(j),Edgeset(j+1)-1
           Edge(:,3,i) = Vmid
        END DO
     END DO

   END SUBROUTINE CubeCut


   SUBROUTINE Edges(E0,Cube,Edge,Ut,dt,dx,dy,Nedge,Outo,Patho)

     IMPLICIT NONE
     REAL(KIND(0D0)) :: E0,dx,dy,dt
     REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1) :: Cube
     REAL(KIND(0D0)),DIMENSION(0:1,0:1)  :: Square
     REAL(KIND(0D0)),DIMENSION(0:2,3,12) :: Edge
     REAL(KIND(0D0)),DIMENSION(0:2,12)   :: Ut
     REAL(KIND(0D0)),DIMENSION(2,4)      :: Cut
     REAL(KIND(0D0)),DIMENSION(2,2)      :: Out
     INTEGER         :: Ncuts,Nedge, i
     LOGICAL         :: Outo, Patho

     DO i = 0,1  ! t = 0 or 1
        Square = Cube(i,:,:)
        CALL FindEdge(E0,Square,Cut,Out,dx,dy,Ncuts,Patho)
        IF (Patho) RETURN
        CALL StoreEdge(Ncuts,Cut,Out,i*dt,0,(/1,2/),Nedge,Edge,Ut,Outo)
     END DO

     DO i = 0,1  ! x = 0 or 1
        Square = Cube(:,i,:)
        CALL FindEdge(E0,Square,Cut,Out,dt,dy,Ncuts,Patho)
        IF (Patho) RETURN
        CALL StoreEdge(Ncuts,Cut,Out,i*dx,1,(/0,2/),Nedge,Edge,Ut,Outo)
     END DO

     DO i = 0,1  ! y = 0 or 1
        Square = Cube(:,:,i)
        CALL FindEdge(E0,Square,Cut,Out,dt,dx,Ncuts,Patho)
        IF (Patho) RETURN
        CALL StoreEdge(Ncuts,Cut,Out,i*dy,2,(/0,1/),Nedge,Edge,Ut,Outo)
     END DO

   END SUBROUTINE Edges


   SUBROUTINE FindEdge(E0,Square,Cut,Out,dx,dy,Ncuts,Patho)

     IMPLICIT NONE

     REAL(KIND(0D0)) :: E0,dx,dy
     REAL(KIND(0D0)),DIMENSION(0:1,0:1) :: Square
     REAL(KIND(0D0)),DIMENSION(2,4)     :: Cut
     REAL(KIND(0D0)),DIMENSION(2,2)     :: Out
     INTEGER :: Ncuts
     LOGICAL :: Patho

     CALL EndsOfEdge(E0,Square,Cut,dx,dy,Ncuts)

     IF (Ncuts .gt. 0) CALL FindOutside(Ncuts,E0,Square,Cut,Out,dx,dy)

     IF ((Ncuts .eq. 3).or.(Ncuts .eq. 1)) THEN
        WRITE(*,*) 'Error in FindEdge, too many (i.e.',Ncuts,') cuts.'
        WRITE(*,*) 'Noncontinuous surface. E0 =',E0
        WRITE(*,*) 'Eps(0,0) =',Square(0,0),' Eps(1,0) =',Square(1,0)
        WRITE(*,*) 'Eps(0,1) =',Square(0,1),' Eps(1,1) =',Square(1,1)
        WRITE(*,*)
        Patho = .true.
     END IF

   END SUBROUTINE FindEdge


   SUBROUTINE EndsOfEdge(E0,Square,Cut,dx,dy,Ncuts)     

     IMPLICIT NONE
     REAL(KIND(0D0)) :: E0,dx,dy
     REAL(KIND(0D0)),DIMENSION(0:1,0:1) :: Square
     REAL(KIND(0D0)),DIMENSION(2,4)     :: Cut
     INTEGER :: Ncuts

     Ncuts = 0
     IF (((Square(0,0)-E0)*(Square(1,0)-E0)) .lt. 0D0) THEN
        Ncuts = Ncuts+1
        Cut(1,Ncuts) = (Square(0,0)-E0)/(Square(0,0)-Square(1,0))*dx
        Cut(2,Ncuts) = 0D0
     ELSE
        IF ((Square(0,0).eq.E0) .or. (Square(1,0).eq.E0))            &
             CALL EndsAtCorner(Square(0,0),Square(1,0),E0,Ncuts,     &
                               Cut(1,Ncuts+1),Cut(2,Ncuts+1),dx,0D0)
     END IF

     IF (((Square(0,0)-E0)*(Square(0,1)-E0)) .lt. 0D0) THEN
        Ncuts = Ncuts + 1
        Cut(1,Ncuts) = 0D0
        Cut(2,Ncuts) = (Square(0,0)-E0)/(Square(0,0)-Square(0,1))*dy
     ELSE
        IF ((Square(0,0).eq.E0) .or. (Square(0,1).eq.E0))            &
             CALL EndsAtCorner(Square(0,0),Square(0,1),E0,Ncuts,     &
                               Cut(2,Ncuts+1),Cut(1,Ncuts+1),dy,0D0)
     END IF

     IF (((Square(1,0)-E0)*(Square(1,1)-E0)) .lt. 0D0) THEN
        Ncuts = Ncuts+1
        Cut(1,Ncuts) = dx
        Cut(2,Ncuts) = (Square(1,0)-E0)/(Square(1,0)-Square(1,1))*dy
     ELSE
        IF ((Square(1,0).eq.E0) .or. (Square(1,1).eq.E0))            &
             CALL EndsAtCorner(Square(1,0),Square(1,1),E0,Ncuts,     &
                               Cut(2,Ncuts+1),Cut(1,Ncuts+1),dy,dx)
     END IF

     IF (((Square(0,1)-E0)*(Square(1,1)-E0)) .lt. 0D0) THEN
        Ncuts = Ncuts+1
        Cut(1,Ncuts) = (Square(0,1)-E0)/(Square(0,1)-Square(1,1))*dx
        Cut(2,Ncuts) = dy
      ELSE
        IF ((Square(0,1).eq.E0) .or. (Square(1,1).eq.E0))            &
             CALL EndsAtCorner(Square(0,1),Square(1,1),E0,Ncuts,     &
                               Cut(1,Ncuts+1),Cut(2,Ncuts+1),dx,dy)
    END IF

   END SUBROUTINE EndsOfEdge


   SUBROUTINE EndsAtCorner(A,B,E0,Ncuts,C1,C2,d1,d2)

     IMPLICIT NONE
     REAL(KIND(0D0)) :: A,B,E0,C1,C2,d1,d2
     INTEGER :: Ncuts

     IF ((A .eq. E0).and.(B .lt. E0)) THEN
        Ncuts = Ncuts+1
        C1 = 1D-9*d1
        C2 = d2
     END IF
     IF ((A .lt. E0).and.(B .eq. E0)) THEN
        Ncuts = Ncuts+1
        C1 = (1D0-1D-9)*d1
        C2 = d2
     END IF

   END SUBROUTINE EndsAtCorner


!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE FindOutside(Ncuts,E0,Square,Cut,Out,dx,dy)

! Finds a point outside the freeze-out surface and sorts ambiguous surfaces
! with two edges on one face of the cube. The rule is to interpolate the value
! at the center of the face, see if it is below or above freeze-out criterion
! and set the surface accordingly.

     IMPLICIT NONE
     INTEGER :: Ncuts
     REAL(KIND(0D0)) :: E0,dx,dy
     REAL(KIND(0D0)),DIMENSION(0:1,0:1) :: Square
     REAL(KIND(0D0)),DIMENSION(2,4)     :: Cut
     REAL(KIND(0D0)),DIMENSION(2,2)     :: Out

     INTEGER :: i, j, Nout
     REAL(KIND(0D0)) :: Eave

     IF (Ncuts .eq. 4) THEN  ! Ambiguous surface check, interpolate the center
        Eave = 2.5D-1*SUM(Square)
!        IF ((Square(0,0)-E0)*(Eave-E0) .ge. 0D0) THEN
! Improved July 26, 2011 to handle the corner at E0 case better
!        IF (    ((Square(0,0).ne.E0).and.((Square(0,0)-E0)*(Eave-E0).ge.0D0)) &
!            .or.((Square(0,0).eq.E0).and.(Eave.ge.E0)) ) THEN
! Improved Aug 6, 2011 to handle the center at E0 case better
        IF (    ((Square(0,0).lt.E0).and.(Eave.lt.E0)) &
            .or.((Square(0,0).ge.E0).and.(Eave.ge.E0))) THEN
           Out(:,1) = Cut(:,2)
           Cut(:,2) = Cut(:,3)
           Cut(:,3) = Out(:,1)
        END IF
        IF ((Eave-E0) .lt. 0D0) THEN ! Outward direction at the center
           Out(1,:) = 5D-1*dx
           Out(2,:) = 5D-1*dy
        ELSE
           IF ((Square(0,0)-E0) .lt. 0D0) THEN
              Out(:,1) = 0D0
              Out(1,2) = dx
              Out(2,2) = dy
           ELSE
              Out(1,1) = dx
              Out(2,1) = 0D0
              Out(1,2) = 0D0
              Out(2,2) = dy
           END IF
        END IF
     ELSE         ! Normal case, only one edge cutting the face of the cube
        Out = 0D0 ! Find the direction outwards (to lower value)
        Nout = 0
        DO i = 0,1
           DO j = 0,1
              IF (Square(i,j) .lt. E0) THEN
                 Out(1,1) = Out(1,1) + i*dx
                 Out(2,1) = Out(2,1) + j*dy
                 Nout = Nout + 1
              END IF
           END DO
        END DO
        IF (Nout .gt. 0) Out = Out/Nout
     END IF

   END SUBROUTINE FindOutside


   SUBROUTINE StoreEdge(Ncuts,Cut,Out,dKnown,Trivial,nonTrivial,       &
                        Nedge,Edge,Ut,Outo)
     IMPLICIT NONE

     LOGICAL              :: Outo
     INTEGER              :: Ncuts, Nedge, Trivial
     INTEGER,DIMENSION(2) :: nonTrivial
     REAL(KIND(0D0)),DIMENSION(2,4) :: Cut
     REAL(KIND(0D0)),DIMENSION(2,2) :: Out
     REAL(KIND(0D0))                :: dKnown
     REAL(KIND(0D0)),DIMENSION(0:2,3,12) :: Edge
     REAL(KIND(0D0)),DIMENSION(0:2,12)   :: Ut

     IF ((Ncuts .eq. 2).or.(Ncuts .eq. 4)) THEN
        Nedge = Nedge + 1
        Edge(Trivial,1:2,Nedge) = dKnown
        Edge(nonTrivial,1:2,Nedge) = Cut(:,1:2)
        Ut(Trivial,Nedge) = dKnown
        Ut(nonTrivial,Nedge) = Out(:,1)
     END IF
     IF (Ncuts .eq. 4) THEN
        Nedge = Nedge + 1
        Edge(Trivial,1:2,Nedge) = dKnown
        Edge(nonTrivial,1:2,Nedge) = Cut(:,3:4)
        Ut(Trivial,Nedge) = dKnown
        Ut(nonTrivial,Nedge) = Out(:,2)
        Outo = .true.
     END IF

   END SUBROUTINE StoreEdge


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE Disconnected2D(N,Edge,Ut,Nsurfs,EdgeSet)

! Subroutine to check whether the surface cuts the cube once or several times

     IMPLICIT NONE
     INTEGER :: N, Nsurfs
     INTEGER,DIMENSION(5) :: EdgeSet
     REAL(KIND(0D0)),DIMENSION(0:2,3,12) :: Edge
     REAL(KIND(0D0)),DIMENSION(0:2,2,12) :: Side
     REAL(KIND(0D0)),DIMENSION(0:2,12)   :: Ut, Ut2
     INTEGER :: i, j, Nside, Nedge

     Nsurfs = 1
     EdgeSet(1) = 1
     Side(:,:,1) = Edge(:,1:2,N)
     Ut2(:,1) = Ut(:,N)
     Nside = 1
     Nedge = N-1

     DO WHILE (Nedge .gt. 0)
        i = 1
        DO WHILE (ANY(Side(:,1,Nside).ne.Edge(:,1,i)).and.               &
                  ANY(Side(:,1,Nside).ne.Edge(:,2,i)).and.(i.le.Nedge))
           i = i+1
        END DO
        IF (i.le.Nedge) THEN
           IF (ALL(Side(:,1,Nside).eq.Edge(:,1,i))) THEN
              Side(:,2,Nside+1) = Edge(:,1,i)
              Side(:,1,Nside+1) = Edge(:,2,i)
           ELSE
              Side(:,:,Nside+1) = Edge(:,1:2,i)
           END IF
           Ut2(:,Nside+1) = Ut(:,i)
           Edge(:,1:2,i) = Edge(:,1:2,Nedge)
           Ut(:,i) = Ut(:,Nedge)
           EdgeSet(Nsurfs+1) = Nside+2
        ELSE
           Nsurfs = Nsurfs+1
           Side(:,:,Nside+1) = Edge(:,1:2,Nedge)
           Ut2(:,Nside+1) = Ut(:,Nedge)
        END IF
        Nside = Nside+1
        Nedge = Nedge-1
     END DO
     Edge(:,1:2,1:N) = Side(:,:,1:N)
     Ut(:,1:N) = Ut2(:,1:N)

   END SUBROUTINE Disconnected2D


!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE StoreCube(Edge,Nedge,Out,CubEdge,NcEdge,Ut,dKnown,Trivial,  &
                        nonTrivial)

     IMPLICIT NONE
     REAL(KIND(0D0)),DIMENSION(0:3,3,96) :: Edge
     REAL(KIND(0D0)),DIMENSION(0:2,3,12) :: CubEdge
     REAL(KIND(0D0)),DIMENSION(0:3,96)   :: Out
     REAL(KIND(0D0)),DIMENSION(0:2,12)   :: Ut
     REAL(KIND(0D0))      :: dKnown
     INTEGER              :: Trivial,Nedge, NcEdge, i
     INTEGER,DIMENSION(3) :: nonTrivial

     DO i = 1,NcEdge
        Nedge = Nedge + 1
        Edge(Trivial,:,Nedge) = dKnown
        Edge(nonTrivial,:,Nedge) = CubEdge(:,:,i)
        Out(Trivial,Nedge) = dKnown
        Out(nonTrivial,Nedge) = Ut(:,i)
     END DO

   END SUBROUTINE StoreCube


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   FUNCTION Suspicious(Nedge,Ngp)

     IMPLICIT NONE
     LOGICAL :: Suspicious
     INTEGER :: Nedge, Ngp, N

     IF (Ngp .gt. 8) THEN 
        N = 16-Ngp
     ELSE
        N = Ngp
     END IF

     Suspicious = ((Nedge .ge. 42).or.                                   &
                   ((Nedge .ge. 36).and.((N .eq. 5).or.(N .eq. 4))).or.  &
                   ((Nedge .ge. 30).and.(N .eq. 3)).or.                  &
                   ((Nedge .ge. 24).and.(N .eq. 2)))

   END FUNCTION Suspicious


   SUBROUTINE Disconnected3D(N,Edge,Ut,Nsurfs,EdgeSet)

! Subroutine to check whether the hypersurface cuts the hypercube 
! once or several times, i.e. is the hypersurface connected or disconnected

     IMPLICIT NONE
     INTEGER :: N, Nsurfs
     REAL(KIND(0D0)),DIMENSION(0:3,3,96) :: Edge, Edge2
     REAL(KIND(0D0)),DIMENSION(0:3,96) :: Ut,Ut2
     REAL(KIND(0D0)),DIMENSION(0:3) :: Point
     INTEGER,DIMENSION(9) :: EdgeSet
     INTEGER          :: Nedge, n1, n2, n3
     LOGICAL,EXTERNAL :: SameEdgeDifferentPolygon
     LOGICAL,EXTERNAL :: NextEdgeSamePolygon

     Edge2(:,:,1:N) = Edge(:,:,1:N)
     Ut2(:,1:N) = Ut(:,1:N)
     n2 = 0; n3 = 1
     Edge(:,:,1) = Edge2(:,:,N)
     Ut(:,1) = Ut2(:,N)
     EdgeSet(1) = 1
     Nedge = N-1
     Nsurfs = 1

     DO WHILE (Nedge .gt. 0)
        n1 = n2+1
        n2 = n3
        CALL FindNextEdges(SameEdgeDifferentPolygon,n1,n2,n3,        &
                           Edge,Edge2,Ut,Ut2,EdgeSet,Nedge,Nsurfs)
        n2 = n3
        CALL FindNextEdges(NextEdgeSamePolygon,n1,n2,n3,             &
                           Edge,Edge2,Ut,Ut2,EdgeSet,Nedge,Nsurfs)
        IF ((n2 .eq. n3).and.(Nedge .gt. 0)) THEN
           Nsurfs = Nsurfs+1
           n3 = n3+1
           Edge(:,:,n3) = Edge2(:,:,Nedge)
           Ut(:,n3) = Ut2(:,Nedge)
           Nedge = Nedge-1
        END IF
     END DO

   END SUBROUTINE Disconnected3D


   SUBROUTINE FindNextEdges(Connects,n1,n2,n3,Edge,Edge2,Ut,Ut2,EdgeSet, &
                            Nedge,Nsurfs)
     IMPLICIT NONE
     INTEGER :: n1, n2, n3, Nedge, Nsurfs
     REAL(KIND(0D0)),DIMENSION(0:3,3,96) :: Edge, Edge2
     REAL(KIND(0D0)),DIMENSION(0:3,96) :: Ut,Ut2
     INTEGER,DIMENSION(9) :: EdgeSet
     LOGICAL,EXTERNAL :: Connects
     INTEGER :: i,j

     DO j = n1,n2
        i = 1
        DO WHILE ((.not.Connects(j,i,Edge,Edge2)).and.       &
                  (i .le. Nedge))
           i = i+1
        END DO
        IF (i .le. Nedge) THEN
           n3 = n3+1
           IF (ALL(Edge(:,1,j).eq.Edge2(:,1,i))) THEN
              Edge(:,1,n3) = Edge2(:,2,i)
              Edge(:,2,n3) = Edge2(:,1,i)
              Edge(:,3,n3) = Edge2(:,3,i)
           ELSE
              Edge(:,:,n3) = Edge2(:,:,i)
           END IF
           Ut(:,n3) = Ut2(:,i)
           Edge2(:,:,i) = Edge2(:,:,Nedge)
           Ut2(:,i) = Ut2(:,Nedge)
           EdgeSet(Nsurfs+1) = n3+1
           Nedge = Nedge-1
        END IF
     END DO

   END SUBROUTINE FindNextEdges


   FUNCTION SameEdgeDifferentPolygon(j,i,Edge,Edge2) RESULT (SedP)

     IMPLICIT NONE
     LOGICAL :: SedP
     INTEGER :: j,i
     REAL(KIND(0D0)),DIMENSION(0:3,3,96) :: Edge, Edge2

     SedP = ALL(Edge(:,1:2,j).eq.Edge2(:,1:2,i)).or.     &
            (ALL(Edge(:,1,j).eq.Edge2(:,2,i)).and.         &
             ALL(Edge(:,2,j).eq.Edge2(:,1,i)))

   END FUNCTION SameEdgeDifferentPolygon

   FUNCTION NextEdgeSamePolygon(j,i,Edge,Edge2) RESULT (NesP)

     IMPLICIT NONE
     LOGICAL :: NesP
     INTEGER :: j,i
     REAL(KIND(0D0)),DIMENSION(0:3,3,96) :: Edge, Edge2

     NesP = ((ALL(Edge(:,1,j).eq.Edge2(:,1,i)).or.       &
              ALL(Edge(:,1,j).eq.Edge2(:,2,i))).and.       &
             ALL(Edge(:,3,j).eq.Edge2(:,3,i)))

   END FUNCTION NextEdgeSamePolygon


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE NormalVector(ns,First,Last,Nedge,Edge,Out,dSigma,Vmids)

     IMPLICIT NONE
     INTEGER :: ns, First, Last, Nedge, i
     REAL(KIND(0D0)),DIMENSION(0:3,3,96) :: Edge    ! Table for ends of edges  
     REAL(KIND(0D0)),DIMENSION(0:3,96) :: Out       ! outside direction
     REAL(KIND(0D0)),DIMENSION(0:3,8) :: dSigma, Vmids

     REAL(KIND(0D0)),DIMENSION(0:3) :: V,Vmid,Vout, A,B,C, SuL
     REAL(KIND(0D0)) :: Volume, V_i

! Calculate the mean vector of intersection points

     V = 0D0
     DO i = First,Last
        V = V + Edge(:,1,i) + Edge(:,2,i)
     END DO
     Vmid = V/(2*Nedge)

! Calculate the centroid, i.e., the center of gravity of the surface element
! (if Nedge=12, the surface element is a tetrahedron and centroid is the mean
!  of corner coordinates)

     IF (Nedge .eq. 12) THEN
        Vmids(:,ns) = Vmid
     ELSE
        V = 0D0
        Volume = 0D0
        DO i = First,Last
           A = Edge(:,1,i) - Vmid
           B = Edge(:,2,i) - Vmid
           C = Edge(:,3,i) - Vmid
           CALL TetraVolume(A,B,C,SuL)
           V_i = Sqrt(SuL(0)**2+SuL(1)**2+SuL(2)**2+SuL(3)**2)
           V = V + V_i*(Edge(:,1,i)+Edge(:,2,i)+Edge(:,3,i)+Vmid)*2.5D-1
           Volume = Volume + V_i
        END DO
        Vmids(:,ns) = V/Volume
     END IF

     dSigma(:,ns) = 0D0
     DO i = First,Last
        A = Edge(:,1,i) - Vmids(:,ns)
        B = Edge(:,2,i) - Vmids(:,ns)
        C = Edge(:,3,i) - Vmids(:,ns)

        CALL TetraVolume(A,B,C,SuL)

        Vout = Out(:,i) - Vmids(:,ns)
        dSigma(:,ns) = dSigma(:,ns) + SIGN(1D0,Dot_Product(Vout,SuL))*SuL
     END DO

   END SUBROUTINE NormalVector


   SUBROUTINE TetraVolume(A,B,C,SuL)

     IMPLICIT NONE
     REAL(KIND(0D0)),DIMENSION(0:3) :: A,B,C,SuL

     SuL(0) =  1/6D0*( A(1)*(B(2)*C(3)-B(3)*C(2)) &
                      -A(2)*(B(1)*C(3)-B(3)*C(1)) &
                      +A(3)*(B(1)*C(2)-B(2)*C(1)) )
     SuL(1) = -1/6D0*( A(0)*(B(2)*C(3)-B(3)*C(2)) &
                      -A(2)*(B(0)*C(3)-B(3)*C(0)) &
                      +A(3)*(B(0)*C(2)-B(2)*C(0)) )
     SuL(2) =  1/6D0*( A(0)*(B(1)*C(3)-B(3)*C(1)) &
                      -A(1)*(B(0)*C(3)-B(3)*C(0)) &
                      +A(3)*(B(0)*C(1)-B(1)*C(0)) )
     SuL(3) = -1/6D0*( A(0)*(B(1)*C(2)-B(2)*C(1)) &
                      -A(1)*(B(0)*C(2)-B(2)*C(0)) &
                      +A(2)*(B(0)*C(1)-B(1)*C(0)) )

   END SUBROUTINE TetraVolume



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE DeadEnd(Cube,E0)

     IMPLICIT NONE
     REAL(KIND(0D0)) :: E0
     REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1,0:1) :: Cube

     WRITE(*,*) 'Unknown surface structure.'
     WRITE(*,*) 'Probably too many gridpoints exactly in freeze-out density.'
     WRITE(*,*) 'Freeze-out value:',E0
     WRITE(*,*) 'Values at HyperCube corners:'
     WRITE(*,*) 'E(0,0,1,0) =',Cube(0,0,1,0),'E(0,1,1,0) =',Cube(0,1,1,0)
     WRITE(*,*) 'E(0,0,0,0) =',Cube(0,0,0,0),'E(0,1,0,0) =',Cube(0,1,0,0)
     WRITE(*,*)
     WRITE(*,*) 'E(0,0,1,1) =',Cube(0,0,1,1),'E(0,1,1,1) =',Cube(0,1,1,1)
     WRITE(*,*) 'E(0,0,0,1) =',Cube(0,0,0,1),'E(0,1,0,1) =',Cube(0,1,0,1)
     WRITE(*,*)
     WRITE(*,*) 'E(1,0,1,0) =',Cube(1,0,1,0),'E(1,1,1,0) =',Cube(1,1,1,0)
     WRITE(*,*) 'E(1,0,0,0) =',Cube(1,0,0,0),'E(1,1,0,0) =',Cube(1,1,0,0)
     WRITE(*,*)
     WRITE(*,*) 'E(1,0,1,1) =',Cube(1,0,1,1),'E(1,1,1,1) =',Cube(1,1,1,1)
     WRITE(*,*) 'E(1,0,0,1) =',Cube(1,0,0,1),'E(1,1,0,1) =',Cube(1,1,0,1)
     STOP

   END SUBROUTINE DeadEnd


   SUBROUTINE StrangeN(Ngp)

     IMPLICIT NONE
     INTEGER :: Ngp

     WRITE(*,*) 'According to Ngp',Ngp,' of the hypercube corners are', &
                ' within the freeze-out surface.'
     WRITE(*,*) 'Strange. Check the value.'
     STOP

   END SUBROUTINE StrangeN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
