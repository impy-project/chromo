   FUNCTION Quadrilinear(Vmid,HyperCube,dt,dx,dy,dz) RESULT (F)

     IMPLICIT NONE
     REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1,0:1) :: HyperCube
     REAL(KIND(0D0)),DIMENSION(0:3) :: Vmid
     REAL(KIND(0D0)) :: dt, dx, dy, dz, F
     REAL(KIND(0D0)) :: t, x, y, z
     INTEGER :: n,i,j,k

     F = 0D0
     t = Vmid(0)/dt
     x = Vmid(1)/dx
     y = Vmid(2)/dy
     z = Vmid(3)/dz

     DO k = 0,1
        DO j = 0,1
           DO i = 0,1
              DO n = 0,1
                 F = F + (n*t+(1-n)*(1-t))*(i*x+(1-i)*(1-x))     &
                        *(j*y+(1-j)*(1-y))*(k*z+(1-k)*(1-z))     &
                        *HyperCube(n,i,j,k)

              END DO
           END DO
        END DO
     END DO

   END FUNCTION Quadrilinear
