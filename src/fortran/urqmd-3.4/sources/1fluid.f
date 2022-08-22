CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Modified 1fhydro code for coupling to UrQMD                        C
C  Changes by H.Petersen and M. Bleicher, 22.11.2007                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c*1fv6.f***********************************D.H.Rischke,September 1996**
c     q                                                                *
c    This program calculates heavy-ion collisions on a 3d grid.       *
c                                                                     *
c    The initial fluid distributions are specified as two nuclei      *
c    colliding at impact parameter b, or read in from a table.        *
c    Optionally, one can initialize a single nucleus.                 *
c                                                                     *
c    The coordinate system is:                                        *
c               z - beam axis,                                        *
c               x - transverse to beam, in reaction plane,            *
c               y - transverse to beam, out of reaction plane.        *
c                                                                     *
c    The underlying hydrodynamical transport algorithm is the         *
c    phoenical SHASTA with (optional) half step and simplified        *
c    source treatment.                                                *
c                                                                     *
c    There is propagation of (baryon) density.                        *
c                                                                     *
c    The maximum grid size is (ngr,ngr,ngr) cells.                    *
c                                                                     *
c    The equation of state is an ultrarelativistic ideal gas eos,     *
c    or is read in from a table.                                      *
c                                                                     *
c    The output requires paw++.                                       *
c                                                                     *
c    The difference to program 1f.f is:                               *
c           -  a reduction of loops in subroutine prop3d,             *
c           -  monitoring conservation laws in each time step,        *
c           -  binning of fluid dN/dy and p_x(y),                     *
c           -  calculation of single inclusive spectra,               *
c                             baryonic p_x(y),                        *
c                             pionic p_x(y)                           *       
c              for freeze-out at constant time.                       *
c                                                                     *
c    The difference to program 1fv2.f is:                             *
c           -  single nucleus initialization option,                  *
c           -  updated comments,                                      *
c           -  variable vacuum cut-off,                               *
c           -  corrected freeze-out in subprogram fileo,              *
c           -  temperature profiles in the reaction plane,            *
c           -  lab energy density in the transverse plane at z=0,     *
c           -  an option to stabilize a cold nucleus,                 *
c           -  mean temperature and (comoving) baryon density         *
c              calculation,                                           *
c           -  calculation of directed transverse momenta,            *
c                                                                     *
c    The difference to program 1fv3.f is:                             *
c           -  calculation of p_x(y).                                 *
c           -  measuring the actual CPU time of calculation,          *
c           -  calculation of the fraction of QGP per cell.           *
c                                                                     *
c    The difference to program 1fv5.f is:                             *
c           -  read initial conditions from a table, write final      *
c              situation on a table.                                  *
c                                                                     *
c**********************************************************************

      subroutine untang(elab,mx,my,mz,r,e,p,n,vx,vy,vz,gamma)
      INCLUDE 'defs.f'
c
c
c
c  This subroutine-subprogram calculates local rest frame quantities
c  from laboratory frame quantities.
c  It is used in subroutine prop3d.
c
c  Type declarations for variables in common-blocks.
c
      real*8 dx,dt,t,vol,cut
c
c  Type declarations for variables in subroutine untang.
c
      real*8 e(ngr),p(ngr),n(ngr)
      real*8 vx(ngr),vy(ngr),vz(ngr),gamma(ngr)
      real*8 elab(ngr),mx(ngr),my(ngr),mz(ngr),r(ngr)
      real*8 v(ngr),m(ngr)
      real*8 velo,press
      integer i
c
c     Common-blocks.
c
      common /gitter/ dx,dt,t,vol,cut
c
c     Start calculation. Each cell is separately treated.
c     Vacuum is assumed if the absolute amount of the 
c     calculational frame field elab is smaller than cut.
c


      do 100 i = 1,ngr
         m(i) = dsqrt(mx(i)*mx(i) + my(i)*my(i) + mz(i)*mz(i))
         if (dabs(elab(i)).le.cut) then
            e(i) = 0d0
            vx(i) = 0d0
            vy(i) = 0d0
            vz(i) = 0d0
            n(i) = 0d0
            gamma(i) = 1d0
         else if (m(i).lt.1d-16) then
            e(i) = elab(i)
            vx(i) = 0d0
            vy(i) = 0d0
            vz(i) = 0d0
            n(i) = r(i)
            gamma(i) = 1d0
         else 
c     
c     The velocity iteration works for the absolute amount
c     of the velocity and the lab momentum.
c     

            e(i) = elab(i)
            vx(i) = 0d0
            vy(i) = 0d0
            vz(i) = 0d0
            gamma(i) = 1d0
            n(i) = r(i)

            v(i) = velo(elab(i),m(i),r(i))
            if (v(i).lt.1d0-1d-16) then
               e(i) = elab(i) - v(i)*m(i)
               vx(i) = v(i)*mx(i)/m(i)
               vy(i) = v(i)*my(i)/m(i)
               vz(i) = v(i)*mz(i)/m(i)
               gamma(i) = 1d0/dsqrt(1d0-v(i)**2)
               n(i) = r(i)/gamma(i)
c     
c     Vacuum is assumed if v .ge. 1 - 10^{-16}.
c     
            else
               e(i) = 0d0
               vx(i) = 0d0
               vy(i) = 0d0
               vz(i) = 0d0
               n(i) = 0d0
               gamma(i) = 1d0 
            end if
         end if
         p(i) = press( e(i), n(i) )
 100  continue
      return
      end
c
c----------------------------------------------------------------------
c
      function velo(el,ml,rl)
      INCLUDE 'defs.f'
c
c  This function-subprogram calculates the velocity as a fixed
c  point of v =  M / ( Elab + p ) .
c  It is used in the subroutines untang and tinit.
c  el,ml,rl are the values of Elab,M,R.
c
c  Type declarations for variables in function velo.
c
      real*8 velo,el,ml,rl
      real*8 v,f,press,root
      integer i
c
c  Fixed-point iteration of the function v = f(v), with
c  f(v) = M / ( Elab + p( Elab-v*M, R*sqrt(1-v*v) ) ).
c  Relative accuracy of the first level of iteration is 1e-10.
c           
      
      i = 0
      f = 0d0
 10   v = f
      if (v.le.1d0) then
         root = dsqrt( 1d0 - v*v )
      else
         root = 0d0
      end if
      f = ml / ( el + press(el-v*ml,rl*root) )
      i = i + 1
c     
c     If more than 1000 iterations: iterate to find v only
c     on the 1 percent relative accuracy level. If that still doesn't
c     work, set v to zero by hand and return.
c     
      if (i.gt.1000) then
         i = 0
         f = 0d0
 20      v = f
         if (v.le.1d0) then
            root = dsqrt( 1d0 - v*v )
         else
            root = 0d0
         end if
         f = ml / ( el + press(el-v*ml,rl*root) )
         i = i + 1
         if (i.gt.3000) then
            write(6,*) ' No root found for Elab=',el,', M=',ml,', R=',rl
c            write(6,*) ' v set to zero by hand.'
             if(f.le.1.0d0)then
             v=(v+f)/2.0
             else
             v = 0d0
             endif
             write(6,*) ' v set by hand: ', v
            return
         end if
         if (dabs((v-f)/f).gt.1d-2) goto 20
         velo = f
         return
      end if
      if (dabs((v-f)/f).gt.1d-7) goto 10
      velo = f
      return
      

      end
c
c----------------------------------------------------------------------
c
      function press(e,n)
      INCLUDE 'defs.f'
c
c  This function-subprogram determines the pressure of the
c  underlying EoS.
c  It is used in subroutines sinit, tinit, untang, velo.
c
c  Type declarations for variables in common-blocks.
c
      real*8 e0,n0,B,hc3,pi2
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)


      integer eos,stabil,anti
c
c  Type declarations for variables used in function press.
c
      real*8 e,n,press,et0
      real*8 de,dn,p1,p2,p3,p4,p13,p24
c      real*8 p12,p34
c
c  Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
      common /eqofst/ eos,stabil,anti
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3


c
c  Fermi energy density in the QGP for given baryon density
c  (in units of e0).
c
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*n*n0*hc3)**(4d0/3d0)+B
      et0 = et0/e0/hc3

      anti = 0


      if (n.lt.0.0)then
         n=-n
         anti = 1
      end if

c
c     Calculate pressure (in units of e0).
c     If flag eos=[else], use tabellized equation of state.
c     
c     
c     arithmetic mean interpolation (yields better values for
c     interpolation at small e,n, but is worse for all other
c     thermodynamic quantities, also at larger e,n)
c     
c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     
      if (eos.eq.3) then
         if (e.le.20d0) then
            de = 0.1d0
            dn = 0.05d0
            p1 = ptab(idint(e/de),idint(n/dn))
            p2 = ptab(idint(e/de)+1,idint(n/dn))
            p3 = ptab(idint(e/de),idint(n/dn)+1)
            p4 = ptab(idint(e/de)+1,idint(n/dn)+1)
            p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
            p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
            p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
            p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
      
            press = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     press = 0.25d0*(p13+p24+p12+p34)
         else if (e.ge.20d0) then
            press = (e-4d0*B/e0/hc3)/3d0
         end if
         
      else if ((eos.eq.4).or.(eos.eq.2).or.(eos.eq.5)) then
         if (e.le.1000.0d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = ptab3(idint(e/de),idint(n/dn))
               p2 = ptab3(idint(e/de)+1,idint(n/dn))
               p3 = ptab3(idint(e/de),idint(n/dn)+1)
               p4 = ptab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2 
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               press = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     write(99,*) p1, ptab3(0,0) 
            end if
            
            if (((e.lt.10.0d0).and.(n.lt.2.0d0)).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = ptab2(idint(e/de),idint(n/dn))
               p2 = ptab2(idint(e/de)+1,idint(n/dn))
               p3 = ptab2(idint(e/de),idint(n/dn)+1)
               p4 = ptab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2 
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               press = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then

               de = 0.5d0
               dn = 0.1d0
 
               p1 = ptab(idint(e/de),idint(n/dn))
               p2 = ptab(idint(e/de)+1,idint(n/dn))
               p3 = ptab(idint(e/de),idint(n/dn)+1)
               p4 = ptab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2 
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               press = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13  
            end if
         else if (e.ge.1000.0d0) then
c     JS muss bestimmt werden
            press = e/3.0d0
         end if
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1).and.(n.ge.-1d0))
     $     then
         press = 0d0
      end if
      
      if (press.lt.0d0) press=0d0
      
      
      
      
      
      if (anti.eq.1)then
         n=-n
      end if
      return
      end

c     
c----------------------------------------------------------------------
c
      function entro(e,n)
      INCLUDE 'defs.f'
c
c  This function-subprogram determines the entropy density of the
c  underlying EoS.
c  It is used in subroutines fileo, prop3d, and the main program.
c
c  Type declarations for variables in common-blocks.
c
      real*8 e0,n0,B,hc3,pi2
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200) 
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer eos,stabil,antie
      real*8 chem
c
c  Type declarations for variables used in function entro.
c
      real*8 e,n,entro,et0,temp
      real*8 de,dn,p1,p2,p3,p4,p13,p24
c      real*8 p12,p34
c
c  Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
      common /eqofst/ eos,stabil,anti

      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3


c
c  Fermi energy density in the QGP for given baryon density
c  (in units of e0).
c
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*n*n0*hc3)**(4d0/3d0)+B
      et0 = et0/e0/hc3
      
      
      antie = 0

      if (n.lt.0)then
         n=-n
         antie = 1
      end if
      
c
c  Calculate entropy density (in units of n0).
c  If flag eos=1, use ideal gas equation of state.   
c  If flag eos=[else], use tabellized equation of state.
c     
c     arithmetic mean interpolation (yields better values for
c     interpolation at small e,n, but is worse for all other
c     thermodynamic quantities, also at larger e,n)
c     
c     
c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     
      if (eos.eq.3) then
         if (e.le.20d0) then
            de = 0.1d0
            dn = 0.05d0
            p1 = stab(idint(e/de),idint(n/dn))
            p2 = stab(idint(e/de)+1,idint(n/dn))
            p3 = stab(idint(e/de),idint(n/dn)+1)
            p4 = stab(idint(e/de)+1,idint(n/dn)+1)
            p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
            p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
            p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
            p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
            entro = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     entro = 0.25d0*(p13+p24+p12+p34)    
            
          else if (e.gt.20d0) then
             entro = 4d0*(e-B/e0/hc3)*e0/3d0/temp(e,n)/n0-
     $            n*chem(e,n)/temp(e,n)
          else
             entro = 0d0
          end if
          
       else if ((eos.eq.4).or.(eos.eq.5).or.(eos.eq.2)) then
          if (e.le.1000.0d0) then
             if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = stab3(idint(e/de),idint(n/dn))
               p2 = stab3(idint(e/de)+1,idint(n/dn))
               p3 = stab3(idint(e/de),idint(n/dn)+1)
               p4 = stab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3

               entro = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if          
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $             ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = stab2(idint(e/de),idint(n/dn))
               p2 = stab2(idint(e/de)+1,idint(n/dn))
               p3 = stab2(idint(e/de),idint(n/dn)+1)
               p4 = stab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3

               entro = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then
               de = 0.5d0
               dn = 0.1d0
               
               p1 = stab(idint(e/de),idint(n/dn))
               p2 = stab(idint(e/de)+1,idint(n/dn))
               p3 = stab(idint(e/de),idint(n/dn)+1)
               p4 = stab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               entro = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if            
            
         else if (e.ge.1000.0d0) then
c     JS muss bestimmt werden
            entro = 100d0
         end if
         
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1).and.(n.ge.-1d0)) 
     $     then
         entro = 0d0
      end if
      
c     if (entro.lt.0d0) entro=0d0
      
      
      
      if (antie.eq.1)then
         n=-n
      end if
      
      return
      end
c
c----------------------------------------------------------------------
c
      function temp(e,n)
      INCLUDE 'defs.f'
c     
c     This function-subprogram determines the temperature of the
c     underlying EoS.
c     It is used in subroutines fileo and entro.
c
c     Type declarations for variables in common-blocks.
c     
      real*8 e0,n0,B,hc3,pi2
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer eos,stabil,antit
c
c     Type declarations for variables used in function temp.
c     
      real*8 e,n,temp,et0,mu
      real*8 de,dn,p1,p2,p3,p4,p13,p24
c     real*8 p12,p34
c     
c     Common-blocks.
c     
      common /grstate/ e0,n0,B,hc3,pi2
      common /eqofst/ eos,stabil,anti
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3


c     
c     Fermi energy density in the QGP for given baryon density
c     (in units of e0).
c     
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*n*n0*hc3)**(4d0/3d0)+B
      et0 = et0/e0/hc3
      
      
      antit = 0
      if (n.lt.0)then
         n=-n
         antit = 1
      end if
c
c     Calculate temperature (in MeV).   
c     If flag eos=[else], use tabellized equation of state.
c     
c     
c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     
      if (eos.eq.3) then
         if (e.le.20d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = ttab3(idint(e/de),idint(n/dn))
               p2 = ttab3(idint(e/de)+1,idint(n/dn))
               p3 = ttab3(idint(e/de),idint(n/dn)+1)
               p4 = ttab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = ttab2(idint(e/de),idint(n/dn))
               p2 = ttab2(idint(e/de)+1,idint(n/dn))
               p3 = ttab2(idint(e/de),idint(n/dn)+1)
               p4 = ttab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then
               de = 0.5d0
               dn = 0.1d0
               p1 = ttab(idint(e/de),idint(n/dn))
               p2 = ttab(idint(e/de)+1,idint(n/dn))
               p3 = ttab(idint(e/de),idint(n/dn)+1)
               p4 = ttab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
         else if (e.gt.20d0) then
            call findqgp(e,n,temp,mu)
         else
            temp = 0d0
         end if
      else if ((eos.eq.4).or.(eos.eq.5).or.(eos.eq.2)) then      
         if (e.le.1000.0d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = ttab3(idint(e/de),idint(n/dn))
               p2 = ttab3(idint(e/de)+1,idint(n/dn))
               p3 = ttab3(idint(e/de),idint(n/dn)+1)
               p4 = ttab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $              ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = ttab2(idint(e/de),idint(n/dn))
               p2 = ttab2(idint(e/de)+1,idint(n/dn))
               p3 = ttab2(idint(e/de),idint(n/dn)+1)
               p4 = ttab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then
               de = 0.5d0
               dn = 0.1d0

               p1 = ttab(idint(e/de),idint(n/dn))
               p2 = ttab(idint(e/de)+1,idint(n/dn))
               p3 = ttab(idint(e/de),idint(n/dn)+1)
               p4 = ttab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if 
            
         else if (e.ge.1000.0d0) then
c     JS muss bestimmt werden
            Temp = 400d0
c     else
            
c     temp = 0d0
            
         end if
         
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1).and.(n.ge.-1d0))
     $     then
         temp = 0d0
      end if
      
      
      if (antit.eq.1)then
         n=-n
      end if    
      return
      end     
      
      
      
      
        

c
c----------------------------------------------------------------------
c
      function chem(e,n)
      INCLUDE 'defs.f'
c
c  This function-subprogram determines the chemical potential of the
c  underlying EoS.
c  It is used in subroutine fileo.
c
c  Type declarations for variables in common-blocks.
c
      real*8 e0,n0,B,hc3,pi2
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer eos,stabil,antic
c
c  Type declarations for variables used in function temp.
c
      real*8 e,n,chem,et0,temp
      real*8 de,dn,p1,p2,p3,p4,p13,p24
c      real*8 p12,p34
c
c  Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
      common /eqofst/ eos,stabil,anti
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3


c
c  Fermi energy density in the QGP for given baryon density
c  (in units of e0).
c
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*n*n0*hc3)**(4d0/3d0)+B
      et0 = et0/e0/hc3
      

      antic = 0
      if (n.lt.0)then
         n=-n
         antic = 1
      end if
c
c  Calculate chemical potential (in MeV).

c     
c  If flag eos=[else], use tabellized equation of state.
c     

c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     
      if (eos.eq.3) then
         if (e.le.20d0) then
            
            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mutab3(idint(e/de),idint(n/dn))
               p2 = mutab3(idint(e/de)+1,idint(n/dn))
               p3 = mutab3(idint(e/de),idint(n/dn)+1)
               p4 = mutab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mutab2(idint(e/de),idint(n/dn))
               p2 = mutab2(idint(e/de)+1,idint(n/dn))
               p3 = mutab2(idint(e/de),idint(n/dn)+1)
               p4 = mutab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then       
               de = 0.5d0
               dn = 0.1d0
               p1 = mutab(idint(e/de),idint(n/dn))
               p2 = mutab(idint(e/de)+1,idint(n/dn))
               p3 = mutab(idint(e/de),idint(n/dn)+1)
               p4 = mutab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13

            end if
            
         else if (e.ge.20d0) then
            call findqgp(e,n,temp,chem)
         else
            chem = 0d0
         end if

      else if ((eos.eq.4).or.(eos.eq.5).or.(eos.eq.2)) then        
         if (e.lt.1000.0d0) then


            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mutab3(idint(e/de),idint(n/dn))
               p2 = mutab3(idint(e/de)+1,idint(n/dn))
               p3 = mutab3(idint(e/de),idint(n/dn)+1)
               p4 = mutab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mutab2(idint(e/de),idint(n/dn))
               p2 = mutab2(idint(e/de)+1,idint(n/dn))
               p3 = mutab2(idint(e/de),idint(n/dn)+1)
               p4 = mutab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then       
               de = 0.5d0
               dn = 0.1d0

               p1 = mutab(idint(e/de),idint(n/dn))
               p2 = mutab(idint(e/de)+1,idint(n/dn))
               p3 = mutab(idint(e/de),idint(n/dn)+1)
               p4 = mutab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
               
            end if
         else if (e.gt.1000.0d0) then
c     JS muss bestimmt werden
            chem = 0.0d0
c     else
c     chem = 0d0
         end if
         
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1)) then
         chem = e0/n0/3
      end if
      
      
      
      
      if (antic.eq.1)then
         n=-n
         chem = -chem
      end if
      return
      end
c     
c----------------------------------------------------------------------
c
      function schem(e,n)
      INCLUDE 'defs.f'
c
c  This function-subprogram determines the chemical potential of the
c  underlying EoS.
c  It is used in subroutine fileo.
c
c  Type declarations for variables in common-blocks.
c
      real*8 e0,n0,B,hc3,pi2
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer eos,stabil,antic
c
c  Type declarations for variables used in function schem.
c
      real*8 e,n,chem,et0,temp,schem
      real*8 de,dn,p1,p2,p3,p4,p13,p24
c      real*8 p12,p34
c
c  Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
      common /eqofst/ eos,stabil,anti
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3


c
c  Fermi energy density in the QGP for given baryon density
c  (in units of e0).
c
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*n*n0*hc3)**(4d0/3d0)+B
      et0 = et0/e0/hc3
      

      antic = 0
      if (n.lt.0)then
         n=-n
         antic = 1
      end if
c
c  Calculate chemical potential (in MeV).     
c  If flag eos=[else], use tabellized equation of state.
c     

c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     
      if (eos.eq.3) then
         if (e.le.20d0) then
            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mustab3(idint(e/de),idint(n/dn))
               p2 = mustab3(idint(e/de)+1,idint(n/dn))
               p3 = mustab3(idint(e/de),idint(n/dn)+1)
               p4 = mustab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mustab2(idint(e/de),idint(n/dn))
               p2 = mustab2(idint(e/de)+1,idint(n/dn))
               p3 = mustab2(idint(e/de),idint(n/dn)+1)
               p4 = mustab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then        
               
               de = 0.5d0
               dn = 0.1d0
               p1 = mustab(idint(e/de),idint(n/dn))
               p2 = mustab(idint(e/de)+1,idint(n/dn))
               p3 = mustab(idint(e/de),idint(n/dn)+1)
               p4 = mustab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13         
              
            end if
         else if (e.ge.20d0) then
            schem = 0d0
         end if
         
      else if ((eos.eq.4).or.(eos.eq.5).or.(eos.eq.2)) then        
         if (e.lt.1000.0d0) then
            
            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mustab3(idint(e/de),idint(n/dn))
               p2 = mustab3(idint(e/de)+1,idint(n/dn))
               p3 = mustab3(idint(e/de),idint(n/dn)+1)
               p4 = mustab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
                  
               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mustab2(idint(e/de),idint(n/dn))
               p2 = mustab2(idint(e/de)+1,idint(n/dn))
               p3 = mustab2(idint(e/de),idint(n/dn)+1)
               p4 = mustab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               

               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then      
               de = 0.5d0
               dn = 0.1d0

               p1 = mustab(idint(e/de),idint(n/dn))
               p2 = mustab(idint(e/de)+1,idint(n/dn))
               p3 = mustab(idint(e/de),idint(n/dn)+1)
               p4 = mustab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     chem = 0.25d0*(p13+p24+p12+p34)    
            end if
         else if (e.gt.1000.0d0) then
c     JS muss bestimmt werden
            schem = 0d0
c     else
c     chem = 0d0
         end if
         
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1)) then
         schem = 0
      end if
      
      
      
      
      if (antic.eq.1)then
         n=-n
         schem = schem
      end if
      return
      end
c 
c----------------------------------------------------------------------
c
      function cs(e,n)
      INCLUDE 'defs.f'
c
c  This function-subprogram determines the chemical potential of the
c  underlying EoS.
c  It is used in subroutine fileo.
c
c  Type declarations for variables in common-blocks.
c
      real*8 e0,n0,B,hc3,pi2
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer eos,stabil,antic
c
c  Type declarations for variables used in function schem.
c
      real*8 e,n,chem,et0,temp,schem,cs
      real*8 de,dn,p1,p2,p3,p4,p13,p24
c      real*8 p12,p34
c
c  Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
      common /eqofst/ eos,stabil,anti
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3
c
c  Fermi energy density in the QGP for given baryon density
c  (in units of e0).
c
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*n*n0*hc3)**(4d0/3d0)+B
      et0 = et0/e0/hc3
      

      antic = 0
      if (n.lt.0)then
         n=-n
         antic = 1
      end if
c
c  Calculate chemical potential (in MeV).
c     
c  If flag eos=[else], use tabellized equation of state.
c     

c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     
      if (eos.eq.3) then
         if (e.le.20d0) then
            de = 0.2d0
            dn = 0.1d0
            p1 = cstab(idint(e/de),idint(n/dn))
            p2 = cstab(idint(e/de)+1,idint(n/dn))
            p3 = cstab(idint(e/de),idint(n/dn)+1)
            p4 = cstab(idint(e/de)+1,idint(n/dn)+1)
            p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
            p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
            p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
            p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3

            cs = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c            chem = 0.25d0*(p13+p24+p12+p34)         
            
         else if (e.ge.20d0) then
            call findqgp(e,n,temp,chem)
         else
            cs = 0d0
         end if

      else if ((eos.eq.4).or.(eos.eq.5).or.(eos.eq.2)) then        
         if (e.lt.1000.0d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = cstab3(idint(e/de),idint(n/dn))
               p2 = cstab3(idint(e/de)+1,idint(n/dn))
               p3 = cstab3(idint(e/de),idint(n/dn)+1)
               p4 = cstab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3


               cs = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = cstab2(idint(e/de),idint(n/dn))
               p2 = cstab2(idint(e/de)+1,idint(n/dn))
               p3 = cstab2(idint(e/de),idint(n/dn)+1)
               p4 = cstab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            


               cs = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then 
               de = 0.5d0
               dn = 0.1d0
               
               p1 = cstab(idint(e/de),idint(n/dn))
               p2 = cstab(idint(e/de)+1,idint(n/dn))
               p3 = cstab(idint(e/de),idint(n/dn)+1)
               p4 = cstab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               cs = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     chem = 0.25d0*(p13+p24+p12+p34)    
            end if
         else if (e.gt.1000.0d0) then
c     JS muss bestimmt werden
            cs = 1/sqrt(3d0)
c     else
c     chem = 0d0
         end if
         
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1)) then
         cs = 0
      end if
     
 
      
      
      
      if (antic.eq.1)then
         n=-n
         cs = cs
      end if
      return
      end
c     
c----------------------------------------------------------------------    
c----------------------------------------------------------------------
c
      function lambda(e,n)
      INCLUDE 'defs.f'
c
c  This function-subprogram determines the fraction of QGP of the
c  underlying EoS.
c  It is used in subroutine prop3d and the main program.
c
c  Type declarations for variables in common-blocks.
c
      real*8 e0,n0,B,hc3,pi2
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400)  
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer eos,stabil,antil
c
c  Type declarations for variables used in function lambda.
c
      real*8 e,n,et0,lambda
      real*8 de,dn,p1,p2,p3,p4,p13,p24
c      real*8 p12,p34
c
c  Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
      common /eqofst/ eos,stabil,anti
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3
c
c  Fermi energy density in the QGP for given baryon density
c  (in units of e0).
c
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*n*n0*hc3)**(4d0/3d0)+B
      et0 = et0/e0/hc3
      

      antil = 0

      if (n.lt.0)then
         n=-n
         antil = 1
      end if
      
c     
c     Calculate fraction of QGP.
c     If eos = 1 or eos = 2, lambda is zero.
c     
      if ((eos.eq.2).or.(eos.eq.4)) then
         lambda = 0d0
c     
c     If eos=[else], use tabellized equation of state.
c     
      else
         if (e.le.20d0) then
            de = 0.1d0
            dn = 0.05d0
            p1 = lamtab(idint(e/de),idint(n/dn))
            p2 = lamtab(idint(e/de)+1,idint(n/dn))
            p3 = lamtab(idint(e/de),idint(n/dn)+1)
            p4 = lamtab(idint(e/de)+1,idint(n/dn)+1)
            p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
            p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
c     

c     
c     Linear interpolation.
c     
            lambda = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     
c     arithmetic mean interpolation (yields better values for
c     interpolation at small e,n, but is worse for all other
c     thermodynamic quantities, also at larger e,n)
c     
c     lambda = 0.25d0*(p13+p24+p12+p34) 
c     
         else if (e.gt.20d0) then
            lambda = 1d0
         else
            lambda = 0d0
         end if
      end if
      if (lambda.gt.1d0) lambda=1d0
      
      if (antil.eq.1)then
         n=-n
      end if
      return
      end
c     
c------------------------------------------------------------------------
c
      subroutine findqgp(energ,dens,t,mu)
      INCLUDE 'defs.f'
c
c  This subroutine-subprogram finds T, mu in the QGP for given e,n.
c  It is used in subroutines temp and chem.
c
      real*8 e0,n0,B,hc3,pi2
      real*8 energ,dens,t,mu
      real*8 pi
      real*8 a1,a2,a3,a4
      real*8 et0
      integer l
c     
c     Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
c     
c  Define constants.
c     
      pi = dsqrt(pi2)
      l = 0
c     
c     Convert energy density and baryon density, 
c     so far given in units of e0 and n0, into MeV**4 
c     resp. MeV**3, so that T and mu emerge in units of MeV.
c     
      energ = energ*e0*hc3
      dens = dens*n0*hc3
c     
c     If e < e(T=0) (the Fermi energy density in the QGP for given
c     baryon density), set everything to zero.
c     
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*dens)**(4d0/3d0)+B
      if (energ.lt.et0) then
         t = 0d0
        mu = 0d0
      else
c     
c     If density is larger than zero, complicated sixth
c     order polynomial in mu has to be solved...
c     
         if (dens.gt.0d0) then
c     
c     Calculate coefficients of sixth order polynomial.
c     
          a1 = 4d0/1215d0/pi2
          a2 = - 4d0/15d0*dens
          a3 = energ - B
          a4 = - 999d0*pi2*dens*dens/40d0
          call froot(mu,a1,a2,a3,a4,1d-7,l)
c
c     After having found mu, calculate T.
c     
          t = dsqrt(4.5d0*dens/mu - mu*mu/9d0/pi2)
c     
c     ...if density is zero, mu = 0 and t is fourth root
c     of e-B (up to a constant). 
c     
        else
           mu = 0d0
           if(energ.ge.B) then
              t = ((energ-B)*30d0/37d0/pi2)**0.25d0
           else
            t=0d0
         end if
      end if
      end if 
c     
c     Finally, convert e, n back into dimensionless
c     units, so that the original quantities will
c     be returned to the main program.
c     
      energ = energ/e0/hc3
      dens = dens/n0/hc3
      return
      end 
c     
c------------------------------------------------------------------------
c
      subroutine froot(mu,a1,a2,a3,a4,acc,l)
      INCLUDE 'defs.f'
c     
c  This subroutine-subprogram finds the root of the sixth order 
c  polynomial in mu.
c  It is used in subprogram findqgp.
c
      real*8 e0,n0,B,hc3,pi2
      real*8 mu,a1,a2,a3,a4,acc
      real*8 mua,mub,fa,fb,fc
      integer i,l
c
c  Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
c
      i = 0
c
c  Lower boundary for interval search in mu is mua=0.
c
      mua = 0d0
      fa = a4
c
c  Upper boundary for interval search in mu corresponds
c  to the maximum possible mu, ie, that mu which
c  corresponds to T=0 and baryon density 12*n0.
c
      mub = (40.5d0*pi2*12d0*n0*hc3)**(1d0/3d0)
      fb = a4 + a3*mub*mub + a2*mub**3 + a1*mub**6
 22   mu = (mua + mub)*0.5d0
      fc = a4 + a3*mu*mu + a2*mu**3 + a1*mu**6
      i = i+1
      if (i.gt.1000) then
        write(6,*) ' froot: more than 1000 iterations to find mu'
        write(6,*) '  fc/e0= ',fc/e0/hc3
        goto 23
      end if
      if (dabs(fa)/e0/hc3.lt.acc) then
        mu = mua
        goto 23
      end if
      if (dabs(fb)/e0/hc3.lt.acc) then
        mu = mub
        goto 23
      end if  
      if (fb*fc.lt.0d0) then
        mua = mu
        fa = fc
      else if (fa*fc.lt.0d0) then
        mub = mu
        fb = fc
      else
c        write(6,*) ' froot: 2 or no zeros in search interval'   
        l = 1
        return
      end if
      if (dabs(fc)/e0/hc3.gt.acc) goto 22
 23   continue
      l = 0
      return
      end 
c
c----------------------------------------------------------------------
c
      function signum(a)
c
c  This function-subprogram calculates the sign.
c  It is used in subroutines prop and fileo.
c
c  Type-declarations for variables used in function-subprogram signum.
c
      real*8 a,signum
      if (dabs(a).gt.0d0) then
        signum = a/dabs(a)
      else
        signum = 0d0
      end if
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine prop(elab,mx,my,mz,r,vx,vy,vz,p)
      INCLUDE 'defs.f'
c
c  This subroutine-subprogram propagates the laboratory frame
c  quantities Elab, M_x, M_y, M_z, R. It is used in subroutine prop3d.
c
c  Propagation direction is ALWAYS the x-direction!
c
c  Type declarations for variables in common-blocks.
c
      real*8 dx,dt,t,vol,cut
      real*8 diff
      integer kp(ngr),km(ngr)
c
c  Type declarations for variables in subroutine prop.
c
      real*8 elab(ngr),mx(ngr),my(ngr),mz(ngr),r(ngr)
      real*8 vx(ngr),vy(ngr),vz(ngr),p(ngr)
      real*8 eps(ngr),qp(ngr),qm(ngr)
      real*8 et(ngr),mxt(ngr),myt(ngr),mzt(ngr),rt(ngr)
      real*8 dele(ngr),delmx(ngr),delmy(ngr),delmz(ngr),delr(ngr)
      real*8 delte(ngr),deltmx(ngr),deltmy(ngr),deltmz(ngr),deltr(ngr)
      real*8 ae(ngr),amx(ngr),amy(ngr),amz(ngr),ar(ngr)
      real*8 ate(ngr),atmx(ngr),atmy(ngr),atmz(ngr),atr(ngr)
      real*8 lam,sgn,mini,mxo,myo,dfe,dfm
      real*8 signum
      integer i
c
c  Common-blocks.
c
      common /gitter/ dx,dt,t,vol,cut
      common /diffus/ diff
      common /indpoi/ kp,km
c
c  Start calculation.
c
c  Step 1: Calculate epsilon and the cell differences.
c
      lam = dt/dx
      do 1999 i = 1,ngr
         eps(i) = vx(i)*lam
         dele(i) = elab(kp(i)) - elab(i)
         delmx(i) = mx(kp(i)) - mx(i)
         delmy(i) = my(kp(i)) - my(i)
         delmz(i) = mz(kp(i)) - mz(i)
         delr(i) = r(kp(i)) - r(i)
 1999 continue  
c
c  Step 2: Calculate transported and diffused quantities,
c
      do 2000 i = 1,ngr
         qp(i) = (0.5d0 - eps(i))/(1d0 + eps(kp(i))-eps(i))
         qm(i) = (0.5d0 + eps(i))/(1d0 - eps(km(i))+eps(i))
         dfe = -0.5d0*(p(kp(i))*vx(kp(i))-p(km(i))*vx(km(i)))
         dfm = -0.5d0*(p(kp(i))-p(km(i)))
c
c  There are NO source terms for M_y, M_z, and R!
c
         et(i) = 0.5d0*qp(i)**2*dele(i)-0.5d0*qm(i)**2*dele(km(i))
     ?        + (qp(i) + qm(i))*elab(i) + dfe*lam
         mxt(i) = 0.5d0*qp(i)**2*delmx(i)-0.5d0*qm(i)**2*delmx(km(i))
     ?        + (qp(i) + qm(i))*mx(i) + dfm*lam
         myt(i) = 0.5d0*qp(i)**2*delmy(i)-0.5d0*qm(i)**2*delmy(km(i))
     ?       + (qp(i) + qm(i))*my(i) 
         mzt(i) = 0.5d0*qp(i)**2*delmz(i)-0.5d0*qm(i)**2*delmz(km(i))
     ?        + (qp(i) + qm(i))*mz(i) 
         rt(i) = 0.5d0*qp(i)**2*delr(i)-0.5d0*qm(i)**2*delr(km(i))
     ?        + (qp(i) + qm(i))*r(i) 
 2000 continue
c
c  Step 3: Calculate differences between the transported and diffused
c          quantities and naive phoenical antidiffusion fluxes.
c
      do 2001 i = 1,ngr
         delte(i) = et(kp(i)) - et(i)
         deltmx(i) = mxt(kp(i)) - mxt(i)
         deltmy(i) = myt(kp(i)) - myt(i)
         deltmz(i) = mzt(kp(i)) - mzt(i)
         deltr(i) = rt(kp(i)) - rt(i)
         ae(i) = diff * 0.125d0*( delte(i) - 
     1       0.125d0*(dele(kp(i))-2d0*dele(i)+dele(km(i))) )
         amx(i) = diff * 0.125d0*( deltmx(i) - 
     1        0.125d0*(delmx(kp(i))-2d0*delmx(i)+delmx(km(i))) )
        amy(i) = diff * 0.125d0*( deltmy(i) - 
     1        0.125d0*(delmy(kp(i))-2d0*delmy(i)+delmy(km(i))) )
        amz(i) = diff * 0.125d0*( deltmz(i) - 
     1       0.125d0*(delmz(kp(i))-2d0*delmz(i)+delmz(km(i))) )
        ar(i) = diff * 0.125d0*( deltr(i) - 
     1       0.125d0*(delr(kp(i))-2d0*delr(i)+delr(km(i))) )
 2001 continue
c     
c     Step 4: Calculate Flux Corrected Antidiffusion Fluxes.
c
      do 2004 i = 1,ngr
         sgn = signum(ae(i))
         mini = dmin1(sgn*delte(km(i)),dabs(ae(i)),sgn*delte(kp(i)))
        ate(i) = sgn * dmax1(0d0,mini)
        sgn = signum(amx(i))
        mini = dmin1(sgn*deltmx(km(i)),dabs(amx(i)),sgn*deltmx(kp(i)))
        atmx(i) = sgn * dmax1(0d0,mini)
        sgn = signum(amy(i))
        mini = dmin1(sgn*deltmy(km(i)),dabs(amy(i)),sgn*deltmy(kp(i)))
        atmy(i) = sgn * dmax1(0d0,mini)
        sgn = signum(amz(i))
        mini = dmin1(sgn*deltmz(km(i)),dabs(amz(i)),sgn*deltmz(kp(i)))
        atmz(i) = sgn * dmax1(0d0,mini)
        sgn = signum(ar(i))
        mini = dmin1(sgn*deltr(km(i)),dabs(ar(i)),sgn*deltr(kp(i)))
        atr(i) = sgn * dmax1(0d0,mini)
 2004 continue
c
c  Step 5: Calculate final densities.
c
      do 2005 i = 1,ngr
         elab(i) = et(i) - ate(i) + ate(km(i)) 
         mx(i) = mxt(i) - atmx(i) + atmx(km(i))
         my(i) = myt(i) - atmy(i) + atmy(km(i))
         mz(i) = mzt(i) - atmz(i) + atmz(km(i))
         r(i) = rt(i) - atr(i) + atr(km(i))
2005  continue
c
c     Step 6: Finally, a consistency check removes acausalities
c     produced in the propagation.
c
      do 2016 i = 1,ngr
         if (elab(i).lt.dsqrt(mx(i)**2+my(i)**2+mz(i)**2)) then
            if (dsqrt(mx(i)**2+my(i)**2+mz(i)**2).gt.1d-16) then
               mxo = mx(i)
               myo = my(i)
               mx(i) = elab(i)*(1d0-1d-12)*mxo
     ?              /dsqrt(mxo**2+myo**2+mz(i)**2)
               my(i) = elab(i)*(1d0-1d-12)*myo
     ?              /dsqrt(mxo**2+myo**2+mz(i)**2)
               mz(i) = elab(i)*(1d0-1d-12)*mz(i)
     ?              /dsqrt(mxo**2+myo**2+mz(i)**2)
c              write(6,*) '  M(',i,') changed to Elab at Step 6'
            else
               mx(i) = 0d0
               my(i) = 0d0
               mz(i) = 0d0
            end if
         end if
 2016 continue
c
c     Step 7: Cells with smaller lab energy density than a cut-off
c     are regarded as vacuum.
c
      do 2017 i = 1,ngr
         if (elab(i).le.cut) then
            elab(i) = 0d0
            mx(i) = 0d0
            my(i) = 0d0
            mz(i) = 0d0
            r(i) = 0d0
         end if
 2017 continue
      return
      end
c
c----------------------------------------------------------------------
c

      SUBROUTINE onefluid(testflag,ithydro,thydro_start,ifr,jfr,kfr,tf,
     &  dstfr,dsxfr,dsyfr,dszfr,tfr,muqfr,musfr,vxfr,vyfr,vzfr,mmax)
      implicit none

c     main program
c     mark
c     Type declarations for variables in common-blocks.
c     
      include 'options.f'
      INCLUDE 'defs.f'
      
      real*8 e0,n0,Bag,hc3,pi2
      real*8 dx,dt,t,vol,cut
      real*8 total0(6)
      integer eos,stabil,avout 
c     
c     Type declarations for variables and function subprograms used in main.
c
      real*4 dxfr(ngr*ngr*ngr)
      real*8 e(ngr,ngr,ngr),p(ngr,ngr,ngr),n(ngr,ngr,ngr)
      real*8 vx(ngr,ngr,ngr),vy(ngr,ngr,ngr),vz(ngr,ngr,ngr)
      real*8 gamma(ngr,ngr,ngr),lam(ngr,ngr,ngr)
      real*8 elab(ngr,ngr,ngr),r(ngr,ngr,ngr)
      real*8 mx(ngr,ngr,ngr),my(ngr,ngr,ngr),mz(ngr,ngr,ngr)
      real*8 s,entro,lambda,temp,chem,press,pr,schem,cs
      real*8 a1,a2,b,v,etot1
      integer itmax,iniflg,stophydro
      integer it,i,j,k
      integer kp(ngr),km(ngr),testflag,ithydro
      real*8 etothy,pxtothy,pytothy,pztothy,bartothy,entrotota
      real*8 qgpfrac,qgp_overall,enorm,h_etot,h_btot,entrotote
      real*8 diff,m,velo,etot

c new variables for OSCAR output    
      integer step,number,vistime,itend,hypertime,hyperstep      

c     New variables for isoerg  freezeout


      real*4 ifr(ngr*ngr*ngr),jfr(ngr*ngr*ngr),kfr(ngr*ngr*ngr)
      integer frozen(ngr),frozenb(ngr),mmax,mfr
      real*4 tf(ngr*ngr*ngr) ,tfr(ngr*ngr*ngr), muqfr(ngr*ngr*ngr)
      real*4 musfr(ngr*ngr*ngr)
      real*8 vxfr(ngr*ngr*ngr), vyfr(ngr*ngr*ngr)
      real*8 vzfr(ngr*ngr*ngr),gfr(ngr*ngr*ngr)
      real*8 sfr,efr,rfr
      real*8 dstfr(ngr*ngr*ngr),dsxfr(ngr*ngr*ngr)
      real*8 dsyfr(ngr*ngr*ngr),dszfr(ngr*ngr*ngr)       


      real*8 thydro_start,entrotot
      integer anti

chp new arrays for Cornelius hypersurface finding
      real*8 elast(ngr,ngr,ngr),nlast(ngr,ngr,ngr)
      real*8 vxlast(ngr,ngr,ngr),vylast(ngr,ngr,ngr),vzlast(ngr,ngr,ngr)
      REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1,0:1) :: Hypervx,Hypervy
      REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1,0:1) :: Hypervz,Hypern
      REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1,0:1) :: HyperCube
      INTEGER         :: Ngp
      REAL(KIND(0D0)),DIMENSION(0:3,8) :: dSigma
      INTEGER         :: Nsurf
      REAL(KIND(0D0)),DIMENSION(0:3,8) :: Vmid
      INTEGER         :: Nambi, Ndisc

      REAL(KIND(0D0)) :: Quadrilinear, Vxmid, Vymid, Vzmid ! Interpolation
      real*8 Tmid, muqmid,musmid,nmid,emid,pressmid,gammid,dxh
      real*8 t0,x0,y0,z0, Hypert,Hyperx,Hypery,Hyperz,dthyper
      integer ii, over, under
c total energy and net baryon number on iso-e hypersurface
      real*8 bflow,umudsigmu,ene,tnn,tnx,tny,tnz
      real*8 negbflow,posbflow,negene,posene 
c
c
c     Common-blocks.
c     
      common /grstate/ e0,n0,Bag,hc3,pi2
      common /gitter/ dx,dt,t,vol,cut
      common /eqofst/ eos,stabil,anti
      common /indpoi/ kp,km
      common /diffus/ diff
      common /propag/ step
                  
      
chp the grid, we get it outside
      common /onefgrid/ elab,mx,my,mz,r
 
chp common block for energy conservation check      
      common /h_etot/ h_etot,h_btot
      
 

c
c     Define constants.
c     
      pi2 = dacos(-1d0)**2
      hc3 = 197.327053d0**3d0 
      e0 = 146.51751415742d0 
      n0 = 0.15891d0 
      Bag = 235d0**4 

c initialization
c     
c     Calculate index pointer (used in subroutine prop).
c
      do 115 i = 1,ngr
        kp(i) = i+1
        if (i.eq.ngr) kp(i) = ngr
        km(i) = i-1
        if (i.eq.1) km(i) = 1
 115  continue

c.. init flags etc.
      itmax=1000 ! number of time steps
c.. is set in uhmerge.f      dx=0.2d0  ! cell size (fm)
      dt=0.4d0*dx ! time step (fm)
      diff=1d0  ! diffusion coefficient in fractions of 1/8
      step=2    ! half step prop. (1=full tep)
      cut=1d-12 ! Energy density (in e/e0) regarded as vacuum
      stabil=0  ! no stabilization of nucleus by hand
      EoS=CTOption(47)
      if(EoS.eq.1)EoS=5
c EoS: 1=5, 2= HG, 3=Bag model w/PT, 4=chiral w/PT, 5=chiral+HG
c
c     Read in equations of state. eos=2: hadronic EoS, 
c     eos=3: EoS with p.t. to QGP (bag modell mu = mu_b) eos=4:
c     chiral EoS with critical endpoint (mu = mu_q)
c     eos=5: chiral EoS but with hadronic T and mu_q.
c    
     
chp file number for OSCAR output
       number=21      
chp vistime*dt gives timesteps to write output
       vistime=8
chp hypertime is a parameter that defines at which timesteps hypersurface is checked
       hypertime=CTParam(71)
       hyperstep=hypertime
       itend=0

      mfr=1
      efr=0d0
      sfr=0d0
      rfr=0d0
CCCC test case do N timesteps ini hydro with HG EoS
      if (testflag.gt.0) then
       itmax=testflag
      end if
CCCCC

      if (eos.eq.2) then
         call readeos1()
      else if (eos.eq.3) then
         call readeos2()
      else if (eos.eq.4) then
         call readeos3()
      else if (eos.eq.5) then
         call readeos3()
      end if
      
c     
c     Zeroth time step of time evolution.
c     
      t = 0.d0 


c     
c     Calculate volume of a cell.
c     

      vol = dx*dx*dx
      entrotota=0d0
      entrotote=0d0
      etot1=0d0
      ene=0d0
      bflow=0d0
      posbflow=0.0d0
      negbflow=0.0d0
      posene=0.0d0
      negene=0.0d0 
c
c
c     Do a special untangle to obtain rest frame quantities.
c

      do 443 k = 1,ngr
         do 442 j = 1,ngr
            do 441 i = 1,ngr
               
               m = dsqrt(mx(i,j,k)*mx(i,j,k) + my(i,j,k)*my(i,j,k)
     $              + mz(i,j,k)*mz(i,j,k))
               frozen(i)=0
               frozenb(i)=0
               if (elab(i,j,k).le.1d-12) then
                  elab(i,j,k) = 0d0
                  mx(i,j,k) = 0d0
                  my(i,j,k) = 0d0
                  mz(i,j,k) = 0d0
                  r(i,j,k) = 0d0
                  e(i,j,k) = 0d0
                  vx(i,j,k) = 0d0
                  vy(i,j,k) = 0d0
                  vz(i,j,k) = 0d0
                  n(i,j,k) = 0d0
                  gamma(i,j,k) = 1d0
               else if (m.lt.1d-16) then
                  e(i,j,k) = elab(i,j,k)
                  vx(i,j,k) = 0d0
                  vy(i,j,k) = 0d0
                  vz(i,j,k) = 0d0
                  n(i,j,k) = r(i,j,k)
                  gamma(i,j,k) = 1d0
               else
                  v = velo(elab(i,j,k),m,r(i,j,k))
                  if (v.lt.1d0-1d-16) then
                     e(i,j,k) = elab(i,j,k) - v*m
                     vx(i,j,k) = v*mx(i,j,k)/m
                     vy(i,j,k) = v*my(i,j,k)/m
                     vz(i,j,k) = v*mz(i,j,k)/m
                     gamma(i,j,k) = 1d0/dsqrt(1d0-v*v)
                     n(i,j,k) = r(i,j,k)/gamma(i,j,k)
c     Vacuum is assumed if v .ge. 1 - 10^{-16}.
                  else
                     elab(i,j,k) = 0d0
                     mx(i,j,k) = 0d0
                     my(i,j,k) = 0d0
                     mz(i,j,k) = 0d0
                     r(i,j,k) = 0d0
                     e(i,j,k) = 0d0
                     vx(i,j,k) = 0d0
                     vy(i,j,k) = 0d0
                     vz(i,j,k) = 0d0
                     n(i,j,k) = 0d0
                     gamma(i,j,k) = 1d0
                  end if
               end if
        
               p(i,j,k) = press(e(i,j,k),n(i,j,k))
               if(elab(i,j,k).gt.1d-12) then
                  etot1=etot1+elab(i,j,k)*vol*e0
               entrotota=entrotota+entro(e(i,j,k),n(i,j,k))*vol
     $                 *n0*gamma(i,j,k)
               endif
 441        continue
 442     continue
 443  continue
      
      
chp write header for OSCAR output
      if(CTOption(54).gt.0)then  
       call oschydro_header(thydro_start,itmax,number) 
      end if

chp collect cells below switching criterion before evolution
      do i=1,ngr
       do j=1,ngr
        do k=1,ngr
         if(elab(i,j,k).gt.1d-8)then
          if(CTOption(52).eq.2)then
           if(e(i,j,k).le.CTParam(64))then
            tf(mfr)=t
            ifr(mfr)=(float(i)-float(ngr)*0.5d0-0.5d0)*sngl(dx)
            jfr(mfr)=(float(j)-float(ngr)*0.5d0-0.5d0)*sngl(dx)
            kfr(mfr)=(float(k)-float(ngr)*0.5d0-0.5d0)*sngl(dx)    
            dstfr(mfr)=dx*dx*dx
            dsxfr(mfr)=0.0d0
            dsyfr(mfr)=0.0d0
            dszfr(mfr)=0.0d0  
            tfr(mfr)=Temp(e(i,j,k),n(i,j,k))
            muqfr(mfr)=chem(e(i,j,k),n(i,j,k))
            musfr(mfr)=schem(e(i,j,k),n(i,j,k))
            vxfr(mfr)= vx(i,j,k)
            vyfr(mfr)= vy(i,j,k)
            vzfr(mfr)= vz(i,j,k)
            gfr(mfr)= gamma(i,j,k)
            dxfr(mfr)=dx
            
            mmax=mfr
            efr=efr+elab(i,j,k)*vol*e0/1000d0
            rfr=rfr+r(i,j,k)*vol*n0
            sfr=sfr+entro(e(i,j,k),n(i,j,k))*vol*n0
     $        *gamma(i,j,k)
            bflow=bflow+n(i,j,k)*vol*n0
            ene=ene+e(i,j,k)*vol*e0/1.d3    
           if(CTOption(54).eq.2)then 
            call oschyper(number,tf(mfr),ifr(mfr),jfr(mfr),kfr(mfr),
     &         e(i,j,k),press(e(i,j,k),n(i,j,k)),n(i,j,k),vxfr(mfr),
     &         vyfr(mfr),vzfr(mfr),tfr(mfr),muqfr(mfr),musfr(mfr),
     &         dstfr(mfr),dsxfr(mfr),dsyfr(mfr),dszfr(mfr))
            end if
c
            mfr=mfr+1
           end if
          end if
         end if 
        end do
       end do
      end do 

      posbflow=bflow
      posene=ene

      bflow=0.d0
      ene=0.d0
c     Start time evolution.
c     
      stophydro=0
      
      do 99 it = 1,itmax
         t = t + dt
         ithydro=it
         etot=0d0
         entrotot=0d0
c         write(*,*) 'thydro=', t

         do 463 k = 1,ngr
            do 462 j = 1,ngr
               do 461 i = 1,ngr
                  if(elab(i,j,k).gt.1d-12)then
                     etot=etot+ elab(i,j,k)*vol*e0                                
                  endif
 461           continue
 462        continue
 463     continue
         




c     first step energy is conserved by definition
         
         if (it.eq.1) then
            etot=etot1
         endif
       
c         write(*,*) t,etot, ((etot1-etot)/etot1)

c     freezing out planes
       if(CTOption(52).lt.2) then   
         do 777 k=1,ngr
            if(frozen(k).eq.0) then
               stophydro=0   
               do 778 j=1,ngr
                  do 779 i=1,ngr
                     if (e(i,j,k).gt.CTParam(64)) then         
                        stophydro = 1
                     endif
 779              continue
 778           continue
               if ((stophydro.eq.0).and.(frozenb(k).eq.1)) then
c                  write(*,*) "Freezing out plane Nr.: ", k 
                  do 458 j=1,ngr
                     do 459 i=1,ngr
                        if (elab(i,j,k).gt.1d-8)then
                           tf(mfr)=t-dt
                           ifr(mfr)=(float(i)-float(ngr)*.5-.5)*sngl(dx)
                           jfr(mfr)=(float(j)-float(ngr)*.5-.5)*sngl(dx)
                           kfr(mfr)=(float(k)-float(ngr)*.5-.5)*sngl(dx)
                           dstfr(mfr)=dx*dx*dx
                           dsxfr(mfr)=0.0d0
                           dsyfr(mfr)=0.0d0
                           dszfr(mfr)=0.0d0  
                           tfr(mfr)=Temp(e(i,j,k),n(i,j,k))
                           muqfr(mfr)=chem(e(i,j,k),n(i,j,k))
                           musfr(mfr)=schem(e(i,j,k),n(i,j,k))
                           vxfr(mfr)= vx(i,j,k)
                           vyfr(mfr)= vy(i,j,k)
                           vzfr(mfr)= vz(i,j,k)
                           gfr(mfr)= gamma(i,j,k)
                           dxfr(mfr)=dx
                           if(CTOption(54).eq.2 .and. CTOption(52).eq.0)
     &                      then 
                             call oschyper(number,tf(mfr),ifr(mfr),
     &                               jfr(mfr),kfr(mfr),e(i,j,k),
     &                               press(e(i,j,k),n(i,j,k)),n(i,j,k),
     &                               vxfr(mfr),vyfr(mfr),vzfr(mfr),
     &                               tfr(mfr),muqfr(mfr),musfr(mfr),
     &                               dstfr(mfr),dsxfr(mfr),dsyfr(mfr),
     &                               dszfr(mfr))
                           end if
                           mfr=mfr+1
                           mmax=mfr
                           frozen(k)=1
                           efr=efr+elab(i,j,k)*vol*e0/1000d0
                           rfr=rfr+r(i,j,k)*vol*n0
                           sfr=sfr+entro(e(i,j,k),n(i,j,k))*vol*n0
     $                          *gamma(i,j,k)
                        endif
 459                 continue
 458              continue
               end if
               if (stophydro.eq.1)then
                  frozenb(k)=1
               end if
               stophydro=1
            end if
 777     continue 
         endif  
 


        
c     if energy loss is greater than 0.1 % then increase dx
         
         if(((etot1-etot)/etot1).gt.0.001d0)then
            write(*,*)'Energy loss:',((etot1-etot)/etot1)*100d0,
     $           '% Changing grid !'            
            
            call changegrid(e,p,n,vx,vy,vz,gamma,efr,rfr,sfr,frozen,
     $     frozenb,mfr,etot1,ifr,jfr,kfr,tf,dstfr,dsxfr,dsyfr,dszfr,
     $     tfr,muqfr,musfr,vxfr,vyfr,vzfr,mmax,elast,vxlast,vylast,
     $     vzlast,nlast)
            number=number+3
            vistime=vistime/2
            itend=it
           if(CTOption(54).eq.1)then
            call oschydro_header(thydro_start,itmax,number)
           endif 
         endif
            
        
chp store grid information of old timestep
	   if((CTOption(52).eq.2).and.
     &     ((it.eq.1).or.(mod(it-1,hypertime).eq.0)))then
          do k=1,ngr
           do j=1,ngr
            do i=1,ngr
              elast(i,j,k) = e(i,j,k)
              vxlast(i,j,k) = vx(i,j,k)
              vylast(i,j,k) = vy(i,j,k)
              vzlast(i,j,k) = vz(i,j,k)
              nlast(i,j,k) = n(i,j,k)
            end do
           end do
          end do                       
         endif      

         

        
         
c     JS Criterium for aborting the evolution when eos=4 or 2:
         if (((eos.eq.4).or.(eos.eq.2).or.(eos.eq.5)))then
            stophydro=0
            do 666 k = 1,ngr
               do 667 j = 1,ngr
                  do 668 i = 1,ngr            
c     hp ctp 64 defines the freezeout criterium, default is 5.0     
                     if (e(i,j,k).gt.CTParam(64)) then
                        stophydro = 1
                     endif
 668              continue
 667           continue
 666        continue
         endif
         
c     hp Criterium for aborting the evolution when eos=3:
         if (eos.eq.3)then
            stophydro=0
            do k = 1,ngr
               do j = 1,ngr
                  do i = 1,ngr                 
                     qgpfrac=qgpfrac+lam(i,j,k)*e(i,j,k)
                     enorm=enorm+e(i,j,k)                
c     hp try what happens with the same value as for the HG
                     if (e(i,j,k).gt.CTParam(64)) then
                        stophydro = 1
                     endif
                  end do
               end do 
            end do
            qgp_overall=qgpfrac/enorm
c     hp freeze-out criterium for bag model (lt 20% QGP) 
c     if (qgp_overall.gt.0.20d0) then
c     stophydro = 1
c     endif
            if(stophydro.eq.0)then
               write(*,*) 'qgp_overall after bag hydro:',qgp_overall
            end if  
         endif
         
         
         
         if (stophydro.eq.0) goto 300           
c     
c     Call 3d-propagation routine.
c     

       if(CTOption(54).eq.1)then
        if(mod(it,vistime).eq.0)then
           call oschydro_event(number,e,n,lam,vx,vy,vz,it,itend)
        end if
       end if 


         call prop3d(elab,mx,my,mz,r,e,p,n,vx,vy,vz,gamma,lam,it)
         


 300     continue
         
         
chp check after the evolution which cells have crossed the hypersurface
         if((CTOption(52).eq.2).and.(mod(it,hypertime).eq.0))then
          dthyper=dt*hypertime   
          dxh=hyperstep*dx 
          do k=hyperstep,ngr-hyperstep,hyperstep
           do j=hyperstep,ngr-hyperstep,hyperstep
            do i=hyperstep,ngr-hyperstep,hyperstep                        
              over=0
              under=0  
              if(elab(i,j,k).gt.1.d-12)then
                if(e(i,j,k).gt.CTParam(64)
     &             .or.e(i,j,k+hyperstep).gt.CTParam(64)
     &             .or.e(i,j+hyperstep,k).gt.CTParam(64)
     &             .or.e(i,j+hyperstep,k+hyperstep).gt.CTParam(64)
     &             .or.e(i+hyperstep,j,k).gt.CTParam(64)
     &             .or.e(i+hyperstep,j,k+hyperstep).gt.CTParam(64)
     &             .or.e(i+hyperstep,j+hyperstep,k).gt.CTParam(64)
     &             .or.e(i+hyperstep,j+hyperstep,k+hyperstep)
     &                 .gt.CTParam(64)
     &             .or.elast(i,j,k).gt.CTParam(64)
     &             .or.elast(i,j,k+hyperstep).gt.CTParam(64)
     &             .or.elast(i,j+hyperstep,k).gt.CTParam(64)
     &             .or.elast(i,j+hyperstep,k+hyperstep).gt.CTParam(64)
     &             .or.elast(i+hyperstep,j,k).gt.CTParam(64)
     &             .or.elast(i+hyperstep,j,k+hyperstep).gt.CTParam(64)
     &             .or.elast(i+hyperstep,j+hyperstep,k).gt.CTParam(64)
     &             .or.elast(i+hyperstep,j+hyperstep,k+hyperstep)
     &                 .gt.CTParam(64))then
                over=1
               end if
               if(e(i,j,k).le.CTParam(64)
     &             .or.e(i,j,k+hyperstep).le.CTParam(64)
     &             .or.e(i,j+hyperstep,k).le.CTParam(64)
     &             .or.e(i,j+hyperstep,k+hyperstep).le.CTParam(64)
     &             .or.e(i+hyperstep,j,k).le.CTParam(64)
     &             .or.e(i+hyperstep,j,k+hyperstep).le.CTParam(64)
     &             .or.e(i+hyperstep,j+hyperstep,k).le.CTParam(64)
     &             .or.e(i+hyperstep,j+hyperstep,k+hyperstep)
     &                 .le.CTParam(64)
     &             .or.elast(i,j,k).le.CTParam(64)
     &             .or.elast(i,j,k+hyperstep).le.CTParam(64)
     &             .or.elast(i,j+hyperstep,k).le.CTParam(64)
     &             .or.elast(i,j+hyperstep,k+hyperstep).le.CTParam(64)
     &             .or.elast(i+hyperstep,j,k).le.CTParam(64)
     &             .or.elast(i+hyperstep,j,k+hyperstep).le.CTParam(64)
     &             .or.elast(i+hyperstep,j+hyperstep,k).le.CTParam(64)
     &             .or.elast(i+hyperstep,j+hyperstep,k+hyperstep)
     &                 .le.CTParam(64))then
                  under=1
               end if
              end if

              if(over.eq.1 .and.under.eq.1)then
 
              call fillcube(i,j,k,hyperstep,e,elast,HyperCube)
                 
              Ngp = COUNT(HyperCube.ge.CTParam(64))

       call Cornelius(CTParam(64),HyperCube,Ngp,dSigma,Nsurf,Vmid,
     &        dthyper,dxh,dxh,dxh,Nambi,Ndisc) 

chp interpolate values on hypersurface and write them to hypersurface array

              call fillcube(i,j,k,hyperstep,n,nlast,Hypern)
              call fillcube(i,j,k,hyperstep,vx,vxlast,Hypervx)
              call fillcube(i,j,k,hyperstep,vy,vylast,Hypervy)
              call fillcube(i,j,k,hyperstep,vz,vzlast,Hypervz)     

              t0=t
              x0=(float(i)-float(ngr)*.5-.5)*sngl(dx)
              y0=(float(j)-float(ngr)*.5-.5)*sngl(dx)   
              z0=(float(k)-float(ngr)*.5-.5)*sngl(dx) 
              
              DO ii=1,Nsurf
               Hypert=t0+Vmid(0,ii)  
               Hyperx=x0+Vmid(1,ii)
               Hypery=y0+Vmid(2,ii)
               Hyperz=z0+Vmid(3,ii) 

            
         emid= Quadrilinear(Vmid(:,ii),HyperCube,dthyper,dxh,dxh,dxh)    
         nmid = Quadrilinear(Vmid(:,ii),Hypern,dthyper,dxh,dxh,dxh)
         Vxmid = Quadrilinear(Vmid(:,ii),Hypervx,dthyper,dxh,dxh,dxh) 
         Vymid = Quadrilinear(Vmid(:,ii),Hypervy,dthyper,dxh,dxh,dxh)
         Vzmid = Quadrilinear(Vmid(:,ii),Hypervz,dthyper,dxh,dxh,dxh)
         
               Tmid=Temp(emid,nmid)
               muqmid=chem(emid,nmid) 
               musmid=schem(emid,nmid)
               pressmid=press(emid,nmid)
c calculate net baryon density flow through surface
               gammid=1.d0/sqrt(1.d0-Vxmid**2-Vymid**2-Vzmid**2)
               umudsigmu=gammid*(dSigma(0,ii)+Vxmid*dSigma(1,ii)
     &             +Vymid*dSigma(2,ii)+Vzmid*dSigma(3,ii))
               bflow=nmid*umudsigmu*n0
            tnn=(emid+pressmid)*gammid**2-pressmid
            tnx=(emid+pressmid)*Vxmid*gammid**2 
            tny=(emid+pressmid)*Vymid*gammid**2 
            tnz=(emid+pressmid)*Vzmid*gammid**2

            ene=(dSigma(0,ii)*tnn+dSigma(1,ii)*tnx
     &         +dSigma(2,ii)*tny+dSigma(3,ii)*tnz)*e0/1.d3

             if((bflow).lt.0)then
              negbflow=negbflow+bflow
             else
              posbflow=posbflow+bflow
             endif    

             if((ene).lt.0)then
              negene=negene+ene
             else
              posene=posene+ene
             endif    


               tf(mfr)=Hypert
               ifr(mfr)=Hyperx
               jfr(mfr)=Hypery
               kfr(mfr)=Hyperz
               tfr(mfr)=Tmid
               muqfr(mfr)=muqmid
               musfr(mfr)=musmid
               vxfr(mfr)= Vxmid
               vyfr(mfr)= Vymid
               vzfr(mfr)= Vzmid
               dstfr(mfr)=dSigma(0,ii)
               dsxfr(mfr)=dSigma(1,ii)
               dsyfr(mfr)=dSigma(2,ii)
               dszfr(mfr)=dSigma(3,ii)  
               gfr(mfr)= 1.d0/sqrt(1-Vxmid**2-Vymid**2-Vzmid**2)
               dxfr(mfr)=dxh
             
               mmax=mfr
               efr=efr+elab(i,j,k)*vol*e0/1000d0
               rfr=rfr+r(i,j,k)*vol*n0
               sfr=sfr+entro(e(i,j,k),n(i,j,k))*vol*n0
     $                          *gamma(i,j,k)
             if(CTOption(54).eq.2)then 
          call oschyper(number,tf(mfr),ifr(mfr),jfr(mfr),kfr(mfr),
     &         e(i,j,k),press(e(i,j,k),n(i,j,k)),n(i,j,k),vxfr(mfr),
     &         vyfr(mfr),vzfr(mfr),tfr(mfr),muqfr(mfr),musfr(mfr),
     &         dstfr(mfr),dsxfr(mfr),dsyfr(mfr),dszfr(mfr))
           end if
                mfr=mfr+1 
              END DO
                
             endif  
            end do
           end do
          end do  
          endif

         if (stophydro.eq.0.or.testflag.eq.1.or.it.eq.itmax)then
            if(it.eq.itmax.and.testflag.eq.0)then
               write(*,*) 'maximum number of timesteps reached:',itmax
               write(*,*) 'before freezeout criterium is fulfilled'
            end if  
            etothy=0.0d0
            pxtothy=0.0d0
            pytothy=0.0d0
            pztothy=0.0d0
            bartothy=0.0d0
            
            if(CTOption(52).eq.1)then
               mfr=0
               efr=0d0
               rfr=0d0
               sfr=0d0
            endif
            
            do k = 1,ngr
               do j = 1,ngr
                  do i = 1,ngr                 
                     if (elab(i,j,k).gt.1d-12)then
                        
                        etothy=etothy+elab(i,j,k)
                        pxtothy=pxtothy+abs(mx(i,j,k))
                        pytothy=pytothy+abs(my(i,j,k))
                        pztothy=pztothy+abs(mz(i,j,k))
                        bartothy=bartothy+r(i,j,k)
                        entrotote=entrotote+entro(e(i,j,k),n(i,j,k))
     $                       *n0*vol*gamma(i,j,k)
                        
                        if((CTOption(52).eq.1)
     $                         .or.
     $                      ((CTOption(52).eq.0).and.(frozen(k).eq.0))
     $                     )then  
                           tf(mfr)=t-dt
                           ifr(mfr)=(float(i)-float(ngr)*.5-.5)*sngl(dx)
                           jfr(mfr)=(float(j)-float(ngr)*.5-.5)*sngl(dx)
                           kfr(mfr)=(float(k)-float(ngr)*.5-.5)*sngl(dx)
                           dstfr(mfr)=dx*dx*dx
                           dsxfr(mfr)=0.0d0
                           dsyfr(mfr)=0.0d0
                           dszfr(mfr)=0.0d0
                           tfr(mfr)=Temp(e(i,j,k),n(i,j,k))
                           muqfr(mfr)=chem(e(i,j,k),n(i,j,k))
                           musfr(mfr)=schem(e(i,j,k),n(i,j,k))
                           vxfr(mfr)= vx(i,j,k)
                           vyfr(mfr)= vy(i,j,k)
                           vzfr(mfr)= vz(i,j,k)
                           gfr(mfr)= gamma(i,j,k)
                           dxfr(mfr)=dx
                           if(CTOption(54).eq.2)then 
                             call oschyper(number,tf(mfr),ifr(mfr),
     &                               jfr(mfr),kfr(mfr),e(i,j,k),
     &                               press(e(i,j,k),n(i,j,k)),n(i,j,k),
     &                               vxfr(mfr),vyfr(mfr),vzfr(mfr),
     &                               tfr(mfr),muqfr(mfr),musfr(mfr),
     &                               dstfr(mfr),dsxfr(mfr),dsyfr(mfr),
     &                               dszfr(mfr))
                           end if
                           mfr=mfr+1
                           mmax=mfr
                           efr=efr+elab(i,j,k)*vol*e0/1000d0
                           rfr=rfr+r(i,j,k)*vol*n0
                           sfr=sfr+entro(e(i,j,k),n(i,j,k))*vol*n0
     $                          *gamma(i,j,k)
                        endif
                     endif
                  end do
               end do
            end do
            
            write(*,301) 'hydroevolution','start','end','diff'
            write(*,302) 'energy:           ',h_etot,etothy*e0*vol/1.d3,
     &           ((etothy*e0*vol/1.d3)-h_etot)/h_etot*100
            write(*,302) 'entropy:          ', entrotota,entrotote,
     &           -(entrotota-entrotote)/entrotota*100
            write(*,302) 'baryons:          ',h_btot,
     &          bartothy*n0*vol,(bartothy*n0*vol-h_btot)/h_btot*100
            write(*,301) 'iso-eps-hyper','total','pos','neg'
            write(*,303) 'energy:         ',posene+negene,posene,negene
            write(*,303) 'netbflow:       ',posbflow+negbflow,
     &         posbflow,negbflow

            if (abs(((etothy*e0*vol/1.d3)-h_etot)/h_etot).gt.0.1d0)then
               write(*,*) 'energy difference is greater than 10 %'
c               write(*,*) 'evolution is stopped'
            end if   
            if(abs((entrotota-entrotote)/entrotota).gt.0.1d0)then
               write(*,*)'entropy difference greater than 10%!!'
               stop 137
            endif
 301  format(a14,x,a5,x,a7,x,a6,x,'(conservation check)')
 302  format(a12,x,f7.1,x,f7.1,x,f5.1,'%')
 303  format(a12,x,f7.1,x,f7.1,x,f7.1)
c            write(*,*) 'total energy on hypersurface', efr
c            write(*,*) 'total baryonnumber on hypersurface', rfr
c            write(*,*) 'total enropy on hypersurface', sfr
            return
         end if   
 99   continue
      return

      
      end


c.......................................................................

      subroutine oschydro_header(thydro_start,itmax,number)

      implicit none

      include 'defs.f'      
      include 'options.f'

      integer Ap, At, Zp, Zt, npart, nbar, nmes, ctag
      integer nsteps,ranseed,event,eos,dectag,uid_cnt
      integer NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm
c 7 integer

      common /sys/ npart, nbar, nmes, ctag,nsteps,uid_cnt,
     +             ranseed,event,Ap,At,Zp,Zt,eos,dectag,
     +             NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm


      character*12 d12_1,d12_2,d12_3
      character*80 outputline
      real*8 e0,n0,Bag,hc3,pi2
      real*8 dx,dt,t,vol,cut  
      real*8 thydro_start
      integer itmax
      integer number

      common /grstate/ e0,n0,Bag,hc3,pi2
      common /gitter/ dx,dt,t,vol,cut

chp header for OSCAR hydro output
 920  format(3a12)
 921  format(a80)
 922  format(i6,6i4)
 923  format(8f10.3) 
 924  format(a,f5.1,a)
 925  format(a,f9.4,a,f9.4,a,f9.4)
 926  format(a,f9.4,a,f9.4)
 927  format(a,f9.4,a,f9.4,a,f9.4,a,f9.4)
 928  format(a,2i5,a,2i5,a,f8.2,a)
        
         d12_1='OSCAR2008H '
         d12_2='ideal      '
        if(CTOption(54).eq.1)then
         d12_3='history    '
        elseif(CTOption(54).eq.2)then
         d12_3='final_hs'
        endif  
      
         write(number,920)d12_1,d12_2,d12_3
c INIT
        write(outputline,928) 'INIT: UrQMD,',Ap,Zp,' on',At,Zt,
     &              ' at E_CM= ',ecm,' AGEV'
     
        write(number,921) outputline

         write(outputline,925)'INIT: b= ',bimp,' n_B,0= ',n0,' e_0= ',
     &                       e0*1.d-3
         write(number,921) outputline

        if(CTOption(52).eq.0)then   
         write(outputline,926)'INIT: start time',thydro_start,
     &     ' grad FO crit.',CTParam(64)*e0*1.d-3
        elseif(CTOption(52).eq.1)then
         write(outputline,926)'INIT: start time',thydro_start,
     &     ' isochronous FO crit. ',CTParam(64)*e0*1.d-3
        elseif(CTOption(52).eq.2)then
         write(outputline,926)'INIT: start time',thydro_start,
     &     ' iso-energy density FO crit.',CTParam(64)*e0*1.d-3
   
        endif
        

         write(number,921) outputline
          

c EOS
         if(CTOption(47).eq.3) then
            outputline='EOS: Bag Model QGP (B=(235*1.d-3)**4) + HRG'
         elseif(CTOption(47).eq.2) then
            outputline='EOS: Hadron Gas, UrQMD degrees of freedom'
         elseif(CTOption(47).eq.5) then
            outputline='EOS: Lattice Eos from chiral model'  
         endif

         write(number,921) outputline
c CHARGES
         outputline='CHARGES: baryon'
         write(number,921) outputline
         
         if(CTOption(54).eq.1) then
           outputline='HYPER: full evolution'
         elseif(CTOption(54).eq.2) then
           if(CTOption(52).eq.0) then
             write(outputline,926)'HYPER: grad FO crit.',
     &                             CTParam(64)*e0*1.d-3
           elseif(CTOption(52).eq.1) then
             write(outputline,926)'HYPER: isochronous FO crit. ',
     &                             CTParam(64)*e0*1.d-3
           elseif(CTOption(52).eq.2) then
             write(outputline,926)'HYPER: iso-energy density FO crit.',
     &                             CTParam(64)*e0*1.d-3
           endif
         endif
         write(number,921) outputline

c GEOM
         outputline='GEOM: 3d-cart'
         write(number,921) outputline
c GRID
         outputline='GRID: Euler'
         write(number,921) outputline
         write(number,922) itmax,int(ngr),int(ngr),
     &    int(ngr),1,0,0
         write(number,923) thydro_start+t,thydro_start+t+itmax*dt,
     &                 (float(0)-float(ngr)*.5-.5)*sngl(dx),
     &                 (float(ngr)-float(ngr)*.5-.5)*sngl(dx),
     &                 (float(0)-float(ngr)*.5-.5)*sngl(dx),
     &                 (float(ngr)-float(ngr)*.5-.5)*sngl(dx),
     &                 (float(0)-float(ngr)*.5-.5)*sngl(dx),
     &                 (float(ngr)-float(ngr)*.5-.5)*sngl(dx)


c VISCOSITY
         outputline='VISCOSITY: none'
         write(number,921) outputline

c COMM
         outputline='COMM: SHASTA by D. Rischke' 
         write(number,921) outputline
         outputline='COMM: modified by H. Petersen & M. Bleicher'
         write(number,921) outputline
         outputline='COMM: axis order, z-y-x' 
         write(number,921) outputline
c END
         outputline='END_OF_HEADER'
         write(number,921) outputline


      return
      end

c.......................................................................
      subroutine oschydro_event(number,e,n,lam,vx,vy,vz,it,itend)
          
      implicit none      
    
      include 'defs.f'

      integer number,i,j,k,ncells,it,itend
      real*8 e(ngr,ngr,ngr),p(ngr,ngr,ngr),n(ngr,ngr,ngr)
      real*8 vx(ngr,ngr,ngr),vy(ngr,ngr,ngr),vz(ngr,ngr,ngr)
      real*8 lam(ngr,ngr,ngr)
      real*8 temp,press,chem
      real*8 e0,n0,Bag,hc3,pi2

      common /grstate/ e0,n0,Bag,hc3,pi2
chp output for visualization
 995  format(4i5,2F11.4,2F9.4,3E24.16,2F11.6)
 996  format(i5,i10)     
       ncells=0 

       do i=1,ngr
        do j=1,ngr
         do k=1,ngr
           write(number,995) it-itend+1,i,j,k,e(i,j,k)*e0*1.d-3,
     &       press(e(i,j,k),n(i,j,k))*1.d-3,
     &       temp(e(i,j,k),n(i,j,k))*1.d-3,        
     &       lam(i,j,k),vx(i,j,k),vy(i,j,k),
     &       0.5d0*log((1+vz(i,j,k))/(1-vz(i,j,k))),
     &       n(i,j,k)*n0, chem(e(i,j,k),n(i,j,k))*1.d-3   
 
         end do
        end do
       end do

      return
      end 
c..................................................................      


       subroutine oschyper(number,t,x,y,z,e,p,n,vx,vy,vz,
     &                       temp,chem,schem,dst,dsx,dsy,dsz)

       integer number
       real*4 t,x,y,z,temp,chem,schem
       real*8 dst,dsx,dsy,dsz
       real*8 vx,vy,vz,e,p,n

       real*8 e0,n0,Bag,hc3,pi2
       common /grstate/ e0,n0,Bag,hc3,pi2

 997   format(6F11.4,2F11.6,3E24.16,3F11.6,4F11.6)    
     

        write(number,997) t,x,y,z,e*e0*1.d-3,p*e0*1.d-3,
     &       temp*1.d-3,0.0d0,vx,vy,vz,n*n0,chem*1.d-3,schem*1.d-3,
     &       dst,dsx,dsy,dsz
    
      return
      end 
      
      
      
c-----------------------------------------------------------------

      subroutine fillcube(i,j,k,hyperstep,e,elast,HyperCube)
      include 'defs.f'      
 
      real*8 e(ngr,ngr,ngr),elast(ngr,ngr,ngr)
      REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1,0:1) :: HyperCube
      integer hyperstep  

      HyperCube(0,0,0,0)=elast(i,j,k)
      HyperCube(0,0,0,1)=elast(i,j,k+hyperstep)
      HyperCube(0,0,1,0)=elast(i,j+hyperstep,k) 
      HyperCube(0,0,1,1)=elast(i,j+hyperstep,k+hyperstep) 
      HyperCube(0,1,0,0)=elast(i+hyperstep,j,k) 
      HyperCube(0,1,0,1)=elast(i+hyperstep,j,k+hyperstep) 
      HyperCube(0,1,1,0)=elast(i+hyperstep,j+hyperstep,k) 
      HyperCube(0,1,1,1)=elast(i+hyperstep,j+hyperstep,k+hyperstep) 
      HyperCube(1,0,0,0)=e(i,j,k)
      HyperCube(1,0,0,1)=e(i,j,k+hyperstep)
      HyperCube(1,0,1,0)=e(i,j+hyperstep,k) 
      HyperCube(1,0,1,1)=e(i,j+hyperstep,k+hyperstep) 
      HyperCube(1,1,0,0)=e(i+hyperstep,j,k) 
      HyperCube(1,1,0,1)=e(i+hyperstep,j,k+hyperstep) 
      HyperCube(1,1,1,0)=e(i+hyperstep,j+hyperstep,k) 
      HyperCube(1,1,1,1)=e(i+hyperstep,j+hyperstep,k+hyperstep) 




      return 
      end

c----------------------------------------------------------------------

c------------------------------------------------------------------------
c     JS
c     Subroutine which changes the grid spacing when Energy conservation
c     is violated
c
      subroutine changegrid(e,p,n,vx,vy,vz,gamma,efr,rfr,sfr,frozen,
     $     frozenb,mfr,etot1,ifr,jfr,kfr,tf,dstfr,dsxfr,dsyfr,dszfr,
     $     tfr,muqfr,musfr,vxfr,vyfr,vzfr,mmax,elast,vxlast,vylast,
     $     vzlast,nlast)
      INCLUDE 'defs.f'
      

      include 'options.f' 

  

      real*8 e0,n0,Bag,hc3,pi2
      real*8 dx,dt,t,vol,cut
      real*4 dxfr(ngr*ngr*ngr)
      real*8 e(ngr,ngr,ngr),p(ngr,ngr,ngr),n(ngr,ngr,ngr)
      real*8 vx(ngr,ngr,ngr),vy(ngr,ngr,ngr),vz(ngr,ngr,ngr)
      real*8 gamma(ngr,ngr,ngr)
      real*8 elab(ngr,ngr,ngr),r(ngr,ngr,ngr)
      real*8 mx(ngr,ngr,ngr),my(ngr,ngr,ngr),mz(ngr,ngr,ngr)
      real*8 elabsmall(ngr/2,ngr/2,ngr/2),rsmall(ngr/2,ngr/2,ngr/2)
      real*8 mxsmall(ngr/2,ngr/2,ngr/2),mysmall(ngr/2,ngr/2,ngr/2)
      real*8 mzsmall(ngr/2,ngr/2,ngr/2)
      integer frozensmall(ngr/2),frozenbsmall(ngr/2)

      real*8 entro,temp,chem,press,schem,velo


chp properties of last timestep for hypersurface finding
      real*8  elast(ngr,ngr,ngr),vxlast(ngr,ngr,ngr),vylast(ngr,ngr,ngr)
      real*8  vzlast(ngr,ngr,ngr),nlast(ngr,ngr,ngr)
      real*8  elastsm(ngr/2,ngr/2,ngr/2),vxlastsm(ngr/2,ngr/2,ngr/2)
      real*8  vylastsm(ngr/2,ngr/2,ngr/2)
      real*8  vzlastsm(ngr/2,ngr/2,ngr/2),nlastsm(ngr/2,ngr/2,ngr/2)


      real*4 ifr(ngr*ngr*ngr),jfr(ngr*ngr*ngr),kfr(ngr*ngr*ngr)
      integer frozen(ngr),frozenb(ngr),mmax,mfr
      real*4 tf(ngr*ngr*ngr),tfr(ngr*ngr*ngr),muqfr(ngr*ngr*ngr)
      real*4 musfr(ngr*ngr*ngr)
      real*8 vyfr(ngr*ngr*ngr),vxfr(ngr*ngr*ngr)
      real*8 vzfr(ngr*ngr*ngr),gfr(ngr*ngr*ngr)
      real*8 sfr,efr,rfr
      real*8 dstfr(ngr*ngr*ngr),dsxfr(ngr*ngr*ngr)
      real*8 dsyfr(ngr*ngr*ngr),dszfr(ngr*ngr*ngr)       

 
      common /onefgrid/ elab,mx,my,mz,r      
      common /gitter/ dx,dt,t,vol,cut
      common /grstate/ e0,n0,Bag,hc3,pi2



      integer i,j,k
      real*8 m,etot1

c     if energy loss is greater than 0.1 % then increase dx
      
      if(CTOption(52).lt.2)then   
      do k=1,ngr-1
         if((frozen(k)*frozen(k+1)*2.ne.frozen(k)+frozen(k+1))
     $        .and.(k/2.ne.(k+1)/2)) then
            
            do 4583 j=1,ngr
               do 4593 i=1,ngr
                  if (elab(i,j,k).gt.1d-8)then
                     if (frozen(k).eq.1) then                            
c                        write(*,*) "Freezing out plane Nr.: ", k+1  
                        tf(mfr)=t-dt
                        ifr(mfr)=(float(i)-float(ngr)*.5-.5)*sngl(dx)
                        jfr(mfr)=(float(j)-float(ngr)*.5-.5)*sngl(dx)
                        kfr(mfr)=(float(k+1)-float(ngr)*.5-.5)*sngl(dx)
                        dstfr(mfr)=dx*dx*dx
                        dsxfr(mfr)=0.0d0
                        dsyfr(mfr)=0.0d0
                        dszfr(mfr)=0.0d0  
                        tfr(mfr)=Temp(e(i,j,k+1),n(i,j,k+1))
                        muqfr(mfr)=chem(e(i,j,k+1),n(i,j,k+1))
                        musfr(mfr)=schem(e(i,j,k+1),n(i,j,k+1))
                        vxfr(mfr)= vx(i,j,k+1)
                        vyfr(mfr)= vy(i,j,k+1)
                        vzfr(mfr)= vz(i,j,k+1)
                        gfr(mfr)= gamma(i,j,k+1)
                        dxfr(mfr)=dx
                        mfr=mfr+1
                        mmax=mfr
                        frozen(k+1)=1
                        efr=efr+elab(i,j,k+1)*vol*e0/1000d0
                        rfr=rfr+r(i,j,k+1)*vol*n0
                        sfr=sfr+entro(e(i,j,k+1),n(i,j,k+1))*vol*n0
     $                       *gamma(i,j,k+1)
                     else
c     write(*,*) "Freezing out plane Nr.: ", k 
                        tf(mfr)=t-dt
                        ifr(mfr)=(float(i)-float(ngr)*.5-.5)*sngl(dx)
                        jfr(mfr)=(float(j)-float(ngr)*.5-.5)*sngl(dx)
                        kfr(mfr)=(float(k)-float(ngr)*.5-.5)*sngl(dx)
                        dstfr(mfr)=dx*dx*dx
                        dsxfr(mfr)=0.0d0
                        dsyfr(mfr)=0.0d0
                        dszfr(mfr)=0.0d0  
                        tfr(mfr)=Temp(e(i,j,k),n(i,j,k))
                        muqfr(mfr)=chem(e(i,j,k),n(i,j,k))
                        musfr(mfr)=schem(e(i,j,k),n(i,j,k))
                        vxfr(mfr)= vx(i,j,k)
                        vyfr(mfr)= vy(i,j,k)
                        vzfr(mfr)= vz(i,j,k)
                        gfr(mfr)= gamma(i,j,k)
                        dxfr(mfr)=dx
                        mfr=mfr+1
                        mmax=mfr
                        frozen(k)=1
                        efr=efr+elab(i,j,k)*vol*e0/1000d0
                        rfr=rfr+r(i,j,k)*vol*n0
                        sfr=sfr+entro(e(i,j,k),n(i,j,k))*vol*n0
     $                       *gamma(i,j,k)
                     endif
                  endif
 4593          continue
 4583       continue
         endif
      end do
      end if
      
      
      dx= 2d0*dx
      dt=0.4d0*dx
      vol = dx*dx*dx
      etot1=0d0
c     now copy everything on a smaller grid
      
      do 4631 k = 1,ngr/2
         frozensmall(k)=frozen(2*k)
         frozenbsmall(k)=frozenb(2*k)
         do 4621 j = 1,ngr/2
            do 4611 i = 1,ngr/2
               
c     edensity          
               elabsmall(i,j,k)=(elab(2*i-1,2*j-1,2*k-1)+
     $              elab(2*i,2*j-1,2*k-1)+elab(2*i-1,2*j,2*k-1)+
     $              elab(2*i-1,2*j-1,2*k)+elab(2*i,2*j,2*k-1)+
     $              elab(2*i,2*j-1,2*k)+elab(2*i-1,2*j,2*k)+
     $              elab(2*i,2*j,2*k))/8d0
c     barydensity
               rsmall(i,j,k)=(r(2*i-1,2*j-1,2*k-1)+
     $              r(2*i,2*j-1,2*k-1)+r(2*i-1,2*j,2*k-1)+
     $              r(2*i-1,2*j-1,2*k)+r(2*i,2*j,2*k-1)+
     $              r(2*i,2*j-1,2*k)+r(2*i-1,2*j,2*k)+
     $              r(2*i,2*j,2*k))/8d0
               
               
c     momdensities
               mxsmall(i,j,k)=(mx(2*i-1,2*j-1,2*k-1)+
     $              mx(2*i,2*j-1,2*k-1)+mx(2*i-1,2*j,2*k-1)+
     $              mx(2*i-1,2*j-1,2*k)+mx(2*i,2*j,2*k-1)+
     $              mx(2*i,2*j-1,2*k)+mx(2*i-1,2*j,2*k)+
     $              mx(2*i,2*j,2*k))/8d0
               
               mysmall(i,j,k)=(my(2*i-1,2*j-1,2*k-1)+
     $              my(2*i,2*j-1,2*k-1)+my(2*i-1,2*j,2*k-1)+
     $              my(2*i-1,2*j-1,2*k)+my(2*i,2*j,2*k-1)+
     $              my(2*i,2*j-1,2*k)+my(2*i-1,2*j,2*k)+
     $              my(2*i,2*j,2*k))/8d0
               
               mzsmall(i,j,k)=(mz(2*i-1,2*j-1,2*k-1)+
     $              mz(2*i,2*j-1,2*k-1)+mz(2*i-1,2*j,2*k-1)+
     $              mz(2*i-1,2*j-1,2*k)+mz(2*i,2*j,2*k-1)+
     $              mz(2*i,2*j-1,2*k)+mz(2*i-1,2*j,2*k)+
     $              mz(2*i,2*j,2*k))/8d0
 
chp hypersurface quantities
               elastsm(i,j,k)=(elast(2*i-1,2*j-1,2*k-1)+
     $              elast(2*i,2*j-1,2*k-1)+elast(2*i-1,2*j,2*k-1)+
     $              elast(2*i-1,2*j-1,2*k)+elast(2*i,2*j,2*k-1)+
     $              elast(2*i,2*j-1,2*k)+elast(2*i-1,2*j,2*k)+
     $              elast(2*i,2*j,2*k))/8d0           
               vxlastsm(i,j,k)=(vxlast(2*i-1,2*j-1,2*k-1)+
     $              vxlast(2*i,2*j-1,2*k-1)+vxlast(2*i-1,2*j,2*k-1)+
     $              vxlast(2*i-1,2*j-1,2*k)+vxlast(2*i,2*j,2*k-1)+
     $              vxlast(2*i,2*j-1,2*k)+vxlast(2*i-1,2*j,2*k)+
     $              vxlast(2*i,2*j,2*k))/8d0           
               vylastsm(i,j,k)=(vylast(2*i-1,2*j-1,2*k-1)+
     $              vylast(2*i,2*j-1,2*k-1)+vylast(2*i-1,2*j,2*k-1)+
     $              vylast(2*i-1,2*j-1,2*k)+vylast(2*i,2*j,2*k-1)+
     $              vylast(2*i,2*j-1,2*k)+vylast(2*i-1,2*j,2*k)+
     $              vylast(2*i,2*j,2*k))/8d0           
               vzlastsm(i,j,k)=(vzlast(2*i-1,2*j-1,2*k-1)+
     $              vzlast(2*i,2*j-1,2*k-1)+vzlast(2*i-1,2*j,2*k-1)+
     $              vzlast(2*i-1,2*j-1,2*k)+vzlast(2*i,2*j,2*k-1)+
     $              vzlast(2*i,2*j-1,2*k)+vzlast(2*i-1,2*j,2*k)+
     $              vzlast(2*i,2*j,2*k))/8d0           
               nlastsm(i,j,k)=(nlast(2*i-1,2*j-1,2*k-1)+
     $              nlast(2*i,2*j-1,2*k-1)+nlast(2*i-1,2*j,2*k-1)+
     $              nlast(2*i-1,2*j-1,2*k)+nlast(2*i,2*j,2*k-1)+
     $              nlast(2*i,2*j-1,2*k)+nlast(2*i-1,2*j,2*k)+
     $              nlast(2*i,2*j,2*k))/8d0           
                  
               
                     
 4611       continue
 4621    continue
 4631 continue
      
c     Now empty the big grid
      
      do 5631 k = 1,ngr
         frozen(k)=0
         frozenb(k)=0
         do 5621 j = 1,ngr
            do 5611 i = 1,ngr
               elab(i,j,k)=0d0
               r(i,j,k)=0d0
               mx(i,j,k)=0d0
               my(i,j,k)=0d0
               mz(i,j,k)=0d0
               elast(i,j,k)=0.d0
               vxlast(i,j,k)=0.d0
               vylast(i,j,k)=0.d0
               vzlast(i,j,k)=0.d0
               nlast(i,j,k)=0.d0
 5611       continue
 5621    continue
 5631 continue
      

c     then put the stuff back in
      
      do 5831 k = 1,ngr/2
         frozen(k+ngr/4)=frozensmall(k)
         frozenb(k+ngr/4)=frozenbsmall(k)
         do 5821 j = 1,ngr/2
            do 5811 i = 1,ngr/2      
               elab(i+ngr/4,j+ngr/4,k+ngr/4)=elabsmall(i,j,k)
               r(i+ngr/4,j+ngr/4,k+ngr/4)=rsmall(i,j,k)
               mx(i+ngr/4,j+ngr/4,k+ngr/4)=mxsmall(i,j,k)
               my(i+ngr/4,j+ngr/4,k+ngr/4)=mysmall(i,j,k)
               mz(i+ngr/4,j+ngr/4,k+ngr/4)=mzsmall(i,j,k)
               elast(i+ngr/4,j+ngr/4,k+ngr/4)=elastsm(i,j,k)
               vxlast(i+ngr/4,j+ngr/4,k+ngr/4)=vxlastsm(i,j,k)
               vylast(i+ngr/4,j+ngr/4,k+ngr/4)=vylastsm(i,j,k)
               vzlast(i+ngr/4,j+ngr/4,k+ngr/4)=vzlastsm(i,j,k)
               nlast(i+ngr/4,j+ngr/4,k+ngr/4)=nlastsm(i,j,k)
 5811       continue
 5821    continue
 5831 continue
      
      

      
      
      
c     do an untangle to obtain rf quantities
      
            
      do 243 k = 1,ngr
         do 242 j = 1,ngr
            do 241 i = 1,ngr
               
               m = dsqrt(mx(i,j,k)*mx(i,j,k) + my(i,j,k)*my(i,j,k)
     $              + mz(i,j,k)*mz(i,j,k))
               if (elab(i,j,k).le.1d-12) then
                  elab(i,j,k) = 0d0
                  mx(i,j,k) = 0d0
                  my(i,j,k) = 0d0
                  mz(i,j,k) = 0d0
                  r(i,j,k) = 0d0
                  e(i,j,k) = 0d0
                  vx(i,j,k) = 0d0
                  vy(i,j,k) = 0d0
                  vz(i,j,k) = 0d0
                  n(i,j,k) = 0d0
                  gamma(i,j,k) = 1d0
               else if (m.lt.1d-16) then
                  e(i,j,k) = elab(i,j,k)
                  vx(i,j,k) = 0d0
                  vy(i,j,k) = 0d0
                  vz(i,j,k) = 0d0
                  n(i,j,k) = r(i,j,k)
                  gamma(i,j,k) = 1d0
               else
                  v = velo(elab(i,j,k),m,r(i,j,k))
                  if (v.lt.1d0-1d-16) then
                     e(i,j,k) = elab(i,j,k) - v*m
                     vx(i,j,k) = v*mx(i,j,k)/m
                     vy(i,j,k) = v*my(i,j,k)/m
                     vz(i,j,k) = v*mz(i,j,k)/m
                     gamma(i,j,k) = 1d0/dsqrt(1d0-v*v)
                     n(i,j,k) = r(i,j,k)/gamma(i,j,k)
c     Vacuum is assumed if v .ge. 1 - 10^{-16}.
                  else
                     elab(i,j,k) = 0d0
                     mx(i,j,k) = 0d0
                     my(i,j,k) = 0d0
                     mz(i,j,k) = 0d0
                     r(i,j,k) = 0d0
                     e(i,j,k) = 0d0
                     vx(i,j,k) = 0d0
                     vy(i,j,k) = 0d0
                     vz(i,j,k) = 0d0
                     n(i,j,k) = 0d0
                     gamma(i,j,k) = 1d0
                  end if
               end if
               
               p(i,j,k) = press(e(i,j,k),n(i,j,k))
               etot1=etot1+elab(i,j,k)*vol*e0
 241        continue
 242     continue
 243  continue
      
      
      
      
      return
      end
c     
c------------------------------------------------------------------------
c     
c     subroutine-subprograms which read in EoS matrices.
c     
c     readeos1 reads in pure hadronic EoS, readeos2 reads in the 
c     EoS with phase transition to the QGP.
c     Both are used in the main program.
c
      subroutine readeos1()
      real*8 t,mu,e,p,n,s,mstar,mus,lam
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400),msttab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer in, ie, j
c     
c     Common-Blocks.
c     
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3
      
 


      open(unit=53,
     $     file='eosfiles/hadgas_eos.dat')
      open(unit=54,
     $     file='eosfiles/hg_eos_small.dat')     
      open(unit=55,
     $     file='eosfiles/hg_eos_mini.dat') 

c  Read in EoS.
c
      do 1169 in = 0,400,1
         j = 53
         do 1116 ie = 0,2000,1
            read(j,7787) t,mu,e,p,n,s,mus,lam
            ptab(ie,in) = p
            ttab(ie,in) = t
            mutab(ie,in) = mu
            stab(ie,in) = s
            mustab(ie,in) = mus
            cstab(ie,in) = lam

 1116    continue
 1169 continue
      close(53)

      do 1269 in = 0,200,1
         j = 54
         do 1216 ie = 0,200,1
            read(j,7787) t,mu,e,p,n,s,mus,lam
            ptab2(ie,in) = p
            ttab2(ie,in) = t
            mutab2(ie,in) = mu
            stab2(ie,in) = s
            mustab2(ie,in) = mus
            cstab2(ie,in) = lam

 1216    continue
 1269 continue
      close(54)      

      do 1569 in = 0,200,1
         j = 55
         do 1516 ie = 0,200,1
            read(j,7787) t,mu,e,p,n,s,mus,lam
            ptab3(ie,in) = p
            ttab3(ie,in) = t
            mutab3(ie,in) = mu
            stab3(ie,in) = s
            mustab3(ie,in) = mus
            cstab3(ie,in) = lam
 1516    continue
 1569 continue
      close(55)


      return
 7787 format(2(1x,f8.3),6(1x,e15.7))

      end
c
c-------------------------------------------------------------------
c     
      subroutine readeos2()
      real*8 t,mu,e,p,n,s,mstar,lam,mus
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400),msttab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer in, ie, j
c
c     Common-blocks.
c     

      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3      
      

c     
c     Open files of EoS table.
c     

      open(unit=47,
     $     file='eosfiles/hadgas_eos.dat')
      open(unit=48,
     $     file='eosfiles/hg_eos_small.dat')     
      open(unit=49,
     $     file='eosfiles/hg_eos_mini.dat') 

      open(unit=51,
     $     file='eosfiles/qgpeos/qgpn00.dat')
      open(unit=52,
     $     file='eosfiles/qgpeos/qgpn05.dat')
      open(unit=53,
     $     file='eosfiles/qgpeos/qgpn10.dat')
      open(unit=54,
     $     file='eosfiles/qgpeos/qgpn15.dat')
      open(unit=55,
     $     file='eosfiles/qgpeos/qgpn20.dat')
      open(unit=56,
     $     file='eosfiles/qgpeos/qgpn25.dat')
      open(unit=57,
     $     file='eosfiles/qgpeos/qgpn30.dat')
      open(unit=58,
     $     file='eosfiles/qgpeos/qgpn35.dat')
      open(unit=59,
     $     file='eosfiles/qgpeos/qgpn40.dat')
      open(unit=60,
     $     file='eosfiles/qgpeos/qgpn45.dat')
      open(unit=61,
     $     file='eosfiles/qgpeos/qgpn50.dat')
      open(unit=62,
     $     file='eosfiles/qgpeos/qgpn55.dat')
      open(unit=63,
     $     file='eosfiles/qgpeos/qgpn60.dat')
      open(unit=64,
     $     file='eosfiles/qgpeos/qgpn65.dat')
      open(unit=65,
     $     file='eosfiles/qgpeos/qgpn70.dat')
      open(unit=66,
     $     file='eosfiles/qgpeos/qgpn75.dat')
      open(unit=67,
     $     file='eosfiles/qgpeos/qgpn80.dat')
      open(unit=68,
     $     file='eosfiles/qgpeos/qgpn85.dat')
      open(unit=69,
     $     file='eosfiles/qgpeos/qgpn90.dat')
      open(unit=70,
     $     file='eosfiles/qgpeos/qgpn95.dat')
      open(unit=71,
     $     file='eosfiles/qgpeos/qgpn100.dat')
      open(unit=72,
     $     file='eosfiles/qgpeos/qgpn105.dat')
      open(unit=73,
     $     file='eosfiles/qgpeos/qgpn110.dat')
      open(unit=74,
     $     file='eosfiles/qgpeos/qgpn115.dat')
c
c     Read in EoS.
c
      do 1009 in = 0,239,1
         j = 51 + in/10
         do 1010 ie = 0,200,1
            read(j,7777) t,mu,e,p,n,s,mstar,lam
            ptab(ie,in) = p
            ttab(ie,in) = t
            mutab(ie,in) = mu
            stab(ie,in) = s
            msttab(ie,in) = mstar
            lamtab(ie,in) = lam
 1010   continue
 1009 continue

         do 4669 in = 0,400,1
            j = 47
            do 4616 ie = 0,2000,1
               read(j,7778) t,mu,e,p,n,s,mus,lam
               ttab(ie,in) = t
               mutab(ie,in) = mu
               mustab(ie,in) = mus
 4616       continue
 4669    continue
         
         do 4769 in = 0,200,1
            j = 48
            do 4716 ie = 0,200,1
               read(j,7778) t,mu,e,p,n,s,mus,lam
               ttab2(ie,in) = t
               mutab2(ie,in) = mu
               mustab2(ie,in) = mus
 4716       continue
 4769    continue
 
           do 4869 in = 0,200,1
            j = 49
            do 4816 ie = 0,200,1
               read(j,7778) t,mu,e,p,n,s,mus,lam
               ttab3(ie,in) = t
               mutab3(ie,in) = mu
               mustab3(ie,in) = mus
 4816       continue
 4869    continue
  

      close(47)
      close(48) 
      close(49)

      close(51)
      close(52) 
      close(53)
      close(54)
      close(55)
      close(56)
      close(57)
      close(58)
      close(59)
      close(60)
      close(61)
      close(62)
      close(63)
      close(64)
      close(65)
      close(66)
      close(67)
      close(68)
      close(69)
      close(70)
      close(71)
      close(72)
      close(73)
      close(74)
      return
 7777 format(2(1x,f8.3),4(1x,e15.7),2(1x,f6.3))
 7778 format(2(1x,f8.3),6(1x,e15.7))
      end
c     
c------------------------------------------------------------------------
c

      subroutine readeos3()
      real*8 t,mu,e,p,n,s,mstar,lam,mus
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400),msttab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200) 
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer in, ie, j, eos

c     
c     Common-blocks.
c
      common /eqofst/ eos,stabil,anti
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3
c     
c     Open files of EoS table.
c     
      open(unit=51,
     $     file='eosfiles/chiraleos.dat')
      open(unit=52,
     $     file='eosfiles/chiralsmall.dat')     
      open(unit=56,
     $     file='eosfiles/chiralmini.dat')  






c  Read in EoS.
c
         do 1109 in = 0,400,1
            j = 51
            do 1110 ie = 0,2000,1
            read(j,7777) t,mu,e,p,n,s,mstar,lam
            ptab(ie,in) = p
            ttab(ie,in) = t
            mutab(ie,in) = mu
            stab(ie,in) = s
            msttab(ie,in) = mstar
            cstab(ie,in) = lam
 1110    continue
 1109 continue
      close(51)

      do 1209 in = 0,200,1
         j = 52
         do 1210 ie = 0,200,1
            read(j,7777) t,mu,e,p,n,s,mstar,lam
            ptab2(ie,in) = p
            ttab2(ie,in) = t
            mutab2(ie,in) = mu
            stab2(ie,in) = s
            msttab2(ie,in) = mstar
            cstab2(ie,in) = lam
 1210    continue
 1209 continue
      close(52)

      do 1509 in = 0,200,1
         j = 56
         do 1510 ie = 0,200,1
            read(j,7777) t,mu,e,p,n,s,mstar,lam
            ptab3(ie,in) = p
            ttab3(ie,in) = t
            mutab3(ie,in) = mu
            stab3(ie,in) = s
            msttab3(ie,in) = mstar
            cstab3(ie,in) = lam
 1510    continue
 1509 continue
      close(56)
      
      if(eos.eq.5) then
         open(unit=53,
     $        file='eosfiles/hadgas_eos.dat')
         open(unit=54,
     $        file='eosfiles/hg_eos_small.dat')     
         open(unit=55,
     $        file='eosfiles/hg_eos_mini.dat')     
                
         
         do 1669 in = 0,400,1
            j = 53
            do 1616 ie = 0,2000,1
               read(j,7777) t,mu,e,p,n,s,mus,lam
               ttab(ie,in) = t
               mutab(ie,in) = mu
               mustab(ie,in) = mus
 1616       continue
 1669    continue
         close(53)
         
         do 1769 in = 0,200,1
            j = 54
            do 1716 ie = 0,200,1
               read(j,7777) t,mu,e,p,n,s,mus,lam
               ttab2(ie,in) = t
               mutab2(ie,in) = mu
               mustab2(ie,in) = mus
 1716       continue
 1769    continue
         close(54)
 
           do 1869 in = 0,200,1
            j = 55
            do 1816 ie = 0,200,1
               read(j,7777) t,mu,e,p,n,s,mus,lam
               ttab3(ie,in) = t
               mutab3(ie,in) = mu
               mustab3(ie,in) = mus
 1816       continue
 1869    continue
         close(55)   
      endif
         
      return

 7777 format(2(1x,f8.3),6(1x,e15.7))
      
      end
c     
c------------------------------------------------------------------------
c     
     
      subroutine prop3d(elab,mx,my,mz,r,e,p,n,vx,vy,vz,gamma,lam,it)
      INCLUDE 'defs.f'
c    
c    This subroutine-subprogram does the 3-d propagation of
c    hydrodynamic fields according to the operator splitting
c    procedure.
c     
c    It is used in the main program.
c     
c    Type declarations for variables in common-blocks.
c    
      real*8 e0,n0,B,hc3,pi2
      real*8 dx,dt,t,vol,cut
      integer step
c
c     Type declarations for variables in subroutine prop3d.
c     
c     Rest frame quantities and velocities.
c     
      real*8 e(ngr,ngr,ngr),p(ngr,ngr,ngr),n(ngr,ngr,ngr)
      real*8 vx(ngr,ngr,ngr),vy(ngr,ngr,ngr),vz(ngr,ngr,ngr)
      real*8 gamma(ngr,ngr,ngr),lam(ngr,ngr,ngr)
c     
c     Lab frame quantities.
c
      real*8 elab(ngr,ngr,ngr),r(ngr,ngr,ngr)
      real*8 mx(ngr,ngr,ngr),my(ngr,ngr,ngr),mz(ngr,ngr,ngr)
c     
c     Auxiliary fields.
c
      real*8 elaba(ngr),mxa(ngr),mya(ngr),mza(ngr),ra(ngr)
      real*8 ea(ngr),pa(ngr),na(ngr)
      real*8 vxa(ngr),vya(ngr),vza(ngr),gammaa(ngr)
      real*8 elabh(ngr),mxh(ngr),myh(ngr),mzh(ngr),rh(ngr)
      real*8 s,entro,lambda
      integer it,ii,i,j,k
c     
c     Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
      common /gitter/ dx,dt,t,vol,cut
      common /propag/ step
    
  


c     
c     The propagation order is  xyz, xzy, zxy, zyx, yzx, yxz .
c     
c     (i) xyz:
c     
      if (mod(it,6).eq.1) then
c     
c     For each yz-plane, put all field variables into
c     auxiliary vectors and then call prop.
c     
         do 131 k = 1,ngr
            do 121 j = 1,ngr
               do 111 i = 1,ngr
                  elaba(i) = elab(i,j,k)
                  mxa(i) = mx(i,j,k) 
                  mya(i) = my(i,j,k)
                  mza(i) = mz(i,j,k) 
                  ra(i) = r(i,j,k)
                  ea(i) = e(i,j,k)
                  pa(i) = p(i,j,k)
                  na(i) = n(i,j,k)
                  vxa(i) = vx(i,j,k)
                  vya(i) = vy(i,j,k)
                  vza(i) = vz(i,j,k) 
                  gammaa(i) = gamma(i,j,k)
 111           continue
               if (step.eq.2) then
                  dt = 0.5d0*dt
                  do 1111 i = 1,ngr
                     elabh(i) = elaba(i)
                     mxh(i) = mxa(i)
                     myh(i) = mya(i)
                     mzh(i) = mza(i)
                     rh(i) = ra(i)
 1111             continue
                  call prop(elabh,mxh,myh,mzh,rh,vxa,vya,vza,pa)
                  call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                 vxa,vya,vza,gammaa)
                  dt = 2d0*dt
               end if
               call prop(elaba,mxa,mya,mza,ra,vxa,vya,vza,pa)
               call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1              vxa,vya,vza,gammaa)
c     
c     Update propagated vectors.
c     
               do 1112 i = 1,ngr
                  elab(i,j,k) = elaba(i)
                  mx(i,j,k) = mxa(i)
                  my(i,j,k) = mya(i)
                  mz(i,j,k) = mza(i)
                  r(i,j,k) = ra(i)
                  e(i,j,k) = ea(i)
                  p(i,j,k) = pa(i)
                  n(i,j,k) = na(i)
                  vx(i,j,k) = vxa(i)
                  vy(i,j,k) = vya(i)
                  vz(i,j,k) = vza(i)
                  gamma(i,j,k) = gammaa(i)
 1112       continue
 121     continue
c     
c     Afterwards, propagate in y-direction.
c     
          do 112 i = 1,ngr
             do 122 j = 1,ngr
                elaba(j) = elab(i,j,k)
                mxa(j) = mx(i,j,k) 
                mya(j) = my(i,j,k) 
                mza(j) = mz(i,j,k) 
                ra(j) = r(i,j,k) 
                ea(j) = e(i,j,k)
                pa(j) = p(i,j,k)
                na(j) = n(i,j,k)
                vxa(j) = vx(i,j,k)
                vya(j) = vy(i,j,k)
                vza(j) = vz(i,j,k)
                gammaa(j) = gamma(i,j,k)
 122         continue  
             if (step.eq.2) then
                dt = 0.5d0*dt
                do 1221 j = 1,ngr
                   elabh(j) = elaba(j)
                   mxh(j) = mxa(j)
                   myh(j) = mya(j)
                   mzh(j) = mza(j)
                   rh(j) = ra(j)
 1221           continue
                call prop(elabh,myh,mxh,mzh,rh,vya,vxa,vza,pa)
                call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1             vxa,vya,vza,gammaa)
                dt = 2d0*dt
             end if
             call prop(elaba,mya,mxa,mza,ra,vya,vxa,vza,pa)
             call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c     
c     Update propagated vectors.
c     
            do 1222 j = 1,ngr
               elab(i,j,k) = elaba(j)
               mx(i,j,k) = mxa(j) 
               my(i,j,k) = mya(j)
               mz(i,j,k) = mza(j)
               r(i,j,k) = ra(j)
               e(i,j,k) = ea(j)
               p(i,j,k) = pa(j)
               n(i,j,k) = na(j)
               vx(i,j,k) = vxa(j)
               vy(i,j,k) = vya(j)
               vz(i,j,k) = vza(j)
               gamma(i,j,k) = gammaa(j)
 1222       continue
 112     continue
 131  continue
c     
c     Finally, propagate in z-direction.
c     
      do 123 j = 1,ngr
         do 113 i = 1,ngr
            do 133 k = 1,ngr
               elaba(k) = elab(i,j,k)
               mxa(k) = mx(i,j,k) 
               mya(k) = my(i,j,k) 
               mza(k) = mz(i,j,k) 
               ra(k) = r(i,j,k) 
               ea(k) = e(i,j,k)
               pa(k) = p(i,j,k)
               na(k) = n(i,j,k)
               vxa(k) = vx(i,j,k)
               vya(k) = vy(i,j,k)
               vza(k) = vz(i,j,k)
               gammaa(k) = gamma(i,j,k)
 133        continue  
            if (step.eq.2) then
               dt = 0.5d0*dt
               do 1331 k = 1,ngr
                  elabh(k) = elaba(k)
                  mxh(k) = mxa(k)
                  myh(k) = mya(k)
                  mzh(k) = mza(k)
                  rh(k) = ra(k)
 1331          continue
               call prop(elabh,mzh,mxh,myh,rh,vza,vxa,vya,pa)
               call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1              vxa,vya,vza,gammaa)
               dt = 2d0*dt
            end if
            call prop(elaba,mza,mxa,mya,ra,vza,vxa,vya,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c     
c     Update propagated vectors and calculate total energy, momentum,
c     baryon number, and entropy.
c     
            do 1332 k = 1,ngr
               elab(i,j,k) = elaba(k)
               mx(i,j,k) = mxa(k) 
               my(i,j,k) = mya(k)
               mz(i,j,k) = mza(k)
               r(i,j,k) = ra(k)
               e(i,j,k) = ea(k)
               p(i,j,k) = pa(k)
               n(i,j,k) = na(k)
               vx(i,j,k) = vxa(k)
               vy(i,j,k) = vya(k)
               vz(i,j,k) = vza(k)
               gamma(i,j,k) = gammaa(k)
               if (e(i,j,k).gt.cut) then
                  s = entro(e(i,j,k),n(i,j,k))
               end if
               lam(i,j,k) = lambda(e(i,j,k),n(i,j,k))
 1332       continue
 113     continue
 123  continue
c     
c     (ii) xzy:
c     
      else if (mod(it,6).eq.2) then
c     
c     For each zy-plane, put all field variables into
c     auxiliary vectors and then call prop.
c     
        do 221 j = 1,ngr
           do 231 k = 1,ngr
              do 211 i = 1,ngr
                 elaba(i) = elab(i,j,k)
                 mxa(i) = mx(i,j,k) 
                 mya(i) = my(i,j,k)
                 mza(i) = mz(i,j,k) 
                 ra(i) = r(i,j,k)
                 ea(i) = e(i,j,k)
                 pa(i) = p(i,j,k)
                 na(i) = n(i,j,k)
                 vxa(i) = vx(i,j,k)
                 vya(i) = vy(i,j,k)
                 vza(i) = vz(i,j,k) 
                 gammaa(i) = gamma(i,j,k)
 211          continue
              if (step.eq.2) then
                 dt = 0.5d0*dt
                 do 2111 i = 1,ngr
                    elabh(i) = elaba(i)
                    mxh(i) = mxa(i)
                    myh(i) = mya(i)
                    mzh(i) = mza(i)
                    rh(i) = ra(i)
 2111            continue
                 call prop(elabh,mxh,myh,mzh,rh,vxa,vya,vza,pa)
                 call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
                 dt = 2d0*dt
              end if
              call prop(elaba,mxa,mya,mza,ra,vxa,vya,vza,pa)
              call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1             vxa,vya,vza,gammaa)
c     
c  Update propagated vectors.
c     
              do 2112 i = 1,ngr
                 elab(i,j,k) = elaba(i)
                 mx(i,j,k) = mxa(i)
                 my(i,j,k) = mya(i)
                 mz(i,j,k) = mza(i)
                 r(i,j,k) = ra(i)
                 e(i,j,k) = ea(i)
                 p(i,j,k) = pa(i)
                 n(i,j,k) = na(i)
                 vx(i,j,k) = vxa(i)
                 vy(i,j,k) = vya(i)
                 vz(i,j,k) = vza(i)
                 gamma(i,j,k) = gammaa(i)
 2112         continue
 231       continue
c     
c  Afterwards, propagate in z-direction.
c     
           do 212 i = 1,ngr
              do 232 k = 1,ngr
                 elaba(k) = elab(i,j,k)
                 mxa(k) = mx(i,j,k) 
                 mya(k) = my(i,j,k) 
                 mza(k) = mz(i,j,k) 
                 ra(k) = r(i,j,k) 
                 ea(k) = e(i,j,k)
                 pa(k) = p(i,j,k)
                 na(k) = n(i,j,k)
                 vxa(k) = vx(i,j,k)
                 vya(k) = vy(i,j,k)
                 vza(k) = vz(i,j,k)
                 gammaa(k) = gamma(i,j,k)
 232        continue  
            if (step.eq.2) then
               dt = 0.5d0*dt
               do 2321 k = 1,ngr
                  elabh(k) = elaba(k)
                  mxh(k) = mxa(k)
                  myh(k) = mya(k)
                  mzh(k) = mza(k)
                  rh(k) = ra(k)
 2321          continue
               call prop(elabh,mzh,mxh,myh,rh,vza,vxa,vya,pa)
               call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1              vxa,vya,vza,gammaa)
               dt = 2d0*dt
            end if
            call prop(elaba,mza,mxa,mya,ra,vza,vxa,vya,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1           vxa,vya,vza,gammaa)
c     
c     Update propagated vectors.
c     
            do 2322 k = 1,ngr
               elab(i,j,k) = elaba(k)
               mx(i,j,k) = mxa(k) 
               my(i,j,k) = mya(k)
               mz(i,j,k) = mza(k)
               r(i,j,k) = ra(k)
               e(i,j,k) = ea(k)
               p(i,j,k) = pa(k)
               n(i,j,k) = na(k)
               vx(i,j,k) = vxa(k)
               vy(i,j,k) = vya(k)
               vz(i,j,k) = vza(k)
               gamma(i,j,k) = gammaa(k)
 2322       continue
 212     continue
 221  continue
c     
c     Finally, propagate in y-direction.
c     
      do 233 k = 1,ngr
         do 213 i = 1,ngr
            do 223 j = 1,ngr
               elaba(j) = elab(i,j,k)
               mxa(j) = mx(i,j,k) 
               mya(j) = my(i,j,k) 
               mza(j) = mz(i,j,k) 
               ra(j) = r(i,j,k) 
               ea(j) = e(i,j,k)
               pa(j) = p(i,j,k)
               na(j) = n(i,j,k)
               vxa(j) = vx(i,j,k)
               vya(j) = vy(i,j,k)
               vza(j) = vz(i,j,k)
               gammaa(j) = gamma(i,j,k)
 223        continue  
            if (step.eq.2) then
               dt = 0.5d0*dt
               do 2231 j = 1,ngr
                elabh(j) = elaba(j)
                mxh(j) = mxa(j)
                myh(j) = mya(j)
                mzh(j) = mza(j)
                rh(j) = ra(j)
 2231        continue
             call prop(elabh,myh,mxh,mzh,rh,vya,vxa,vza,pa)
             call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
             dt = 2d0*dt
          end if
          call prop(elaba,mya,mxa,mza,ra,vya,vxa,vza,pa)
          call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c     
c     Update propagated vectors and calculate total energy, momentum,
c     baryon number, and entropy.
c
          do 2232 j = 1,ngr
             elab(i,j,k) = elaba(j)
             mx(i,j,k) = mxa(j) 
             my(i,j,k) = mya(j)
             mz(i,j,k) = mza(j)
             r(i,j,k) = ra(j)
             e(i,j,k) = ea(j)
             p(i,j,k) = pa(j)
             n(i,j,k) = na(j)
             vx(i,j,k) = vxa(j)
             vy(i,j,k) = vya(j)
             vz(i,j,k) = vza(j)
             gamma(i,j,k) = gammaa(j)
            if (e(i,j,k).gt.cut) then
                s = entro(e(i,j,k),n(i,j,k))
             end if
             lam(i,j,k) = lambda(e(i,j,k),n(i,j,k))
 2232     continue
 213   continue
 233  continue
c     
c     (iii) zxy:
c     
      else if (mod(it,6).eq.3) then
c     
c     For each xy-plane, put all field variables into
c     auxiliary vectors and then call prop.
c     
         do 321 j = 1,ngr
            do 311 i = 1,ngr
               do 331 k = 1,ngr
                  elaba(k) = elab(i,j,k)
                  mxa(k) = mx(i,j,k) 
                  mya(k) = my(i,j,k)
                  mza(k) = mz(i,j,k) 
                  ra(k) = r(i,j,k)
                  ea(k) = e(i,j,k)
                  pa(k) = p(i,j,k)
                  na(k) = n(i,j,k)
                  vxa(k) = vx(i,j,k)
                  vya(k) = vy(i,j,k)
                  vza(k) = vz(i,j,k) 
                  gammaa(k) = gamma(i,j,k)
 331           continue
               if (step.eq.2) then
              dt = 0.5d0*dt
              do 3311 k = 1,ngr
                 elabh(k) = elaba(k)
                 mxh(k) = mxa(k)
                 myh(k) = mya(k)
                 mzh(k) = mza(k)
                rh(k) = ra(k)
 3311        continue
             call prop(elabh,mzh,mxh,myh,rh,vza,vxa,vya,pa)
             call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1            vxa,vya,vza,gammaa)
              dt = 2d0*dt
           end if
           call prop(elaba,mza,mxa,mya,ra,vza,vxa,vya,pa)
           call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1          vxa,vya,vza,gammaa)
c     
c     Update propagated vectors.
c     
           do 3312 k = 1,ngr
              elab(i,j,k) = elaba(k)
              mx(i,j,k) = mxa(k)
              my(i,j,k) = mya(k)
              mz(i,j,k) = mza(k)
              r(i,j,k) = ra(k)
              e(i,j,k) = ea(k)
              p(i,j,k) = pa(k)
              n(i,j,k) = na(k)
              vx(i,j,k) = vxa(k)
              vy(i,j,k) = vya(k)
              vz(i,j,k) = vza(k)
              gamma(i,j,k) = gammaa(k)
 3312      continue
 311    continue
c     
c     Afterwards, propagate in x-direction.
c
        do 332 k = 1,ngr
           do 312 i = 1,ngr
              elaba(i) = elab(i,j,k)
              mxa(i) = mx(i,j,k) 
              mya(i) = my(i,j,k) 
              mza(i) = mz(i,j,k) 
              ra(i) = r(i,j,k) 
              ea(i) = e(i,j,k)
              pa(i) = p(i,j,k)
              na(i) = n(i,j,k)
              vxa(i) = vx(i,j,k)
              vya(i) = vy(i,j,k)
              vza(i) = vz(i,j,k)
              gammaa(i) = gamma(i,j,k)
 312        continue  
            if (step.eq.2) then
              dt = 0.5d0*dt
              do 3121 i = 1,ngr
                elabh(i) = elaba(i)
                mxh(i) = mxa(i)
                myh(i) = mya(i)
                mzh(i) = mza(i)
                rh(i) = ra(i)
 3121         continue
              call prop(elabh,mxh,myh,mzh,rh,vxa,vya,vza,pa)
              call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
              dt = 2d0*dt
            end if
            call prop(elaba,mxa,mya,mza,ra,vxa,vya,vza,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c
c  Update propagated vectors.
c
            do 3122 i = 1,ngr
              elab(i,j,k) = elaba(i)
              mx(i,j,k) = mxa(i) 
              my(i,j,k) = mya(i)
              mz(i,j,k) = mza(i)
              r(i,j,k) = ra(i)
              e(i,j,k) = ea(i)
              p(i,j,k) = pa(i)
              n(i,j,k) = na(i)
              vx(i,j,k) = vxa(i)
              vy(i,j,k) = vya(i)
              vz(i,j,k) = vza(i)
              gamma(i,j,k) = gammaa(i)
 3122       continue
 332      continue
 321    continue
c
c  Finally, propagate in y-direction.
c
        do 333 k = 1,ngr
          do 313 i = 1,ngr
            do 323 j = 1,ngr
              elaba(j) = elab(i,j,k)
              mxa(j) = mx(i,j,k) 
              mya(j) = my(i,j,k) 
              mza(j) = mz(i,j,k) 
              ra(j) = r(i,j,k) 
              ea(j) = e(i,j,k)
              pa(j) = p(i,j,k)
              na(j) = n(i,j,k)
              vxa(j) = vx(i,j,k)
              vya(j) = vy(i,j,k)
              vza(j) = vz(i,j,k)
              gammaa(j) = gamma(i,j,k)
 323        continue  
            if (step.eq.2) then
              dt = 0.5d0*dt
              do 3231 j = 1,ngr
                elabh(j) = elaba(j)
                mxh(j) = mxa(j)
                myh(j) = mya(j)
                mzh(j) = mza(j)
                rh(j) = ra(j)
 3231         continue
              call prop(elabh,myh,mxh,mzh,rh,vya,vxa,vza,pa)
              call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
              dt = 2d0*dt
            end if
            call prop(elaba,mya,mxa,mza,ra,vya,vxa,vza,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c
c  Update propagated vectors and calculate total energy, momentum,
c  baryon number, and entropy.
c
            do 3232 j = 1,ngr
              elab(i,j,k) = elaba(j)
              mx(i,j,k) = mxa(j) 
              my(i,j,k) = mya(j)
              mz(i,j,k) = mza(j)
              r(i,j,k) = ra(j)
              e(i,j,k) = ea(j)
              p(i,j,k) = pa(j)
              n(i,j,k) = na(j)
              vx(i,j,k) = vxa(j)
              vy(i,j,k) = vya(j)
              vz(i,j,k) = vza(j)
              gamma(i,j,k) = gammaa(j)
              if (e(i,j,k).gt.cut) then
                s = entro(e(i,j,k),n(i,j,k))
              end if
              lam(i,j,k) = lambda(e(i,j,k),n(i,j,k))
 3232       continue
 313      continue
 333    continue
c
c  (iv) zyx:
c
      else if (mod(it,6).eq.4) then
c
c  For each yx-plane, put all field variables into
c  auxiliary vectors and then call prop.
c
        do 411 i = 1,ngr
          do 421 j = 1,ngr
            do 431 k = 1,ngr
              elaba(k) = elab(i,j,k)
              mxa(k) = mx(i,j,k) 
              mya(k) = my(i,j,k)
              mza(k) = mz(i,j,k) 
              ra(k) = r(i,j,k)
              ea(k) = e(i,j,k)
              pa(k) = p(i,j,k)
              na(k) = n(i,j,k)
              vxa(k) = vx(i,j,k)
              vya(k) = vy(i,j,k)
              vza(k) = vz(i,j,k) 
              gammaa(k) = gamma(i,j,k)
 431        continue
            if (step.eq.2) then
              dt = 0.5d0*dt
              do 4311 k = 1,ngr
                elabh(k) = elaba(k)
                mxh(k) = mxa(k)
                myh(k) = mya(k)
                mzh(k) = mza(k)
                rh(k) = ra(k)
 4311         continue
              call prop(elabh,mzh,mxh,myh,rh,vza,vxa,vya,pa)
              call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
              dt = 2d0*dt
            end if
            call prop(elaba,mza,mxa,mya,ra,vza,vxa,vya,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c
c  Update propagated vectors.
c
            do 4312 k = 1,ngr
              elab(i,j,k) = elaba(k)
              mx(i,j,k) = mxa(k)
              my(i,j,k) = mya(k)
              mz(i,j,k) = mza(k)
              r(i,j,k) = ra(k)
              e(i,j,k) = ea(k)
              p(i,j,k) = pa(k)
              n(i,j,k) = na(k)
              vx(i,j,k) = vxa(k)
              vy(i,j,k) = vya(k)
              vz(i,j,k) = vza(k)
              gamma(i,j,k) = gammaa(k)
 4312       continue
 421      continue
c
c  Afterwards, propagate in y-direction.
c
          do 432 k = 1,ngr
            do 422 j = 1,ngr
              elaba(j) = elab(i,j,k)
              mxa(j) = mx(i,j,k) 
              mya(j) = my(i,j,k) 
              mza(j) = mz(i,j,k) 
              ra(j) = r(i,j,k) 
              ea(j) = e(i,j,k)
              pa(j) = p(i,j,k)
              na(j) = n(i,j,k)
              vxa(j) = vx(i,j,k)
              vya(j) = vy(i,j,k)
              vza(j) = vz(i,j,k)
              gammaa(j) = gamma(i,j,k)
 422        continue  
            if (step.eq.2) then
              dt = 0.5d0*dt
              do 4221 j = 1,ngr
                elabh(j) = elaba(j)
                mxh(j) = mxa(j)
                myh(j) = mya(j)
                mzh(j) = mza(j)
                rh(j) = ra(j)
 4221         continue
              call prop(elabh,myh,mxh,mzh,rh,vya,vxa,vza,pa)
              call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
              dt = 2d0*dt
            end if
            call prop(elaba,mya,mxa,mza,ra,vya,vxa,vza,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c
c  Update propagated vectors.
c
            do 4222 j = 1,ngr
              elab(i,j,k) = elaba(j)
              mx(i,j,k) = mxa(j) 
              my(i,j,k) = mya(j)
              mz(i,j,k) = mza(j)
              r(i,j,k) = ra(j)
              e(i,j,k) = ea(j)
              p(i,j,k) = pa(j)
              n(i,j,k) = na(j)
              vx(i,j,k) = vxa(j)
              vy(i,j,k) = vya(j)
              vz(i,j,k) = vza(j)
              gamma(i,j,k) = gammaa(j)
 4222       continue
 432      continue
 411    continue
c
c  Finally, propagate in x-direction.
c
        do 433 k = 1,ngr
          do 423 j = 1,ngr
            do 413 i = 1,ngr
              elaba(i) = elab(i,j,k)
              mxa(i) = mx(i,j,k) 
              mya(i) = my(i,j,k) 
              mza(i) = mz(i,j,k) 
              ra(i) = r(i,j,k) 
              ea(i) = e(i,j,k)
              pa(i) = p(i,j,k)
              na(i) = n(i,j,k)
              vxa(i) = vx(i,j,k)
              vya(i) = vy(i,j,k)
              vza(i) = vz(i,j,k)
              gammaa(i) = gamma(i,j,k)
 413        continue  
            if (step.eq.2) then
              dt = 0.5d0*dt
              do 4131 i = 1,ngr
                elabh(i) = elaba(i)
                mxh(i) = mxa(i)
                myh(i) = mya(i)
                mzh(i) = mza(i)
                rh(i) = ra(i)
 4131         continue
              call prop(elabh,mxh,myh,mzh,rh,vxa,vya,vza,pa)
              call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
              dt = 2d0*dt
            end if
            call prop(elaba,mxa,mya,mza,ra,vxa,vya,vza,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c
c  Update propagated vectors and calculate total energy, momentum,
c  baryon number, and entropy.
c
            do 4132 i = 1,ngr
              elab(i,j,k) = elaba(i)
              mx(i,j,k) = mxa(i) 
              my(i,j,k) = mya(i)
              mz(i,j,k) = mza(i)
              r(i,j,k) = ra(i)
              e(i,j,k) = ea(i)
              p(i,j,k) = pa(i)
              n(i,j,k) = na(i)
              vx(i,j,k) = vxa(i)
              vy(i,j,k) = vya(i)
              vz(i,j,k) = vza(i)
              gamma(i,j,k) = gammaa(i)
              if (e(i,j,k).gt.cut) then
                s = entro(e(i,j,k),n(i,j,k))
              end if
              lam(i,j,k) = lambda(e(i,j,k),n(i,j,k))
 4132       continue
 423      continue
 433    continue
c
c  (v) yzx:
c
      else if (mod(it,6).eq.5) then
c
c  For each zx-plane, put all field variables into
c  auxiliary vectors and then call prop.
c
        do 511 i = 1,ngr
          do 531 k = 1,ngr
            do 521 j = 1,ngr
              elaba(j) = elab(i,j,k)
              mxa(j) = mx(i,j,k) 
              mya(j) = my(i,j,k)
              mza(j) = mz(i,j,k) 
              ra(j) = r(i,j,k)
              ea(j) = e(i,j,k)
              pa(j) = p(i,j,k)
              na(j) = n(i,j,k)
              vxa(j) = vx(i,j,k)
              vya(j) = vy(i,j,k)
              vza(j) = vz(i,j,k) 
              gammaa(j) = gamma(i,j,k)
 521        continue
            if (step.eq.2) then
              dt = 0.5d0*dt
              do 5211 j = 1,ngr
                elabh(j) = elaba(j)
                mxh(j) = mxa(j)
                myh(j) = mya(j)
                mzh(j) = mza(j)
                rh(j) = ra(j)
 5211         continue
              call prop(elabh,myh,mxh,mzh,rh,vya,vxa,vza,pa)
              call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
              dt = 2d0*dt
            end if
            call prop(elaba,mya,mxa,mza,ra,vya,vxa,vza,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c
c  Update propagated vectors.
c
            do 5212 j = 1,ngr
              elab(i,j,k) = elaba(j)
              mx(i,j,k) = mxa(j)
              my(i,j,k) = mya(j)
              mz(i,j,k) = mza(j)
              r(i,j,k) = ra(j)
              e(i,j,k) = ea(j)
              p(i,j,k) = pa(j)
              n(i,j,k) = na(j)
              vx(i,j,k) = vxa(j)
              vy(i,j,k) = vya(j)
              vz(i,j,k) = vza(j)
              gamma(i,j,k) = gammaa(j)
 5212       continue
 531      continue
c
c  Afterwards, propagate in z-direction.
c
          do 522 j = 1,ngr
            do 532 k = 1,ngr
              elaba(k) = elab(i,j,k)
              mxa(k) = mx(i,j,k) 
              mya(k) = my(i,j,k) 
              mza(k) = mz(i,j,k) 
              ra(k) = r(i,j,k) 
              ea(k) = e(i,j,k)
              pa(k) = p(i,j,k)
              na(k) = n(i,j,k)
              vxa(k) = vx(i,j,k)
              vya(k) = vy(i,j,k)
              vza(k) = vz(i,j,k)
              gammaa(k) = gamma(i,j,k)
 532        continue  
            if (step.eq.2) then
              dt = 0.5d0*dt
              do 5321 k = 1,ngr
                elabh(k) = elaba(k)
                mxh(k) = mxa(k)
                myh(k) = mya(k)
                mzh(k) = mza(k)
                rh(k) = ra(k)
 5321         continue
              call prop(elabh,mzh,mxh,myh,rh,vza,vxa,vya,pa)
              call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
              dt = 2d0*dt
            end if
            call prop(elaba,mza,mxa,mya,ra,vza,vxa,vya,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c
c  Update propagated vectors.
c
            do 5322 k = 1,ngr
              elab(i,j,k) = elaba(k)
              mx(i,j,k) = mxa(k) 
              my(i,j,k) = mya(k)
              mz(i,j,k) = mza(k)
              r(i,j,k) = ra(k)
              e(i,j,k) = ea(k)
              p(i,j,k) = pa(k)
              n(i,j,k) = na(k)
              vx(i,j,k) = vxa(k)
              vy(i,j,k) = vya(k)
              vz(i,j,k) = vza(k)
              gamma(i,j,k) = gammaa(k)
 5322       continue
 522      continue
 511    continue
c
c  Finally, propagate in x-direction.
c
        do 533 k = 1,ngr
          do 523 j = 1,ngr
            do 513 i = 1,ngr
              elaba(i) = elab(i,j,k)
              mxa(i) = mx(i,j,k) 
              mya(i) = my(i,j,k) 
              mza(i) = mz(i,j,k) 
              ra(i) = r(i,j,k) 
              ea(i) = e(i,j,k)
              pa(i) = p(i,j,k)
              na(i) = n(i,j,k)
              vxa(i) = vx(i,j,k)
              vya(i) = vy(i,j,k)
              vza(i) = vz(i,j,k)
              gammaa(i) = gamma(i,j,k)
 513        continue  
            if (step.eq.2) then
              dt = 0.5d0*dt
              do 5131 i = 1,ngr
                elabh(i) = elaba(i)
                mxh(i) = mxa(i)
                myh(i) = mya(i)
                mzh(i) = mza(i)
                rh(i) = ra(i)
 5131         continue
              call prop(elabh,mxh,myh,mzh,rh,vxa,vya,vza,pa)
              call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
              dt = 2d0*dt
            end if
            call prop(elaba,mxa,mya,mza,ra,vxa,vya,vza,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c
c  Update propagated vectors and calculate total energy, momentum,
c  baryon number, and entropy.
c
            do 5132 i = 1,ngr
              elab(i,j,k) = elaba(i)
              mx(i,j,k) = mxa(i) 
              my(i,j,k) = mya(i)
              mz(i,j,k) = mza(i)
              r(i,j,k) = ra(i)
              e(i,j,k) = ea(i)
              p(i,j,k) = pa(i)
              n(i,j,k) = na(i)
              vx(i,j,k) = vxa(i)
              vy(i,j,k) = vya(i)
              vz(i,j,k) = vza(i)
              gamma(i,j,k) = gammaa(i)
              if (e(i,j,k).gt.cut) then
                s = entro(e(i,j,k),n(i,j,k))
              end if
              lam(i,j,k) = lambda(e(i,j,k),n(i,j,k))
 5132       continue
 523      continue
 533   continue

c
c  (vi) yxz:
c
      else 
c
c  For each xz-plane, put all field variables into
c  auxiliary vectors and then call prop.
c
        do 631 k = 1,ngr
          do 611 i = 1,ngr
            do 621 j = 1,ngr
              elaba(j) = elab(i,j,k)
              mxa(j) = mx(i,j,k) 
              mya(j) = my(i,j,k)
              mza(j) = mz(i,j,k) 
              ra(j) = r(i,j,k)
              ea(j) = e(i,j,k)
              pa(j) = p(i,j,k)
              na(j) = n(i,j,k)
              vxa(j) = vx(i,j,k)
              vya(j) = vy(i,j,k)
              vza(j) = vz(i,j,k) 
              gammaa(j) = gamma(i,j,k)
 621        continue
            if (step.eq.2) then
              dt = 0.5d0*dt
              do 6211 j = 1,ngr
                elabh(j) = elaba(j)
                mxh(j) = mxa(j)
                myh(j) = mya(j)
                mzh(j) = mza(j)
                rh(j) = ra(j)
 6211         continue
              call prop(elabh,myh,mxh,mzh,rh,vya,vxa,vza,pa)
              call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
              dt = 2d0*dt
            end if
            call prop(elaba,mya,mxa,mza,ra,vya,vxa,vza,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c
c  Update propagated vectors.
c
            do 6212 j = 1,ngr
              elab(i,j,k) = elaba(j)
              mx(i,j,k) = mxa(j)
              my(i,j,k) = mya(j)
              mz(i,j,k) = mza(j)
              r(i,j,k) = ra(j)
              e(i,j,k) = ea(j)
              p(i,j,k) = pa(j)
              n(i,j,k) = na(j)
              vx(i,j,k) = vxa(j)
              vy(i,j,k) = vya(j)
              vz(i,j,k) = vza(j)
              gamma(i,j,k) = gammaa(j)
 6212       continue
 611      continue
c
c  Afterwards, propagate in x-direction.
c
          do 622 j = 1,ngr
            do 612 i = 1,ngr
              elaba(i) = elab(i,j,k)
              mxa(i) = mx(i,j,k) 
              mya(i) = my(i,j,k) 
              mza(i) = mz(i,j,k) 
              ra(i) = r(i,j,k) 
              ea(i) = e(i,j,k)
              pa(i) = p(i,j,k)
              na(i) = n(i,j,k)
              vxa(i) = vx(i,j,k)
              vya(i) = vy(i,j,k)
              vza(i) = vz(i,j,k)
              gammaa(i) = gamma(i,j,k)
 612        continue  
            if (step.eq.2) then
              dt = 0.5d0*dt
              do 6121 i = 1,ngr
                elabh(i) = elaba(i)
                mxh(i) = mxa(i)
                myh(i) = mya(i)
                mzh(i) = mza(i)
                rh(i) = ra(i)
 6121         continue
              call prop(elabh,mxh,myh,mzh,rh,vxa,vya,vza,pa)
              call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1                    vxa,vya,vza,gammaa)
              dt = 2d0*dt
            end if
            call prop(elaba,mxa,mya,mza,ra,vxa,vya,vza,pa)
            call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1                  vxa,vya,vza,gammaa)
c
c  Update propagated vectors.
c
            do 6122 i = 1,ngr
              elab(i,j,k) = elaba(i)
              mx(i,j,k) = mxa(i) 
              my(i,j,k) = mya(i)
              mz(i,j,k) = mza(i)
              r(i,j,k) = ra(i)
              e(i,j,k) = ea(i)
              p(i,j,k) = pa(i)
              n(i,j,k) = na(i)
              vx(i,j,k) = vxa(i)
              vy(i,j,k) = vya(i)
              vz(i,j,k) = vza(i)
              gamma(i,j,k) = gammaa(i)
 6122       continue
 622      continue
 631    continue
c
c  Finally, propagate in z-direction.
c
        do 623 j = 1,ngr
           do 613 i = 1,ngr
              do 633 k = 1,ngr
                 elaba(k) = elab(i,j,k)
                 mxa(k) = mx(i,j,k) 
                 mya(k) = my(i,j,k) 
                 mza(k) = mz(i,j,k) 
                 ra(k) = r(i,j,k) 
                 ea(k) = e(i,j,k)
                 pa(k) = p(i,j,k)
                 na(k) = n(i,j,k)
                 vxa(k) = vx(i,j,k)
                 vya(k) = vy(i,j,k)
                 vza(k) = vz(i,j,k)
              gammaa(k) = gamma(i,j,k)
 633       continue  
           if (step.eq.2) then
              dt = 0.5d0*dt
              do 6331 k = 1,ngr
                 elabh(k) = elaba(k)
                 mxh(k) = mxa(k)
                 myh(k) = mya(k)
                 mzh(k) = mza(k)
                 rh(k) = ra(k)
 6331        continue
             call prop(elabh,mzh,mxh,myh,rh,vza,vxa,vya,pa)
             call untang(elabh,mxh,myh,mzh,rh,ea,pa,na,
     1            vxa,vya,vza,gammaa)
             dt = 2d0*dt
          end if
          call prop(elaba,mza,mxa,mya,ra,vza,vxa,vya,pa)
          call untang(elaba,mxa,mya,mza,ra,ea,pa,na,
     1         vxa,vya,vza,gammaa)
c     
c     Update propagated vectors and calculate total energy, momentum,
c     baryon number, and entropy.
c     
          do 6332 k = 1,ngr
             elab(i,j,k) = elaba(k)
             mx(i,j,k) = mxa(k) 
             my(i,j,k) = mya(k)
             mz(i,j,k) = mza(k)
             r(i,j,k) = ra(k)
             e(i,j,k) = ea(k)
             p(i,j,k) = pa(k)
             n(i,j,k) = na(k)
             vx(i,j,k) = vxa(k)
             vy(i,j,k) = vya(k)
             vz(i,j,k) = vza(k)
             gamma(i,j,k) = gammaa(k)
            if (e(i,j,k).gt.cut) then
                s = entro(e(i,j,k),n(i,j,k))
             end if
             lam(i,j,k) = lambda(e(i,j,k),n(i,j,k))
 6332     continue
 613   continue
 623  continue
      end if
      return
 7777 format(7(1x,f13.5))
      end
