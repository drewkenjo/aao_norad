c*******************Valery Model****************************************
      subroutine DVMPW(COST,W,Q2, Phi_g,E,heli,MESONMASS,
     1                 sigma0,S_T,S_L,S_TT,S_LT,S_LTP)
C
c  dsigma/dQ2 dphi dCosTheta* dW for ep-->ep pi0
C
C  exc pi0 x-section
c
cinput:
c COST = Cos(Theta*), Theta* is the angle (pi0-gamma* or  p-p')in the CM system
c W is W
c xb,Q2 x and Q^2 (GeV^2)
c Phi_g angle in the photon frame (radians)
c E energy of the electron in GeV
c heli electron helicity -1.0 or +1.0
c
cdel2=t (negative GeV^2)           NEGATIVE !!!!
C
C MESONMASS is the mass of the pi0 or eta.
C The actual masses that will be used for the calculations are in pi0eta.par file


      IMPLICIT NONE
      REAL    COST,W,del2,xb,Q2, Phi_g,E,heli, mesonmass
      REAL    Mp, mele, pi
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)

      REAL     sigma0,S_T, S_L, S_LT, S_TT,S_LTP
      REAL     SIGMA_T, SIGMA_L, SIGMA_LT, SIGMA_TT,SIGMA_LTP
      REAL     EPS, EPSILON, FLUXW, SIGMA_TOT
      EXTERNAL EPSILON,SIGMA_T,SIGMA_L, SIGMA_LT, SIGMA_TT
      EXTERNAL SIGMA_LTP
      REAL     JACG,JACR
      LOGICAL  CHECK_KINE
c      REAL DVMPX
      INCLUDE 'pi0eta.par'
c
      
      CALL XSINIT(MESONMASS)

c     
      S_T = 0.0
      S_L = 0.0
      S_LT = 0.0
      S_TT = 0.0
      S_LTP = 0.0
      IF(.NOT.CHECK_KINE(COST,W,Q2,E,del2,xb,JACG,JACR)) RETURN 
       EPS=EPSILON(XB,q2,e)
       S_T  =  SIGMA_T   (COST,W,Q2,E)
       S_L  =  SIGMA_L   (COST,W,Q2,E)
       S_LT =  2.0*SIGMA_LT  (COST,W,Q2,E)
       S_TT =  SIGMA_TT  (COST,W,Q2,E)
       S_LTP = 2.0*SIGMA_LTP (COST,W,Q2,E)

c
       sigma0 =FLUXW(xb,Q2,E)/(2.*PI)*JACG*(
     *                       S_T   + EPS*S_L          + 
     * EPS                  *S_TT  * COS(2*PHI_G)       + 
     * SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G)     +
     * HELI*SQRT(2.*EPS*(1.-EPS))*S_LTP * SIN(2*PHI_G)     
     * )
      
       IF (sigma0.LT.0.0) then
         sigma0=0.0
         S_T = 0.0 
         S_L = 0.0 
         S_LT = 0.0 
         S_TT = 0.0 
         S_LTP = 0.0 
       ENDIF
c       PRINT *,DVMPW,COST,W,Q2,del2,xb,phi_g,JACG,JACR
c      PRINT *,DVMPX(del2,xb,Q2, Phi_g,E,heli,MESONMASS)*JACG*JACR
c      PRINT *,DVMPX(del2,xb,Q2, Phi_g,E,heli,MESONMASS)
c      PRINT *,del2,xb,Q2, Phi_g,E,heli,MESONMASS

c       PRINT *,DVMPW,COST,W,Q2,phi_g,eps
c       print *,'DVMP',sigma0, S_T,S_L,S_TT,S_LT
c       print *,'S_TT  ',  S_TT, EPS     *S_TT  * COS(2*PHI_G) 
c       print *, 'S_LT ',S_LT,SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G) 

      RETURN
      END

c-------------------------------------------------------------------------

      SUBROUTINE XSINIT(MESONMASS)
C
      IMPLICIT NONE
      REAL MESONMASS
      INCLUDE 'pi0eta.par'

      IF(MESONMASS.GT.0.140) THEN
           K=2
      ELSE
           K=1
      ENDIF

      RETURN
      END


C====================================================================================c

      REAL FUNCTION SIGMA_T(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T,JACG,JACR
      LOGICAL       CHECK_KINE
      EXTERNAL      XSIGMA_T
      REAL          XSIGMA_T
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_T=JACR*XSIGMA_T(del2,x,Q2,E)
      ELSE
        SIGMA_T=0.0
      ENDIF
      RETURN
      END
C=============================================================      
      REAL FUNCTION SIGMA_L(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T,JACG,JACR
      LOGICAL CHECK_KINE
      EXTERNAL      XSIGMA_L
      REAL          XSIGMA_L
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_L=JACR*XSIGMA_L(del2,x,Q2,E)
      ELSE
        SIGMA_L=0.0
      ENDIF
      RETURN
      END
      
C=============================================================      
      REAL FUNCTION SIGMA_TT(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T,JACG,JACR
      LOGICAL CHECK_KINE
      EXTERNAL      XSIGMA_TT
      REAL          XSIGMA_TT
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_TT=JACR*XSIGMA_TT(del2,x,Q2,E)
      ELSE
        SIGMA_TT=0.0
      ENDIF
      RETURN
      END
      
C=============================================================      
      REAL FUNCTION SIGMA_LT(COST,W,Q2,E)
      IMPLICIT NONE
      REAL          COST,W,del2,x,Q2,E
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T,JACG,JACR
      LOGICAL CHECK_KINE
      EXTERNAL      XSIGMA_LT
      REAL          XSIGMA_LT
      
      IF(CHECK_KINE(COST,W,Q2,E,del2,x,JACG,JACR)) THEN
        SIGMA_LT=JACR*XSIGMA_LT(del2,x,Q2,E)
      ELSE
        SIGMA_LT=0.0
      ENDIF
      RETURN
      END
      
C====================================================================================c

       REAL FUNCTION SIGMA_LTP(COST,W,Q2,E)
       IMPLICIT NONE
       REAL COST,W,del2,x,Q2,E
       REAL EPS, EPSILON,A2_PHI,FLUXW,SIGMA_TOT
       EXTERNAL EPSILON,FLUXW

       SIGMA_LTP=0.

       RETURN
       END

C====================================================================================c

      LOGICAL FUNCTION CHECK_KINE(COST,W,Q2,E,del2,xb,JACG,JACR)
      IMPLICIT NONE
      REAL COST,W,del2,xb,Q2,E,JACG,JACR
      real Mp, mele, pi,mpi0
      real fluxw, rxs
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)
           

      real nu,W2,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      real  xmin1,xmax1
      real y, e1, epsilon
      INCLUDE 'pi0eta.par'

      MPI0=AM(K)

      CHECK_KINE=.FALSE.
      W2=W*W
      XB=Q2/(W2+Q2-Mp*Mp)
      nu  = (W2+Q2-Mp*Mp)/(2*Mp)
      qmod = sqrt(nu**2 + Q2)
      IF(W         .LT. Mp+Mpi0)              RETURN
      IF(W         .GT. -Q2+2*Mp*E+Mp*Mp)     RETURN
      IF(Q2/(4*E*(E-NU)) .GE. 1.0)            RETURN
      IF(XB        .LT. Q2/(2*Mp*E))          RETURN
      IF(ABS(COST) .GT. 1.0)                  RETURN
      xmin1 = Q2/(2.0*Mp*E)
      xmax1 = 1.0

      E1cm = Mp*(Mp + nu)/W
      E2cm = (W2 + Mp**2-Mpi0**2)/(2D0*W)
      IF(E1cm.LE.Mp .OR. E2cm.LE. Mp) RETURN
      P1cm = Mp*qmod/W
      P2cm = SQRT(E2CM**2 - Mp**2)
      del2max = 2.0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2.0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)
      del2=del2min-2.*P1cm*P2cm*(1-COST)
      IF( xb.le.xmin1 .or. xb.gt.xmax1 )         return   !    x  out of range
      IF( del2.ge.del2min .or. del2.le.del2max ) return   ! delta out of range

      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)**2/Q2
      EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)

      IF(EPSILON.LT.0. .OR.EPSILON .GT.1.)         RETURN

C      jacobian
C      JAC=4*P1cm*P2CM*W*XB/(W2+Q2-Mp*Mp)
      JACG=2*W*XB/(W2+Q2-Mp*Mp)
      JACR=2*P1cm*P2CM

      CHECK_KINE=.TRUE.

      RETURN
      END

C=======================================================================================C


      REAL FUNCTION DVMPX(del2,xb,Q2, Phi_g,E,heli,MESONMASS)
C
c  dsigma/dQ2 dX dt dphi for ep-->ep pi0
C
C exc pi0 x-section
c
cinput:
cdel2=t (negative GeV^2)           NEGATIVE !!!!
cxb,Q2 x and Q^2 (GeV^2)
c Phi_g angle in the photon frame (radians)
c E energy of the electron in GeV
c heli electron helicity -1.0 or +1.0
c MESONMASS is the mass of th epi0 or eta

      IMPLICIT NONE
      REAL del2,xb,Q2, Phi_g,E,heli,MESONMASS
      real Mp, mele, pi
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)

      REAL     S_T, S_L, S_LT, S_TT,S_LTP
      REAL     XSIGMA_T, XSIGMA_L, XSIGMA_LT, XSIGMA_TT,XSIGMA_LTP
      REAL     EPS, EPSILON, FLUXW, SIGMA_TOT
      EXTERNAL EPSILON,XSIGMA_T,XSIGMA_L, XSIGMA_LT, XSIGMA_TT
      EXTERNAL XSIGMA_LTP
      REAL     FL
      LOGICAL  XCHECK_KINE
      
      INCLUDE 'pi0eta.par'
c
      CALL XSINIT(MESONMASS)

      DVMPX = 0.0
      
      IF(.NOT.XCHECK_KINE(del2,xb,Q2,E)) RETURN
       EPS  =  EPSILON(XB,q2,e)
       S_T  =  XSIGMA_T   (del2,xb,Q2,E)
       S_L  =  XSIGMA_L   (del2,xb,Q2,E)
       S_LT =  XSIGMA_LT  (del2,xb,Q2,E)
       S_TT =  XSIGMA_TT  (del2,xb,Q2,E)
       S_LTP = XSIGMA_LTP (del2,xb,Q2,E)

       DVMPX =FLUXW(xb,Q2,E)/(2.*PI)*(
     *                       S_T   + EPS*S_L          + 
     * EPS                  *S_TT  * COS(2*PHI_G)       + 
     * SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G)     +
     * HELI*SQRT(2.*EPS*(1.-EPS))*S_LTP * SIN(2*PHI_G)     
     * ) 
      if(DVMPX.lt.0.) DVMPX=0.

      RETURN
      END
C=============================================================================
      REAL FUNCTION XRXSB(del2,x,Q2, Phi_g,E,heli,MESONMASS)
      IMPLICIT NONE
      REAL       del2,x,Q2,Phi_g,E,heli,MESONMASS
      REAL       DVMPX,fluxw, pi
      PARAMETER (PI=3.1415926536)

      INCLUDE 'pi0eta.par'

      CALL XSINIT(MESONMASS)

      XRXSB=DVMPX(del2,x,Q2,Phi_g,E,heli,MESONMASS)/FLUXW(x,Q2,E)
     *      *(2.*PI)

      RETURN
      END

C====================================================================================c

      REAL FUNCTION XSIGMA_T(del2,x,Q2,E)
      IMPLICIT NONE
      REAL          del2,x,Q2,E
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T
      LOGICAL       XCHECK_KINE
      INCLUDE 'pi0eta.par'
      
      IF(XCHECK_KINE(del2,x,Q2,E)) THEN
        T0=tminq(Q2,X)
        SLOPE = 2.*1.1*ALOG(X)
        T=-DEL2
        XSIGMA_T=(P(1,k)+P(2,k)*SQRT(T-T0))*EXP(SLOPE*T*P(3,k))
     *           *X**P(11,k)*(1-X)**P(12,k)/(Q2+P(13,k))**P(14,k)
      ELSE
         XSIGMA_T=0.0
      ENDIF
      RETURN
      END
      
C====================================================================================c

      REAL FUNCTION XSIGMA_L(del2,x,Q2,E)
      IMPLICIT NONE
      REAL          del2,x,Q2,E
      EXTERNAL      tminq
      REAL          tminq,T0,SLOPE,T
      LOGICAL XCHECK_KINE
      INCLUDE 'pi0eta.par'
      
      IF(XCHECK_KINE(del2,x,Q2,E)) THEN
        T0=tminq(Q2,X)
        SLOPE = 2.*1.1*ALOG(X)
        T=-DEL2
        XSIGMA_L=Q2*(P(4,k)+P(5,k)*(T-T0))*EXP(SLOPE*T*P(6,k))
     *           *X**P(11,k)*(1-X)**P(12,k)/(Q2+P(13,k))**P(14,k)
      ELSE
         XSIGMA_L=0.0
      ENDIF
       RETURN
       END

C====================================================================================c

       REAL FUNCTION XSIGMA_TT(del2,x,Q2,E)
       IMPLICIT NONE
       REAL          del2,x,Q2,E
       EXTERNAL      tminq
       REAL          tminq,T0,SLOPE,T
       LOGICAL XCHECK_KINE
       INCLUDE 'pi0eta.par'

       IF(XCHECK_KINE(del2,x,Q2,E)) THEN
         T0=tminq(Q2,X)
         SLOPE = 2.*1.1*ALOG(X)
         T=-DEL2
         XSIGMA_TT=P(7,k)*(T-T0)*EXP(SLOPE*T*P(8,k))
     *           *X**P(11,k)*(1-X)**P(12,k)/(Q2+P(13,k))**P(14,k)
       ELSE
         XSIGMA_TT=0.0
       ENDIF
       RETURN
       END

C====================================================================================c

       REAL FUNCTION XSIGMA_LT(del2,x,Q2,E)
       IMPLICIT NONE
       REAL          del2,x,Q2,E
       EXTERNAL      tminq
       REAL          tminq,T0,SLOPE,T
       LOGICAL XCHECK_KINE
       INCLUDE 'pi0eta.par'

       IF(XCHECK_KINE(del2,x,Q2,E)) THEN
         T0=tminq(Q2,X)
         SLOPE = 2.*1.1*ALOG(X)
         T=-DEL2
         XSIGMA_LT=P(9,k)*SQRT(T-T0)*EXP(SLOPE*T*P(10,k))
     *           *X**P(11,k)*(1-X)**P(12,k)/(Q2+P(13,k))**P(14,k)
       ELSE
         XSIGMA_LT=0.0
       ENDIF
       RETURN
       END


C====================================================================================c

       REAL FUNCTION XSIGMA_LTP(del2,x,Q2,E)
       IMPLICIT NONE
       REAL del2,x,Q2,E
       REAL EPS, EPSILON,A2_PHI,FLUXW,SIGMA_TOT
       EXTERNAL EPSILON,FLUXW
       INCLUDE 'pi0eta.par'

       XSIGMA_LTP=0.

       RETURN
       END

C====================================================================================c

      LOGICAL FUNCTION XCHECK_KINE(del2,xb,Q2,E)
      IMPLICIT NONE
      REAL del2,xb,Q2,E
      real Mp, mele, pi,mpi0
      real fluxw, rxs
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)
           

      real nu,W2,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      real  xmin1,xmax1
      real y, e1, epsilon

       INCLUDE 'pi0eta.par'

      MPI0=AM(K)
      XCHECK_KINE=.FALSE.
      xmin1 = Q2/(2.0*Mp*E)
      xmax1 = 1.0
      nu  = Q2/(2D0*Mp*xb)
      W2  = Mp**2 + 2.0*Mp*nu - Q2
      IF(W2.LT.(Mp+Mpi0)**2)                               RETURN
      W   = sqrt(W2)
      qmod = sqrt(nu**2 + Q2)

      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp**2-Mpi0**2)/(2.*W)
      IF(E2cm.LE.Mp)                             RETURN
      P2cm = SQRT(E2CM**2 - Mp**2)
      del2max = 2.0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2.0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)

      IF( xb.le.xmin1 .or. xb.gt.xmax1 )         return   !    x  out of range
      IF( del2.ge.del2min .or. del2.le.del2max ) return   ! delta out of range

      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)**2/Q2
      EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)
  
      IF(EPSILON.LT.0. .OR.EPSILON .GT.1.)         RETURN

      XCHECK_KINE=.TRUE.

      RETURN
      END

C====================================================================================c

      real function fluxw(x,Q2,E)
      implicit none
      real x,Q2,E,y,eps,e1
      real alpha,Mp,PI
      parameter (alpha=1.0/137.036,Mp=0.93827231,PI=3.14151926)

c
      y=Q2/(2*Mp*x*E)
      e1=(y*x*Mp)**2/Q2
      eps=(1.0-y-e1)/(1-y+y**2/2+e1)
      fluxw=alpha/(2*PI)*y*y/(1-eps)*(1-x)/x/Q2
      return
      end

C====================================================================================c

      real function EPSILON(x,Q2,E)
      implicit none
      real x,Q2,E,y,eps,e1
      real alpha,Mp,PI
      parameter (alpha=1.0/137.036,Mp=0.93827231,PI=3.14151926)
c
      y=Q2/(2*Mp*x*E)
      e1=(y*x*Mp)**2/Q2
      EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)
      return
      end


C======================================================================================c

      REAL FUNCTION tminq(Q2,X)
      IMPLICIT NONE
      REAL Q2,X
      REAL W2,W,S,E1CM,E3CM,P1CM,P3CM,TMAX
      INTEGER I
      real alpha,Mp,PI,Me
      REAL MPI0

       INCLUDE 'pi0eta.par'
c      print *,x,q2,e,fluxw

      MPI0=AM(K)
      tminq=0.

      Mp=0.9382723
      Me=0.0051
      PI=3.14151926

      IF(X.LE.0. .OR. X.GE.1.)        RETURN
      W2 = Q2*(1./X-1.)+Mp**2
      W=SQRT(W2)
      IF(W.LT.Mp+Mpi0)                RETURN

      E1CM=(W2+Q2+Mp**2)/(2*W)
      P1CM=SQRT(E1CM**2-MP**2)

      E3CM=(W2-MPI0**2+Mp**2)/(2*W)
      P3CM=SQRT(E3CM**2-MP**2)
      
      TMINQ=-((Q2+MPI0**2)**2/4./W2-(P1CM-P3CM)**2)
C      TMAXQ=-(Q2+MPI0**2)**2/4./W2-(P1CM+P3CM)**2
      RETURN
      END


c*****************The end of Valery Model******************************
