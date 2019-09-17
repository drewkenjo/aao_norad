C========================================================================
      
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
c MESONMASS is the mass of the pi0 or eta

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
     * EPS                  *S_TT  * COS(2*PHI_G)     + 
     * SQRT(2.*EPS*(1.+EPS))*S_LT  * COS(PHI_G)       +
     * HELI*SQRT(2.*EPS*(1.-EPS))*S_LTP * SIN(2*PHI_G)     
     * ) 
      IF(DVMPX.lt.0.) DVMPX=0.

      RETURN
      END
c-------------------------------------------------------------------------

      SUBROUTINE XSINIT(MESONMASS)
C
      IMPLICIT NONE
      REAL MESONMASS

      REAL PAR
      COMMON /PAWPAR/ PAR(12)

      INTEGER I
C
C------Parametrization based on the pi0/eta CLAS12 data
C
C      H_T and E_T convolutions depends on Q2, xB and t.
C      See the HT and ET functions for more details
C      Valery Kubarovsky
C
C------------------------------------  23 April 2019
       REAL AM(2)
       DATA AM/0.134976,0.547300/
       INTEGER K
       COMMON/PIETA/K
       REAL P(12,2)
       DATA P/
     *   6.2624,  0.8032, -0.1152,  1.7431,  0.0,
     *  92.9727,  3.4071, -1.9099, -0.1199,  0.0, 16.8603,  2.1722,
     *   6.9439,  1.7523, -1.2875,  0.6822,  0.0,
     *  17.0423,  1.1264,  0.0491,  1.6594,  0.0, 21.6379,  3.9873/

C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C---------------------------------------------------------------------------------------------------
C       KUMAC file: fitxs.kumac
C
C                     pi0:  With systematic errors
C       6.2624  0.8032 -0.1152  1.7431  0.0    92.9727  3.4071 -1.9099 -0.1199  0.0 16.8603  2.1722
C  +/-  0.1308  0.0097  0.0167  0.0582  0.0     2.4662  0.0205  0.0234  0.0999  0.0  1.4065  0.1689

C                     pi0:  statistical errors ONLY errors
C       6.2657  0.7996 -0.1395  1.6899  0.0    94.1568  3.4664 -1.9103 -0.0380  0.0 17.2384  2.1922
C  +/-  0.1038  0.0079  0.0141  0.0487  0.0     1.9730  0.0174  0.0207  0.0821  0.0  1.1245  0.1342

C                     eta:  With systematic errors
C       6.9439  1.7523 -1.2875  0.6822  0.0    17.0423  1.1264  0.0491  1.6594  0.0 21.6379  3.9875
C  +/-  0.2615  0.0239  0.0229  0.1260  0.0     1.4616  0.0288  0.0663  0.3117  0.0  5.2832  0.8186

C                     eta:  statistical errors ONLY errors
C       7.2312  1.7912 -1.3468  0.5201  0.0    15.3211  1.0784  0.1371  1.9489  0.0 11.8179  2.2773
C  +/-  0.2270  0.0200  0.0191  0.1070  0.0     1.1253  0.0242  0.0541  0.2607  0.0  2.4073  0.4269

C---------------------------------------------------------------------------------------------------

      

      IF(MESONMASS.GT.0.140) THEN
         K=2
      ELSE
         K=1
      ENDIF

      DO I=1,12
         PAR(I)=P(I,K)
      ENDDO

      RETURN
      END


C=============================================================================
      REAL FUNCTION XRXSB(del2,x,Q2, Phi_g,E,heli,MESONMASS)
      IMPLICIT NONE
      REAL       del2,x,Q2,Phi_g,E,heli,MESONMASS
      REAL       DVMPX,fluxw, pi
      PARAMETER (PI=3.1415926536)

      CALL XSINIT(MESONMASS)

      XRXSB=DVMPX(del2,x,Q2,Phi_g,E,heli,MESONMASS)/FLUXW(x,Q2,E)

      RETURN
      END

C====================================================================================c
      REAL FUNCTION HT(del2,X,Q2)
      IMPLICIT NONE
      REAL          del2,x,Q2
      REAL          SLOPE,T
      REAL P1,P2,P3,P4,P5
      REAL P
      REAL W2
      COMMON /PAWPAR/ P(12)

      T=-DEL2
      P1=P(1)
      P2=P(2)
      P3=P(3)
      P4=P(4)
      P5=P(5)
      SLOPE=P2+P3*(ALOG(X)-ALOG(0.15))
      W2=Q2*(1./X-1.)+0.93827**2
      HT=P1 * EXP(-SLOPE*T) * Q2**(P4/2.)
      RETURN
      END      

C====================================================================================c
      REAL FUNCTION ET(del2,X,Q2)
      IMPLICIT NONE
      REAL          del2,x,Q2
      REAL          SLOPE,T
      REAL P1,P2,P3,P4,P5
      REAL P
      REAL W2
      COMMON /PAWPAR/ P(12)

      T=-DEL2
      P1=P(6)
      P2=P(7)
      P3=P(8)
      P4=P(9)
      P5=P(10)
      SLOPE=P2+P3*(ALOG(X)-ALOG(0.15))
      W2=Q2*(1./X-1.)+0.93827**2
      ET=P1*EXP(-SLOPE*T) * Q2**(P4/2.)
      RETURN
      END
C====================================================================================c
      REAL FUNCTION HTEBAR(del2,X,Q2)
      IMPLICIT NONE
      REAL          del2,x,Q2
      REAL T
      REAL P11,P12
      REAL P
      COMMON /PAWPAR/ P(12)

      T=-DEL2
      P11=P(11)
      P12=P(12)
      HTEBAR=P11*EXP(-P12*T)

      RETURN
      END
C====================================================================================c

      REAL FUNCTION XSIGMA_U(del2,x,Q2,E)
      IMPLICIT NONE
      REAL          del2,x,Q2,E
      EXTERNAL      XSIGMA_T,XSIGMA_L,EPSILON
      REAL          XSIGMA_T,XSIGMA_L,EPSILON
      LOGICAL       XCHECK_KINE

      IF(XCHECK_KINE(del2,x,Q2,E)) THEN
        XSIGMA_U=XSIGMA_T(DEL2,X,Q2,E)+
     *  EPSILON(x,Q2,E)*XSIGMA_L(DEL2,X,Q2,E)
      ELSE
         XSIGMA_U=0.0
      ENDIF
      RETURN
      END

C====================================================================================c

      REAL FUNCTION XSIGMA_T(del2,x,Q2,E)
      IMPLICIT NONE
      REAL          del2,x,Q2,E,ALPHA,KSI
      EXTERNAL      tminq,PHASE,HT,ET
      REAL          tminq,T0,SLOPE,T,PI,PHASE
      REAL          HT,ET,MP,HC2
      LOGICAL       XCHECK_KINE
      DATA ALPHA/0.00729927/
      DATA PI/3.14159265/
      DATA MP/0.938272/
      DATA HC2/389379.36/
      
      IF(XCHECK_KINE(del2,x,Q2,E)) THEN
        T0=tminq(Q2,X)
        T=-DEL2
        KSI=X/(2-X)*(1+MP**2/Q2)
        XSIGMA_T=4.*PI*ALPHA/2./PHASE(Q2,X)/Q2**2*(
     *          (1-KSI**2)*    HT(del2,X,Q2)**2
     *       +  (T-T0)/8/MP**2*ET(del2,X,Q2)**2
     *   )
      ELSE
         XSIGMA_T=0.0
      ENDIF
      XSIGMA_T=XSIGMA_T*HC2

      RETURN
      END
      

C====================================================================================c

       REAL FUNCTION XSIGMA_TT(del2,x,Q2,E)
       IMPLICIT NONE
       REAL          del2,x,Q2,E
       EXTERNAL      tminq
       REAL          tminq,T0,SLOPE,T,HC2
       REAL ALPHA,PI,MP,PHASE,ET
       LOGICAL XCHECK_KINE
       DATA ALPHA/0.00729927/
       DATA PI/3.14159265/
       DATA MP/0.938272/
       DATA HC2/389379.36/

       IF(XCHECK_KINE(del2,x,Q2,E)) THEN
         T0=tminq(Q2,X)
         T=-DEL2
        XSIGMA_TT=-4.*PI*ALPHA/2./PHASE(Q2,X)/Q2**2
     *         *  (T-T0)/8/MP**2*ET(del2,X,Q2)**2
      ELSE
         XSIGMA_TT=0.0
      ENDIF

      XSIGMA_TT=XSIGMA_TT*HC2
      
      RETURN
      END

C====================================================================================c

      REAL FUNCTION XSIGMA_L(del2,x,Q2,E)
      IMPLICIT NONE
      REAL          del2,x,Q2,E

      XSIGMA_L=0.0

      RETURN
      END
C====================================================================================c

       REAL FUNCTION XSIGMA_LT(del2,x,Q2,E)
       IMPLICIT NONE
       REAL          del2,x,Q2,E
       EXTERNAL      tminq
       REAL          tminq,T0,SLOPE,T,HC2
       REAL          ALPHA,PI,MP,PHASE,KSI
       REAL          HTEBAR
       LOGICAL       XCHECK_KINE
       DATA          ALPHA/0.00729927/
       DATA          PI/3.14159265/
       DATA          MP/0.938272/
       DATA          HC2/389379.36/

       IF(XCHECK_KINE(del2,x,Q2,E)) THEN
       T0=tminq(Q2,X)
       T=-DEL2
       KSI=X/(2-X)*(1+MP**2/Q2)
       XSIGMA_LT=4.*PI*ALPHA/SQRT(2.)/PHASE(Q2,X)/Q2**(3./2.)
     *         * KSI*SQRT(1.0-KSI**2)
     *         * SQRT((T-T0))/2/MP * HTEBAR(del2,X,Q2)**2
      ELSE
         XSIGMA_LT=0.0
      ENDIF

      XSIGMA_LT=XSIGMA_LT*HC2
      RETURN
      END
C====================================================================================c

      REAL FUNCTION XSIGMA_LTP(del2,x,Q2,E)
      IMPLICIT NONE
      REAL          del2,x,Q2,E

      XSIGMA_LTP=0.

      RETURN
      END

C====================================================================================c

      LOGICAL FUNCTION XCHECK_KINE(del2,xb,Q2,E)
      IMPLICIT NONE
      REAL del2,xb,Q2,E
      REAL Mp, mele, pi, mpi0
      REAL fluxw, rxs
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)
           

      REAL nu,W2,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      REAL  xmin1,xmax1
      REAL y, e1, epsilon,tminq

      REAL AM(2)
      DATA AM/0.134976,0.547300/
      INTEGER K
      COMMON/PIETA/K

      MPI0=AM(K)
      XCHECK_KINE=.FALSE.
      xmin1 = Q2/(2.0*Mp*E)
      xmax1 = 1.0
      nu  = Q2/(2D0*Mp*xb)
      W2  = Mp**2 + 2.0*Mp*nu - Q2
C      print *,'W2=',W2
      IF(W2.LT.(Mp+Mpi0)**2)                  RETURN
      W   = sqrt(W2)
      qmod = sqrt(nu**2 + Q2)

      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp**2-Mpi0**2)/(2.*W)
C      print *,'E2cm=',E2cm,Mp
      IF(E2cm.LE.Mp)                          RETURN
      P2cm = SQRT(E2CM**2 - Mp**2)
      del2max = 2.0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2.0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)
C      print *,'xB=',xb,xmin1,xmax1
      IF( xb.le.xmin1 .or. xb.gt.xmax1 )     RETURN !    x  out of range
C      print *,'t=',del2,del2min,del2max,tminq(Q2,xb)
      IF( del2.ge.del2min .or. del2.le.del2max ) return   ! delta out of range

      y=Q2/(2*Mp*xb*E)
      e1=(y*xb*Mp)**2/Q2
      EPSILON=(1.0-y-e1)/(1-y+y**2/2+e1)
C      print *,'Epsilon',epsilon
      IF(EPSILON.LT.0. .OR.EPSILON .GT.1.)   RETURN

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

      REAL AM(2)
      DATA AM/0.134976,0.547300/
      INTEGER K
      COMMON/PIETA/K
      
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

      RETURN
      END

C====================================================================================c
      REAL FUNCTION PHASE(Q2,X)
      IMPLICIT NONE
      REAL Q2,X
      REAL W,W2,W4,Q4,MP,MP2,MP4,PI,LAMBDA
      DATA MP/0.938272/
      DATA PI/3.14159265/
      PHASE=0.0
      
      IF(X.GE.1.) RETURN
      W2=Q2*(1/X-1.)+MP**2
      W4=W2*W2
      MP2=MP*MP
      MP4=MP2*MP2

      LAMBDA=W4+Q4+MP4+2*W2*Q2-2*W2*MP2+2*Q2*MP2
      IF(LAMBDA.LE.0) RETURN

      PHASE=16.*PI*(W2-MP2)*SQRT(LAMBDA)
      RETURN
      END
      
c*****************The end of Valery Model******************************
c
      REAL FUNCTION SIGMA_T_TT(del2,xB,Q2,E)
      SIGMA_T_TT=XSIGMA_T(del2,xB,Q2,E)+XSIGMA_TT(del2,xB,Q2,E)
      RETURN
      END
C----------------------------------------
      SUBROUTINE PRINTCOM
      COMMON /PAWPAR/F(12)
      PRINT 100, F
 100  FORMAT('/PAWPAR/',12F8.4)
      RETURN
      END
      
C====================================================================================c
      REAL FUNCTION ET2(DEL2,X,Q2)
C
C     Coefficient to calculate GFF^2: |ET|^2 = ET2(del2,x,Q2)*sigma_TT
C
      IMPLICIT NONE
      REAL          del2,x,Q2,E,ALPHA,KSI,MP,HC2
      EXTERNAL      tminq,PHASE
      REAL          tminq,T0,SLOPE,T,PI,PHASE
      DATA ALPHA/0.00729927/
      DATA PI/3.14159265/
      DATA MP/0.938272/
      DATA HC2/389379.36/

      T0=tminq(Q2,X)
      T=-DEL2
      ET2=-1.0/(HC2*
     * 4.*PI*ALPHA/2./PHASE(Q2,X)/Q2**2*(T-T0)/8/MP**2
     *     )
C      print *, 'ET2= ',del2,x,q2,T,T0,ET2

      RETURN
      END
C====================================================================================c
      REAL FUNCTION HT2(DEL2,X,Q2)
C
C     Coefficient to calculate GFF^2: |HT|^2 = HT2(del2,x,Q2)*(sigma_T + sigma_TT)
C
      IMPLICIT NONE
      REAL          del2,x,Q2,E,ALPHA,KSI,MP,HC2
      EXTERNAL      tminq,PHASE
      REAL          tminq,T0,SLOPE,T,PI,PHASE
      DATA ALPHA/0.00729927/
      DATA PI/3.14159265/
      DATA MP/0.938272/
      DATA HC2/389379.36/

      T0=tminq(Q2,X)
      T=-DEL2
      KSI=X/(2-X)*(1+MP**2/Q2)
      HT2=1.0/(HC2*
     * 4.*PI*ALPHA/2./PHASE(Q2,X)/Q2**2*(1-ksi**2)
     * )

C      print *, 'HT2= ',del2,x,q2,HT2

      RETURN
      END
C====================================================================================c
      REAL FUNCTION HTE2(DEL2,X,Q2)
C
C     Coefficient to calculate Re(HT*Ebar): HTETBAR = HTE(del2,x,Q2)*sigma_LT
C
      IMPLICIT NONE
      REAL          del2,x,Q2,E,ALPHA,KSI,MP,HC2
      EXTERNAL      tminq,PHASE
      REAL          tminq,T0,SLOPE,T,PI,PHASE
      DATA ALPHA/0.00729927/
      DATA PI/3.14159265/
      DATA MP/0.938272/
      DATA HC2/389379.36/

      T0=tminq(Q2,X)
      T=-DEL2
      KSI=X/(2-X)*(1+MP**2/Q2)
      HTE2=1.0/(HC2*
     *  4.*PI*ALPHA/SQRT(2.)/PHASE(Q2,X)/Q2**(3./2.)
     *         * KSI*SQRT(1-KSI**2)
     *         * SQRT((T-T0))/2/MP)
C      print *, 'HTE2= ',del2,x,q2,T,T0,HTE2

      RETURN
      END


