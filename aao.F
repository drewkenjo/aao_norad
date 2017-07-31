      subroutine aao(epq2,epw,epeps,epcos,epphi,epirea,
     *sigma0,sigu,sigt,sigl,sigi)
     
      implicit none
      
      real pi
      real epq2,epw,epeps,epcos,epphi,sigma0,sigu,sigt,sigl,sigi
      real tandd,fkt,sig0,sigmao
      real xmprot

***************************************************************************
*      AUTHOR:        V. BURKERT AND Z. LI
*      FIRST VERSION: SUMMER, 1991.  RECENT UPDATE: MAY.1993
* AO IS WRITTEN BASED ON THE ORIGINAL PROGRAM A_AND_O FROM V.BURKERT
* THIS PROGRAM SHOULD BE LINKED TO EITHER QKMC OR QKXM FOR THE CALCULATION
* OF QUARK MODEL.  
* QKMC AND QKXM USE DIFFERENT FORM FACTORS:
* QKMC USES THE TREATMENT OF FOSTER ET. AL.
* QKXM USES THE DIPOLE FORM FACTOR 
* QKMC IS THE DEFAULT CHOICE IN THIS PACKAGE.
* AO CAN BE USED TO EXAM: Q2-DEPENDENCE OF HELICITY AMP. AT RES. POSITION (GO1)
*                         OUTPUT:GDH.TOP: A B CA(A_1/2) CB(A_3/2)
*                         GDH SUM RULE (GO2):OUTPUT GDH.TOP
*                         OBERVABLES  (RETURN):OUTPUT TEST.TOP
**************************************************************************
*           AO.FOR
* AO CONTAINS THE FOLLOWING SUBROUTINES/FUNCTIONS:
* SIGMA--CALCULATES OBSERVABLES
* EPRES--CALCULATES BREIT-WIGNER RESONANCE AMPLITUDES
* EXPA --CALCULATES THE HELICITY AMPLITUTES FROM EXPT (V.BURKERT ORIGINAL)
* BORNT--CALCULATES THE BORN TERM CONTRIBUTIONS
* BACK --CALCULATES BORN AND NON-BORN BACK GROUND TERMSS
* QKMA --CALCULATES THE HELICITY AMPLITUTES FROM QUARK MODEL
* RAMPF --CALCULATES THE Q2-DEPENDENCE OF THE HELICITY AMP. AT RES. POSITIONS
* HAMP --CALCULATES THE ENERGY-DEPENDECE OF THE RESONANCE HELICITY AMPLITUDES
* QKMC --CALCULATES THE COUPLING CONSTANTS FROM QUARK MODEL
************************************************************************
**************************************************************************
* UPDTATED: NOV., 1992
*         * WITH NEW OPTIONS TO TURN THE BORN TERM ON AND OFF
*         * BORN TERM ARE MODIFIED WITH A CUT OFF FACTOR AT HIGHER 
*           ENERGIES (WCM>1.3 GEV)
*         * SOME CORRECTIONS HAVE BEEN MADE TO THE EXPA SUBROUTINE
***************************************************************************  

      COMMON /IP/IT,IB,NF,IBORN,CUT,IP11
      
c     COMMON /CA/CS11P,CP11P,CP333,CP331
C* CA HERE USED FOR ADJUSTING THE COUPLING CONSTANT CALCULATED
C* FROM QUARK MODEL 

C* IT:		DETERMINE WHICH MODEL TO USE, 
C* IB=1,0:	TURNS THE NON-BORN BACKGROUND ON AND OFF
C* IBORN=1,0:	TURNS THE BORN  TERM ON AND OFF
C* NF=1,3:	DETERMINE THE FORM OF P11_1440
C* IP11=1,4:	DETERMINES THE USE OF P11_1440 AMPLITUDE IN EXPT

      real cut
      integer it,ib,nf,iborn,ip11,epirea

       real
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710, 
     * WP13_1720, WD13_1520, WD13_1700, 
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600, 
     * WS31_1620, WS31_1900, WP31_1910,   
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710, 
     * LP13_1720, LD13_1520, LD13_1700, 
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600, 
     * LS31_1620, LS31_1900, LP31_1910,  
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710, 
     * PIBP13_1720, PIBD13_1520, PIBD13_1700, 
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600, 
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,  
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975

       COMMON/RCONST/      
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710, 
     * WP13_1720, WD13_1520, WD13_1700, 
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600, 
     * WS31_1620, WS31_1900, WP31_1910,   
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710, 
     * LP13_1720, LD13_1520, LD13_1700, 
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600, 
     * LS31_1620, LS31_1900, LP31_1910,  
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710, 
     * PIBP13_1720, PIBD13_1520, PIBD13_1700, 
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600, 
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,  
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975
C
C      Branching ratios into N-pion channel
C
      DATA PIBS11_1535      /0.40    / 
      DATA PIBS11_1650      /0.60    / 
      DATA PIBP11_1440      /0.60    / 
      DATA PIBP11_1710      /0.15    /
      DATA PIBP13_1720      /0.15    /
      DATA PIBD13_1520      /0.55    / 
      DATA PIBD13_1700      /0.10    / 
      DATA PIBD15_1675      /0.40    /
      DATA PIBF15_1680      /0.60    /
      DATA PIBG17_2190      /0.14    /
      DATA PIBG19_2250      /0.10    /
      DATA PIBH19_2220      /0.18    /
      DATA PIBI111_2600     /0.05    /
 
      DATA PIBS31_1620      /0.30    / 
      DATA PIBS31_1900      /0.10    / 
      DATA PIBP31_1910      /0.20    /
      DATA PIBP33_1232      /0.994   /
      DATA PIBP33_1920      /0.20    / 
      DATA PIBD33_1700      /0.15    /
      DATA PIBD35_1930      /0.10    /
      DATA PIBF35_1905      /0.10    /
      DATA PIBF37_1950      /0.40    /
      DATA PIBH311_2420     /0.10    / 
 
      DATA PIBP33_1600      /0.19    /
      DATA PIBF17_1990      /0.05    /
      DATA PIBF15_2000      /0.12    /
      DATA PIBP11_2100      /0.11    /
      DATA PIBF35_2000      /0.02    /
      
      DATA PIBP13_1870      /0.02   /
      DATA PIBP31_1925      /0.02   /      
      DATA PIBP13_1980      /0.02   /
      DATA PIBF15_1955      /0.02   /      
      DATA PIBP13_1955      /0.02   /
      DATA PIBP33_1975      /0.02   /      

      DATA LS11_1535      /1.532    / 
      DATA LS11_1650      /1.650    / 
      DATA LP11_1440      /1.440    / 
      DATA LP11_1710      /1.710    /
      DATA LP13_1720      /1.740    /
      DATA LD13_1520      /1.520    / 
      DATA LD13_1700      /1.700    / 
      DATA LD15_1675      /1.675    /
      DATA LF15_1680      /1.688    /
      DATA LG17_2190      /2.175    /
      DATA LG19_2250      /2.250    /
      DATA LH19_2220      /2.220    /
      DATA LI111_2600     /2.640    / 
      DATA LS31_1620      /1.620    / 
      DATA LS31_1900      /1.925    / 
      DATA LP31_1910      /1.905    /
      DATA LP33_1232      /1.232    /
      DATA LP33_1920      /1.920    / 
      DATA LD33_1700      /1.670    /
      DATA LD35_1930      /1.930    /
      DATA LF35_1905      /1.905    /
      DATA LF37_1950      /1.950    /
      DATA LH311_2420     /2.420    / 
      DATA LP33_1600      /1.600   /
      DATA LF17_1990      /1.990   /
      DATA LF15_2000      /2.000   /
      DATA LP11_2100      /2.100   /
      DATA LF35_2000      /2.000   /
      DATA LP13_1870      /1.870   /
      DATA LP31_1925      /1.925   /
      DATA LP13_1980      /1.980   /
      DATA LF15_1955      /1.955   /
      DATA LP13_1955      /1.955   /
      DATA LP33_1975      /1.975   /

      DATA WS11_1535      /.150     / 
      DATA WS11_1650      /.150     / 
      DATA WP11_1440      /.300     / 
      DATA WP11_1710      /.110     /
      DATA WP13_1720      /.200     /
      DATA WD13_1520      /.125     / 
      DATA WD13_1700      /.100     /
      DATA WD15_1675      /.155     /
      DATA WF15_1680      /.125     / 
      DATA WG17_2190      /.350     /
      DATA WG19_2250      /.300     /
      DATA WH19_2220      /.400     /
      DATA WI111_2600     /.400     / 
      DATA WS31_1620      /.140     / 
      DATA WS31_1900      /.150     / 
      DATA WP31_1910      /.220     /
      DATA WP33_1232      /.115     / 
      DATA WP33_1920      /.250     / 
      DATA WD33_1700      /.250     / 
      DATA WD35_1930      /.250     /
      DATA WF35_1905      /.300     /
      DATA WF37_1950      /.240     /
      DATA WH311_2420     /.300     /
      DATA WP33_1600      /.230     /
      DATA WF17_1990      /.300     /
      DATA WF15_2000      /.135     /
      DATA WP11_2100      /.230     /
      DATA WF35_2000      /.300     /
                                                                        
      DATA WP13_1870      /.300     /
      DATA WP31_1925      /.300     /
      DATA WP13_1980      /.300     /
      DATA WF15_1955      /.300     /
      DATA WP13_1955      /.300     /
      DATA WP33_1975      /.300     /
      DATA XMPROT     /.938     /     
c      DATA XMPIP      /.1395    /     
c      DATA XMPI0      /.1349    /     
c      DATA XMETA      /.5488    /     

      pi	= 4.*atan(1.)
      nf	= 1
      it	= 1
      ip11	= 1
      ib	= 0
      iborn	= 1
      cut	= 0.4       
      sigma0	= 0.


      sigma0 = SIGMAO(EPW,EPQ2,EPEPS,EPCOS,EPPHI,EPIREA,
     *   0.,0.,0.,0.,sig0,sigu,sigt,sigl,sigi,fkt)
     
      end
C**************************************************************************
         COMPLEX FUNCTION EPRES(IREA,XMP,XMPI,W,W0,GAM,X,L,J,ISO,     
     *   PIMP,GENER,TSTS11)                                             
         COMPLEX QETA,GAMMAT,Q0ETA                                      
         INTEGER L,J,ISO,IREA                                           
         REAL XMP,XMPI,W,W0,GAM,X,PIMP0                                 
         REAL PENER0,GENER0,XMETA,GAMMA                                 
         LOGICAL TSTS11                                                 
         COMPLEX IM                                                     
         DATA IM/(0.,1.)/                                               
         DATA XMETA/.5488/                                              
                   
         PENER0=(W0**2-XMP**2+XMPI**2)/(2.*W0)                          
         PIMP0 = SQRT ( PENER0**2 - XMPI**2 )                           
         GENER0 = (W0**2-XMP**2)/(2.*W0)                                
         if (PIMP0 .le. 0.)then
            write(6,*)' epres: pimp0 =',PIMP0
            PIMP0=abs(PIMP0)
         endif
                                                                       
         GAMMA=GAM * ((PIMP/PIMP0)**(2.*L+1))*                          
     *   (((PIMP0**2+X**2)/(PIMP**2+X**2))**L)                          
                                                                       
         GAMMAG=GAM * ((GENER/GENER0)**(2.*J))*                         
     *   (((GENER0**2+X**2)/(GENER**2+X**2))**J)                        
                                                                       
         IF(TSTS11)THEN                                                 
           XETA=(W**2-XMP**2+XMETA**2)/(2.*W)                           
           XETA=XETA**2-XMETA**2                                        
           IF(XETA.GE.0.)THEN                                           
          QETA=CMPLX(SQRT(XETA),0.)                                  
           ELSE                                                         
             QETA=CMPLX(0.,SQRT(-XETA))                                 
           END IF                                                       
          Q0ETA=(W0**2-XMP**2+XMETA**2)/(2.*W0)                        
          Q0ETA=SQRT(Q0ETA**2-XMETA**2)                                
          GAMMAT=GAM*(.55*(QETA/Q0ETA)+.4*(PIMP/PIMP0)+.05)            
C            GAMMAT=GAM*(PIMP/PIMP0)
           ELSE                                                           
            GAMMAT=CMPLX(GAMMA,0.)                                       
          END IF                                                         
         EPRES=SQRT((GENER0*PIMP0)/(GENER*PIMP))                        
         EPRES=EPRES*W0*SQRT(GAMMA*GAMMAG)                              
         EPRES=EPRES/(W0**2-W**2-IM*W0*GAMMAT)                          
         IF(IREA.EQ.1.OR.IREA.EQ.2)THEN                                              
           IF(ISO.EQ.1)THEN                                             
             EPRES=-EPRES/SQRT(2.)                                      
           ELSE                                                         
             EPRES=EPRES*SQRT(2.)                                       
           END IF                                                       
         END IF                                                         
99        RETURN                                                         
         END      

***************************************************************************
      SUBROUTINE EXPA(EPIREA,EPQ2)
      implicit none
c
c     Amplitudes from Experiment
c

      COMMON /IP/IT,IB,NF,IBORN,CUT,IP11

      real cut
      real pi,xmprot,xmpip,xmpi0,xmpi,epw
      real elchar,epq2,qme11,qmm11,qmm12,qme22,qmm21,qmm22,qmm23
      real fact,cpin,pvec_gam,pvec_pi,pvec_gam0,q2evf,q2evf0,dip_evf
      real q0_gam,ffac
      real cosdd,sindd

      integer it,ib,nf,iborn,ip11
      INTEGER EPIREA


       real
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710, 
     * WP13_1720, WD13_1520, WD13_1700, 
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600, 
     * WS31_1620, WS31_1900, WP31_1910,   
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710, 
     * LP13_1720, LD13_1520, LD13_1700, 
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600, 
     * LS31_1620, LS31_1900, LP31_1910,  
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710, 
     * PIBP13_1720, PIBD13_1520, PIBD13_1700, 
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600, 
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,  
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975

       COMMON/RCONST/      
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710, 
     * WP13_1720, WD13_1520, WD13_1700, 
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600, 
     * WS31_1620, WS31_1900, WP31_1910,   
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710, 
     * LP13_1720, LD13_1520, LD13_1700, 
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600, 
     * LS31_1620, LS31_1900, LP31_1910,  
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710, 
     * PIBP13_1720, PIBD13_1520, PIBD13_1700, 
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600, 
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,  
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975

         COMMON/RAMP/
     * RAS11_1535,             RCS11_1535, 
     * RAS11_1650,             RCS11_1650, 
     * RAP11_1440,             RCP11_1440, 
     * RAP11_1710,             RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,         
     * RAP13_1910, RBP13_1910, RCP13_1910,         
     * RAD13_1520, RBD13_1520, RCD13_1520, 
     * RAD13_1700, RBD13_1700, RCD13_1700, 
     * RAD15_1675, RBD15_1675, RCD15_1675, 
     * RAF15_1680, RBF15_1680, RCF15_1680, 
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,             RCS31_1620, 
     * RAS31_1900,             RCS31_1900, 
     * RAP31_1910,             RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232, 
     * RAP33_1920, RBP33_1920, RCP33_1920, 
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * RAP33_1600, RBP33_1600, RCP33_1600,
     * RAF17_1990, RBF17_1990, RCF17_1990,
     * RAF15_2000, RBF15_2000, RCF15_2000,
     * RAP11_2100, RBP11_2100, RCP11_2100,
     * RAF35_2000, RBF35_2000, RCF35_2000,
     * RAP13_1870, RBP13_1870, RCP13_1870,
     * RAP31_1925, RBP31_1975, RCP31_1975,
     * RAP13_1980, RBP13_1980, RCP13_1980,
     * RAF15_1955, RBF15_1955, RCF15_1955,
     * RAP13_1955, RBP13_1955, RCP13_1955,
     * RAP33_1975, RBP33_1975, RCP33_1975

         real
     * RAS11_1535,             RCS11_1535, 
     * RAS11_1650,             RCS11_1650, 
     * RAP11_1440,             RCP11_1440, 
     * RAP11_1710,             RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,         
     * RAP13_1910, RBP13_1910, RCP13_1910,         
     * RAD13_1520, RBD13_1520, RCD13_1520, 
     * RAD13_1700, RBD13_1700, RCD13_1700, 
     * RAD15_1675, RBD15_1675, RCD15_1675, 
     * RAF15_1680, RBF15_1680, RCF15_1680, 
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,             RCS31_1620, 
     * RAS31_1900,             RCS31_1900, 
     * RAP31_1910,             RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232, 
     * RAP33_1920, RBP33_1920, RCP33_1920, 
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * RAP33_1600, RBP33_1600, RCP33_1600,
     * RAF17_1990, RBF17_1990, RCF17_1990,
     * RAF15_2000, RBF15_2000, RCF15_2000,
     * RAP11_2100, RBP11_2100, RCP11_2100,
     * RAF35_2000, RBF35_2000, RCF35_2000,
     * RAP13_1870, RBP13_1870, RCP13_1870,
     * RAP31_1925, RBP31_1975, RCP31_1975,
     * RAP13_1980, RBP13_1980, RCP13_1980,
     * RAF15_1955, RBF15_1955, RCF15_1955,
     * RAP13_1955, RBP13_1955, RCP13_1955,
     * RAP33_1975, RBP33_1975, RCP33_1975


      REAL A12S11_1535,              S12S11_1535
      REAL A12S11_1650,              S12S11_1650
      REAL A12P11_1440,              S12P11_1440
      REAL A12P11_1710,              S12P11_1710
      REAL A12P13_1720, A32P13_1720, S12P13_1720
      REAL A12P13_1910, A32P13_1910, S12P13_1910
      REAL A12D13_1520, A32D13_1520, S12D13_1520
      REAL A12D13_1700, A32D13_1700, S12D13_1700
      REAL A12D15_1675, A32D15_1675, S12D15_1675
      REAL A12F15_1680, A32F15_1680, S12F15_1680
      REAL A12G17_2190, A32G17_2190, S12G17_2190
      REAL A12G19_2250, A32G19_2250, S12G19_2250
      REAL A12H19_2220, A32H19_2220, S12H19_2220
      REAL A12I111_2600,A32I111_2600,S12I111_2600
      REAL A12S31_1620,              S12S31_1620
      REAL A12S31_1900,              S12S31_1900
      REAL A12P31_1910,              S12P31_1910
      REAL A12P33_1232, A32P33_1232, S12P33_1232
      REAL A12P33_1920, A32P33_1920, S12P33_1920
      REAL A12D33_1700, A32D33_1700, S12D33_1700
      REAL A12D35_1930, A32D35_1930, S12D35_1930
      REAL A12F35_1905, A32F35_1905, S12F35_1905
      REAL A12F37_1950, A32F37_1950, S12F37_1950
      REAL A12H311_2420,A32H311_2420,S12H311_2420
      REAL THETA_S11
c      REAL THETA_D13
c      REAL THETA_D15

      DATA A12S11_1535,              S12S11_1535  /0.051  ,0.   / 
      DATA A12S11_1650,              S12S11_1650  /0.061  ,0.   / 
      DATA A12P11_1440,              S12P11_1440  /-0.071 ,0. / 
      DATA A12P11_1710,              S12P11_1710  /-0.039  ,0.   /
      DATA A12P13_1720, A32P13_1720, S12P13_1720  /0.029  ,-0.034 ,0./         
      DATA A12P13_1910, A32P13_1910, S12P13_1910  /0.  ,0.   ,0./         
      DATA A12D13_1520, A32D13_1520, S12D13_1520  /-0.021  ,0.164 ,0./ 
      DATA A12D13_1700, A32D13_1700, S12D13_1700  /0.  ,0.   ,0./ 
      DATA A12D15_1675, A32D15_1675, S12D15_1675  /0.004  ,0.025  ,0./ 
      DATA A12F15_1680, A32F15_1680, S12F15_1680  /-0.016  ,0.146 ,0./ 
      DATA A12G17_2190, A32G17_2190, S12G17_2190  /0.  ,0.   ,0./
      DATA A12G19_2250, A32G19_2250, S12G19_2250  /0.  ,0.   ,0./
      DATA A12H19_2220, A32H19_2220, S12H19_2220  /0.  ,0.   ,0./
      DATA A12I111_2600,A32I111_2600,S12I111_2600 /0.  ,0.   ,0./
      DATA A12S31_1620,              S12S31_1620  /0.040  ,0.   / 
      DATA A12S31_1900,              S12S31_1900  /0.  ,0.   / 
      DATA A12P31_1910,              S12P31_1910  /0.033  ,0.   /
      DATA A12P33_1232, A32P33_1232, S12P33_1232  /-0.143,-0.259,0./ 
      DATA A12P33_1920, A32P33_1920, S12P33_1920  /0.  ,0.   ,0.   /    
      DATA A12D33_1700, A32D33_1700, S12D33_1700  /0.110,0.093 ,0. /
      DATA A12D35_1930, A32D35_1930, S12D35_1930  /-0.021,0.004,0. /
      DATA A12F35_1905, A32F35_1905, S12F35_1905  /0.049 ,-0.011,0./
      DATA A12F37_1950, A32F37_1950, S12F37_1950  /-0.090,-0.113,0./
      DATA A12H311_2420,A32H311_2420,S12H311_2420 /0.  ,0.   ,0.   /
                                                                        
      DATA THETA_S11 / 38. /
c      DATA THETA_D13 /  0. /
c      DATA THETA_D15 /  0. /


c      DATA SQ2        /1.41421  /     
      DATA PI         /3.14159  /     
      DATA XMPROT     /.938     /     
      DATA XMPIP      /.1395    /     
      DATA XMPI0      /.1349    /     
c      DATA XMETA      /.5488    /     




C       Now calculate the resonant helicity amplitudes 
C
         IF(EPIREA.EQ.1.OR.EPIREA.EQ.2)XMPI=XMPI0
         IF(EPIREA.EQ.3.OR.EPIREA.EQ.4)XMPI=XMPIP
         ELCHAR = SQRT(4.*PI/137.)

        IF (EPIREA.EQ.1.OR.EPIREA.EQ.3) THEN
        CALL QKEM(LS11_1535,EPQ2,QME11,QMM11,QMM12)
        A12S11_1535 = ELCHAR/2./SQRT(XMPROT*(LS11_1535**2-XMPROT**2))* 
     *       (SQRT(1./3.)*QME11 - SQRT(2./3.)*QMM11)*COSDD(THETA_S11)
        CALL QKEM(LS11_1650,EPQ2,QME11,QMM11,QMM12)
        A12S11_1650 = ELCHAR/2./SQRT(XMPROT*(LS11_1650**2-XMPROT**2)) * 
     *        (SQRT(1./3.)*QME11 - SQRT(2./3.)*QMM11)*SINDD(THETA_S11)
  
C
C        Assume no mixing for D13 and D15
C
        CALL QKEM(LD13_1520,EPQ2,QME11,QMM11,QMM12)
        A12D13_1520=ELCHAR/2./SQRT(XMPROT*(LD13_1520**2-XMPROT**2))*
     *   (SQRT(1./6.)*QME11 + SQRT(1./12.)*QMM11 - SQRT(3./4.)*QMM12)
        A32D13_1520=ELCHAR/2./SQRT(XMPROT*(LD13_1520**2-XMPROT**2))*
     *           (SQRT(1./2.)*QME11 + 0.5*QMM11 + 0.5*QMM12)

         A12D13_1700 = 0. 
         A32D13_1700 = 0. 

         A12D15_1675 = 0. 
         A32D15_1675 = 0. 

        CALL QKEM(LD33_1700,EPQ2,QME11,QMM11,QMM12)
        A12D33_1700 = -ELCHAR/2./SQRT(XMPROT*(LD33_1700**2-XMPROT**2))*
     *   (-SQRT(1./12.) * QMM12 - SQRT(1./6.) * QME11 +
     *    0.5*SQRT(1./27.)*QMM11)
        A32D33_1700 = -ELCHAR/2./SQRT(XMPROT*(LD33_1700**2-XMPROT**2))*
     *        (1./6.*QMM12 - SQRT(1./2.)*QME11 + 1./6.*QMM11)

        CALL QKEM(LS31_1620,EPQ2,QME11,QMM11,QMM12)
        A12S31_1620 = ELCHAR/2./SQRT(XMPROT*(LS31_1620**2-XMPROT**2))*
     *        ( SQRT(1./3.)*QME11 + SQRT(2./27.)*QMM11)

C*  SET ALL TO ZERO EXCEPT S11(1535) AND D13(1520)
C*         A12S11_1650=0.
C*         A12S31_1620=0.
C*         A12D33_1700=0.
C*         A32D33_1700=0.
        

C
C         Now do the {56,0+} ---> {56,2+} transition.
C         Here we use the quark multipole moments of Cottingham & Dunbar
C         directly. Note thet quark multipole QMM21 is assumed = 0.
C         This assumption has no justification (data do not allow its
C         determination).
C                  
          CALL QKEM2(LF15_1680,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LF15_1680**2 - XMPROT**2))
          A12F15_1680 = FACT*(-SQRT(2./3.)*QMM23 + SQRT(1./5.)*QME22+
     *                  SQRT(2./15.)*QMM22)
          A32F15_1680 = FACT*(SQRT(1./3.)*QMM23 + SQRT(2./5.)*QME22+
     *                  SQRT(4./15.)*QMM22)

          CALL QKEM2(LP13_1720,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LP13_1720**2 - XMPROT**2))
          A12P13_1720 = -FACT*(SQRT(3./4.)*QMM21 + SQRT(1./10.)*QME22-
     *                   SQRT(3./20.)*QMM22)
          A32P13_1720 =  FACT*(0.5*QMM21 - SQRT(3./10.)*QME22 + 
     *                   SQRT(9./20.)*QMM22)

          CALL QKEM2(LF35_1905,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LF35_1905**2 - XMPROT**2))
          A12F35_1905 =  FACT*(SQRT(16./189.)*QMM23 + 
     *                  SQRT(28./135.)*QMM22)
          A32F35_1905 = -FACT*(SQRT(8./189.)*QMM23 - 
     *                  SQRT(56./135.)*QMM22)

          CALL QKEM2(LP33_1920,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LP33_1920**2 - XMPROT**2))
          A12P33_1920 = -FACT*(SQRT(3./15.)*QMM22 + 1./3.*QMM21)
          A32P33_1920 =  FACT*(SQRT(1./15.)*QMM22 - SQRT(1./3.)*QMM21)

          CALL QKEM2(LF37_1950,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LF37_1950**2 - XMPROT**2))
          A12F37_1950 = -FACT*SQRT(2./7.)*QMM23
          A32F37_1950 = -FACT*SQRT(10./21.)*QMM23
cccc
cccc        multiply A12F37_1950 and A32F37_1950 by factor 2 (temporarily).
cccc
          A12F37_1950 = A12F37_1950 * 2.
          A32F37_1950 = A32F37_1950 * 2.
cccc
cccc
cccc
          CALL QKEM2(LP31_1910,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LP31_1910**2 - XMPROT**2))
          A12P31_1910 = -FACT*2./3.*QMM21

C*   SET ALL STATES ZERO EXCEPT F15(1680)
C          A12P13_1720=0.
C          A32P13_1720=0.
C          A12F37_1950=0.
C          A32F37_1950=0.
C          A12F35_1905=0.
C          A32F35_1905=0.
C          A12P33_1920=0.
C          A32P33_1920=0.
C          A12P31_1910=0.
C
          ELSE
C* THE FOLLOWING IS FOR NEUTRON TARGET
        A12P11_1440=0.056

        CALL QKEM(LS11_1535,EPQ2,QME11,QMM11,QMM12)
        A12S11_1535 = ELCHAR/2./SQRT(XMPROT*(LS11_1535**2-XMPROT**2))* 
     *       ((-SQRT(1./3.)*QME11 +SQRT(2./27.)*QMM11)*COSDD(THETA_S11)
     *      +SINDD(THETA_S11)*(-SQRT(2./27.)*QMM11))
        CALL QKEM(LS11_1650,EPQ2,QME11,QMM11,QMM12)
        A12S11_1650 = ELCHAR/2./SQRT(XMPROT*(LS11_1650**2-XMPROT**2)) *
     *        ((-SQRT(1./3.)*QME11 +SQRT(2./27.)*QMM11)*
     *        (SINDD(THETA_S11))
     *      +COSDD(THETA_S11)*(-SQRT(2./27.)*QMM11))
C
C        Assume no mixing for D13 and D15
C
        CALL QKEM(LD13_1520,EPQ2,QME11,QMM11,QMM12)
        A12D13_1520=ELCHAR/2./SQRT(XMPROT*(LD13_1520**2-XMPROT**2))*
     *   (-SQRT(1./6.)*QME11-SQRT(1./108.)*QMM11 +SQRT(1./12.)*QMM12)
        A32D13_1520=ELCHAR/2./SQRT(XMPROT*(LD13_1520**2-XMPROT**2))*
     *           (-SQRT(1./2.)*QME11 - 1./6.*QMM11 -1./6.*QMM12)


        CALL QKEM(LD13_1700,EPQ2,QME11,QMM11,QMM12)
         A12D13_1700 =ELCHAR/2./SQRT(XMPROT*(LD13_1700**2-XMPROT**2))* 
     *    (-SQRT(5./54.)*QMM11-SQRT(1./120.)*QMM12)
         A32D13_1700 =ELCHAR/2./SQRT(XMPROT*(LD13_1700**2-XMPROT**2))*
     *    (-SQRT(5./18.)*QMM11+SQRT(1./360.)*QMM12)

        CALL QKEM(LD15_1675,EPQ2,QME11,QMM11,QMM12)
         A12D15_1675 =ELCHAR/2./SQRT(XMPROT*(LD15_1675**2-XMPROT**2))* 
     *    (-SQRT(2./15.)*QMM12)
         A32D15_1675 =ELCHAR/2./SQRT(XMPROT*(LD15_1675**2-XMPROT**2))*
     *     (-SQRT(4./15.)*QMM12)
        CALL QKEM(LD33_1700,EPQ2,QME11,QMM11,QMM12)
        A12D33_1700 = -ELCHAR/2./SQRT(XMPROT*(LD33_1700**2-XMPROT**2))*
     *   (-SQRT(1./12.) * QMM12 - SQRT(1./6.) * QME11 +
     *    0.5*SQRT(1./27.)*QMM11)
        A32D33_1700 = -ELCHAR/2./SQRT(XMPROT*(LD33_1700**2-XMPROT**2))*
     *        (1./6.*QMM12 - SQRT(1./2.)*QME11 + 1./6.*QMM11)

        CALL QKEM(LS31_1620,EPQ2,QME11,QMM11,QMM12)
        A12S31_1620 = ELCHAR/2./SQRT(XMPROT*(LS31_1620**2-XMPROT**2))*
     *        ( SQRT(1./3.)*QME11 + SQRT(2./27.)*QMM11)

C*  SET ALL TO ZERO EXCEPT S11(1535) AND D13(1520)
C*         A12S11_1650=0.
C*         A12S31_1620=0.
C*         A12D33_1700=0.
C*         A32D33_1700=0.
        

C
C         Now do the {56,0+} ---> {56,2+} transition.
C         Here we use the quark multipole moments of Cottingham & Dunbar
C         directly. Note thet quark multipole QMM21 is assumed = 0.
C         This assumption has no justification (data do not allow its
C         determination).
C                  

          CALL QKEM2(LF15_1680,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LF15_1680**2 - XMPROT**2))
          A12F15_1680 = FACT*(SQRT(8./27.)*QMM23-SQRT(8./135.)*QMM22)
          A32F15_1680 = FACT*(-SQRT(4./27.)*QMM23-SQRT(16./135.)*QMM22)

          CALL QKEM2(LP13_1720,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LP13_1720**2 - XMPROT**2))
          A12P13_1720 = -FACT*(-SQRT(1./3.)*QMM21 + SQRT(1./15.)*QMM22)
          A32P13_1720 = -FACT*(1./3.*QMM21+SQRT(1./5.)*QMM22)  

          CALL QKEM2(LF35_1905,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LF35_1905**2 - XMPROT**2))
          A12F35_1905 =  FACT*(SQRT(16./189.)*QMM23 + 
     *                  SQRT(28./135.)*QMM22)
          A32F35_1905 = -FACT*(SQRT(8./189.)*QMM23 - 
     *                  SQRT(56./135.)*QMM22)

          CALL QKEM2(LP33_1920,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LP33_1920**2 - XMPROT**2))
          A12P33_1920 = -FACT*(SQRT(3./15.)*QMM22 + 1./3.*QMM21)
          A32P33_1920 =  FACT*(SQRT(1./15.)*QMM22 - SQRT(1./3.)*QMM21)

          CALL QKEM2(LF37_1950,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LF37_1950**2 - XMPROT**2))
          A12F37_1950 = -FACT*SQRT(2./7.)*QMM23
          A32F37_1950 = -FACT*SQRT(10./21.)*QMM23
cccc
cccc        multiply A12F37_1950 and A32F37_1950 by factor 2 (temporarily !!!)
cccc        to bring them in quantitative agreement with photoproduction
cccc        analysis
cccc
          A12F37_1950 = A12F37_1950 * 2.
          A32F37_1950 = A32F37_1950 * 2.
cccc
cccc
cccc
          CALL QKEM2(LP31_1910,EPQ2,QME22,QMM21,QMM22,QMM23)
          FACT = ELCHAR/2./SQRT(XMPROT*(LP31_1910**2 - XMPROT**2))
          A12P31_1910 = -FACT*2./3.*QMM21
      
          END IF

C
C       Now calculate the resonant partial wave helicity elements
C

          CPIN=SQRT(2./3.)
          IF(EPIREA.EQ.1.OR.EPIREA.EQ.3)CPIN=-CPIN
98         PVEC_GAM = SQRT((EPQ2 + (LS11_1535+XMPROT)**2)*
     *   (EPQ2+(LS11_1535-XMPROT)**2))/(2.*LS11_1535)
         PVEC_PI= SQRT((LS11_1535**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LS11_1535**2) - XMPI**2)
         FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LS11_1535*PIBS11_1535/WS11_1535)
         RAS11_1535 =-FACT*A12S11_1535*19.73*CPIN
         RCS11_1535 = 0.

         PVEC_GAM = SQRT((EPQ2 + (LS11_1650+XMPROT)**2)*
     *   (EPQ2+(LS11_1650-XMPROT)**2))/(2.*LS11_1650)
         PVEC_PI= SQRT((LS11_1650**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LS11_1650**2) - XMPI**2)
         FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LS11_1650*PIBS11_1650/WS11_1650)
         RAS11_1650 = -FACT*A12S11_1650*19.73*CPIN
         RCS11_1650 = 0.

         PVEC_GAM = SQRT((EPQ2 + (LP11_1440+XMPROT)**2)*
     *   (EPQ2+(LP11_1440-XMPROT)**2))/(2.*LP11_1440)
         PVEC_GAM0 = SQRT(((LP11_1440+XMPROT)**2)*
     *   ((LP11_1440-XMPROT)**2))/(2.*LP11_1440)
         PVEC_PI= SQRT((LP11_1440**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LP11_1440**2) - XMPI**2)
         FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LP11_1440*PIBP11_1440/WP11_1440)
         EPW=LP11_1440
         Q2EVF=((EPW**2-XMPROT**2)**2+EPQ2*(EPW+XMPROT)**2)/
     *   (4.*EPW*XMPROT)
         Q2EVF0=((EPW**2-XMPROT**2)**2)/(4.*EPW*XMPROT)
         DIP_EVF = (1.+ Q2EVF/0.71)**(-2.)
         
         RAP11_1440 = FACT*A12P11_1440*(Q2EVF/Q2EVF0)*
     *   DIP_EVF*19.73*CPIN
         RCP11_1440 =               0.0000D+00

         RAP11_1710 =               0.0000D+00
         RCP11_1710 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LP13_1720+XMPROT)**2)*
     *   (EPQ2+(LP13_1720-XMPROT)**2))/(2.*LP13_1720)
         PVEC_PI= SQRT((LP13_1720**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LP13_1720**2) - XMPI**2)
         FACT = SQRT(1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LP13_1720*PIBP13_1720/WP13_1720)
         RAP13_1720 = - FACT*A12P13_1720*19.73*CPIN
         RBP13_1720 = FACT*SQRT(16./12.)*A32P13_1720*19.73*
     *                 CPIN
         RCP13_1720 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LD13_1520+XMPROT)**2)*
     *   (EPQ2+(LD13_1520-XMPROT)**2))/(2.*LD13_1520)
         PVEC_PI= SQRT((LD13_1520**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LD13_1520**2) - XMPI**2)
         FACT = SQRT(1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LD13_1520*PIBD13_1520/WD13_1520)
         RAD13_1520 = FACT*A12D13_1520*19.73*CPIN
         RBD13_1520 =  -FACT*SQRT(16./12.)*A32D13_1520*19.73*
     *                 CPIN
         RCD13_1520 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LD13_1700+XMPROT)**2)*
     *   (EPQ2+(LD13_1700-XMPROT)**2))/(2.*LD13_1700)
         PVEC_PI= SQRT((LD13_1700**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LD13_1700**2) - XMPI**2)
         FACT = SQRT(1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LD13_1700*PIBD13_1700/WD13_1700)
         RAD13_1700 = FACT*A12D13_1700*19.73*CPIN
         RBD13_1700 = -FACT*SQRT(16./12.)*A32D13_1700*19.73*CPIN
         RCD13_1700 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LD15_1675+XMPROT)**2)*
     *   (EPQ2+(LD15_1675-XMPROT)**2))/(2.*LD15_1675)
         PVEC_PI= SQRT((LD15_1675**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LD15_1675**2) - XMPI**2)
         FACT = SQRT(1./6./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LD15_1675*PIBD15_1675/WD15_1675)
         RAD15_1675 =-FACT*A12D15_1675*19.73*CPIN
         RBD15_1675 = FACT*SQRT(16./32.)*A32D15_1675*19.73*CPIN
         RCD15_1675 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LF15_1680+XMPROT)**2)*
     *   (EPQ2+(LF15_1680-XMPROT)**2))/(2.*LF15_1680)
         PVEC_PI= SQRT((LF15_1680**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LF15_1680**2) - XMPI**2)
         FACT = SQRT(1./6./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LF15_1680*PIBF15_1680/WF15_1680)
         RAF15_1680 = FACT*A12F15_1680*19.73*CPIN
         RBF15_1680 = -FACT*SQRT(16./32.)*A32F15_1680*19.73*
     *                 CPIN

         RCF15_1680 =              0.0000D+00

         RAG17_2190 =              0.0000D+00
         RBG17_2190 =              0.0000D+00
         RCG17_2190 =              0.0000D+00

         RAG19_2250 =              0.0000D+00
         RBG19_2250 =              0.0000D+00
         RCG19_2250 =              0.0000D+00

         RAH19_2220 =              0.0000D+00
         RBH19_2220 =              0.0000D+00
         RCH19_2220 =              0.0000D+00

         RAI111_2600 =             0.0000D+00
         RBI111_2600 =             0.0000D+00
         RCI111_2600 =             0.0000D+00
C          
         PVEC_GAM = SQRT((EPQ2 + (LS31_1620+XMPROT)**2)*
     *   (EPQ2+(LS31_1620-XMPROT)**2))/(2.*LS31_1620)
         PVEC_PI= SQRT((LS31_1620**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LS31_1620**2) - XMPI**2)
         FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LS31_1620*PIBS31_1620/WS31_1620)
         RAS31_1620 =  FACT*A12S31_1620*19.73*SQRT(1./3.)
         RCS31_1620 =              0.0000D+00

         RAS31_1900 =              0.0000D+00
         RCS31_1900 =              0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LP31_1910+XMPROT)**2)*
     *   (EPQ2+(LP31_1910-XMPROT)**2))/(2.*LP31_1910)
         PVEC_PI= SQRT((LP31_1910**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LP31_1910**2) - XMPI**2)
         FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LP31_1910*PIBP31_1910/WP31_1910)
         RAP31_1910 = -FACT*A12P31_1910*19.73*SQRT(1./3.)
         RCP31_1910 =              0.0000D+00

         
         PVEC_GAM = SQRT((EPQ2 + (LP33_1232+XMPROT)**2)*
     *   (EPQ2+(LP33_1232-XMPROT)**2))/(2.*LP33_1232)
         PVEC_GAM0 = SQRT(((LP33_1232+XMPROT)**2)*
     *   ((LP33_1232-XMPROT)**2))/(2.*LP33_1232)
         Q0_GAM=(LP33_1232**2 - XMPROT**2+EPQ2)/(2.*LP33_1232)
         PVEC_PI= SQRT((LP33_1232**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LP33_1232**2) - XMPI**2)
           FACT = SQRT(1./137.*LP33_1232*PIBP33_1232/6./XMPROT/PVEC_PI/
     *              XMPROT**2/WP33_1232)*PVEC_GAM
         FFAC = 1./(1+EPQ2/0.71)**2*EXP(-0.21*EPQ2)
         RAP33_1232 =  -0.5*3.1*FFAC*FACT*SQRT(1./3.)*19.73
         RBP33_1232 =       3.1*FFAC*FACT*SQRT(1./3.)*19.73
         
C*  SET THE AMPLITUDE OF P33 TO ZERO
ccc          RAP33_1232=0.
ccc          RBP33_1232=0.
c         RCP33_1232 =  -0.045*0.07*2.*SQRT(EPQ2)/Q0_GAM*3.*FFAC*
c     *                  FACT*SQRT(1./3.)*19.73
 
         
         PVEC_GAM = SQRT(EPQ2+ (LP33_1920+XMPROT)**2)*
     *   (EPQ2+(LP33_1920-XMPROT)**2)/(2.*LP33_1920)
         PVEC_PI= SQRT(((LP33_1920**2+XMPI**2-XMPROT**2)/
     *   (4*LP33_1920**2) - XMPI**2))
         FACT = SQRT((1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LP33_1920*PIBP33_1920/WP33_1920))
         RAP33_1920 =  FACT*A12P33_1920*SQRT(1./3.)
         RBP33_1920 = -FACT*SQRT(16./12.)*A32P33_1920*SQRT(1./3.)
         RCP33_1920 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LD33_1700+XMPROT)**2)*
     *   (EPQ2+(LD33_1700-XMPROT)**2))/(2.*LD33_1700)
         PVEC_PI= SQRT((LD33_1700**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LD33_1700**2) - XMPI**2)
         FACT = SQRT(1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LD33_1700*PIBD33_1700/WD33_1700)
         RAD33_1700 = -FACT*A12D33_1700*19.73*SQRT(1./3.)
         RBD33_1700 =  FACT*SQRT(16./12.)*A32D33_1700*19.73*
     *                SQRT(1./3.)
         RCD33_1700 =            0.0000D+00
         
         RAD35_1930 =               0.0000D+00
         RBD35_1930 =               0.0000D+00
         RCD35_1930 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LF35_1905+XMPROT)**2)*
     *   (EPQ2+(LF35_1905-XMPROT)**2))/(2.*LF35_1905)
         PVEC_PI= SQRT((LF35_1905**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LF35_1905**2) - XMPI**2)
         FACT = SQRT(1./6./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LF35_1905*PIBF35_1905/WF35_1905)
         RAF35_1905 = -FACT*A12F35_1905*19.73*SQRT(1./3.)
         RBF35_1905 =  FACT*SQRT(16./32.)*A32F35_1905*19.73*
     *                 SQRT(1./3.)
         RCF35_1905 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LF37_1950+XMPROT)**2)*
     *   (EPQ2+(LF37_1950-XMPROT)**2))/(2.*LF37_1950)
         PVEC_PI= SQRT((LF37_1950**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LF37_1950**2) - XMPI**2)
         FACT = SQRT(1./8./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LF37_1950*PIBF37_1950/WF37_1950)
         RAF37_1950 =  FACT*A12F37_1950*19.73*SQRT(1./3.)
         RBF37_1950 = -FACT*SQRT(16./60.)*A32F37_1950*19.73*
     *                 SQRT(1./3.)
	 RCF37_1950 =               0.0000D+00

         RAH311_2420 =              0.0000D+00
         RBH311_2420 =              0.0000D+00
         RCH311_2420 =              0.0000D+00
c      rm1m=-0.19
c      RAP11_1440=RM1M/sqrt(0.60)
c      RE0P=-0.11
c      RAS11_1535=RE0P/sqrt(0.40)
c      RE0P=0.090
c      RAS11_1650=RE0P/sqrt(0.60)
c      RE2M=-0.55
c      RM2M=-0.16
c      RAD13_1520=(3./2.*RM2M-1./2.*RE2M)/sqrt(0.55)
c      RBD13_1520=(RE2M+RM2M)/sqrt(0.55) 
c      RE2M=-0.034
c      RM2M=-0.016   
c      RAD13_1700=(3./2.*RM2M-1./2.*RE2M)/sqrt(0.10)
c      RBD13_1700=(RE2M+RM2M)/sqrt(0.10)
c      RE3M=0.039
c      RM3M=0.0
c      RAF15_1680=(2.*RM3M-RE3M)/sqrt(0.60)
c      RBF15_1680=(RE3M+RM3M)/sqrt(0.60)
c      RE3M=0.23
c      RM3M=-0.17
c      RAD33_1700=(3./2.*RM2M-1./2.*RE2M)/sqrt(0.15)
c      RBD33_1700=(RE2M+RM2M)/sqrt(0.15)

         RETURN
         END
***************************************************************************
***************************************************************************
* DEVELOPED IN END OF SEPT. 1992 Z.LI
* With IFLAG and VFLAG both equal to 0, we calculate electric born only
* Both neutron and proton form factor took from original program hborn.
* The axial form factor may not correct.
* A cut off factor is added for the born terms when the center of
* mass energy is larger than 1.3 GEV:
* Fcut=cut**2/(cut**2-(w-1.3)**2), cut is set to 0.4 
**************************************************************************
       SUBROUTINE BORNT(Q2,W,SINX,COSX,SINX2,COSX2,BORN,IR,CUT)
       implicit none
       real q2,w,sinx,cosx,sinx2,cosx2,cut
       integer ir
       real xmpip,xmpi0,xmpi,xmp,up,un,ec,gn
       real BORN(6),GAMMA(6,3)
       real A(6,3),AR(6),FR(6),F1(3),F2(3),FP(3),F2A(3)
       real C(6),B(6,6),ETA(6),EPS(3)
       real pi
       real s,q0,e1,e2,p,p0,t,u,sq2
       real gg,gep,cap,f1p,f1n,f2p,f2n,f1v,f1s,f2v,f2s
       real q,z1,z2,fc,pq,fcut
       integer iflag
       real vflag
       integer i,j

       DATA XMPIP/0.139563/,XMPI0/0.1349630/,XMP/0.93827/
       DATA UP,UN/1.793,-1.913/
       DATA EC/.302862/,GN/13.5/
       DATA ETA/1.,1.,-1.,1.,-1.,-1./,EPS/1.,1.,-1./
       DATA IFLAG/0/,VFLAG/0./
C* IFLAG=0, VFLAG=0. WOULD GET RIDE OF MAGNETIC AND PV EXTRA TERM 
C* RESPECTIVELY.
       IF(IR.EQ.1.OR.IR.EQ.2)XMPI=XMPI0
       IF(IR.EQ.3.OR.IR.EQ.4)XMPI=XMPIP
       PI=4.0*ATAN(1.0)
       GG=SQRT(4.*PI*GN)
C*  THE PI-N COUPLING CONSTANT
       GEP=1./(1.+3.04*Q2+1.54*Q2*Q2+0.068*Q2*Q2*Q2)
       CAP=Q2/(4.*XMP*XMP)
       F1P=EC*GEP*(1.+(1.+UP)*CAP)/(1.+CAP)
       F1N=EC*GEP*UN*CAP/(1.+CAP)
       IF(IFLAG.EQ.0) THEN
         F2P=0.
         F2N=0.
       ELSE
         F2P=EC/(2.*XMP)*GEP*UP/(1.+CAP)
         F2N=EC/(2.*XMP)*GEP*UN/(1.+CAP)
       END IF
       F1V=F1P-F1N
       F1S=F1P+F1N
       F2V=F2P-F2N
       F2S=F2P+F2N
       F1(1)=F1V
       F1(2)=F1S
       F1(3)=F1V
       F2(1)=F2V
       F2(2)=F2S
       F2(3)=F2V
C* 1,3,2, CORESSPONDING TO THE ISOVECTOR (+,-) AND ISOSCALAR (0).
       FP(3)=EC/(1.+Q2/0.5)
       FP(1)=0.
       FP(2)=0.
C* PION FORM FACTOR.
       S=W**2
       Q0=(S-XMP**2-Q2)/(2.*W)
       Q=SQRT(Q2+Q0**2)
       E1=(S+XMP**2+Q2)/(2.*W)
       E2=(S+XMP**2-XMPI**2)/(2.*W)
       P=SQRT(E2**2-XMP**2)
       P0=SQRT(P**2+XMPI**2)
       T=2.*Q*P*COSX-2.*Q0*P0+XMPI**2-Q2
       U=-2.*Q*P*COSX-2.*Q0*E2+XMP**2-Q2
       DO 15,J=1,3
         GAMMA(1,J)=0.5*GG*F1(J)
         GAMMA(2,J)=-GG*F1(J)/(T-XMPI**2)
         GAMMA(3,J)=-0.5*GG*F2(J)
         GAMMA(4,J)=GAMMA(3,J)
         GAMMA(5,J)=0.5*GAMMA(2,J)
         GAMMA(6,J)=0.
15     CONTINUE   
       DO 25,I=1,4
       DO 35,J=1,3
         A(I,J)=(1./(S-XMP**2)+EPS(J)*ETA(I)/(U-XMP**2))*GAMMA(I,J)
35     CONTINUE
25     CONTINUE
       F2A(1)=F2(1)
       F2A(2)=F2(2)
       A(1,1)=VFLAG*0.5*GG/XMP*F2A(1)+A(1,1)
       A(1,2)=VFLAG*0.5*GG/XMP*F2A(2)+A(1,2)
c* Add in the pseudovector extra term, here we simply used F2(1,2) same
c* as in non-axial term.
       DO 45,J=1,3
       IF (Q2.EQ.0.) THEN
         A(5,J)=0.
         A(6,J)=0.
       ELSE
         A(5,J)=(1./(S-XMP**2)+ETA(5)/(U-XMP**2))*GAMMA(5,J)
     C   +0.5*(1.-EPS(J))*2.*GG/Q2*(FP(J)-F1(J))/(T-XMPI**2)
         A(6,J)=(1./(S-XMP**2)+ETA(6)/(U-XMP**2))*GAMMA(6,J)
       END IF
45     CONTINUE
       SQ2=SQRT(2.0)
       DO 50, I=1,6
         IF(IR.EQ.1) AR(I)=(A(I,2)+A(I,1))
         IF(IR.EQ.2) AR(I)=(-A(I,2)+A(I,1))
         IF(IR.EQ.3) AR(I)=-SQ2*(A(I,2)+A(I,3))
         IF(IR.EQ.4) AR(I)=SQ2*(A(I,2)-A(I,3))
50     CONTINUE
       Z1=SQRT(E1+XMP)
       Z2=SQRT(E2+XMP)
       FC=1./(8.*PI*W)
       C(1)=FC*(W-XMP)*Z1*Z2
       C(2)=FC*(W+XMP)*P*Q/Z1/Z2
       C(3)=FC*(W+XMP)*P*Q*Z2/Z1
       C(4)=FC*(W-XMP)*P**2*Z1/Z2
       C(5)=FC*Z1/Z2*P
       C(6)=FC*Q*Z2/Z1  

c     The following is not used:
c       XX=C(1)/C(2)

       PQ=(T-XMPI**2+Q2)/2.
       B(1,1)=1.
       B(1,2)=0.
       B(1,3)=-PQ/(W-XMP)
       B(1,4)=W-XMP+PQ/(W-XMP)
       B(1,5)=0.
       B(1,6)=Q2/(W-XMP)
       B(2,1)=-1.
       B(2,2)=0.
       B(2,3)=-PQ/(W+XMP)
       B(2,4)=W+XMP+PQ/(W+XMP)
       B(2,5)=0.
       B(2,6)=Q2/(W+XMP)
       B(3,1)=0.
       B(3,2)=(S-XMP**2+Q2/2.)/(W+XMP)
       B(3,3)=1.
       B(3,4)=-1.
       B(3,5)=-Q2/(W+XMP)
       B(3,6)=0.
       B(4,1)=0.
       B(4,2)=-(S-XMP**2+Q2/2.)/(W-XMP)
       B(4,3)=1.
       B(4,4)=-1.
       B(4,5)=Q2/(W-XMP)
       B(4,6)=0.
       B(5,1)=-E1+XMP
       B(5,2)=1./(2.*Q0)*(Q**2*(3.*PQ+2.*Q0*W)-(PQ+P0*Q0)*
     C   (2.*(S-XMP**2)+Q2))
       B(5,3)=PQ+P0*(W-XMP)
       B(5,4)=-PQ-P0*(W+XMP)+(E1-XMP)*(W+XMP)
       B(5,5)=P0*Q2-Q0*PQ
       B(5,6)=-(E1-XMP)*(W+XMP)
       B(6,1)=E1+XMP
       B(6,2)=-B(5,2)
       B(6,3)=PQ+P0*(W+XMP)
       B(6,4)=-PQ-P0*(W+XMP)+(E1+XMP)*(W-XMP)
       B(6,5)=-B(5,5)
       B(6,6)=-(E1+XMP)*(W-XMP)
       DO 55, I=1,6
       FR(I)=0.
       DO 65, J=1,6
         FR(I)=FR(I)+C(I)*B(I,J)*AR(J)*19.732857
65     CONTINUE
55     CONTINUE
       IF(W.LE.1.3) FCUT=1.
       IF(W.GT.1.3) FCUT=CUT**2/(CUT**2+(W-1.3)**2)
c       fcut=1
       BORN(1)=-SINX*COSX2/SQ2*(FR(3)+FR(4))
       BORN(2)=-SQ2*COSX2*(FR(1)-FR(2))+SINX2*SINX*(FR(3)-FR(4))/SQ2
       BORN(3)=SINX2*SINX/SQ2*(FR(3)-FR(4))
       BORN(4)=SQ2*SINX2*(FR(1)+FR(2))+COSX2*SINX*(FR(3)+FR(4))/SQ2 
       BORN(5)=COSX2*(FR(5)+FR(6))*Q0/Q
       BORN(6)=SINX2*(FR(5)-FR(6))*Q0/Q
       DO 75, I=1,6
       BORN(I)=BORN(I)*FCUT
75     CONTINUE
       RETURN
       END   
c*********************************************************************
         SUBROUTINE QKEM(EPW,EPQ2,QME11,QMM11,QMM12)        
         implicit none
         real epw,epq2,qme11,qmm11,qmm12,xmprot
         real q2evf,q2evf0,dip_evf,dip_evf0
         DATA XMPROT/0.938/
         Q2EVF=((EPW**2-XMPROT**2)**2+EPQ2*(EPW+XMPROT)**2)/
     *   (4.*EPW*XMPROT)
         Q2EVF0=((EPW**2-XMPROT**2)**2)/(4.*EPW*XMPROT)
         DIP_EVF = (1.+ Q2EVF/0.71)**(-2.)
         dip_evf0 = (1.+ Q2EVF0/0.71)**(-2.)
C* PUT THE FOLLOWING TO USE THE COUPLING CONSTANT GIVEN IN DATA STATEMENT
C          GO TO 98
                  QME11 = 3.  * DIP_EVF
                  QMM11 =  3.8 * (0.4 - Q2EVF) * DIP_EVF
         IF(Q2EVF.LE.1.2)THEN
                  QMM12 = 5.* Q2EVF * DIP_EVF
          ELSE
                  QMM12 = ( 7.2 - 1.0*Q2EVF ) * DIP_EVF
         ENDIF
         RETURN
         END

         SUBROUTINE QKEM2(EPW,EPQ2,QME22,QMM21,QMM22,QMM23)        
         implicit none
         real epw,epq2,qme22,qmm21,qmm22,qmm23
         real xmprot,q2evf,q2evf0,dip_evf,dip_evf0

         DATA XMPROT/0.938/
         Q2EVF=((EPW**2-XMPROT**2)**2+EPQ2*(EPW+XMPROT)**2)/
     *   (4.*EPW*XMPROT)
         Q2EVF0=((EPW**2-XMPROT**2)**2)/(4.*EPW*XMPROT)
         DIP_EVF = (1.+ Q2EVF/0.71)**(-2.)
         DIP_EVF0 = (1.+ Q2EVF0/0.71)**(-2)
CCC
CCC       NOTE THAT QME22 HAS TO BE RENORMALIZED TO REPRODUCE DATA
CCC           REASON IS UNCLEAR !!!!
CCC
         QME22 = 0.99*SQRT(Q2EVF)*DIP_EVF
     *            /SQRT(Q2EVF0)/DIP_EVF0
C
          QMM21 = 0.* SQRT(Q2EVF)*DIP_EVF
          QMM22 = (0.75*SQRT(5.)-1.5*SQRT(5.)*Q2EVF)* 
     *            SQRT(Q2EVF)*DIP_EVF
          IF(Q2EVF.LE.1.2)THEN
          QMM23 = 5.*Q2EVF*SQRT(Q2EVF)*DIP_EVF
          ELSE
          QMM23 = (7.2 - Q2EVF)*SQRT(Q2EVF)*DIP_EVF
          ENDIF
          RETURN
          END
*************************************************************************                                                            
         REAL FUNCTION SIGMAO(EPW,EPQ2,EPEPS,EPCOS,EPPHI,EPIREA,
     *   POL_ELEC,POL_TARG,POL_TARG_THETA,POL_TARG_PHI,sig0
     *   ,sigu,sigt,sigl,sigi,fkt)

         implicit none
         real epw,epq2,epeps,epcos,EPPHI,pol_elec,pol_targ
         real pol_targ_theta,pol_targ_phi,sig0,sigu,sigt,sigl,sigi
         real fkt
         real sindd,cosdd,pol_x,pol_y,pol_z,sigpt,siget
         real xila,sige,theta_gam,enue,theta_elec,e_elec_scatt
         real e_elec_prim
         real c2,c3
         real q2,sinx,cosx,sinx2,cosx2
         real eph1b,eph2b,eph3b,eph4b,eph5b,eph6b
         real xmpot
         real sq2,pi,xmprot,xmpip,xmpi0,xmpi,dum
         real epsin,epsin2,epcos2,pener,pimp,gener,gimp,gamk
         real p1s,p2s,p3s,p4s,p5s,p6s,p7s,p2s2,p3s2,p4s2,p5s2
         real p6s2,p7s2
         real x,x0,f
         real bpipy0,bpipy1,bpipy2,bpipy3
         integer irea
c	The folowing were not declared in the original AO (RCM)
c	I think they are supposed to be complex
         complex xa4m,xb4m,xc4m,xa4p,xb4p,xc4p,xa5m,xb5m,xc5m,xa5p
         complex xb5p,xc5p,xa6m,xb6m,xc6m

c      CHARACTER*3 GON
c      LOGICAL TEST
      COMMON /IP/IT,IB,NF,IBORN,CUT,IP11
c      COMMON /GOCOM/TEST,GON
      real cut
      integer it,ib,nf,iborn,ip11

         COMMON /BACKG/      
     * BGA0P0,BGA1P0,BGA2P0,BGA3P0,  
     * BGA1M0,BGA2M0,BGA3M0,
     * BGB1P0,BGB2P0,BGB3P0, 
     * BGB2M0,BGB3M0,  
     * BGC0P0,BGC1P0,BGC2P0,BGC3P0,
     * BGC1M0,BGC2M0,BGC3M0,
     * BGA0PP,BGA1PP,BGA2PP,BGA3PP,
     * BGA1MP,BGA2MP,BGA3MP,
     * BGB1PP,BGB2PP,BGB3PP,
     * BGB2MP,BGB3MP,
     * BGC0PP,BGC1PP,BGC2PP,BGC3PP,
     * BGC1MP,BGC2MP,BGC3MP
       COMMON/RCONST/      
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710, 
     * WP13_1720, WD13_1520, WD13_1700, 
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600, 
     * WS31_1620, WS31_1900, WP31_1910,   
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710, 
     * LP13_1720, LD13_1520, LD13_1700, 
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600, 
     * LS31_1620, LS31_1900, LP31_1910,  
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710, 
     * PIBP13_1720, PIBD13_1520, PIBD13_1700, 
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600, 
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,  
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975
         COMMON/RAMP/
     * RAS11_1535,             RCS11_1535, 
     * RAS11_1650,             RCS11_1650, 
     * RAP11_1440,             RCP11_1440, 
     * RAP11_1710,             RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,         
     * RAP13_1910, RBP13_1910, RCP13_1910,         
     * RAD13_1520, RBD13_1520, RCD13_1520, 
     * RAD13_1700, RBD13_1700, RCD13_1700, 
     * RAD15_1675, RBD15_1675, RCD15_1675, 
     * RAF15_1680, RBF15_1680, RCF15_1680, 
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,             RCS31_1620, 
     * RAS31_1900,             RCS31_1900, 
     * RAP31_1910,             RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232, 
     * RAP33_1920, RBP33_1920, RCP33_1920, 
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * RAP33_1600, RBP33_1600, RCP33_1600,
     * RAF17_1990, RBF17_1990, RCF17_1990,
     * RAF15_2000, RBF15_2000, RCF15_2000,
     * RAP11_2100, RBP11_2100, RCP11_2100,
     * RAF35_2000, RBF35_2000, RCF35_2000,
     * RAP13_1870, RBP13_1870, RCP13_1870,
     * RAP31_1925, RBP31_1925, RCP31_1925,
     * RAP13_1980, RBP13_1980, RCP13_1980,
     * RAF15_1955, RBF15_1955, RCF15_1955,
     * RAP13_1955, RBP13_1955, RCP13_1955,
     * RAP33_1975, RBP33_1975, RCP33_1975

         real
     * BGA0P0,BGA1P0,BGA2P0,BGA3P0,  
     * BGA1M0,BGA2M0,BGA3M0,
     * BGB1P0,BGB2P0,BGB3P0, 
     * BGB2M0,BGB3M0,  
     * BGC0P0,BGC1P0,BGC2P0,BGC3P0,
     * BGC1M0,BGC2M0,BGC3M0,
     * BGA0PP,BGA1PP,BGA2PP,BGA3PP,
     * BGA1MP,BGA2MP,BGA3MP,
     * BGB1PP,BGB2PP,BGB3PP,
     * BGB2MP,BGB3MP,
     * BGC0PP,BGC1PP,BGC2PP,BGC3PP,
     * BGC1MP,BGC2MP,BGC3MP

       real
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710, 
     * WP13_1720, WD13_1520, WD13_1700, 
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600, 
     * WS31_1620, WS31_1900, WP31_1910,   
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * WP33_1600,WF17_1990,WF15_2000,WP11_2100,WF35_2000,
     * WP13_1870,WP31_1925,WP13_1980,WF15_1955,WP13_1955,WP33_1975,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710, 
     * LP13_1720, LD13_1520, LD13_1700, 
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600, 
     * LS31_1620, LS31_1900, LP31_1910,  
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420,
     * LP33_1600,LF17_1990,LF15_2000,LP11_2100,LF35_2000,
     * LP13_1870,LP31_1925,LP13_1980,LF15_1955,LP13_1955,LP33_1975,
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710, 
     * PIBP13_1720, PIBD13_1520, PIBD13_1700, 
     * PIBD15_1675, PIBF15_1680,
     * PIBG17_2190,
     * PIBG19_2250, PIBH19_2220,
     * PIBI111_2600, 
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,  
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420,
     * PIBP33_1600,PIBF17_1990,PIBF15_2000,PIBP11_2100,PIBF35_2000,
     * PIBP13_1870,PIBP31_1925,PIBP13_1980,PIBF15_1955,PIBP13_1955,
     * PIBP33_1975

       real
     * RAS11_1535,             RCS11_1535, 
     * RAS11_1650,             RCS11_1650, 
     * RAP11_1440,             RCP11_1440, 
     * RAP11_1710,             RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,         
     * RAP13_1910, RBP13_1910, RCP13_1910,         
     * RAD13_1520, RBD13_1520, RCD13_1520, 
     * RAD13_1700, RBD13_1700, RCD13_1700, 
     * RAD15_1675, RBD15_1675, RCD15_1675, 
     * RAF15_1680, RBF15_1680, RCF15_1680, 
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,             RCS31_1620, 
     * RAS31_1900,             RCS31_1900, 
     * RAP31_1910,             RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232, 
     * RAP33_1920, RBP33_1920, RCP33_1920, 
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * RAP33_1600, RBP33_1600, RCP33_1600,
     * RAF17_1990, RBF17_1990, RCF17_1990,
     * RAF15_2000, RBF15_2000, RCF15_2000,
     * RAP11_2100, RBP11_2100, RCP11_2100,
     * RAF35_2000, RBF35_2000, RCF35_2000,
     * RAP13_1870, RBP13_1870, RCP13_1870,
     * RAP31_1925, RBP31_1925, RCP31_1925,
     * RAP13_1980, RBP13_1980, RCP13_1980,
     * RAF15_1955, RBF15_1955, RCF15_1955,
     * RAP13_1955, RBP13_1955, RCP13_1955,
     * RAP33_1975, RBP33_1975, RCP33_1975

    
      COMPLEX                                                           
     * QAS11_1535, QCS11_1535, 
     * QAS11_1650, QCS11_1650, 
     * QAP11_1440, QCP11_1440,
     * QAP11_1710, QCP11_1710,
     * QAP13_1720, qbp13_1720, QCP13_1720,
     * QAD13_1520, QBD13_1520, QCD13_1520, 
     * QAD13_1700, QBD13_1700, QCD13_1700, 
     * QAD15_1675, QBD15_1675, QCD15_1675, 
     * QAF15_1680, QBF15_1680, QCF15_1680,
     * QAG17_2190, QBG17_2190, QCG17_2190, 
     * QAG19_2250, QBG19_2250, QCG19_2250, 
     * QAH19_2220, QBH19_2220, QCH19_2220,
     * QAI111_2600,QBI111_2600,QCI111_2600, 
     * QAS31_1620, QCS31_1620,
     * QAP31_1910, QCP31_1910,
     * QAP33_1232, QBP33_1232, QCP33_1232,
     * QAP33_1920, QBP33_1920, QCP33_1920,
     * QAD33_1700, QBD33_1700, QCD33_1700,
     * QAD35_1930, QBD35_1930, QCD35_1930,
     * QAF35_1905, QBF35_1905, QCF35_1905,
     * QAF37_1950, QBF37_1950, QCF37_1950,
     * QAH311_2420,QBH311_2420,QCH311_2420,
     * QAP33_1600, QBP33_1600, QCP33_1600,
     * QAF17_1990, QBF17_1990, QCF17_1990,
     * QAF15_2000, QBF15_2000, QCF15_2000,
     * QAP11_2100, QBP11_2100, QCP11_2100,
     * QAF35_2000, QBF35_2000, QCF35_2000,
     * QAP13_1870, QBP13_1870, QCP13_1870,
     * QAP31_1925, QBP31_1925, QCP31_1925,
     * QAP13_1980, QBP13_1980, QCP13_1980,
     * QAF15_1955, QBF15_1955, QCF15_1955,
     * QAP13_1955, QBP13_1955, QCP13_1955,
     * QAP33_1975, QBP33_1975, QCP33_1975,
     * qas31_1900, qcs31_1900,
     * XA0P  ,XA1P  ,XA2P  ,XA3P  ,                                     
     * XA1M  ,XA2M  ,XA3M  ,                                            
     * XB1P  ,XB2P  ,XB3P  ,                                            
     * XB2M  ,XB3M  ,                                                   
     * XC0P  ,XC1P  ,XC2P  ,XC3P  ,                                     
     * XC1M  ,XC2M  ,XC3M                                               
      COMPLEX                                                           
     * RS11_1535, RS11_1650, RP11_1440, RP11_1710, 
     * RP13_1720, RD13_1520, RD13_1700, 
     * RD15_1675, RF15_1680,
     * RG17_2190,
     * RG19_2250, RH19_2220,
     * RI111_2600, 
     * RS31_1620, RS31_1900, RP31_1910,  
     * RP33_1232, RP33_1920, RD33_1700,
     * RD35_1930, RF35_1905,
     * RF37_1950,
     * RH311_2420,
     * RP33_1600, RF17_1990, RF15_2000, RP11_2100, RF35_2000,
     * RP13_1870,RP31_1925,RP13_1980,RF15_1955,RP13_1955,RP33_1975
     
      COMPLEX EPH1,EPH2,EPH3,EPH4,EPH5,EPH6                             
      COMPLEX ADDA0,ADDA1,ADDA2,ADDA3,ADDA4,ADDA5
      COMPLEX ADDB1,ADDB2,ADDB3,ADDB4,ADDB5
      COMPLEX ADDC0,ADDC1,ADDC2,ADDC3,ADDC4,ADDC5
      COMPLEX SUBA0,SUBA1,SUBA2,SUBA3,SUBA4,SUBA5
      COMPLEX SUBB1,SUBB2,SUBB3,SUBB4,SUBB5
      COMPLEX SUBC0,SUBC1,SUBC2,SUBC3,SUBC4,SUBC5
      COMPLEX HNP,HNM,HFP,HFM,HN0,HF0                                   
      COMPLEX X1,X2,Y1,Y2,Y3,Y4,Z1,Z2                                   
      COMPLEX R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11
      REAL HB(6)
                                      
         COMPLEX ZERO                 
         INTEGER EPIREA               
         INTEGER IR                   
         COMPLEX EPRES                
         REAL YMPI,W,PIM,GK           
                                      
         DATA ZERO/(0.,0.)/           
      DATA SQ2        /1.41421  /     
      DATA PI         /3.14159  /     
      DATA XMPROT     /.938     /     
      DATA XMPIP      /.1395    /     
      DATA XMPI0      /.1349    /     
c      DATA XMETA      /.5488    /     
        

        sig0=0.

         DUM=ACOS(EPCOS)                                                
         EPSIN=SIN(DUM)                                                 
         EPSIN2=SIN(DUM/2.)                                             
         EPCOS2=COS(DUM/2.)                                             
         IF(EPIREA.EQ.1.OR.EPIREA.EQ.2)XMPI=XMPI0                                      
         IF(EPIREA.EQ.3.OR.EPIREA.EQ.4)XMPI=XMPIP                                      
                                                                                 
           PENER  =(EPW**2-XMPROT**2+XMPI   **2)/(2.*EPW   )            
           PIMP   =SQRT(PENER   **2 - XMPI   **2)                       
           GENER  =(EPW   **2-XMPROT**2-EPQ2   )/(2.*EPW   )            
           GIMP   =SQRT(GENER**2 + EPQ2   )                             
           GAMK   =(EPW**2-XMPROT**2)/(2.*EPW   )                       
           X=EPCOS                                                      
           P1S  =1.   
           P2S  =3.*X 
           P3S  =(15.*X*X-3.)/2.
           P4S  =(35.*X*X*X-15.*X)/2.
           P5S  =(315.*X*X*X*X-210.*X*X+15.)/8. 
           P6S  =(693.*X*X*X*X*X-630.*X*X*X)/8.
           P7S  =(3003.*X*X*X*X*X*X-3465.*X*X*X*X)/16.
           P2S2 =3.
           P3S2 =15.*X 
           P4S2 =(105.*X*X-15.)/2.
           P5S2 =(315.*X*X*X-105.*X)/2.
           P6S2 =(3465.*X*X*X*X-1890.*X*X)/8.                             
           P7S2 =(18018.*X*X*X*X*X-13860.*X*X*X)/16.
           X0=1./SQRT((PIMP   **2)+(.35**2))                            
           F=(1.+(PIMP   **2)/.71)                                      
           BPIPY0   =1./F                                               
           BPIPY1   =X0*PIMP/F                                       
           BPIPY2   =X0*X0*PIMP*PIMP/F                            
           BPIPY3   =X0*X0*X0*PIMP*PIMP*PIMP/F  

           if (gamk .gt. 0.)then
              FKT=PIMP/GAMK
           else
              FKT=0.
           endif


           Q2=EPQ2                                                      
           W =EPW                                                       
           SINX =EPSIN                                                  
           COSX =EPCOS                                                  
           SINX2=EPSIN2                                                 
           COSX2=EPCOS2                                                 
           IREA=EPIREA   
           
c           IBORN=0                                                                                              
           IF(IBORN.EQ.1) THEN
            CALL BORNT(Q2,W,SINX,COSX,SINX2,COSX2,HB,IREA,CUT)               
            EPH1B   =HB(1)
            EPH2B   =HB(2)
            EPH3B   =HB(3)
            EPH4B   =HB(4)
            EPH5B   =HB(5)
            EPH6B   =HB(6)                                           
           ELSE
            EPH1B=0.
            EPH2B=0.
            EPH3B=0.
            EPH4B=0.
            EPH5B=0.
            EPH6B=0.
           END IF
                                                                        
           IR=EPIREA         
           YMPI  =  XMPI                                                
           W     =  EPW                                                 
           PIM   =  PIMP                                                
           GK    =  GAMK                                                
      
      RS11_1535 = EPRES(IR,XMPROT,YMPI,W,LS11_1535,WS11_1535,
     * .35,0,1,1,PIM,GK,.TRUE.)
       RS11_1650  = EPRES(IR,XMPROT,YMPI,W,LS11_1650,WS11_1650,
     * .35,0,1,1,PIM,GK,.FALSE.)
       RP11_1440 = EPRES(IR,XMPROT,YMPI,W,LP11_1440,WP11_1440,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP11_1710 = EPRES(IR,XMPROT,YMPI,W,LP11_1710,WP11_1710,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP13_1720 = EPRES(IR,XMPROT,YMPI,W,LP13_1720,WP13_1720,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RD13_1520 = EPRES(IR,XMPROT,YMPI,W,LD13_1520,WD13_1520,
     * .35,2,1,1,PIM,GK,.FALSE.)
       RD13_1700 = EPRES(IR,XMPROT,YMPI,W,LD13_1700,WD13_1700,
     * .35,2,1,1,PIM,GK,.FALSE.)
       RD15_1675 = EPRES(IR,XMPROT,YMPI,W,LD15_1675,WD15_1675,
     * .35,2,2,1,PIM,GK,.FALSE.)
       RF15_1680 = EPRES(IR,XMPROT,YMPI,W,LF15_1680,WF15_1680,
     * .35,3,2,1,PIM,GK,.FALSE.)
       RG17_2190 = EPRES(IR,XMPROT,YMPI,W,LG17_2190,WG17_2190,
     * .35,4,3,1,PIM,GK,.FALSE.)
       RG19_2250 = EPRES(IR,XMPROT,YMPI,W,LG19_2250,WG19_2250,
     * .35,4,4,1,PIM,GK,.FALSE.)
       RH19_2220 = EPRES(IR,XMPROT,YMPI,W,LH19_2220,WH19_2220,
     * .35,5,4,1,PIM,GK,.FALSE.)
       RI111_2600 = EPRES(IR,XMPROT,YMPI,W,LI111_2600,WI111_2600,
     * .35,6,5,1,PIM,GK,.FALSE.)


       RS31_1620 = EPRES(IR,XMPROT,YMPI,W,LS31_1620,WS31_1620,
     * .35,0,1,3,PIM,GK,.FALSE.)
       RS31_1900 = EPRES(IR,XMPROT,YMPI,W,LS31_1900,WS31_1900,
     * .35,0,1,3,PIM,GK,.FALSE.)
       RP31_1910 = EPRES(IR,XMPROT,YMPI,W,LP31_1910,WP31_1910,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RP33_1232 = EPRES(IR,XMPROT,YMPI,W,LP33_1232,WP33_1232,
     * .185,1,1,3,PIM,GK,.FALSE.)
       RP33_1920 = EPRES(IR,XMPROT,YMPI,W,LP33_1920,WP33_1920,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RD33_1700 = EPRES(IR,XMPROT,YMPI,W,LD33_1700,WD33_1700,
     * .35,2,1,3,PIM,GK,.FALSE.)
       RD35_1930 = EPRES(IR,XMPROT,YMPI,W,LD35_1930,WD35_1930,
     * .35,2,2,3,PIM,GK,.FALSE.)
       RF35_1905 = EPRES(IR,XMPROT,YMPI,W,LF35_1905,WF35_1905,
     * .35,3,2,3,PIM,GK,.FALSE.)
       RF37_1950 = EPRES(IR,XMPROT,YMPI,W,LF37_1950,WF37_1950,
     * .35,3,3,3,PIM,GK,.FALSE.)
       RH311_2420 = EPRES(IR,XMPROT,YMPI,W,LH311_2420,WH311_2420,
     * .35,5,5,3,PIM,GK,.FALSE.)
                                                                         
C* THE FOLOWING ARE ONE OR TWO STAR RESONANCES NOT INCLUDED IN THE ABOVE
       RP33_1600 = EPRES(IR,XMPROT,YMPI,W,LP33_1600,WP33_1600,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RF17_1990 = EPRES(IR,XMPOT,YMPI,W,LF17_1990,WF17_1990,
     * .35,3,3,1,PIM,GK,.FALSE.)
       RF15_2000 = EPRES(IR,XMPOT,YMPI,W,LF15_2000,WF15_2000,
     * .35,3,2,1,PIM,GK,.FALSE.)
       RP11_2100 = EPRES(IR,XMPROT,YMPI,W,LP11_2100,WP11_2100,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RF35_2000 = EPRES(IR,XMPROT,YMPI,W,LF35_2000,WF35_2000,
     * .35,3,2,3,PIM,GK,.FALSE.)

C*  THE FOLOWING ARE RESONANCES PREDICTED BY THE QUARK MODEL ONLY
       RP13_1870 = EPRES(IR,XMPROT,YMPI,W,LP13_1870,WP13_1870,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP31_1925 = EPRES(IR,XMPROT,YMPI,W,LP31_1925,WP31_1925,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RP13_1980 = EPRES(IR,XMPROT,YMPI,W,LP13_1980,WP13_1980,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RF15_1955 = EPRES(IR,XMPROT,YMPI,W,LF15_1955,WF15_1955,
     * .35,3,2,1,PIM,GK,.FALSE.)
       RP13_1955 = EPRES(IR,XMPROT,YMPI,W,LP13_1955,WP13_1955,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP33_1975 = EPRES(IR,XMPROT,YMPI,W,LP33_1975,WP33_1975,
     * .35,1,1,3,PIM,GK,.FALSE.)
                                                     
c           CALL BACK(IB)
                                                                        
           IF(EPIREA   .EQ.1)THEN                                       
             XA0P    = BGA0P0*BPIPY0   +ZERO                            
             XA1P    = BGA1P0*BPIPY1   +ZERO                            
             XA2P    = BGA2P0*BPIPY2   +ZERO                            
	     XA3P    = BGA3P0*BPIPY3   +ZERO
	     xa4p    = zero
	     xa5p    = zero
	     
             XB1P    = BGB1P0*BPIPY1   +ZERO 
             XB2P    = BGB2P0*BPIPY2   +ZERO 
             XB3P    = BGB3P0*BPIPY3   +ZERO
	     xb4p    = zero
	     xb5p    = zero
	     
             XC0P    = BGC0P0*BPIPY0   +ZERO                            
             XC1P    = BGC1P0*BPIPY1   +ZERO                            
             XC2P    = BGC2P0*BPIPY2   +ZERO
             XC3P    = BGC3P0*BPIPY3   +ZERO
             xc4p     =  zero
             xc5p     = zero
                                            
             XA1M    = BGA1M0*BPIPY1   +ZERO
             XA2M    = BGA2M0*BPIPY2   +ZERO
             XA3M    = BGA3M0*BPIPY3   +ZERO
	     xa4m    = zero
	     xa5m     = zero
	     xa6m     = zero
	     	
             XB2M    = BGB2M0*BPIPY2   +ZERO
             XB3M    = BGB3M0*BPIPY3   +ZERO                            
	     xb4m     = zero
	     xb5m    = zero
	     xb6m     = zero
	     
             XC1M    = BGC1M0*BPIPY1   +ZERO                            
             XC2M    = BGC2M0*BPIPY2   +ZERO
             XC3M    = BGC3M0*BPIPY3   +ZERO
             xc4m    =  zero
             xc5m     = zero
             xc6m    = zero
                            
           ELSE 

             XA0P    = BGA0PP*BPIPY0   +ZERO 
             XA1P    = BGA1PP*BPIPY1   +ZERO 
             XA2P    = BGA2PP*BPIPY2   +ZERO 
             XA3P    = BGA3PP*BPIPY3   +ZERO 

             XB1P    = BGB1PP*BPIPY1   +ZERO 
             XB2P    = BGB2PP*BPIPY2   +ZERO 
             XB3P    = BGB3PP*BPIPY3   +ZERO 

             XC0P    = BGC0PP*BPIPY0   +ZERO 
             XC1P    = BGC1PP*BPIPY1   +ZERO 
             XC2P    = BGC2PP*BPIPY2   +ZERO 
             XC3P    = BGC3PP*BPIPY3   +ZERO 
                                    
             XA1M    = BGA1MP*BPIPY1   +ZERO 
             XA2M    = BGA2MP*BPIPY2   +ZERO 
             XA3M    = BGA3MP*BPIPY3   +ZERO 

             XB2M    = BGB2MP*BPIPY2   +ZERO 
             XB3M    = BGB3MP*BPIPY3   +ZERO 

             XC1M    = BGC1MP*BPIPY1   +ZERO 
             XC2M    = BGC2MP*BPIPY2   +ZERO 
             XC3M    = BGC3MP*BPIPY3   +ZERO 

           END IF
 
C* IT=1 FOR EXPT ONLY, IT=2 FOR EXPT+ONE OR TWO STAR FROM QKM
C* IT=2 FOR EXPT+ALL THE OTHER FROM QKM
C* IT=4 FOR ALL FROM QKM
C* C2,C3 ARE USED TO TURN THE QKM PREDICTION OFF AND ON

         IF(IT.EQ.1) THEN
           CALL EXPA(EPIREA,EPQ2)
         IF(IP11.EQ.2) THEN
           RAP11_1440=0.0
         ENDIF
         IF(IP11.EQ.3.OR.IP11.EQ.4) THEN
c           CALL QKMA(EPIREA,EPQ2)
         ENDIF
         C2=0.
         C3=0.
         GO TO 35
         ENDIF

         IF(IT.EQ.2) THEN
         C2=1.
         C3=0.
c         CALL QKMA(EPIREA,EPQ2)
         CALL EXPA(EPIREA,EPQ2)
         GO TO 25
         ENDIF

         IF (IT.EQ.3) THEN
c         CALL QKMA(EPIREA,EPQ2)
         CALL EXPA(EPIREA,EPQ2)
         C2=1.
         C3=1.
         ENDIF

         IF(IT.EQ.4) THEN
         C2=0.
         C3=0.
C* here I get ride of the terms which is not used in
c* the EXPT
c         CALL QKMA(EPIREA,EPQ2)
C*         GO TO 35
         ENDIF

C* THIS FOLLOWING RESONANCES ARE PREDICTED BY QK MODEL ONLY
15       QAP13_1870   =RP13_1870   * RAP13_1870
         QBP13_1870   =RP13_1870   * RBP13_1870
         QCP13_1870   =RP13_1870   * 0.0

         QAP31_1925   =RP31_1925   * RAP31_1925
         QBP31_1925   =RP31_1925   * RBP31_1925
         QCP31_1925   =RP31_1925   * 0.0
 
         QAP13_1980   =RP13_1980   * RAP13_1980
         QBP13_1980   =RP13_1980   * RBP13_1980
         QCP13_1980   =RP13_1980   * 0.0
        
         QAF15_1955   =RF15_1955   * RAF15_1955
         QBF15_1955   =RF15_1955   * RBF15_1955
         QCF15_1955   =RF15_1955   * 0.0
         
         QAP13_1955   =RP13_1955   * RAP13_1955
         QBP13_1955   =RP13_1955   * RBP13_1955
         QCP13_1955   =RP13_1955   * 0.0
     
         QAP33_1975   =RP33_1975   * RAP33_1975
         QBP33_1975   =RP33_1975   * RBP33_1975
         QCP33_1975   =RP33_1975   * 0.0
         

C*  THE FOLLOWING RESONANCES ARE SEEN AS ONE OR TWO STAR
25       QAP33_1600   =RP33_1600   * RAP33_1600
         QBP33_1600   =RP33_1600   * RBP33_1600
         QCP33_1600   =RP33_1600   * 0.0

         QAF17_1990   =RF17_1990   * RAF17_1990
         QBF17_1990   =RF17_1990   * RBF17_1990
         QCF17_1990   =RF17_1990   * 0.0
        
         QAF15_2000   =RF15_2000   * RAF15_2000
         QBF15_2000   =RF15_2000   * RBF15_2000
         QCF15_2000   =RF15_2000   * 0.0
    
         QAP11_2100   =RP11_2100   * RAP11_2100
         QBP11_2100   =RP11_2100   * RBP11_2100
         QCP11_2100   =RP11_2100   * 0.0

         QAF35_2000   =RF35_2000   * RAF35_2000
         QBF35_2000   =RF35_2000   * RBF35_2000
         QCF35_2000   =RF35_2000   * RCF35_2000

         
C* THE FOLLWING RESONANCES ARE SEEN BY EXPERIMENT                                                                        
35       QAS11_1535   =RS11_1535   * RAS11_1535 
         QCS11_1535   =RS11_1535   * RCS11_1535

         QAS11_1650   =RS11_1650   * RAS11_1650
         QCS11_1650   =RS11_1650   * RCS11_1650

         QAP11_1440   =RP11_1440   * RAP11_1440  
         QCP11_1440   =RP11_1440   * RCP11_1440  

         QAP13_1720   =RP13_1720   * RAP13_1720
         QBP13_1720   =RP13_1720   * RBP13_1720
         QCP13_1720   =RP13_1720   * RCP13_1720

         QAD13_1520   =RD13_1520   * RAD13_1520
         QBD13_1520   =RD13_1520   * RBD13_1520
         QCD13_1520   =RD13_1520   * RCD13_1520

         QAD13_1700   =RD13_1700   * RAD13_1700
         QBD13_1700   =RD13_1700   * RBD13_1700
         QCD13_1700   =RD13_1700   * RCD13_1700

         QAD15_1675   =RD15_1675   * RAD15_1675
         QBD15_1675   =RD15_1675   * RBD15_1675
         QCD15_1675   =RD15_1675   * RCD15_1675

         QAF15_1680   =RF15_1680   * RAF15_1680  
         QBF15_1680   =RF15_1680   * RBF15_1680
         QCF15_1680   =RF15_1680   * RCF15_1680

         QAG17_2190   =RG17_2190   * RAG17_2190
         QBG17_2190   =RG17_2190   * RBG17_2190
         QCG17_2190   =RG17_2190   * RCG17_2190

         QAG19_2250   =RG19_2250   * RAG19_2250
         QBG19_2250   =RG19_2250   * RBG19_2250
         QCG19_2250   =RG19_2250   * RCG19_2250

         QAH19_2220   =RH19_2220   * RAH19_2220
         QBH19_2220   =RH19_2220   * RBH19_2220
         QCH19_2220   =RH19_2220   * RCH19_2220

         QAI111_2600  =RI111_2600  * RAI111_2600
         QBI111_2600  =RI111_2600  * RBI111_2600
         QCI111_2600  =RI111_2600  * RCI111_2600

         QAS31_1620   =RS31_1620   * RAS31_1620  
         QCS31_1620   =RS31_1620   * RCS31_1620  

         QAS31_1900   =RS31_1900   * RAS31_1900  
         QCS31_1900   =RS31_1900   * RCS31_1900  

         QAP31_1910   =RP31_1910   * RAP31_1910  
         QCP31_1910   =RP31_1910   * RCP31_1910  

         QAP33_1232   =RP33_1232   * RAP33_1232  
         QBP33_1232   =RP33_1232   * RBP33_1232  
         QCP33_1232   =RP33_1232   * RCP33_1232

         QAP33_1920   =RP33_1920   * RAP33_1920  
         QBP33_1920   =RP33_1920   * RBP33_1920  
         QCP33_1920   =RP33_1920   * RCP33_1920  


         QAD33_1700   =RD33_1700   * RAD33_1700  
         QBD33_1700   =RD33_1700   * RBD33_1700  
         QCD33_1700   =RD33_1700   * RCD33_1700  

         QAD35_1930   =RD35_1930   * RAD35_1930
         QBD35_1930   =RD35_1930   * RBD35_1930
         QCD35_1930   =RD35_1930   * RCD35_1930

         QAF35_1905   =RF35_1905   * RAF35_1905
         QBF35_1905   =RF35_1905   * RBF35_1905
         QCF35_1905   =RF35_1905   * RCF35_1905

         QAF37_1950   =RF37_1950   * RAF37_1950
         QBF37_1950   =RF37_1950   * RBF37_1950
         QCF37_1950   =RF37_1950   * RCF37_1950

         QAH311_2420  =RH311_2420  * RAH311_2420
         QBH311_2420  =RH311_2420  * RBH311_2420
         QCH311_2420  =RH311_2420  * RCH311_2420

c	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c	The following were not initialized in my copy of AO (RCM)
	qap11_1710=zero
    	qcp11_1710=zero
c	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
    				
         XA0P=XA0P + QAS11_1535 + QAS11_1650 + QAS31_1620 + QAS31_1900
         XC0P=XC0P + QCS11_1535 + QCS11_1650 + QCS31_1620 + QCS31_1900
                                                       
         XA1M=XA1M + QAP11_1440 + QAP11_1710 + QAP31_1910
     *       +C2*(QAP11_2100)+C3*(QAP31_1925)
         XC1M=XC1M + QCP11_1440 + QCP11_1710 + QCP31_1910
     *       +C2*(QCP11_2100)+C3*(QCP31_1925)
                                                                   
         XA1P=XA1P + QAP13_1720 + QAP33_1232 + QAP33_1920
     *+C2*(QAP33_1600)+C3*(QAP13_1870+QAP13_1980+QAP13_1955+QAP33_1975)
         XB1P=XB1P + QBP13_1720 + QBP33_1232 + QBP33_1920
     *+C2*(QBP33_1600)+C3*(QBP13_1870+QBP13_1980+QBP13_1955+QBP33_1975)
         XC1P=XC1P + QCP13_1720 + QCP33_1232 + QCP33_1920
     *+C2*(QCP33_1600)+C3*(QCP13_1870+QCP13_1980+QCP13_1955+QCP33_1975)
                                                                        
         XA2M = XA2M + QAD13_1520 + QAD13_1700 + QAD33_1700
         XB2M = XB2M + QBD13_1520 + QBD13_1700 + QBD33_1700
         XC2M = XC2M + QCD13_1520 + QCD13_1700 + QCD33_1700
         
         XA2P = XA2P + QAD15_1675 + QAD35_1930
         XB2P = XB2P + QBD15_1675 + QBD35_1930
         XC2P = XC2P + QCD15_1675 + QCD35_1930

         XA3M = XA3M + QAF15_1680 + QAF35_1905
     *        +C2*(QAF15_2000+QAF35_2000)+C3*(QAF15_1955)    
         XB3M = XB3M + QBF15_1680 + QBF35_1905
     *        +C2*(QBF15_2000+QBF35_2000)+C3*(QBF15_1955)    
         XC3M = XC3M + QCF15_1680 + QCF35_1905
     *        +C2*(QCF15_2000+QCF35_2000)+C3*(QCF15_1955)    
                                      
         XA3P = XA3P + QAF37_1950
     *        +C2*(QAF17_1990)
         XB3P = XB3P + QBF37_1950
     *        +C2*(QBF17_1990)
         XC3P = XC3P + QCF37_1950
     *        +C2*(QCF17_1990)

         XA4M = XA4M + QAG17_2190
         XB4M = XB4M + QBG17_2190
         XC4M = XC4M + QCG17_2190
                                  
         XA4P = XA4P + QAG19_2250
         XB4P = XB4P + QBG19_2250
         XC4P = XC4P + QCG19_2250

         XA5M = XA5M + QAH19_2220
         XB5M = XB5M + QBH19_2220
         XC5M = XC5M + QCH19_2220

         XA5P = XA5P + QAH311_2420
         XB5P = XB5P + QBH311_2420
         XC5P = XC5P + QCH311_2420

         XA6M = XA6M + QAI111_2600
         XB6M = XB6M + QBI111_2600
         XC6M = XC6M + QCI111_2600

         SUBA0=XA0P-XA1M 
         SUBA1=XA1P-XA2M 
         SUBA2=XA2P-XA3M 
         SUBA3=XA3P-XA4M 
         SUBA4=XA4P-XA5M 
         SUBA5=XA5P-XA6M 
                         
         SUBC0=XC0P-XC1M 
         SUBC1=XC1P-XC2M 
         SUBC2=XC2P-XC3M 
         SUBC3=XC3P-XC4M 
         SUBC4=XC4P-XC5M 
         SUBC5=XC5P-XC6M 
                         
         SUBB1=XB1P-XB2M 
         SUBB2=XB2P-XB3M 
         SUBB3=XB3P-XB4M 
         SUBB4=XB4P-XB5M 
         SUBB5=XB5P-XB6M 
                         
         ADDA0=XA0P+XA1M 
         ADDA1=XA1P+XA2M 
         ADDA2=XA2P+XA3M 
         ADDA3=XA3P+XA4M 
         ADDA4=XA4P+XA5M 
         ADDA5=XA5P+XA6M 
                         
         ADDC0=XC0P+XC1M 
         ADDC1=XC1P+XC2M 
         ADDC2=XC2P+XC3M 
         ADDC3=XC3P+XC4M 
         ADDC4=XC4P+XC5M 
         ADDC5=XC5P+XC6M 
                         
         ADDB1=XB1P+XB2M 
         ADDB2=XB2P+XB3M 
         ADDB3=XB3P+XB4M 
         ADDB4=XB4P+XB5M 
         ADDB5=XB5P+XB6M 
c
C   Now do the partial wave expansion for the Walker amplitudes
C
       EPH1   =SUBB1*(       -P2S2   )                                  
     *       + SUBB2*(P2S2   -P3S2   )                                  
     *       + SUBB3*(P3S2   -P4S2   )                                  
     *       + SUBB4*(P4S2   -P5S2   )                                  
     *       + SUBB5*(P5S2   -P6S2   )
       EPH1   = EPH1   * EPSIN   * EPCOS2   /SQ2                           
                                                                        
                                                                        
       EPH2  = SUBA0   *(       -P1S   )                                
     *       + SUBA1   *(P1S    -P2S   )                                
     *       + SUBA2   *(P2S    -P3S   )                                
     *       + SUBA3   *(P3S    -P4S   )
     *       + SUBA4   *(P4S    -P5S   )
     *       + SUBA5   *(P5S    -P6S   )
       EPH2  = EPH2    * EPCOS2   * SQ2

                                                                        
       EPH3   = ADDB1    *(       +P2S2   )                             
     *       +  ADDB2    *(P2S2   +P3S2   )                             
     *       +  ADDB3    *(P3S2   +P4S2   )
     *       +  ADDB4    *(P4S2   +P5S2   )
     *       +  ADDB5    *(P5S2   +P6S2   )
       EPH3   = EPH3   *  EPSIN   * EPSIN2   /SQ2     
                                                                        
                                                                        
       EPH4  = ADDA0     * (       +P1S   )                              
     *       + ADDA1     * (P1S    +P2S   )                             
     *       + ADDA2     * (P2S    +P3S   )                             
     *       + ADDA3     * (P3S    +P4S   )
     *       + ADDA4     * (P4S    +P5S   )
     *       + ADDA5     * (P5S    +P6S   )
       EPH4   = EPH4   * EPSIN2   * SQ2              
                                                                        
                                                                        
       EPH5   =SUBC0     *(       -P1S   )                              
     *       + SUBC1     *(P1S    -P2S   )                              
     *       + SUBC2     *(P2S    -P3S   )                              
     *       + SUBC3     *(P3S    -P4S   )
     *       + SUBC4     *(P4S    -P5S   )
     *       + SUBC5     *(P5S    -P6S   )
       EPH5   = EPH5   * EPCOS2   * SQ2                
                                                                        
                                                                        
       EPH6  = ADDC0     * (       +P1S   )                              
     *       + ADDC1     * (P1S    +P2S   )                             
     *       + ADDC2     * (P2S    +P3S   )                             
     *       + ADDC3     * (P3S    +P4S   )
     *       + ADDC4     * (P4S    +P5S   )
     *       + ADDC5     * (P5S    +P6S   )
        EPH6   =  EPH6   * EPSIN2   * SQ2              
                                                                        
C
C      Now add the Born term contributions. Remember, these are real functions!
C                                                                        
         EPH1   =EPH1   +CMPLX(EPH1B   ,0.)                             
         EPH2   =EPH2   +CMPLX(EPH2B   ,0.)                             
         EPH3   =EPH3   +CMPLX(EPH3B   ,0.)                             
         EPH4   =EPH4   +CMPLX(EPH4B   ,0.)                             
         EPH5   =EPH5   +CMPLX(EPH5B   ,0.)                             
         EPH6   =EPH6   +CMPLX(EPH6B   ,0.)               
   
         
C
C       Now we are ready to calculate all the observables
C
C          First the unpolarized transverse term
C              
                                                                        
           SIGU   =(CABS(EPH1   )**2                                    
     *            +CABS(EPH2   )**2                                     
     *            +CABS(EPH3   )**2                                     
     *            +CABS(EPH4   )**2)/2                                  
C
C            Next the transverse polarized term
C                                                                       
           SIGT   =REAL(EPH2   *CONJG(EPH3   )                          
     *            -     EPH1   *CONJG(EPH4   ))                         
C
C             And the longitudinal term
C                                                                        
           SIGL   =CABS(EPH5   )**2                                     
     *            +CABS(EPH6   )**2                                     
C
C             Finally, the longitudinal-transverse interference term
C
           SIGI   =SQ2*REAL(EPH5   *CONJG(EPH1   -EPH4   )              
     *            +         EPH6   *CONJG(EPH3   +EPH2   ))             
C
C              Sum up all the unpolarized cross section terms
C
                                                                        
             sig0  = SIGU       + EPEPS*SIGL
           SIGMAO   = sig0
     *             +EPEPS       *SIGT       *COS(EPPHI   *PI/90.)       
     *             +SQRT(EPEPS       *(1+EPEPS       )/2.)*SIGI         
     *             *COS(EPPHI   *PI/180.)                               



C
         if(POL_ELEC .eq. 0. .and. POL_TARG .eq. 0.)then
            sigmao=sigmao*FKT

            return
         endif



C
C              Convert Walker amplitudes to Bartl & Majerotto spin-flip and non spin-flip amplitudes
C
           HNP=(EPH4+EPH1)/SQ2                                          
           HNM=(EPH4-EPH1)/SQ2                                          
           HFP=(EPH3-EPH2)/SQ2                                          
           HFM=(EPH3+EPH2)/SQ2                                          
           HN0=EPH5                                                     
           HF0=EPH6                                                     
C
C               Define the appreviations for the various polarized cross section terms
C               Polarized beam, polarized target, polarized beam - polarized target
C
           X1=HF0*CONJG(HNP) + HN0*CONJG(HFP)                           
           X2=HFM*CONJG(HNP) + HNM*CONJG(HFP)                           
           Y1=HNP*CONJG(HFP) + HNM*CONJG(HFM)                           
           Y2=HNM*CONJG(HFM) - HNP*CONJG(HFP)                           
           Y3=HN0*CONJG(HF0)                                            
           Y4=HN0*CONJG(HFM) - HF0*CONJG(HNM)                           
	   Z1=HN0*CONJG(HNP) - HF0*CONJG(HFP)
           Z2=HNM*CONJG(HNP) - HFM*CONJG(HFP)                           
C
C               These are the terms for the recoil polarizations
C
	   R1=CONJG(HN0)*HF0
           R2=CONJG(HNP)*HFP - CONJG(HNM)*HFM
           R3=CONJG(HN0)*HFM + CONJG(HF0)*HNP
           R4=CONJG(HF0)*HNM - CONJG(HN0)*HFM
           R5=CONJG(HN0)*HNP + CONJG(HF0)*HFP
           R6=CONJG(HNP)*HFP + CONJG(HNM)*HFM
           R7=CONJG(HFM)*HNP + CONJG(HFP)*HNM
           R8=CONJG(HFP)*HFM + CONJG(HNP)*HNM
           R9=CONJG(HNP)*HFM - CONJG(HNM)*HFP
           R10=CONJG(HNP)*HNM + CONJG(HFP)*HFM
           R11=CONJG(HN0)*HNM + CONJG(HF0)*HFM
           

           XILA=SQRT(2.*EPEPS*(EPEPS+1))                                
           SIGE=-POL_ELEC*SQRT(2*EPEPS*(1-EPEPS))*SIN(EPPHI*PI/180.)*         
     *     AIMAG(HN0*CONJG(HNM)+HF0*CONJG(HFM))

                         
C           SIGPT=PX*(-XILA*SIN(EPPHI*PI/180)*AIMAG(X1)-EPEPS*            
C     *         SIN(EPPHI*PI/90.)*AIMAG(X2))                             
C     *       -PY*(AIMAG(Y1)+EPEPS*COS(EPPHI*PI/90.)*AIMAG(Y2)           
C     *          +2.*EPEPS*AIMAG(Y3)+XILA*COS(EPPHI*PI/180)*AIMAG(Y4))   
C     *         +PZ*(EPEPS*SIN(EPPHI*PI/90.)*AIMAG(Z2)+                  
C     *              XILA*SIN(EPPHI*PI/180.)*AIMAG(Z1))                  
C       SIGET=-POL_ELEC*(-PX*(SQRT(2.*EPEPS*(1-EPEPS))*COS(EPPHI*PI/180)  
C     *               *REAL(X1)+SQRT(1-EPEPS**2)*REAL(X2))               
C     *       +PY*SQRT(2.*EPEPS*(1-EPEPS))*SIN(EPPHI*PI/180.)*REAL(Y4)   
C     *       +PZ*(SQRT(1-EPEPS**2)*REAL(Z2)+SQRT(2*EPEPS*(1-EPEPS))*    
C     *        COS(EPPHI*PI/180.)*REAL(Z1)))                             
C           SIGMA = SIGMA + SIGE + SIGPT + SIGET
C           SIGMA=SIGMA*FKT                                              
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C         The target is polarized along the beam axis.
c                                       PS = -1 anti-parallel to beam
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
           IF(EPQ2.LT.1.E-06.OR.EPEPS.LT.1.E-06)THEN
                THETA_GAM = 0.
           ELSE
           ENUE = (EPW**2 - XMPROT**2 + EPQ2)/2./XMPROT
           THETA_ELEC = 2.*ATAN(SQRT(EPQ2/2./(EPQ2+ENUE**2)*
     *                  (1.-EPEPS)/EPEPS))
           E_ELEC_SCATT = -ENUE/2 + SQRT((ENUE/2.)**2 + 
     *                   EPQ2/4./SIN(THETA_ELEC/2.)**2)
           E_ELEC_PRIM = E_ELEC_SCATT + ENUE
           THETA_GAM = ATAN(SIN(THETA_ELEC)/
     *                 (E_ELEC_PRIM/E_ELEC_SCATT - COS(THETA_ELEC)))
           THETA_GAM = THETA_GAM*180./PI
           ENDIF

           POL_X = POL_TARG*SINDD(THETA_GAM)*COSDD(EPPHI)
           POL_Y = POL_TARG*SINDD(THETA_GAM)*SINDD(EPPHI)
           POL_Z = POL_TARG*COSDD(THETA_GAM)

           SIGPT = - POL_X*(SINDD(EPPHI)*
     1  SQRT(2.*EPEPS*(1. + EPEPS))*AIMAG(X1) + 
     2  EPEPS*SINDD(2.*EPPHI)*AIMAG(X2)) - 
     3  POL_Y*(AIMAG(Y1) + 
     4  EPEPS*COSDD(2.*EPPHI)*AIMAG(Y2) + 
     5  2.*EPEPS*AIMAG(Y3) + 
     6  SQRT(2.*EPEPS*(1.+EPEPS))*COSDD(EPPHI)*AIMAG(Y4)) + 
     7  POL_Z*(EPEPS*SINDD(2.*EPPHI)*AIMAG(Z2) + 
     8  SQRT(2.*EPEPS*(1. + EPEPS))*SINDD(EPPHI)*AIMAG(Z1))

        SIGET = POL_ELEC*POL_TARG*(SINDD(THETA_GAM)*(COSDD(EPPHI))**2*
     1  (REAL(Y4) + SQRT(2.*EPEPS*(1.-EPEPS))*REAL(X1))+
     2  COSDD(EPPHI)*(SINDD(THETA_GAM)*SQRT(1.-EPEPS**2)*REAL(X2) - 
     3  COSDD(THETA_GAM)*SQRT(2.*EPEPS*(1.-EPEPS))*REAL(Z1)) -
     4  SINDD(THETA_GAM)*REAL(Y4) - COSDD(THETA_GAM)*
     5  SQRT(1.-EPEPS**2)*REAL(Z2))
   
           SIGMAO = SIGMAO + SIGE + SIGPT + SIGET 

 99        SIGMAO = SIGMAO*FKT
           
           RETURN                                                       
           END                                                          

*****************************************************************************

      real function sindd(x)
      implicit none

      real pi,x
      pi=4.*atan(1.)
      sindd=sin(pi*x/180.)
      return
      end
*****************************************************************************
*****************************************************************************


*****************************************************************************
*****************************************************************************

      real function cosdd(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      cosdd=cos(pi*x/180.)
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************

      real function tandd(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      tandd=tan(pi*x/180.)
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************
      
      real function asind(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      asind=180.*asin(x)/pi
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************
      
      real function acosd(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      acosd=180.*acos(x)/pi
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************
      
      real function atand(x)
      implicit none
      real pi,x

      pi=4.*atan(1.)
      atand=180.*atan(x)/pi
      return
      end
*****************************************************************************
