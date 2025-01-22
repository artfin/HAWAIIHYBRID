      SUBROUTINE POTV(V,R1,R2,XCOS2) bind(C, name="potv")
C
C     Written by Thomas Bondo Pedersen, november 2001.
C     - minor modifications jan. 2002 (tbp).
C     - 12-C 16-O and 13-C 16-O vibr. av. included
C       may 2002 (tbp).
C
C     Driver routine for calculating the CO-Ar
C     2D, 3D, or vibrationally averaged potential
C     (including coupling matrix elements if needed).
C     Can also be used for full 3-D calculations, or
C     1-D CO calculations.
C
C     Input:
C
C        R1    - CO distance in bohr
C        R2    - distance from CO center-of-mass to Ar in bohr
C        XCOS2 - cos to the intermolecular angle
C                (angle =   0 degrees corresponds to linear CO---Ar,
C                 angle = 180 degrees corresponds to linear Ar---CO)
C
C     Output:
C
C        V - interaction energy in hartree
C
C     Details of use:
C
C        The form of the potential is based on that used by
C        Toczylowski and Cybulski, J. Chem. Phys. 112, 4604 (2000) 
C
C        Which potential parameters are used depends upon
C        the parameter IPOT transferred in the common block POTDEF
C
C        IPOT = 1 : Use the parameters of Toczylowski and Cybulski's 
C                   CCSD(T)/aug-cc-pVTZ-33221 2D surface for R1 = RE =
C                   2.132 bohr. In this case, R1 is a dummy.
C
C        IPOT = 2 : Use the parameters of the CCSD(T)/aug-cc-pVQZ-33211
C                   2D surface for R1 = RE = 2.132 bohr. In this case,
C                   R1 is a dummy.
C
C        IPOT = 3 : Use the parameters of the CCSD(T)/aug-cc-pVQZ-33211
C                   3D surface at R1. This surface is probably not
C                   trustworthy beyond the interval (this is NOT tested!)
C                      1.898 bohr <= R1 <= 2.234 bohr
C                   due to the single-reference nature of CCSD(T).
C                   Thus, intended only for 2-D calculations at
C                   a given CO distance in the above trust interval.
C
C        IPOT = 4 : Use a matrix element of the potential in the
C                   CO vibrational basis based on a Taylor expansion
C                   through order 8 of the 3D potential surface. The
C                   CO vibrational matrix elements, <IV1|(r - re)^k|IV2>
C                   for k = 0,1,2,...,8, were obtained using MOLCAS and
C                   the CO potential curve of Huxley and Murrel, J. Chem.
C                   Soc. Faraday Trans. II 79, 323 (1983) with
C                   re = 2.132 bohr.
C                   The CO vibrational states are specified through
C                   variables IV1 and IV2.
C                   The isotope numbers of C and O are stored in
C                   common ISODEF.
C
C        IPOT = 5 : Use the potential of 3 plus the CO potential
C                   of Huxley and Maurray, J. Chem. Soc. Faraday Trans.
C                   II 79, 323 (1983) with re = 2.132 bohr. For full
C                   3-D calculations.
C
C        IPOT = 6 : Use the empirical CO potential of Huxley and Murray,
C                   J. Chem. Soc. Faraday Trans. II 79, 323 (1983) with
C                   re = 2.132 bohr. For diatomic calculations.
C
C        IPOT = 7 : As 2, but with R1 = 1.898 bohr.
C
C        IPOT = 8 : As 2, but with R1 = 2.234 bohr.
C
C        IPOT = 9 : As 5, except that the intermolecular part is kept
C                   fixed at re = 2.132 bohr to mimic the usual assumption
C                   that the intermolecular dynamics is affected
C                   by the CO vibrations only through the value of the
C                   CO rotational constant.
C
C        IPOT = 10: Use a matrix element of the potential in the
C                   CO vibrational basis based on a Taylor expansion
C                   through order 8 of the 3D potential surface. The
C                   CO vibrational matrix elements, <IV1|(r - re)^k|IV2>
C                   for k = 0,1,2,...,8, were obtained using MOLCAS and
C                   the CO potential curve on the interval
C                   1.898 bohr <= r <= 2.234 bohr calculated at the
C                   CCSD(T)/aug-cc-pVQZ level with all orbitals included.
C                   re = 2.132 bohr (NOT the CCSD(T) re !!!).
C                   The CO vibrational states are specified through
C                   variables IV1 and IV2.
C                   The isotope numbers of C and O are stored in
C                   common ISODEF.
C
C     ISTAT = 1: write calculated potential point to unit LUSTAT
C                (LUSTAT MUST BE OPEN ON ENTRY!!!!)
C
C     The routine uses coordinate transformations depending on the masses
C     used for C and O. Thus, subroutine COORD_TRF2 needs the common block
C     MASS containing info for this transformation. The structure of the
C     common block is adapted to TRIATOM.
C
C
      use, intrinsic :: iso_c_binding, only : c_double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     Common block(s) defining the potential to be calculated.
C     ========================================================
C

      COMMON / POTDEF / IPOT, IV1, IV2, ISTAT, LUSTAT
      COMMON / ISODEF / ISOC, ISOO

C
C     Parameter arrays.
C     =================
C
      DIMENSION B(6), D(6), G0(6), G1(6), G2(6), G3(6)
C
C     Name of subroutine.
C     ===================
C
      CHARACTER*4 SECNAM
      PARAMETER (SECNAM = 'POTV')
C
C     Some variables.
C     ===============
C
      PARAMETER (RE = 2.132D0)
      PARAMETER (BOHR = 0.52917725D0)
      PARAMETER (AUTOCM = 219474.625D0)
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)
C
C     Start calculation according to IPOT.
C     ====================================

C     XMASS(3) is assumed to contain the masses of
C     Ar, C, and O in that order!
C     --------------------------------------------
      DIMENSION XMASS(3)
      COMMON /MASS/ XMASS
      
      XMASS(1) = 39.962384
      XMASS(2) = 12.000000
      XMASS(3) = 15.994915 
      
      IPOT = 3
C
      IF (IPOT .EQ. 1) THEN
C
C        CCSD(T)/atz-33221 surface (Cybulski's) at RE = 2.132 bohr.
C        ----------------------------------------------------------
C
         R    = R2
         XCOS = XCOS2
         CALL PARINIT(B,D,G0,G1,G2,G3,C60,C62,C71,C73)
         CALL POTCYB(V,R,XCOS,B,D,G0,G1,G2,G3,C60,C62,C71,C73)

      ELSE IF (IPOT .EQ. 2) THEN
C
C        CCSD(T)/aqz-33211 surface at RE = 2.132 bohr.
C        ---------------------------------------------
C
         RCO = RE
         CALL PARINIT2(B,D,G0,G1,G2,G3,C60,C62,C71,C73)
         CALL COORD_TRF2(RCO,R2,XCOS2,R,XCOS)
         CALL POTCYB(V,R,XCOS,B,D,G0,G1,G2,G3,C60,C62,C71,C73)

      ELSE IF (IPOT .EQ. 3) THEN
C
C        CCSD(T)/aqz-33211 surface at R1.
C        --------------------------------
C
         RCO = R1
         CALL PARINIT3(B,D,G0,G1,G2,G3,C60,C62,C71,C73,RCO)
         CALL COORD_TRF2(RCO,R2,XCOS2,R,XCOS)
         CALL POTCYB(V,R,XCOS,B,D,G0,G1,G2,G3,C60,C62,C71,C73)

      ELSE IF (IPOT .EQ. 4) THEN
C
C        CO vibrationally averaged CCSD(T)/aqz-33211 surface.
C        CO vibrational states are empirical.
C        ----------------------------------------------------
C
         RCO = R1
         CALL COORD_TRF2(RCO,R2,XCOS2,R,XCOS)
         CALL POTVIB(V,RCO,R,XCOS,IV1,IV2)

      ELSE IF (IPOT .EQ. 5) THEN
C     
C        CO empirical + 3-D CCSD(T)/aug-cc-pVQZ-33211 surface.
C        -----------------------------------------------------
C
         RCO = R1
         CALL PARINIT3(B,D,G0,G1,G2,G3,C60,C62,C71,C73,RCO)
         CALL COORD_TRF2(RCO,R2,XCOS2,R,XCOS)
         CALL POTCYB(VI,R,XCOS,B,D,G0,G1,G2,G3,C60,C62,C71,C73)
         CALL POTCO(VCO,RCO)
         V = VI + VCO

      ELSE IF (IPOT .EQ. 6) THEN
C     
C        CO empirical.
C        -------------
C
         RCO  = R2
         R    = ZERO
         XCOS = ONE
         CALL POTCO(V,RCO)

      ELSE IF (IPOT .EQ. 7) THEN
C
C        CCSD(T)/aqz-33211 surface at RCO = 1.898 bohr.
C        ---------------------------------------------
C
         RCO = 1.898D0
         CALL PARINIT3(B,D,G0,G1,G2,G3,C60,C62,C71,C73,RCO)
         CALL COORD_TRF2(RCO,R2,XCOS2,R,XCOS)
         CALL POTCYB(V,R,XCOS,B,D,G0,G1,G2,G3,C60,C62,C71,C73)

      ELSE IF (IPOT .EQ. 8) THEN
C
C        CCSD(T)/aqz-33211 surface at RCO = 2.234 bohr.
C        ---------------------------------------------
C
         RCO = 2.234D0
         CALL PARINIT3(B,D,G0,G1,G2,G3,C60,C62,C71,C73,RCO)
         CALL COORD_TRF2(RCO,R2,XCOS2,R,XCOS)
         CALL POTCYB(V,R,XCOS,B,D,G0,G1,G2,G3,C60,C62,C71,C73)

      ELSE IF (IPOT .EQ. 9) THEN
C     
C        CO empirical + 2-D CCSD(T)/aug-cc-pVQZ-33211 surface
C                       at RCO = 2.132 bohr.
C        ----------------------------------------------------
C
         RCO = RE
         CALL PARINIT2(B,D,G0,G1,G2,G3,C60,C62,C71,C73)
         CALL COORD_TRF2(RCO,R2,XCOS2,R,XCOS)
         CALL POTCYB(VI,R,XCOS,B,D,G0,G1,G2,G3,C60,C62,C71,C73)
         CALL POTCO(VCO,RCO)
         V = VI + VCO

      ELSE IF (IPOT .EQ. 10) THEN
C
C        CO vibrationally averaged CCSD(T)/aqz-33211 surface.
C        The CO vibrational states are from CCSD(T)/aqz.
C        ----------------------------------------------------
C
         RCO = R1
         CALL COORD_TRF2(RCO,R2,XCOS2,R,XCOS)
         CALL POTVIB(V,RCO,R,XCOS,IV1,IV2)

      ELSE
C
C        Undefined surface: abort.
C        -------------------------
C
         WRITE(6,'(//,5X,A,A,I6,/,5X,A)')
     &   SECNAM,': Undefined potential, IPOT = ',IPOT,
     &   ' - program will be aborted !!!'
         STOP

      ENDIF

      IF (ISTAT .EQ. 1) THEN
         PI     = DACOS(-1.00D0)
         TODEG  = 180.000D0/PI
         THETA2 = DACOS(XCOS2)*TODEG 
         THETA  = DACOS(XCOS)*TODEG 
         WRITE(LUSTAT,'(/,A,1X,F15.7,1X,F15.7,1X,F15.7)')
     &   'Input  coord.:',R1*BOHR,R2*BOHR,THETA2
         WRITE(LUSTAT,'(A,1X,F15.7,1X,F15.7,1X,F15.7)')
     &   'Actual coord.:',RCO*BOHR,R*BOHR,THETA
         WRITE(LUSTAT,'(A,1X,F15.7)')
     &   'Potential    :',V*AUTOCM
         IF ((IPOT.EQ.5) .OR. (IPOT.EQ.9)) THEN
            WRITE(LUSTAT,'(A,1X,F15.7,1X,F15.7)')
     &      'V(CO), VI    :',VCO*AUTOCM,VI*AUTOCM
         ENDIF
      ENDIF

      RETURN
      END
      SUBROUTINE POTCO(V,R)
C
C     Evaluate the empirical CO potential of Huxley and Murray,
C     J. Chem. Soc. Faraday Trans. II 79, 323 (1983), with
C     RE = 2.132 bohr.
C
C     Input:
C
C        R - in bohr
C
C     Output:
C
C        V - in hartree
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION A(3), RHOP(0:3)
C
C     The next 3 lines define the potential parameters.
C     =================================================
C
      DATA A / 3.897D0, 2.305D0, 1.898D0 /
      PARAMETER (DE = 11.226D0)
      PARAMETER (RE = 2.132D0)

      PARAMETER (BOHR = 0.52917725D0)
      PARAMETER (AUTOEV = 27.2113834D0)
      PARAMETER (ONE = 1.00D0)
C
C     Calculate V in eV.
C     ==================
C
      RHO = (R - RE)*BOHR
      V   = ONE
      CALL CALPOW(RHOP,RHO,3)
      DO I = 1,3
         V = V + A(I)*RHOP(I)
      ENDDO
      XPO = -A(1)*RHOP(1)
      V   = -DE*V*DEXP(XPO)
C
C     Convert V to hartree.
C     =====================
C
      V = V/AUTOEV

      RETURN
      END
      SUBROUTINE POTCYB(V,R2,XCOS,BC,DC,G0,G1,G2,G3,C60,C62,C71,C73)
C
C     Evaluate the potential form from JCP 112, 4604 (2000)
C
C        - Changes relative to this paper:
C          [they are due to B*R being negative here]
C
C           Eq. (2): D - B*R    ---> D + B*R
C           Eq. (6): In sum:  x ---> |x|
C                    In exp: -x ---> +x
C           (Tab. II: b and d parameters should be interchanged.)
C
C     Input:
C
C        R2 in bohr, XCOS (dimensionless)
C        Parameters BC, DC, etc.
C
C     Output:
C
C        V in hartree
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION BC(6), DC(6), G0(6), G1(6), G2(6), G3(6)

      DIMENSION PL(6)

      PARAMETER (BOHR = 0.52917725D0, XMILLI = 1.00D-3)
C
C     Convert R2 to Angstrom.
C     =======================
C
      R = R2*BOHR
C
C     Get the lowest Legendre functions.
C     ----------------------------------
C
      CALL CYB_LEGENDRE(PL,XCOS,6)
C
C     Short-range part.
C     =================
C
      CALL CYB_X(B,PL,BC)
      CALL CYB_X(D,PL,DC)
      CALL CYB_G(G,R,PL,G0,G1,G2,G3)

      BR   = B*R
      VSH  = G*DEXP(D + BR)

C
C     Asymptotic part.
C     ----------------
C
      CALL CYB_TT(F6,BR,6)
      CALL CYB_TT(F7,BR,7)

      RQ  = R*R
      R6  = RQ*RQ*RQ
      R7  = R6*R
      AS6 = (PL(1)*C60 + PL(3)*C62)/R6
      AS7 = (PL(2)*C71 + PL(4)*C73)/R7
      VAS = F6*AS6 + F7*AS7
C
C     Total, convert to hartree (from mEh).
C     -------------------------------------
C
      V = VSH + VAS
      V = V*XMILLI

      RETURN
      END
      SUBROUTINE POTVIB(V,R1,R2,XCOS,IV1,IV2)
C
C     Evaluate the potential matrix element between CO vibrational
C     states IV1 and IV2 for the aQZ-33211 3D potential
C     (defined through IPOT = 3 below).
C
C     Input:
C
C        R1 in bohr (only needed for debugging)
C        R2 in bohr, XCOS (dimensionless)
C        IV1, IV2 = 0,1,2 represent the CO vibrational states.
C
C     Output:
C
C        V in hartree
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (NORDR = 8, IPOT = 3)
      DIMENSION DVSH(0:NORDR), DVAS(0:NORDR)
      DIMENSION FACT(0:NORDR)
      DIMENSION VIBR(0:2,0:2,0:NORDR)

      PARAMETER (BOHR = 0.52917725D0)
      PARAMETER (RE = 2.132D0)
      PARAMETER (ZERO = 0.000D0, ONE  = 1.000D0)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*6 SECNAM
      PARAMETER (SECNAM = 'POTVIB')
C
C     Debug: echo input.
C     ==================
C
      IF (LOCDBG) THEN
         WRITE(6,'(//,A,A)')
     &   SECNAM,': input parameters:'
         WRITE(6,'(A,/,3F15.7)')
     &   'R1, R2, XCOS (bohr, bohr, dimensionless):',
     &   R1,R2,XCOS
         WRITE(6,'(A,/,2I6)')
     &   'Vibrational element calculated for states:',
     &   IV1,IV2
      ENDIF
C
C     Check IV1 and IV2.
C     ==================
C
      IF ((IV1.LT.0) .OR. (IV1.GT.2) .OR. (IV2.LT.0) .OR. (IV2.GT.2))
     & THEN
         WRITE(6,'(//,5X,A,A,/,5X,A,I4,A,I4,/,5X,A)')
     &   SECNAM,': State indices must be 0, 1, or 2.',
     &   'IV1 = ',IV1,' IV2 = ',IV2,
     &   ' - program will be aborted !!!'
         STOP
      ENDIF
C
C     Get vibrational matrix elements of (r - re)^k.
C     ==============================================
C
      CALL VIBMAT(VIBR,NORDR)
      CONV = ONE
      DO K = 1,NORDR
         CONV = CONV*BOHR
         VIBR(IV1,IV2,K) = VIBR(IV1,IV2,K)*CONV
      ENDDO

      IF (LOCDBG) THEN
         WRITE(6,'(A)')
     &   'Matrix elements of (r - re)^k: (k, element in Angstrom^k)'
         DO K = 0,NORDR
            WRITE(6,'(I6,2X,1P,D15.6)') K,VIBR(IV1,IV2,K)
         ENDDO
      ENDIF
C
C     Calculate Taylor expansion of V through order NORDR.
C     ====================================================
C
      CALL VTAYLOR(DVSH,DVAS,R1,R2,XCOS,NORDR,IPOT)
C
C     Calculate matrix element.
C     =========================
C
      CALL CALFACT(FACT,NORDR)
      V = ZERO
      DO K = 0,NORDR
         FCT = ONE/FACT(K)
         V   = V
     &       + DVSH(K)*VIBR(IV1,IV2,K)*FCT
     &       + DVAS(K)*VIBR(IV1,IV2,K)*FCT
      ENDDO

      RETURN
      END
      SUBROUTINE VTAYLOR(DVSH,DVAS,R1,R2,XCOS2,N,IPOT)
C
C     Purpose: Calculate the derivatives
C
C        DVSH(k) = (d^k VSH/dR1^k)_R1=RE
C
C        DVAS(k) = (d^k VAS/dR1^k)_R1=RE
C
C        for k=0,1,2,...,N<=NMAX (NMAX defined as parameter).
C
C        NOTICE: N CANNOT BE LESS THAN 1 !!!!
C
C        The value of RE is defined implicitly by the
C        parameters, its precise value being irrelevant here.
C
C        R1 and R2 must be supplied in units of bohr, although
C        the parameters are expected to be appropriate for these
C        in Angstrom !!!
C        Unless debugging is turned on, R1 is dummy.
C
C        The derivatives are returned in units of hartree/Angstrom^k
C
C        The potential has the Cybulski-form
C        [JCP 112, 4604 (2000)]
C
C        VSH = G(R1,R2,XCOS2) * EXP[D(R1,XCOS2) + B(R1,XCOS2)*R2]
C
C        VAS = F6(B(R1,XCOS2)*R2) * C6(R1,XCOS2)/R2**6
C            + F7(B(R1,XCOS2)*R2) * C7(R1,XCOS2)/R2**7
C
C        The parameters are obtained according to IPOT,
C
C        IPOT = 3 : Get parameters from routine GETPAR3.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION DVSH(0:N), DVAS(0:N)

      PARAMETER (NMAX = 50)

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'VTAYLOR')

      DIMENSION GPAR(0:5,0:3,0:2)
      DIMENSION DPAR(0:5,0:2), BPAR(0:5,0:2)
      DIMENSION C60PAR(0:2), C62PAR(0:2)
      DIMENSION C71PAR(0:2), C73PAR(0:2)

      DIMENSION G(0:2), D(0:2), BR(0:2)
      DIMENSION C6(0:2), C7(0:2)
      DIMENSION XPO(0:NMAX), DPBR(0:2)
      DIMENSION F6(0:NMAX), F7(0:NMAX)
      DIMENSION SUM(0:7,0:NMAX)

      DIMENSION PLEG(0:5), RPOWER(0:7)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0, TWO = 2.00D0)
      PARAMETER (BOHR = 0.52917725D0)
      PARAMETER (AUTOCM = 219474.625D0)
      PARAMETER (XMILLI = 1.00D-3)
C
C     Abort if N < 1.
C     ===============
C
      IF (N .LT. 1) THEN
         WRITE(6,'(//,5X,A,A,/,5X,A)')
     &   SECNAM,' cannot handle orders less than 1 !!',
     &   ' - aborting execution'
         STOP ' !!! TAYLOR EXPANSION TOO SHORT !!! '
      ELSE IF (N .GT. NMAX) THEN
         WRITE(6,'(//,5X,A,A,I6,A,/,5X,A)')
     &   SECNAM,' cannot handle orders larger than ',NMAX,' !!',
     &   ' - aborting execution'
         STOP ' !!! TAYLOR EXPANSION TOO LONG !!! '
      ENDIF
C
C     Get the parameters according to IPOT.
C     =====================================
C
      IF (IPOT .EQ. 3) THEN
         CALL GETPAR3(BPAR,DPAR,GPAR,C60PAR,C62PAR,C71PAR,C73PAR)
      ELSE
         WRITE(6,'(//,5X,A,I9,/,5X,A,A)')
     &   'Parameters not defined by IPOT = ',IPOT,
     &   ' - aborting execution in subroutine ',SECNAM
         STOP ' !!! UNDEFINED PARAMETERS !!! '
      ENDIF
C
C     Convert R2 to Angstrom.
C     =======================
C
      R   = R2*BOHR
C
C     Get Legendre polynomials through order 5.
C     =========================================
C
      CALL CYB_LEGENDRE(PLEG,XCOS2,6)
C
C     Set up powers of R through order 7.
C     ===================================
C
      CALL CALPOW(RPOWER,R,7)
C
C     Calculate 0th, 1st, and 2nd order parameter functions:
C     par = par(0) + par(1)*(r - re) + par(2)*(r - re)^2.
C     ======================================================
C
      DO K = 0,2
         BR(K) = ZERO
         D(K)  = ZERO
         G(K)  = ZERO
         DO L = 0,5
            BR(K) = BR(K) + BPAR(L,K)*PLEG(L)
            D(K)  = D(K)  + DPAR(L,K)*PLEG(L)
         ENDDO
         DO M = 0,3
            DO L = 0,5
               G(K) = G(K) + GPAR(L,M,K)*RPOWER(M)*PLEG(L)
            ENDDO
         ENDDO
         BR(K)   = BR(K)*R
         D(K)    = D(K)
         G(K)    = G(K)
         DPBR(K) = D(K) + BR(K)
         C6(K)   = (C60PAR(K)*PLEG(0) + C62PAR(K)*PLEG(2))/RPOWER(6)
         C7(K)   = (C71PAR(K)*PLEG(1) + C73PAR(K)*PLEG(3))/RPOWER(7)
      ENDDO
C
C     Debug: print parameters.
C     ========================
C
      IF (LOCDBG) THEN
         WRITE(6,'(A,A,/,A,A)')
     &   'Order:        0               1               2        ',
     &   '   Total value',
     &   '*******************************************************',
     &   '***************'
         RHO1 = (R1 - 2.132D0)*BOHR
         RHO2 = RHO1*RHO1
         VAL1 = BR(0) + BR(1)*RHO1 + BR(2)*RHO2
         VAL2 = D(0)  + D(1)*RHO1  + D(2)*RHO2
         VAL3 = G(0)  + G(1)*RHO1  + G(2)*RHO2
         VAL4 = C6(0) + C6(1)*RHO1 + C6(2)*RHO2
         VAL5 = C7(0) + C7(1)*RHO1 + C7(2)*RHO2
         WRITE(6,1111) 'B*R  :',(BR(K),K=0,2),VAL1
         WRITE(6,1111) 'D    :',(D(K),K=0,2),VAL2
         WRITE(6,1111) 'G    :',(G(K),K=0,2),VAL3
         WRITE(6,1111) 'C6   :',(C6(K),K=0,2),VAL4
         WRITE(6,1111) 'C7   :',(C7(K),K=0,2),VAL5
         WRITE(6,'(A,A)')
     &   '*******************************************************',
     &   '***************'
 1111    FORMAT(A6,1X,F15.7,1X,F15.7,1X,F15.7,1X,F15.7)
         XVSH = VAL3*DEXP(VAL2+VAL1)
         CALL CYB_TT(TT6,VAL1,6)
         CALL CYB_TT(TT7,VAL1,7)
         XVAS = TT6*VAL4 + TT7*VAL5
         XV   = XVSH + XVAS
         XSC  = XMILLI*AUTOCM
         WRITE(6,'(A)') 'RCO, R, THETA (A, A, Deg.):'
         WRITE(6,'(F15.7,1X,F15.7,1X,F15.7)')
     &   R1*BOHR,R2*BOHR,180.00D0*DACOS(XCOS2)/DACOS(-1.00D0)
         WRITE(6,'(A)') 'Exact potential values:  VSH, VAS, V (cm-1)'
         WRITE(6,'(F15.7,1X,F15.7,1X,F15.7)')
     &   XVSH*XSC,XVAS*XSC,XV*XSC
      ENDIF
C
C     Calculate derivatives of exp(D+B*R) through order N.
C     ====================================================
C
      CALL CALXPO(XPO,DPBR,N)
C
C     Calculate DVSH through order N.
C     ===============================
C
      DO K = 0,N
         DVSH(K) = ZERO
      ENDDO
      CALL FINCON(DVSH,G,XPO,N)
C
C     Calculate F6 and F7 through order N.
C     ====================================
C
      CALL CYB_DTT(F6,F7,BR,SUM,N)
C
C     Calculate DVAS through order N.
C     ===============================
C
      DO K = 0,N
         DVAS(K) = ZERO
      ENDDO
      CALL FINCON(DVAS,C6,F6,N)
      CALL FINCON(DVAS,C7,F7,N)
C
C     Finally convert to Eh from mEh.
C     ===============================
C
      DO K = 0,N
         DVSH(K) = DVSH(K)*XMILLI
         DVAS(K) = DVAS(K)*XMILLI
      ENDDO
C
C     Debug: Print energies.
C     ======================
C
      IF (LOCDBG) THEN
         TSH = DVSH(0)
         TAS = DVAS(0)
         FAC = ONE
         RHO = (R1 - 2.132D0)*BOHR
         RHP = ONE
         DO K = 1,N
            FAC = FAC*DFLOAT(K)
            RHP = RHP*RHO
            SCL = RHP/FAC
            TSH = TSH + DVSH(K)*SCL
            TAS = TAS + DVAS(K)*SCL
         ENDDO
         TT = TSH + TAS
         WRITE(6,'(A)') 'Taylor potential values:  VSH, VAS, V (cm-1)'
         WRITE(6,'(F15.7,1X,F15.7,1X,F15.7,//)')
     &   TSH*AUTOCM,TAS*AUTOCM,TT*AUTOCM
      ENDIF

      RETURN
      END
      SUBROUTINE CYB_X(X,PL,XC)
C
C     Evaluate the sum:
C
C        X = sum_(i=1,6) XC(i)*PL(i)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION PL(6), XC(6)
C
      PARAMETER (ZERO = 0.00D0)
C
      X = ZERO
      DO I = 1,6
         X = X + XC(I)*PL(I)
      ENDDO
C
      RETURN
      END
      SUBROUTINE CYB_G(G,R,PL,G0,G1,G2,G3)
C
C     Evaluate the sum:
C
C        G = sum_(i=1,6) [ G0(i)     + G1(i)*R
C                        + G2(i)*R*R + G3(i)*R*R*R]
C                      * PL(i)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION PL(6), G0(6), G1(6), G2(6), G3(6)
C
      PARAMETER (ZERO = 0.00D0)
C
      R2 = R*R
      R3 = R2*R
      G  = ZERO
      DO I = 1,6
         FAC = G0(I) + G1(I)*R + G2(I)*R2 + G3(I)*R3
         G   = G + FAC*PL(I)
      ENDDO
C
      RETURN
      END
      SUBROUTINE CYB_TT(FN,X,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     Evaluate the Tang-Toennis function:
C
C        FN = 1 - exp(X) sum_(k=1,N) |X|^k/k!
C
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)
C
      XSAV = X
      X    = DABS(XSAV)
C
      DCOUN = ZERO
      FACN  = ONE
      XN    = ONE
      SUM   = ONE
      DO I = 1,N
         DCOUN = DCOUN + ONE
         FACN  = FACN*DCOUN
         XN    = XN*X
         SUM   = SUM + XN/FACN
      ENDDO
C
      X  = XSAV
      FN = ONE - DEXP(X)*SUM
C
      RETURN
      END
      SUBROUTINE CYB_LEGENDRE(PL,X,N)
C
C     Calculate (recursively) the Legendre polymials from
C     order 0 to N-1 for a given value of their argument
C     -1 <= X <= 1 (THIS IS NOT TESTED !).
C
C     Store in PL(i=1,N), N >= 1.
C
C     Note: the Legendre polynomials are NOT normalized.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION PL(*)
C
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0, TWO = 2.00D0)
C
C     Check N.
C     --------
C
      IF (N .LE. 0) RETURN
C
C     0th order = 1 (constant)
C     ------------------------
C
      PL(1) = ONE
      IF (N .EQ. 1) GO TO 101
C
C     1st order = X
C     -------------
C
      PL(2) = X
C
C     Higher orders recursively
C     -------------------------
C
      DN = ONE
      DO I = 3,N
         DN = DN + ONE
         PL(I) = (TWO*DN - ONE)*X*PL(I-1)
     &         - (DN - ONE)*PL(I-2)
         PL(I) = PL(I)/DN
      ENDDO
C
  101 RETURN
      END
      SUBROUTINE COORD_TRF2(R1,R2,XCOS,R2P,XCOSP)
C
C     PURPOSE:
C     ========
C              Transform the internal coordinates
C              R1,R2,XCOS to coordinates appropriate
C              for the potential energy evaluation,
C              returned in R2P,XCOSP (R1 is unchanged).
C
C     NOTES:
C     ======
C           1) The routine gets the triatom masses from
C              COMMON /MASS/
C           2) It is assumed that the coordinate system
C              is defined such that Theta = 0 Degr.
C              corresponds to the linear Ar---O=C.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     XMASS(3) is assumed to contain the masses of
C     Ar, C, and O in that order!
C     --------------------------------------------
C
      COMMON /MASS/ XMASS(3),G1,G2
C
C     These are the masses used in the CO-Ar calculations
C     for the electronic interaction energy, and hence
C     the ones that POTV refers to.
C     ---------------------------------------------------
C
      PARAMETER(XMO = 15.9994D0, XMC = 12.0112D0)
C
C     Some parameters.
C     ----------------
C
      PARAMETER (XMONE = -1.00D0, ZERO = 0.00D0, HALF = 0.50D0)
      PARAMETER (TOL = 1.00D-15)
C
C     Calculate pi/2.
C     ---------------
C
      PIHF = HALF*DACOS(XMONE)
C
C     Calculate Cartesian coordinates of Ar.
C     --------------------------------------
C
      THETA = DACOS(XCOS)
      GAMMA = PIHF - THETA
C
      XAR = R2*DCOS(GAMMA)
      YAR = R2*DSIN(GAMMA)
C
C     Calculate center of mass of CO using XMC and XMO.
C     -------------------------------------------------
C
      DNOM = (XMC + XMO)*(XMASS(2) + XMASS(3))
C
      YCM = R1*(XMASS(2)*XMO - XMC*XMASS(3))/DNOM
C
C     Translate origin to YCM.
C     ------------------------
C
      YARP = YAR - YCM
C
C     Calculate new coordinates.
C     --------------------------
C
      GAMMA = DATAN2(YARP,XAR)
      THETA = PIHF - GAMMA
C
      R2P   = DSQRT(XAR*XAR + YARP*YARP)
      XCOSP = DCOS(THETA)
C
      RETURN
      END
      SUBROUTINE CYB_DTT(F6,F7,BR,SUM,N)
C
C     Purpose: Calculate the derivatives of the Tang-Toennis
C              damping functions.
C
C     NOTE: N MUST BE AT LEAST 1. (NOT TESTED!)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F6(0:N), F7(0:N)
      DIMENSION BR(0:2)
      DIMENSION SUM(0:7,0:N)
      DIMENSION BRPOW(0:7), FACT(0:7)
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0, TWO = 2.00D0)
C
C     Calculate powers of BR(0) and the factorial.
C     ============================================
C
      CALL CALPOW(BRPOW,BR(0),7)
      CALL CALFACT(FACT,7)
C
C     Calculate 0th order sums.
C     =========================
C
      SUM(0,0) = DEXP(BR(0))
      DO J = 1,7
         SUM(J,0) = BRPOW(J)*SUM(0,0)/FACT(J)
      ENDDO
      DO J = 1,7,2
         SUM(J,0) = -SUM(J,0)
      ENDDO
C
C     Calculate 1st order sums.
C     =========================
C
      SUM(0,1) = BR(1)*SUM(0,0)
      DO J = 1,7
         SUM(J,1) = BR(1)*(SUM(J,0) - SUM(J-1,0))
      ENDDO
C
C     Calculate sums through order N.
C     ===============================
C
      DO K = 2,N
         SCALBR2  = TWO*DFLOAT(K-1)*BR(2)
         SUM(0,K) = SCALBR2*SUM(0,K-2) + BR(1)*SUM(0,K-1)
         DO J = 1,7
            SUM(J,K) = SCALBR2*(SUM(J,K-2) - SUM(J-1,K-2))
     &               +   BR(1)*(SUM(J,K-1) - SUM(J-1,K-1))
         ENDDO
      ENDDO
C
C     Calculate the Tang-Toennis derivatives.
C     =======================================
C
      F6(0) = ONE
      DO K = 1,N
         F6(K) = ZERO
      ENDDO
      DO K = 0,N
         DO J = 0,6
            F6(K) = F6(K) - SUM(J,K)
         ENDDO
         F7(K) = F6(K) - SUM(7,K)
      ENDDO

      RETURN
      END
      SUBROUTINE FINCON(DV,F1,F2,N)
C
C     Purpose: accumulate the final contributions to the potential derivatives:
C
C        DV(k) = DV(k)
C              + F1(0) * F2(k)
C              + k * F1(1) * F2(k-1)
C              + k*(k-1) * F1(2) * F2(k-2)
C
C     for k = 0,1,2,...,N. N must be at least 1 (not tested).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DV(0:N), F1(0:2), F2(0:N)
      DIMENSION SCAL(0:2)
      PARAMETER (ONE = 1.00D0)

      DV(0) = DV(0) + F1(0)*F2(0)
      DV(1) = DV(1) + F1(0)*F2(1) + F1(1)*F2(0)
      SCAL(0) = ONE
      DO K = 2,N
         SCAL(1) = DFLOAT(K)
         SCAL(2) = DFLOAT(K*(K-1))
         DO M = 0,2
            DV(K) = DV(K) + SCAL(M)*F1(M)*F2(K-M)
         ENDDO
      ENDDO

      RETURN
      END
      SUBROUTINE CALXPO(RES,FN,N)
C
C     Purpose: Calculate derivatives of the exponential function,
C              given the derivatives of the exponent through order
C              2,
C
C              F = FN(0) + FN(1)*RHO + FN(2)*RHO^2
C
C              where RHO = (r - re).
C
C     NOTE: N MUST BE AT LEAST 1 (NOT TESTED!!!)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RES(0:N), FN(0:2)
      PARAMETER (TWO = 2.00D0)

      RES(0) = DEXP(FN(0))
      RES(1) = FN(1)*RES(0)
      DO K = 2,N
         SCAL   = TWO*DFLOAT(K-1)
         RES(K) = SCAL*FN(2)*RES(K-2) + FN(1)*RES(K-1)
      ENDDO

      RETURN
      END
      SUBROUTINE CALPOW(XTOK,X,N)
C
C     Purpose: Calculate powers of X
C
C        XTOK(k) = X**k
C
C        for k=0,1,2,...,N
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XTOK(0:N)
      PARAMETER (ONE = 1.00D0)

      IF (N .LT. 0) RETURN

      XTOK(0) = ONE
      DO K = 1,N
         XTOK(K) = XTOK(K-1)*X
      ENDDO

      RETURN
      END
      SUBROUTINE CALFACT(FACT,N)
C
C     Purpose: Calculate the factorial.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FACT(0:N)
      PARAMETER (ONE = 1.00D0)

      IF (N .LT. 0) RETURN

      FACT(0) = ONE
      DO K = 1,N
         DK      = DFLOAT(K)
         FACT(K) = FACT(K-1)*DK
      ENDDO

      RETURN
      END
      SUBROUTINE VIBMAT(XMAT,N)
C
C     CO vibrational matrix elements from MOLCAS.
C
C        <v1(CO)| (r - re)^k |v2(CO)>
C
C     for v1(CO), v2(CO) = 0,1,2
C     and k              = 0,1,2,3,...,N
C
C     Unit: bohr^k
C
C     NOTE: N can be at most NMAX (defined as parameter).
C
C     The isotope numbers of C and O are taken from common ISODEF.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XMAT(0:2,0:2,0:N)
      PARAMETER (NMAX = 8)
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)

      LOGICAL C12O16, C13O16

      COMMON / POTDEF / IPOT, IV1, IV2, ISTAT, LUSTAT
      COMMON / ISODEF / ISOC, ISOO
C
C     Check N.
C     ========
C
      IF ((N.LT.0) .OR. (N.GT.NMAX)) GO TO 102
C
C     Check isotopes.
C     ===============
C
      IF ((ISOC.EQ.12) .AND. (ISOO.EQ.16)) THEN
         C12O16 = .TRUE.
         C13O16 = .FALSE.
      ELSE IF ((ISOC.EQ.13) .AND. (ISOO.EQ.16)) THEN
         C12O16 = .FALSE.
         C13O16 = .TRUE.
      ELSE
         C12O16 = .FALSE.
         C13O16 = .FALSE.
      ENDIF
C
C     (r - re)^0 = 1
C     ==============
C
      XMAT(0,0,0) =  ONE
      XMAT(1,0,0) =  ZERO
      XMAT(1,1,0) =  ONE
      XMAT(2,0,0) =  ZERO
      XMAT(2,1,0) =  ZERO
      XMAT(2,2,0) =  ONE
      XMAT(0,1,0) =  XMAT(1,0,0)
      XMAT(0,2,0) =  XMAT(2,0,0)
      XMAT(1,2,0) =  XMAT(2,1,0)

      IF (N .EQ. 0) GO TO 101
C
C     (r - re)^1
C     ==========
C
      IF (IPOT .EQ. 4) THEN
         IF (C12O16) THEN
            XMAT(0,0,1) =  0.007708D0
            XMAT(1,0,1) = -0.063792D0
            XMAT(1,1,1) =  0.023255D0
            XMAT(2,0,1) = -0.003643D0
            XMAT(2,1,1) = -0.090479D0
            XMAT(2,2,1) =  0.038996D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,1) =  0.007535D0
            XMAT(1,0,1) = -0.063073D0
            XMAT(1,1,1) =  0.022731D0
            XMAT(2,0,1) = -0.003561D0
            XMAT(2,1,1) = -0.089453D0
            XMAT(2,2,1) =  0.038113D0
         ELSE
            GO TO 104
         ENDIF
      ELSE IF (IPOT .EQ. 10) THEN
         IF (C12O16) THEN
            XMAT(0,0,1) =  0.012613D0
            XMAT(1,0,1) = -0.065180D0
            XMAT(1,1,1) =  0.047140D0
            XMAT(2,0,1) = -0.008267D0
            XMAT(2,1,1) = -0.096926D0
            XMAT(2,2,1) =  0.095363D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,1) =  0.012289D0
            XMAT(1,0,1) = -0.064396D0
            XMAT(1,1,1) =  0.045668D0
            XMAT(2,0,1) = -0.008000D0
            XMAT(2,1,1) = -0.095663D0
            XMAT(2,2,1) =  0.092819D0
         ELSE
            GO TO 104
         ENDIF
      ELSE
         GO TO 103
      ENDIF

      XMAT(0,1,1) =  XMAT(1,0,1)
      XMAT(0,2,1) =  XMAT(2,0,1)
      XMAT(1,2,1) =  XMAT(2,1,1)

      IF (N .EQ. 1) GO TO 101

C
C     (r - re)^2
C     ==========
C
      IF (IPOT .EQ. 4) THEN
         IF (C12O16) THEN
            XMAT(0,0,2) =  0.004142D0
            XMAT(1,0,2) = -0.001643D0
            XMAT(1,1,2) =  0.012837D0
            XMAT(2,0,2) =  0.005641D0
            XMAT(2,1,2) = -0.004689D0
            XMAT(2,2,2) =  0.022154D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,2) =  0.004048D0
            XMAT(1,0,2) = -0.001588D0
            XMAT(1,1,2) =  0.012536D0
            XMAT(2,0,2) =  0.005517D0
            XMAT(2,1,2) = -0.004531D0
            XMAT(2,2,2) =  0.021617D0
         ELSE
            GO TO 104
         ENDIF
      ELSE IF (IPOT .EQ. 10) THEN
         IF (C12O16) THEN
            XMAT(0,0,2) =  0.004477D0
            XMAT(1,0,2) = -0.003077D0
            XMAT(1,1,2) =  0.016188D0
            XMAT(2,0,2) =  0.005541D0
            XMAT(2,1,2) = -0.011029D0
            XMAT(2,2,2) =  0.034173D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,2) =  0.004363D0
            XMAT(1,0,2) = -0.002952D0
            XMAT(1,1,2) =  0.015691D0
            XMAT(2,0,2) =  0.005425D0
            XMAT(2,1,2) = -0.010572D0
            XMAT(2,2,2) =  0.033034D0
         ELSE
            GO TO 104
         ENDIF
      ELSE
         GO TO 103
      ENDIF

      XMAT(0,1,2) =  XMAT(1,0,2)
      XMAT(0,2,2) =  XMAT(2,0,2)
      XMAT(1,2,2) =  XMAT(2,1,2)

      IF (N .EQ. 2) GO TO 101

C
C     (r - re)^3
C     ==========
C
      IF (IPOT .EQ. 4) THEN
         IF (C12O16) THEN
            XMAT(0,0,3) =  0.000116D0
            XMAT(1,0,3) = -0.000818D0
            XMAT(1,1,3) =  0.000765D0
            XMAT(2,0,3) =  0.000264D0
            XMAT(2,1,3) = -0.002428D0
            XMAT(2,2,3) =  0.002109D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,3) =  0.000111D0
            XMAT(1,0,3) = -0.000790D0
            XMAT(1,1,3) =  0.000731D0
            XMAT(2,0,3) =  0.000253D0
            XMAT(2,1,3) = -0.002342D0
            XMAT(2,2,3) =  0.002013D0
         ELSE
            GO TO 104
         ENDIF
      ELSE IF (IPOT .EQ. 10) THEN
         IF (C12O16) THEN
            XMAT(0,0,3) =  0.000209D0
            XMAT(1,0,3) = -0.001011D0
            XMAT(1,1,3) =  0.001853D0
            XMAT(2,0,3) =  0.000530D0
            XMAT(2,1,3) = -0.003838D0
            XMAT(2,2,3) =  0.006797D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,3) =  0.000199D0
            XMAT(1,0,3) = -0.000970D0
            XMAT(1,1,3) =  0.001747D0
            XMAT(2,0,3) =  0.000504D0
            XMAT(2,1,3) = -0.003659D0
            XMAT(2,2,3) =  0.006416D0
         ELSE
            GO TO 104
         ENDIF
      ELSE
         GO TO 103
      ENDIF

      XMAT(0,1,3) =  XMAT(1,0,3)
      XMAT(0,2,3) =  XMAT(2,0,3)
      XMAT(1,2,3) =  XMAT(2,1,3)

      IF (N .EQ. 3) GO TO 101

C
C     (r - re)^4
C     ==========
C
      IF (IPOT .EQ. 4) THEN
         IF (C12O16) THEN
            XMAT(0,0,4) =  0.000052D0
            XMAT(1,0,4) = -0.000046D0
            XMAT(1,1,4) =  0.000286D0
            XMAT(2,0,4) =  0.000151D0
            XMAT(2,1,4) = -0.000235D0
            XMAT(2,2,4) =  0.000810D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,4) =  0.000050D0
            XMAT(1,0,4) = -0.000044D0
            XMAT(1,1,4) =  0.000272D0
            XMAT(2,0,4) =  0.000144D0
            XMAT(2,1,4) = -0.000222D0
            XMAT(2,2,4) =  0.000770D0
         ELSE
            GO TO 104
         ENDIF
      ELSE IF (IPOT .EQ. 10) THEN
         IF (C12O16) THEN
            XMAT(0,0,4) =  0.000065D0
            XMAT(1,0,4) = -0.000104D0
            XMAT(1,1,4) =  0.000498D0
            XMAT(2,0,4) =  0.000210D0
            XMAT(2,1,4) = -0.000726D0
            XMAT(2,2,4) =  0.002069D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,4) =  0.000061D0
            XMAT(1,0,4) = -0.000097D0
            XMAT(1,1,4) =  0.000467D0
            XMAT(2,0,4) =  0.000199D0
            XMAT(2,1,4) = -0.000676D0
            XMAT(2,2,4) =  0.001931D0
         ELSE
            GO TO 104
         ENDIF
      ELSE
         GO TO 103
      ENDIF
      XMAT(0,1,4) =  XMAT(1,0,4)
      XMAT(0,2,4) =  XMAT(2,0,4)
      XMAT(1,2,4) =  XMAT(2,1,4)

      IF (N .EQ. 4) GO TO 101

C
C     (r - re)^5
C     ==========
C
      IF (IPOT .EQ. 4) THEN
         IF (C12O16) THEN
            XMAT(0,0,5) =  0.000003D0
            XMAT(1,0,5) = -0.000018D0
            XMAT(1,1,5) =  0.000028D0
            XMAT(2,0,5) =  0.000013D0
            XMAT(2,1,5) = -0.000085D0
            XMAT(2,2,5) =  0.000116D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,5) =  0.000003D0
            XMAT(1,0,5) = -0.000017D0
            XMAT(1,1,5) =  0.000026D0
            XMAT(2,0,5) =  0.000013D0
            XMAT(2,1,5) = -0.000080D0
            XMAT(2,2,5) =  0.000108D0
         ELSE
            GO TO 104
         ENDIF
      ELSE IF (IPOT .EQ. 10) THEN
         IF (C12O16) THEN
            XMAT(0,0,5) =  0.000006D0
            XMAT(1,0,5) = -0.000029D0
            XMAT(1,1,5) =  0.000087D0
            XMAT(2,0,5) =  0.000036D0
            XMAT(2,1,5) = -0.000211D0
            XMAT(2,2,5) =  0.000535D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,5) =  0.000005D0
            XMAT(1,0,5) = -0.000027D0
            XMAT(1,1,5) =  0.000080D0
            XMAT(2,0,5) =  0.000033D0
            XMAT(2,1,5) = -0.000194D0
            XMAT(2,2,5) =  0.000490D0
         ELSE
            GO TO 104
         ENDIF
      ELSE
         GO TO 103
      ENDIF

      XMAT(0,1,5) =  XMAT(1,0,5)
      XMAT(0,2,5) =  XMAT(2,0,5)
      XMAT(1,2,5) =  XMAT(2,1,5)

      IF (N .EQ. 5) GO TO 101

C
C     (r - re)^6
C     ==========
C
      IF (IPOT .EQ. 4) THEN
         IF (C12O16) THEN
            XMAT(0,0,6) =  0.000001D0
            XMAT(1,0,6) = -0.000002D0
            XMAT(1,1,6) =  0.000009D0
            XMAT(2,0,6) =  0.000005D0
            XMAT(2,1,6) = -0.000012D0
            XMAT(2,2,6) =  0.000038D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,6) =  0.000001D0
            XMAT(1,0,6) = -0.000002D0
            XMAT(1,1,6) =  0.000008D0
            XMAT(2,0,6) =  0.000005D0
            XMAT(2,1,6) = -0.000011D0
            XMAT(2,2,6) =  0.000035D0
         ELSE
            GO TO 104
         ENDIF
      ELSE IF (IPOT .EQ. 10) THEN
         IF (C12O16) THEN
            XMAT(0,0,6) =  0.000002D0
            XMAT(1,0,6) = -0.000005D0
            XMAT(1,1,6) =  0.000023D0
            XMAT(2,0,6) =  0.000011D0
            XMAT(2,1,6) = -0.000051D0
            XMAT(2,2,6) =  0.000165D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,6) =  0.000002D0
            XMAT(1,0,6) = -0.000004D0
            XMAT(1,1,6) =  0.000021D0
            XMAT(2,0,6) =  0.000010D0
            XMAT(2,1,6) = -0.000046D0
            XMAT(2,2,6) =  0.000148D0
         ELSE
            GO TO 104
         ENDIF
      ELSE
         GO TO 103
      ENDIF

      XMAT(0,1,6) =  XMAT(1,0,6)
      XMAT(0,2,6) =  XMAT(2,0,6)
      XMAT(1,2,6) =  XMAT(2,1,6)

      IF (N .EQ. 6) GO TO 101

C
C     (r - re)^7
C     ==========
C
      IF (IPOT .EQ. 4) THEN
         IF (C12O16) THEN
            XMAT(0,0,7) =  0.000000D0
            XMAT(1,0,7) = -0.000001D0
            XMAT(1,1,7) =  0.000001D0
            XMAT(2,0,7) =  0.000001D0
            XMAT(2,1,7) = -0.000004D0
            XMAT(2,2,7) =  0.000007D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,7) =  0.000000D0
            XMAT(1,0,7) = -0.000001D0
            XMAT(1,1,7) =  0.000001D0
            XMAT(2,0,7) =  0.000001D0
            XMAT(2,1,7) = -0.000003D0
            XMAT(2,2,7) =  0.000006D0
         ELSE
            GO TO 104
         ENDIF
      ELSE IF (IPOT .EQ. 10) THEN
         IF (C12O16) THEN
            XMAT(0,0,7) =  0.000000D0
            XMAT(1,0,7) = -0.000001D0
            XMAT(1,1,7) =  0.000005D0
            XMAT(2,0,7) =  0.000002D0
            XMAT(2,1,7) = -0.000015D0
            XMAT(2,2,7) =  0.000049D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,7) =  0.000000D0
            XMAT(1,0,7) = -0.000001D0
            XMAT(1,1,7) =  0.000005D0
            XMAT(2,0,7) =  0.000002D0
            XMAT(2,1,7) = -0.000014D0
            XMAT(2,2,7) =  0.000044D0
         ELSE
            GO TO 104
         ENDIF
      ELSE
         GO TO 103
      ENDIF

      XMAT(0,1,7) =  XMAT(1,0,7)
      XMAT(0,2,7) =  XMAT(2,0,7)
      XMAT(1,2,7) =  XMAT(2,1,7)

      IF (N .EQ. 7) GO TO 101

C
C     (r - re)^8
C     ==========
C
      IF (IPOT .EQ. 4) THEN
         IF (C12O16) THEN
            XMAT(0,0,8) =  0.000000D0
            XMAT(1,0,8) =  0.000000D0
            XMAT(1,1,8) =  0.000000D0
            XMAT(2,0,8) =  0.000000D0
            XMAT(2,1,8) = -0.000001D0
            XMAT(2,2,8) =  0.000002D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,8) =  0.000000D0
            XMAT(1,0,8) =  0.000000D0
            XMAT(1,1,8) =  0.000000D0
            XMAT(2,0,8) =  0.000000D0
            XMAT(2,1,8) = -0.000001D0
            XMAT(2,2,8) =  0.000002D0
         ELSE
            GO TO 104
         ENDIF
      ELSE IF (IPOT .EQ. 10) THEN
         IF (C12O16) THEN
            XMAT(0,0,8) =  0.000000D0
            XMAT(1,0,8) =  0.000000D0
            XMAT(1,1,8) =  0.000001D0
            XMAT(2,0,8) =  0.000001D0
            XMAT(2,1,8) = -0.000004D0
            XMAT(2,2,8) =  0.000016D0
         ELSE IF (C13O16) THEN
            XMAT(0,0,8) =  0.000000D0
            XMAT(1,0,8) =  0.000000D0
            XMAT(1,1,8) =  0.000001D0
            XMAT(2,0,8) =  0.000001D0
            XMAT(2,1,8) = -0.000004D0
            XMAT(2,2,8) =  0.000014D0
         ELSE
            GO TO 104
         ENDIF
      ELSE
         GO TO 103
      ENDIF

      XMAT(0,1,8) =  XMAT(1,0,8)
      XMAT(0,2,8) =  XMAT(2,0,8)
      XMAT(1,2,8) =  XMAT(2,1,8)

  101 RETURN

  102 WRITE(6,'(//,A,/,A,I6,/,A,I6)')
     & 'The expansion order is out of bounds in VIBMAT:',
     & ' - requested orders up to and including ',N,
     & ' - minimum order is 0, maximum order is ',NMAX
      STOP ' !!! Order mismatch in VIBMAT !!! '

  103 WRITE(6,'(//,A,/,A,I6,/,A,/,A)')
     & 'Potential definition is out of bounds in VIBMAT:',
     & ' - IPOT is given as (from COMMON/POTDEF/) ',IPOT,
     & ' - IPOT must be either 4 (empirical CO states)',
     & '   or 10 (CCSD(T)/aug-cc-pVQZ CO states)'
      STOP ' !!! IPOT error in VIBMAT !!! '

  104 WRITE(6,'(//,A,/,A,I6,/,A,I6,/,A)')
     & 'Isotope combination not implemented in VIBMAT:',
     & ' - ISOC is given as (from COMMON/ISODEF/) ',ISOC,
     & ' - ISOO is given as (from COMMON/ISODEF/) ',ISOO,
     & ' PROGRAM WILL BE ABORTED !!! '
      STOP ' !!! Isotope mismatch in VIBMAT !!! '
      END
      SUBROUTINE VIBMAT_ORIG(XMAT,N)
C
C     CO vibrational matrix elements from MOLCAS.
C
C        <v1(CO)| (r - re)^k |v2(CO)>
C
C     for v1(CO), v2(CO) = 0,1,2
C     and k              = 0,1,2,3,...,N
C
C     Unit: bohr^k
C
C     NOTE: N can be at most NMAX (defined as parameter).
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XMAT(0:2,0:2,0:N)
      PARAMETER (NMAX = 8)
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)
C
C     Check N.
C     ========
C
      IF ((N.LT.0) .OR. (N.GT.NMAX)) GO TO 102
C
C     (r - re)^0 = 1
C     ==============
C
      XMAT(0,0,0) =  ONE
      XMAT(1,0,0) =  ZERO
      XMAT(1,1,0) =  ONE
      XMAT(2,0,0) =  ZERO
      XMAT(2,1,0) =  ZERO
      XMAT(2,2,0) =  ONE
      XMAT(0,1,0) =  XMAT(1,0,0)
      XMAT(0,2,0) =  XMAT(2,0,0)
      XMAT(1,2,0) =  XMAT(2,1,0)

      IF (N .EQ. 0) GO TO 101
C
C     (r - re)^1
C     ==========
C
      XMAT(0,0,1) =  0.007708D0
      XMAT(1,0,1) = -0.063792D0
      XMAT(1,1,1) =  0.023255D0
      XMAT(2,0,1) = -0.003643D0
      XMAT(2,1,1) = -0.090479D0
      XMAT(2,2,1) =  0.038996D0
      XMAT(0,1,1) =  XMAT(1,0,1)
      XMAT(0,2,1) =  XMAT(2,0,1)
      XMAT(1,2,1) =  XMAT(2,1,1)

      IF (N .EQ. 1) GO TO 101

C
C     (r - re)^2
C     ==========
C
      XMAT(0,0,2) =  0.004142D0
      XMAT(1,0,2) = -0.001643D0
      XMAT(1,1,2) =  0.012837D0
      XMAT(2,0,2) =  0.005641D0
      XMAT(2,1,2) = -0.004689D0
      XMAT(2,2,2) =  0.022154D0
      XMAT(0,1,2) =  XMAT(1,0,2)
      XMAT(0,2,2) =  XMAT(2,0,2)
      XMAT(1,2,2) =  XMAT(2,1,2)

      IF (N .EQ. 2) GO TO 101

C
C     (r - re)^3
C     ==========
C
      XMAT(0,0,3) =  0.000116D0
      XMAT(1,0,3) = -0.000818D0
      XMAT(1,1,3) =  0.000765D0
      XMAT(2,0,3) =  0.000264D0
      XMAT(2,1,3) = -0.002428D0
      XMAT(2,2,3) =  0.002109D0
      XMAT(0,1,3) =  XMAT(1,0,3)
      XMAT(0,2,3) =  XMAT(2,0,3)
      XMAT(1,2,3) =  XMAT(2,1,3)

      IF (N .EQ. 3) GO TO 101

C
C     (r - re)^4
C     ==========
C
      XMAT(0,0,4) =  0.000052D0
      XMAT(1,0,4) = -0.000046D0
      XMAT(1,1,4) =  0.000286D0
      XMAT(2,0,4) =  0.000151D0
      XMAT(2,1,4) = -0.000235D0
      XMAT(2,2,4) =  0.000810D0
      XMAT(0,1,4) =  XMAT(1,0,4)
      XMAT(0,2,4) =  XMAT(2,0,4)
      XMAT(1,2,4) =  XMAT(2,1,4)

      IF (N .EQ. 4) GO TO 101

C
C     (r - re)^5
C     ==========
C
      XMAT(0,0,5) =  0.000003D0
      XMAT(1,0,5) = -0.000018D0
      XMAT(1,1,5) =  0.000028D0
      XMAT(2,0,5) =  0.000013D0
      XMAT(2,1,5) = -0.000085D0
      XMAT(2,2,5) =  0.000116D0
      XMAT(0,1,5) =  XMAT(1,0,5)
      XMAT(0,2,5) =  XMAT(2,0,5)
      XMAT(1,2,5) =  XMAT(2,1,5)

      IF (N .EQ. 5) GO TO 101

C
C     (r - re)^6
C     ==========
C
      XMAT(0,0,6) =  0.000001D0
      XMAT(1,0,6) = -0.000002D0
      XMAT(1,1,6) =  0.000009D0
      XMAT(2,0,6) =  0.000005D0
      XMAT(2,1,6) = -0.000012D0
      XMAT(2,2,6) =  0.000038D0
      XMAT(0,1,6) =  XMAT(1,0,6)
      XMAT(0,2,6) =  XMAT(2,0,6)
      XMAT(1,2,6) =  XMAT(2,1,6)

      IF (N .EQ. 6) GO TO 101

C
C     (r - re)^7
C     ==========
C
      XMAT(0,0,7) =  0.000000D0
      XMAT(1,0,7) = -0.000001D0
      XMAT(1,1,7) =  0.000001D0
      XMAT(2,0,7) =  0.000001D0
      XMAT(2,1,7) = -0.000004D0
      XMAT(2,2,7) =  0.000007D0
      XMAT(0,1,7) =  XMAT(1,0,7)
      XMAT(0,2,7) =  XMAT(2,0,7)
      XMAT(1,2,7) =  XMAT(2,1,7)

      IF (N .EQ. 7) GO TO 101

C
C     (r - re)^8
C     ==========
C
      XMAT(0,0,8) =  0.000000D0
      XMAT(1,0,8) =  0.000000D0
      XMAT(1,1,8) =  0.000000D0
      XMAT(2,0,8) =  0.000000D0
      XMAT(2,1,8) = -0.000001D0
      XMAT(2,2,8) =  0.000002D0
      XMAT(0,1,8) =  XMAT(1,0,8)
      XMAT(0,2,8) =  XMAT(2,0,8)
      XMAT(1,2,8) =  XMAT(2,1,8)

  101 RETURN

  102 WRITE(6,'(//,A,/,A,I6,/,A,I6)')
     & 'The expansion order is out of bounds in VIBMAT:',
     & ' - requested orders up to and including ',N,
     & ' - minimum order is 0, maximum order is ',NMAX
      STOP ' !!! Order mismatch in VIBMAT !!! '
      END
      SUBROUTINE GETPAR3(BPAR,DPAR,GPAR,C60PAR,C62PAR,C71PAR,C73PAR)
C
C     Our "full" aug-cc-pVQZ-33211 surface.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION BPAR(0:5,0:2)
      DIMENSION DPAR(0:5,0:2)
      DIMENSION GPAR(0:5,0:3,0:2)
      DIMENSION C60PAR(0:2), C62PAR(0:2)
      DIMENSION C71PAR(0:2), C73PAR(0:2)
C
      BPAR(0,0) =  -2.773531400530731D0
      BPAR(1,0) =   0.355978611060635D0
      BPAR(2,0) =  -0.038461094313790D0
      BPAR(3,0) =  -0.065597122410946D0
      BPAR(4,0) =  -0.017415924916096D0
      BPAR(5,0) =   0.009165061618361D0
C
      BPAR(0,1) =   0.383748086062260D0
      BPAR(1,1) =  -0.267512471179551D0
      BPAR(2,1) =   0.237733994757126D0
      BPAR(3,1) =  -0.466599098913902D0
      BPAR(4,1) =   0.222943449032861D0
      BPAR(5,1) =  -0.019012147801577D0
C
      BPAR(0,2) =   4.896202058362745D0
      BPAR(1,2) =  -9.446664416716875D0
      BPAR(2,2) =  -9.207034029053412D0
      BPAR(3,2) =   6.011576154886567D0
      BPAR(4,2) =  -0.191562479666016D0
      BPAR(5,2) =  -0.365735693068647D0
C
      DPAR(0,0) =  12.447246080504405D0
      DPAR(1,0) =  -2.612275572240643D0
      DPAR(2,0) =   0.622683086317908D0
      DPAR(3,0) =   0.252752167976621D0
      DPAR(4,0) =   0.090218496365782D0
      DPAR(5,0) =  -0.048238070953660D0
C
      DPAR(0,1) =   1.837823678601632D0
      DPAR(1,1) =  -1.236773331233324D0
      DPAR(2,1) =  -0.342754115413572D0
      DPAR(3,1) =   1.438028441494077D0
      DPAR(4,1) =  -0.369337751726034D0
      DPAR(5,1) =  -0.081322869855238D0
C
      DPAR(0,2) = -24.074189169057551D0
      DPAR(1,2) =  32.523132561232011D0
      DPAR(2,2) =  27.834718735481708D0
      DPAR(3,2) = -25.905392215443936D0
      DPAR(4,2) =   4.031779135514781D0
      DPAR(5,2) =   0.487956040833263D0
C
      GPAR(0,0,0) =   1.339631498881794D0
      GPAR(1,0,0) =   2.948722284496854D0
      GPAR(2,0,0) =   2.031907693472959D0
      GPAR(3,0,0) =   0.696650687381108D0
      GPAR(4,0,0) =   0.121157437425799D0
      GPAR(5,0,0) =  -0.058540074115915D0
      GPAR(0,1,0) =  -0.716696943349055D0
      GPAR(1,1,0) =  -1.799379832289937D0
      GPAR(2,1,0) =  -1.272374773587762D0
      GPAR(3,1,0) =  -0.422537827515598D0
      GPAR(4,1,0) =  -0.066726893855753D0
      GPAR(5,1,0) =   0.046805922626871D0
      GPAR(0,2,0) =   0.133466525438003D0
      GPAR(1,2,0) =   0.363354836255877D0
      GPAR(2,2,0) =   0.278274964310456D0
      GPAR(3,2,0) =   0.089554042766288D0
      GPAR(4,2,0) =   0.013086477976977D0
      GPAR(5,2,0) =  -0.011042103563629D0
      GPAR(0,3,0) =  -0.009525087912154D0
      GPAR(1,3,0) =  -0.025258865602849D0
      GPAR(2,3,0) =  -0.021827843859677D0
      GPAR(3,3,0) =  -0.006926659963151D0
      GPAR(4,3,0) =  -0.000880022720031D0
      GPAR(5,3,0) =   0.000834993057544D0
C
      GPAR(0,0,1) =   0.265157869877913D0
      GPAR(1,0,1) =  -0.623695350637235D0
      GPAR(2,0,1) =   1.674633318254092D0
      GPAR(3,0,1) =   1.832031291880343D0
      GPAR(4,0,1) =  -0.630177209092439D0
      GPAR(5,0,1) =  -0.514436144698261D0
      GPAR(0,1,1) =  -0.767216076396932D0
      GPAR(1,1,1) =   0.204506170017565D0
      GPAR(2,1,1) =  -1.026962739406349D0
      GPAR(3,1,1) =  -1.499270438208511D0
      GPAR(4,1,1) =   0.617151174631303D0
      GPAR(5,1,1) =   0.374467312887415D0
      GPAR(0,2,1) =   0.291554027998303D0
      GPAR(1,2,1) =   0.073924141379195D0
      GPAR(2,2,1) =   0.183942919605771D0
      GPAR(3,2,1) =   0.416043391498683D0
      GPAR(4,2,1) =  -0.190246683084211D0
      GPAR(5,2,1) =  -0.090246903296188D0
      GPAR(0,3,1) =  -0.030968444075873D0
      GPAR(1,3,1) =  -0.020411320334507D0
      GPAR(2,3,1) =  -0.009669711535675D0
      GPAR(3,3,1) =  -0.038366340553356D0
      GPAR(4,3,1) =   0.018653689407513D0
      GPAR(5,3,1) =   0.007343003097835D0
C
      GPAR(0,0,2) =   8.785016107647921D0
      GPAR(1,0,2) =  -4.895140952257445D0
      GPAR(2,0,2) = -19.071047814522810D0
      GPAR(3,0,2) =  -0.651246137053238D0
      GPAR(4,0,2) =   8.433900748345344D0
      GPAR(5,0,2) =  -7.514200564005597D0
      GPAR(0,1,2) =  -2.600511428982215D0
      GPAR(1,1,2) =   4.439659662444800D0
      GPAR(2,1,2) =  14.877518045476110D0
      GPAR(3,1,2) =   4.915499523879274D0
      GPAR(4,1,2) =  -7.167713194561587D0
      GPAR(5,1,2) =   6.754486599741076D0
      GPAR(0,2,2) =  -0.033889744561149D0
      GPAR(1,2,2) =  -1.097798816676007D0
      GPAR(2,2,2) =  -3.182774594628603D0
      GPAR(3,2,2) =  -1.975283365814364D0
      GPAR(4,2,2) =   2.039158636810153D0
      GPAR(5,2,2) =  -1.870762601953492D0
      GPAR(0,3,2) =   0.036634516644001D0
      GPAR(1,3,2) =   0.065423196677937D0
      GPAR(2,3,2) =   0.170669053301343D0
      GPAR(3,3,2) =   0.190098586427090D0
      GPAR(4,3,2) =  -0.194641294138869D0
      GPAR(5,3,2) =   0.163741529410431D0
C
      C60PAR(0) = -1996.575537892978900D0
      C60PAR(1) =    -0.307468414487437D0
      C60PAR(2) =     5.743607568908502D0
C
      C62PAR(0) =  -388.206699666093300D0
      C62PAR(1) =    -0.014941745228257D0
      C62PAR(2) =    -0.235609028363907D0
C
      C71PAR(0) =  4059.088160207028977D0
      C71PAR(1) =     0.021629800103675D0
      C71PAR(2) =    -0.101783055744035D0
C
      C73PAR(0) =   558.322859066182559D0
      C73PAR(1) =     0.017816070337537D0
      C73PAR(2) =    -0.120764317323077D0
C
      RETURN
      END
      SUBROUTINE PARINIT(BC,DC,G0,G1,G2,G3,C60,C62,C71,C73)
C
C     Cybulski's aug-cc-pVTZ-33221 surface for rCO = re = 2.132 bohr.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION BC(6), DC(6), G0(6)
      DIMENSION G1(6), G2(6), G3(6)
C
      BC(1) = -2.81869D0
      BC(2) =  0.41714D0
      BC(3) = -0.09493D0
      BC(4) = -0.07197D0
      BC(5) = -0.03237D0
      BC(6) =  0.02129D0
C
      DC(1) = 12.42564D0
      DC(2) = -2.62075D0
      DC(3) =  0.69161D0
      DC(4) =  0.24167D0
      DC(5) =  0.13941D0
      DC(6) = -0.07974D0
C
      G0(1) =  1.36083D0
      G0(2) =  2.88520D0
      G0(3) =  2.17234D0
      G0(4) =  0.58897D0
      G0(5) =  0.17074D0
      G0(6) = -0.09171D0
C
      G1(1) = -0.71415D0
      G1(2) = -1.68790D0
      G1(3) = -1.38119D0
      G1(4) = -0.30349D0
      G1(5) = -0.10148D0
      G1(6) =  0.07914D0
C
      G2(1) =  0.13467D0
      G2(2) =  0.31669D0
      G2(3) =  0.31489D0
      G2(4) =  0.05272D0
      G2(5) =  0.02339D0
      G2(6) = -0.02005D0
C
      G3(1) = -0.01027D0
      G3(2) = -0.01989D0
      G3(3) = -0.02619D0
      G3(4) = -0.00360D0
      G3(5) = -0.00203D0
      G3(6) =  0.00157D0
C
      C60   = -1996.5746D0
      C62   =  -388.20555D0
C
      C71   =  4059.0878D0
      C73   =   558.32279D0
C
      RETURN
      END
      SUBROUTINE PARINIT2(BC,DC,G0,G1,G2,G3,C60,C62,C71,C73)
C
C     Our aug-cc-pVQZ-33211 surface for rCO = re = 2.132 bohr.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION BC(6), DC(6), G0(6)
      DIMENSION G1(6), G2(6), G3(6)
C
      BC(1) = -2.77201D0
      BC(2) =  0.3548550
      BC(3) = -0.0386127D0
      BC(4) = -0.065697D0
      BC(5) = -0.0166002D0
      BC(6) =  0.00926849D0
C
      DC(1) = 12.4464D0
      DC(2) = -2.61161D0
      DC(3) =  0.622713D0
      DC(4) =  0.253088D0
      DC(5) =  0.0879602D0
      DC(6) = -0.0485837D0
C
      G0(1) =  1.33986D0
      G0(2) =  2.94754D0
      G0(3) =  2.03162D0
      G0(4) =  0.699013D0
      G0(5) =  0.121969D0
      G0(6) = -0.0580699D0
C
      G1(1) = -0.717386D0
      G1(2) = -1.79928D0
      G1(3) = -1.27189D0
      G1(4) = -0.424176D0
      G1(5) = -0.0669338D0
      G1(6) =  0.0466011D0
C
      G2(1) =  0.133656D0
      G2(2) =  0.363559D0
      G2(3) =  0.27804D0
      G2(4) =  0.089929D0
      G2(5) =  0.0130131D0
      G2(6) = -0.0110455D0
C
      G3(1) = -0.00953441D0
      G3(2) = -0.0252936D0
      G3(3) = -0.0217967D0
      G3(4) = -0.00695421D0
      G3(5) = -0.000862247D0
      G3(6) =  0.000841174D0
C
      C60   = -1996.58D0
      C62   =  -388.207D0
C
      C71   =  4059.09D0
      C73   =   558.323D0
C
      RETURN
      END
      SUBROUTINE PARINIT3(BC,DC,G0,G1,G2,G3,C60,C62,C71,C73,RCO)
C
C     Our "full" aug-cc-pVQZ-33211 surface.
C     Parameters are calculated for a given CO distance, RCO (in bohr).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION BC(6), DC(6), G0(6)
      DIMENSION G1(6), G2(6), G3(6)
C
      DIMENSION BCPAR(6,3), DCPAR(6,3), G0PAR(6,3)
      DIMENSION G1PAR(6,3), G2PAR(6,3), G3PAR(6,3)
      DIMENSION C60PAR(3), C62PAR(3)
      DIMENSION C71PAR(3), C73PAR(3)
      DIMENSION RD(3)
C
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)
      PARAMETER (BOHR = 0.52917725D0)
      PARAMETER (RE   = 2.13200000D0)
C
      BCPAR(1,1) =  -2.773531400530731D0
      BCPAR(2,1) =   0.355978611060635D0
      BCPAR(3,1) =  -0.038461094313790D0
      BCPAR(4,1) =  -0.065597122410946D0
      BCPAR(5,1) =  -0.017415924916096D0
      BCPAR(6,1) =   0.009165061618361D0
      BCPAR(1,2) =   0.383748086062260D0
      BCPAR(2,2) =  -0.267512471179551D0
      BCPAR(3,2) =   0.237733994757126D0
      BCPAR(4,2) =  -0.466599098913902D0
      BCPAR(5,2) =   0.222943449032861D0
      BCPAR(6,2) =  -0.019012147801577D0
      BCPAR(1,3) =   4.896202058362745D0
      BCPAR(2,3) =  -9.446664416716875D0
      BCPAR(3,3) =  -9.207034029053412D0
      BCPAR(4,3) =   6.011576154886567D0
      BCPAR(5,3) =  -0.191562479666016D0
      BCPAR(6,3) =  -0.365735693068647D0
C
      DCPAR(1,1) =  12.447246080504405D0
      DCPAR(2,1) =  -2.612275572240643D0
      DCPAR(3,1) =   0.622683086317908D0
      DCPAR(4,1) =   0.252752167976621D0
      DCPAR(5,1) =   0.090218496365782D0
      DCPAR(6,1) =  -0.048238070953660D0
      DCPAR(1,2) =   1.837823678601632D0
      DCPAR(2,2) =  -1.236773331233324D0
      DCPAR(3,2) =  -0.342754115413572D0
      DCPAR(4,2) =   1.438028441494077D0
      DCPAR(5,2) =  -0.369337751726034D0
      DCPAR(6,2) =  -0.081322869855238D0
      DCPAR(1,3) = -24.074189169057551D0
      DCPAR(2,3) =  32.523132561232011D0
      DCPAR(3,3) =  27.834718735481708D0
      DCPAR(4,3) = -25.905392215443936D0
      DCPAR(5,3) =   4.031779135514781D0
      DCPAR(6,3) =   0.487956040833263D0
C
      G0PAR(1,1) =   1.339631498881794D0
      G0PAR(2,1) =   2.948722284496854D0
      G0PAR(3,1) =   2.031907693472959D0
      G0PAR(4,1) =   0.696650687381108D0
      G0PAR(5,1) =   0.121157437425799D0
      G0PAR(6,1) =  -0.058540074115915D0
      G0PAR(1,2) =   0.265157869877913D0
      G0PAR(2,2) =  -0.623695350637235D0
      G0PAR(3,2) =   1.674633318254092D0
      G0PAR(4,2) =   1.832031291880343D0
      G0PAR(5,2) =  -0.630177209092439D0
      G0PAR(6,2) =  -0.514436144698261D0
      G0PAR(1,3) =   8.785016107647921D0
      G0PAR(2,3) =  -4.895140952257445D0
      G0PAR(3,3) = -19.071047814522810D0
      G0PAR(4,3) =  -0.651246137053238D0
      G0PAR(5,3) =   8.433900748345344D0
      G0PAR(6,3) =  -7.514200564005597D0
C
      G1PAR(1,1) =  -0.716696943349055D0
      G1PAR(2,1) =  -1.799379832289937D0
      G1PAR(3,1) =  -1.272374773587762D0
      G1PAR(4,1) =  -0.422537827515598D0
      G1PAR(5,1) =  -0.066726893855753D0
      G1PAR(6,1) =   0.046805922626871D0
      G1PAR(1,2) =  -0.767216076396932D0
      G1PAR(2,2) =   0.204506170017565D0
      G1PAR(3,2) =  -1.026962739406349D0
      G1PAR(4,2) =  -1.499270438208511D0
      G1PAR(5,2) =   0.617151174631303D0
      G1PAR(6,2) =   0.374467312887415D0
      G1PAR(1,3) =  -2.600511428982215D0
      G1PAR(2,3) =   4.439659662444800D0
      G1PAR(3,3) =  14.877518045476110D0
      G1PAR(4,3) =   4.915499523879274D0
      G1PAR(5,3) =  -7.167713194561587D0
      G1PAR(6,3) =   6.754486599741076D0
C
      G2PAR(1,1) =   0.133466525438003D0
      G2PAR(2,1) =   0.363354836255877D0
      G2PAR(3,1) =   0.278274964310456D0
      G2PAR(4,1) =   0.089554042766288D0
      G2PAR(5,1) =   0.013086477976977D0
      G2PAR(6,1) =  -0.011042103563629D0
      G2PAR(1,2) =   0.291554027998303D0
      G2PAR(2,2) =   0.073924141379195D0
      G2PAR(3,2) =   0.183942919605771D0
      G2PAR(4,2) =   0.416043391498683D0
      G2PAR(5,2) =  -0.190246683084211D0
      G2PAR(6,2) =  -0.090246903296188D0
      G2PAR(1,3) =  -0.033889744561149D0
      G2PAR(2,3) =  -1.097798816676007D0
      G2PAR(3,3) =  -3.182774594628603D0
      G2PAR(4,3) =  -1.975283365814364D0
      G2PAR(5,3) =   2.039158636810153D0
      G2PAR(6,3) =  -1.870762601953492D0
C
      G3PAR(1,1) =  -0.009525087912154D0
      G3PAR(2,1) =  -0.025258865602849D0
      G3PAR(3,1) =  -0.021827843859677D0
      G3PAR(4,1) =  -0.006926659963151D0
      G3PAR(5,1) =  -0.000880022720031D0
      G3PAR(6,1) =   0.000834993057544D0
      G3PAR(1,2) =  -0.030968444075873D0
      G3PAR(2,2) =  -0.020411320334507D0
      G3PAR(3,2) =  -0.009669711535675D0
      G3PAR(4,2) =  -0.038366340553356D0
      G3PAR(5,2) =   0.018653689407513D0
      G3PAR(6,2) =   0.007343003097835D0
      G3PAR(1,3) =   0.036634516644001D0
      G3PAR(2,3) =   0.065423196677937D0
      G3PAR(3,3) =   0.170669053301343D0
      G3PAR(4,3) =   0.190098586427090D0
      G3PAR(5,3) =  -0.194641294138869D0
      G3PAR(6,3) =   0.163741529410431D0
C
      C60PAR(1) = -1996.575537892978900D0
      C60PAR(2) =    -0.307468414487437D0
      C60PAR(3) =     5.743607568908502D0
C
      C62PAR(1) =  -388.206699666093300D0
      C62PAR(2) =    -0.014941745228257D0
      C62PAR(3) =    -0.235609028363907D0
C
      C71PAR(1) =  4059.088160207028977D0
      C71PAR(2) =     0.021629800103675D0
      C71PAR(3) =    -0.101783055744035D0
C
      C73PAR(1) =   558.322859066182559D0
      C73PAR(2) =     0.017816070337537D0
      C73PAR(3) =    -0.120764317323077D0
C
      DIFF  = (RCO - RE)*BOHR
      RD(1) = ONE
      RD(2) = DIFF
      RD(3) = DIFF*DIFF
C
      C60 = ZERO
      C62 = ZERO
      C71 = ZERO
      C73 = ZERO
      DO I = 1,6
         BC(I) = ZERO
         DC(I) = ZERO
         G0(I) = ZERO
         G1(I) = ZERO
         G2(I) = ZERO
         G3(I) = ZERO
      ENDDO

      DO J = 1,3
         DO I = 1,6
            BC(I) = BC(I) + BCPAR(I,J)*RD(J)
            DC(I) = DC(I) + DCPAR(I,J)*RD(J)
            G0(I) = G0(I) + G0PAR(I,J)*RD(J)
            G1(I) = G1(I) + G1PAR(I,J)*RD(J)
            G2(I) = G2(I) + G2PAR(I,J)*RD(J)
            G3(I) = G3(I) + G3PAR(I,J)*RD(J)
         ENDDO
         C60 = C60 + C60PAR(J)*RD(J)
         C62 = C62 + C62PAR(J)*RD(J)
         C71 = C71 + C71PAR(J)*RD(J)
         C73 = C73 + C73PAR(J)*RD(J)
      ENDDO
C
      RETURN
      END
      SUBROUTINE POTINP
C
C     Read potential definitions from file and
C     store in common block POTDEF.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON / POTDEF / IPOT, IV1, IV2, ISTAT, LUSTAT
      COMMON / ISODEF / ISOC, ISOO

      CHARACTER*6 SECNAM
      PARAMETER (SECNAM = 'POTINP')

      CHARACTER*80 TEXT

      LOGICAL ERR1, ERR2, WARN

      CHARACTER*7 FINP
      PARAMETER (LUINP = 7, FINP = 'pot.inp')
      CHARACTER*7 FSTAT
      PARAMETER (FSTAT = 'potstat')

      PARAMETER (MXREF = 5)
      CHARACTER*50 POTTYP,METHOD
      CHARACTER*50 REF(MXREF)

      LUSTAT = 37

      WRITE(6,'(//,15X,A,A,/,15X,A,/)')
     & SECNAM,': INFORMATION ABOUT THE POTENTIAL',
     &  '======================================='

      OPEN(LUINP,STATUS='OLD',FORM='FORMATTED',FILE=FINP,ERR=101)

      READ(LUINP,*) NTEXT
      DO ITEXT = 1,NTEXT
         READ(LUINP,'(A80)') TEXT
      ENDDO
      READ(LUINP,*) IPOT,IV1,IV2,ISTAT
      READ(LUINP,*) ISOC, ISOO

      WARN = .FALSE.
      IF (IPOT .EQ. 1) THEN
         POTTYP = '2-DIMENSIONAL INTERMOLECULAR                      '
         METHOD = 'CCSD(T)/AUG-CC-PVTZ-33221 AT R(CO) = 2.132 BOHR   '
         NREF   = 1
         REF(1) = 'J. CHEM. PHYS. 112, 4604 (2000)                   '
      ELSE IF (IPOT .EQ. 2) THEN
         POTTYP = '2-DIMENSIONAL INTERMOLECULAR                      '
         METHOD = 'CCSD(T)/AUG-CC-PVQZ-33211 AT R(CO) = 2.132 BOHR   '
         NREF   = 0
      ELSE IF (IPOT .EQ. 3) THEN
         POTTYP = '3-DIMENSIONAL INTERMOLECULAR                      '
         METHOD = 'CCSD(T)/AUG-CC-PVQZ-33211                         '
         NREF   = 0
      ELSE IF (IPOT .EQ. 4) THEN
         POTTYP = '3-DIMENSIONAL INTERMOLECULAR                      '
         METHOD = 'CCSD(T)/AUG-CC-PVQZ-33211 (empirical CO states)   '
         NREF   = 1
         REF(1) = 'CO POTENTIAL: J. CHEM. SOC. FAR. II 79, 323 (1983)'
         ERR1 = (IV1.LT.0) .OR. (IV1.GT.2)
         ERR2 = (IV2.LT.0) .OR. (IV2.GT.2)
         IF (ERR1 .OR. ERR2) WARN = .TRUE.
      ELSE IF (IPOT .EQ. 5) THEN
         POTTYP = '3-DIMENSIONAL INTERMOLECULAR PLUS DIATOMIC        '
         METHOD = 'CCSD(T)/AUG-CC-PVQZ-33211 AND EMPIRICAL           '
         NREF   = 1
         REF(1) = 'CO POTENTIAL: J. CHEM. SOC. FAR. II 79, 323 (1983)'
      ELSE IF (IPOT .EQ. 6) THEN
         POTTYP = '1-DIMENSIONAL EMPIRICAL CO POTENTIAL              '
         METHOD = 'EMPIRICAL                                         '
         NREF   = 1
         REF(1) = 'J. CHEM. SOC. FARADAY TRANS. II 79, 323 (1983)    '
      ELSE IF (IPOT .EQ. 7) THEN
         POTTYP = '2-DIMENSIONAL INTERMOLECULAR                      '
         METHOD = 'CCSD(T)/AUG-CC-PVQZ-33211 AT R(CO) = 1.898 BOHR   '
         NREF   = 0
      ELSE IF (IPOT .EQ. 8) THEN
         POTTYP = '2-DIMENSIONAL INTERMOLECULAR                      '
         METHOD = 'CCSD(T)/AUG-CC-PVQZ-33211 AT R(CO) = 2.234 BOHR   '
         NREF   = 0
      ELSE IF (IPOT .EQ. 9) THEN
         POTTYP = '2-D INTERMOLECULAR R(CO) = 2.132 BOHR + DIATOMIC  '
         METHOD = 'CCSD(T)/AUG-CC-PVQZ-33211 AND EMPIRICAL           '
         NREF   = 1
         REF(1) = 'CO POTENTIAL: J. CHEM. SOC. FAR. II 79, 323 (1983)'
      ELSE IF (IPOT .EQ. 10) THEN
         POTTYP = '3-DIMENSIONAL INTERMOLECULAR                      '
         METHOD = 'CCSD(T)/AUG-CC-PVQZ-33211 (CCSD(T)/aQZ CO states) '
         NREF   = 0
         ERR1 = (IV1.LT.0) .OR. (IV1.GT.2)
         ERR2 = (IV2.LT.0) .OR. (IV2.GT.2)
         IF (ERR1 .OR. ERR2) WARN = .TRUE.
      ELSE
         POTTYP = 'UNKNOWN !!!                                       '
         METHOD = 'UNKNOWN !!!                                       '
         WARN   = .TRUE.
         NREF   = 0
      ENDIF

      WRITE(6,1002) POTTYP
      WRITE(6,1003) METHOD
      IF ((IPOT.EQ.4) .OR. (IPOT.EQ.10)) THEN
         IF (IV1 .EQ. IV2) THEN
            WRITE(6,1004) 'AVERAGE TAKEN FOR CO STATE ',IV1
         ELSE
            WRITE(6,1005) 'ELEMENT TAKEN BETWEEN CO STATES ',IV1,
     &                    ' AND ',IV2
         ENDIF
         WRITE(6,1008) 'CARBON : ',ISOC
         WRITE(6,1008) 'OXYGEN : ',ISOO
      ENDIF
      IF (NREF .GT. 0) THEN
         WRITE(6,1006) REF(1)
         DO IREF = 2,NREF
            WRITE(6,1007) REF(IREF)
         ENDDO
      ENDIF
      IF (WARN) THEN
         WRITE(6,'(/,5X,A,A,A,/,5X,A,I2,A,I2,A,I2,A,I2)')
     &   '*** WARNING *** INPUT WAS NOT UNDERSTOOD BY ',SECNAM,':',
     &   ' IPOT  = ',IPOT,
     &   ' IV1   = ',IV1,
     &   ' IV2   = ',IV2,
     &   ' ISTAT = ',ISTAT
         WRITE(6,'(5X,A)')
     &   ' - PROGRAM CONTINUES NEVERTHELESS !!!'
      ENDIF
      WRITE(6,'(//)')

      IF (ISTAT .EQ. 1) THEN
         OPEN(LUSTAT,STATUS='UNKNOWN',FORM='FORMATTED',FILE=FSTAT)
         WRITE(LUSTAT,1002) POTTYP
         WRITE(LUSTAT,1003) METHOD
         IF ((IPOT.EQ.4) .OR. (IPOT.EQ.10)) THEN
            IF (IV1 .EQ. IV2) THEN
               WRITE(LUSTAT,1004) 'AVERAGE TAKEN FOR CO STATE ',IV1
            ELSE
               WRITE(LUSTAT,1005) 'ELEMENT TAKEN BETWEEN CO STATES ',
     &                             IV1,' AND ',IV2
            ENDIF
            WRITE(LUSTAT,1008) 'CARBON : ',ISOC
            WRITE(LUSTAT,1008) 'OXYGEN : ',ISOO
         ENDIF
         IF (NREF .GT. 0) THEN
            WRITE(LUSTAT,1006) REF(1)
            DO IREF = 2,NREF
               WRITE(LUSTAT,1007) REF(IREF)
            ENDDO
         ENDIF
         WRITE(LUSTAT,'(/)')
      ENDIF

      CLOSE(LUINP,STATUS='KEEP',ERR=102)

      RETURN

  101 WRITE(6,1001)
     & SECNAM,': Error opening unit ',LUINP,' named ',FINP
      STOP

  102 WRITE(6,1001)
     & SECNAM,': Error closing unit ',LUINP,' named ',FINP
      STOP
     
 1001 FORMAT(//,5X,A,A,I2,A)
 1002 FORMAT(5X,'POTENTIAL TYPE    :  ',A)
 1003 FORMAT(5X,'THEORETICAL METHOD:  ',A)
 1004 FORMAT(5X,'VIBRATIONAL MATRIX:  ',A,I2)
 1005 FORMAT(5X,'VIBRATIONAL MATRIX:  ',A,I2,A,I2)
 1006 FORMAT(5X,'REFERENCES        :  ',A)
 1007 FORMAT(5X,'                     ',A)
 1008 FORMAT(5X,'ISOTOPE NUMBER FOR   ',A,I2)
      END
