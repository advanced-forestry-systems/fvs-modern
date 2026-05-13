!== last modified  4-9-2002
  SUBROUTINE SF_ZERO(X1, F1, X2, F2, MODE, TOL)
!           Finds the zero of a univariate single precision function  SF_ZERO
!           Adopted from SZERO found in MATH77 Subr Library of
!               Language Systems   Sterling, VA.
!
!           Code is Copyright by Language Systems.
!           Source, object and executable MAY be distributed with SF_TEST.

!>> 1996       SZERO  JWF     Removed calls AMACH , MESS, SMESS and UMESS
!>> 1993-04-27 SZERO  Krogh   Additions for Conversion to C.
!>> 1993-04-13 SZERO  Krogh   Minor change for new MESS.
!>> 1992-04-08 SZERO  Krogh   Unused label 400 removed.
!>> 1992-01-09 SZERO  Krogh   Moved calc. of XXMXO up (for error msg.)
!>> 1991-11-26 SZERO  Krogh   Converted to new error processor.
!>> 1988-08-14 SZERO  Krogh   Labels runumbered.
!>> 1988-03-07 SZERO  Krogh   Initial code.
! ========== This program has been specialized ==========
! SUBROUTINE TO FIND A BOUNDED ZERO
!
! analysis and coding by Fred T.Krogh at the Jet Propulsion
! Laboratory, Pasadena, Calif.  April 25, 1972.
! Modified for portability, April 1984 by Krogh.
! Algorithmic changes, vars. added to save stmt., Sept. 1987 by Krogh
!
! Parameters in the calling sequence are defined as follows:
!
!  X1  = independent variable
!  F1  = dependent variable --  initially   F1=F(X1).
!        When MODE=1 (or 5) the user is to compute F(X1) given X1
!  X2  = second value of independent variable
!  F2  = F(X2) on the initial entry.  When MODE = 2-4, F2=F(X2) and
!        F1*F2 .le. 0.
!  MODE  is a parameter used for communication between this
!        subroutine and the user. (The user should set MODE
!        only to initialize it to 0 before the first call)
!      =1  compute F(X1) and call $ZERO
!      =2  F(X1) is approximately 0, the iteration is finished
!          and the error criterion is satisfied.
!      =3  same as MODE=2, except the error criterion can
!          not be satisfied.
!      =4  apparently the function has a discontinuity
!          between X1 and X2 -- No zero can be found
!      =5  F1*F2 was greater than zero on the first call, and an attempt
!          to bound the zero on both sides have failed.
!      =6  fatal error -- $ZERO was called after mode was set .ge.2.
!          If $ZERO is called again, the program will be stopped.
!          (Unless MODE is set to 0)
!      <0  If MODE is set <0 and $ZERO is called, no action is taken
!          except that print is turned on for -MODE calls to $ZERO.
!          This print gives all values of X and F used in the iteration.
!  TOL    is the error tolerance
!     TOL.GT.0  Iterate until values of X1 and X2 are known
!              for which abs(X1-X2) .le. tol and F1*F2 .le. 0.
!     TOL.LT.0  Iterate until a value of X1 is found for which
!              abs(F1) .le. abs(TOL).
!     TOL  = 0  Iterate until the zero is determined as
!              precisely as possible.  MODE = 3 is impossible
!              in this case.
!
! Parameters in the calling sequence have the following types
!
  integer MODE
  real             X1, X2, F1, F2, TOL
!
! Usage is as follows (of course, variations are possible.)
!         Somehow one has available X1, F1, X2, and F2 such
!         that F1 = F(X1), F2 = F(X2) and F1*F2 .le. 0.
!         In addition, one should assign a value to TOL.
!     MODE = 0
!***  In the statement below, $ is replaced by an 'S' for single
!***  precision and a 'D' for double.
! XXX call $ZERO(X1,F1,X2,F2,MODE,TOL)
!     go to  (N1,N2,N3,N4,N5,N6), MODE
!  N1 COMPUTE  F1=F(X1)
!     go to XXX
!
!  N4 continue
!  N5 continue
!  N6 stop
!  N3 If you want to -- print results to note that error
!                       is too big.
!  N2 zero is found, do whatever you want to with it.
!
! End of comments explaining usage.
!
! ************************* Usage of internal variables ****************
!
! C0     Parameter = 0.
! C1     Parameter = 1.
! C1P01  Parameter = 1.01
! C1P25  Parameter = 1.25
! C2     Parameter = 2.
! C4     Parameter = 4.
! CP01   Parameter = 0.01
! CP125  Parameter = 1.25
! CP25   Parameter = 0.25
! CP5    Parameter = 0.5
! CP75   Parameter = 0.75
! CP99   Parameter = 0.99
! R1MACH Gets constants associated with floating point arithmetic.
! DFDXXX = (XXMXL/FFMFL) * (est. deriv. of f w.r.t. x at x = XX).  All
!   derivatives are base on a second degree polynonial that interpolates
!   the last three points generated.
! DFDXXX = (XXMXO/FFMFL) * (est. deriv. of f w.r.t. x at x = X0).
! DIV    If artificial subdivision of the interval is used, determines
!   the amount of the sudivision.  (-XXMXOO * DIV / (1. + DIV))
! SMESS  Prints error messages.
! DXDFFF = (FFMFL/XXMXL) * (est. deriv. of x w.r.t. f at f = FF).
! DXDFFO = (FFMFO/XXMXL) * (est. deriv. of x w.r.t. f at f = FO).
! F1     (formal arg.) The last value of F computed, on return the value
!   of F(X1).
! F2     (formal arg.) The other initial value provided for F.  Set to
!   the value of F(X2) on returns.
! FDAT   Temporary storage for floating point values for messages.
! FF     Value of F1 after F is computed.
! FFDFO  FF / FO
! FFMFB  FFMFL + FLMFB = FF - value of FF 2 iterations back.
! FFMFL  FF - FL
! FL     Value of FF from the previous iteration.
! FLMFB  Value of FFMFL from the previous iteration
! FO     F(XO)
! I      Comments for LCHNG define how I is set.
! IDAT   Temporary storage for integer values for messages.
! INDIC  Internal value to be assigned to MODE on exit.  Equivalenced to
!   MACT(3) and set to 0 initially.
! J      This is 1 if FF .le. 0., and is 2 if FF > 0.
! KNKP   KNKP(J) (see J above) is set to 0 whenever there are signs of
!   decent convergence.  It is counted up when convergence is slow.
! KS     =-1  initially,
!        = 0  whenever F changes sign, otherwise
!        = number of times the sign of F has remained the same
! KTYP   = 1 if interpolation was used to get the last iterate, = 0 if
!   an artificial subdivision was used.
! LCHG  the J-th continuation in the data statement for LCHG below gives
! new states for the case when the state number is J-1.  State 0 is the
! initial state.  The I-th entry on a row gives the state for case on I
! as follows:  (TP is the ratio (new f) / (last f of same sign)
!    I = 1   TP < 0.01
!    I = 2   .01 <= TP < 1
!    I = 3   TP = 1
!    I = 4   1 < TP <= 4
!    I = 5   TP > 4.
! States are as follows:
!    0   initial state, or big increase, or small increase in state 0
!    1   after big decrease, perhaps followed by small decreases
!    2   after one or more small decreases from state 0
!    3   one or more small increases from state 2
!    4   one or more small decreases from state 3
!    5   decision made that noise is a problem on this side of zero.
! LINIT  - the maximum number of iterations that can be taken with out
!   getting a sign change in F.
! LMODE  The value of MODE the last time in this routine.
! LNLP
! LTXTxx Names of this form are used in setting up data statements for
!   error messages.  These names are generated automatically by PMESS,
!   the program that makes up these messages.
! MACT   This array difines the actions to be taken by the error message
!   program.  See comments in MESS for details.
! MACT1  As for MACT except used for the diagnostic print.
! MExxxx Parameters defining constants used for interaction with the
!   error message program MESS.  See comments there for definitions.
! MLOC   Contains locations in MTXTAA for error messages.
! MODE   (formal) See comments above.
! MTXTAA Text for error messages.
! MTXTAB Text for diagnostic message.
! MTXTAC Text for diagnostic message.
! NP     If > 0, gives number of iterations till diagnostic print stops.
! QFM
! QXM
! RND    Largest relative difference between succesive floating point
!   numbers.
! SMALL  .5 / (RND * largest floating point number)
! TOL    (Formal) See description above.
! TOLX   Actually tolerance required for accuracy in X.  Usually =
!   max(TOL, XRND).  It can be cut by a factor of 2 for use in setting
!   bounds on an acceptable interval.
! TP     Ordinarily the ratio (FF / prev. FF of the same sign.
! TP1    Used for temporary storage.
! X1     (Formal) Value of x where F is to be computed, and value
!   returned for the zero after convergence.
! X2     (Formal) Initially other value of x where F is given.  After
!   convergence gives the other closest x which gives an F of opposite
!   sign from that given by x1.
! XL     Value of XX from the previous iteration.
! XLMXB  Value of XXMXL from the previous iteration.
! XO     Value of x on opposite side of the zero from the current x.
! XRND   Best accuracy that one could hope for based on the finite
!   precision of floating point numbers.
! XX     Current x, the last value of X1 where F was computed.
! XXMXL  XX - XL
! XXMXO  XX - XO = length of interval in which 0 lies.
! XXMXOL Value of XXMXO from a previous iteration.
!
  parameter (LINIT = -40)
  integer KNKP(2), LCHG(30), LMODE, LNLP(2), NP
  real             XX, XO, XL, FF, FO, FL, FFDFO
  real             DIV, QFM, QXM, TP, TP1, XXMXO, XXMXOL
  real             RND, XRND, SMALL, TOLX
  real             XXMXL, XLMXB, FFMFL, FFMFB, FLMFB
  real             DXDFFF, DXDFFO, DFDXXX, DFDXXO
  real             C0, C1, C2, C4, CP125, CP25, CP5, CP75, C1P25
  real             C8, CP01, CP99, C1P01, CP001, C1P031
  real             R1MACH
!
  parameter (C0 = 0.E0, C1 = 1.E0, C2 = 2.E0, C4 = 4.E0)
  parameter (C8 = 8.E0)
  parameter (CP125 = 0.125E0, CP25 = 0.25E0, CP75 = 0.75E0)
  parameter (CP5 = 0.5E0)
  parameter (C1P25 = 1.25E0)
  parameter (CP01 = 0.01E0)
  parameter (CP001 = 0.001E0)
  parameter (CP99 = 0.99E0)
  parameter (C1P01 = 1.01E0)
  parameter (C1P031 = 1.03125E0)
!
!                      Declarations for error message processing.
!
  real             FDAT(4)
  integer MACT(5), MLOC(4), IDAT(2)
  equivalence (INDIC, MACT(3))
  save DIV, FL, FLMFB, FO, KNKP, KS, KTYP, LCHG, LMODE, &
     LNLP, MACT, NP, RND, SMALL, XL, XLMXB, XO, XX, XXMXOL
  parameter (MERET  =51)
  parameter (MEEMES =52)
  parameter (METEXT =53)
!
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA SZERO$B
!AB Best bound for zero is [$F, $F], but tolerance is $F.$E
!AC Apparent discontinuity in function near X = $F.$E
!AD Can not find a sign change: X1=$F, X2=$F, F1=$F, F2=$F$E
!AE Called with MODE = $I.$E
!   $
!AF In SZERO -- X1=$F F1=$(E11.4) KTYP=$I DIV=$G KS=$I$E
!   $
!AG             X2=$F F2=$G$E
  integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE,LTXTAF,LTXTAG
  parameter (LTXTAA=  1,LTXTAB=  8,LTXTAC= 63,LTXTAD=112,LTXTAE=168, &
   LTXTAF=  1,LTXTAG=  1)
!      character MTXTAA(1) * (191)
!      character MTXTAB(1) * (52)
!      character MTXTAC(1) * (25)
!      data MTXTAA/'SZERO$BBest bound for zero is [$F, $F], but tolerance
!     * is $F.$EApparent discontinuity in function near X = $F.$ECan not
!     *find a sign change: X1=$F, X2=$F, F1=$F, F2=$F$ECalled with MODE =
!     * $I.$E'/
!      data MTXTAB/'In SZERO -- X1=$F F1=$(E11.4) KTYP=$I DIV=$G KS=$I$E'
!     * /
!      data MTXTAC/'            X2=$F F2=$G$E'/
!
!                      1  2  3  4      5
  data MACT / MEEMES, 0, 0, 0, MERET /
!      data MACT1 / METEXT, MERET /
  data MLOC / LTXTAB, LTXTAC, LTXTAD, LTXTAE /
!
  data RND / C0 /
  data KS, KTYP, LMODE, DIV / 0, 0, 2, C0 /
  data LCHG / &
     1, 2, 0, 0, 0, &
     1, 1, 4, 5, 0, &
     1, 2, 3, 3, 0, &
     1, 4, 4, 3, 0, &
     1, 4, 5, 5, 0, &
     1, 5, 5, 5, 0 /
  data NP / 0 /

!
! INITIALIZE
!
  if (MODE .lt. 0) then
     NP = -1 - MODE
     return
  ENDIF
  if (NP .gt. 0) then
     NP = NP - 1
     FDAT(1) = X1
     FDAT(2) = F1
     FDAT(3) = DIV
     IDAT(1) = KTYP
     IDAT(2) = KS
     GO TO 1000
!         call SMESS(MACT1, MTXTAB, IDAT, FDAT)
!         if (MODE .ne. 0) if (LMODE - 1) 70, 80, 450
!         FDAT(1) = X2
!         FDAT(2) = F2
!         call SMESS(MACT1, MTXTAC, IDAT, FDAT)
  else if (MODE .ne. 0) then
     if (LMODE - 1) 70, 80, 450
  ENDIF
!
  if (RND .eq. C0) then
     RND = R1MACH(4)
     SMALL = CP5 / (RND * R1MACH(2))
  ENDIF
  XL = X2
  FL = F2
30 TP = C1
  MODE = 1
  INDIC = 2
  XXMXOL = C0
  KNKP(1) = 0
  KNKP(2) = 0
  LNLP(1) = 0
  LNLP(2) = 0
  KS = -1
  XX = X1
  FF = F1
  if (FL) 40, 75, 50
40 if (FF) 60, 230, 100
50 if (FF) 100, 230, 60
60 LMODE = 0
!             Take care of points on same side of zero.
70 FF = F1
  XX = X1
  TP = FF / FL
  if (TP .lt. C0) go to 30
  LMODE = LMODE - 1
  if (LMODE .lt. LINIT) then
     INDIC = 5
     FDAT(1) = XX
     FDAT(2) = XL
     FDAT(3) = FF
     FDAT(4) = FL
     go to 250
  ENDIF
  if (TP .gt. C1) then
     FF = FL
     XX = XL
     FL = F1
     XL = X1
  ENDIF
  if (abs(FF) .ge. C8 * abs(FL-FF)) then
     TP = C8
  else
     TP = max(-CP25*real(LMODE), FF / (FL - FF))
  ENDIF
  FL = FF
  XO = XL
  XL = XX
  if (XX .eq. XO) XO = C1P031 * XX + sign(CP001, XX)
  XX = XX + TP * (XX - XO)
  X1 = XX
  MODE = 1
  return
!
75 X1 = XL
  F1 = FL
  go to 250
! END OF INITIALIZATION
!
!
! ENTRY AFTER COMPUTING F FOR THE LAST ITERATE
80 FF = F1
  TP = FF / FL
  if (TP) 90, 230, 110
90 TP = FF / FO
  KS = 0
100 FO = FL
  XO = XL
  go to 120
110 KS = KS + 1
120 J = 1
  if (FF .gt. C0) J = 2
  if (TP - C1) 150, 140, 130
130 I = 4
  if (TP .gt. C4) I = 5
  go to 160
140 I = 3
  go to 160
150 I = 2
  if (TP .lt. CP01) I = 1
  if (TP .lt. CP99) go to 170
160 KNKP(J) = KNKP(J) + 1
  go to 180
170 KNKP(J) = 0
180 XXMXO = XX - XO
  LNLP(J) = LCHG(5*LNLP(J) + I)
  if (LNLP(J) .ge. 4) then
     if (LNLP(3 - J) .ge. 4) go to 210
  ENDIF
! XXMXO GIVES THE LENGTH OF THE INTERVAL INSIDE WHICH
! THE ZERO IS KNOWN TO LIE.
  if (C2 * abs(XXMXO) .lt. abs(XXMXOL)) then
     KNKP(1) = 0
     KNKP(2) = 0
  ENDIF
  XXMXOL = XXMXO
  XRND = RND * (abs(XX) + abs(XO) + SMALL)
!
! TEST FOR CONVERGENCE
  if (TOL) 190, 200, 200
190 continue
  if (abs(FF) .le. abs(TOL)) go to 220
200 continue
  TOLX = max(TOL, XRND)
  if (abs(XXMXO) .gt. TOLX) go to 310
!
! CONVERGENCE -- PREPARE FOR FINAL EXIT
210 if ((abs(XXMXO) .gt. TOL) .and. (TOL .ne. C0)) then
     INDIC = 3
     FDAT(3) = TOL
     if (XXMXO .gt. 0) then
        FDAT(2) = XX
        FDAT(1) = XO
     else
        FDAT(1) = XX
        FDAT(2) = XO
     ENDIF
  ENDIF
! SET FINAL VALUES FOR X1,F1,X2,AND F2
220 continue
  if (abs(FF) .le. abs(FO)) go to 240
  F1 = FO
  X1 = XO
230 FO = FF
  XO = XX
240 X2 = XO
  F2 = FO
! TEST FOR DISCONTINUITY
  if ((KNKP(1) .gt. 5) .or. (KNKP(2) .gt. 5)) then
    INDIC = 4
    FDAT(1) = XX
  ENDIF
250 MODE = INDIC
  if (INDIC - 2) 420, 420, 430
! END OF CODE FOR FINAL EXIT
!
! F NOT DECREASING (OR THE FIRST ITERATE)
! PREPARE TO DIVIDE THE INTERVAL
260 TP = C1
  if (KS) 370, 280, 270
270 if (KTYP .eq. 0) go to 290
280 DIV = C2
290 continue
  DIV = max(DIV, FFDFO)
! KTYP=0 IF AND ONLY IF THE INTERVAL WAS DIVIDED (USING DIV)
! ON THE LAST ITERATION
  if (KTYP .eq. 0) DIV = DIV * (C1P25 / (C1P25 - TP))
! DIVIDE THE INTERVAL AS SPECIFIED BY DIV
300 TP1 = -XXMXO * (DIV/(DIV+C1))
  KTYP = 0
  go to 410
!
310 continue
  XXMXL = XX - XL
  FFMFL = FF - FL
  FFDFO = abs(FF / FO)
  TOLX = CP5 * TOLX
  if (TP .ge. C1) go to 260
! DIVIDE THE INTERVAL IF F HAS HAD THE SAME SIGN FOR
! FOUR OR MORE TIMES IN SUCCESSION
  if (KS - 4) 320, 340, 290
320 continue
  if (FLMFB .eq. C0) go to 340
! BEGINNING OF CODE TO DETERMINE IF INVERSE QUADRATIC
! INTERPOLATION IS TO BE USED.
  FFMFB = FFMFL + FLMFB
  if (FFMFB .eq. C0) go to 330
  QFM = C1 - (FFMFL / FLMFB) * (XLMXB / XXMXL)
  QXM = C1 - (XXMXL / XLMXB) * (FLMFB / FFMFL)
  DXDFFF = C1 + (FFMFL / FFMFB) * QFM
  DXDFFO = DXDFFF + C2 * ((FO - FF) / FFMFB) * QFM
  TP1 = XXMXL + XLMXB
  DFDXXX = C1 + (XXMXL / TP1) * QXM
  DFDXXO = DFDXXX + C2 * ((XO - XX) / TP1) * QXM
  TP1 = DXDFFF * DFDXXX
  if ((TP1 .le. CP25) .or. (TP1 .ge. C4)) go to 330
  TP1 = DXDFFO * DFDXXO
  if ((TP1 .gt. CP25) .and. (TP1 .lt. C4)) go to 380
!
! DERIVATIVES DO NOT MATCH WELL ENOUGH
330 continue
  if (KS .eq. 0) if (FFDFO - C1) 350, 370, 360
340 continue
  if ((KTYP .eq. 0) .and. (TP .ge. CP75)) go to 290
  continue
  TP = C1 - TP
  if (TP .le. FFDFO) go to 280
  FFDFO = FFDFO / TP
  DIV = CP125
  go to 290
350 continue
  DIV = CP5 * max(max(CP25, FFDFO), TP / (C1P25 - min(TP, C1)))
  go to 300
360 continue
  DIV = min(C4, CP5 * FFDFO)
  go to 300
! INTERPOLATE WITH SECANT METHOD
370 TP1 = -XXMXL
  go to 390
!
! DERIVATIVES MATCH UP PRETTY WELL.
380 continue
! INTERPOLATE USING THE INVERSE QUADRATIC
  TP1 = XXMXL * (QFM * (FL / FFMFB) - C1)
390 TP1 = (FF/FFMFL) * TP1
  KTYP = 1
!
! EXIT TO GET F(X)
410 continue
  FL = FF
  FLMFB = FFMFL
  XLMXB = XXMXL
  XL = XX
! COMPUTE X1, INSURING THAT IT IS NOT TOO CLOSE TO THE
! ENDS OF THE INTERVAL
  XX = min(max(XL + TP1, min(XL, XO) + TOLX), max(XL, XO) - TOLX)
  X1 = XX
420 LMODE = MODE
  return
!
430 MACT(2) = 11*INDIC  - 9
440 MACT(4) = MLOC(INDIC-2)
  GO TO 1000
!     call SMESS(MACT, MTXTAA, IDAT, FDAT)
!      go to 420
!
! A CALL TO THE SUBROUTINE HAS BEEN MADE WITH MODE.NE.1
450 IDAT(1) = MODE
  INDIC = 6
  MODE = 6
  if (LMODE .ne. 6) go to 430
  MACT(2) = 99
  go to 440
!
1000 MODE=6
  RETURN

!$$   END OF -ZERO SUBROUTINE
  end
  subroutine AMACH(MODE, I, I1, R1, D1)
!>> 1992-04-07 AMACH  Oken    Removed ^Z at EOF (error found by VAX compile)
!>> 1992-02-20 AMACH  Snyder  Added Cray-YMP stuff, q.v.
!>> 1990-06-11 AMACH  Snyder  Added Apollo DN-10000 stuff, q.v.
!>> 1990-12-14 AMACH  Lawson  Changed to eliminate ENTRY statements.
!>> 1990-08-21 AMACH  Krogh   No test was getting done for bad machine.
!>> 1990-02-28 AMACH  Krogh   Correct missing DOUBLE PRECISION AMSUB1
!>> 1989-08-14 AMACH  Krogh   Parameterized everything -- Massive change
!>> 1989-03-30 AMACH  Snyder  Correct missing "/" line 921
!>> 1989-01-30 AMACH  Snyder  Incorporate more constants from NETLIB.
!>> 1988-05-19 AMACH  Lawson  Initial code.
! File AMACH.FOR contains user-callable functions I1MACH, D1MACH, and
! R1MACH, plus second-level subroutines AMACH, AMTEST, and AMSUB1.
! Appropriate lines must be switched between comment and non-comment
! status when this code is moved to a different computer system.
!     These changes can be done with any text editor, however the "c++"
! lines permit automation of the change using the MARVEL processor.
! Note that when the MARVEL processor activates a line it shifts
! Columns 2-72 to 1-71 and puts a blank in Column 72.  When it inactiv-
! ates a line it shifts Columns 1-71 to 2-72 and puts a C in Column 1.
!     The possible choices using MARVEL are:
!              c++ SET SYS = IEEE
!              c++ SET SYS = AMDAHL
!              c++ SET SYS = APOLLO-10000
!              c++ SET SYS = BUR1700
!              c++ SET SYS = BUR5700
!              c++ SET SYS = BUR67-7700
!              c++ SET SYS = CDC60-7000
!              c++ SET SYS = CONVEXC-1
!              c++ SET SYS = CRAY1
!              c++ SET SYS = CRAY1-SD (Sngl prec.arith. used for dble.)
!              c++ SET SYS = CRAY1-64 (64 bit integers)
!              c++ SET SYS = CRAY1-SD-64 (64 bit int, SP used for DP)
!              c++ SET SYS = CRAY-YMP
!              c++ SET SYS = DG-S2000
!              c++ SET SYS = HARRIS220
!              c++ SET SYS = HON600-6000
!              c++ SET SYS = HON-DPS-8-70
!              c++ SET SYS = IBM360-370
!              c++ SET SYS = INTERDATA-8-32
!              c++ SET SYS = PDP10-KA
!              c++ SET SYS = PDP10-KB
!              c++ SET SYS = PDP11
!              c++ SET SYS = PRIME50
!              c++ SET SYS = SEQ-BAL-8000
!              c++ SET SYS = UNIVAC
!              c++ SET SYS = VAX
!     The current choice is:
!++ SET SYS = IEEE
!
!  I/O UNIT NUMBERS:
!
!    IM1 = I1MACH( 1) = THE STANDARD INPUT UNIT.
!    IM2 = I1MACH( 2) = THE STANDARD OUTPUT UNIT.
!    IM3 = I1MACH( 3) = THE STANDARD PUNCH UNIT.
!    IM4 = I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
!
!  WORDS:
!
!    IM5 = I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
!    IM6 = I1MACH( 6) = THE NUMBER OF CHARACTERS/INTEGER STORAGE UNIT.
!
!  INTEGERS:
!
!    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
!
!               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
!
!    IM7 = I1MACH( 7) = A, THE BASE.
!    IM8 = I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
!    IM9 = I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
!
!  FLOATING-POINT NUMBERS:
!
!    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
!    BASE-B FORM
!
!               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
!               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
!
!    IM10 = I1MACH(10) = B, THE BASE.
!
!  SINGLE-PRECISION:
!
!    IM11 = I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
!    IM12 = I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
!    IM13 = I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
!
!  DOUBLE-PRECISION:
!
!    IM14 = I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
!    IM15 = I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
!    IM16 = I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
!
!  CONVERSION FROM FUNCTIONAL TO STRUCTURAL FLOATING POINT CONSTANTS
!
!    IM17 = CONSTANT SUCH THAT IM14 + IM17 = ACTUAL NUMBER OF BASE-B
!           DIGITS IN DOUBLE PRECISION, USED FOR CHECKING THAT CORRECT
!           VERSION OF THIS PROGRAM IS INSTALLED.  (SEE DEFINITION OF
!           DM6, AND THE USE OF DM6 IN CALLING AMTEST.)
!
!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
!  THE DESIRED SET OF PARAMETER STATEMENTS SHOULD BE ACTIVATED BY
!  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
!  IM1 - IM4 SHOULD BE CHECKED FOR CONSISTENCY
!  WITH THE LOCAL OPERATING SYSTEM.
!     -----------------------------------------------------------------
!     Original design and code due to P. A. Fox, A. D. Hall, and
!     N. L. Schryer, Bell Laboratories.  See ACM TOMS, 4,(1978),177-188.
!     Adapted to Univac 1100 by Kris Stewart, JPL, 7/30/81.
!     Adapted for the JPL MATH77 library by C. L. Lawson and F. T. Krogh
!     Sept, 1987.
!     1989-08-14 AMACH  Krogh   Parameterized everything. Major changes.
!     1990 Dec. CLL reorganized code to avoid using ENTRY statements
!     for functions of different types.  Also added save statements.
!     -----------------------------------------------------------------
!     On the first call to this function, tests are done to verify that
!     IM10 and IM14 are not grossly wrong for the host environment.
!     This gives some protection against using the wrong version of this
!     subprogram.
!     -----------------------------------------------------------------
  integer MODE, I, I1
  real R1
  double precision D1, TEST
!
  integer IMACH(16)
  integer IM1, IM2, IM3, IM4, IM5, IM6, IM7, IM8, IM9, IM10, IM11, &
     IM12, IM13, IM14, IM15, IM16, IM17
  real             RMACH(5), RM1, RM2, RM3, RM4, RM5, &
                   RMA, RMB, RBASE
  double precision DMACH(5), DM1, DM2, DM3, DM4, DM5, DM6, &
                   DMA, DMB, DBASE
  save TEST, IMACH, RMACH, DMACH
!     -----------------------------------------------------------------
!     Machine constants for IEEE standard binary floating-point
!     processors.  This includes PC's and work-stations using the
!     Intel 8087, 80287, 80387, ... processors or the
!     Motorola 68881, 68882, ... processors.
!     Note:  We are setting the "most negative exponent" (IMACH(12) and
!     IMACH(15)) to be the exponent of the smallest normalized number.
!     An IEEE processor actually handles smaller numbers before
!     underflowing, however these "unnormalized" numbers have
!     diminished precision.
!
!++ Code for SYS = IEEE is ACTIVE
   PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
   PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
   PARAMETER (IM10 =2, IM11 =24, IM12 =-125, IM13 =128)
   PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17=0)
!     -----------------------------------------------------------------
!++ Code for SYS = AMDAHL is INACTIVE
!C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
!      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
!      -----------------------------------------------------------------
!++ Code for SYS = APOLLO-10000 is INACTIVE
!c     MACHINE CONSTANTS FOR APOLLO DN-10000 MACHINES.
!c     The only difference from IEEE is IM13.  This difference has
!c     nothing to do with the arithmetic or representation used by the
!c     machine.  It is caused by a bug in the compiler:  The right-hand
!c     side of RM2 (below) is apparently evaluated in double precision.
!c     When the compiler is ready to store the resulting value into its
!c     internal data structures, it compares it to an incorrect value
!c     of the overflow limit.  It appears the incorrect value has the
!c     correct exponent, but the fraction is 1.5 instead of 2-2**(-p),
!c     where p is the precision in bits.  You can get the correct result
!c     by changing IM13 to 128, changing RM2 from a parameter to a
!c     variable, and changing the parameter statement that assigns a
!c     value to RM2 into an ordinary assignment statement.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-125, IM13 =127)
!      PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17 =0)
!C     -----------------------------------------------------------------
!++ Code for SYS = BUR1700 is INACTIVE
!C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
!C
!      PARAMETER (IM1 =7, IM2 =2, IM3 =2, IM4 =2)
!      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =33)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-256, IM13 =255)
!      PARAMETER (IM14 =60, IM15 =-256, IM16 =255, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = BUR5700 is INACTIVE
!C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =48, IM6 =6, IM7 =2, IM8 =39)
!      PARAMETER (IM10 =8, IM11 =13, IM12 =-50, IM13 =76)
!      PARAMETER (IM14 =26, IM15 =-50, IM16 =76, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = BUR67-7700 is INACTIVE
!C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =48, IM6 =6, IM7 =2, IM8 =39)
!      PARAMETER (IM10 =8, IM11 =13, IM12 =-50, IM13 =76)
!      PARAMETER (IM14 =26, IM15 =-32754, IM16 =32780, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = CDC60-7000 is INACTIVE
!C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =60, IM6 =10, IM7 =2, IM8 =48)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-929, IM13 =1070)
!      PARAMETER (IM14 =94, IM15 =-929, IM16 =1069, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = CONVEXC-1 is INACTIVE
!C     MACHINE CONSTANTS FOR CONVEX C-1.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-128, IM13 =127)
!      PARAMETER (IM14 =53, IM15 =-1024, IM16 =1023, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY1 is INACTIVE
!C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
!      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
!      PARAMETER (IM14 =94, IM15 =-8099, IM16 =8190, IM17=2)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY-YMP is INACTIVE
!C     MACHINE CONSTANTS FOR THE CRAY YMP
!C     Cray claims the overflow exponent (IM13 and IM16) is 8189, and
!C     the underflow exponent (IM12 and IM15) is -8189, but these values
!C     don't seem to work in cf77:  the underflow limit underflows, and
!C     the overflow limit overflows when using Cray's values.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
!      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-8188, IM13 =8189)
!      PARAMETER (IM14 =94, IM15 =-8188, IM16 =8189, IM17=2)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY1-SD is INACTIVE
!C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3, WHEN DOUBLE
!C     PRECISION IS TO USE SINGLE PRECISION ARITHMETIC.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
!      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
!      PARAMETER (IM14 =47, IM15 =-8189, IM16 =8190, IM17=1)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY1-64 is INACTIVE
!C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
!      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =63)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
!      PARAMETER (IM14 =94, IM15 =-8099, IM16 =8190, IM17=2)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY1-SD-64 is INACTIVE
!C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3, WHEN DOUBLE
!C     PRECISION IS TO USE SINGLE PRECISION ARITHMETIC.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
!      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =63)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
!      PARAMETER (IM14 =47, IM15 =-8189, IM16 =8190, IM17=1)
!C     -----------------------------------------------------------------
!++ Code for SYS = DG-S2000 is INACTIVE
!C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!C
!      PARAMETER (IM1 =11, IM2 =12, IM3 =8, IM4 =10)
!      PARAMETER (IM5 =16, IM6 =2, IM7 =2, IM8 =15)
!      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
!      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = HARRIS220 is INACTIVE
!C     MACHINE CONSTANTS FOR THE HARRIS 220, SLASH 6, SLASH 7.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =0, IM4 =6)
!      PARAMETER (IM5 =24, IM6 =3, IM7 =2, IM8 =23)
!      PARAMETER (IM10 =2, IM11 =23, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =38, IM15 =-127, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = HON600-6000 is INACTIVE
!C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =43, IM4 =6)
!      PARAMETER (IM5 =36, IM6 =6, IM7 =2, IM8 =35)
!      PARAMETER (IM10 =2, IM11 =27, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =63, IM15 =-127, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = HON-DPS-8-70 is INACTIVE
!C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =43, IM4 =6)
!      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =35)
!      PARAMETER (IM10 =2, IM11 =27, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =63, IM15 =-127, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = IBM360-370 is INACTIVE
!C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
!
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
!      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = INTERDATA-8-32 is INACTIVE
!C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
!C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =6, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =62)
!      PARAMETER (IM14 =14, IM15 =-64, IM16 =62, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = PDP10-KA is INACTIVE
!C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =36, IM6 =5, IM7 =2, IM8 =35)
!      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
!      PARAMETER (IM14 =54, IM15 =-101, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = PDP10-KB is INACTIVE
!C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =36, IM6 =5, IM7 =2, IM8 =35)
!      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
!      PARAMETER (IM14 =62, IM15 =-128, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = PDP11 is INACTIVE
!C     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING
!C     16-BIT INTEGER ARITHMETIC.
!
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =16, IM6 =2, IM7 =2, IM8 =15)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =56, IM15 =-127, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = PRIME50 is INACTIVE
!C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
!C     WITH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
!C     SUPPLIED BY IGOR BRAY.
!
!      PARAMETER (IM1 =1, IM2 =1, IM3 =2, IM4 =1)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =23, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =47, IM15 =-32895, IM16 =32637, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = SEQ-BAL-8000 is INACTIVE
!C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
!C
!      PARAMETER (IM1 =0, IM2 =0, IM3 =7, IM4 =0)
!      PARAMETER (IM5 =32, IM6 =1, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-125, IM13 =128)
!      PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = UNIVAC is INACTIVE
!C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!C
!C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 1
!C     WHICH IS APPROPRIATE FOR THE UNIVAC-FTN SYSTEM.
!C     IF YOU HAVE THE UNIVAC-FOR SYSTEM, SET IT TO 7.
!C     IM6 = 4 for FTN (4 chars per word), 6 for FOR (6 chars per word).
!c
!      PARAMETER (IM1 =5, IM2 =6, IM3 =1, IM4 =6)
!      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =35)
!      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
!      PARAMETER (IM14 =60, IM15 =-1024, IM16 =1023, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = VAX is INACTIVE
!     MACHINE CONSTANTS for the VAX/VMS
!     and for PDP-11 FORTRAN SUPPORTING 32-BIT INTEGER ARITHMETIC.
!
!     PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!     PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!     PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
!     PARAMETER (IM14 =56, IM15 =-127, IM16 =127, IM17=0)
!++ end
!     -----------------------------------------------------------------
!
  PARAMETER (IM9 = 2 * (2**(IM8-1) - 1) + 1)
!
! Real parameters
!
!  RM1 = R1MACH(1) = B**(EMIN-1), The smallest positive number, i.e.,
!                    the underflow limit.
!  RM2 = R1MACH(2) = B**EMAX*(1 - B**(-T)), The largest number, i.e.,
!                    the overflow limit.
!  RM3 = R1MACH(3) = B**(-T), The smallest relative spacing, i.e., the
!                    difference between 1.0 and the next smaller number.
!  RM4 = R1MACH(4) = B**(1-T), The largest relative spacing, i.e., the
!                     difference between 1.0 and the next larger number.
!  RM5 = R1MACH(5) = LOG10(B).  When B = 2 this value is
!              Log10(2) = 0.30102_99956_63981_19521_37388_94724
!
! Parameter RMA and RMB are selected so that for values of the base =
! 2, 8, 16, 10, RMA has the values 1, 3, 4, 0, and RMB has the values 0,
! 0, 0, 1.  These values are used in computing RM5.
! $$$$ Note that if other bases are to be supported, the calculation of
! $$$$ RMA and RMB will have to be generalized.
!
  PARAMETER (RMA = ((IM10 - 10) * (-3 + ((IM10 - 2) * (-77 + &
      12 * (IM10 - 8))) / 14)) / 24)
  PARAMETER (RMB = ((IM10 - 2) * (IM10 - 8) * (16 - IM10)) / 96)
  PARAMETER (RBASE = IM10)
!
!     Weird subterfuges below are NECESSARY to compute DM1 and DM2 on
!     some systems.  DON'T SIMPLIFY THEM.  We compute RM1 and RM2 using
!     these subterfuges so it will be clear we're computing the REAL and
!     DOUBLE PRECISION characteristics in the same way.
  PARAMETER (RM1 = (RBASE**(IM12/2)) * (RBASE**(IM12-IM12/2-1)))
  PARAMETER (RM2 = RBASE**(IM13-IM11) * ((RBASE**IM11 - RBASE) &
                 + (RBASE - 1.0E0)))
  PARAMETER (RM3 = RBASE**(-IM11))
  PARAMETER (RM4 = RBASE**(1-IM11))
!     PARAMETER (RM5 = RMA*0.30102 99956 63981 19521 37388 94724E0+RMB)
  PARAMETER (RM5 = RMA*0.301029995663981195213738894724E0+RMB)
!
! Double precision paramters -- (Defined like the real ones.)
!
  PARAMETER (DMA = ((IM10 - 10) * (-3 + ((IM10 - 2) * (-77 + &
      12 * (IM10 - 8))) / 14)) / 24)
  PARAMETER (DMB = ((IM10 - 2) * (IM10 - 8) * (16 - IM10)) / 96)
  PARAMETER (DBASE = IM10)
!
!     Weird subterfuges below are NECESSARY to compute DM1 and DM2 on
!     some systems.  DON'T SIMPLIFY THEM.
  PARAMETER (DM1 = (DBASE**(IM15/2)) * (DBASE**(IM15-IM15/2-1)))
  PARAMETER (DM2 = DBASE**(IM16-IM14) * ((DBASE**IM14 - DBASE) &
                 + (DBASE - 1.0D0)))
  PARAMETER (DM3 = DBASE**(-IM14))
  PARAMETER (DM4 = DBASE**(1-IM14))
!     PARAMETER (DM5 = DMA*0.30102 99956 63981 19521 37388 94724D0+DMB)
  PARAMETER (DM5 = DMA*0.301029995663981195213738894724D0+DMB)
! DM6 and TEST are used in checking that the correct constants have been
! selected.
  PARAMETER (DM6 = DBASE**(-IM14-IM17))
  data TEST / 0.D0 /
!
!     DATA IMACH / IM1, IM2, IM3, IM4, IM5, IM6, IM7, IM8, IM9, IM10,
!    1   IM11, IM12, IM13, IM14, IM15, IM16 /
!     DATA RMACH / RM1, RM2, RM3, RM4, RM5 /
!     DATA DMACH / DM1, DM2, DM3, DM4, DM5 /
!     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (TEST .eq. 0.0D0) then
     IMACH(1) = IM1
     IMACH(2) = IM2
     IMACH(3) = IM3
     IMACH(4) = IM4
     IMACH(5) = IM5
     IMACH(6) = IM6
     IMACH(7) = IM7
     IMACH(8) = IM8
     IMACH(9) = IM9
     IMACH(10) = IM10
     IMACH(11) = IM11
     IMACH(12) = IM12
     IMACH(13) = IM13
     IMACH(14) = IM14
     IMACH(15) = IM15
     IMACH(16) = IM16
     RMACH(1) = RM1
     RMACH(2) = RM2
     RMACH(3) = RM3
     RMACH(4) = RM4
     RMACH(5) = RM5
     DMACH(1) = DM1
     DMACH(2) = DM2
     DMACH(3) = DM3
     DMACH(4) = DM4
     DMACH(5) = DM5
     CALL AMTEST (TEST, DM6)
  ENDIF

  if (MODE .eq. 0) then
     I1=IMACH(I)
  else if (MODE .eq. 1) then
     R1=RMACH(I)
!                                  Here we assume MODE = 2.
  else
     D1=DMACH(I)
  endif
  return
  end
!     ==================================================================
  integer function I1MACH(I)
  integer I, I1
  real R1
  double precision D1
!      IF (I .LT. 1  .OR.  I .GT. 16) THEN
!         PRINT*,'I1MACH.. Bad argument: I =',I
!         STOP 'I1MACH error'
!      ENDIF
  call AMACH (0, I, I1, R1, D1)
  I1MACH = I1
  return
  end
!     ==================================================================
!
  real function R1MACH(I)
  integer I, I1
  real R1
  double precision D1
!      IF (I .lt. 1  .or.  I .gt. 5) THEN
!         print*,'R1MACH.. Bad argument: I = ',I
!         stop 'R1MACH error'
!      ENDIF
  call AMACH (1, I, I1, R1, D1)
  R1MACH = R1
  RETURN
  end
!     ==================================================================
!
  double precision function D1MACH(I)
  integer I, I1
  real R1
  double precision D1
!      IF (I .lt. 1  .or.  I .gt. 5) THEN
!         print*,'D1MACH.. Bad argument: I = ',I
!         stop 'D1MACH error'
!      ENDIF
  call AMACH (2, I, I1, R1, D1)
  D1MACH = D1
  RETURN
  END
!     ==================================================================
!
  SUBROUTINE AMTEST (TEST, D6)
! Verifies that D6 is an appropriate values for DM6.
  DOUBLE PRECISION AMSUB1, D6, TEST
  TEST = AMSUB1(1.D0 + D6)
!
! The comparison with 1.875E0*D6 in the line below is to guard
! against the possibility that TEST is > 0 as a result of rounding
! up in the addition of D6 to 1.
!
  IF ((TEST .eq. 0.D0) .or. (TEST .gt. 1.875D0*D6)) THEN
     TEST = (D6 + D6) + 1.D0
     IF (AMSUB1(TEST) .ne. 0.D0) RETURN
  ENDIF
!      print*,'AMACH has bad parameters for current environment.'
!      stop
  RETURN
  END
!     ==================================================================
!
  DOUBLE PRECISION FUNCTION AMSUB1 (TEST1)
  DOUBLE PRECISION TEST1
!     Returns the value of TEST1 - 1.
  AMSUB1 = TEST1 - 1.0D0
  RETURN
  END
