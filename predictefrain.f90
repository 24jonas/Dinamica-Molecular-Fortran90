subroutine predictefrain(NP,X,Y,Z,XR,YR,ZR,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,X5,Y5,Z5,FX,FY,FZ,AL)
    !****************************************************************
    !
    ! Este programa predice la siguiente configuración, utilizando
    ! una expansion de Taylor de acuerdo con el metodo de Gear
    ! predictor-corrector.
    !
    ! Unidades: Reducidas Sigma-Beta
    !
    ! Autores: Jonas Valenzuela y Dr. Efrain Urrutia Bañuelos
    !
    !****************************************************************
    implicit none

    integer :: NP
    real :: AL
    real, dimension(NP) :: X,Y,Z,XR,YR,ZR,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,X5,Y5,Z5,FX,FY,FZ

    integer :: I7

    DO 200 I7=1,NP
        X(I7)=X(I7)+X1(I7)+X2(I7)+X3(I7)+X4(I7)+X5(I7)
        XR(I7)=XR(I7)+X1(I7)+X2(I7)+X3(I7)+X4(I7)+X5(I7)
        Y(I7)=Y(I7)+Y1(I7)+Y2(I7)+Y3(I7)+Y4(I7)+Y5(I7)
        YR(I7)=YR(I7)+X1(I7)+X2(I7)+X3(I7)+X4(I7)+X5(I7)
        Z(I7)=Z(I7)+Z1(I7)+Z2(I7)+Z3(I7)+Z4(I7)+Z5(I7)
        ZR(I7)=ZR(I7)+X1(I7)+X2(I7)+X3(I7)+X4(I7)+X5(I7)
        X1(I7)=X1(I7)+2.0*X2(I7)+3.0*X3(I7)+4.0*X4(I7)+5.0*X5(I7)
        Y1(I7)=Y1(I7)+2.0*Y2(I7)+3.0*Y3(I7)+4.0*Y4(I7)+5.0*Y5(I7)
        Z1(I7)=Z1(I7)+2.0*Z2(I7)+3.0*Z3(I7)+4.0*Z4(I7)+5.0*Z5(I7)
        X2(I7)=X2(I7)+3.0*X3(I7)+6.0*X4(I7)+10.0*X5(I7)
        Y2(I7)=Y2(I7)+3.0*Y3(I7)+6.0*Y4(I7)+10.0*Y5(I7)
        Z2(I7)=Z2(I7)+3.0*Z3(I7)+6.0*Z4(I7)+10.0*Z5(I7)
        X3(I7)=X3(I7)+4.0*X4(I7)+10.0*X5(I7)
        Y3(I7)=Y3(I7)+4.0*Y4(I7)+10.0*Y5(I7)
        Z3(I7)=Z3(I7)+4.0*Z4(I7)+10.0*Z5(I7)
        X4(I7)=X4(I7)+5.0*X5(I7)
        Y4(I7)=Y4(I7)+5.0*Y5(I7)
        Z4(I7)=Z4(I7)+5.0*Z5(I7)
!         CONDICIONES PERIODICAS
        X(I7)=X(I7)-AL*NINT(X(I7)/AL)
        Y(I7)=Y(I7)-AL*NINT(Y(I7)/AL)
        Z(I7)=Z(I7)-AL*NINT(Z(I7)/AL)
        FX(I7)=0.00
        FY(I7)=0.00
        FZ(I7)=0.00
 200  CONTINUE
      RETURN



end subroutine
