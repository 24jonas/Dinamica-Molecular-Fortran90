subroutine correctefrain(NP,X,Y,Z,XR,YR,ZR,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,X5,Y5,Z5,FX,FY,FZ,AL,dt,AHEAT,etot)
    !****************************************************************
    !
    ! Este programa corrige la configuración actual, utilizando
    ! las fuerzas medidas de acuerdo con el metodo de Gear
    ! predictor-corrector.
    !
    ! Unidades: Reducidas Sigma-Beta
    !
    ! Autores: Jonas Valenzuela y Dr. Efrain Urrutia Bañuelos
    !
    !****************************************************************

      integer :: NP
      real, dimension(NP) :: X,Y,Z,XR,YR,ZR,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,X5,Y5,Z5,FX,FY,FZ
      real :: AL, dt, AHEAT
      real :: etot

      real :: ALFA0, ALFA1, ALFA3, ALFA4, ALFA5, STSQH, vlcosq
      real :: XERROR, YERROR, ZERROR
      integer :: I8






      vlcosq=0.0
      ALFA0=3.0/16.0
      ALFA1=251.0/360.0
      ALFA3=11.0/18.0
      ALFA4=1.0/6.0
      ALFA5=1.0/60.0
      STSQH=dt*dt*0.50
      DO 500 I8=1,NP
        XERROR=STSQH*FX(I8)-X2(I8)
        YERROR=STSQH*FY(I8)-Y2(I8)
        ZERROR=STSQH*FZ(I8)-Z2(I8)
        X(I8)=X(I8)+XERROR*ALFA0
        XR(I8)=XR(I8)+XERROR*ALFA0
        X1(I8)=X1(I8)+XERROR*ALFA1
        X2(I8)=X2(I8)+XERROR
        X3(I8)=X3(I8)+XERROR*ALFA3
        X4(I8)=X4(I8)+XERROR*ALFA4
        X5(I8)=X5(I8)+XERROR*ALFA5
        Y(I8)=Y(I8)+YERROR*ALFA0
        YR(I8)=YR(I8)+YERROR*ALFA0
        Y1(I8)=Y1(I8)+YERROR*ALFA1
        Y2(I8)=Y2(I8)+YERROR
        Y3(I8)=Y3(I8)+YERROR*ALFA3
        Y4(I8)=Y4(I8)+YERROR*ALFA4
        Y5(I8)=Y5(I8)+YERROR*ALFA5
        Z(I8)=Z(I8)+ZERROR*ALFA0
        ZR(I8)=ZR(I8)+ZERROR*ALFA0
        Z1(I8)=Z1(I8)+ZERROR*ALFA1
        Z2(I8)=Z2(I8)+ZERROR
        Z3(I8)=Z3(I8)+ZERROR*ALFA3
        Z4(I8)=Z4(I8)+ZERROR*ALFA4
        Z5(I8)=Z5(I8)+ZERROR*ALFA5
!         CONDICIONES PERIODICAS
        X(I8)=X(I8)-AL*NINT(X(I8)/AL)
        Y(I8)=Y(I8)-AL*NINT(Y(I8)/AL)
        Z(I8)=Z(I8)-AL*NINT(Z(I8)/AL)
! SUMA DE VELOCIDADES PARA ESCALARLAS DURANTE LA EQUILIBRACION
        VLCOSQ=VLCOSQ+X1(I8)**2+Y1(I8)**2+Z1(I8)**2
 500  CONTINUE
        !write(*,*) VLCOSQ, "     Potencial abajo"
       FACTOR=SQRT(AHEAT/VLCOSQ)
       !if(l1.eq.nc)then
       !factor=1.00
       !endif
        do ij8=1,np
        !X1R(IJ8)=X1(IJ8)
        !Y1R(IJ8)=Y1(IJ8)
        !Z1R(IJ8)=Z1(IJ8)
        X1(Ij8)=X1(Ij8)*FACTOR
        Y1(Ij8)=Y1(Ij8)*FACTOR
        Z1(Ij8)=Z1(Ij8)*FACTOR
        enddo
        etot = 0.5*VLCOSQ/(dt*dt*NP) !Suma de velocidades (unidad pos)
      RETURN


end subroutine
