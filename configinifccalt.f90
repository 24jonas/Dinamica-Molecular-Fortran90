subroutine configinifccalt(np,X,Y,Z,al)
    !****************************************************************
    !
    ! Este programa crea una configuracion inicial aleatoria, donde
    ! los puntos estan distanciados aleatoriamente en X, Y y Z dentro
    ! de un cubo, con maria luisa.
    !
    ! Unidades: Reducidas Sigma-Beta
    !
    ! Autores: Dr. Efrain Urrutia Bañuelos y Dr. Raúl Garibay Alonso
    !
    !****************************************************************
      implicit none

      integer :: np
      real, dimension(np) :: X, Y, Z
      real :: al

      real :: nunit, ncheck, dist
      integer :: i, j, k, ij, kct, m


      NUNIT=(real(np)/4.0)**(1.0/3.0)+.10
 5    NCHECK=4*(NUNIT**3)
      IF (NCHECK.LT.NP) THEN
        NUNIT=NUNIT+1
        GOTO 5
      ENDIF
      DIST =0.50*al/NUNIT
      X(1)=0.
      Y(1)=0.
      Z(1)=0.
      X(2)=0.
      Y(2)=DIST
      Z(2)=DIST
      X(3)=DIST
      Y(3)=0.
      Z(3)=DIST
      X(4)=DIST
      Y(4)=DIST
      Z(4)=0.
      M=0
      KCT=0
      DO 12 I=1,int(NUNIT)
       DO 12 J=1,int(NUNIT)
        DO 12 K=1,int(NUNIT)
         DO 10 IJ=1,4
          IF (KCT.LT.NP) THEN
            X(IJ+M)= X(IJ)+2.0*DIST*(K-1)
            Y(IJ+M)= Y(IJ)+2.0*DIST*(J-1)
            Z(IJ+M)= Z(IJ)+2.0*DIST*(I-1)
          ENDIF
          KCT=KCT+1
 10      CONTINUE
         M=M+4
 12     CONTINUE



    open(31,file="puntosinifcc.txt", status="unknown")
    do i=1, np
        write(31,*) X(i), Y(i), Z(i)
    end do
    close(31)

end subroutine
