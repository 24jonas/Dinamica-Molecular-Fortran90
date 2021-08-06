subroutine veliniale(N,TEMP,dt,X,Y,Z,XM,YM,ZM,VX,VY,VZ)
    !************************************************
    ! Esta subrutina asigna velocidades aleatorias
    ! a cada particula y les da una escala apropiada
    ! a la temperatura del sistema
    !************************************************

    implicit none
    integer, intent(in) :: N
    real, intent(in) :: TEMP, dt
    real, dimension(N) :: X, Y, Z, XM, YM, ZM
    real, dimension(N), intent(out) :: VX, VY, VZ
    real :: sumvx, sumvy, sumvz, sumv2
    real :: fs
    integer :: i

    sumvx = 0.
    sumvy = 0.
    sumvz = 0.
    sumv2 = 0.



    do i=1, N
        VX(i) = rand() - 0.5
        VY(i) = rand() - 0.5
        VZ(i) = rand() - 0.5
        sumvx = sumvx + VX(i)
        sumvy = sumvy + VY(i)
        sumvz = sumvz + VZ(i)

    end do
    sumvx = sumvx/N
    sumvy = sumvy/N
    sumvz = sumvz/N

    do i=1, N
        VX(i) = (VX(i)-sumvx)
        VY(i) = (VY(i)-sumvy)
        VZ(i) = (VZ(i)-sumvz)
        sumv2 = sumv2 + VX(i)**2 + VY(i)**2 + VZ(i)**2
    end do

    fs = sqrt(3.0*dt*dt*real(N)*TEMP/sumv2) !Para version efrain gear (velocidad tiene unidades de posicion)

    do i=1, N
        VX(i) = VX(i)*fs
        VY(i) = VY(i)*fs
        VZ(i) = VZ(i)*fs
        XM(i) = X(i) - VX(i)*dt
        YM(i) = Y(i) - VY(i)*dt
        ZM(i) = Z(i) - VZ(i)*dt
    end do

end subroutine
