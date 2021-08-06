subroutine fuerzas(N,BOXL,X,Y,Z,FX,FY,FZ,ENPOT)
    !*******************************************
    ! Esta subrutina calcula la fuerza y energia
    ! sobre cada particula.
    !*******************************************

    implicit none
    integer, intent(in) :: N
    real, intent(in) :: BOXL
    real, dimension(N), intent(IN) :: X, Y, Z

    real :: ENPOT
    real, dimension(N), intent(out) :: FX, FY, FZ

    integer :: i, j
    real :: r, rx, ry, rz
    real :: RCUT

    !Variables LJ
    !real :: U2, ECUT

    !Variables gupta
    real :: A, XI, P, Q, U0
    real :: FPIJ, FQIJ
    real :: sexq, vq, vp
    real, dimension(N) :: SQEXQ

    do i=1, N
        FX(i) = 0
        FY(i) = 0
        FZ(i) = 0
        SQEXQ(i)= 0
    end do
    ENPOT = 0
    RCUT = BOXL/2

    !Para Lennard Jones
    !ECUT = 4*((1/RCUT**12)-(1/RCUT**6))

    !Para Gupta, metal Ni
    A=0.0376
    XI=1.07
    P=16.999
    Q=1.189
    U0=4.44

    !Para Gupta, metal Cu
    !A=0.0855
    !XI=1.224
    !P=10.96
    !Q=2.278
    !U0=3.54

    !Para Gupta, metal Pt
    !A=0.2975
    !XI=2.695
    !P=10.612
    !Q=4.004
    !U0=5.85

    !Para Gupta, metal Ag
    !A=0.1028
    !XI=1.178
    !P=10.928
    !Q=3.139
    !U0=2.96





    !Potencial gupta (se necesita correr antes de fuerza gupta para sqexp)
    do i=1, N
        sexq = 0
        vq = 0
        vp = 0
        do j=1, N
            if (i == j) then
                cycle
            end if
            rx = X(i) - X(j)
            rx = rx - BOXL*anint(rx/BOXL) !Condicion periodica
            ry = Y(i) - Y(j)
            ry = ry - BOXL*anint(ry/BOXL)
            rz = Z(i) - Z(j)
            rz = rz - BOXL*anint(rz/BOXL)

            r = sqrt(rx**2 + ry**2 + rz**2)

            if (r < RCUT) then  !Radio de corte
                sexq=sexq+exp(-2.0*q*r)         !El programa principal hace doble calculo
                vq=vq+exp(-2.0*q*(r-1.0))
                vp=vp+exp(-p*(r-1.0))
            end if

        end do
        enpot=enpot+A*vp-XI*sqrt(vq)
        SQEXQ(i)=1.0/sqrt(sexq)
    end do



    do i=1, N
        do j=i+1, N
            !if (i == j) then
            !    cycle
            !end if
            rx = X(i) - X(j)
            rx = rx - BOXL*anint(rx/BOXL) !Condicion periodica
            ry = Y(i) - Y(j)
            ry = ry - BOXL*anint(ry/BOXL)
            rz = Z(i) - Z(j)
            rz = rz - BOXL*anint(rz/BOXL)

            r = sqrt(rx**2 + ry**2 + rz**2)

            if (r < RCUT) then  !Radio de corte
                !Fuerza Lennard Jones
                !U2 = 48*(1/r**14) - 24*(1/r**8)

                !FX(i) = FX(i) + U2*rx
                !FY(i) = FY(i) + U2*ry
                !FZ(i) = FZ(i) + U2*rz

                !FX(j) = FX(j) - U2*rx
                !FY(j) = FY(j) - U2*ry
                !FZ(j) = FZ(j) - U2*rz

                !Potencial Lennard Jones
                !ENPOT = ENPOT + 4*((1/r**12)-(1/r**6)) - ECUT !Pot LJ


                !Fuerza Gupta
                FPIJ=P*A*exp(-p*(r-1.0))/U0
                FQIJ=XI*Q*EXP(-Q*(2.0*r-1.0))*SQEXQ(i)/U0

                FX(i) = FX(i) + (FPIJ-FQIJ)*rx/r
                FY(i) = FY(i) + (FPIJ-FQIJ)*ry/r
                FZ(i) = FZ(i) + (FPIJ-FQIJ)*rz/r

                FX(j) = FX(j) - (FPIJ-FQIJ)*rx/r
                FY(j) = FY(j) - (FPIJ-FQIJ)*ry/r
                FZ(j) = FZ(j) - (FPIJ-FQIJ)*rz/r

            end if

        end do
    end do

    ENPOT = ENPOT/N


end subroutine

