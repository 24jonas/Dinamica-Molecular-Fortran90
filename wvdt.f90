subroutine wvdt(CXR,CYR,CZR,CVX,CVY,CVZ,NN2,DT,NFREC,N)
    !*********************************************************
    !
    ! Esta subrutina se encarga de calcular el desplazamiento
    ! cuadratico medio y coeficiente de difusion dependiente
    ! del tiempo.
    !
    !*********************************************************



    !**************************
    ! Declaracion de variables
    !**************************

    implicit none

    !Variables de entrada
    integer, intent(in) :: NN2, NFREC, N
    real, intent(in) :: DT
    real, dimension(N,NN2) :: CXR, CYR, CZR, CVX, CVY, CVZ

    !Variables auxiliares
    real :: TIM, TIME, DIF
    integer :: i, j, k, l
    integer :: NTMAX
    real :: WT, WTX, WTY, WTZ, VT, VTX, VTY, VTZ

    real :: VT0, VT0M !? Para replicar Efrain


    !************************
    ! Calculos
    !************************


    open(5,file="wt.txt",status="unknown")


    !Tiempo entre configuraciones
    TIM = real(nfrec)*DT*0.1776 !en ps para Ag

    do i=1, NN2-1 !Cadencia en el barrido temporal
        NTMAX = NN2-i

        VT0 = 0
        DO l=1,N !?
            VT0=VT0+CVX(l,i)**2+CVY(l,i)**2+CVZ(l,i)**2
        END DO
        VT0M=VT0/REAL(N)

        WTX=0.0 !Desplazamiento cuadrático medio en x
        WTY=0.0 !en y
        WTZ=0.0 !en z
        WT=0.0
        VTX=0   !volume velocity en componentes
        VTY=0
        VTZ=0
        VT=0

        do j=1, N   !Barrido de particulas
            do k=1, NTMAX !Barrido de tiempo
                WTX = WTX + (CXR(j,i+k)-CXR(j,k))**2
                WTY = WTY + (CYR(j,i+k)-CYR(j,k))**2
                WTZ = WTZ + (CZR(j,i+k)-CZR(j,k))**2
                VTX = VTX + CVX(j,i+k)*CVX(j,k)
                VTY = VTY + CVY(j,i+k)*CVY(j,k)
                VTZ = VTZ + CVZ(j,i+k)*CVZ(j,k)
            end do
        end do

        TIME = TIM*i
        WT = (WTX + WTY + WTZ)/(real(NTMAX)*real(N)) !Quitar el 6 corresponde a frenkel, ponerlo corresponde a Laura
        VT = (VTX + VTY + VTZ)/(real(NTMAX)*real(N)*VT0M)

        DIF = WT/(6.0*TIME) !forma reducida pag.304 haile figura 7.13 !Y poner el 6 aqui corresponde a Efrain
        !Se multiplica por r0 para tener unidades de A^2/ps

        write(5,*) TIME, WT, DIF, VT
    end do

    close(5)

end subroutine
