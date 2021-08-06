subroutine instgdr(N,BOXL,DENS,X,Y,Z,conf)
    !***************************************************
    !
    ! Esta subrutina se encarga de calcular la gdr
    ! de la configuracion actual, es util para calcular
    ! propiedades instantaneas que dependen de la gdr,
    ! como el exceso de entropia S2.
    !
    !***************************************************

    !**************************
    ! Declaracion de variables
    !**************************
    implicit none

    integer, intent(in) :: N                        !Variables entrada
    real, intent(in) :: BOXL, DENS
    real, dimension(N), intent(in) :: X, Y, Z
    integer :: conf

    integer :: i, j, l, capas                       !Variables auxiliares
    real, parameter :: dr = 0.025
    real, parameter :: pi = 3.14159265
    real :: RX, RY, RZ, rij, C1, C2


    integer, dimension(:), allocatable :: Hist      !Variables salida
    real, dimension(:), allocatable :: G
    real, dimension(:), allocatable :: fdx          !Para integracion entropia exceso
    real :: S2                                      !Exceso de entropia
    real :: sum2, sum3




    !***************************
    ! Calculo de histograma
    !***************************
    capas = int(BOXL/(2.0*dr))
    C1 = 4.0*pi*DENS/3.0
    allocate(Hist(50000))
    allocate(G(capas))
    allocate(fdx(capas))
    Hist=0

    do i=1, N
        do j=1, N
            if (i /= j) then
                RX = X(i)-X(j)
                RY = Y(i)-Y(j)
                RZ = Z(i)-Z(j)

                !Distancia y condicion de imagen minima
                RX = RX - BOXL*ANINT(RX/BOXL)
                RY = RY - BOXL*ANINT(RY/BOXL)
                RZ = RZ - BOXL*ANINT(RZ/BOXL)
                Rij=sqrt(RX**2+RY**2+RZ**2)

                !do l=1, capas
                !    if (Rij < dr*l) then
                !        Hist(l) = Hist(l) + 1
                !        exit
                !    end if
                !end do

                l = int(Rij/dr)+1
                Hist(l) = Hist(l) + 1
            end if
        end do
    end do


    !**************************
    ! Calculo de g(r)
    !**************************
    G=0

    do l=1, capas
        G(l) = G(l) + Hist(l)
    end do

    do l=1, capas
        C2=C1*((l*dr)**3-((l-1)*dr)**3)
        G(l) = G(l)/(N*C2)
    end do



    !*******************************
    ! Calculos finales con gdr inst
    !*******************************
    !Calcular argumento de integral / integrando
    do l=1, capas
        fdx(l) = 0
        if (G(l) /= 0.0) then
            fdx(l) = (G(l)*log(G(l)) - (G(l)-1))*((l-0.5)*dr)**2
        end if
    end do

    !Calcular integral con simpson compuesta 1/3
    sum2=0
    do l=2, capas-2, 2
        sum2 = sum2 + 2.0*fdx(l)
    end do

    sum3=0
    do l=3, capas-1, 2
        sum3 = sum3 + 4.0*fdx(l)
    end do

    S2 = -2.0*pi*DENS*(dr/3.0)*(fdx(1)+sum2+sum3+fdx(capas))

    !write(*,*) "Exceso de entropia: ", S2



    !******************************
    ! Exportar resultados
    !******************************

    open(60,file="instagdr.txt",status="unknown")

    do l=1, capas
        write(60,*) dr*(l-1.0/2.0), G(l)
    end do
    close(60)

    open(61,file="instaS2.txt",status="unknown")
    write(61,*) conf, S2
    close(61)











    deallocate(G,Hist,fdx)





end subroutine
