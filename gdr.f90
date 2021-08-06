subroutine gdr(CX,CY,CZ,N,NN2,BOXL,DENS)
    !****************************************************
    !
    ! Esta subrutina calcula la funcion de destribucion
    ! g(r), para que a partir de esta se puedan calcular
    ! magnitudes termodinamicas.
    !
    ! Autores: Jonas Valenzuela
    ! Agradecimientos: Dra. Laura Laura Lorenia Yeomans Reyna
    !
    !****************************************************

    !**************************
    ! Declaracion de variables
    !**************************
    implicit none

    integer, intent(in) :: N                           !Variables de entrada
    integer, intent(in) :: NN2
    real, intent(in) :: BOXL, DENS
    real, dimension(N,NN2), intent(in) :: CX, CY, CZ

    integer :: capas                                       !Variables auxiliares
    real :: rij, RX, RY, RZ, C1, C2
    real, parameter :: dr = 0.025
    real, parameter :: pi = 3.14159265
    integer :: i, j, k, l
    character(24) :: nombre

    integer, dimension(:), allocatable :: Hist                   !Variables de salida
    real, dimension(:), allocatable :: G


    !***************************
    ! Calculo de histograma
    !***************************
    capas = int(BOXL/(2.0*dr))
    C1 = 4.0*pi*DENS/3.0
    allocate(Hist(50000))
    allocate(G(capas))
    Hist=0

    do k=1, NN2
        do i=1, N
            do j=1, N
                if (i /= j) then
                    RX = CX(i,k)-CX(j,k)
                    RY = CY(i,k)-CY(j,k)
                    RZ = CZ(i,k)-CZ(j,k)

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
        G(l) = G(l)/(N*NN2*C2)
    end do


    !**************************
    ! Exportamos resultados
    !**************************
    write(nombre, '(a,f4.2,a)') 'gdr',DENS,'.txt'
    open(4,file=nombre,status="unknown")

    do l=1, capas
        write(4,*) dr*(l-1.0/2.0), G(l)
        !write(4,*) dr*(l-1.0/2.0), G(l)*2.49 !En angstrom, se quita reduccion con diametro ref
    end do
    close(4)

    !write(nombre, '(a,i5,a)') 'gconact',int(DENS*1000),'.txt'
    !open(68,file=nombre,status="unknown")
    !l=0
    !do
    !    l=l+1
    !    if (G(l) > 0) then
    !        write(68,*) DENS, G(l)
    !        exit
    !    end if
    !end do
    !close(68)


    deallocate(G,Hist)

end subroutine
