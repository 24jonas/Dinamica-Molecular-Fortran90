subroutine configiniale(BOXL,N,X,Y,Z,SEED,DENS)
    !****************************************************************
    !
    ! Este programa crea una configuracion inicial aleatoria, donde
    ! los puntos estan distanciados aleatoriamente en X, Y y Z dentro
    ! de un cubo, con maria luisa.
    !
    ! Unidades: Reducidas Sigma-Beta
    !
    ! Autores: Jonas Valenzuela
    ! Agradecimientos: Dra. Laura Laura Lorenia Yeomans Reyna
    !
    !****************************************************************

    !***************************
    ! Declaracion de variables
    !***************************
    implicit none
    integer, intent(in) :: N, SEED                          !Variables de entrada
    real, intent(in) :: BOXL, DENS

    integer :: i, j                                         !Variables auxiliares
	real :: RMAX, Rij
	character(24) :: nombre

    real, dimension(N), intent(out) :: X, Y, Z              !Variables de salida


    !***************************
    ! Creacion de Configuracion
    !***************************
	RMAX = (BOXL-1.0)/2.0 !Maria luisa
	call srand(SEED)

	X(1) = (2*rand()-1.0)*RMAX
	Y(1) = (2*rand()-1.0)*RMAX
	Z(1) = (2*rand()-1.0)*RMAX

	do i=2, N
		!Posicion aleatoria
4		X(i) = (2*rand()-1.0)*RMAX
		Y(i) = (2*rand()-1.0)*RMAX
		Z(i) = (2*rand()-1.0)*RMAX

		!Sin traslape
		do j=i-1, 1, -1
			Rij = sqrt((X(i)-X(j))**2+(Y(i)-Y(j))**2+(Z(i)-Z(j))**2)
			if (Rij < 1.0) then
                !write(*,*) "Traslape", i, j
				go to 4
			end if
		end do
	end do


    write(nombre, '(a,i5,a)') 'puntosini',int(DENS*1000),'.txt'
    open(31,file=nombre, status="unknown")
    do i=1, N
        write(31,*) X(i), Y(i), Z(i)
    end do
    close(31)

end subroutine
