subroutine configinireg(BOXL,N,X,Y,Z,DENS)
    !****************************************************************
    !
    ! Este programa crea una configuracion inicial regular, donde
    ! los puntos estan distanciados uniformemente en X, Y y Z dentro
    ! de un cubo.
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
    integer, intent(in) :: N                                !Variables de entrada
    real, intent(in) :: BOXL, DENS

    integer :: i, j, k, m                                   !Variables auxiliares
    integer :: nlin
    real, dimension(:), allocatable :: Xpos, Ypos, Zpos
    character(24) :: nombre

    real, dimension(N), intent(out) :: X, Y, Z              !Variables de salida


    !**************************
    ! Calculos iniciales
    !**************************
    nlin = int(real(N)**(1.0/3.0))
    allocate(Xpos(nlin),Ypos(nlin),Zpos(nlin))


    !***************************
    ! Creacion de Configuracion
    !***************************
    call DistLineaML(BOXL,nlin,1.0,Xpos)
    call DistLineaML(BOXL,nlin,1.0,Ypos)
    call DistLineaML(BOXL,nlin,1.0,Zpos)
    m=1
    do i=1, nlin
        do j=1, nlin
            do k=1, nlin
                X(m)=Xpos(i)
                Y(m)=Ypos(j)
                Z(m)=Zpos(k)
                m=m+1
            end do
        end do
    end do
    deallocate(Xpos)
    deallocate(Ypos)
    deallocate(Zpos)

    write(nombre, '(a,i5,a)') 'puntosini',int(DENS*1000),'.txt'
    open(31,file=nombre, status="unknown")
    do i=1, N
        write(31,*) X(i), Y(i), Z(i)
    end do
    close(31)

end subroutine

