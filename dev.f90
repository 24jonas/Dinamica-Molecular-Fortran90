subroutine dev(col,rows,maxcol)
    !**********************************************************
    !
    ! Esta subrutina toma un archivo de nombre "name", lee la
    ! variable de la columna "col", y saca un promedio y
    ! desviación utilizando "row" datos.
    !
    !**********************************************************

    !***************************
    ! Declaracion de variables
    !***************************
    implicit none

    character(24) :: name = "energiastempconf.txt"          !Variables entrada
    integer, intent(in) :: col, rows, maxcol

    real, dimension(maxcol,rows) :: Var          !Variables auxiliar
    integer :: i, j

    real :: prom, desv                           !Variables salida
    character(24) :: newname = "PromDesvPotencial.txt"



    !**************************
    ! Lectura de datos
    !**************************

    open(62,file=name,status="old")
    do j=1, rows
        read(62,*) ( Var(i,j), i=1,maxcol)
    end do
    close(62)



    !*************************
    ! Calculos
    !*************************
    prom = 0
    do j=1, rows
        prom = prom + Var(col,j)
    end do
    prom = prom/rows

    desv = 0
    do j=1, rows
        desv = desv + (Var(col,j)-prom)**2
    end do
    desv = sqrt(desv/rows)



    !*************************
    ! Exportar resultados
    !*************************
    open(63,file=newname,status="unknown")

    write(63,*) prom
    write(63,*) desv

    close(63)


end subroutine
