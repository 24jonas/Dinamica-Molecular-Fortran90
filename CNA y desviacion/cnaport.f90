program cnaport
    !******************************************************************************
    !
    ! Programa CNA portatil. Analiza 1 configuración con el metodo de vecinos
    ! comunes. Puede utilizarse radio de corte fijo o adaptativo. El formato del
    ! archivo de entrada debe ser X, Y, Z. Se debe indicar en las variables de
    ! entrada el numero de particulas, el nombre del archivo de configuracion,
    ! si se desea utilizar radio de corte fijo o adaptativo, y concentracion
    ! reducida.
    !
    ! Para analizar los primeros vecinos, rcut debe ser la distancia del primer
    ! minimo de la gdr.
    !
    !******************************************************************************



    !Declaracion de variables
    implicit none

    !Variables de entrada
    integer, parameter :: N = 500
    character(24) :: name = "puntosfin.txt"
    logical, parameter :: RcutAdaptativo = .false.
    real, parameter :: RcutFijo = 1.3375 !1.3637
    real, parameter :: DENS = 1.312
    real, dimension(N) :: X, Y, Z

    !Variables auxiliares
    real :: BOXL
    real, dimension(N) :: dist
    integer, dimension(N) :: ind
    real :: RX, RY, RZ
    real, dimension(N-1) :: Rij
    integer, dimension(N-1) :: indaux
    integer :: i, j, k, m, l, s, r1, r2, r2a, s3, r3, r3a, r4, s4
    real :: fcc, hcp, bcc, icos, no
    real :: norm
    integer :: icna, jcna, kcna, f421, f311, f422, f221, f441, f321, f444, f666
    integer :: f555, f544, f433, f665
    integer, dimension(24) :: jind
    integer, dimension(24,2) :: indk
    real, dimension(24) :: distint
    logical :: encontre1, encontre2

    real :: rcut



    !Inicializacion

    open(62,file=name,status="old") !Leer datos
    do j=1, N
        read(62,*) X(j), Y(j), Z(j)
    end do
    close(62)

    BOXL = ((real(N)/DENS)**(1./3.))        !Longitud de el contenedor cubico

    fcc = 0
    hcp = 0
    bcc = 0
    icos = 0
    no = 0

    norm = 0

    if (RcutAdaptativo) then
        rcut = 0            !Obtenemos rcut
        do i=1, N
            m=0
            do j=1, N
                if (i /= j) then
                    m=m+1
                    RX = X(i)-X(j)
                    RY = Y(i)-Y(j)
                    RZ = Z(i)-Z(j)

                    !Distancia y condicion de imagen minima
                    RX = RX - BOXL*ANINT(RX/BOXL)
                    RY = RY - BOXL*ANINT(RY/BOXL)
                    RZ = RZ - BOXL*ANINT(RZ/BOXL)
                    Rij(m)=sqrt(RX**2+RY**2+RZ**2)
                end if
            end do

            call sort(N-1,Rij,0.0000001,indaux)

            do l=1, 6
                rcut = rcut + Rij(l)
            end do

        end do
        rcut = ((1.0+sqrt(2.0))/2.0)*rcut/(6.0*N)
        write(*,*) "rcut adaptativo: ", rcut

    else
        rcut = RcutFijo
    end if



    !Ciclo de calculos

    fcc = 0
    hcp = 0
    bcc = 0
    icos = 0
    no = 0

    norm = 0

    f421 = 0
    f311 = 0
    f422 = 0
    f221 = 0
    f441 = 0
    f321 = 0
    f444 = 0
    f666 = 0
    f665 = 0
    f544 = 0
    f433 = 0
    f555 = 0

    open(73,file="estporparticula.txt",status="unknown")
    write(73,*) "  i", "  221", "  311", "  321", "  421", "  422", "  433", "  441", "  444", "  544", "  555", "  665", "  666"
    write(73,*) " "

    do i=1, N !Cambio de particula de referencia
        !write(*,*) "... Progreso: particula", i, " ..."

        !Obtener vecindad de particula de referencia i
        do j=1, N
            if (i /= j) then
                RX = X(i)-X(j)
                RY = Y(i)-Y(j)
                RZ = Z(i)-Z(j)

                !Distancia y condicion de imagen minima
                RX = RX - BOXL*ANINT(RX/BOXL)
                RY = RY - BOXL*ANINT(RY/BOXL)
                RZ = RZ - BOXL*ANINT(RZ/BOXL)
                dist(j)=sqrt(RX**2+RY**2+RZ**2) !Se guardan las distancias
            end if
        end do
        ind(1) = 0
        call sort(N,dist,0.0000001,ind) !Se ordenan de menor a mayor

        do j=1, N
            if (dist(j) < rcut) then
                cycle
            else
                exit
            end if
        end do
        j=j-1 !Se obtiene el numero de vecinos de la referencia

        !f421 = 0
        !f311 = 0
        !f422 = 0
        !f221 = 0
        !f441 = 0
        !f321 = 0
        !f444 = 0
        !f666 = 0

        !Analisis de vecindad con 1 (referencia) + j (vecindad) particulas
        do m=1, j
            icna = 0
            jcna = 0
            kcna = 0
            jind = 0

            do l=1, j
                !Identificar indice i para particula m de vecindad
                if (m /= l) then
                    RX = X(ind(m))-X(ind(l))
                    RY = Y(ind(m))-Y(ind(l))
                    RZ = Z(ind(m))-Z(ind(l))
                    RX = RX - BOXL*ANINT(RX/BOXL) !Distancia y condicion de imagen minima
                    RY = RY - BOXL*ANINT(RY/BOXL)
                    RZ = RZ - BOXL*ANINT(RZ/BOXL)
                    distint(l)=sqrt(RX**2+RY**2+RZ**2)
                    if (distint(l) < rcut) then
                        icna = icna + 1
                        jind(icna) = ind(l) !Guarda para calculo de j
                    end if
                end if
            end do

            !Identificar indice j para particula m de vecindad
            do l=1, icna-1
                do s=l+1, icna
                    if (s /= l) then
                        RX = X(jind(l))-X(jind(s))
                        RY = Y(jind(l))-Y(jind(s))
                        RZ = Z(jind(l))-Z(jind(s))
                        RX = RX - BOXL*ANINT(RX/BOXL) !Distancia y condicion de imagen minima
                        RY = RY - BOXL*ANINT(RY/BOXL)
                        RZ = RZ - BOXL*ANINT(RZ/BOXL)
                        dist(l)=sqrt(RX**2+RY**2+RZ**2)
                        if (dist(l) < rcut) then
                            jcna = jcna + 1
                            indk(jcna,1) = jind(l)
                            indk(jcna,2) = jind(s)
                            kcna = 1
                        end if
                    end if
                end do
            end do

            !Identificar indice k para particula m de vecindad
            do l=1, jcna-1 !Verificar camino de long. 2
                do r1=1, 2
                    do s=l+1, jcna
                        do r2=1, 2
                            if (indk(l,r1) == indk(s,r2)) then
                                if (kcna < 2) then
                                    kcna = 2
                                end if



                                if (r2 == 2) then !Verificar camino de long. 3, lado1
                                    r2a = 1
                                else
                                    r2a = 2
                                end if
                                encontre1 = .false.
                                encontre2 = .false.
                                do s3=1, jcna
                                    if (s3 /= l .and. s3 /= s) then
                                        do r3=1, 2
                                            if (indk(s,r2a) == indk(s3,r3)) then
                                                if (kcna < 3) then
                                                    kcna = 3
                                                end if
                                                encontre1 = .true.

                                                if (r3 == 2) then !Verificar camino de long4
                                                    r3a = 1
                                                else
                                                    r3a = 2
                                                end if
                                                do s4=1, jcna
                                                    if (s4 /= l .and. s4 /= s .and. s4 /= s3) then
                                                        do r4=1, 2
                                                            if (indk(s3,r3a) == indk(s4,r4)) then
                                                                if (kcna < 4) then
                                                                    kcna = 4
                                                                end if
                                                                encontre2 = .true.
                                                            end if
                                                        end do
                                                    end if
                                                end do


                                            end if
                                        end do
                                    end if
                                end do
                                if (r1 == 2) then !Verificar camino de long. 3, lado2, y si existe cola para long4 o 5
                                    r2a = 1
                                else
                                    r2a = 2
                                end if
                                do s3=1, jcna
                                    if (s3 /= l .and. s3 /= s) then
                                        do r3=1, 2
                                            if (indk(l,r2a) == indk(s3,r3)) then
                                                if (kcna < 3) then
                                                    kcna = 3
                                                end if
                                                if (encontre1 .and. k < 4) then
                                                    kcna = 4
                                                end if
                                                if (encontre2 .and. k < 5) then
                                                    kcna = 5
                                                end if

                                                if (r3 == 2) then !Verificar camino de long4, y si existe cola para long 5 o 6
                                                    r3a = 1
                                                else
                                                    r3a = 2
                                                end if
                                                do s4=1, jcna
                                                    if (s4 /= l .and. s4 /= s .and. s4 /= s3) then
                                                        do r4=1, 2
                                                            if (indk(s3,r3a) == indk(s4,r4)) then
                                                                if (kcna < 4) then
                                                                    kcna = 4
                                                                end if
                                                                if (encontre1 .and. k < 5) then
                                                                    kcna = 5
                                                                end if
                                                                if (encontre2 .and. k < 6) then
                                                                    kcna = 6
                                                                end if
                                                            end if
                                                        end do
                                                    end if
                                                end do

                                            end if
                                        end do
                                    end if
                                end do

                            end if
                        end do
                    end do
                end do
            end do

            if (kcna > jcna) then !Evitar contar de mas en caminos cerrados
                kcna = jcna
            end if

            !Guardar estructura
            if (icna == 4 .and. jcna == 2 .and. kcna == 1) then
                fcc = fcc + 1
            else if (icna == 4 .and. jcna == 2 .and. kcna == 2) then
                fcc = fcc + 1
            else if (icna == 4 .and. jcna == 4 .and. kcna == 4) then
                bcc = bcc + 1
            else if (icna == 6 .and. jcna == 6 .and. kcna == 5) then
                bcc = bcc + 1
            else if (icna == 4 .and. jcna == 3 .and. kcna == 3) then
                icos = icos + 1
            else if (icna == 5 .and. jcna == 4 .and. kcna == 4) then
                icos = icos + 1
            else if (icna == 5 .and. jcna == 5 .and. kcna == 5) then
                icos = icos + 1
            else
            no = no + 1
            end if

            !Escribir cuenta para comparar
            if (icna == 4 .and. jcna == 2 .and. kcna == 2) then
                f422 = f422 + 1
            else if (icna == 4 .and. jcna == 2 .and. kcna == 1) then
                f421 = f421 + 1
            else if (icna == 3 .and. jcna == 1 .and. kcna == 1) then
                f311 = f311 + 1
            else if (icna == 2 .and. jcna == 2 .and. kcna == 1) then
                f221 = f221 + 1
            else if (icna == 4 .and. jcna == 4 .and. kcna == 1) then
                f441 = f441 + 1
            else if (icna == 3 .and. jcna == 2 .and. kcna == 1) then
                f321 = f321 + 1
            else if (icna == 6 .and. jcna == 6 .and. kcna == 6) then
                f666 = f666 + 1
            else if (icna == 4 .and. jcna == 4 .and. kcna == 4) then
                f444 = f444 + 1
            else if (icna == 5 .and. jcna == 5 .and. kcna == 5) then
                f555 = f555 + 1
            else if (icna == 6 .and. jcna == 6 .and. kcna == 5) then
                f665 = f665 + 1
            else if (icna == 5 .and. jcna == 4 .and. kcna == 4) then
                f544 = f544 + 1
            else if (icna == 4 .and. jcna == 3 .and. kcna == 3) then
                f433 = f433 + 1
            end if




            !write(*,*) icna, jcna, kcna
            !write(48,*) f422, f421, f311, f221, f441, f321, f666, f444
            !write(48,*) icna, jcna, kcna


        end do

        !Identificar estructura
    !    if (f311 == 6 .and. f421 == 3) then
    !        fcc111 = fcc111 + 1
    !    else if (f421 == 12) then
    !        fcc = fcc + 1
    !    else if (f421 == 6 .and. f422 == 6) then
    !        hcp = hcp + 1
    !    else if (f221 == 4 .and. f441 == 4) then
    !        fcc100f = fcc100f + 1
    !    else if (f321 == 6 .and. f441 == 3) then
    !        fcc111f = fcc111f + 1
    !    else if (f666 == 8 .and. f444 == 6) then
    !        bcc = bcc + 1
    !    end if

    write(73,*) i, f221, f311, f321, f421, f422, f433, f441, f444, f544, f555, f665, f666

    end do


    norm = no + fcc + bcc + icos

    fcc = fcc/norm
    bcc = bcc/norm
    icos = icos/norm
    !hcp = hcp/N

    !fcc111 = fcc111/N
    !fcc111f = fcc111f/N
    !fcc100f = fcc100f/N


    open(48,file="estructura.txt",status="unknown")
    open(70,file="cuentaest.txt",status="unknown")


    write(48,*) fcc, bcc, icos

    write(70,*) " 221 ", f221, f221/norm
    write(70,*) " 311 ", f311, f311/norm
    write(70,*) " 321 ", f321, f321/norm
    write(70,*) " 421 ", f421, f421/norm
    write(70,*) " 422 ", f422, f422/norm
    write(70,*) " 433 ", f433, f433/norm
    write(70,*) " 441 ", f441, f441/norm
    write(70,*) " 444 ", f444, f444/norm
    write(70,*) " 544 ", f544, f544/norm
    write(70,*) " 555 ", f555, f555/norm
    write(70,*) " 665 ", f665, f665/norm
    write(70,*) " 666 ", f666, f666/norm


    close(48)
    close(70)
    close(73)





end program









! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! Adapted from flib/hpsort_eps
  !---------------------------------------------------------------------
  subroutine sort(n, ra, eps, ind)
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps".
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  implicit none
  !-input/output variables
  integer, intent(in)   :: n
  real, intent(in)  :: eps
  integer :: ind (n)
  real :: ra (n)
  !-local variables
  integer :: i, ir, j, l, iind
  real :: rra
!
  ! initialize index array
  IF (ind (1) .eq.0) then
     DO i = 1, n
        ind (i) = i
     ENDDO
  ENDIF

  ! nothing to order
  IF (n.lt.2) return
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1

  ir = n

  sorting: do

    ! still in hiring phase
    IF ( l .gt. 1 ) then
       l    = l - 1
       rra  = ra (l)
       iind = ind (l)
       ! in retirement-promotion phase.
    ELSE
       ! clear a space at the end of the array
       rra  = ra (ir)
       !
       iind = ind (ir)
       ! retire the top of the heap into it
       ra (ir) = ra (1)
       !
       ind (ir) = ind (1)
       ! decrease the size of the corporation
       ir = ir - 1
       ! done with the last promotion
       IF ( ir .eq. 1 ) then
          ! the least competent worker at all !
          ra (1)  = rra
          !
          ind (1) = iind
          exit sorting
       ENDIF
    ENDIF
    ! wheter in hiring or promotion phase, we
    i = l
    ! set up to place rra in its proper level
    j = l + l
    !
    DO while ( j .le. ir )
       IF ( j .lt. ir ) then
          ! compare to better underling
          IF ( hslt( ra (j),  ra (j + 1) ) ) then
             j = j + 1
          !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
           !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
          ENDIF
       ENDIF
       ! demote rra
       IF ( hslt( rra, ra (j) ) ) then
          ra (i) = ra (j)
          ind (i) = ind (j)
          i = j
          j = j + j
       !else if ( .not. hslt ( ra(j) , rra ) ) then
          !this means rra == ra(j) within tolerance
          ! demote rra
         ! if (iind.lt.ind (j) ) then
         !    ra (i) = ra (j)
         !    ind (i) = ind (j)
         !    i = j
         !    j = j + j
         ! else
             ! set j to terminate do-while loop
         !    j = ir + 1
         ! endif
          ! this is the right place for rra
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1
       ENDIF
    ENDDO
    ra (i) = rra
    ind (i) = iind

  END DO sorting
contains

  !  internal function
  !  compare two real number and return the result

  logical function hslt( a, b )
    REAL :: a, b
    IF( abs(a-b) <  eps ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
  end function hslt

  !
end subroutine
