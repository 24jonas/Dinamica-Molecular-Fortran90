program main
    !********************************************************************************************
    !
    ! Programa principal que controla simulación de dinámica
    ! molecular clásica.
    !
    ! Estructura del programa:
    !  1- Declaración de variables
    !  2- Lectura de parámetros físicos y de simulación
    !  3- Cálculos preliminares y condición inicial
    !  4- Simulación
    !  5- Cálculo de propiedades
    !  6- Mensajes finales
    !
    ! Notas:
    !  - Se utilizan unidades reducidas sigma=1
    !
    ! Autor: Jonás Valenzuela Terán
    ! Agradecimientos: Dra. Laura Laura Lorenia Yeomans Reyna, Dr. Efrain Urrutia Bañuelos y
    !                  clase de Licenciatura en Física, Unison 2016-2020
    !
    !********************************************************************************************


    !****************************
    ! 1- Declaración de variables
    !****************************

    implicit none

    integer :: N                                    !Numero de particulas
    integer :: NSTEP                                !Numero de ciclos
    integer :: NENER                                !Configuracion donde se alcanza el equilibrio
    integer :: ISAVE                                !Frecuencia de guardado de configuraciones
    integer :: SEED                                 !Semillas para procesos aleatorios
    real :: DENS                                    !Concentracion reducida
    real :: TEMP                                    !Temperatura reducida
    real :: ENPOT                                   !Energia potencial

    integer :: NTRAZ                                !Numero id de particula trazadora
    real :: BOXL                                    !Longitud de arista de cubo de simulación
    real :: RCUT                                    !Radio de corte para optimizacion
    integer :: NN2                                  !Numero de configuraciones por guardar
    integer :: NN1                                  !Num de conf para variables dinamicas
    real :: etot                                    !Energia cinetica
    real, dimension(:), allocatable :: X, Y, Z      !Posiciones de cada particula (N)
    real, dimension(:), allocatable :: XR, YR, ZR   !Conf. ini origen
    real, dimension(:), allocatable :: XM, YM, ZM   !Posiciones del paso anterior (N)
    real, dimension(:), allocatable :: VX, VY, VZ   !Velocidades de cada particula (N)
    real, dimension(:), allocatable :: FX, FY, FZ   !Fuerzas de cada particula (N)
    real, dimension(:,:), allocatable :: CX, CY, CZ !Configuraciones guardadas (N,NN2)
    real, dimension(:,:), allocatable :: CXR, CYR, CZR !Conf. para difusion (N, NN2)
    real, dimension(:,:), allocatable :: CVX, CVY, CVZ !Velocidades guardadas (N, NN2)

    real, dimension(:), allocatable :: AX, AY, AZ   !Aceleraciones de cada particula
    real, dimension(:), allocatable :: X3, Y3, Z3   !Tercera derivada de posicion
    real, dimension(:), allocatable :: X4, Y4, Z4   !Cuarta derivada de posicion
    real, dimension(:), allocatable :: X5, Y5, Z5   !Quinta derivada de posicion

    real :: dt                                      !Tiempo transcurrido y paso de tiempo

    integer :: i, j, k, l                           !Indices de barrido para ciclos
    real, parameter :: PI = 3.14159265

    integer, external :: fccsize                    !Calcula N de particulas en arreglo fcc

    integer :: TEMPi                                !Generar configuracion inicial (1), o leerla (2)
    real :: TEMPred                                 !Factor de referencia para definir temperatura reducida
    real :: AHEAT                                   !Factor para escalar velocidades


    !**************************
    ! 2- Lectura de parámetros
    !**************************

    write(*,*) "******************************************************************"
    write(*,*) "                 Simulacion de dinamica molecular                 "
    write(*,*) "******************************************************************"
    write(*,*) " "

    write(*,*) "De el numero de particulas"
    !read(*,*) N
    N=500
    write(*,*) " "

    write(*,*) "De el numero de ciclos"
    !read(*,*) NSTEP
    NSTEP=50000
    write(*,*) " "

    write(*,*) "Describa el numero de configuracion donde se alcanza el equilibrio"
    !read(*,*) NENER
    NENER=30000
    write(*,*) " "

    write(*,*) "Cada cuantos pasos guardo configuraciones?"
    !read(*,*) ISAVE
    ISAVE=100
    write(*,*) " "

    write(*,*) "Introduzca un numero semilla para seleccion de trazadora"
    !read(*,*) SEED
    SEED=43
    write(*,*) " "

    write(*,*) "Concentración reducida"
    !read(*,*) DENS
    DENS= 1.312 !Para Ni           !1.403 !Para Cu                    !1.38 !Para Ag                   !1.366 !Para Pt
    write(*,*) " "

    TEMPred = 1.9416e-5 !con U0 de Ni          !2.4343e-5 !con U0 de Cu                      !2.91232E-5 !con U0 de Ag             !1.4723e-5 !con U0 de Pt
    TEMP = 1728.0*TEMPred !k boltzman eV- entre U0 de gupta
    TEMPi = 1 !Elegir si realizar configuracion inicial (1), o leerla (2) (Solo debe cambiarse TEMP)
    dt = 0.008

    !Temperaturas de fusion experimentales:
    ! Ni 1728k
    ! Cu 1356k
    ! Pt 2047k
    ! Ag 1234k


    !**********************************************
    ! 3- Cálculos preliminares y condición inicial
    !**********************************************

    !Calculo de parámetros dependientes
    NN2 = int((NSTEP-NENER)/ISAVE)          !Numero de estados en equilibrio muestreados
    NN1 = int(NSTEP/ISAVE)                  !Numero de estados muestreados para variables dinámicas (incluye no equilibrio)
    BOXL = ((real(N)/DENS)**(1./3.))        !Longitud de el contenedor cubico
    RCUT = BOXL/2.0                         !Radio de corte
    AHEAT = 3*Real(N)*dt*dt*TEMP            !Escala de velocidades para cumplir con equipartición de la energía

    !Asignar tamaño a arreglos que se utilizarán en la simulación
    allocate(X(N),Y(N),Z(N))
    allocate(XR(N),YR(N),ZR(N))
    allocate(VX(N),VY(N),VZ(N))
    allocate(AX(N),AY(N),AZ(N),X3(N),Y3(N),Z3(N),X4(N),Y4(N),Z4(N),X5(N),Y5(N),Z5(N)) !Para predictor-corrector
    allocate(XM(N),YM(N),ZM(N)) !Para verlet
    allocate(FX(N),FY(N),FZ(N))
    allocate(CX(N,NN2),CY(N,NN2),CZ(N,NN2))
    allocate(CXR(N,NN1),CYR(N,NN1),CZR(N,NN1))
    allocate(CVX(N,NN1),CVY(N,NN1),CVZ(N,NN1))

    !Generar o leer configuración inicial
    if (TEMPi == 1) then

        write(*,*) "...Creando configuracion inicial..."

        call configinifccalt(N,X,Y,Z,BOXL) !Funciona solo con N=500 y casos especiales

        !Propuesta general de posiciones iniciales
        !if (DENS <= 0.63) then
        !    call configiniale(BOXL,N,X,Y,Z,SEED,DENS)
        !    write(*,*) "Aleatoria"
        !else
        !    call configinireg(BOXL,N,X,Y,Z,DENS)
        !    write(*,*) "Regular, N=",N
        !end if

        call srand(SEED) !Seleccionar particula trazadora al azar
        NTRAZ = int(rand()*N)
        if (NTRAZ == 0) NTRAZ = 1

        !Velocidades iniciales
        call veliniale(N,TEMP,dt,X,Y,Z,XM,YM,ZM,VX,VY,VZ) !Nota: Unidades reducidas de posicion a velocidades

        XR = 0
        YR = 0
        ZR = 0

        !Aceleraciones iniciales
        call fuerzas(N,BOXL,X,Y,Z,FX,FY,FZ,ENPOT)
        do i=1, N
            AX(i) = FX(i)*dt*dt/2.0
            AY(i) = FY(i)*dt*dt/2.0
            AZ(i) = FZ(i)*dt*dt/2.0
        end do

        X3 = 0
        Y3 = 0
        Z3 = 0

        X4 = 0
        Y4 = 0
        Z4 = 0

        X5 = 0
        Y5 = 0
        Z5 = 0

    else

        !Lectura de configuración incial
        open(22,file="puntosfin.txt", status="old")
        open(51,file="velfin.txt", status="old")
        open(52,file="acefin.txt", status="old")
        open(53,file="x3fin.txt", status="old")
        open(54,file="x4fin.txt", status="old")
        open(55,file="x5fin.txt", status="old")
        open(56,file="xrfin.txt", status="old")
        do i=1, N
            read(22,*) X(i), Y(i), Z(i)
            read(51,*) VX(i), VY(i), VZ(i)
            read(52,*) AX(i), AY(i), AZ(i)
            read(53,*) X3(i), Y3(i), Z3(i)
            read(54,*) X4(i), Y4(i), Z4(i)
            read(55,*) X5(i), Y5(i), Z5(i)
            read(56,*) XR(i), YR(i), ZR(i)
        end do
        close(22)
        close(51)
        close(52)
        close(53)
        close(54)
        close(55)
        close(56)

    end if


    !****************
    ! 4- Simulación
    !****************

    write(*,*) "...Moviendo a las particulas..."
    write(*,*) "Conf., E total, Temperatura"

    open(26,file="energiastempconf.txt", status="unknown")
    open(30,file="traza.txt", status="unknown")
    k=0
    l=0

    do i=1, NSTEP
        !Predecir siguente configuracion Gear
        call predictefrain(N,X,Y,Z,XR,YR,ZR,VX,VY,VZ,AX,AY,AZ,X3,Y3,Z3,X4,Y4,Z4,X5,Y5,Z5,FX,FY,FZ,BOXL) !Usa unidades reducidas de posicion para X1, X2, X3, X4, X5, ...

        !Medir fuerzas
        call fuerzas(N,BOXL,X,Y,Z,FX,FY,FZ,ENPOT)

        !Corregir configuracion predicha Gear
        call correctefrain(N,X,Y,Z,XR,YR,ZR,VX,VY,VZ,AX,AY,AZ,X3,Y3,Z3,X4,Y4,Z4,X5,Y5,Z5,FX,FY,FZ,BOXL,dt,AHEAT,etot) !Usa unidades reducidas de posicion para X1, X2, X3, X4, X5, ...

        write(30,*) X(NTRAZ), Y(NTRAZ), Z(NTRAZ)
        !write(*,*) i, ENPOT, etot, ENPOT + etot, TEMP  !Monitoreo

        !Guardado de configuraciones
        if (mod(i,ISAVE) == 0) then
            l=l+1
            do j=1, N
                CXR(j,l) = XR(j)
                CYR(j,l) = YR(j)
                CZR(j,l) = ZR(j)

                CVX(j,l) = VX(j)/dt !/dt para unidades de velocidad en Efrain
                CVY(j,l) = VY(j)/dt
                CVZ(j,l) = VZ(j)/dt
            end do

            if (i > NENER) then
                k=k+1
                do j=1, N
                    CX(j,k) = X(j)
                    CY(j,k) = Y(j)
                    CZ(j,k) = Z(j)
                end do

                write(*,*) i, ENPOT, etot, ENPOT + etot, TEMP
                write(26,*) i, ENPOT, etot, ENPOT + etot, TEMP

                if (i == 30000) then
                    call instgdr(N,BOXL,DENS,X,Y,Z,i) !Util para propiedades dependientes de g(r)
                end if
            end if

        end if

    end do

    close(26)
    close(30)



    !***************************
    ! 5- Calculo de Propiedades
    !***************************
    write(*,*) "...Guardando configuracion final..."
    open(22,file="puntosfin.txt", status="unknown")
    open(51,file="velfin.txt", status="unknown")
    open(52,file="acefin.txt", status="unknown")
    open(53,file="x3fin.txt", status="unknown")
    open(54,file="x4fin.txt", status="unknown")
    open(55,file="x5fin.txt", status="unknown")
    open(56,file="xrfin.txt", status="unknown")
    do i=1, N
        write(22,*) CX(i,NN2), CY(i,NN2), CZ(i,NN2)
        write(51,*) CVX(i,NN2)*dt, CVY(i,NN2)*dt, CVZ(i,NN2)*dt !Para efrain
        write(52,*) AX(i), AY(i), AZ(i)
        write(53,*) X3(i), Y3(i), Z3(i)
        write(54,*) X4(i), Y4(i), Z4(i)
        write(55,*) X5(i), Y5(i), Z5(i)
        write(56,*) XR(i), YR(i), ZR(i)
    end do
    close(22)
    close(51)
    close(52)
    close(53)
    close(54)
    close(55)
    close(56)

    write(*,*) "...Calculando Gdr..."
    call gdr(CX,CY,CZ,N,NN2,BOXL,DENS)

    write(*,*) "...Calculando Var dinamicas..."
    call wvdt(CXR,CYR,CZR,CVX,CVY,CVZ,NN1,DT,ISAVE,N)

    write(*,*) "...Calculando Promedios y Desviaciones..."
    call dev(2,NN2,5) !Cambiar parámetros dependiendo de archivo y columna que se desea analizar

    !Limpiar memoria
    deallocate(X,Y,Z)
    deallocate(XR,YR,ZR)
    deallocate(VX,VY,VZ)
    deallocate(AX,AY,AZ,X3,Y3,Z3,X4,Y4,Z4,X5,Y5,Z5)
    deallocate(XM,YM,ZM)
    deallocate(FX,FY,FZ)
    deallocate(CX,CY,CZ)
    deallocate(CXR,CYR,CZR)
    deallocate(CVX,CVY,CVZ)


    !*********************
    ! 6- Mensajes finales
    !*********************
    write(*,*) " "
    write(*,*) "Finalizado"
    write(*,*) " "








end program
