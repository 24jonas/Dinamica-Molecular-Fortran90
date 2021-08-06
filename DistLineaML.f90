subroutine DistLineaML(L,N,sigma,X)

    !**Declaracion de variables
    implicit none
    REAL, intent(in) :: L, sigma
    INTEGER, intent(in) :: N
    INTEGER :: i
    REAL, DIMENSION(N), intent(out):: X


    !**Calculos**
    DO i=1, N
        X(i) = -(L/2.0) + (sigma/2.0) + (L-sigma)*(real(i)-1.0)/(real(N)-1.0)
    END DO

end subroutine
