module Initialization
    use Diagonalization
    use PotentialModule
    implicit none

    public :: Eigenvalue3Point
    private

    contains

    function Eigenvalue3Point(Grid) result(Eigenvalues)
        real*8, intent(in) :: Grid(:)
        real*8, allocatable :: Matrix_S(:,:), Matrix_V(:,:), Matrix_L(:,:), Eigenvalues(:)
        real*8 :: Length, NeighbDistance
        integer :: i

        allocate(Matrix_S(size(Grid),size(Grid)))
        allocate(Matrix_V(size(Grid),size(Grid)))
        allocate(Matrix_L(size(Grid),size(Grid)))
        Matrix_S = 0.
        Matrix_V = 0.
        Matrix_L = 0.
        allocate(Eigenvalues(size(Grid)))

        Length = Grid(size(Grid)) - Grid(1)
        NeighbDistance = Length/size(Grid)

        do i=1, size(Grid)
            Matrix_S(i,i) = -2.
            Matrix_V(i,i) = Potential(Grid(i),Length)
            if (i-1 > 0) then
                Matrix_S(i-1,i) = 1.
            endif
            if (i+1 <= size(Grid)) then
                Matrix_S(i+1,i) = 1.
            endif
        enddo

        Matrix_L = - Matrix_S/(2*NeighbDistance**2) + Matrix_V

        call diagonalize(Matrix_L, eigenvalues=Eigenvalues)
    end
end module Initialization