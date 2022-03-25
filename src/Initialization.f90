module Initialization
    use input
    use Diagonalization
    use PotentialModule
    implicit none

    public :: Eigenvalue3Point
    private

    contains

    function Eigenvalue3Point(Grid) result(Eigenvalues)
        real*8, intent(in) :: Grid(:)
        real*8, allocatable :: Matrix_S(:,:), Matrix_V(:,:), Matrix_L(:,:)
        real*8 :: Eigenvalues(In_GridPoints)
        real*8 :: NeighbDistance
        integer :: i

        allocate(Matrix_S(In_GridPoints,In_GridPoints))
        allocate(Matrix_V(In_GridPoints,In_GridPoints))
        allocate(Matrix_L(In_GridPoints,In_GridPoints))
        Matrix_S = 0.
        Matrix_V = 0.
        Matrix_L = 0.
        Eigenvalues = 0.

        NeighbDistance = In_Length/In_GridPoints

        ! Fill in Matrix V and S
        do i=1, In_GridPoints
            Matrix_S(i,i) = -2.
            Matrix_V(i,i) = Potential(Grid(i))
            if (i-1 > 0) then
                Matrix_S(i-1,i) = 1.
            endif
            if (i+1 <= In_GridPoints) then
                Matrix_S(i+1,i) = 1.
            endif
        enddo

        ! Combine V and S in Matrix L
        Matrix_L = - Matrix_S/(2*NeighbDistance**2) + Matrix_V

        ! Call module diagonalization to calculate eigenvalues for this matrix
        call diagonalize(Matrix_L, eigenvalues=Eigenvalues)
    end
end module Initialization