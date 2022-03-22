module GridSetup
    implicit none

    public :: CreateGrid
    private

    contains

    function CreateGrid(GridPoints, LowerBound, UpperBound) result(Grid)
        integer, intent(in) :: GridPoints
        real, intent(in) :: LowerBound, UpperBound
        real, allocatable :: Grid(:)
        real :: NeighbDistance
        integer :: i

        NeighbDistance = (UpperBound-LowerBound)/real(GridPoints-1)
        
        allocate(Grid(GridPoints))

        do i=1, GridPoints
            Grid(i) = LowerBound + (i-1)*NeighbDistance
        enddo
    end
end module GridSetup