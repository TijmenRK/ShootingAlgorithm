module GridSetup
    implicit none

    public :: CreateGrid
    private

    contains

    function CreateGrid(GridPoints, LowerBound, UpperBound) result(Grid)
        integer, intent(in) :: GridPoints
        real*8, intent(in) :: LowerBound, UpperBound
        real*8, allocatable :: Grid(:)
        real*8 :: NeighbDistance
        integer :: i

        NeighbDistance = (UpperBound-LowerBound)/(GridPoints-1)
        
        allocate(Grid(GridPoints))

        do i=1, GridPoints
            Grid(i) = LowerBound + (i-1)*NeighbDistance
        enddo
    end
end module GridSetup