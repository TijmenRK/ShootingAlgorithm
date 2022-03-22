program main
    use GridSetup
    implicit none
    real, allocatable :: grid(:)
    real :: Length, LowerBound, UpperBound
    integer :: GridPoints
    
    Length = 2.
    GridPoints = 20
    allocate(grid(GridPoints))

    LowerBound = -Length/2.
    UpperBound = Length/2.

    grid = CreateGrid(GridPoints,LowerBound,UpperBound)

    
end program main