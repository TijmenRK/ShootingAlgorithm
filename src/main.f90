program main
    use GridSetup
    use Initialization
    use ShootingMethod

    implicit none
    real*8, allocatable :: grid(:)
    real*8 :: boundarycond(4)
    real*8, allocatable :: eigenvalues(:), EigenFunctions(:)
    real*8 :: Length, LowerBound, UpperBound
    character(len=30) :: FileName
    integer :: GridPoints, iu, i
    
    Length = 2.
    GridPoints = 1000

    allocate(grid(GridPoints))
    allocate(eigenvalues(GridPoints))

    boundarycond(1) = 0
    boundarycond(2) = 0.001
    boundarycond(3) = 0.001
    boundarycond(4) = 0

    LowerBound = -Length/2.
    UpperBound = Length/2.

    grid = CreateGrid(GridPoints,LowerBound,UpperBound)
    ! write (*,'(*(f8.4))') grid
    eigenvalues = Eigenvalue3Point(grid)
    call Shooting(eigenvalues(6), boundarycond, grid, EigenFunctions)

    FileName = 'data.txt'
    open(newunit=iu,file=FileName,action="write",status="replace")
    do i=1, GridPoints
        write (iu,*) grid(i), EigenFunctions(i)
    enddo
    
end program main