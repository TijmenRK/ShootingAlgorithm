program main
    use GridSetup
    use Initialization
    use ShootingMethod

    implicit none
    real*8, allocatable :: grid(:)
    real*8 :: boundarycond(4)
    real*8, allocatable :: eigenvalues(:), EigenFunctions(:,:)
    real*8 :: Length, LowerBound, UpperBound
    character(len=30) :: FileName
    integer :: GridPoints, iu, i, AmountofStates
    
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
    eigenvalues = Eigenvalue3Point(grid)

    AmountofStates = 10
    allocate(EigenFunctions(AmountofStates,GridPoints))
    do i=1, AmountofStates
        EigenFunctions(i,:) = Shooting(eigenvalues(i), boundarycond, grid)
    enddo

    FileName = 'data.txt'
    open(newunit=iu,file=FileName,action="write",status="replace")
    write (iu, '(i3,x,a,t30,a,i10)') AmountofStates, 'Eigenvalues', 'Potential out of well:', 1000
    write (iu, *) '                         ', eigenvalues(1:AmountofStates)
    write (iu, *) ''
    do i=1, GridPoints
        write (iu,*) grid(i), EigenFunctions(:,i)
    enddo
    
end program main