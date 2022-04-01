program main
    use input
    use GridSetup
    use Initialization
    use ShootingMethod

    implicit none
    real*8, allocatable :: grid(:)
    real*8, allocatable :: eigenvalues(:), EigenFunctions(:,:)
    real*8 :: LowerBound, UpperBound
    integer :: iu, i

    allocate(grid(In_GridPoints))
    allocate(eigenvalues(In_GridPoints))

    ! Calculate upper and lower values (X) for the grid
    LowerBound = -In_Length/2.
    UpperBound = In_Length/2.
    
    ! Create grid with module GridSetup
    grid = CreateGrid(In_GridPoints,LowerBound,UpperBound)

    ! Calculate trialeigenvalues for the grid with module Initialization
    eigenvalues = Eigenvalue3Point(grid)

    ! Fill matrix EigenFunctions with module ShootingMethod
    allocate(EigenFunctions(In_AmountofStates,In_GridPoints))
    do i=1, In_AmountofStates
        EigenFunctions(i,:) = ShootingAlgorithm(eigenvalues(i), grid)
    enddo

    ! Print matrix EigenFunctions in output file
    open(newunit=iu,file=In_FileName,action="write",status="replace")
    write (iu, '(i3,x,a,t30,a,es8.1)') In_AmountofStates, 'Trialeigenvalues'
    write (iu, *) '                         ', eigenvalues(1:In_AmountofStates)
    write (iu, *) ''
    do i=1, In_GridPoints
        write (iu,*) grid(i), EigenFunctions(:,i)
    enddo
end program main