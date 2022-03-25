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

    LowerBound = -In_Length/2.
    UpperBound = In_Length/2.

    
    grid = CreateGrid(In_GridPoints,LowerBound,UpperBound)
    eigenvalues = Eigenvalue3Point(grid)

    allocate(EigenFunctions(In_AmountofStates,In_GridPoints))
    do i=1, In_AmountofStates
        EigenFunctions(i,:) = Shooting(eigenvalues(i), grid)
    enddo

    open(newunit=iu,file=In_FileName,action="write",status="replace")
    write (iu, '(i3,x,a,t30,a,f8.1)') In_AmountofStates, 'Eigenvalues', 'Potential out of well:', In_PotentialConstOutofWell
    write (iu, *) '                         ', eigenvalues(1:In_AmountofStates)
    write (iu, *) ''
    do i=1, In_GridPoints
        write (iu,*) grid(i), EigenFunctions(:,i)
    enddo
    
end program main