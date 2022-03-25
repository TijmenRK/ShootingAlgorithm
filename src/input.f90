module input

    public
    
    ! Define parameters

    integer, parameter :: In_GridPoints = 1000                         ! Amount of Gridpoints (Diagonalization with >1000 Gridpoints takes a long time)
    real*8, parameter :: In_Length = 4.                                ! Grid with X values from -Length/2 to Length/2

    real*8, parameter :: In_BoundaryCond(4) = (/0.,0.001,0.001,0./)    ! Boundaries of eigenfunctions approach 0
    
    real*8, parameter :: In_PotentialLengthofTotal = 0.9               ! Fraction of the X values within the potential well (two-sided around 0)
    real*8, parameter :: In_PotentialConstOutofWell = 1000             ! Define the constant within the potential well (Constant within well = 0)

    ! Output file
    character(30), parameter :: In_FileName = 'output.txt'             ! Name of the output file
    integer, parameter :: In_AmountofStates = 10                       ! Number of Eigenstates printed to the output file (from 1 to ...)

end