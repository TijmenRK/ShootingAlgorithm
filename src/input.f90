module input

    public
    
    ! Define parameters

    integer, parameter :: In_GridPoints = 1000                         ! Amount of Gridpoints (Diagonalization with >1000 Gridpoints takes a long time)
    real*8, parameter :: In_Length = 4.                                ! Grid with X values from -Length/2 to Length/2

    real*8, parameter :: In_BoundaryCond(4) = (/0.,0.001,0.001,0./)    ! Boundaries of eigenfunctions approach 0
    

    character(len=20), parameter :: In_PotentialType = 'GaussianBox'       ! Choose 'InfiniteBox', 'FiniteBox, 'GaussianWell', 'GaussianBox'

    ! Parameters for PotentialType 'InfiniteBox'
    real*8, parameter :: In_PotentialLengthofTotal = 0.9                ! Fraction of the X values within the potential well (two-sided around 0)
    ! Extra Parameter for PotentialType 'FiniteBox'
    real*8, parameter :: In_PotentialConstOutofWell = 100               ! Define the constant out of the potential well (Constant within well = 0)

    ! Parameters for PotentialType 'GaussianWell'
    real*8, parameter :: In_V0 = 50
    real*8, parameter :: In_Alpha = 0.5
    ! Extra parameters for PotentialType 'GaussianBox'
    real*8, parameter :: InPotentialGaussianBoxLoT = 0.9                ! Fraction of the X values within the potential well (two-sided around 0)

    ! Output file
    character(30), parameter :: In_FileName = 'output.txt'             ! Name of the output file
    integer, parameter :: In_AmountofStates = 5                        ! Number of Eigenstates printed to the output file (from 1 to ...)

end