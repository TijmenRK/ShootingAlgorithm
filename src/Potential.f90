module PotentialModule
    use input
    implicit none

    public :: Potential
    private 

    contains
    
    real*8 function Potential(GridPoint)
    real*8, intent(in) :: GridPoint
    real*8 :: PotentialLength
    PotentialLength = In_Length*In_PotentialLengthofTotal
    if (abs(GridPoint) <= PotentialLength/2) then
        Potential = 0
    else
        Potential = In_PotentialConstOutofWell
    endif
    end
end module PotentialModule