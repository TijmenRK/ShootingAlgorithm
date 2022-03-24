module PotentialModule
    implicit none

    public :: Potential
    private 

    contains
    
    real function Potential(GridPoint, Length)
    real*8, intent(in) :: GridPoint, Length
    real*8 :: PotentialLength
    PotentialLength = Length*0.9
    if (abs(GridPoint) <= PotentialLength/2) then
        Potential = 0
    else
        Potential = 10000
    endif
    end
end module PotentialModule