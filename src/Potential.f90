module PotentialModule
    use input

    public :: Potential
    private 

    contains
    
    real*8 function Potential(GridPoint)
        real*8, intent(in) :: GridPoint
        real*8 :: PotentialLength

        if (In_PotentialType .eq. 'InfiniteBox' .or. In_PotentialType .eq. 'FiniteBox') then
            PotentialLength = In_Length*In_PotentialLengthofTotal
            if (abs(GridPoint) <= PotentialLength/2) then
                Potential = 0
            else
                if (In_PotentialType .eq. 'InfiniteBox') then
                    Potential = 1000
                elseif (In_PotentialType .eq. 'FiniteBox') then
                    Potential = In_PotentialConstOutofWell
                endif
            endif

        elseif (In_PotentialType .eq. 'GaussianWell' .or. In_PotentialType .eq. 'GaussianBox') then
            if (In_PotentialType .eq. 'GaussianWell') then
                Potential = -In_V0*exp((-(GridPoint**2))/(2.*In_Alpha**2))
            elseif (In_PotentialType .eq. 'GaussianBox') then
                PotentialLength = In_Length*InPotentialGaussianBoxLoT
                if (abs(GridPoint) <= PotentialLength/2) then
                    Potential = -In_V0*exp((-(GridPoint**2))/(2.*In_Alpha**2))
                else
                    Potential = 1000
                endif
            endif            
        else
            print *, 'No valid PotentialType'
        endif

    end
end module PotentialModule