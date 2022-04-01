module ShootingMethod
    use input
    use PotentialModule
    implicit none

    public :: ShootingAlgorithm
    private

    contains

    function ShootingAlgorithm(ApproxEigenvalue, Grid) result(Eigenfunctions)
        real*8, intent(in) :: ApproxEigenvalue, Grid(:)
        real*8 :: Eigenfunctions(In_GridPoints)
        real*8, allocatable :: EigenfunctionsIn(:), EigenfunctionsOut(:)
        real*8 :: TrialEigenvalue, TempNewEigenvalue, EigenvalueCorrection
        real*8 :: NeighbDistance
        integer :: MatchingPoint, i, MPIn
        
        ! MatchingPoint placed in the middle of the grid
        if (modulo(size(Grid),2) == 0) then
            MatchingPoint = size(Grid)/2
            ! Position of MatchingPoint in the EigenfunctionsIn-array
            MPIn = 2
        else
            MatchingPoint = (size(Grid)+1)/2
            MPIn = 3
        endif

        allocate(EigenfunctionsIn(MatchingPoint+2))
        allocate(EigenfunctionsOut(MatchingPoint+2))

        NeighbDistance = In_Length/In_GridPoints

        ! Boundary filled in from input module (once, because never overwritten)
        call AssignBoundary(In_BoundaryCond, EigenfunctionsIn)
        call AssignBoundary(In_BoundaryCond, EigenfunctionsOut)

        TrialEigenvalue = ApproxEigenvalue
        do
            ! Calculate eigenfunctions from two sides with TrialEigenvalue
            do i=3, MatchingPoint+2
                EigenfunctionsOut(i) = -EigenfunctionsOut(i-2) + 2.*NeighbDistance**2.*(Potential(Grid(i-1))- &
                TrialEigenvalue+1./NeighbDistance**2.)*EigenfunctionsOut(i-1)

                EigenfunctionsIn(size(EigenfunctionsIn)-i+1) = -EigenfunctionsIn(size(EigenfunctionsIn)-i+3) + 2.* &
                    NeighbDistance**2.*(Potential(Grid(In_GridPoints-i+2))-TrialEigenvalue+1./NeighbDistance**2.)* &
                    EigenfunctionsIn(size(EigenfunctionsIn)-i+2)
            enddo

            call Normalize(EigenfunctionsIn)
            call Normalize(EigenfunctionsOut)

            ! Calculate a correction for the TrialEigenvalue
            EigenvalueCorrection = 1./2.*( &
                (Derivative(EigenfunctionsIn,MPIn)/EigenfunctionsIn(MPIn))- &
                (Derivative(EigenfunctionsOut,MatchingPoint)/EigenfunctionsOut(MatchingPoint)))* &
                1/(sum(EigenfunctionsOut(1:MatchingPoint)**2)/EigenfunctionsOut(MatchingPoint)**2 &
                + sum(EigenfunctionsIn(MPIn:size(EigenfunctionsIn))**2)/EigenfunctionsIn(MPIn)**2)

            ! Exit loop if correction is small enough, else try again with corrected eigenvalue
            TempNewEigenvalue = TrialEigenvalue - EigenvalueCorrection
            if (abs(TrialEigenvalue - TempNewEigenvalue) < 1E-13) then
                exit
            else
                TrialEigenvalue = TempNewEigenvalue
            endif
        enddo

        deallocate(EigenfunctionsIn)
        deallocate(EigenfunctionsOut)

        ! Start shooting algorithm from left to right (similar to EigenfunctionsOut) with final eigenvalue
        call AssignBoundary(In_BoundaryCond, Eigenfunctions)
        do i=3, In_GridPoints
            Eigenfunctions(i) = -Eigenfunctions(i-2) + 2.*NeighbDistance**2.*(Potential(Grid(i-1))- &
            TrialEigenvalue+1./NeighbDistance**2.)*Eigenfunctions(i-1)
        enddo
        call Normalize(Eigenfunctions)
    end


    subroutine AssignBoundary(BoundaryCond, Array)
        real*8, intent(in) :: BoundaryCond(4)
        real*8, intent(inout) :: Array(:)

        Array(1) = BoundaryCond(1)
        Array(2) = BoundaryCond(2)
        Array(size(Array)-1) = BoundaryCond(3)
        Array(size(Array)) = BoundaryCond(4)
    end subroutine


    subroutine Normalize(Array)
        real*8, intent(inout) :: Array(:)
        real*8 :: HighestValue
        integer :: i

        if (abs(minval(Array)) > maxval(Array)) then
            HighestValue = abs(minval(Array))
        else
            HighestValue = maxval(Array)
        endif

        do i=1, size(Array)
            Array(i) = Array(i)/HighestValue
        enddo
    end subroutine

    real*8 function Derivative(Array, Index)
        real*8, intent(in) :: Array(:)
        integer, intent(in) :: Index
        real*8 :: NeighbDistance
        NeighbDistance = In_Length/In_GridPoints
        Derivative = (Array(Index+1)-Array(Index-1))/(2*NeighbDistance)
    end
end module ShootingMethod