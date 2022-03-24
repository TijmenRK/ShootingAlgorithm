module ShootingMethod
    use PotentialModule
    implicit none

    public :: Shooting
    private

    contains

    subroutine Shooting(ApproxEigenvalue, BoundaryCond, Grid, Eigenfunctions)
        real*8, intent(in) :: ApproxEigenvalue, Grid(:)
        real*8, intent(in) :: BoundaryCond(4)
        real*8, allocatable, intent(out) :: Eigenfunctions(:)
        real*8, allocatable :: EigenfunctionsIn(:), EigenfunctionsOut(:)
        real*8 :: TrialEigenvalue, TempNewEigenvalue, EigenvalueCorrection
        real*8 :: Length, NeighbDistance
        integer :: MatchingPoint, i, GridSize, MPIn
        
        if (modulo(size(Grid),2) == 0) then
            MatchingPoint = size(Grid)/2
            MPIn = 2
        else
            MatchingPoint = (size(Grid)+1)/2
            MPIn = 3
        endif

        allocate(EigenfunctionsIn(MatchingPoint+2))
        allocate(EigenfunctionsOut(MatchingPoint+2))

        GridSize = size(Grid)
        Length = Grid(Gridsize) - Grid(1)
        NeighbDistance = Length/GridSize

        call AssignBoundary(BoundaryCond, EigenfunctionsIn)
        call AssignBoundary(BoundaryCond, EigenfunctionsOut)

        TrialEigenvalue = ApproxEigenvalue
        do
            do i=3, MatchingPoint+2
                EigenfunctionsOut(i) = -EigenfunctionsOut(i-2) + 2.*NeighbDistance**2.*(Potential(Grid(i-1),Length)- &
                TrialEigenvalue+1./NeighbDistance**2.)*EigenfunctionsOut(i-1)

                EigenfunctionsIn(size(EigenfunctionsIn)-i+1) = -EigenfunctionsIn(size(EigenfunctionsIn)-i+3) + 2.* &
                    NeighbDistance**2.*(Potential(Grid(GridSize-i+2),Length)-TrialEigenvalue+1./NeighbDistance**2.)* &
                    EigenfunctionsIn(size(EigenfunctionsIn)-i+2)
            enddo

            call Normalize(EigenfunctionsIn)
            call Normalize(EigenfunctionsOut)

            EigenvalueCorrection = 1./2.*( &
                (Derivative(EigenfunctionsIn,NeighbDistance,MPIn)/EigenfunctionsIn(MPIn))- &
                (Derivative(EigenfunctionsOut,NeighbDistance,MatchingPoint)/EigenfunctionsOut(MatchingPoint)))* &
                1/(sum(EigenfunctionsOut(1:MatchingPoint)**2)/EigenfunctionsOut(MatchingPoint)**2 &
                + sum(EigenfunctionsIn(MPIn:size(EigenfunctionsIn))**2)/EigenfunctionsIn(MPIn)**2)

            TempNewEigenvalue = TrialEigenvalue - EigenvalueCorrection
            if (abs(TrialEigenvalue - TempNewEigenvalue) < 1E-13) then
                exit
            else
                TrialEigenvalue = TempNewEigenvalue
            endif
        enddo

        deallocate(EigenfunctionsIn)
        deallocate(EigenfunctionsOut)

        allocate(Eigenfunctions(GridSize))
        call AssignBoundary(BoundaryCond, Eigenfunctions)

        do i=3, GridSize
            Eigenfunctions(i) = -Eigenfunctions(i-2) + 2.*NeighbDistance**2.*(Potential(Grid(i-1),Length)- &
            TrialEigenvalue+1./NeighbDistance**2.)*Eigenfunctions(i-1)
        enddo
        call Normalize(Eigenfunctions)
    end subroutine


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

    real*8 function Derivative(Array, NeighbDistance, Index)
        real*8, intent(in) :: Array(:), NeighbDistance
        integer, intent(in) :: Index
        Derivative = (Array(Index+1)-Array(Index-1))/(2*NeighbDistance)
    end
end module ShootingMethod