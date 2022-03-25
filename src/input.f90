module input

    public
    
    integer :: In_GridPoints = 1000
    real*8 :: In_Length = 4.
    real*8 :: In_BoundaryCond(4) = (/0.,0.001,0.001,0./)

    
    real*8 :: In_PotentialLengthofTotal = 0.9
    real*8 :: In_PotentialConstOutofWell = 1000

    
    integer :: In_AmountofStates = 10
    character(30) :: In_FileName = 'output.txt'

end