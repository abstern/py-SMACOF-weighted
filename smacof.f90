subroutine smacof_embed(numvertices, target_dimension, weights, lengths, Vplus, maximum_iterations, X)
  implicit none

  ! input
  
  ! dimension
  integer, intent(in) :: numvertices
  ! iterations cap
  integer, intent(in) :: maximum_iterations
  ! target dimension
  integer, intent(in) :: target_dimension
  ! the edge lengths
  real, dimension (numvertices, numvertices), intent(in) :: lengths
  ! the edge weights 
  real, dimension (numvertices, numvertices), intent(in) :: weights

  ! storage during calculation and output
  integer step

  ! V_+ is th emoore-penrose inverse of the V matrix (the 'hessian') that holds the transformed edges, weights
  ! calculated (once) using numpy for convenience
  real, dimension (numvertices, numvertices) :: Vplus

  real :: stress

  ! previous step
  real, dimension (numvertices, target_dimension) :: Z
  ! B = B(Z) 
  real, dimension (numvertices, numvertices) :: B

  ! output

  ! array of resulting coordinates
  real, dimension (numvertices, target_dimension), intent(out) :: X

  call RANDOM_NUMBER(X)

  do step=1,maximum_iterations
     Z = X
     call BTransform(numvertices, target_dimension, weights, lengths, Z, B)
     X = matmul(Vplus, matmul(B, Z))
     call FindStress(numvertices, target_dimension, weights, lengths, X, stress)
     print*, "step: ",step,"stress: ",stress
  end do
     
end subroutine smacof_embed

subroutine FindStress(numvertices, target_dimension, weights, lengths, Z, stress)
  implicit none

  ! input
  integer, intent(in) :: target_dimension
  integer, intent(in) :: numvertices
  real, intent(in), dimension (numvertices, target_dimension) :: Z
  real, dimension (numvertices, numvertices), intent(in) :: lengths
  real, dimension (numvertices, numvertices), intent(in) :: weights

  ! intermediate
  real :: distance
  integer i, j, k 
  
  ! output
  real, intent(out) :: stress

  stress = 0.0

  ! i != j
  do j=2,numvertices
     do i=1,j-1
        distance = 0.0
        do k=1,target_dimension
           distance = distance + (Z(i, k) - Z(j,k))**2
        end do
        distance = sqrt(distance)
        stress = stress + (distance - lengths(i, j))**2 * weights(i, j)
     end do
  end do

end subroutine FindStress


subroutine BTransform(numvertices, target_dimension, weights, lengths, Z, B) 
  implicit none

  ! input
  integer, intent(in) :: target_dimension
  integer, intent(in) :: numvertices
  real, intent(in), dimension (numvertices, target_dimension) :: Z
  real, dimension (numvertices, numvertices), intent(in) :: lengths
  real, dimension (numvertices, numvertices), intent(in) :: weights

  ! intermediate
  real :: distance
  integer i, j, k
  
  ! output
  real, dimension (numvertices, numvertices), intent(out) :: B

  ! initialize
  B = 0.0

  ! i != j
  do j=2,numvertices
     do i=1,j-1
        distance = 0.0
        do k=1,target_dimension
           distance = distance + (Z(i, k) - Z(j,k))**2
        end do
        distance = sqrt(distance)
        if (distance > 0) then
           B(i, j) = - weights(i, j) * lengths(i, j) / distance
           B(j, i) = - weights(j, i) * lengths(j, i) / distance
        else
           B(i, j) = 0
           B(j, i) = 0
        end if
        
     end do
  end do

  ! i = j
  do i=1, numvertices
     do j=1, i-1
        B(i, i) = B(i, i) - B(i, j)
     end do
     do j=i+1, numvertices
        B(i, i) = B(i, i) - B(i, j)
     end do
  end do

end subroutine BTransform
