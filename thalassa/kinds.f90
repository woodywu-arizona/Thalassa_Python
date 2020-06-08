module KINDS
! Description:
!    Contains values for the kind parameters used in Thalassa. To work in
!    quadruple precision, set dk = qk = 33.
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    davideamato@email.arizona.edu
!
! ==============================================================================
implicit none



integer,parameter  ::  ik = selected_int_kind(9)
integer,parameter  ::  dk = selected_real_kind(15)

TYPE structOut
integer, pointer  :: npts;
real(dk), pointer :: orb_point(:,:);
real(dk), pointer :: cart_point(:,:);
END TYPE structOut

end module KINDS
