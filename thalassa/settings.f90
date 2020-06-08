module SETTINGS
! Description:
!    Settings for Thalassa.
!
! Author:
!    Davide Amato
!    Space Dynamics Group - Technical University of Madrid
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! Revisions:
!    180531: Add THALASSA version.
!    180608: Add F107 flag.
!    190113: add Moon collision flag.
!
! ==============================================================================

! MODULES
use KINDS, only: dk

implicit none
! Settings file id and path
integer,parameter      ::  id_set = 13
character(len=4096)    ::  input_path
! Physical model
integer  ::  insgrav            ! Non-spherical gravity field flag.
integer  ::  isun               ! 0 = no Sun perturbation, >1 otherwise.
integer  ::  imoon              ! 0 = no Moon perturbation, >1 otherwise.
integer  ::  idrag              ! 0 = no atmospheric drag, 1 = Wertz model, 2 = US76 (PATRIUS), 3 = J77 (Carrara - INPE), 4 = NRLMSISE-00 (Picone - NRL)
integer  ::  iF107              ! 0 = constant F10.7 flux, 1 = variable F10.7 flux
integer  ::  iSRP               ! 0 = no SRP, 1 = otherwise.
integer  ::  iephem             ! Ephemerides source. 1 = DE431 ephemerides. 2 = Simpl. Meeus & Brown
integer  ::  gdeg,gord          ! Gravitational potential - maximum degree and order
integer  ::  Mord,Sord          ! Order of the Legendre expansion for the Moon and the Sun
! Integrator settings
integer  ::  mxstep             ! Max. number of integration/output steps.
real(dk) ::  tol                ! Integrator tolerance.
integer  ::  imcoll             ! 0 = do not check for collisions with the Moon, 1 = otherwise.
! Equations of motion settings
integer  ::  eqs                ! Equations of motion type. 1 = Cowell,
!                                 2 = EDromo(t), 3 = EDromo(c), 4 = EDromo(l)
! Output settings
character(len=512)  ::  outpath
integer  ::  verb

! THALASSA version
character(len=20),parameter   ::  THALASSA_ver = 'v1.2'

contains


subroutine READ_SETTINGS(val1,val2,val3,val4,val5,val6,val7, &
                         val8, val9, val10, val11, val12, val13)
! Description:
!    Reads the settings file.
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! Revisions:
!     180608: add iF107 flag.
!     190113: add Moon collision flag.
! 
! ==============================================================================

implicit none
integer,intent(in) ::  val1,val2,val3,val4,val5,val6,val7,val8,val9,val11,val12,val13
real(dk),intent(in) ::  val10
integer,parameter    ::  hlines = 13  ! <- Check this when modifying input.txt
integer  :: i
character(len=4096)  ::  dummy
real(dk)  ::  rmxstep

insgrav = val1
isun = val2
imoon = val3
idrag = val4
iF107 = val5
iSRP = val6
iephem = val7
gdeg = val8
gord = val9
tol = val10
rmxstep = val11
imcoll = val12
eqs = val13

mxstep = int(rmxstep)

! If either of isun and imoon > 1, interpret it as the order of the Legendre expansion
Sord = isun
Mord = imoon

end subroutine READ_SETTINGS

end module SETTINGS
