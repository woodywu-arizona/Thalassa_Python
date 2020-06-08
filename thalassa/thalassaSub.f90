subroutine thalassasub(MJD0, COE0, tspan, tstep, insgrav, isun, imoon, idrag, &
                       iF107, iSRP, iephem, gdeg, gord, rmxstep, tol, imcoll, &
                       eqs, SCMass, ADrag, ASRP, CD, CR, mxpts, npts, tag, exitcode, &
                        cart_target, orb_target)
!!
!! Thalassa fortran assambled subroutine. Generate trajectory for given conditions.
!! Notice, this subroutine is created to compile python thalassa module.
!! Input:   MJD0, COE0, tspan, tstep, insgrav, isun, imoon, idrag, iF107, iSRP, iephem, gdeg, gord,             rmxstep, tol, imcoll, eqs, SCMass, ADrag, ASRP, CD, CR, mxpts, npts, tag, exitcode
!!          Check VARIABLES Assignment Stage for detailed declaration and definition
!! Output:  cart_target, matrix, size 7 x mxpts, each column follows [MJD,X,Y,Z,Vx,Vy,Vz]
!!          orb_target, matrix, size 7 x mxpts, each column follows [MJD,SMA,ECC,INC,RAAN,AOP,M]
!!

! MODULES
USE ISO_C_BINDING
use KINDS,       only: dk
use CART_COE,    only: COE2CART,CART2COE
use PHYS_CONST,  only: READ_PHYS
use NSGRAV,      only: INITIALIZE_NSGRAV, INITIALIZE_EOP
use PROPAGATE,   only: DPROP_REGULAR
use SUN_MOON,    only: INITIALIZE_LEGENDRE
use IO,          only: READ_IC,object_path
use SETTINGS,    only: READ_SETTINGS
use PHYS_CONST,  only: GE,d2r,r2d,secsPerDay,twopi
use SUN_MOON,    only: GslSun,GslMoon
use AUXILIARIES, only: set_mjd

implicit none

!! VARIABLES Assignment Stage
double precision  ::  COE0(1:6),COE0_rad(1:6)
real(dk)  ::  R0(1:3),V0(1:3)
! Integration span and dt [solar days]
double precision,intent(in)  ::  tspan,tstep
! Trajectory
integer, intent(in)           :: mxpts, tag
integer, intent(inout)          ::  npts
integer                       ::  ipt, start, finish
real(dk)                      ::  R(1:3),V(1:3)
real(dk)                      ::  A2M_DRAG, A2M_SRP
double precision,intent(in)           ::  SCMass, ADrag, ASRP, CD, CR

! Measurement of CPU time, diagnostics
integer, intent(inout)  ::  exitcode
! Function calls and integration steps
integer  ::  int_steps,tot_calls
! Paths
character(len=512) :: earth_path,phys_path, eop_path
! READ_IC
!real(dk),intent(in)  ::  MJD0
double precision,intent(in)  ::  MJD0
! READ_SETTINGS
! Physical model
integer,intent(in)  ::  insgrav            ! Non-spherical gravity field flag.
integer,intent(in)  ::  isun               ! 0 = no Sun perturbation, >1 otherwise.
integer,intent(in)  ::  imoon              ! 0 = no Moon perturbation, >1 otherwise.
integer,intent(in)  ::  idrag              ! 0 = no atmospheric drag, 1 = Wertz model, 2 = US76 (PATRIUS), 3 = J77 (Carrara - INPE), 4 = NRLMSISE-00 (Picone - NRL)
integer,intent(in)  ::  iF107              ! 0 = constant F10.7 flux, 1 = variable F10.7 flux
integer,intent(in)  ::  iSRP               ! 0 = no SRP, 1 = otherwise.
integer,intent(in)  ::  iephem             ! Ephemerides source. 1 = DE431 ephemerides. 2 = Simpl. Meeus & Brown
integer,intent(in)  ::  gdeg,gord          ! Gravitational potential - maximum degree and order
! Integrator settings
!integer             ::  mxstep             ! Max. number of integration/output steps.
integer,intent(in)  ::  rmxstep             ! Max. number of integration/output steps.
!real(dk),intent(in) ::  tol                ! Integrator tolerance.
double precision,intent(in)  ::  tol
integer,intent(in)  ::  imcoll             ! 0 = do not check for collisions with the Moon, 1 = otherwise.
! Equations of motion settings
integer,intent(in)  ::  eqs                ! Equations of motion type. 1 = Cowell,
! 2 = EDromo(t), 3 = EDromo(c), 4 = EDromo(l)
! OUTPUTS
! format is each column is in the order of [MJD,X,Y,Z,VX,VY,VZ]
double precision, intent(out)  ::  cart_target(7,mxpts)
! format is each column is in the order of [MJD,SMA,ECC,INC,RAAN,AOP,M]
double precision, intent(out)  ::  orb_target(7,mxpts)
integer  :: mjdStart

!! Printing to debug
!print *, MJD0,tspan,tstep, insgrav, isun, imoon, idrag, iF107, iSRP, iephem, gdeg, gord
!write(*,*) "cart size:",size(cart_target)

! compute area to mass ratios
A2M_Drag = ADrag/SCMass
A2M_SRP  = ASRP/SCmass

! paths
earth_path = './thalassa/data/earth_potential/GRIM5-S1.txt'
phys_path = './thalassa/data/physical_constants.txt'
eop_path = './thalassa/data/10_FINALS.DATA_IAU2000_V2013_0110.txt'

!! Initialization Stage
! Read initial conditions, settings and physical model data.
call set_mjd(MJD0)
call READ_SETTINGS(insgrav,isun,imoon,idrag,iF107,iSRP,iephem,gdeg,gord,tol,rmxstep,imcoll,eqs)
call READ_IC(A2M_Drag,A2M_SRP,CD,CR)
call READ_PHYS(phys_path)

! Initialize Earth data
call INITIALIZE_NSGRAV(earth_path)

! Initialize earth orientation data
call INITIALIZE_EOP(eop_path)

! Initialize Legendre coefficients, if needed
if (isun > 1) then
call INITIALIZE_LEGENDRE(isun,GslSun)
end if

if (imoon > 1) then
call INITIALIZE_LEGENDRE(imoon,GslMoon)
end if

! Load SPICE kernels
call FURNSH('./thalassa/data/kernels_to_load.furnsh')

! Convert to Cartesian coordinates
COE0_rad = [COE0(1:2),COE0(3:6)*real(d2r,dk)]
call COE2CART(COE0_rad,R0,V0,GE)

!! Real Calculation Stage
! Davide Propogation with Regulation Magic
call DPROP_REGULAR(R0,V0,MJD0,tspan,tstep,tol,eqs,cart_target,int_steps,tot_calls,exitcode,mxpts,npts)

!! Unloading and Converting
! unload Spice kernels
call UNLOAD('./thalassa/data/kernels_to_load.furnsh')

! debugging output
!write(*,*) "cart_target shape:",shape(cart_target)

! Convert to orbital elements.
do ipt=1,npts
    mjdStart = 1
    start = 2
    finish = 7
    R = cart_target(start:start + 2, ipt)
    V = cart_target(finish - 2:finish, ipt)
    orb_target(mjdStart, ipt)     = cart_target(mjdStart,ipt)
    call CART2COE(R,V,orb_target(start:finish,ipt) ,GE)
    orb_target(start + 2: finish,ipt) = orb_target(start + 2: finish,ipt)/d2r;
!    print *, orb_target(:,ipt)
end do

end subroutine thalassasub

