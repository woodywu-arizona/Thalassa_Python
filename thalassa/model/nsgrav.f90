module NSGRAV
! Description:
!    Contains procedures for the calculation of perturbations due to the non-
!    sphericity of the main body. INITIALIZE_NSGRAV reads main body data from a
!    data file, and initializes coefficient matrices.
!    The calculation of the perturbing potential, perturbing acceleration, and
!    the time derivative of the potential in the body-fixed frame (the latter is
!    needed in regularized formulations) takes place in PINES_NSG.
!    NORMFACT gives the normalization factor for the gravitational coefficients.
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! Revisions:
!    180806: Overhaul and implementation of Pines method.
!    181204: Use GMST (IERS 2006 conventions) rather than ERA. Consider Earth
!            Rotation Rate derived from IERS 2006 conventions.
! 
! ==============================================================================

! MODULES
use KINDS,      only: dk
use SETTINGS,   only: gdeg,gord
use IO,         only: id_earth, id_eop
use PHYS_CONST, only: qk,GE,RE,flatt,omegaE,secsPerDay,secsPerSidDay,twopi,&
&ERR_constant
implicit none

! VARIABLES
! Spherical harmonics (unnormalized) 
integer               ::  maxDeg, maxOrd
real(dk),allocatable  ::  Cnm(:,:),Snm(:,:)

! Pines algorithm arrays
real(dk),allocatable  ::  Anm(:,:),Dnm(:,:),Enm(:,:),Fnm(:,:),Gnm(:,:)
real(dk),allocatable  ::  Rm(:),Im(:),Pn(:)
real(dk),allocatable  ::  Aux1(:),Aux2(:),Aux3(:),Aux4(:,:)

! earth orientation parameters array
real(dk),allocatable  ::  eop(:,:)

contains

subroutine INITIALIZE_EOP(eopfile)

character(len=*),intent(in)  ::  eopfile
integer :: io, nlines, i
character(200) :: line
character(200) :: format
real(dk) :: mjd, pmx, pmy, UT1_UTC, LOD, dx00,dy00
nlines = 0

! count lines
open(unit=id_eop,file=trim(eopfile),status='old',action='read')
do
  read(id_eop,*,iostat=io)
  if (io/=0) exit
  nlines = nlines + 1
end do
close(id_eop)

! rewind 
rewind(id_eop)

! allocate
if (allocated(eop)) deallocate(eop)
allocate(eop(1:nlines, 1:7))

! read in eop data
open(unit=id_eop,file=trim(eopfile),status='old',action='read')
do i = 1,nlines

  ! read line
  read(id_eop, '(7X,F8.2,3X,F9.6,10X,F9.6,12X,F10.7,11X,F7.4,11X,F9.3,10X,F9.3,67X)')mjd,pmx,pmy,UT1_UTC,LOD,dx00,dy00

  ! debug 
  ! write(*,*)mjd,pmx,pmy,UT1_UTC,LOD,dx00,dy00

  ! add to array
  eop(i,1) = mjd
  eop(i,2) = pmx
  eop(i,3) = pmy
  eop(i,4) = UT1_UTC
  eop(i,5) = LOD
  eop(i,6) = dx00
  eop(i,7) = dy00
  ! write(*,*)mjd,pmx,pmy,UT1_UTC,LOD,dx00,dy00
enddo
close(id_eop)

! index based on fact that every entry is a new mjd


end subroutine INITIALIZE_EOP





subroutine INITIALIZE_NSGRAV(earthFile)
! Description:
!    Reads Earth gravity data from a text file. Initializes the gravitational
!    parameter, equatorial radius, flattening, rotational velocity, spherical
!    harmonics coefficients, and auxiliary matrices for Pines' algorithm. The
!    latter computes the Earth potential in Cartesian coordinates, avoiding
!    singularities due to the spherical harmonics.
!    Part of this subroutine is due to Hodei Urrutxua (Universidad Rey Juan
!    Carlos, Madrid, Spain) and Claudio Bombardelli (Universidad Politécnica de
!    Madrid, Madrid, Spain).
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! Revisions:
!    180806: Subroutine created from part of READ_PHYS().
! 
! ==============================================================================

! Arguments
character(len=*),intent(in)  ::  earthFile

! Locals
character(len=72)  ::  dummy
real(dk)           ::  invFlatt
integer            ::  i,j,l,m,n


! ==============================================================================

open(unit=id_earth,file=trim(earthFile),status='old',action='read')
read(id_earth,'(a)') (dummy, i=1,4)
read(id_earth,'(a37,i3)') dummy, maxDeg
read(id_earth,'(a37,i3)') dummy, maxOrd
read(id_earth,'(a36,e22.15)') dummy, GE
read(id_earth,'(a36,e22.15)') dummy, RE
read(id_earth,'(a36,e22.15)') dummy, invFlatt
read(id_earth,'(a36,e22.15)') dummy, omegaE

flatt    = 1._dk/invFlatt

! Read and de-normalize spherical harmonics coefficients and auxiliary matrices
! for Pines' algorithm.
! Initialize
l = 1; m = 0;
if (allocated(Cnm)) deallocate(Cnm)
allocate(Cnm(1:maxDeg,0:maxDeg)); Cnm = 0._dk
if (allocated(Snm)) deallocate(Snm)
allocate(Snm(1:maxDeg,0:maxDeg)); Snm = 0._dk
if (allocated(Anm)) deallocate(Anm)
allocate(Anm(0:maxDeg+2, 0:maxDeg+2))
if (allocated(Dnm)) deallocate(Dnm)
allocate(Dnm(1:maxDeg,   0:maxDeg))
if (allocated(Gnm)) deallocate(Gnm)
allocate(Gnm(1:maxDeg,   0:maxDeg))
if (allocated(Enm)) deallocate(Enm)
allocate(Enm(1:maxDeg,   1:maxDeg))
if (allocated(Fnm)) deallocate(Fnm)
allocate(Fnm(1:maxDeg,   1:maxDeg))
if (allocated(Rm)) deallocate(Rm)
allocate(Rm(0:maxDeg))
if (allocated(Im)) deallocate(Im)
allocate(Im(0:maxDeg))
if (allocated(Pn)) deallocate(Pn)
allocate(Pn(0:maxDeg + 1))
if (allocated(Aux1)) deallocate(Aux1)
allocate(Aux1(1:maxDeg+1))
if (allocated(Aux2)) deallocate(Aux2)
allocate(Aux2(1:maxDeg+1))
if (allocated(Aux3)) deallocate(Aux3)
allocate(Aux3(1:maxDeg+1))
if (allocated(Aux4)) deallocate(Aux4)
allocate(Aux4(1:maxDeg+1, 0:maxDeg+1))

read(id_earth,'(a)') (dummy, i=1,2)
do i=1,maxDeg
  do j=0,minval([i,maxOrd])
    read(id_earth,'(2(1x,i2),2(1x,e24.17))') l,m,Cnm(i,j),Snm(i,j)
    Cnm(i,j) = Cnm(i,j)/NORMFACT(i,j)
    Snm(i,j) = Snm(i,j)/NORMFACT(i,j)
  end do
end do

close(id_earth)

! Fill coefficient arrays for Pines algorithm
Anm(:,:) = 0._dk
Dnm(:,:) = 0._dk
Gnm(:,:) = 0._dk
Enm(:,:) = 0._dk
Fnm(:,:) = 0._dk
Anm(0,0) = 1._dk
Anm(1,1) = 1._dk
do n = 1,maxDeg + 1
  Aux1(n) = 2._dk*n + 1._dk
  Aux2(n) = Aux1(n) / (n+1._dk)
  Aux3(n) = n / (n+1._dk)
  do m = 0, n-1
    Aux4(n,m) = (n+m+1._dk)
  end do
end do

! Earth Rotation Rate (revolutions per tropical day)
secsPerSidDay = twopi/omegaE
ERR_constant = secsPerDay/secsPerSidDay

end subroutine INITIALIZE_NSGRAV




function NORMFACT(l,m)
! Description:
!    Normalization factor for the spherical harmonics coefficients and
!    associated Legendre functions:
! 
!    sqrt( (l + m)! / ( (2 - delta_{0,m}) * (2n  + 1) * (n - m)! ) )
! 
! Reference:
!    [1] Montenbruck, O., Gill, E., "Satellite Orbits", p. 58, Springer, 2000.
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! ==============================================================================

! Arguments and function definition
integer,intent(in)  ::  l,m
real(qk)            ::  NORMFACT
! Locals
real(qk)            ::  lr,mr
real(qk)            ::  kron
real(qk)            ::  numer,denom

! ==============================================================================
lr = real(l,qk)
mr = real(m,qk)

numer = gamma(lr + mr + 1._qk)

if (m == 0) then
	kron = 1._qk
else
	kron = 0._qk
end if

denom = (2._qk - kron) * (2._qk*lr + 1._qk) * gamma(lr - mr + 1._qk)

NORMFACT = sqrt(numer/denom)

end function NORMFACT

subroutine GCRFtoITRF(MJD_UTC, Rgcrf, Ritrf, TU, Vgcrf, Vitrf, Agcrf, Aitrf, gcrf_itrf_out)

  ! external subroutines
  use PHYS_CONST, only: UTC2TT

  ! no implicit vars
  implicit none

  ! external functions
  real(dk), external :: iau_SP00
  real(dk), external :: iau_ERA00

  ! Formals
  real(dk), intent(in)            ::  MJD_UTC
  real(dk), intent(in)            ::  Rgcrf(3)
  real(dk), intent(in), optional  ::  Vgcrf(3)
  real(dk), intent(in), optional  ::  Agcrf(3)
  real(dk), intent(out)           ::  Ritrf(3)
  real(dk), intent(out), optional ::  Vitrf(3)
  real(dk), intent(out), optional ::  Aitrf(3)
  real(dk), intent(in)            :: TU
  real(dk), intent(out), optional :: gcrf_itrf_out(3,3)


  ! Locals
  integer   :: eopIdx
  integer   :: imjd, imjd0, imjdcurr, imjdDiff
  real(dk)  :: UTCjd1, UTCjd2, TTjd1, TTjd2, MJD_TT   
  real(dk)  :: UT1_UTC, UT1jd1, UT1jd2
  integer   :: J
  real(dk)  :: Xp, Yp, dx00, dy00, LOD
  real(dk)  :: OMEGAEARTH, omegaE
  real(dk)  :: X, Y, S
  real(dk)  :: BPN_T(3,3), R_T(3,3), W_T(3,3), tempmat1(3,3), tempmat2(3,3), gcrf_itrf(3,3)
  real(dk)  :: era
  real(dk)  :: Rtirs(3), Vtirs(3), tempvec1(3), tempvec2(3), tempvec3(3), tempvec4(3), tempvec5(3), tempvec6(3)
  real(dk)  :: omegavec(3), omegavec2(3)
  real(dk)  :: DAS2R, DMAS2R

  ! constants
  DAS2R  = 4.848136811095359935899141e-6_dk
  DMAS2R = DAS2R / 1e3_dk

  ! get UTC and TT
  MJD_TT  = UTC2TT(MJD_UTC)
  UTCjd1  = 2400000.5_dk
  UTCjd2  = MJD_UTC
  TTjd1   = 2400000.5_dk
  TTjd2   = MJD_TT

  ! get eop data index
  imjdCurr = MJD_UTC
  imjd0 = eop(1,1)

  ! check for unprovided date
  if (imjdCurr < imjd0) then
    write(*,*)imjdCurr,imjd0,'Requested date is before start of Earth Orientation Data'
  endif

  ! index of requested mjd
  imjdDiff = imjdCurr - imjd0
  eopIdx = 1 + imjdDIff

  ! debug
  ! write(*,*)imjd0,imjdCurr,eopIdx,eop(eopIdx,1)

  ! polar motion
  Xp = eop(eopIdx, 2) * DAS2R
  Yp = eop(eopIdx, 3) * DAS2R

  ! CIP offsets wrt IAU 20000
  dx00 = eop(eopIdx, 5) * DMAS2R
  dy00 = eop(eopIdx, 6) * DMAS2R

  ! length of day (convert microseconds to seconds)
  LOD = eop(eopIdx, 7)/1000_dk

  ! compute earth rotation rate
  OMEGAEARTH = 7.292115146706979e-5_dk
  omegaE     = OMEGAEARTH*(1_dk - (LOD/86400_dk))

  ! nondimensionalize rotation rate
  omegaE     = omegaE/TU 

  ! convert UTC to UT1
  ! write(*,*)'eopIdx',eopIdx,eop(eopIdx,1:7)
  UT1_UTC = eop(eopIdx, 4) 
  ! write(*,*)'UT1_UTC',UT1_UTC
  call iau_UTCUT1(UTCjd1, UTCjd2, UT1_UTC, UT1jd1, UT1Jd2, J)

  ! CIP and CIO, IAU 2000A (Represents the PN matrix in Vallado)
  call iau_XYS00B(TTjd1, TTjd2, X, Y, S)

  ! add CIP corrections
  X = X + dx00;
  Y = Y + dy00;

  ! GCRS to CIRS matrix (celestial to intermediate matrix [BPN]' = [N]'[P]'[B]')
  call iau_C2IXYS(X, Y, S, BPN_T)

  ! Earth rotation angle
  era = iau_ERA00(UT1jd1, UT1jd2)

  ! rotation matrix about pole
  R_T(1,1) = 1
  R_T(1,2) = 0
  R_T(1,3) = 0
  R_T(2,1) = 0
  R_T(2,2) = 1
  R_T(2,3) = 0
  R_T(3,1) = 0
  R_T(3,2) = 0
  R_T(3,3) = 1
  call iau_RZ(era, R_T)

  ! Polar motion matrix (TIRS->ITRS, IERS 2003) (W_T matrix)
  call iau_POM00(Xp, Yp, iau_SP00(TTjd1, TTjd2), W_T)

  ! multiply [R]',[BPN]',and [W]'
  call iau_RXR(R_T, BPN_T, tempmat1)
  call iau_RXR(W_T, tempmat1, gcrf_itrf)

  ! rotate position vector
  call iau_RXP(gcrf_itrf, Rgcrf, Ritrf)

  ! output rotation matrix if present
  if (present(gcrf_itrf_out)) then
    gcrf_itrf_out = gcrf_itrf
  endif

  if (present(Vgcrf)) then

    ! compute position and velocity in TIRS frame
    call iau_RXR(R_T, BPN_T, tempmat1)
    call iau_RXP(tempmat1, Rgcrf, Rtirs)

    ! define angular velocity vector
    omegavec(1) = 0
    omegavec(2) = 0
    omegavec(3) = omegaE
    
    ! perform transformation
    call iau_PXP(omegavec, Rtirs, tempvec1)
    call iau_PMP(Vtirs, tempvec1, tempvec2)
    call iau_RXP(W_T, tempvec2, Vitrf)

    if (present(Agcrf)) then

      ! define two times angular velocity vector
      omegavec2(1) = 0
      omegavec2(2) = 0
      omegavec2(3) = 2_dk*omegaE

      ! vector-matrix operations from Vallado

      ! R'*N'*P'*B'
      call iau_RXR(R_T, BPN_T, tempmat1)

      ! R'*N'*P'*B'*Agcrf
      call iau_RXP(tempmat1, Agcrf, tempvec1)

      ! omega x omega
      call iau_PXP(omegavec, omegavec, tempvec2)

      ! omega x omega x Rtirs
      call iau_PXP(tempvec2, Rtirs, tempvec3)

      ! 2*omega x Vtirs
      call iau_PXP(omegavec2, Vtirs, tempvec4)

      ! (omega x omega x Rtirs + 2*omega x Vtirs)
      call iau_PPP(tempvec3, tempvec4, tempvec5)

      ! R'*N'*P'*B'*Agcrf - (omega x omega x Rtirs + 2*omega x Vtirs)
      call iau_PMP(tempvec1, tempvec5, tempvec6)

      ! Aitrf = W'*(R'*N'*P'*B'*Agcrf - (omega x omega x Rtirs + 2*omega x Vtirs))
      call iau_RXP(W_T, tempvec6, Aitrf)

    endif

  endif

end subroutine GCRFtoITRF

subroutine ITRFtoGCRF(MJD_UTC, Ritrf, Rgcrf, TU, Vitrf, Vgcrf, Aitrf, Agcrf, itrf_gcrf_out)

  ! external subroutines
  use PHYS_CONST, only: UTC2TT

  ! no implicit vars
  implicit none

  ! external functions
  real(dk), external :: iau_SP00
  real(dk), external :: iau_ERA00

  ! Formals
  real(dk), intent(in)            ::  MJD_UTC
  real(dk), intent(out)           ::  Rgcrf(3)
  real(dk), intent(out), optional ::  Vgcrf(3)
  real(dk), intent(out), optional ::  Agcrf(3)
  real(dk), intent(in)            ::  Ritrf(3)
  real(dk), intent(in), optional  ::  Vitrf(3)
  real(dk), intent(in), optional  ::  Aitrf(3)
  real(dk), intent(in)            :: TU
  real(dk), intent(out), optional :: itrf_gcrf_out(3,3)


  ! Locals
  integer   :: eopIdx
  integer   :: imjd, imjd0, imjdcurr, imjdDiff
  real(dk)  :: UTCjd1, UTCjd2, TTjd1, TTjd2, MJD_TT   
  real(dk)  :: UT1_UTC, UT1jd1, UT1jd2
  integer   :: J
  real(dk)  :: Xp, Yp, dx00, dy00, LOD
  real(dk)  :: OMEGAEARTH, omegaE
  real(dk)  :: X, Y, S
  real(dk)  :: BPN_T(3,3), R_T(3,3), W_T(3,3), tempmat1(3,3), tempmat2(3,3), itrf_gcrf(3,3)
  real(dk)  :: BPN(3,3), R(3,3), W(3,3)
  real(dk)  :: era
  real(dk)  :: Rtirs(3), Vtirs(3), tempvec1(3), tempvec2(3), tempvec3(3), tempvec4(3), tempvec5(3), tempvec6(3)
  real(dk)  :: omegavec(3), omegavec2(3)
  real(dk)  :: DAS2R, DMAS2R

  ! constants
  DAS2R  = 4.848136811095359935899141e-6_dk
  DMAS2R = DAS2R / 1e3_dk
  ! write(*,*)DAS2R,DMAS2R

  ! get UTC and TT
  MJD_TT  = UTC2TT(MJD_UTC)
  UTCjd1  = 2400000.5_dk
  UTCjd2  = MJD_UTC
  TTjd1   = 2400000.5_dk
  TTjd2   = MJD_TT

  ! get eop data index
  imjdCurr = MJD_UTC
  imjd0 = eop(1,1)

  ! check for unprovided date
  if (imjdCurr < imjd0) then
    write(*,*)imjdCurr,imjd0,'Requested date is before start of Earth Orientation Data'
  endif

  ! index of requested mjd
  imjdDiff = imjdCurr - imjd0
  eopIdx = 1 + imjdDIff

  ! debug
  ! write(*,*)imjd0,imjdCurr,eopIdx,eop(eopIdx,1)

  ! polar motion
  Xp = eop(eopIdx, 2) * DAS2R
  Yp = eop(eopIdx, 3) * DAS2R

  ! CIP offsets wrt IAU 20000
  dx00 = eop(eopIdx, 5) * DMAS2R
  dy00 = eop(eopIdx, 6) * DMAS2R

  ! length of day (convert microseconds to seconds)
  LOD = eop(eopIdx, 7)/1000_dk

  ! compute earth rotation rate
  OMEGAEARTH = 7.292115146706979e-5_dk
  omegaE     = OMEGAEARTH*(1_dk - (LOD/86400_dk))

  ! nondimensionalize rotation rate
  omegaE     = omegaE/TU 

  ! convert UTC to UT1
  UT1_UTC = eop(eopIdx, 4) 
  call iau_UTCUT1(UTCjd1, UTCjd2, UT1_UTC, UT1jd1, UT1Jd2, J)

  ! CIP and CIO, IAU 2000A (Represents the PN matrix in Vallado)
  call iau_XYS00B(TTjd1, TTjd2, X, Y, S)

  ! add CIP corrections
  X = X + dx00;
  Y = Y + dy00;

  ! GCRS to CIRS matrix (celestial to intermediate matrix [BPN]' = [N]'[P]'[B]')
  call iau_C2IXYS(X, Y, S, BPN_T)

  ! transpose
  call iau_TR(BPN_T, BPN)

  ! Earth rotation angle
  era = iau_ERA00(UT1jd1, UT1jd2)

  ! rotation matrix about pole
  R_T(1,1) = 1
  R_T(1,2) = 0
  R_T(1,3) = 0
  R_T(2,1) = 0
  R_T(2,2) = 1
  R_T(2,3) = 0
  R_T(3,1) = 0
  R_T(3,2) = 0
  R_T(3,3) = 1
  call iau_RZ(era, R_T)

  ! transpose
  call iau_TR(R_T, R)

  ! Polar motion matrix (TIRS->ITRS, IERS 2003) (W_T matrix)
  call iau_POM00(Xp, Yp, iau_SP00(TTjd1, TTjd2), W_T)

  ! transpose
  call iau_TR(W_T, W)

  ! multiply [BPN], [R], and [W]
  call iau_RXR(R, W, tempmat1)
  call iau_RXR(BPN, tempmat1, itrf_gcrf)

  ! rotate position vector
  call iau_RXP(itrf_gcrf, Ritrf, Rgcrf)
  ! write(*,*)Rgcrf

  ! output rotation matrix if present
  if (present(itrf_gcrf_out)) then
    itrf_gcrf_out = itrf_gcrf
  endif

  if (present(Vitrf)) then

    ! compute position and velocity in TIRS frame
    call iau_RXR(R_T, BPN_T, tempmat1)
    call iau_RXP(tempmat1, Rgcrf, Rtirs)

    ! define angular velocity vector
    omegavec(1) = 0
    omegavec(2) = 0
    omegavec(3) = omegaE
    
    ! perform transformation
    call iau_PXP(omegavec, Rtirs, tempvec1)
    call iau_PPP(Vtirs, tempvec1, tempvec2)
    call iau_RXR(BPN, R, tempmat1)
    call iau_RXP(tempmat1, tempvec2, Vgcrf)
    ! write(*,*)Vgcrf

    if (present(Aitrf)) then

      ! define two times angular velocity vector
      omegavec2(1) = 0
      omegavec2(2) = 0
      omegavec2(3) = 2_dk*omegaE

      ! vector-matrix operations from Vallado

      ! B*P*N*R
      call iau_RXR(BPN, R, tempmat1)

      ! W*Aitrf
      call iau_RXP(W, Aitrf, tempvec1)

      ! omega x omega
      call iau_PXP(omegavec, omegavec, tempvec2)

      ! omega x omega x Rtirs
      call iau_PXP(tempvec2, Rtirs, tempvec3)

      ! 2*omega x Vtirs
      call iau_PXP(omegavec2, Vtirs, tempvec4)

      ! (omega x omega x Rtirs + 2*omega x Vtirs)
      call iau_PPP(tempvec3, tempvec4, tempvec5)

      ! W*Aitrf + omega x omega x Rtirs + 2*omega x Vtirs
      call iau_PPP(tempvec1, tempvec5, tempvec6)

      ! Aitrf = W'*(R'*N'*P'*B'*Agcrf - (omega x omega x Rtirs + 2*omega x Vtirs))
      call iau_RXP(tempmat1, tempvec6, Agcrf)

    endif

  endif

end subroutine ITRFtoGCRF

subroutine PINES_NSG(GM,RE,rIn,tau,FIn,pot,dPot,vIn)
! Description:
!    Compute the perturbing acceleration, perturbing potential (optional), and
!    time derivative of the perturbing potential in the body-fixed frame
!    (optional), given the position vector wrt the non-spherical body, its
!    gravitational parameter and its potential coefficients.
!    Uses the method described in the Reference to perform the calculation in
!    Cartesian coordinates, thereby avoiding associated Legendre functions and
!    polynomials. The formulation is non-singular everywhere for r > RE, where
!    RE is the radius of the non-spherical body.
! 
!    Adapted from code developed by Hodei Urrutxua (Universidad Rey Juan Carlos,
!    Madrid, Spain) and Claudio Bombardelli (Universidad Politécnica de Madrid,
!    Madrid, Spain).
! 
! Author:
!    Davide Amato
!    The University of Arizona
!    davideamato@email.arizona.edu
! 
! Reference:
!    S. Pines, "Uniform Representation of the Gravitational Potential and its
!    derivatives," AIAA J. 11 (11), pp. 1508-1511, 1973.
! 
! Revisions:
!    180806: First working version of subroutine.
!    181204: Use GMST (IERS 2006 conventions) rather than ERA. Consider Earth
!            Rotation Rate derived from IERS 2006 conventions.
!
! ==============================================================================

! Use associations
use AUXILIARIES, only: MJD0, TU, T2MJD
use PHYS_CONST,  only: delta_JD_MJD
use PHYS_CONST,  only: ERR_IAU06, UTC2TT

! Arguments
real(dk),intent(in)            ::  GM             ! Gravitational parameter
real(dk),intent(in)            ::  RE             ! Equatorial radius
real(dk),intent(in)            ::  rIn(1:3)       ! Position in the inertial frame
real(dk),intent(in), optional  ::  vIn(1:3)       ! Velocity in the inertial frame
real(dk),intent(in)            ::  tau            ! Physical time
real(dk),intent(out),optional  ::  FIn(1:3)  ! Perturbing acceleration in the inertial frame
real(dk),intent(out),optional  ::  pot   ! Perturbing potential
real(dk),intent(out),optional  ::  dPot  ! Time derivative of the potential in body-fixed frame
! Locals
real(dk)  ::  rNorm,rNormSq  ! Norm of position vector and square
real(dk)  ::  s,t,u  ! Direction cosines in the body-fixed frame
real(dk)  ::  rho    ! = equatorial radius / r 
real(dk)  ::  F(1:3),a1,a2,a3,a4  ! Acceleration and its components 
integer   ::  n,m    ! Harmonic indices
logical   ::  skip_EFG
! GMST-related quantities
real(dk)  ::  MJD_UTC, MJD_TT      ! UTC and TT dates
real(dk)  ::  GMST,cosGMST,sinGMST ! GMST and its trig functions
real(dk)  ::  ERR, ERR_nd          ! Earth Rotation Rate [rad/s, -]            

! SOFA routines
real(dk) :: iau_GMST06

! gcrf -> itrf conversion
real(dk) :: Rgcrf(3), Ritrf(3)
real(dk) :: Vgcrf(3), Vitrf(3)
real(dk) :: Agcrf(3), Aitrf(3)
real(dk) :: dummy
real(dk) :: gcrf_itrf(3,3)
real(dk) :: itrf_gcrf(3,3)


! ==============================================================================

rNormSq = dot_product(rIn,rIn)
rNorm = sqrt(rNormSq)

! ==============================================================================
! 01. Transform from inertial to body-fixed frame (gcrf -> itrf using IAU2000 CIO conventions)
! ==============================================================================

! get UTC and TT
MJD_UTC = T2MJD(tau)
MJD_TT  = UTC2TT(MJD_UTC)


! required formals:
! UTCjd1, UTCjd2, dx00, dy00, UT1_UTC, TT_UT1, xp, yp

! convert gcrf to itrf
Rgcrf = rIn
dummy = 1
call GCRFtoITRF(MJD_UTC, Rgcrf, Ritrf, TU, gcrf_itrf_out=gcrf_itrf)
u = Ritrf(3)/rNorm
s = Ritrf(1)/rNorm
t = Ritrf(2)/rNorm
! write(*,*)Rgcrf,Ritrf

! transpose
call iau_TR(gcrf_itrf, itrf_gcrf)

!  --------------------- gmst rot only ---------------------
! u = rIn(3)/rNorm

! ! Greenwich Mean Sidereal Time (IAU 2006 conventions)
! GMST = iau_GMST06 ( delta_JD_MJD, MJD_UTC, delta_JD_MJD, MJD_TT )
! cosGMST = cos(GMST); sinGMST = sin(GMST)

! ! Rotate equatorial components of rIn to get direction cosines in the body-fixed frame
! s = (rIn(1) * cosGMST + rIn(2) * sinGMST)/rNorm
! t = (rIn(2) * cosGMST - rIn(1) * sinGMST)/rNorm
!  --------------------- gmst rot only ---------------------

! write(*,*)'old:',u,s,t

! write(*,*)tau,MJD_UTC

rho = RE/rNorm

! ==============================================================================
! 02. Fill in coefficient matrices and auxiliary vectors
! ==============================================================================

! Fill A Matrix
Anm(0,0) = 1._dk
Anm(1,1) = 1._dk
Anm(1,0) = u
do n = 1, gdeg + 1
  Anm(n+1,n+1) = Aux1(n) * Anm(n,n) ! Fill the diagonal
  Anm(n+1,0) = Aux2(n) * u * Anm(n,0) - Aux3(n) * Anm(n-1,0) ! Fill the 1st column
  Anm(n+1,n) = u * Anm(n+1,n+1)    ! Fill the subdiagonal

end do
do n = 2, gdeg + 1 ! Fill remaining elements
  do m = 0, n - 2
    Anm(n+1,m+1) = Aux4(n,m) * Anm(n,m) + u * Anm(n,m+1)
  end do

end do

! Fill R, I, and P vectors
Rm(0) = 1._dk
Im(0) = 0._dk
Pn(0) = GM / rNorm
Pn(1) = rho * Pn(0)
do n = 1, gdeg
  Rm(n)   = s  * Rm(n-1) - t * Im(n-1)
  Im(n)   = s  * Im(n-1) + t * Rm(n-1)
  Pn(n+1) = rho * Pn(n)

end do

! Fill D, E, and F matrices
skip_EFG = present(pot) .and. .not.(present(FIn)) .and. .not.(present(dPot))
do m = 1, gord
  do n = m, gdeg
    Dnm(n,m) = Cnm(n,m)*Rm(m)   + Snm(n,m)*Im(m)
    if (.not.(skip_EFG)) then
      Enm(n,m) = Cnm(n,m)*Rm(m-1) + Snm(n,m)*Im(m-1)
      Fnm(n,m) = Snm(n,m)*Rm(m-1) - Cnm(n,m)*Im(m-1)

    end if
  end do

end do
do n = 1, gdeg
  Dnm(n,0) = Cnm(n,0)*Rm(0)  !+ S(n,0)*I(0) = 0
end do

! ==============================================================================
! 03. Perturbing potential
! ==============================================================================

if (present(pot)) then
  pot = 0._dk
  do m = 1, gord
    do n = m, gdeg
      pot = pot + Pn(n) * Anm(n,m) * Dnm(n,m)
		
    end do
	
  end do
  do n = 1, gdeg
  pot = pot + Pn(n) * Anm(n,0) * Dnm(n,0)
	
  end do
  pot = -pot  ! Change the sign to get the potential
end if

! ==============================================================================
! 04. Perturbing acceleration
! ==============================================================================

if (present(FIn)) then
  a1 = 0._dk; a2 = 0._dk; a3 = 0._dk; a4 = 0._dk
  do m = 1, gord
    do n = m, gdeg
      a1 = a1 + Pn(n+1) * Anm(n,m) * m * Enm(n,m)
      a2 = a2 + Pn(n+1) * Anm(n,m) * m * Fnm(n,m)
      a3 = a3 + Pn(n+1) * Anm(n,m+1)   * Dnm(n,m)
      a4 = a4 - Pn(n+1) * Anm(n+1,m+1) * Dnm(n,m)
    end do
  end do
  do n = 1, gdeg
    a3 = a3 + Pn(n+1) * Anm(n,1)     * Dnm(n,0)
    a4 = a4 - Pn(n+1) * Anm(n+1,1)   * Dnm(n,0)

  end do
  F = [a1, a2, a3] + [s, t, u] * a4
  F = F / RE

  ! Transform to inertial frame
  Aitrf = F
  ! write(*,*)Aitrf,Vitrf

  ! call ITRFtoGCRF(MJD_UTC, Aitrf, Agcrf, TU)

  call iau_RXP(itrf_gcrf, Aitrf, Agcrf)
  Fin(1) = Agcrf(1)
  Fin(2) = Agcrf(2)
  Fin(3) = Agcrf(3)
  ! write(*,*)Fin(1),Fin(2),Fin(3)

  !  --------------------- gmst rot only ---------------------
  ! FIn(1) = F(1)*cosGMST - F(2)*sinGMST
  ! FIn(2) = F(1)*sinGMST + F(2)*cosGMST
  ! FIn(3) = F(3)
  !  --------------------- gmst rot only ---------------------
end if


! ==============================================================================
! 05. Time derivative of potential in body-fixed frame
! ==============================================================================

if(present(dPot)) then
  dPot = 0._dk
  ERR = ERR_IAU06(0._dk, MJD_TT)
  ERR_nd = ERR / TU
  do m = 1, gord
    do n = m, gdeg
      Gnm(n,m) = m * ( t * Enm(n,m) - s * Fnm(n,m) )
      Gnm(n,m) = ERR_nd * Gnm(n,m)
      dPot = dPot + Pn(n) * Anm(n,m) * Gnm(n,m)
      
    end do
    
  end do
  dPot = -dPot
end if

end subroutine PINES_NSG




end module NSGRAV
