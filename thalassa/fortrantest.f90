program reductionTest

    ! external modules
    use KINDS,  only: dk
    use NSGRAV, only: INITIALIZE_EOP, ITRFtoGCRF, GCRFtoITRF

    ! no implict variables
    implicit none

    ! local variables
    real(dk) :: MJD_UTC
    real(dk) :: result_t_to_g(3,3)
    real(dk) :: result_g_to_t(3,3)
    real(dk) :: RitrfIn(3), RgcrfIn(3)
    real(dk) :: RitrfOut(3), RgcrfOut(3)
    real(dk) :: gcrf_itrf(3,3), itrf_gcrf(3,3)
    real(dk) :: gcrf_itrf_err(3,3), itrf_gcrf_err(3,3)
    real(dk) :: TU = 1_dk

    ! MJD UTC
    MJD_UTC = 54195.500000

    ! test vals for position vector
    RitrfIn(1) = 5000
    RitrfIn(2) = 5000
    RitrfIn(3) = 5000
    RgcrfIn(1) = 5000
    RgcrfIn(2) = 5000
    RgcrfIn(3) = 5000

    ! expected results for itrf  to gcrf
    result_g_to_t(1,1) = 0.973104317697512_dk
    result_g_to_t(1,2) = 0.230363826239227_dk
    result_g_to_t(1,3) = -0.000703163482268_dk
    result_g_to_t(2,1) = -0.230363800456136_dk
    result_g_to_t(2,2) = 0.973104570632777_dk
    result_g_to_t(2,3) = 0.000118545366806_dk
    result_g_to_t(3,1) = 0.000711560162777_dk
    result_g_to_t(3,2) = 0.000046626403835_dk
    result_g_to_t(3,3) = 0.999999745754024_dk

    ! transpose
    call iau_TR(result_g_to_t, result_t_to_g)

    ! call subroutines
    call INITIALIZE_EOP('./data/10_FINALS.DATA_IAU2000_V2013_0110.txt')
    call GCRFtoITRF(MJD_UTC, RgcrfIn, RitrfOut, TU, gcrf_itrf_out=gcrf_itrf)
    call ITRFtoGCRF(MJD_UTC, RitrfIn, RgcrfOut, TU, itrf_gcrf_out=itrf_gcrf)

    ! calculate gcrf to itrf errors
    gcrf_itrf_err(1,1) = result_g_to_t(1,1) - gcrf_itrf(1,1)
    gcrf_itrf_err(1,2) = result_g_to_t(1,2) - gcrf_itrf(1,2)
    gcrf_itrf_err(1,3) = result_g_to_t(1,3) - gcrf_itrf(1,3)
    gcrf_itrf_err(2,1) = result_g_to_t(2,1) - gcrf_itrf(2,1)
    gcrf_itrf_err(2,2) = result_g_to_t(2,2) - gcrf_itrf(2,2)
    gcrf_itrf_err(2,3) = result_g_to_t(2,3) - gcrf_itrf(2,3)
    gcrf_itrf_err(3,1) = result_g_to_t(3,1) - gcrf_itrf(3,1)
    gcrf_itrf_err(3,2) = result_g_to_t(3,2) - gcrf_itrf(3,2)
    gcrf_itrf_err(3,3) = result_g_to_t(3,3) - gcrf_itrf(3,3)

    ! calculate itrf to gcrf errors
    itrf_gcrf_err(1,1) = result_t_to_g(1,1) - itrf_gcrf(1,1)
    itrf_gcrf_err(1,2) = result_t_to_g(1,2) - itrf_gcrf(1,2)
    itrf_gcrf_err(1,3) = result_t_to_g(1,3) - itrf_gcrf(1,3)
    itrf_gcrf_err(2,1) = result_t_to_g(2,1) - itrf_gcrf(2,1)
    itrf_gcrf_err(2,2) = result_t_to_g(2,2) - itrf_gcrf(2,2)
    itrf_gcrf_err(2,3) = result_t_to_g(2,3) - itrf_gcrf(2,3)
    itrf_gcrf_err(3,1) = result_t_to_g(3,1) - itrf_gcrf(3,1)
    itrf_gcrf_err(3,2) = result_t_to_g(3,2) - itrf_gcrf(3,2)
    itrf_gcrf_err(3,3) = result_t_to_g(3,3) - itrf_gcrf(3,3)

    ! output
    write(*,*)'GCRF to ITRF Errors'
    write(*,*)gcrf_itrf_err(1,1:3)
    write(*,*)gcrf_itrf_err(2,1:3)
    write(*,*)gcrf_itrf_err(3,1:3)
    write(*,*)'ITRF to GCRF Errors'
    write(*,*)itrf_gcrf_err(1,1:3)
    write(*,*)itrf_gcrf_err(2,1:3)
    write(*,*)itrf_gcrf_err(3,1:3)

end program reductionTest