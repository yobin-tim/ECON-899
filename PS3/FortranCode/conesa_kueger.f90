! ------------------------------------------------------------------------
! module : params_grid
!
! Description : This module will form the foudation for our program. In it
! we will allocate space for all paramaters used in this program and set up the
! grids to create a discritized state space
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

module params_grid

    implicit none
    
    ! -----------------------------------------------------------------------
    ! *******************DECLARATION OF PARAMETERS AND VARIABLES*************
    ! -----------------------------------------------------------------------
    
    ! Model Parameters
    INTEGER, PARAMETER          :: N_final     = 66               ! Lifespan of the agents
    INTEGER, PARAMETER          :: J_R   = 46               ! Retirement age
    DOUBLE PRECISION, PARAMETER :: n     = 0.011d0          ! Population growth rate
    DOUBLE PRECISION, PARAMETER :: a_1   = 0d0              ! Initial assets holding for newborns
    DOUBLE PRECISION, PARAMETER :: THETA = 0.11d0           ! Labor income tax rate
    DOUBLE PRECISION, PARAMETER :: z_H   = 1.0d0            ! Idiosyncratic productivity High
    DOUBLE PRECISION, PARAMETER :: z_L   = 0.5d0            ! Idiosyncratic productivity Low
    DOUBLE PRECISION, PARAMETER :: GAMMA = 0.42d0           ! Utillity weight on consumption
    DOUBLE PRECISION, PARAMETER :: SIGMA = 2.0d0            ! coefcient of relative risk aversion
    DOUBLE PRECISION, PARAMETER :: ALPHA = 0.36d0           ! Capital share in production
    DOUBLE PRECISION, PARAMETER :: DELTA = 0.06d0           ! Capital depreciation rate
    DOUBLE PRECISION, PARAMETER :: BETA  = 0.97d0           ! Discount factor

    ! Model Probabilities
    DOUBLE PRECISION, PARAMETER :: p_H   = 0.2037d0         ! Probability of z_H at birth
    DOUBLE PRECISION, PARAMETER :: p_L   = 0.7963d0         ! Probability of z_L at birth
    DOUBLE PRECISION, PARAMETER :: PI_HH = 0.9261d0         ! Probability of transision from z_H to z_H
    DOUBLE PRECISION, PARAMETER :: PI_HL = 1.0d0 - PI_HH    ! Probability of transision from z_H to z_L
    DOUBLE PRECISION, PARAMETER :: PI_LL = 0.0739d0         ! Probability of transision from z_L to z_L
    DOUBLE PRECISION, PARAMETER :: PI_LH = 1.0d0 - PI_LL    ! Probability of transision from z_L to z_H
    
    ! Model Grids
    DOUBLE PRECISION           :: ez(J_R-1)                ! Age eficiency profiles
    ! Model Variables
    
    DOUBLE PRECISION          :: b                !  Pension benefits
    DOUBLE PRECISION          :: r                !  Interest rate
    DOUBLE PRECISION          :: w                !  Wage
    
    ! -----------------------------------------------------------------------
	! ****************************GRID SET**********************************
	! -----------------------------------------------------------------------
	! Set up for discritizing the state space (Capital Grid)
    
    ! Set up for discritizing the state space (Employment Grid)

    ! Set up for discritizing the state space (Age Grid)

    end module params_grid ! end of module

program conesa_krueger
    use params_grid
    implicit none

    INTEGER :: i

    call read_ez()

    do i = 1, J_R-1
        write(*,*) "ez(",i,") = ", ez(i)
    end do

end program

! ------------------------------------------------------------------------
! subroutine : read_ez
!
! Description : This module reads ez.txt the file containing the deterministic
! age efficiency profile.
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine read_ez()
    use params_grid
    implicit none

    integer        :: fu, rc, i

    open(action='read', file='/home/mitch/Work/ECON-899/PS3/FortranCode/ef.csv', iostat=rc, newunit=fu)
    
    if (rc /= 0) stop 'Unable to open file ez.txt'

    do i = 1, J_R - 1
        read(fu, *, iostat=rc) ez(i)
        if (rc /= 0) stop 'Error reading file ez.txt'
    end do

    close(fu)

    end subroutine read_ez ! end of subroutine