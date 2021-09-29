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
    ! ---------------------------------------------------------------------
    ! Model Parameters
    INTEGER, PARAMETER          :: nJ = 66               ! Lifespan of the agents
    INTEGER                     :: j                          ! Index of age group
    INTEGER, PARAMETER          :: J_R     = 46               ! Retirement age
    DOUBLE PRECISION, PARAMETER :: n       = 0.011d0          ! Population growth rate
    DOUBLE PRECISION, PARAMETER :: a_1     = 0d0              ! Initial assets holding for newborns
    DOUBLE PRECISION, PARAMETER :: THETA   = 0.11d0           ! Labor income tax rate
    DOUBLE PRECISION, PARAMETER :: GAMMA   = 0.42d0           ! Utillity weight on consumption
    DOUBLE PRECISION, PARAMETER :: SIGMA   = 2.0d0            ! coefcient of relative risk aversion
    DOUBLE PRECISION, PARAMETER :: ALPHA   = 0.36d0           ! Capital share in production
    DOUBLE PRECISION, PARAMETER :: DELTA   = 0.06d0           ! Capital depreciation rate
    DOUBLE PRECISION, PARAMETER :: BETA    = 0.97d0           ! Discount fact
    DOUBLE PRECISION            :: ez(J_R-1)                  ! Age eficiency profiles
    ! -----------------------------------------------------------------------
	! ****************************GRID SET**********************************
	! -----------------------------------------------------------------------
    ! Parameters regarding stochastic processe
    ! Set up for discritizing the state space (Productivity Grid)
    INTEGER                     :: i_z                        ! Indexes for productivity                 ! Number of points in the grid
    DOUBLE PRECISION, PARAMETER :: z_H     = 3.0d0            ! Idiosyncratic productivity High
    DOUBLE PRECISION, PARAMETER :: z_L     = 0.5d0            ! Idiosyncratic productivity Low
    INTEGER, PARAMETER          :: nZ      = 2                ! Number of points in the grid for z
    DOUBLE PRECISION            :: grid_Z(nZ)                 ! Grid for productivity
    DOUBLE PRECISION, PARAMETER :: p_H     = 0.2037d0         ! Probability of z_H at birth
    DOUBLE PRECISION, PARAMETER :: p_L     = 0.7963d0         ! Probability of z_L at birth
    DOUBLE PRECISION, PARAMETER :: PI_HH   = 0.9261d0         ! Probability of transision from z_H to z_H
    DOUBLE PRECISION, PARAMETER :: PI_HL   = 1.0d0 - PI_HH    ! Probability of transision from z_H to z_L
    DOUBLE PRECISION, PARAMETER :: PI_LL   = 0.0739d0         ! Probability of transision from z_L to z_L
    DOUBLE PRECISION, PARAMETER :: PI_LH   = 1.0d0 - PI_LL    ! Probability of transision from z_L to z
	! Set up for discritizing the state space (Asset Grid)
    INTEGER                     :: i_A, i_Anext               ! Indexes for the asset grid
    INTEGER, PARAMETER          :: nA      = 1000             ! Number of points in the grid for assets
    DOUBLE PRECISION            :: grid_A(nA)                 ! Grid for assets
    DOUBLE PRECISION, PARAMETER :: A_min   = 0d0              ! Minimum asset level
    DOUBLE PRECISION, PARAMETER :: A_max   = 75d0             ! Maximum asset level
    DOUBLE PRECISION, PARAMETER :: sA      = dble(nA) - 1d0   ! Asset grid spacing factor
    DOUBLE PRECISION, PARAMETER :: step_A  = (A_max-A_min)/sA ! Asset grid step
    ! Model Variables
    DOUBLE PRECISION            :: b                          !  Pension benefits
    DOUBLE PRECISION            :: r                          !  Interest rate
    DOUBLE PRECISION            :: w                          !  Wage
    ! Global parameters for optimization
    ! Variables for dynamic programming
    DOUBLE PRECISION            :: cand_max                   ! Initialize candite maximum
    DOUBLE PRECISION            :: v_current                  ! Initialize current value of the value fucntion
    DOUBLE PRECISION            :: c_current                  ! Initialize current value of the consumption
    DOUBLE PRECISION            :: a_current                  ! Initialize current value of the asset
    DOUBLE PRECISION            :: a_next                     ! Initialize next value of the asset
    DOUBLE PRECISION            :: a_maxim                    ! Initialize canditate maximizer asset holding level
    INTEGER                     :: a_ind_maxim                ! Initialize canditate maximizer asset holding index
    DOUBLE PRECISION            :: l_opt                      ! Initialize optimal labor supply
    DOUBLE PRECISION            :: v_next                     ! Initialize continuation value of the value function
    DOUBLE PRECISION            :: u_current                  ! Initialize current value of utility
    DOUBLE PRECISION            :: Pr(nZ)                     ! Initialize probability of productivity transition     
    DOUBLE PRECISION            :: e                          ! Initialize the age efficiency
    ! Allocating space for Policy Functions
    DOUBLE PRECISION            :: pf_c(nA, nZ, nJ)
    DOUBLE PRECISION            :: pf_A(nA, nZ, nJ)
    INTEGER                     :: pf_A_ind(nA, nZ, nJ)
    DOUBLE PRECISION            :: pf_v(nA, nZ, nJ)
    INTEGER                     :: i_stat
    INTEGER                     :: iMaxThreads
    ! Variables for paralellization

    ! -----------------------------------------------------------------------
    end module params_grid ! end of module

program conesa_krueger
    use params_grid
    implicit none

    ! Begin Computational Timer
    INTEGER                     :: beginning, rate! ,end
    ! Variables for reading parameters from comand line
    CHARACTER(100)              :: 				r_char
    CHARACTER(100)              :: 				w_char
    CHARACTER(100)              :: 				b_char

    call system_clock(beginning, rate)

    call read_ez() ! Read in the ez array
    
    ! Read in the parameters from the command line
    IF (COMMAND_ARGUMENT_COUNT() < 1) THEN
        ! If no command line arguments are given, use the default values
        ! default values correspond to the values on Questions 1 and 2
        r = 0.05d0
        w = 1.05d0
        b = 0.2d0
    ELSE IF (COMMAND_ARGUMENT_COUNT() < 3 ) THEN
        ! If less than three command line arguments are given throw an error and stop
        WRITE(*,*) 'ERROR: WAS EXPECTING 3 COMAND LINE ARGUMENTS, ONLY GOT', COMMAND_ARGUMENT_COUNT()
        WRITE(*,*) 'IF YOU WANT TO USE THE DEFAULT VALUES, RUN THE PROGRAM WITHOUT ANY ARGUMENTS'
        WRITE(*,*) 'IF YOU WANT TO SUPPLY ARGUMENTS TO THE PROGRAM USE: ./program $r $w $b'
        WRITE(*,*) 'EXITING...'
        STOP
    END IF

    ! Read User supplied values for r, w, b
    CALL GET_COMMAND_ARGUMENT(1,r_char)
    READ(r_char,*)r

    CALL GET_COMMAND_ARGUMENT(2,w_char)
    READ(w_char,*)w

    CALL GET_COMMAND_ARGUMENT(3,b_char)
    READ(b_char,*)b

    ! write(*,*) 'r = ', r, ' w = ', w, ' b = ', b

    call housekeeping()                 ! Set up the grids and allocate space for the policy functions
    
    call V_Func_Ret()                   ! Solve for the value function for retirees

    call V_Func_Work()                   ! Solve for the value function for working age agents

    call coda()                         ! Write resutls to file
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

    open(action='read', file='./PS3/FortranCode/ef.csv', iostat=rc, newunit=fu)
    
    if (rc /= 0) stop 'Unable to open file ez.csv'

    do i = 1, J_R - 1
        read(fu, *, iostat=rc) ez(i)
        if (rc /= 0) stop 'Error reading file ez.csv'
    end do

    close(fu)

    end subroutine read_ez ! end of subroutine

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : housekeeping
!
! description : Initializes Grids and Policy Functions
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
subroutine housekeeping()

    use params_grid

    implicit none

    ! Discretize the state space (Assets)
    do i_A = 1, nA
        grid_A(i_A) = A_min + step_A * dble(i_A-1)
    end do

    ! Discretize the state space (Productivity)
    do i_Z = 1, nZ
        if (i_Z == 1) then
            grid_z(i_Z) = z_H
        else
            grid_z(i_Z) = z_L
        end if
    end do
    
    ! Initialize Policy Functions
    do i_A = 1,nA
        do i_Z = 1,nZ
            do j = 1, nJ
                if ( j < nJ ) then
                    pf_c(i_A, i_Z, j) = 0d0
                    pf_A(i_A, i_Z, j) = 0d0
                    pf_A_ind(i_A, i_Z, j) = 0
                    pf_v(i_A, i_Z, j) = 0d0
                else
                    pf_c(i_A, i_Z, j) = grid_A(i_A)*(1 + r) + b
                    pf_A(i_A, i_Z, j) = 0d0
                    pf_A_ind(i_A, i_Z, j) = 1
                    pf_v(i_A, i_Z, j) = pf_c(i_A, i_Z, j)**(GAMMA * (1 - SIGMA)) / (1 - SIGMA)
                end if 
            end do
        end do
    end do

    end subroutine housekeeping ! end of subroutine

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : V_Func_Ret
!
! description : This subroutine computes the value function for 
! retirement using backward induction.
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
subroutine V_Func_Ret()
    use params_grid

    implicit none

    do j = nJ-1, J_R, -1 ! Loop backwards over the age groups until the retirement age
        do i_A = 1, nA   ! Loop over the asset grid
            cand_max = -100000000d0 ! A large negative number for candidate to maximum
            a_current = grid_A(i_A)
            do i_Anext = 1, nA ! Loop over all possible choises of the next asset level
                a_next = grid_A(i_Anext)
                c_current = (1 + r) * a_current + b - a_next
                if (c_current > 0d0) then ! Check if consumption is positive
                    u_current = c_current**(1 - GAMMA) / (1 - GAMMA)
                    v_current = u_current + BETA * pf_v(i_Anext, 1, j+1)
                    if (v_current > cand_max) then  ! Update the candidate maximizer
                        cand_max = v_current
                        a_maxim = a_next
                        a_ind_maxim = i_Anext
                    end if
                else
                    exit ! If consumption is negative, then stop searching
                end if
            end do ! End of loop over i_Anext
            pf_v(i_A, 1, j) = cand_max
            pf_v(i_A, 2, j) = cand_max
            pf_A(i_A, 1, j) = a_max
            pf_A(i_A, 2, j) = a_max
            pf_c(i_A, 1, j) = c_current
            pf_c(i_A, 2, j) = c_current
            pf_A_ind(i_A, 1, j) = a_ind_maxim
            pf_A_ind(i_A, 2, j) = a_ind_maxim
        end do ! End of loop over i_A
    end do ! End of loop over j
    
end subroutine V_Func_Ret ! end of subroutine
    
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : V_Func_Work
!
! description : This subroutine computes the value function for 
! working age using backward induction.
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
subroutine V_Func_Work()
    use params_grid
    

    implicit none

    do j = J_R - 1, 1, -1 ! Loop backwards over the age groups from last working age to the birth age
        do i_Z = 1, nZ   ! Loop over the productivity grid
            e = ez(j) * grid_z(i_Z) ! Compute the effective productivity for this age group 
            ! *********************
            ! Define today's shock
            ! *********************
            if (i_Z == 1) then
                Pr(1) = PI_HH
                Pr(2) = PI_HL
            else
                Pr(1) = PI_LH
                Pr(2) = PI_LL
            end if
            do i_A = 1, nA   ! Loop over the asset grid
                a_ind_maxim = -100
                a_maxim = -100d0 ! Initialize the value function for this age group
                cand_max = -100d0 ! A large negative number for candidate to maximum
                ! a_maxim  = -100000000d0 ! A large negative number for the asset level that maximizes the value function
                a_current = grid_A(i_A) ! Compute the current asset level
                do i_Anext = 1, nA ! Loop over all possible choises of the next asset level

                    ! *********************
                    ! Compute agent choices
                    ! *********************
                    a_next = grid_A(i_Anext)
                    l_opt = (GAMMA*(1-THETA)*e*w-(1-GAMMA)*((1+r)*a_current-a_next))/((1-THETA)*e*w) ! Compute the labor supply
                    c_current = w*(1-THETA)*e*l_opt+(1+r)*a_current-a_next ! Compute the consumption
                    ! *******************************************************
                    ! Compute the fesibility of labor and consumption choices
                    ! *******************************************************
                    if (c_current < 0d0) then ! If consumption is negative, then go to the next iteration
                        CYCLE
                    end if 
                    if (l_opt < 0d0) then ! If labor supply is negative, then go to the next iteration
                        CYCLE
                    end if
                    if (l_opt > 1.0d0) then ! If labor supply is greater than the maximum, then go to the next iteration
                        CYCLE
                    end if
                    ! **************************
                    ! Compute the value function
                    ! **************************
                    u_current = ( (c_current**GAMMA ) * ( (1 - l_opt) ** (1 - GAMMA)) ** (1 - SIGMA)) / (1 - SIGMA) !  Compute the utility
                    v_next = pf_v(i_Anext, 1, j+1) * Pr(1) + pf_v(i_Anext, 2, j+1) * Pr(2) ! Compute the expected value of next period
                    v_current = u_current + BETA * v_next

                    ! Check if the current value is greater than the previous one and update the value function
                    if (v_current > cand_max) then 
                        cand_max = v_current
                        a_maxim = a_next
                        a_ind_maxim = i_Anext
                    end if
            
                    pf_v(i_A, i_Z, j) = cand_max
                    pf_A(i_A, i_Z, j) = a_maxim
                    pf_A_ind(i_A, i_Z, j) = a_ind_maxim
                    pf_c(i_A, i_Z, j) = c_current
                end do ! End of loop over i_Anext 
            end do ! End of loop over i_A
        end do ! End of loop over i_Z
    end do ! End of loop over j

    end subroutine V_Func_Work ! end of subroutine

! ************************************************************************
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : coda
!
! description : Writes results to .csv file
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
subroutine coda()

    use params_grid
    INTEGER                         :: rc ! Identifier for error checking

    open(unit=2, file='./PS3/FortranCode/results.csv', status='replace', action='write', iostat=rc)
    if (rc /= 0) stop 'Unable to open file results.csv'
    200 format(i2,2x, f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,i4,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x)
    
    do j = 1, nJ
        do i_Z = 1, nZ
            do i_A = 1, nA
                write(2,200)j,grid_Z(i_Z),grid_A(i_A),pf_c(i_A, i_Z, j),pf_A(i_A, i_Z, j),pf_A_ind(i_A, i_Z, j),pf_v(i_A, i_Z, j)
            end do
        end do
    end do

return
end subroutine coda ! end of subroutine