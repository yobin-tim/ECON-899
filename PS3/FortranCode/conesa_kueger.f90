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
    ! ------------------------------------------------------------------------
    ! First we define a variable to determine if we are solving the whole problem or just the Value Function
    INTEGER                     :: solve_whole_problem        ! 1 = solve whole problem, 0 = solve value function
    ! Model Parameters
    INTEGER, PARAMETER          :: nJ      = 66               ! Lifespan of the agents
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
    DOUBLE PRECISION, PARAMETER :: PI_LL   = 0.9811d0         ! Probability of transision from z_L to z_L
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
    INTEGER                     :: min_ind_search             ! Index of the minimum value in the search space
    DOUBLE PRECISION            :: cand_max                   ! Initialize candite maximum
    DOUBLE PRECISION            :: v_current                  ! Initialize current value of the value fucntion
    DOUBLE PRECISION            :: c_current                  ! Initialize current value of the consumption
    DOUBLE PRECISION            :: c_cand                     ! Initialize candidate value of the consumption
    DOUBLE PRECISION            :: a_current                  ! Initialize current value of the asset
    DOUBLE PRECISION            :: a_next                     ! Initialize next value of the asset
    DOUBLE PRECISION            :: a_maxim                    ! Initialize canditate maximizer asset holding level
    INTEGER                     :: a_ind_maxim                ! Initialize canditate maximizer asset holding index
    DOUBLE PRECISION            :: l_opt                      ! Initialize optimal labor supply
    DOUBLE PRECISION            :: l_opt_cand                 ! Initialize optimal labor supply candidate
    DOUBLE PRECISION            :: v_next                     ! Initialize continuation value of the value function
    DOUBLE PRECISION            :: u_current                  ! Initialize current value of utility
    DOUBLE PRECISION            :: Pr(nZ)                     ! Initialize probability of productivity transition     
    DOUBLE PRECISION            :: e                          ! Initialize the age efficiency
    ! Allocating space for Policy Functions
    DOUBLE PRECISION            :: pf_c(nA, nZ, nJ)           ! Policy function for consumption
    DOUBLE PRECISION            :: pf_A(nA, nZ, nJ)           ! Policy function for asset
    INTEGER                     :: pf_A_ind(nA, nZ, nJ)       ! Policy function for asset index
    DOUBLE PRECISION            :: pf_v(nA, nZ, nJ)           ! Policy function for value
    DOUBLE PRECISION            :: pf_l(nA, nZ, nJ)           ! Policy function for labor
    
    ! Allocating space for the steady state distribution
    DOUBLE PRECISION            :: mu(nJ)
    DOUBLE PRECISION            :: F_SS(nA, nZ, nJ)

    ! Allocating space for the aggreate capital and labor levels
    DOUBLE PRECISION            :: K_SS = 0                   ! Steady state aggregate capital
    DOUBLE PRECISION            :: L_SS = 0                   ! Steady state aggregate labor
    DOUBLE PRECISION            :: K                          ! Initialize initial guess of capital
    DOUBLE PRECISION            :: L                          ! Initialize initial guess of labor
    
    ! Parameters for market clearing iteration
    DOUBLE PRECISION            :: LAMBDA     = 0.7d0         ! Adjustment parameter for the market clearing iteration
    DOUBLE PRECISION            :: ERR        = 100d0         ! Initialize error
    DOUBLE PRECISION            :: TOL        = 1.0d-2        ! Tolerance for the market clearing iteration
    INTEGER                     :: MAX_ITER   = 200          ! Maximum number of iterations for the market clearing iteration
    INTEGER                     :: COUNT_ITER = 0             ! Maximum number of iterations for the market clearing iteration
    

    INTEGER                     :: i_stat
    INTEGER                     :: iMaxThreads
    ! Variables for paralellization

    ! -----------------------------------------------------------------------
    end module params_grid ! end of module

program conesa_krueger
    use params_grid
    implicit none

    ! Begin Computational Timer
    INTEGER                     ::  beginning, rate ,end ! for timing
    DOUBLE PRECISION            ::  m = 0 ! To store mass of retirment
    ! Variables for reading parameters from comand line
    CHARACTER(100)              ::  solve_whole_problem_char      
    CHARACTER(100)              ::  r_char
    CHARACTER(100)              ::  w_char
    CHARACTER(100)              ::  b_char
    CHARACTER(100)              ::  K_char
    CHARACTER(100)              ::  L_char

    call system_clock(beginning, rate)

    call read_ez() ! Read in the ez array
    
    ! Read in the parameters from the command line
    IF (COMMAND_ARGUMENT_COUNT() < 1) THEN
        ! If no command line arguments are given, use the default values
        ! default values correspond to the values on Questions 1 and 2
        solve_whole_problem = 0
        r = 0.05d0
        w = 1.05d0
        b = 0.2d0
    ELSE 
        CALL GET_COMMAND_ARGUMENT(1,solve_whole_problem_char)
        READ(solve_whole_problem_char,*)solve_whole_problem
        IF (solve_whole_problem == 0 .AND. COMMAND_ARGUMENT_COUNT() < 4 ) THEN
            ! If less than 4 command line arguments are given throw an error and stop
            WRITE(*,*) 'ERROR: WAS EXPECTING 4 COMAND LINE ARGUMENTS, ONLY GOT', COMMAND_ARGUMENT_COUNT()
            WRITE(*,*) 'IF YOU WANT TO USE THE DEFAULT VALUES, RUN THE PROGRAM WITHOUT ANY ARGUMENTS'
            WRITE(*,*) 'IF YOU WANT TO SUPPLY ARGUMENTS TO THE PROGRAM FIRST:'
            WRITE(*,*) 'DECIDE WHETHER TO SOLVE THE WHOLE PROBLEM (I=1) OR JUST THE VALUE FUNCTION (I=0)'
            WRITE(*,*) 'IF I = 1 THEN USE:    ./program $I $K $L    TO SUPPLY INITIAL GUESSES FOR K AND L'
            WRITE(*,*) 'IF I = 0 THEN USE:    ./program $I $r $w $b TO SUPPLY VALUES FOR r, w AND b'
            WRITE(*,*) 'EXITING...'
            STOP
        ELSE IF (solve_whole_problem == 0) THEN
            
            ! Read User supplied values for r, w, b
            CALL GET_COMMAND_ARGUMENT(2,r_char)
            READ(r_char,*)r
    
            CALL GET_COMMAND_ARGUMENT(3,w_char)
            READ(w_char,*)w
        
            CALL GET_COMMAND_ARGUMENT(4,b_char)
            READ(b_char,*)b
        ELSE IF (solve_whole_problem == 1 .AND. COMMAND_ARGUMENT_COUNT() < 3 ) THEN
            ! If less than 4 command line arguments are given throw an error and stop
            WRITE(*,*) 'ERROR: WAS EXPECTING 3 COMAND LINE ARGUMENTS, ONLY GOT', COMMAND_ARGUMENT_COUNT()
            WRITE(*,*) 'IF YOU WANT TO USE THE DEFAULT VALUES, RUN THE PROGRAM WITHOUT ANY ARGUMENTS'
            WRITE(*,*) 'IF YOU WANT TO SUPPLY ARGUMENTS TO THE PROGRAM FIRST:'
            WRITE(*,*) 'DECIDE WHETHER TO SOLVE THE WHOLE PROBLEM (I=1) OR JUST THE VALUE FUNCTION (I=0)'
            WRITE(*,*) 'IF I = 1 THEN USE:    ./program $I $K $L    TO SUPPLY INITIAL GUESSES FOR K AND L'
            WRITE(*,*) 'IF I = 0 THEN USE:    ./program $I $r $w $b TO SUPPLY VALUES FOR r, w AND b'
            WRITE(*,*) 'EXITING...'
            STOP
        ELSE
            ! Read User supplied values for K and L
            CALL GET_COMMAND_ARGUMENT(2,K_char)
            READ(K_char,*)K
        
            CALL GET_COMMAND_ARGUMENT(3,L_char)
            READ(L_char,*)L
    
            ! Using the values for K and L, solve for r, w
        END IF
    END IF
    
    
    DO WHILE (ERR > TOL .AND. COUNT_ITER  < MAX_ITER)
        ! Set up the grids and allocate space for the policy functions
        call housekeeping()                
        ! To get b we use the the pupulation distribution 
        IF (solve_whole_problem == 1) THEN
            do j = J_R, nJ
                m = m + mu(j)
            end do
            w = (1d0 - ALPHA) * ( ( K ** ALPHA) * ( L ** (-ALPHA) ) )
            r = ALPHA * ( ( K ** (ALPHA - 1)) * ( L **(1 - ALPHA) ) ) - DELTA
            b = THETA * w * L / m
        END IF
        
        call V_Func_Ret()                   ! Solve for the value function for retirees
        
        call V_Func_Work()                   ! Solve for the value function for working age agents
        
        call steady_state_dist()            ! Solve for the steady state distribution
        !  We calculate the steady state aggregate capital and labor
        ! WRITE(*,*) 'r = ', r, ' w = ', w, ' b = ', b, 'K = ', K, 'L = ', L, 'SUM(F) = ', SUM(F_SS)
        
        K_SS = 0d0
        L_SS = 0d0
        do j = 1, nJ
            do i_Z = 1, nZ
                do i_A = 1, nA
                    K_SS = K_SS + F_SS(i_A, i_Z, j)*grid_A(i_A)
                    L_SS = L_SS + F_SS(i_A, i_Z, j)*pf_l(i_A, i_Z, j)
                end do
            end do

            ! WRITE(*,*) "Age = ", j, "K_SS = ", K_SS, "L_SS = ", L_SS
        end do 

        IF (solve_whole_problem == 0) THEN
            COUNT_ITER = MAX_ITER + 100
        ELSE
            ! Calculate Error
            ERR = MAX(ABS(K_SS - K), ABS(L_SS - L))

            ! Tune adjustment parameter
            IF (err > tol*5 .AND. LAMBDA <= 0.85) THEN
                LAMBDA = 0.85
            ELSE IF (err > tol*1.01 .AND. LAMBDA <= 0.90)  THEN
                LAMBDA = 0.90
            ELSE IF (LAMBDA<= 0.95)  THEN
                LAMBDA = 0.95
            END IF
            ! Update Guesses and the values of r, w and b
            K = (1 - LAMBDA) * K_SS + LAMBDA * K
            L = (1 - LAMBDA) * L_SS + LAMBDA * L
            r = ALPHA * ( K ** (-ALPHA)) * ( L ** ALPHA )
            w = (1d0 - ALPHA) * ( K ** ALPHA) * ( L ** (-ALPHA) )
            b = THETA * w * L / m
        END IF
        IF (solve_whole_problem == 1) THEN
            WRITE(*,*)'Iteration:',COUNT_ITER,'ERR=',ERR,'K=', K ,'L=', L,'LAMBDA=',LAMBDA
        END IF
        ! Increment iteration counter
        COUNT_ITER = COUNT_ITER + 1
        ! Clean the Steady State values of K and L
    END DO
    ! End Computational Timer
    call system_clock(end)
    IF (solve_whole_problem == 1) THEN
        write(*,*) "Total elapsed time = ", real(end - beginning) / real(rate)," seconds"
    END IF
    call coda()                         ! Write resutls to file(s)
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

    DOUBLE PRECISION :: sum_mu=0.0d0

    if (COUNT_ITER < 1) then
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
        
        ! initialize the population distribution
        mu(1) = 1.0d0
        sum_mu = mu(1)
        do j = 2, nJ
            mu(j) = mu(j-1)/(1 + n)
        end do
        
        sum_mu = sum(mu, dim=1)
        do j = 1, nJ
            mu(j) = mu(j)/sum_mu
        end do
    end if

    ! Initialize Policy Functions
    do j = 1, nJ
        do i_A = 1,nA
            do i_Z = 1,nZ
                if ( j < nJ ) then
                    pf_c(i_A, i_Z, j) = 0d0
                    pf_A(i_A, i_Z, j) = 0d0
                    pf_A_ind(i_A, i_Z, j) = 0
                    pf_v(i_A, i_Z, j) = 0d0
                    pf_l(i_A, i_Z, j) = 0d0
                else
                    pf_c(i_A, i_Z, j) = grid_A(i_A)*(1 + r) + b
                    pf_A(i_A, i_Z, j) = 0d0
                    pf_A_ind(i_A, i_Z, j) = 1
                    pf_v(i_A, i_Z, j) = pf_c(i_A, i_Z, j)**(GAMMA * (1 - SIGMA)) / (1 - SIGMA)
                    pf_l(i_A, i_Z, j) = 0d0
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
    ! INTEGER :: count = 0
    do j = nJ-1, J_R, -1 ! Loop backwards over the age groups until the retirement age
        min_ind_search = 1
        do i_A = 1, nA   ! Loop over the asset grid
            cand_max = -100000000d0 ! A large negative number for candidate to maximum
            a_current = grid_A(i_A)
            do i_Anext = min_ind_search, nA ! Loop over all possible choises of the next asset level
                a_next = grid_A(i_Anext)
                c_current = (1 + r) * a_current + b - a_next
                if (c_current >= 0d0) then ! Check if consumption is positive
                    u_current = (c_current**( GAMMA * (1 - SIGMA) )) / (1 - SIGMA)
                    v_current = u_current + BETA * pf_v(i_Anext, 1, j+1)
                    if (v_current > cand_max) then  ! Update the candidate maximizer
                        cand_max = v_current
                        a_maxim = a_next
                        a_ind_maxim = i_Anext
                        pf_v(i_A, 1, j) = cand_max
                        pf_v(i_A, 2, j) = cand_max
                        pf_A(i_A, 1, j) = a_maxim
                        pf_A(i_A, 2, j) = a_maxim
                        pf_c(i_A, 1, j) = c_current
                        pf_c(i_A, 2, j) = c_current
                        pf_A_ind(i_A, 1, j) = a_ind_maxim
                        pf_A_ind(i_A, 2, j) = a_ind_maxim
                    end if
                else
                    exit ! If consumption is negative, then stop searching
                end if
                min_ind_search = a_ind_maxim
            end do ! End of loop over i_Anext
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
    use omp_lib
    implicit none
    iMaxThreads = omp_get_max_threads()
    call omp_set_num_threads(iMaxThreads) ! Set the number of threads to use
    
    do j = J_R-1, 1, -1 ! Loop backwards over the age groups from last working age to the birth age
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE( i_Z , i_A, i_Anext, a_current, a_next, Pr, l_opt, c_current, u_current, v_next, &
        !$OMP 																 v_current, cand_max, a_maxim , a_ind_maxim, l_opt_cand)
        !$OMP DO
        do i_Z = 1, nZ   ! Loop over the productivity grid 
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
            e = ez(j) * grid_z(i_Z) ! Compute the effective productivity for this age group             
            min_ind_search = 1 ! Exploiting monotonicity in the policy function
            do i_A = 1, nA   ! Loop over the asset grid
                a_current = grid_A(i_A) ! Compute the current asset holdings level
                cand_max = -500000d0 ! A large negative number for candidate to maximum
                l_opt_cand = -500000d0
                do i_Anext = min_ind_search, nA ! Loop over all possible choises of the next asset level
                    a_next = grid_A(i_Anext)
                    ! *********************
                    ! Compute agent choices
                    ! *********************
                    ! Calculate optimal labor supply
                    l_opt = ( GAMMA * (1-THETA) * e * w - (1-GAMMA)*((1+r) * a_current - a_next) ) / ( (1 - THETA) * w * e ) 
                    if (l_opt < 0d0) then ! Check if labor is negative
                        l_opt = 0d0
                    end if
                    if (l_opt > 1d0) then ! Check if labor is greater than one
                        l_opt = 1d0 ! If greater than one, set to one
                    end if 
                    ! Calculate consumption
                    c_current = w * (1-THETA) * e * l_opt + (1+r) * a_current - a_next
                    if ( c_current < 0d0 ) then ! Check if consumption is negative
                        exit ! If negative, stop searching
                    end if  
                    ! **************************
                    ! Compute the value function
                    ! **************************
                    u_current = ( ( (c_current**GAMMA ) * ( (1 - l_opt) ** (1 - GAMMA) ) ) ** (1 - SIGMA) ) / (1 - SIGMA) !  Compute the utility
                    v_next = pf_v(i_Anext, 1, j+1) * Pr(1) + pf_v(i_Anext, 2, j+1) * Pr(2) ! Compute the expected value function for the next period 
                    v_current = u_current + BETA * v_next ! Compute the value function
                    if (v_current > cand_max) then  ! Update the candidate maximizer
                        cand_max = v_current
                        a_maxim = a_next
                        a_ind_maxim = i_Anext
                        l_opt_cand = l_opt * e
                        c_cand = c_current
                    !     if (j == J_R-3 .and. i_z == 1 .and. i_A == 350) then
                    ! WRITE(*,*) "started form:", min_ind_search, "a_next = ", a_next, "l_opt = ", l_opt, "c_current = ", c_current
                    !         ! WRITE(*,*) "v_current = ", cand_max, "c_current = ", c_cand, "a_maxim = ", a_maxim
                    !     end if
                    end if
                end do ! End of loop over i_Anext 
                pf_v(i_A, i_Z, j) = cand_max
                pf_A(i_A, i_Z, j) = a_maxim
                pf_A_ind(i_A, i_Z, j) = a_ind_maxim
                pf_c(i_A, i_Z, j) = c_cand
                pf_l(i_A, i_Z, j) = l_opt_cand
                min_ind_search = a_ind_maxim
            end do ! End of loop over i_A
        end do ! End of loop over i_Z
        !$OMP END DO
        !$OMP END PARALLEL
        
    end do ! End of loop over j
    
    end subroutine V_Func_Work ! end of subroutine

! ************************************************************************
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : steady_state_dist
!
! description : Gives the steady state distribution of the model
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
subroutine steady_state_dist()
    use params_grid
    ! use omp_lib
    implicit none

    ! Start by cleaning up the distriburtion
    do j = 1, nJ
        do i_A = 1,nA
            do i_Z = 1,nZ
                F_SS(i_A, i_Z, j) = 0.0d0
            end do
        end do
    end do
    
    F_SS(1, 2, 1) = mu(1) * p_L
    F_SS(1, 1, 1) = mu(1) * p_H

    do j = 2, nJ
        do i_z = 1, nZ
            if (i_Z == 1) then
                Pr(1) = PI_HH
                Pr(2) = PI_HL
            else
                Pr(1) = PI_LH
                Pr(2) = PI_LL
            end if
            ! WRITE(*,*) "Pr(1) + Pr(2) = ", Pr(1) + Pr(2) 
            do i_A = 1, nA
                a_ind_maxim = pf_A_ind(i_A, i_z, j-1)
                if (a_ind_maxim == 0d0) then ! Level not reached
                    CYCLE
                end if
                F_SS(a_ind_maxim, 1, j) = F_SS(a_ind_maxim, 1, j) + F_SS(i_A, i_z, j-1) * Pr(1) * (mu(j)/mu(j-1))
                F_SS(a_ind_maxim, 2, j) = F_SS(a_ind_maxim, 2, j) + F_SS(i_A, i_z, j-1) * Pr(2) * (mu(j)/mu(j-1))
            end do ! End of loop over i_A
        end do ! End of loop over i_z
        
    end do ! End of loop over j
    end subroutine steady_state_dist ! end of subroutine
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
    INTEGER                         :: rc1, rc2 ! Identifier for error checking

    open(unit=2, file='./PS3/FortranCode/results.csv', status='replace', action='write', iostat=rc1)
    if (rc1 /= 0) stop 'Unable to open file results.csv'
    200 format(i2,2x, f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,i4,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x)
    
    do j = 1, nJ
        do i_Z = 1, nZ
            do i_A = 1, nA
write(2,200)j,grid_Z(i_Z),grid_A(i_A),pf_l(i_A, i_Z, j),pf_c(i_A, i_Z, j),pf_A(i_A, i_Z, j),pf_A_ind(i_A, i_Z, j),pf_v(i_A, i_Z, j)
            end do
        end do
    end do

    open(unit=3, file='./PS3/FortranCode/agg_results.csv', status='replace', action='write', iostat=rc2)
    if (rc2/=0) stop 'Unable to open file agg_results.csv'
    300 format(f25.15,2x,f25.15)
    write(3, 300) K_SS, L_SS

    

return
end subroutine coda ! end of subroutine