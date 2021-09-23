! ************************************************************************
! Filename :T_operator_paralle.f90
!
! Author : Mitchell Valdes-Bobes (adapted from a code from Philip Coyle)
!
! Date Created : September 17th, 2021
!
! Description : This program will use dynamic programming techniques to solve
! a simple neoclassical growth model with a two state markov productivity shock.
!
! Routine:
! cd "to/directory/where/code/is/located"
! gfortran -fopenmp -O2 -o name_of_compiled_file T_operator_parallel.f90
! ifort -qopenmp -O2 -o name_of_compiled_file neogrowth_parallel.f90
! ./name_of_compiled_file
! ************************************************************************

! ************************************************************************
! ------------------------------------------------------------------------

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

	! Price as a parameter 
	! TODO: make this a variable
	DOUBLE PRECISION			:: 				q ! Price 
	INTEGER 					::				n_glob_iter = 0
	double precision, parameter :: 				cBET 				= 0.9932d0
	double precision, parameter :: 				cALPHA 				= 1.5d0

	! Model Probabilities
	double precision, parameter :: 				Pgg 				= 0.97d0
	double precision, parameter :: 				Pbg 				= 1d0 - Pgg
	double precision, parameter :: 				Pbb 				= 0.5d0
	double precision, parameter :: 				Pgb 				= 1 - Pbb
	double precision	 		:: 				Pr(2)


	! Tolerance level for convergence and max itations
	double precision, parameter :: 				tol			 		= 1d-5!
	integer, parameter 			:: 				max_it 				= 10000
	integer 					::				it		 			= 1 !itation counter
	integer 					:: 				converged			= 0


	! -----------------------------------------------------------------------
	! ****************************GRID SET**********************************
	! -----------------------------------------------------------------------
	! Set up for discritizing the state space (Capital Grid)
	integer						:: 				i_A, i_Apr
	integer, parameter 			:: 				n_A 			= 1000
	double precision 			:: 				grid_A(n_A)
	double precision, parameter :: 				min_A 			= -2d0 !1d0
	double precision, parameter :: 				max_A 			= 5d0 !75d0
	double precision, parameter :: 				step_A 			= (max_A - min_A)/(dble(n_A) - 1d0)
	double precision			:: 				A_today
	double precision			:: 				A_tomorrow

	! Set up for discritizing the state space (Employment Grid)
	integer						:: 				i_S
	double precision, parameter :: 				Sg 					= 1d0
	double precision, parameter :: 				Sb 					= 0.5d0
	integer, parameter 			:: 				n_S 				= 2
	double precision 			:: 				grid_S(n_S)
	double precision			:: 				S_today

	integer 					:: 				i_state
	integer, parameter  		:: 				n_state 		= n_A*n_S ! Writing grid in one whole loop

	! Global variables for Dynamic Progamming
	double precision 			:: 				c_today
	double precision 			:: 				c_today_temp
	double precision 			:: 				y_today
	double precision 			:: 				A_tomorrow_max
	double precision 			:: 				v_today
	double precision 			:: 				v_today_temp
	double precision 			:: 				v_tomorrow
	double precision 			:: 				v_yesterday


	! Allocating space for Policy Functions
	double precision 			:: 				pf_c(n_A, n_S)
	double precision 			:: 				pf_A(n_A, n_S)
	double precision 			:: 				pf_v(n_A, n_S)

	integer 					::			  i_stat
	integer 					::			  iMaxThreads

end module params_grid


! ************************************************************************
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! program : neogrowth
!
! Description : This program will use dynamic programming techniques to solve
! a simple neoclassical growth model with a two state markov productivity shock
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

program neogrowth

use params_grid
use omp_lib


implicit none





! Begin Computational Timer
integer 						::				beginning, end, rate
CHARACTER(100) 					:: 				q_char, iter_char
DOUBLE PRECISION 				:: 				q_double

call system_clock(beginning, rate)

IF(COMMAND_ARGUMENT_COUNT() < 1)THEN
	WRITE(*,*)'ERROR, PRICE SHOULD BE PASSED AS COMMAND-LINE ARGUMENT, STOPPING'
	STOP
ENDIF

CALL GET_COMMAND_ARGUMENT(1,q_char)
READ(q_char,*)q

IF(COMMAND_ARGUMENT_COUNT() > 1)THEN
	CALL GET_COMMAND_ARGUMENT(2,iter_char)
	READ(iter_char,*)n_glob_iter
ENDIF



! write(*,*) ""
! write(*,*) "******************************************************"
! write(*,*) "******************Start OF T operator*****************"
! write(*,*) "******************************************************"
! write(*,*) ""

! Initialize Grids and Policy Functions
call housekeeping()

! Do Value Function Iteration
call bellman()

call system_clock(end)
! write(*,*) ""
! write(*,*) "******************************************************"
! write(*,*) "Total Iterations   = ", it
! write(*,*) "Total elapsed time = ", real(end - beginning) / real(rate)," seconds"
! write(*,*) "******************************************************"
! write(*,*) ""

! Write results
call coda()

! write(*,*) ""
! write(*,*) "******************************************************"
! write(*,*) "*******************END OF T operator******************"
! write(*,*) "******************************************************"
! write(*,*) ""

end program neogrowth

! ************************************************************************
! ************************************************************************
! ************************************************************************
! **************************** SUBROUTINES *******************************
! ************************************************************************
! ************************************************************************
! ************************************************************************


! ************************************************************************


! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : housekeeping
!
! description : Initializes Grids and Policy Functions
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine housekeeping()

use params_grid
use omp_lib

implicit none

! Discretizing the state space (capital)
do i_A = 1,n_A
	grid_A(i_A) = min_A + (dble(i_A) - 1d0)*step_A
end do

! Discretizing the state space (productivity)
do i_S = 1,n_S
	if (i_S == 1) then
		grid_S(i_S) = Sg
	else
		grid_S(i_S) = Sb
	end if
end do

! Setting up Policy Function guesses
if (n_glob_iter == 0) then
	! If this is the first time we run the program, we set the policy functions to 0
	do i_A = 1,n_A
		do i_S = 1,n_S
			pf_c(i_A, i_S) 			= 0d0
			pf_A(i_A, i_S) 		    = 0d0
			pf_v(i_A, i_S)    	    = 0d0
		end do
	end do
else 
	! If not the first time we run the program we use the previous policy functions
	call read_data()
ENDIF


return

end subroutine housekeeping




! ************************************************************************


! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : bellman
!
! description : Solves the dynamic programming problem for policy and value
! functions.
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine bellman()

use params_grid
use omp_lib

implicit none
! allocating space for policy function updates
double precision 						:: 				pf_c_up(n_A, n_S)
double precision 						:: 				pf_A_up(n_A, n_S)
double precision 						:: 				pf_v_up(n_A, n_S)

double precision 						:: 				diff_c
double precision 						:: 				diff_A
double precision 						:: 				diff_v
double precision 						:: 				max_diff


converged = 0
it = 1

! Begin Dynamic Programming Algo
do while (converged == 0 .and. it < max_it)

	iMaxThreads = omp_get_max_threads()
	call omp_set_num_threads(iMaxThreads) ! iMaxThreads = 2*(num of cores)


	! When parallelizing over one grid, need to include state vars in private.
	! Break OMP line up for gfortran compiler
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i_state, i_S, i_A, S_today, Pr, A_today, A_tomorrow, y_today, c_today_temp, &
	!$OMP 																 v_tomorrow, v_today_temp, v_today, A_tomorrow_max, c_today)
	!$OMP DO
	do i_state = 1,n_state

		! ***********************************
		! Define the today's state variables
		! ***********************************
		i_S         = MOD(i_state,n_S) + 1
		i_A         = MOD(CEILING(dble(i_state)/dble(n_S)),n_A) + 1

		S_today = grid_S(i_S)
		A_today = grid_A(i_A)

		! *********************
		! Define today's shock
		! *********************
		if (i_S == 1) then
			Pr(1) = Pgg
			Pr(2) = Pbg
		else
			Pr(1) = Pgb
			Pr(2) = Pbb
		end if

		! ******************************************************
		! Solve for the optimal consumption / capital investment
		! ******************************************************
		v_today = -10000000000d0
		v_yesterday = v_today

		do i_Apr = 1,n_A
			A_tomorrow = grid_A(i_Apr)
			if ( A_tomorrow <= (S_today + A_today)/q ) then

			c_today_temp = S_today + A_today - q * A_tomorrow
			c_today_temp = max(0d0,c_today_temp)

			v_tomorrow = Pr(1)*pf_v(i_Apr,1) + Pr(2)*pf_v(i_Apr,2)

			v_today_temp = ( c_today_temp**(1d0 - cALPHA) - 1d0 )  / (1d0 - cALPHA ) + cBET*v_tomorrow

			if ( v_today_temp < v_yesterday ) then
				exit
			end if 

			if (v_today_temp > v_today) then
				v_today = v_today_temp
				A_tomorrow_max = A_tomorrow
			end if
			else 
				! write(*,*) "a' too high "
			end if

			v_yesterday = v_today
		end do

		! v_today = v_today_temp
		A_tomorrow = A_tomorrow_max
		c_today = S_today + A_today - q * A_tomorrow

		! *******************************
		! ****Update Policy Functions****
		! *******************************
		pf_c_up(i_A, i_S) = c_today
		pf_A_up(i_A, i_S) = A_tomorrow
		pf_v_up(i_A, i_S) = v_today
	end do
	!$OMP END DO
	!$OMP END PARALLEL


	! Find the difference between the policy functions and updates
	! diff_c  = maxval(abs(pf_c - pf_c_up))
	! diff_A  = maxval(abs(pf_A - pf_A_up))
	max_diff  = maxval(abs(pf_v - pf_v_up))

	! max_diff = diff_c + diff_A + diff_v

	! if (mod(it,50) == 0) then
	! 	write(*,*) ""
	! 	write(*,*) "********************************************"
	! 	write(*,*) "At itation = ", it
	! 	write(*,*) "Max Difference = ", max_diff
	! 	write(*,*) "********************************************"
	! end if


	if (it > max_it) then
		converged = 1
	end if

	if (max_diff < tol) then
		converged = 1
		! write(*,*) ""
		! write(*,*) "********************************************"
		! write(*,*) "At iteration = ", it
		! write(*,*) "Max Difference = ", max_diff
		! write(*,*) "********************************************"
	end if

	it = it+1

	pf_c 		= pf_c_up
	pf_A 		= pf_A_up
	pf_v		= pf_v_up
end do

return

end subroutine bellman

! ************************************************************************
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : coda
!
! description : Writes results to .csv file
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine read_data()
	use params_grid

	IMPLICIT NONE
    
    type :: data
        DOUBLE PRECISION  :: S
        DOUBLE PRECISION  :: A
        DOUBLE PRECISION  :: C
        DOUBLE PRECISION  :: A_next
        DOUBLE PRECISION  :: Value
    end type data

    type(data) :: data_prev(2000)
    integer            :: fu, rc, i

    open (action='read', file='value_funct.csv', iostat=rc, newunit=fu)

    if (rc /= 0) stop

    do i = 1, 2000
        read (fu, *, iostat=rc) data_prev(i)
        if (rc /= 0) exit
    end do

    close (fu)

	do i_A = 1,n_A
		do i_S = 1,n_S
			pf_c(i_A, i_S) = data_prev( 1000 * (i_S - 1) +  i_A)%C
			pf_A(i_A, i_S) = data_prev( 1000 * (i_S - 1) +  i_A)%A_next
			pf_v(i_A, i_S) = data_prev( 1000 * (i_S - 1) +  i_A)%Value
		end do
	end do

end

subroutine coda()

use params_grid
use omp_lib

implicit none

! write(*,*) ""
! write (*,*) "Writing PFs to CSV file"
open(unit = 2, file = 'value_funct.csv', status = 'replace', action = 'write', iostat = i_stat)
200 format(f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x)

do i_S = 1, n_S
	do i_A = 1,n_A
	    write(2,200) grid_S(i_S),grid_A(i_A), pf_c(i_A, i_S), pf_A(i_A, i_S), pf_v(i_A, i_S)
	end do
end do

return

end subroutine coda
