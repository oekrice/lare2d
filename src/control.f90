  ! Copyright 2020 University of Warwick

  ! Licensed under the Apache License, Version 2.0 (the "License");
  ! you may not use this file except in compliance with the License.
  ! You may obtain a copy of the License at

  !    http://www.apache.org/licenses/LICENSE-2.0

  ! Unless required by applicable law or agreed to in writing, software
  ! distributed under the License is distributed on an "AS IS" BASIS,
  ! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ! See the License for the specific language governing permissions and
  ! limitations under the License.
  
!******************************************************************************
! Set up the control values required by the core code
!******************************************************************************

MODULE control

  USE shared_data
  USE normalise

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: user_normalisation, control_variables, set_output_dumps

CONTAINS

  !****************************************************************************
  ! Normalisation constants
  !****************************************************************************

  SUBROUTINE user_normalisation

    ! Set the normalising constants for LARE
    ! This is needed to allow the use of some physics modules which are coded
    ! in SI units

    ! Gamma is the ratio of specific heat capacities
    gamma = 5.0_num / 3.0_num

    ! Average mass of an ion in proton masses
    ! The code assumes a single ion species with this mass
    mf = 1.2_num

    ! The equations describing the normalisation in LARE have three free
    ! parameters which must be specified by the end user. These must be the
    ! normalisation used for your initial conditions. Strictly only needed for
    ! non-ideal MHD terms.

    ! Magnetic field normalisation in Tesla
    B_norm = 1.0_num!0.03_num

    ! Length normalisation in m
    L_norm = 1.0_num!180.e3_num

    ! Density normalisation in kg / m^3
    rho_norm = 1.0_num!1.67e-4_num

  END SUBROUTINE user_normalisation



  !****************************************************************************
  ! General control variables. Commented in detail
  !****************************************************************************

  SUBROUTINE control_variables

    REAL(num), DIMENSION(11):: vars

    OPEN(1, FILE = 'variables.txt')
    READ(1, *) vars
    CLOSE(1)

    run_num = vars(1)
    nplots = vars(4)
    energy_variable = vars(5)
    outflow_variable = vars(6)
    bfact_variable = vars(7)
    density_variable = vars(8)
    shearfact_variable = vars(9)
    eta0_variable = vars(10)
    grav_variable = vars(11)

    print*, 'IMPORTED VARIABLES'
    print*, ''
    print*, 'nx', int(vars(2))
    print*, 'Run number', int(run_num)
    print*, 'Number of plots', int(nplots)
    print*, 'Initial Energy', energy_variable
    print*, 'Initial Outflow', outflow_variable
    print*, 'Magnetic Field Factor', bfact_variable
    print*, 'Initial Density', density_variable
    print*, 'Shearing Factor', shearfact_variable
    print*, 'Supergranular Diffusion Factor', eta0_variable
    print*, 'Gravity Factor', grav_variable
    print*, 'Max. time', vars(3), 'units, or', vars(3)*24.2*shearfact_variable, 'days'
    ! Set the number of gridpoints in x and y directions
    nx_global = int(vars(2))
    ny_global = int(vars(2)/2)

    ! Set the maximum number of iterations of the core solver before the code
    ! terminates. If nsteps < 0 then the code will run until t = t_end
    nsteps = -1
    ! The maximum runtime of the code
    t_end = vars(3)

    ! Shock viscosities as detailed in manual - they are dimensionless
    visc1 = 0.1_num
    visc2 = 1.0_num
    ! \nabla^2 v damping 
    ! visc3 is an array set initial conditions
    use_viscous_damping = .TRUE.

    ! Set these constants to manually override the domain decomposition.
    ! If either constant is set to zero then the code will try to automatically
    ! decompose in this direction
    nprocx = 0
    nprocy = 1

    ! The length of the domain in the x direction
    x_min = -1.0_num
    x_max = 1.0_num
    ! Should the x grid be stretched or uniform
    x_stretch = .FALSE.

    ! The length of the domain in the y direction
    y_min = -0.0_num
    y_max = 1.0_num
    ! Should the y grid be stretched or uniform
    y_stretch = .FALSE.

    ! Turn on or off the resistive parts of the MHD equations
    resistive_mhd = .TRUE.

    ! The background resistivity expressed as the inverse Lundquist number
    !eta_background = 0.0005_num   !This value produces nice arcade eruptions
    eta_background = 0.0005
    ! The critical current for triggering anomalous resistivity
    ! and the resistivity when above the critical current.
    ! The resistivity is expressed as the inverse Lundquist number.
    j_max = 1e15_num!100.0_num
    eta0 = 0.01_num

    ! Turn on or off the hall_mhd term in the MHD equations
    ! If true than lambda_i must be set in the initial conditions
    hall_mhd = .FALSE.

    ! Turn on or off the Braginskii thermal conduction term in
    ! the MHD equations
    ! WARNING: this is not robust. It is known to have problems
    ! with steep temperature gradients and very hot regions with
    ! large thermal conductivity. For many problems it is however
    ! fine.
    conduction = .TRUE.
    ! Apply a flux limiter to stop heat flows exceeding free streaming limit
    heat_flux_limiter = .TRUE.
    ! Fraction of free streaming heat flux used if limiter on
    flux_limiter = 0.06_num

    ! Use radiation as specified in SUBROUTINE rad_losses
    ! in src/radiative.f90
    radiation = .FALSE.

    ! Include user specified heating function as specified in 
    ! SUBROUTINE rad_losses user_defined_heating in src/radiative.f90
    coronal_heating = .FALSE.

    ! Remap kinetic energy correction. LARE does not perfectly conserve kinetic
    ! energy during the remap step. This missing energy can be added back into
    ! the simulation as a uniform heating. Setting rke to true turns on this
    ! addition.
    rke = .FALSE.

    ! The code to choose the initial conditions. The valid choices are
    ! IC_NEW     - Use set_initial_conditions in "initial_conditions.f90" to
    !              setup new initial conditions
    ! IC_RESTART - Load the output file with index restart_snapshot and use it
    !              as the initial conditions
    initial = IC_NEW
    restart_snapshot = 1

    ! If cowling_resistivity is true then the code calculates and
    ! applies the Cowling Resistivity to the MHD equations
    ! only possible if not EOS_IDEAL
    ! resistive_mhd must be TRUE for this to actaully be applied
    cowling_resistivity = .FALSE.

    ! Set the boundary conditions on the four edges of the simulation domain
    ! Valid constants are
    ! BC_PERIODIC - Periodic boundary conditions
    ! BC_OPEN     - Riemann far-field characteristic boundary conditions
    ! BC_USER     - User boundary conditions specified in boundary.f90
    xbc_min = BC_USER
    xbc_max = BC_USER
    ybc_min = BC_USER
    ybc_max = BC_USER

    !If any user boundaries are driven set this flag
    driven_boundary = .FALSE.

    ! Control Boris scheme for limiting the Alfven speed
    ! Logical boris to turn on/off
    ! va_max controls the effective mass density and is
    ! the reduced light speed in Boris's paper in Lare normalised units
    boris = .FALSE.
    va_max = 4.7e3_num

    ! Set the equation of state. Valid choices are
    ! EOS_IDEAL - Simple ideal gas for perfectly ionised plasma
    ! EOS_PI    - Simple ideal gas for partially ionised plasma
    ! EOS_ION   - EOS_PI plus the ionisation potential
    ! N.B. read the manual for notes on these choices
    eos_number = EOS_IDEAL
    ! EOS_IDEAL also requires that you specific whether
    ! the gas is ionised or not. Some stratified atmospheres
    ! only work for neutral hydrogen even though using MHD
    ! For fully ionised gas set .FALSE.
    ! For neutral hydrogen set .TRUE.
    ! This flag is ignored for all other EOS choices.
    neutral_gas = .TRUE.

    !An exponential moving average 
    !(https://en.wikipedia.org/wiki/Moving_average#Exponential_moving_average)
    !Tweak this to get a "good" cooling function that doesn't just remove all
    !heating effects
    ! Works for viscosity and first order resistive effects
    cooling_term = .FALSE.
    alpha_av = 0.05_num

  END SUBROUTINE control_variables



  !****************************************************************************
  ! Output controls.
  !****************************************************************************

  SUBROUTINE set_output_dumps

    ! The output directory for the code
    if (run_num < 10) then
      write (data_dir, "(A5, I1)") 'Data/', int(run_num)
    else if (run_num < 100) then
      write (data_dir, "(A5, I2)") 'Data/', int(run_num)
    else
      write (data_dir, "(A5, I3)") 'Data/', int(run_num)
    end if

    ! The interval between output snapshots.
    dt_snapshots = t_end / nplots

    ! Force dt to adjust to output exactly at times set by dt_snapshot
    force_exact_time_outputs = .FALSE.

    ! dump_mask is an array which specifies which quantities the code should
    ! output to disk in a data dump.
    ! The codes are
    ! 1  - rho
    ! 2  - energy
    ! 3  - vx
    ! 4  - vy
    ! 5  - vz
    ! 6  - bx
    ! 7  - by
    ! 8  - bz
    ! 9  - temperature
    ! 10 - pressure
    ! 11 - cs (sound speed)
    ! 12 - parallel_current
    ! 13 - perp_current
    ! 14 - neutral_faction
    ! 15 - eta_perp
    ! 16 - eta
    ! 17 - jx
    ! 18 - jy
    ! 19 - jz
    ! 20 - accumulated viscous and resistive heating
    ! If a given element of dump_mask is true then that field is dumped
    ! If the element is false then the field isn't dumped
    ! N.B. if dump_mask(1:8) not true then the restart will not work
    dump_mask = .FALSE.
    dump_mask(1:8) = .TRUE.
    dump_mask(19) = .TRUE.
    IF (eos_number /= EOS_IDEAL) dump_mask(14) = .TRUE.
    IF (cowling_resistivity) dump_mask(15) = .TRUE.
    IF (resistive_mhd) dump_mask(16) = .TRUE.

  END SUBROUTINE set_output_dumps

END MODULE control
