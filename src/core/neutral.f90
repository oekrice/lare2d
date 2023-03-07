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
  
! All the subroutines in this module are for the partially ionised flux
! emergence simulations; see Leake & Arber, 2006

MODULE neutral

  USE shared_data
  USE boundary

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: perpendicular_resistivity, neutral_fraction, setup_neutral
  PUBLIC :: get_neutral

CONTAINS

  SUBROUTINE setup_neutral

    ! Ion neutral collision cross section(m^2)

    REAL(num) :: sigma_in = 5.0e-19_num
    REAL(num) :: mbar

    ALLOCATE(xi_n(-1:nx+2, -1:ny+2))
    xi_n = 0.0_num

    IF (cowling_resistivity) THEN
      ALLOCATE(eta_perp(-1:nx+2, -1:ny+2))
      ALLOCATE(parallel_current(0:nx, 0:ny))
      ALLOCATE(perp_current(0:nx, 0:ny))
    END IF

    ! Temperature of the photospheric radiation field
    tr = 6420.0_num

    ! Calculate fbar^(2 / 3) in (k^-1 m^-2)
    f_bar = 2.0_num * pi * (me_si / hp_si) * (kb_si / hp_si)
    f_bar = f_bar**1.5_num

    ! Calculate tbar in (K)
    t_bar = ionise_pot_si / kb_si

    ! Calculate rbar in (kg^-1)
    mbar = mh_si * mf
    r_bar = 4.0_num / mbar

    ! Calculate eta_bar in (m^4 / (k s kg^2))
    eta_bar = 2.0_num * mbar / (SQRT(16.0_num * kb_si / (pi * mbar)) * sigma_in)

  END SUBROUTINE setup_neutral



  SUBROUTINE perpendicular_resistivity

    ! This subroutine calculates the cross field resistivity at the current
    ! temperature.

    REAL(num) :: f, xi_v, bxv, byv, bzv, bfieldsq, rho_v, t_v
    INTEGER :: ixp, iyp, ix2, iy2

    DO iy = 0, ny
      iyp = iy + 1
      DO ix = 0, nx
        ixp = ix + 1

        ! Get the vertex density
        rho_v = rho(ix ,iy ) * cv(ix ,iy ) + rho(ixp,iy ) * cv(ixp,iy ) &
            +   rho(ix ,iyp) * cv(ix ,iyp) + rho(ixp,iyp) * cv(ixp,iyp)
        rho_v = rho_v &
            / (cv(ix,iy ) + cv(ixp,iy ) + cv(ix,iyp) + cv(ixp,iyp))

        ! Get the vertex magnetic field
        bxv = bx(ix,iy) + bx(ix ,iyp)
        byv = by(ix,iy) + by(ixp,iy )
        bzv = bz(ix,iy) + bz(ixp,iy ) + bz(ix,iyp) + bz(ixp,iyp)
        bxv = 0.5_num  * bxv
        byv = 0.5_num  * byv
        bzv = 0.25_num * bzv
        bfieldsq = bxv**2 + byv**2 + bzv**2

        ! Get the vertex temperature
        t_v = 0.0_num
        DO iy2 = iy, iy + 1
          DO ix2 = ix, ix + 1
            t_v = t_v + (energy(ix2,iy2) &
                - (1.0_num - xi_n(ix2,iy2)) * ionise_pot) &
                / (2.0_num - xi_n(ix2,iy2))
          END DO
        END DO

        t_v = t_v * (gamma - 1.0_num) / 4.0_num

        xi_v = get_neutral(t_v, rho_v)

        f = MAX(1.0_num - xi_v, none_zero)
        IF (f > 0) THEN
          eta_perp(ix,iy) = eta_bar * xi_v / f * bfieldsq &
              / rho_v**2 / SQRT(t_v)
        ELSE
          eta_perp(ix,iy) = 0.0_num
        END IF
      END DO
    END DO

    eta_perp = MIN(eta_perp, 100.0_num)

  END SUBROUTINE perpendicular_resistivity



  FUNCTION get_neutral(t_v, rho_v)

    REAL(num), INTENT(IN) :: t_v, rho_v
    REAL(num) :: get_neutral
    REAL(num) :: bof, r, t_rad, dilution

    t_rad = tr
    dilution = 0.5_num

    bof = 1.0_num / (dilution * f_bar * t_rad * SQRT(t_v)) &
        * EXP((0.25_num * (t_v / t_rad - 1.0_num) + 1.0_num) * T_bar / t_v)
    r = 0.5_num * (-1.0_num + SQRT(1.0_num + r_bar * rho_v * bof))
    get_neutral = r / (1.0_num + r)

  END FUNCTION  get_neutral



  SUBROUTINE neutral_fraction

    REAL(num) :: T, rho0, e0, dx, x, fa, xi_a
    REAL(num), DIMENSION(2) :: ta
    INTEGER :: loop

    DO iy = -1, ny + 2
      DO ix = -1, nx + 2
        rho0 = rho(ix,iy)
        e0 = energy(ix,iy)
        ta = (gamma - 1.0_num) &
            * (/ MAX((e0 - ionise_pot) / 2.0_num, none_zero), e0 /)

        IF (ta(1) > ta(2)) THEN
          PRINT*, 'Temperature bounds problem', ta
          STOP
        END IF

        dx = ta(2) - ta(1)
        t = ta(1)

        DO loop = 1, 100
          dx = dx / 2.0_num
          x = t  + dx
          xi_a = get_neutral(x, rho0)
          fa = x - (gamma - 1.0_num) &
              * (e0 - (1.0_num - xi_a) * ionise_pot) / (2.0_num - xi_a)
          IF (fa <= 0.0_num) t = x
          IF (ABS(dx) < 1.e-8_num .OR. ABS(fa) < 1.e-8_num) EXIT
        END DO

        xi_n(ix,iy) = get_neutral(x, rho0)
      END DO
    END DO

  END SUBROUTINE neutral_fraction

END MODULE neutral
