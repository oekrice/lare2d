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
  
MODULE zremap

  USE shared_data; USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: remap_z

CONTAINS

  SUBROUTINE remap_z

    REAL(num) :: v_advect, flux1, flux2

    DO iy = -1, ny + 1
      DO ix = -1, nx + 1
        ixm = ix - 1
        iym = iy - 1

        v_advect = (vz1(ix ,iy) + vz1(ix ,iym)) * 0.5_num
        flux1 = v_advect * bx(ix ,iy) * dt

        v_advect = (vz1(ixm,iy) + vz1(ixm,iym)) * 0.5_num
        flux2 = v_advect * bx(ixm,iy) * dt

        bz(ix,iy) = bz(ix,iy) + flux1 - flux2
      END DO
    END DO

    DO iy = -1, ny + 1
      DO ix = -1, nx + 1
        ixm = ix - 1
        iym = iy - 1

        v_advect = (vz1(ix,iy ) + vz1(ixm,iy )) * 0.5_num
        flux1 = v_advect * by(ix,iy ) * dt

        v_advect = (vz1(ix,iym) + vz1(ixm,iym)) * 0.5_num
        flux2 = v_advect * by(ix,iym) * dt

        bz(ix,iy) = bz(ix,iy) + flux1 - flux2
      END DO
    END DO

    CALL bz_bcs

  END SUBROUTINE remap_z

END MODULE zremap
