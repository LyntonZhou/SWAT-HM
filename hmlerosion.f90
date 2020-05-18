      subroutine hmlerosion(iwave)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates h.metal transported with suspended sediment 

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    enratio       |none         |enrichment ratio calculated for day in HRU
!!    hruhml(:)     |none         |h.metal use flag:
!!                                | 0: no h.metals used in HRU
!!                                | 1: h.metals used in HRU
!!    ihru          |none         |HRU number
!!    iwave         |none         |flag to differentiate calculation of HRU and
!!                                |subbasin sediment calculation
!!                                |iwave = 0 for HRU
!!                                |iwave = subbasin # for subbasin
!!    n_hml_mx      |none         |number of different metals used in the simulation
!!    n_hml_no(:)   |none         |array of unique metals used in watershed
!!    hml_enr(:,:)  |none         |h.metal enrichment ratio
!!    sedyld(:)     |metric tons  |daily soil loss caused by water erosion in
!!                                |HRU
!!    sol_hmlkp(:,:)|(kg/kg)/(kg/m3)|h. metal sorption coefficient, Kp; the
!!                                |ratio of the concentration in the solid
!!                                |phase to the concentration in solution
!!    sol_hml_lab(:,:,:)|kg/ha   |amount of h.metal in layer in HRU
!!    sub_hml(:,:)  |kg/ha        |amount of h.metal in layer in subbasin

!!    sol_bd(:,:)  |Mg/m**3       |bulk density of the soil
!!    sol_z(:,:)   |mm            |depth to bottom of soil layer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hml_sed(3,:,:)     |kg/ha   |h.metal loading from HRU sorbed onto
!!                                |sediment
!!    sol_hml_lab(:,:,:)|kg/ha   |amount of  labile h.metal in layer in HRU
!!    sol_hml_lab(:,:,:)|kg/ha   |amount of non-labile h.metal in layer in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    conc        |(kgCd/kgSolid)|concentration of h.metal in soil
!!    soil_weight_in_hur_ha      |Soil mass per ha, (kg/ha)
!!    er          |none          |enrichment ratio for h.metals
!!    j           |none          |HRU number
!!    k           |none          |counter
!!    kk          |none          |h.metal number from database
!!    xx          |kg/ha         |amount of h.metal in soil
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer, intent (in) :: iwave
      integer :: j, k, kk
      real :: xx1, xx2, conc, er

      j = 0
      j = ihru

      if (hruhml(j) == 0) return

      dg = sol_z(1,j)

      do k = 1, n_hml_mx
          kk = 0
          kk = n_hml_no(k)

          if (kk > 0) then
              xx1 = 0.
              xx2 = 0.
              if (iwave <= 0) then
                  xx1 = sol_hml_lab(k,j,1) ! 应根据侵蚀深度计算进入悬浮质的量
                  xx2 = sol_hml_nlab(k,j,1)
              else
                  xx1 = sub_hml(kk,iwave)
                  xx2 = sub_hml(kk,iwave)
              end if

              if (xx1 >= 1.e-8 .or. xx2 >= 1.e-8) then
                  er = 0.
                  if (hml_enr(k,j) > 0.) then
                      er = hml_enr(k,j)  ! 应根据侵蚀深度调整降低er
                      !er = enratio
                  else
                      er = enratio
                  end if

                  ! concentration of  labile h.metal in soil
                  !! HRU sediment calculations
             soil_weight_in_hru_ha = sol_bd(1,j) * dg * 1.e4	! Unit: kg Soil per ha;根据侵蚀深度计算
                  ! Convert unit: Mg/m3 * mm per ha=ton per ha
             conc=sol_hml_lab(k,j,1)/ soil_weight_in_hru_ha ! Unit: kg per kg
                  ! conc unit: Kg
             hml_sed(1,k,j) = 1000. * sedyld(j) * conc * er / hru_ha(j)
             conc=sol_hml_nlab(k,j,1)/ soil_weight_in_hru_ha
                  ! conc unit: Kg
             hml_sed(2,k,j) = 1000. * sedyld(j) * conc * er / hru_ha(j)

                  if (hml_sed(1,k,j) < 0.) hml_sed(1,k,j) = 0.
                  if (hml_sed(2,k,j) < 0.) hml_sed(2,k,j) = 0.
                  if (hml_sed(1,k,j) > xx1) hml_sed(1,k,j) = xx1
                  if (hml_sed(2,k,j) > xx2) hml_sed(2,k,j) = xx2

                  sol_hml_lab(k,j,1) = xx1 - hml_sed(1,k,j)  ! 当侵蚀深度不足sol_z(1)时可以，超过时怎么办？
                  !sol_hml_nlab(k,j,2) = xx2 - hml_sed(2,k,j)
                  sol_hml_nlab(k,j,1) = xx2 - hml_sed(2,k,j)

              end if
          end if
      end do

      return
      end
