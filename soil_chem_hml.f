      subroutine soil_chem_hml

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine initializes soil chemical properties associating heavy metals

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hruhml(:)     |none           |heavy metal involvement flag:
!!                                  | 0: no metals involved in HRU
!!                                  | 1: metals involved in HRU
!!    i             |none           |HRU number
!!    n_hml_mx      |none           |number of different metals used in the simulation
!!    n_hml_no(:)   |none           |array of unique metals used in watershed
!!    conv_wt(:,:)  |none           |factor which converts kg/kg soil to kg/ha
!!    sol_cov(:)    |kg/ha          |amount of residue on soil surface
!!    sol_hum(:,:)  |kg humus/ha    |amount of organic matter in the soil layer classified as humic substances
!!    sol_rsd(:,:)  |kg/ha          |amount of organic matter in the soil layer classified as residue
!!    sol_z(:,:)    |mm             |depth to bottom of soil layer
!!                                 
!!    sol_hmlkp(:)  |(mg/kg)/(kg/m3)|Soil partition coef. between solid and aqueous phase
!!    sol_bd(:,:)   |Mg/m**3        |bulk density of the soil
!!    sol_cbn(:,:)  |%              |percent organic carbon in soil layer
!!    sol_nly(:)    |none           |number of soil layers 
!!    sol_hml_lab(:,:,:)   |mg/kg   |amount of metal stored as  labile species in solid
!!    sol_hml_nlab(:,:,:)  |mg/kg   |amount of metal stored as non-labile species in solid

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    sol_hml_sol(:,:,:)   |kg /ha        |amount of metal stored as free ion in solution
!!    sol_hml_lig(:,:,:)   |kg /ha        |amount of metal stored as ligand-fixed species in solution
!!    sol_hml_lab(:,:,:)   |kg /ha        |amount of metal stored as  labile species in solid
!!    sol_hml_nlab(:,:,:)  |kg /ha        |amount of metal stored as non-labile species in solid


!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dg          |mm            |depth of layer
!!    j           |none          |counter
!!    jj          |none          |dummy variable to hold value
!!    n           |none          |counter
!!    nly         |none          |number of soil layers
!!    soldepth    |mm            |depth from bottom of 1st soil layer to
!!                               |the bottom of the layer of interest
!!    solpst      |mg/kg         |concentration of pesticide in soil
!!    summinp     |kg P/ha       |amount of phosphorus stored in the mineral P
!!                               |pool in the profile
!!    sumno3      |kg N/ha       |amount of nitrogen stored in the nitrate pool
!!                               |in the soil profile
!!    sumorgn     |kg N/ha       |amount of nitrogen stored in the organic N
!!                               |pools in the profile
!!    sumorgp     |kg P/ha       |amount of phosphorus stored in the organic P
!!                               |pools in the profile
!!    wt1         |none          |converts mg/kg (ppm) to kg/ha
!!    xx          |none          |variable to hold value
!!    zdst        |none          |variable to hold value
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: nly, j, jj, n, i_hml
      real :: xx, dg, wt1, zdst, soldepth

      nly = 0
      nly = sol_nly(i)

!!    calculate initial metal contents of layers,converting mg/kg (ppm) to kg/ha
      xx = 0.
      do i_hml = 1, mhml
          do j = 1, nly
              if (sol_hml_lab(i_hml, i, j) <= 0.) then
                  zdst = 0.
                  zdst = Exp(-sol_z(j,i) / 1000.)
                  !zdst = Exp(-sol_z(j,i) / 333.)
                  !zdst = Exp(-sol_z(j,i) / 300.)
                  !zdst = Exp(-sol_z(j,i) / 250.)
                  !zdst = Exp(-sol_z(j,i) / 200.)
                  sol_hml_lab(i_hml, i, j) 
     &            = zdst * sol_hml_lab(i_hml, i, 1)
                  sol_hml_nlab(i_hml, i, j) 
     &            = zdst * sol_hml_nlab(i_hml, i, 1)
              end if
          end do
      end do

      ! converting mg/kg (ppm) to kg/ha
      dg = 0.
      do i_hml=1, mhml
          do j = 1, nly
              if(j == 1) then
                  dg = sol_z(j,i)
              else
                  dg = sol_z(j,i) - sol_z(j-1,i)
              end if
              wt1 = 0.
              wt1 = sol_bd(j,i) * dg / 100.

              sol_hml_lab(i_hml,i,j) = sol_hml_lab (i_hml,i,j) * wt1
              sol_hml_nlab(i_hml,i,j) = sol_hml_nlab (i_hml,i,j) * wt1
              sumhml_exch = sumhml_exch + sol_hml_lab(i_hml,i,j)
              sumhml_nlab = sumhml_nlab + sol_hml_nlab(i_hml,i,j)
          end do
      end do

      return
      end
