      subroutine hmlinput
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine adds heavy metal from atmosphere to the soil profile
!!    and adds heavy metal from agriculture to the soil profile
      
!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    curyr        |none          |current year of simulation
!!    hml_dep_total |g/ha/yr      |atmospheric deposition of total metal
!!    hml_dep_frac  |             |fraction of atmospheric deposition
!!    nyskip       |none          |number of years to skip output summarization
!!                                |and printing
!!    precipday    |mm H2O        |precipitation for the day in HRU
!!    atmo_hml     |kg/ha         |average annual amount of metal added to soil
!!                                |by deposition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: j

      j = 0
      j = ihru
     
      !! calculate metal deposition 
      sol_hml_lab(1,j,1) = sol_hml_lab(1,j,1) + hml_dep_total(j) * hml_dep_frac(j)/365.
      
      sol_hml_nlab(1,j,1) = sol_hml_nlab(1,j,1) + hml_dep_total(j) * (1-hml_dep_frac(j))/365.
     
      !! metal from fertilizers 
      sol_hml_lab(1,j,1) = sol_hml_lab(1,j,1) + hml_agr_total(j)* hml_agr_frac(j)/365.
      
      sol_hml_nlab(1,j,1) = sol_hml_nlab(1,j,1) + hml_agr_total(j)* (1-hml_agr_frac(j))/365.
      
      !! summary calculations
      if (curyr > nyskip) then
        atmo_hml = atmo_hml + hml_dep_total(j) * hru_dafr(j)
      end if

      return
      end