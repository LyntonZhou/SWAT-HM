      subroutine hmlpu
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates plant heavy metal uptake

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units          |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    curyr       |none           |current year of simulation
!!    bioday      |kg             |biomass generated on current day in HRU
!!    hru_dafr(:) |km**2/km**2    |fraction of watershed area in HRU
!!    icr(:)      |none           |sequence number of crop grown within the
!!                                |current year
!!    idplt(:)    |none           |land cover code from crop.dat
!!    ihru        |none           |HRU number
!!    nro(:)      |none           |sequence number of year in rotation
!!    nyskip      |none           |number of years to skip output summarization/
!!                                |printing
!!    phuacc(:)   |none           |fraction of plant heat units accumulated
!!    plt_hml(:)  |kg/ha          |amount of metal stored in plant
!!    sol_nly(:)  |none           |number of soil layers in profile
!!    sol_z(:,:)  |mm             |depth to bottom of soil layer
!!    uobp        |none           |phosphorus uptake normalization parameter
!!                                |This variable normalizes the phosphorus
!!                                |uptake so that the model can easily verify
!!                                |that uptake from the different soil layers
!!                                |sums to 1.0
!!    wshd_pup    |kg P/ha        |average annual amount of plant uptake of 
!!                                |phosphorus 
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    plt_hml(:)  |kg/ha         |amount of metal stored in plant
!!    pltfr_p(:)  |none          |fraction of plant biomass that is phosphorus
!!    pplnt(:)    |kg P/ha       |plant uptake of phosphorus in HRU for the day
!!    sol_solp(:,:)|kg P/ha      |amount of phosohorus stored in solution
!!    wshd_pup    |kg P/ha       |average annual amount of plant uptake of 
!!                               |phosphorus
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
     
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    gx          |mm            |lowest depth in layer from which phosphorus
!!                               |may be removed
!!    icrop       |none          |land cover code
!!    ir          |none          |flag for bottom of root zone
!!    j           |none          |HRU number
!!    l           |none          |counter (soil layers)
!!    uapd        |kg P/ha       |plant demand of phosphorus
!!    uapl        |kg P/ha       |amount of phosphorus removed from layer
!!    up2         |kg P/ha       |optimal plant phosphorus content
!!    upmx        |kg P/ha       |maximum amount of phosphorus that can be
!!                               |removed from the soil layer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Min
!!    SWAT: nuts

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: j, icrop, ly, rdl
      real :: dg, dgf
      real :: HMlLtot,SOLtot, zlf
      real :: plt_hml_dmd, hml_plt_ly

      j = 0
      j = ihru
      
      HMLtot = 0
      SOLtot = 0

      icrop = 0
      icrop = idplt(j)
      
      rdl = 0
      do ly = 1, sol_nly(j)
          
        if (sol_z(ly,j) >= sol_rd) then
 !        if (sol_z(ly,j) >= sol_zmx(j)) then 
          rdl = ly 
          exit
        end if
   
      end do
      
      if (rdl > 1) then
          zlf = 0.
      end if
      
      do ly = 1, rdl
          
         if (ly == 1) then
             yy = 0.
         else
             yy = 0.
             yy = sol_z(ly-1,j)
         end if
         dg = 0.
         dg = sol_z(ly,j) - yy ! mm
                
         HMLtot = HMLtot + sol_hml_lab(1,j,ly) + sol_hml_nlab(1,j,ly) ! kg/ha
         SOLtot = SOLtot + dg/1000 * sol_bd(ly,j) * 1000 * 10000 ! kg/ha 
         
      end do
      
      if (bio_ms(j) > 4000.) then
          zlf = 0.
      end if    
      
      plt_hml_dmd = HMLtot / SOLtot * plt_hmlku(1) * bio_ms(j)
      
      hml_plt(j) = plt_hml_dmd - plt_hml(j) 
      
      if (hml_plt(j) < 0.) then
         hml_plt(j) = 0.
      else
          zlf = 0.
      end if
      
      plt_hml(j) = hml_plt(j) + plt_hml(j)
      
      do ly = 1, rdl
          
         if (ly == 1) then
             yy = 0.
         else
             yy = 0.
             yy = sol_z(ly-1,j)
         end if
         dg = 0.
         dg = sol_z(ly,j) - yy ! mm
         dgf = dg / sol_z(rdl,j)
         
         hml_plt_ly = hml_plt(j) * dgf
         
         if (hml_plt_ly > sol_hml_sol(1,j,ly)) then
             ! if metal uptake is large than dissloved metal 
             hml_plt_ly = sol_hml_sol(1,j,ly)
             hml_plt(j) = hml_plt(j)
     &       - (hml_plt_ly - sol_hml_sol(1,j,ly))
             sol_hml_sol(1,j,ly) = 0.      
         else 
             sol_hml_sol(1,j,ly) = sol_hml_sol(1,j,ly)
     &       - hml_plt_ly         
         
        end if     
         
      end do

!! summary calculations
      if (curyr > nyskip) then
        wshd_pup = wshd_pup + hml_plt(j) * hru_dafr(j)
      end if

      return
      end