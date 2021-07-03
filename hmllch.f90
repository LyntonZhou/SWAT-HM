      subroutine hmllch
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates h.metal leached through each layer,
!!    h. metal transported with lateral subsurface flow, and h. metal
!!    transported with surface runoff

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    iflag         |none			  |Flaging whether dissolved ligand-bound metal 
!!				  				  |be ignored (iflag=1), or dissolved free metal
!!								  |cation be ignored (iflag=-1), or hmllch should
!!                                  |not be ignored.
!!    flat(:,:)     |mm H2O         |lateral flow in soil layer on current day
!!    hruhml(:)     |none           |h. metal use flag:
!!                                  | 0: no h. metals used in HRU
!!                                  | 1: h. metals used in HRU
!!    ihru          |none           |HRU number
!!    n_hml_mx      |none           |number of different metals used in the simulation
!!    n_hml_no(:)   |none           |array of unique metals used in watershed
!!    percop        |none           |h. metal percolation coefficient (0-1)
!!                                  |0: concentration of h. metal in surface
!!                                  |   runoff is zero
!!                                  |1: percolate has same concentration of
!!                                  |   h. metal as surface runoff
!!    hml_wsol(:)   |kg/m3          |solubility of h.metal in water, unit is 1000ppm
!!    sol_bd(:,:)   |Mg/m**3        |bulk density of the soil
!!    sol_hmlkp(:,:)|(kg/kg)/(kg/m3)|h. metal sorption coefficient, Kp; the
!!                  |(kg/kg)/(kg/m3)|ratio of the concentration in the solid
!!                  |(m3/kg)        |phase to the concentration in solution
!!    sol_nly(:)    |none           |number of layers in soil profile
!!    sol_por(:,:)  |none           |total porosity of soil layer expressed as
!!                                  |a fraction of the total volume
!!    sol_prk(:,:)  |mm H2O         |percolation from soil layer on current day
!!    sol_hml_sol(:,:,:)  |kg /ha   |amount of metal stored as free ion in solution
!!    sol_hml_lab(:,:,:)  |kg /ha   |amount of metal stored as  labile species in solid
!!    sol_hml_nlab(:,:,:) |kg /ha   |amount of metal stored as non-labile species in solid
!!    frac_sol_cbn_soluble(:,:)     |Fraction of soil layer carbon that is soluble organic ligand
!! 	  point_gamma(:)                |Ligand binding number for specific soluble metal 
!!    sol_wpmm(:,:)|mm H20          |water content of soil at -1.5 MPa (wilting point)
!!    sol_z(:,:)   |mm              |depth to bottom of soil layer
!!    surfq(:)     |mm H2O          |surface runoff generated on day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hml_lat(:,:) |kg pst/ha      |amount of h. metal in lateral flow in HRU
!!                                 |for the day
!!    hml_surq(:,:)|kg/ha          |amount of h. metal type lost in surface
!!                                 |runoff on current day in HRU
!!    hml_prk(:,:)  |kg/ha         |amount of h. metal type leached from soil
!!                                 |profile on current day
!!    sol_hml_sol(:,:,:)  |kg /ha  |amount of metal stored as free ion in solution
!!    sol_hml_lab(:,:,:)  |kg /ha  |amount of metal stored as  labile species in solid
!!    sol_hml_nlab(:,:,:) |kg /ha  |amount of metal stored as non-labile species in solid
!!    hml_zdb(:,:)     |
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    co          |kg/mm-ha      |concentration of h. metal in water
!!    cocalc      |kg/mm-ha      |calc concentration of h. metal in water
!!    csurf       |kg/mm-ha      |concentration of h. metal in surq and latq
!!    dg          |mm            |depth of soil layer
!!    j           |none          |HRU number
!!    k           |none          |counter
!!    kk          |none          |h. metal number from hml.dat
!!    ly          |none          |counter (soil layers)
!!    qsurf       |mm H2O        |surface runoff for layer
!!    vf          |
!!    xx          |kg/ha         |amount of h. metal removed from soil layer
!!    yy          |
!!    zdb1        |mm
!!    gamma
!!    kappa
!!    sol_in_soil(1:2) |	       |Temperory 2-element array for metal, 1st element for soluable, 2nd for ligand bound metal
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Min

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm
	
      integer :: j, ly, k, kk
      real :: dg, yy, qsurf, vf, zdb1, xx, co, csurf, cocalc, SAT
      real :: kp, kL, L, kappa, gamma
      real :: sol_in_soil(1:2)=0.
      real :: HMLtot, zlf

      percop=max(percop,0.9)
      j = 0
      j = ihru
	
      if (hruhml(j) /= 0) then
          do ly = 1, sol_nly(j)
              if (ly == 1) then
                  yy = 0.
              else
                  yy = 0.
                  yy = sol_z(ly-1,j)
              end if
              dg = 0.
              dg = sol_z(ly,j) - yy
              SAT = sol_ul(ly,j)

              do k_i = 1, n_hml_mx
                  k = n_hml_no(k_i)	 ! ID No. of metal in db
                  !! write(*,*) hmlname(k) ! For test
                  !! Sequence No. of input metal, in fact this sequence number is useless
                  kk = i_hml_no(k) 
                  
                  if (hml_eqn == 2) then
                      kp = sol_hmlkp(k,1)
                  else
                      kp = hmlkd(sol_ph(ly,j), sol_cbn(ly,j),ly) / 1000.
                      kp = (1 + sol_hmlkp(k,1)) * kp
                  endif

                  if (kk > 0) then
                      qsurf = 0.
                      if (ly == 1) then
                          qsurf = surfq(j)
                      else
                          qsurf = 0.
                      endif
                      
                      !! debug by Zhou 20151031
                      if (qsurf > 20) then
                          zlf = 0
                      end if
                      
                      ! Total amount of h.metal in soil layer
                      HMLtot = sol_hml_sol(k,j,ly) + sol_hml_lab(k,j,ly)				

                      call hml2EQspecies(X1,X2,HMLtot,SAT,kp,sol_bd(ly,j),dg)

                      sol_hml_sol(k,j,ly) = X1
                      sol_hml_lab(k,j,ly)= X2

                      zdb1 = 0.
                      zdb1 = sol_ul(ly,j) + 1.e3 * kp * sol_bd(1,j) * dg ! Unit: mm
                      !! units: mm + (m^3/ton)*(ton/m^3)*mm = mm
                      if (ly == 1) hml_zdb(k,j) = zdb1

                      vf = 0.
                      if (ly == 1) then
                          vf = qsurf + sol_prk(ly,j) + flat(ly,j)
                      else
                          vf = sol_prk(ly,j) + flat(ly,j)
                      endif

                      sol_in_soil = sol_hml_sol(k,j,ly:ly+1) ! Transport soluble metal

                      if (sol_in_soil(1) >= 0. .and. vf > 0.) then
                          if (ly == sol_nly(j)) then
                              zlf = 0
                          end if 
                          xx = 0.
                          xx = sol_in_soil(1) * (1. - Exp(-vf / (zdb1 + 1.e-6)))
                          cocalc = 0.
                          co = 0.
                          if (ly == 1) then
                              cocalc = xx / (sol_prk(ly,j) + percop * (qsurf + flat(ly,j))+ 1.e-6)
                              ! kg/(mm ha)==0.1 kg/m3
                          else
                              cocalc = xx / (sol_prk(ly,j) + flat(ly,j) + 1.e-6)
                          end if
                          ! Unit: hml_wsol unit: from kg/m3 to 0.1 kg/m3
                          co = Min(hml_wsol(kk)*10, cocalc)

                          !! calculate concentration of h. metal in surface
                          !! runoff and lateral flow
                          csurf = 0.
                          if (ly == 1) then
                              csurf = percop * co
                          else
                              csurf = co
                          end if

                          !! calculate h. metal leaching
                          xx = 0.
                          xx = co * sol_prk(ly,j)
                          if (xx > sol_in_soil(1)) then
                              xx = sol_in_soil(1)
                          end if

                          sol_in_soil(1) = sol_in_soil(1) - xx

                          if (ly < sol_nly(j)) then
                              sol_in_soil(2) = sol_in_soil(2) + xx
                          else
                              hml_prk(k,j) = xx
                          end if

                          !! calculate h. metal lost in surface runoff
                          if (ly == 1) then
                              yy = 0.
                              yy = csurf * surfq(j)
                              if (yy > sol_in_soil(1)) then
                                  yy = sol_in_soil(1)
                              end if
                              sol_in_soil(1) = sol_in_soil(1) - yy
                              hml_surq(k,j) = yy
                          endif

                          !! calculate h. metal lost in lateral flow
                          yy = 0.
                          yy = csurf * flat(ly,j)
                          if (yy > sol_in_soil(1)) then
                              yy = sol_in_soil(1)
                          end if
                          sol_in_soil(1) = sol_in_soil(1) - yy
                          hml_lat(k,j) = hml_lat(k,j) + yy
                          !! Save to soluble metal
                          sol_hml_sol(k,j,ly:ly+1) = sol_in_soil 

                          !! Calculate eqiulibrium again for current layer
                          !! Total amount of h.metal in soil layer
                          HMLtot = sol_hml_sol(k,j,ly) + sol_hml_lab(k,j,ly)				

                          call hml2EQspecies(X1,X2,HMLtot,SAT,kp,sol_bd(ly,j),dg)

                          sol_hml_sol(k,j,ly) = X1
                          sol_hml_lab(k,j,ly)= X2

                      end if

                  end if
              end do
          end do
      end if

      return
      end
