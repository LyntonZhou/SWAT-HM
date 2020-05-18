      subroutine lakehml

!    ~ ~ ~ PURPOSE ~ ~ ~
!    this subroutine computes the lake hydrologic h.metal balance.
!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!    name          |units         |definition
!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!    inum1         |none          |reservoir number
!    lkhml_conc(:) |mg/m^3        |h.metal concentration in lake water
!    lkhml_koc(:)  |m**3/g        |h.metal partition coefficient between
!                                 |water and sediment in lake water
!    lkhml_mix(:)  |m/day         |mixing velocity (diffusion/dispersion) in
!                                 |lake water for h.metal
!    lkhml_rea(:)  |1/day         |h.metal reaction coefficient in lake water
!    lkhml_rsp(:)  |m/day         |resuspension velocity in lake water for
!                                 |h.metal sorbed to sediment
!    lkhml_stl(:)  |m/day         |settling velocity in lake water for
!                                 |h.metal sorbed to sediment
!    lkhml_vol(:)  |m/day         |h.metal volatilization coefficient in lake
!                                 |water
!    lkshml_act(:) |m             |depth of active sediment layer in lake for
!                                 |for h.metal
!    lkshml_bry(:) |m/day         |h.metal burial velocity in lake bed
!                                 |sediment
!    lkshml_conc(:)|mg/m^3        |h.metal concentration in lake bed sediment
!    lkshml_rea(:) |1/day         |h.metal reaction coefficient in lake bed
!                                 |sediment
!    res_sed(:)    |kg/L (ton/m^3)|amount of sediment in reservoir
!    res_vol(:)    |m^3 H2O       |reservoir volume
!    resflwo       |m^3 H2O       |water leaving reservoir on day
!    ressa         |ha            |surface area of reservoir on day
!    ressedo       |metric tons   |sediment leaving reservoir during time step
!    solhmli      |mg hml        |soluble h.metal entering reservoir
!    sorhmli      |mg hml        |sorbed h.metal entering reservoir
!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!    name          |units         |definition
!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!    bury          |mg hml        |loss of h.metal from active sediment layer
!                                 |by burial
!    difus         |mg hml        |diffusion of h.metal from sediment to lake
!                                 |water
!    lkhml_conc(:) |mg/m^3        |h.metal concentration in lake water
!    lkshml_conc(:)|mg/m^3        |h.metal concentration in lake bed sediment
!    reactw        |mg hml        |amount of h.metal in lake water lost
!                                 |through reactions
!    reactb        |mg hml        |amount of h.metal in sediment that is lost
!                                 |through reactions
!    reshmli      |mg hml        |h.metal entering reservoir on day
!    resushml      |mg hml        |amount of h.metal moving from sediment to
!                                 |lake water due to resuspension
!    setlhml       |mg hml        |amount of h.metal moving from water to
!                                 |sediment due to settling
!    solhmlo      |mg hml        |soluble h.metal in outflow on day
!    sorhmlo      |mg hml        |sorbed h.metal in outflow on day
!    volathml      |mg hml        |amount of h.metal lost from lake water
!                                 |by volatilization
!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!    name        |units         |definition
!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!    dlake       |m             |depth of water in reservoir
!    fd1         |none          |fraction of h.metal in water that is soluble
!    fd2         |none          |fraction of h.metal in sediment that is
!                               |soluble
!    fp1         |none          |fraction of h.metal in water that is sorbed
!    fp2         |none          |fraction of h.metal in sediment that is 
!                               |sorbed
!    jres        |none          |reservoir number
!    res_hml_sol    |kg hml         |amount of h.metal in lake watertable as soluble species
!    sed_hml_sol    |kg hml         |amount of h.metal in lake sediment as soluble species
!    sed_por                        |porosity of river bed sediment
!    sol_pd                         |soil particle density
!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!    Intrinsic: Abs

!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

    use parm

    integer :: jres
    real :: fd1, fp1, fd2, dlake, fp2
    integer :: k
    real :: res_hml_2, sed_hml_2, sedcon
	real :: kp1, kp2
    real ::  sed_por, sol_pd
      

	real :: res_hml_exch(n_hml_mx),res_hml_nlab(n_hml_mx)
	real :: res_hml_sol(n_hml_mx)
	real :: sed_hml_exch(n_hml_mx),sed_hml_nlab(n_hml_mx)
	real :: sed_hml_sol(n_hml_mx)
	real :: x1,x2,y1,y2,difus_hml_tmp
    

      jres = 0
      jres = inum1
      
      sed_por = 0.8 
      sol_pd =2.65e3
      
      !hml_lkwtr = 0.
      !hml_lksed = 0.
      !hml_lkwtr = lkhml_mass(jres)
      !hml_lksed = lkshml_mass(jres)
      !lkshml_act(jres)=lkspst_act(jres)
	  bedvol=lkshml_act(jres) * ressa * 10000.

      if (res_vol(jres) > 1.) then
        ! calculate depth of lake
        dlake = 0.
        dlake = res_vol(jres) / (ressa * 10000.)
	 else
	    return
	 end if
! h.metal transported into lake/reservior during day
	do k = 1, n_hml_mx

	    ! Fundamental parameters for specific metal
		kp1 = sol_hmlkp(k,2)
        kp2 = sol_hmlkp(k,3)
		
     !   fd1 = 0.
     !   fp1 = 0.
     !   fd2 = 0.
     !   fp2 = 0.
     !   fd1 = 1. / (1. + Kp1 * res_sed(jres) * 1.e3)	
     !   ! Note Kp has unit of kg/m3, different from Kp(pest) which has the unit of g/m3
     !   fp1 = 1. - fd1
     !   ! ASSUME POR=0.8; DENSITY=2.6E6, then concsed = 5.2e5; KD2=KD1
     !   fd2 = 1. / (.8 + 5.2e2 * Kp2)	! Note the unit of Kp
	    !fp2 = 1. - fd2

		! Incoming metal
		!!hmlin = 0.
		res_hml_sol_in = hmlroute(1,k,inum2)
        res_hml_exch_in = hmlroute(3,k,inum2)
		res_hml_nlab_in = hmlroute(4,k,inum2)

		! calculate mass of h.metal in lake/reservior
		res_hml_sol(k) = res_hml_sol_in + lkhml_conc(1,k,jres) * res_vol(jres)
		res_hml_exch(k) = res_hml_exch_in + lkhml_conc(3,k,jres) * res_vol(jres)
		res_hml_nlab(k) = res_hml_nlab_in + lkhml_conc(4,k,jres) * res_vol(jres)

		! calculate mass of h.metal in bed sediment
		sed_hml_sol(k) = lksedhml_conc(1,k,jres) * bedvol*sed_por
 		sed_hml_exch(k) = lksedhml_conc(3,k,jres) * bedvol*(1-sed_por)
		sed_hml_nlab(k) = lksedhml_conc(4,k,jres) * bedvol*(1-sed_por)


		if (res_hml_sol(k) + res_hml_exch(k) + res_hml_nlab(k) < 1.e-10 ) then 
			lkhml_conc(:,k,jres) = 0.
		end if
     	    if (sed_hml_sol(k) + sed_hml_exch(k) + sed_hml_nlab(k) < 1.e-10) then
		    lksedhml_conc(:,k,jres) = 0.
		end if


        ! determine amount of h.metal settling to sediment layer
	      !lkhml_stl=lkpst_stl
	      !lkhml_rsp=lkpst_rsp
	      !lkhml_mix=lkpst_mix
	      !lksedhml_bry=lkspst_bry

      setlhml_exch(k)=lkhml_stl(jres)*res_hml_exch(k)/dlake
	  setlhml_nlab(k)=lkhml_stl(jres)*res_hml_nlab(k)/dlake
        if (setlhml_nlab(k) > res_hml_nlab(k)) then
          setlhml_nlab(k) = res_hml_nlab(k)
          res_hml_nlab(k) = 0.
	      setlhml_exch(k) = res_hml_exch(k)
          res_hml_exch(k) = 0
        else
          res_hml_nlab(k) = res_hml_nlab(k) -setlhml_nlab(k)
	      res_hml_exch(k) = res_hml_exch(k) -setlhml_exch(k)
        end if
      sed_hml_exch(k) = sed_hml_exch(k) + setlhml_exch(k)
      sed_hml_nlab(k) = sed_hml_nlab(k) + setlhml_nlab(k)

        ! determine h.metal resuspended into lake water
        rsushml_exch(k)=lkhml_rsp(jres)*sed_hml_exch(k)/lkshml_act(jres)
        rsushml_nlab(k)=lkhml_rsp(jres)*sed_hml_nlab(k)/lkshml_act(jres)
        if (rsushml_nlab(k) > sed_hml_nlab(k)) then
          rsushml_nlab(k) = sed_hml_nlab(k)
          sed_hml_nlab(k) = 0.
	      rsushml_exch(k) = sed_hml_exch(k)
          sed_hml_exch(k) = 0
        else
          sed_hml_nlab(k) = sed_hml_nlab(k) -rsushml_nlab(k)
	      sed_hml_exch(k) = sed_hml_exch(k) -rsushml_exch(k)
        end if
        res_hml_exch(k) = res_hml_exch(k) + rsushml_exch(k)
	    res_hml_nlab(k) = res_hml_nlab(k) + rsushml_nlab(k)

        ! determine h.metal diffusing from sediment to water, including soluble and ligand-bnd species.
        difus_hml(k) = lkhml_mix(jres) *(sed_hml_sol(k)/lkshml_act(jres) - res_hml_sol(k)/dlake )      
        if (difus_hml(k) > 0.) then
          if (difus_hml(k) > sed_hml_sol(k)) then
            difus_hml(k) = sed_hml_sol(k)
            sed_hml_sol(k) = 0.
          else
            sed_hml_sol(k) = sed_hml_sol(k) - Abs(difus_hml(k))
          end if
          res_hml_sol(k) = res_hml_sol(k) + Abs(difus_hml(k))
        else
          if (Abs(difus_hml(k)) > res_hml_sol(k)) then
            difus_hml(k) = -res_hml_sol(k)
            res_hml_sol(k) = 0.
          else
            res_hml_sol(k) = res_hml_sol(k) - Abs(difus_hml(k))
          end if
          sed_hml_sol(k) = sed_hml_sol(k) + Abs(difus_hml(k))
        end if


        ! determine h.metal lost from sediment by burial
	  bury_exch(k)=lksedhml_bry(jres)*sed_hml_exch(k)/lkshml_act(jres)
	  bury_nlab(k)=lksedhml_bry(jres)*sed_hml_nlab(k)/lkshml_act(jres)
        if (bury_nlab(k) > sed_hml_nlab(k)) then
          bury_nlab(k) = sed_hml_nlab(k)
          sed_hml_nlab(k) = 0.
          bury_exch(k) = sed_hml_exch(k)
          sed_hml_exch(k) = 0.
        else
          sed_hml_nlab(k) = sed_hml_nlab(k) - bury_nlab(k)
	      sed_hml_exch(k) = sed_hml_exch(k) - bury_exch(k)
        end if
        
        ! calculate soluble h.metal transported out of reservoir
        !solhmlo = resflwo * fd1 * hml_lkwtr / res_vol(jres)
        !if (solhmlo > hml_lkwtr) then
        !  solhmlo = hml_lkwtr
        !  hml_lkwtr = 0.
        !else
        !  hml_lkwtr = hml_lkwtr - solhmlo
        !end if

        ! calculate sorbed h.metal transported out of reservoir
        !sorhmlo = resflwo * fp1 * hml_lkwtr / res_vol(jres)
        !if (sorhmlo > hml_lkwtr) then
        !  sorhmlo = hml_lkwtr
        !  hml_lkwtr = 0.
        !else
        !  hml_lkwtr = hml_lkwtr - sorhmlo
        !end if

        ! update concentration of h.metal in lake water and sediment

	if (res_hml_exch(k) + res_hml_nlab(k) < 1.e-10) then
		res_hml_exch(k) =0.
		res_hml_nlab(k) =0.
	end if
	if (res_hml_sol(k) < 1.e-10) then
		res_hml_sol(k) =0.
	end if

! Equilibrium striking in terms of water column and sediment bed respectively
	x1 = 0.
	x2 = 0.
	sedcon = res_sed(jres) * 1.e3	! Sediment concentration in lake, in unit of kg/m3
	call hml2EQspecies_reach(x1,x2,Kp1,sedcon)
	
	res_hml_2=res_hml_sol(k) + res_hml_exch(k)
	res_hml_sol (k)=res_hml_2 * x1
	res_hml_exch(k)=res_hml_2 * x2
	
	! Convert mass in reservior water into concentration 
	lkhml_conc(1, k,jres) = res_hml_sol(k) / res_vol(jres)
	lkhml_conc(3, k,jres) = res_hml_exch(k)/ res_vol(jres)
	lkhml_conc(4, k,jres) = res_hml_nlab(k)/ res_vol(jres)


	y1 = 0.
	y2 = 0.
	! ASSUME POR=0.8; DENSITY=2.6E6, then concsed = 5.2e5 g/m3 in lake sediment bed. 
	!sedcon=2.6e3*(1-0.8)	! Sediment concentration in lake bed layer, in unit of kg/m3
	call hml2EQspecies_sed(y1,y2,Kp2,sed_por,sol_pd)
	
	sed_hml_2=sed_hml_sol(k) + sed_hml_exch(k)
	sed_hml_sol (k)=sed_hml_2 * y1
	sed_hml_exch(k)=sed_hml_2 * y2
	
	lksedhml_conc(1,k,jres) = sed_hml_sol(k)/(lkshml_act(jres)* sed_por * ressa * 10000. + 1.)	! Unit: kg/m3
	lksedhml_conc(3,k,jres) = sed_hml_exch(k)/(lkshml_act(jres)*(1-sed_por) * ressa * 10000. + 1.)
    lksedhml_conc(4,k,jres) = sed_hml_nlab(k)/(lkshml_act(jres)*(1-sed_por) * ressa * 10000. + 1.)
      if (rtwtr > 0.001) then
        hml_sol_o(k)  = lkhml_conc(1, k,jres)
        hml_exch_o(k) = lkhml_conc(3, k,jres)
        hml_nlab_o(k) = lkhml_conc(4, k,jres)
      else
        hml_sol_o(k) = 0.
        hml_exch_o(k) = 0.
        hml_nlab_o(k) = 0.
      end if

	end do

      return
      end



