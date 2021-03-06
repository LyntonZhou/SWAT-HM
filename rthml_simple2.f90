    subroutine rthml_simple

    !     ~ ~ ~ PURPOSE ~ ~ ~
    !     this subroutine computes the daily stream h.metal balance
    !     (soluble and sorbed)

    !    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
    !    name          |units         |definition
    !    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !    ch_l2(:)      |km            |length of main channel
    !    ch_w(2,:)     |m             |average width of main channel
    !    chhml_conc(:) |kg/(m**3)     |initial h.metal concentration in reach
    !    chhml_koc(:)  |m**3/g        |h.metal partition coefficient between
    !                                 |water and sediment in reach
    !    chhml_mix(:,:)|m/day         |mixing velocity (diffusion/dispersion) for
    !                                 |h.metal in reach
    !    chhml_rsp(:)  |m/day         |resuspension velocity in reach for h.metal
    !                                 |sorbed to sediment
    !    chhml_stl(:)  |m/day         |settling velocity in reach for h.metal
    !                                 |sorbed to sediment
    !    hru_sub(:)    |none          |subbasin number where reach is located
    !    inum1         |none          |reach number
    !    inum2         |none          |inflow hydrograph storage location number
    !    rchdep        |m             |depth of flow on day
    !    sedrchin      |metric tons   |sediment transported into channel
    !    sedrchout     |metric tons   |sediment transported out of channel
    !    rchwtr        |m^3 H2O       |water stored in reach at beginning of day
    !    rnum1         |none          |fraction of overland flow
    !    rtwtr         |m^3 H2O       |water leaving reach on day
    !    sedhml_act(:) |m             |depth of active sediment layer in reach for
    !                                 |h.metal
    !    sedhml_bry(:) |m/day         |h.metal burial velocity in river bed
    !                                 |sediment
    !    sedhml_rea(:) |1/day         |h.metal reaction coefficient in river bed
    !                                 |sediment
    !    hmlroute(:,:,:)              |h.metal amount in aqueous input
    !    sedhml_conc(:,:,:)|kg/m3     |h.metal amount in sediment
    !    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    !    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
    !    name        |units         |definition
    !    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !    bury           |mg hml     |loss of h.metal from active sediment layer
    !                               |by burial
    !    difus_hml(k)   |mg hml     |diffusion of h.metal from sediment to reach
    !    rsushml(k)     |mg hml     |amount of h.metal moving from sediment to
    !                               |reach due to resuspension
    !    setlhml(k)     |mg hml     |amount of h.metal moving from water to
    !                               |sediment due to settling
    !    hml_sol_o(k)	|kg/m3      |soluble hml concentration in outflow on day
    !    hml_lig_o(k)	|kg/m3      |ligand-bnd hml concentration in outflow on day
    !    hml_exch_o(k)	|kg/m3      |exchangeable hml concentration in outflow on day
    !    hml_nlab_o(k)	|kg/m3      |non-labile hml concentration in outflow on day
    !    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    !    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
    !    name        |units         |definition
    !    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !    bedvol      |m^3           |volume of river bed sediment
    !    depth       |m             |depth of water in reach
    !    frsol       |none          |fraction of h.metal in reach that is soluble
    !    frsrb       |none          |fraction of h.metal in reach that is sorbed
    !    jrch        |none          |reach number
    !    sedcon      |g/m^3         |sediment concentration
    !    rch_hml_sol_in    |kg hml  |soluble h.metal entering reach during
    !                               |time step
    !    rch_hml_exch_in   |kg hml  |sorbed h.metal entering reach during
    !							    |time step
    !    tday        |days          |flow duration
    !    wtrin       |m^3 H2O       |volume of water entering reach during time
    !                               |step
    !    setlhml_exch(:)            |Amount settled
    !    rsushml_exch(:)            |Amount resuspended
    !    bury_exch                  |Amount buried as exchangeable species
    !    bury_nlab                  |Amount buried as nonlabile species
    !    sed_por                    |porosity of river bed sediment
    !    sol_pd                     |soil particle density
    !    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    !    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
    !    Intrinsic: Abs

    !    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

    use parm

    integer :: jrch
    real :: rch_hml_sol_in, rch_hml_exch_in, rch_hml_nlab_in
    real :: depth, chhmlmass, frsol, frsrb
    real :: bedvol, bedmass, wtrin, solmax, sedcon, tday
    real :: rch_hml_exch(n_hml_mx),rch_hml_nlab(n_hml_mx)
    real :: rch_hml_sol(n_hml_mx)
    real :: sed_hml_exch(n_hml_mx),sed_hml_nlab(n_hml_mx)
    real :: sed_hml_sol(n_hml_mx)
    real :: rch_hml_2,sed_hml_2,x1,x2,y1,y2
    real :: kp1, kp2
    real :: sed_por, sol_pd
    real :: rsp, stl
    
    jrch = 0
    jrch = inum1

    setlhml_exch(1:2) = 0.
    setlhml_nlab(1:2) = 0.
    rsushml_exch(1:2) = 0.
    rsushml_nlab(1:2) = 0.
    
    ! assume por=0.5, density=2650 kg/m3 in river bed sediment 
    sol_pd = 2.65e3 
    sed_por = 0.5
    
    ! added by Zhou for debug
    if (jrch == 48) then
        zlf = 0.
    end if
    
    ! initialize depth of water for metal calculations
    depth = 0.
    if (rchdep < 0.1) then
        depth = .1
    else
        depth = rchdep
    endif

    ! calculate volume of active river bed sediment layer
    bedvol = 0.
    bedvol = ch_w(2,jrch) * ch_l2(jrch) * 1000. * sedhml_act(jrch)
    
    ! calculate mass of active river bed sediment layer
    bedmass = 0.
    bedmass = bedvol * sol_pd / 1000 * (1-sed_por)

    ! calculate volume of water entering reach
    wtrin = 0.
    wtrin = varoute(2,inum2) * (1. - rnum1)

    ! metal loop start
    do k = 1, n_hml_mx
        
        ! partition coefficient in water and river bed for specific metal
        kp1 = sol_hmlkp(k,2)
        kp2 = sol_hmlkp(k,3)
        
        ! metal entering reach
        rch_hml_sol_in = hmlroute(1,k,inum2) * (1. -rnum1)
        rch_hml_exch_in = hmlroute(3,k,inum2) * (1. -rnum1)
        rch_hml_nlab_in = hmlroute(4,k,inum2) * (1. -rnum1)
        
        ! add metal point source 
        if (pointsourcehmlihout > 0) then
            if(varoute(19+k,pointsourcehmlihout) > 0) then
                rch_hml_sol_in = rch_hml_sol_in + varoute(19+k,pointsourcehmlihout)
                varoute(19+k,pointsourcehmlihout) = 0
                pointsourcehmlihout = 0
            end if
        end if

        ! calculate mass of metal in reach
        rch_hml_sol(k) = rch_hml_sol_in + chhml_conc(1,k,jrch) * rchwtr
        rch_hml_exch(k) = rch_hml_exch_in + chhml_conc(3,k,jrch) * rchwtr
        rch_hml_nlab(k) = rch_hml_nlab_in + chhml_conc(4,k,jrch) * rchwtr

        ! calculate mass of metal in active sediment layer
        sed_hml_sol(k) = sedhml_conc(1,k,jrch) * bedvol * sed_por
        sed_hml_exch(k) = sedhml_conc(3,k,jrch) * bedvol *(1-sed_por)
        sed_hml_nlab(k) = sedhml_conc(4,k,jrch) * bedvol *(1-sed_por)

        if (rch_hml_sol(k) + rch_hml_exch(k) + rch_hml_nlab(k) < 1.e-10 ) then
            chhml_conc(:,k,jrch) = 0.
        end if
        if (sed_hml_sol(k) + sed_hml_exch(k) + sed_hml_nlab(k) < 1.e-10) then
            sedhml_conc(:,k,jrch) = 0.
        end if

        !in-stream processes
        if (rtwtr / 86400. > 0.002) then ! If outflow exists
            ! calculate sediment concentration
            sedcon = 0.
            sedcon = sedrch / rtwtr * 1.e3  ! Unit: ton/m3 * 1e3 = kg/m3
                        
            ! calculate metal partition in water
            call hml2EQspecies_reach(x1,x2,kp1,sedcon)
            rch_hml_2 = rch_hml_sol(k) + rch_hml_exch(k)
            rch_hml_sol(k) = rch_hml_2 * x1
            rch_hml_exch(k) = rch_hml_2 * x2
            
            ! calculate metal partition in sed
            call hml2EQspecies_sed(y1,y2,kp2,sed_por,sol_pd)        
            sed_hml_2 = sed_hml_sol(k) + sed_hml_exch(k)
            sed_hml_sol(k) = sed_hml_2 * y1
            sed_hml_exch(k) = sed_hml_2 * y2
            
            ! calculate flow duration
            tday = 0.
            tday = rttime / 24.0
            if (tday > 1.0) tday = 1.0
            
            ! sediment processes
            if (sedrchin > sedrchout) then
                ! settling
                setlhml_exch(k)=(sedrchin-sedrchout)/sedrchin*rch_hml_exch(k)
                setlhml_nlab(k)=(sedrchin-sedrchout)/sedrchin*rch_hml_nlab(k)
                ! check settling veolocity
                !stl = (sedrchin-sedrchout)/sedrchin*depth/tday
                !if (stl > 2.5) then ! 10 m/s
                !    zlf = 0.
                !    stl = 10.
                !    setlhml_exch(k)=2.5*rch_hml_exch(k)*tday/depth
                !    setlhml_nlab(k)=2.5*rch_hml_nlab(k)*tday/depth
                !endif
            else
                ! resuspension
                rsushml_exch(k)=(sedrchout-sedrchin)/bedmass*sed_hml_exch(k)
                rsushml_nlab(k)=(sedrchout-sedrchin)/bedmass*sed_hml_nlab(k)
                 ! check resuspension veolocity
                rsp = (sedrchin-sedrchout)/bedmass*sedhml_act(jrch)/tday
                if (rsp > 50) then
                    zlf = 0.
                endif
            endif

            !calculate amount of metal removed from reach by settling
            !setlhml_exch(k)=chhml_stl(k,jrch)*rch_hml_exch(k)*tday/depth
            if (setlhml_exch(k) > rch_hml_exch(k)) then
                setlhml_exch(k) = rch_hml_exch(k)
                rch_hml_exch(k) = 0
            else
                rch_hml_exch(k) = rch_hml_exch(k) - setlhml_exch(k)
            end if
            sed_hml_exch(k) = sed_hml_exch(k) + setlhml_exch(k)

            !setlhml_nlab(k)=chhml_stl(k,jrch)*rch_hml_nlab(k)*tday/depth
            if (setlhml_nlab(k) > rch_hml_nlab(k)) then
                setlhml_nlab(k) = rch_hml_nlab(k)
                rch_hml_nlab(k) = 0.
            else
                rch_hml_nlab(k) = rch_hml_nlab(k) - setlhml_nlab(k)
            end if
            sed_hml_nlab(k) = sed_hml_nlab(k) + setlhml_nlab(k)

            ! calculate resuspension of metal in reach
            !rsushml_exch(k)=chhml_rsp(k,jrch)*sed_hml_exch(k)*tday/ sedhml_act(jrch)
            if (rsushml_exch(k) > sed_hml_exch(k)) then
                rsushml_exch(k) = sed_hml_exch(k)
                sed_hml_exch(k) = 0
            else
                sed_hml_exch(k) = sed_hml_exch(k) - rsushml_exch(k)
            end if
            rch_hml_exch(k) = rch_hml_exch(k) + rsushml_exch(k)

            !!rsushml_nlab(k)=chhml_rsp(k,jrch)*sed_hml_nlab(k)*tday/sedhml_act(jrch)
            if (rsushml_nlab(k) > sed_hml_nlab(k)) then
                rsushml_nlab(k) = sed_hml_nlab(k)
                sed_hml_nlab(k) = 0.
            else
                sed_hml_nlab(k) = sed_hml_nlab(k) - rsushml_nlab(k)
            end if
            rch_hml_nlab(k) = rch_hml_nlab(k) + rsushml_nlab(k)

            ! calculate diffusion of dissolved metal between reach water and active sediment layer
            !difus_hml(k) = chhml_mix(k,jrch) * tday * (sed_hml_sol(k) / sedhml_act(jrch) - rch_hml_sol(k) / depth)
            difus_hml(k) = chhml_mix(k,jrch) * tday * (sed_hml_sol(k) / sedhml_act(jrch) - rch_hml_sol(k) / depth)
            if (difus_hml(k) > 0.) then
                if (difus_hml(k) > sed_hml_sol(k)) then
                    difus_hml(k) = sed_hml_sol(k)
                    sed_hml_sol(k) = 0.
                else
                    sed_hml_sol(k) = sed_hml_sol(k) - Abs(difus_hml(k))
                end if
                rch_hml_sol(k) = rch_hml_sol(k) + Abs(difus_hml(k))
            else
                if (Abs(difus_hml(k)) > rch_hml_sol(k)) then
                    difus_hml(k) = -rch_hml_sol(k)
                    rch_hml_sol(k) = 0.
                else
                    rch_hml_sol(k) = rch_hml_sol(k) - Abs(difus_hml(k))
                end if
                sed_hml_sol(k) = sed_hml_sol(k) + Abs(difus_hml(k))
            end if

            ! calculate removal of metal from active sediment layer by burial
            bury_exch(k)=sedhml_bry(k,jrch)*sed_hml_exch(k)/sedhml_act(jrch)
            if (bury_exch(k) > sed_hml_exch(k)) then
                bury_exch(k) = sed_hml_exch(k)
                sed_hml_exch(k) = 0.
            else
                sed_hml_exch(k) = sed_hml_exch(k) - bury_exch(k)
            end if

            bury_nlab(k)=sedhml_bry(k,jrch)*sed_hml_nlab(k)/sedhml_act(jrch)
            if (bury_nlab(k) > sed_hml_nlab(k)) then
                bury_nlab(k) = sed_hml_nlab(k)
                sed_hml_nlab(k) = 0.
            else
                sed_hml_nlab(k) = sed_hml_nlab(k) - bury_nlab(k)
            end if

            ! verify that water concentration is at or below solubility
            solmax = 0.
            solmax = hml_wsol(k) * (rchwtr + wtrin)	! kg/m3 * m3 = kg
            sol_aqueous = rch_hml_sol(k)
            if (solmax < sol_aqueous) then
                ! those beyond solubility is settled onto solid as labile metal
                sed_hml_exch(k) = sed_hml_exch(k) + (sol_aqueous - solmax)
                rch_hml_sol(k) = rch_hml_sol(k) - (sol_aqueous - solmax)
            end if
            
            rch_hml_2 = rch_hml_sol(k) + rch_hml_exch(k)
            rch_hml_sol(k) = rch_hml_2 * x1
            rch_hml_exch(k) = rch_hml_2 * x2
            
            sed_hml_2 = sed_hml_sol(k) + sed_hml_exch(k)
            sed_hml_sol(k) = sed_hml_2 * y1
            sed_hml_exch(k) = sed_hml_2 * y2

        else	! if insignificant outflow
            ! calculate metal partition in sed
            call hml2EQspecies_sed(y1,y2,kp2,sed_por,sol_pd)
            
            sed_hml_sol(k) = sed_hml_sol(k) + rch_hml_sol(k)
            sed_hml_exch(k) = sed_hml_exch(k) + rch_hml_exch(k)
            sed_hml_nlab(k) = sed_hml_nlab(k) + rch_hml_nlab(k)
            rch_hml_sol(k) = 0.
            rch_hml_exch(k) = 0.
            rch_hml_nlab(k) = 0.

            sed_hml_2 = sed_hml_sol(k) + sed_hml_exch(k)
            sed_hml_sol(k) = sed_hml_2 * y1
            sed_hml_exch(k) = sed_hml_2 * y2

        end if
        
        !! sediment processes
        !! calculate transformation of metal in bed sediment by ageing reaction
        
        
        ! calculate metal concentrations at end of day
        chhml_conc(:,k,jrch) = 0.
        sedhml_conc(:,k,jrch) = 0.
        chhml_conc(1,k,jrch) = rch_hml_sol(k) / (rchwtr + wtrin)
        chhml_conc(3,k,jrch) = rch_hml_exch(k) / (rchwtr + wtrin)
        chhml_conc(4,k,jrch) = rch_hml_nlab(k) / (rchwtr + wtrin)
        sedhml_conc(1,k,jrch) = sed_hml_sol(k) / (bedvol * sed_por)
        sedhml_conc(3,k,jrch) = sed_hml_exch(k) / (bedvol * (1-sed_por))
        sedhml_conc(4,k,jrch) = sed_hml_nlab(k) / (bedvol * (1-sed_por))
        
        ! calculate outflow metal concentration
        if (rtwtr > 0.001) then
            hml_sol_o(k)  = chhml_conc(1,k,jrch)
            hml_exch_o(k) = chhml_conc(3,k,jrch)
            hml_nlab_o(k) = chhml_conc(4,k,jrch)
            
            if (hml_sol_o(k) + hml_exch_o(k) + hml_nlab_o(k) < 0.) then
                hml_sol_o(k) = 0.
                hml_exch_o(k) = 0.
                hml_nlab_o(k) = 0.
            end if

        else
            hml_sol_o(k) = 0.
            hml_exch_o(k) = 0.
            hml_nlab_o(k) = 0.
        end if

    end do ! metal loop end
    
    return
    end
