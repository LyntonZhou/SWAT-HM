      subroutine readhml

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine reads data from the HRU/subbasin heavy metal source input
!!    file (.hml). This file contains initial amounts of heavy metals in the
!!    mining pit/tailing pile. in the first soil layer. (Specifics about the first soil layer are given
!!    in the .sol file.) All data in the .hml file is optional input.

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ihru        |none          |HRU number
!!    mhml        |none          |maximum number of metals involved in
!!                               |watershed
!!    n_hml_mx    |none          |number of different metals used in
!!                               |the simulation
!!    n_hml_no(:) |none          |array of unique heavy metals used in
!!                               |watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hruhml(:)     |none          |heavy metal involvement flag:
!!                                 | 0: no metals involved in HRU
!!                                 | 1: metals involved in HRU
!!    hml_fr(:,:)   |km2/km2       |fraction of mine/tailing/piling area contained in a HRU	
!!    i_hml_no(:)	  |none			 |sequence number of metal in N_hml_no(:)
!!    n_hml_no(:)   |none          |array of unique heavy metals used in watershed
!!    n_hml_mx      |none          |number of different metals used in the simulation
!!    sol_hml(:,:,1)|mg/kg         |metal concentration in soil
!!    sol_hml_lab(:,:,1)  |mg/kg   |concentration of metal stored as  labile species in solid in 1st soil layer
!!    sol_hml_nlab(:,:,1) |mg/kg   |concentration of metal stored as non-labile species in solid in 1st soil layer
!!    hml_enr(:,:)  |none          |heavy metal enrichment ratio
!!    gw_sol_hml(:) |ug/L          |soluble metal concentration in groundwater loading to reach
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    eof         |none          |end of file flag
!!    j           |none          |counter
!!    k           |none          |counter
!!    newhml      |none          |new heavy metal flag
!!    hmlfr		|km2/km2       |fraction of mine/tailing/piling area in HRU	
!!    hmlrock     |kg/ha         |heavy metal in rock to be weathered, kg HML per hac
!!    hmlenr      |none          |heavy metal enrichment ratio
!!    hmlnum      |none          |heavy metal number
!!    solhml_ex   |mg/kg         | labile Metal concentration in 1st layer soil 
!!    solhml_nl   |mg/kg         |Non-labile Metal concentration in 1st layer soil 
!!    titldum     |NA            |title line for .hml file
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


      use parm

      character (len=80) :: titldum
      integer :: j, eof, k, newhml, hmlnum
      real :: hmlfr, hmlrock solhml_ex, solhml_nl, hmlenr 
      real :: hmlagrtotal, hmlagrfrac
      real :: hmldeptotal, hmldepfrac
      real :: gwsolhml

      eof = 0
      
      do
	  read (32,5000,iostat=eof) titldum
        if (eof < 0) exit
        read (32,5000,iostat=eof) titldum
        if (eof < 0) exit

        do j = 1, mhml
        hmlnum = 0
	  hmlfr  = 0.
	  hmlrock = 0.
        solhml_ex = 0.
	  solhml_nl = 0.
        hmlenr = 0.
        hmlagrtotal = 0.
        hmlagrfrac = 0.
        hmldeptotal = 0.
        hmldepfrac = 0.
        gwsolhml = 0.
        read (32,*,iostat=eof) hmlnum
        if (eof < 0) exit
        read (32,*,iostat=eof) hmlfr
        if (eof < 0) exit
        read (32,*,iostat=eof) hmlrock
        if (eof < 0) exit
        read (32,*,iostat=eof) solhml_ex
        if (eof < 0) exit
        read (32,*,iostat=eof) solhml_nl
        if (eof < 0) exit
        read (32,*,iostat=eof) hmlenr
        if (eof < 0) exit
        read (32,*,iostat=eof) hmlagrtotal
        if (eof < 0) exit
        read (32,*,iostat=eof) hmlagrfrac
        if (eof < 0) exit
        read (32,*,iostat=eof) hmldeptotal
        if (eof < 0) exit
        read (32,*,iostat=eof) hmldepfrac
        if (eof < 0) exit
        read (32,*,iostat=eof) gwsolhml
        if (eof < 0) exit
        
        if (hmlnum > 0) then
          hruhml(ihru) = 1
          newhml = 0
          do k = 1, n_hml_mx
            if (hmlnum == n_hml_no(k)) then
              newhml = 1		! Matching an existing metal in metal db.
              exit
            endif
          end do
          if (newhml == 0) then ! Find a new metal not yet included in metal db.
            n_hml_no(n_hml_mx) = hmlnum
            i_hml_no(hmlnum) = n_hml_mx
            n_hml_mx = n_hml_mx + 1
          end if
          if (ihru == 1364) then
              zlf = 0
          endif
          k = 0
          k = i_hml_no(hmlnum)
	    hml_fr(k,ihru) = hmlfr
		hml_rock(k,ihru) = hmlrock
          sol_hml_lab(k,ihru,1) = solhml_ex
	    sol_hml_nlab(k,ihru,1) = solhml_nl
          hml_enr(k,ihru) = hmlenr
          hml_agr_total(ihru) = hmlagrtotal/1000.
          hml_agr_frac(ihru) = hmlagrfrac
          hml_dep_total(ihru) = hmldeptotal/1000.
          hml_dep_frac(ihru) = hmldepfrac
          gw_sol_hml(ihru) = gwsolhml
        end if

      if (eof < 0) exit
      end do
      exit
      end do

      close (32)

      return
 5000 format (a)
 5100 format (27x,10f12.2)
      end
