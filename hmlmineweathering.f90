	subroutine hmlmineweathering

	! Calculation of heavy metals weathering in source (mine/tailing/pile)
	!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hru_ra(:)      |MJ/m^2        |solar radiation for the day in HRU
!!    precipday      |mm H2O        |precipitation for the day in HRU
!!    tmn(:)         |deg C         |minimum temperature for the day in HRU
!!    tmpav(:)       |deg C         |average temperature for the day in HRU
!!    tmx(:)         |deg C         |maximum temperature for the day in HRU
!!    rhd(:)         |none          |relative humidity for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    nhml(:)    |none           |sequence number of h. metal application
!!                               |within the year 
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |HRU number
!!    jj          |none          |subbasin number
!!    k           |none          |sequence number of h. metal in N_hml_NO(:)
!!    kk          |none          |h. metal identification number from metal.dat
!!    xx          |kg/ha         |amount of h. metal input to HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: Erfc

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

    use parm

    integer :: j, kk, k, jj, src, zlf
    real :: ra, t_av, weth_rate, Tadj, f

    j = 0
    j = ihru
    f = 1.0

    !! initialize local variables
    kk = 0
    k = 0
    jj = 0
    t_av =tmpav(j)

    do kk = 1, n_hml_mx
        k = n_hml_no(kk)
        jj = inum1
        weth_rate = hml_weth(k)
        !! adjustment weathering rate by temperature
        tadj = 1 / 293.15 - 1 / (t_av + 273.15)
        tadj = exp(1000 * tadj)
        weth_rate = weth_rate * tadj
        
        if (hml_rock(k,j) > 0) then	
            
            !! only if mine rocks exist in hru
            weth_amt(j) = hml_fr(k,j) * hml_rock(k,j) * weth_rate
            !weth_amt(j) = hml_fr(k,j) * hml_rock(k,j) * 0.95**(curyr-1) * weth_rate
            sol_hml_lab(k,j,1) = sol_hml_lab(k,j,1) + f * weth_amt(j)
            sol_hml_nlab(k,j,1) = sol_hml_nlab(k,j,1) + (1-f) * weth_amt(j)
            
            !! check by zhou 
            if (sol_hml_lab(k,j,1) < weth_amt(j) ) then
                zlf = 0
            end if
            
            !! summary calculations
            if (curyr > nyskip) then
                weth_hml(kk) = weth_hml(kk) + weth_amt(j)
            end if
            
            !! update sequence number for h. metal input for the jth hru
            nhml(j) = nhml(j) + 1
            
        end if
    end do
    
    return
    end