    real function hmlkd(ph,soc,sly)
    !! This function calcualtes Kd1 value using pH and SOC

    use parm
    implicit none
    !! ph	: soil pH
    !! soc	: soil organic carbon

    real :: ph, soc, zlf
    integer :: sly

    if (soc < 0.2) then
        zlf = 0.
    end if

    !if (sly > 1) then
    !    ph = ph + 0.2
    !end if
    if (hml_eqn == 2) then

    if (hmlname(1) == '        Zn') then
        hmlkd = 10**(-1.77  + 0.66 * ph + 0.79 * LOG10(soc))
    else if (hmlname(1) == '        Cd') then
        hmlkd = 10**(-1.04  + 0.55 * ph + 0.70 * LOG10(soc))
    else if (hmlname(1) == '        Pb') then
        hmlkd = 10**(1.32  + 0.40 * ph + 0.50 * LOG10(soc))
    else if (hmlname(1) == '        Cu') then
        hmlkd = 10**(0.45  + 0.34 * ph + 0.65 * LOG10(soc))   
    else if (hmlname(1) == '        Ni') then
        hmlkd = 10**(0.99  + 0.30 * ph)
    else if (hmlname(1) == '        Hg') then
        hmlkd = 1    
    end if

    else

    if (hmlname(1) == '        Zn') then
        hmlkd = 10**(-2.48  + 0.69 * ph + 0.67 * LOG10(soc))
    else if (hmlname(1) == '        Cd') then
        hmlkd = 10**(-1.7  + 0.62 * ph + 0.61 * LOG10(soc))
    else if (hmlname(1) == '        Ni') then
        hmlkd = 10**(-0.58  + 0.44 * ph)
    else if (hmlname(1) == '        Hg') then
        hmlkd = 1    
    end if

    end if

    return

    end function