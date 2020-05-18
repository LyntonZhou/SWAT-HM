    real function hmlkd(ph,soc,sly)
    
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
    
    if (hmlname(1) == '        Zn') then
        
        hmlkd = 10**(-2.48  +0.69 * ph + 0.67 * LOG10(soc))
        
    else if (hmlname(1) == '        Cd') then
        
        hmlkd = 10**(-1.7  +0.62 * ph + 0.61 * LOG10(soc))
    
    else if (hmlname(1) == '        Pb') then
        
        hmlkd = 10**(1.32  +0.40 * ph + 0.50 * LOG10(soc))    
        
    end if

    return

    end function