    subroutine hmlaging
!! This function calcualtes the transformation between
!!  labile and non-labile amount of metals in soil
    !implicit none

    use parm

    real :: exchtonlab, nlabtoexch, t_av
    real :: kx, km, ph, Tadj, pHadj, zlf
    integer :: j, ly, k, kk, k_i

    j = 0
    j = ihru
    kk = 0
    k = 0
    k_i = 0
    kx = 0.
    km = 0.
    exchtonlab = 0.
    nlabtoexch = 0.
    t_av =tmpav(j)
    
    do ly = 1, sol_nly(j)
        do k_i = 1, n_hml_mx
            
            k = n_hml_no(k_i)	 ! ID No. of metal in db 
            kk = i_hml_no(k)   ! Sequence No. of input metal, in fact this sequence number is useles
            kx = sol_hmlkx(k)
            km = sol_hmlkm(k)
            
            !! debug 
            if (hml_rock(k,j) > 0) then
                zlf = 0
            end if
            if (j == 860) then
                zlf = 0
            end if
            
            ! Adjustment aging rate by pH for Zn
            ph = hml_ph(j)
            if (hmlname(kk) == '        Zn') then
              km = km + 4.18e-4 * (7- ph)
              kx = kx - 4.18e-4 * (7- ph)
            else if (hmlname(kk) == '        Cd') then
              km = km + 1.38e-4 * (7- ph)
              kx = kx - 1.38e-4 * (7- ph)
            else if (hmlname(kk) == '        Pb') then
              km = km
              kx = kx
            end if
            
            !! Adjustment aging rate by Temperature
            Tadj = 1 / 293.15 - 1 / (t_av + 273.15)
            Tadj = exp(1000 * Tadj)
            kx = kx * Tadj
            km = km * Tadj
            
            
            exchtonlab = sol_hml_lab(k,j,ly) * kx
            nlabtoexch = sol_hml_nlab(k,j,ly) * km
            
            !! debug by Zhou
            if (abs(nlabtoexch - exchtonlab) > 10 ) then
                zlf = 0
            end if
            
            !sol_hml_nlab(k,j,ly) = sol_hml_nlab(k,j,ly) - nlabtoexch
            !sol_hml_lab(k,j,ly) = sol_hml_lab(k,j,ly) + nlabtoexch
            
            sol_hml_nlab(k,j,ly) = sol_hml_nlab(k,j,ly) - nlabtoexch + exchtonlab
            sol_hml_lab(k,j,ly) = sol_hml_lab(k,j,ly) + nlabtoexch - exchtonlab
            
        end do
    end do

    return
    end 