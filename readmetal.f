      subroutine readmetal

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine reads parameters from the heavy metal database 
!!    (metal.dat)

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    mhml        |none          |maximum number of heavy metals in 
!!                               |database
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    rainmob_ef(:)    |none     |Rain moblization efficiency (0-1)
!!    sol_hmlkp(:,:)   |(mg/kg)/(kg/m3)  |Soil partition coef. between solid and aqueous phase
!!    sol_hmlkx(:)     |per day          |Rate constant for conversion of  labile species to non-labile in soil
!!    sol_hmlkm(:)     |per day          |Rate constant for conversion of non-labile species to  labile in soil
!!    wtr_hmlkl(:)     |per day          |Rate constant for conversion of aqueous ions to metal-ligands in water
!!    wtr_hmlkr(:)     |per day          |Rate constant for conversion of aqueous metal-ligands to ions in water
!!    plt_hmlku(:)     |per day          |Rate constant for uptake of aqueous metal by plants
!!    hml_wsol(:)      |kg/m3 (1000ppm)  |solubility of heavy metal in water
!!                     |k/L			   |Different from the unit of organics etc.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    eof         |none          |end of file flag
!!    ip          |none          |counter which represents the array
!!                               |storage number of the metal data
!!                               |the array storage number is used by the
!!                               |model to access data for a specific
!!                               |metal
!!    ipnum       |none          |number of metal
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: ip, ipnum, eof, eqn
      character (len=10) :: nm
      real :: kx, km, kl, kr, ku, ksol, kwash, kweth
      real, dimension (3) :: kp 
      eof = 0


      do
        !!initialize local variables
        nm = ""
        kp = 0.0
        kx = 0.0
	  km = 0.0
	  kl = 0.0
	  kr = 0.0
	  ku = 0.0
	  ksol = 0.0
	  kwash = 0.0
        kweth = 0.0
	  gamma = 1.0

      !   read (34,5000,iostat=eof) ip, nm, eqn, kp, kx, 
      !&   km, kl, kr, ku, ksol, kwash, kweth, gamma
        read (34,5000,iostat=eof) ip, nm, eqn, kp, kx, 
     &   km, ksol, kweth, ku, gamma, kl, kr, kwash
        if (eof < 0) exit
        
        if (ip == 0) exit

        hmlname(ip) = nm
        hml_eqn  = eqn
        sol_hmlkp(ip,:) = kp
        sol_hmlkx(ip) = kx
        sol_hmlkm(ip) = km
        hml_wsol(ip) = ksol
        hml_weth(ip) = kweth
        plt_hmlku(ip) = ku
        point_gamma(ip)=gamma
        wtr_hmlkl(ip) = kl
        wtr_hmlkr(ip) = kr
        hmlwashout_ef(ip)= kwash
        
      !! set values for heavy metal routed through main channel network
        if (ip == irtpest) then
          hml_wsol_channel = hml_wsol(ip) * 1000.
        end if

      end do

      close (34)

      return
 !5000 format (i3,a10,i3,3f10.3,5e10.2,f10.3,e10.2,e10.2,f10.3)
 !5000format (i3,a10,i3,3f10.3,2e10.2,f10.3,3e10.2,f10.3,2e10.2)
 5000 format (i3,a10,i3,3f10.3,2f10.6,f10.3,3f10.6,f10.3,2f10.6)
      end
