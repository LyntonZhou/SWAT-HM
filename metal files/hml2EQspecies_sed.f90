subroutine hml2EQspecies_sed(y1,y2,Kp,Por,rho)
	
	! This subroutine calculates the equilibrium amounts of h. metals given total amount,
	! partition constants, and other parameters

	!~~~~~~~~~~~~~~~~~~~~~~Input arguments~~~~~~~~~~~~~~~~~
	! Por    			| porosity of bed sediemnt
	! Kp (m3/kg)		| Partition constant between free cation and sorbed  labile species
	! rho	(ton/m3)	| Particle density of soil
	! d	(mm)			| Layer depth
	


	!~~~~~~~~~~~~~~~~~~~~~~Local variables~~~~~~~~~~~~~~~~~
	! y1
	! y2
    
	real Por, Kp,rho,d
	real y1, y2 
    real z1, z2, ztot

	z1=Por
    z2=(1-por)*rho*Kp
	!z2=(rho*d*Kp)/Por*1000
    ztot = z1 + z2
    y1 = z1/ztot
    y2 = z2/ztot
	
	return
	end