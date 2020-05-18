	subroutine hml2EQspecies_reach(x1,x2,Kp,SED)
	
	! This subroutine calculates the equilibrium fraction of h. metals given sediment concentration,
	!  soluble organic concentration,partition constants, and other parameters

	!~~~~~~~~~~~~~~~~~~~~~~Input arguments~~~~~~~~~~~~~~~~~
	! SED				| Sediment concentration in reach flow
	! Kp				| Partition constant between free cation and sorbed  labile species
	! See "表层土壤重金属的在固液两相间的平衡20120610.docx"	/"关于水中三元金属形态的平衡"
	
	!~~~~~~~~~~~~~~~~~~~~~~Output arguments~~~~~~~~~~~~~~~~
	! x1				| Fraction of H.metal as aqueous cation		
	! x2				| Fraction of H.metal as solid-sorbed  labile species

	!~~~~~~~~~~~~~~~~~~~~~~Local variables~~~~~~~~~~~~~~~~~
	! y1
	! y2
	! y3
	real x1,x2
	real SED, Kp
	real z1, z2

	z1=1.
	z2=Kp*SED
	
	Ztot = z1 + z2
	x1 = z1 / Ztot
	x2 = z2 / Ztot
	
	return
	end
