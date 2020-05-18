subroutine hml2EQspecies(X1,X2,X,SAT,Kp,rho,d)
	
	! This subroutine calculates the equilibrium amounts of h. metals given total amount,
	! partition constants, and other parameters

	!~~~~~~~~~~~~~~~~~~~~~~Input arguments~~~~~~~~~~~~~~~~~
	! X	(kg)			| total amount of h.metal in soil
	! SAT (mm)			| Saturation water content in soil
	! Kp (m3/kg)		| Partition constant between free cation and sorbed  labile species
	! KL (m3/g)**gamma	| Partition constant between ligand-bound species and free cation
	! rho	(ton/m3)	| Bulk density of soil
	! d	(mm)			| Layer depth
	! L	(mg/L;g/m3)  	| Ligand organic concentration in aquesou phase
	! gamma				| Ligand bounding number
	! See "表层土壤重金属的在固液两相间的平衡20120610.docx"	/"关于三元金属形态的平衡"
	
	!~~~~~~~~~~~~~~~~~~~~~~Output arguments~~~~~~~~~~~~~~~~
	! X1				| Amount as aqueous cation		
	! X2				| Amount as aqueous ligand bound metal
	! X3				| Amount as solid-sorbed  labile species

	!~~~~~~~~~~~~~~~~~~~~~~Local variables~~~~~~~~~~~~~~~~~
	! y1
	! y2
	! y3
	real X1,X2,X
	real SAT, Kp,rho,d
	real y1, y2, y3
	real z1, z2, z3

	z1=1.
	z2=(rho*1000.*d*Kp)/SAT
	
	Ztot = z1 + z2
	X1 = z1 / Ztot * X
	X2 = z2 / Ztot * X
	
	
	return
	end
