	 subroutine zero_hml

	!!    ~ ~ ~ PURPOSE ~ ~ ~
	!!This subroutine zeros values for array variables involving heavy metals
	use parm

	chhml_conc=0.	    !(mch)
	chhml_kd=0.         !(mch)
	chhml_mix=0.        !(mch)
	chhml_rea=0.        !(mch)
	chhml_rsp=0.        !(mch)
	chhml_stl=0.        !(mch)	
	orig_sedhmlconc=0.  !(mch)	
	sedhml_act=0.       !(mch)
	sedhml_bry=0.       !(mch)
	sedhml_conc=0.      !(mch)
	sedhml_rea=0.       !(mch)
		
	sub_solhml=0.       !(msub)
	sub_sorhml=0.       !(msub)

	sub_hml=0.        !(mhml)
	
	orig_solhml=0.    !(mhml)
	sol_hmlkp=0.      !(mhml)
	sol_hml=0.        !(mhml)
	
	lkhml_conc=0.	    !(mres)
	lkhml_kd=0.         !(mres)
	lkhml_mix=0.        !(mres)
	lkhml_rea=0.        !(mres)
	lkhml_rsp=0.        !(mres)
	lkhml_stl=0.        !(mres)
	lkshml_act=0.       !(mres)
	lkshml_bry=0.       !(mres)
	lkshml_conc=0.      !(mres)
      lksedhml_conc=0.    !(mres)
	lkshml_rea=0.       !(mres)
	
	orig_lkhmlconc=0.     !(mres)
	orig_lksedhmlconc=0.    !(mres)	
	lkhml_mass=0.         !(mres)
	lkshml_mass=0.        !(mres)
		
	ap_ef=0.          !(mhmldb)
	decay_f=0.        !(mhmldb)
	decay_s=0.        !(mhmldb)
	hlife_f=0.        !(mhmldb)
	hlife_s=0.        !(mhmldb)
	pst_wsol=0.       !(mhmldb)
	skoc=0.           !(mhmldb)
	
	hml_dep=0.		!(mnr)
	
	nhml=0.			!(mhru)
      hruhml= 0	!(mhru)
      hruhmla = 0.	!(mhml,4,mhru)
      hruhmlm = 0.	!(mhml,4,mhru)
      hruhmly = 0.	!(mhml,4,mhru)
		
	hml_enr=0.        !(mhml)
	hml_sed=0.        !(mhml)
	hml_surq=0.       !(mhml)
	hml_zdb=0.        !(mhml)
	!hml_ph=0.         !(mhml)
      hml_agr_frac=0.  !(mhml)
      hml_agr_frac=0.  !(mhml)
	hml_lat=0.        !(mhml)
	n_hml_no=0.       !(mhml)
	hml_prk=0.        !(mhml)
	wshd_hmlap=0.     !(mhml)
	wshd_hmldg=0.     !(mhml)

	solhmlcnst=0.      !(mrecc)
	srbhmlcnst=0.      !(mrecc)
	
	solhmlyr=0.        !(mrecy)
	srbhmlyr=0.        !(mrecy)
	
	solhmlmon=0.       !(mrecm)
	srbhmlmon=0.       !(mrecm)
	
	hsolhml=0.        !(24)
	hsorhml=0.        !(24)

	whmlaao=0.        !(mhml)
	whmlmono=0.       !(mhml)
	whmlyro=0.        !(mhml)
	whmldayo=0.       !(mhml)
	
	weth_hml=0.	  !(mhml)
	wshd_hmldg=0.	  !(mhml)
      atmo_hml=0.       !(mhml)
	
	hmlname=""
	n_hml_mx=1	

	resushml =0.0
	setlhml= 0.0
	hml_zdb= 0.

	i_hml_no = 0
	irthml= 0
	percohml =0.
      
      hml_eqn = 0
	sol_hmlkp = 0.
	sol_hmlkx = 0.
	sol_hmlkm = 0.
	wtr_hmlkl = 0.
	wtr_hmlkr = 0.
	plt_hmlku = 0.
      hml_wsol = 0.
	hmlwashout_ef =0.
      hml_weth = 0.
      
      !! added by Zhou
      hml_dep_total = 0.
      hml_dep_frac = 0.
      hml_agr_total = 0.
      hml_agr_frac = 0.
      weth_amt = 0.
      plt_hml = 0.
      hml_plt = 0.

	hml_wsol_channel =0.

      orig_solhml_sol = 0.
	orig_solhml_lig = 0.
	orig_solhml_exch = 0.
	orig_solhml_nlab = 0.
      sol_hml_sol = 0.
	sol_hml_lig = 0.
	sol_hml_lab = 0.
	sol_hml_nlab = 0.

	hmlroute = 0.
	hruhmld = 0.

	sub_hml_sol = 0.
	sub_hml_lig = 0.
	sub_hml_exch = 0.
	sub_hml_nlab = 0.
      sub_hml_dep = 0.
      sub_hml_weth = 0.
      sub_hml_agr = 0.

	hml_fr= 0.
	hml_rock=0.    !20150509
	
	point_gamma= 1.0
	frac_sol_cbn_soluble = 0.1
	frac_rch_cbn_soluble = 0.1
	frac_res_cbn_soluble = 0.1

	hml_sol_o = 0.
	hml_lig_o = 0.
	hml_exch_o = 0.
	hml_nlab_o = 0.

	hmlrchdy = 0.
      hmlrchmono = 0.

	hml_lag = 0.

	return
	end