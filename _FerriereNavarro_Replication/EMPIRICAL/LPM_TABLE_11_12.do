
clear matrix
clear
set more off

****************************************************************************************
*---set up and directories
****************************************************************************************

global sep = "/"  // sep = "/" for mac, and sep = "\" for windows

*--set baseline directory
global baseline_path "/Users/gastonnavarro/Dropbox/Axelle-Gaston/GovernmentSpending/REPLICATION/EMPIRICAL"

global outputdir "$baseline_path${sep}output_tables"
global auxdir    "$baseline_path${sep}auxiliar"


****************************************************************************************
*---load data
****************************************************************************************

import delimited "$baseline_path${sep}DATA_MACRO_FN.csv", delimiter(comma) varnames(1) 


****************************************************************************************
*---LPM set-up
****************************************************************************************
global Hmax = 16 	// Maximum horizon for the LPM  
global TP   = 4     // trend = t^[1, 2 , ... , TP]
global QLAG = 4     // number of lags for the Newey-West correction matrix

*--lhs variable
local yvar "gdp"

*---can switch pvar: pfed (federal taxes only), and psoc (federal + social security taxes)
local pvar "pfed"

*-- n cases: 1) benchmark, 2) no MTR,  3) with ATR, 4) with deficit, 5) with tb3, 6) from 1953
*-           7) to 2015,   8) only bp, 9) only bp from 1953, 10) only rz, 11) only rz from 1953
*-           12) 4 lags,   13) 2 lags, 14) wa 8, 15) wa 16, 16) wb 4, 17) wb 12    
global case = 17
if $case == 1{  // *--benchmark
	*------------------------------------------------------------------------------	
	local cases "benchmark"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}
else if $case == 2{  // *--no MTR
	*------------------------------------------------------------------------------	
	local cases "no_mtr"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}
else if $case == 3{  // *--with ATR
	*------------------------------------------------------------------------------	
	local cases "w_atr"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R atr_L*_P atr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 4{  // *--with rdef
	*------------------------------------------------------------------------------	
	local cases "w_rdef"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	// local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R rdefgdp_L*_P rdefgdp_L*_R c_P c_R trend_*
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R rdef_L*_P rdef_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 5{  // *--with tb3
	*------------------------------------------------------------------------------	
	local cases "w_tb3"
	
	global id   = 1920    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R tb3_L*_P tb3_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 6{  // *--from 1953
	*------------------------------------------------------------------------------	
	local cases "from_1953"
	
	global id   = 1953    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 7{  // *--to 2015
	*------------------------------------------------------------------------------	
	local cases "to_2015"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2015.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 8{  // *--only bp
	*------------------------------------------------------------------------------	
	local cases "only_bp"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 9{  // *--only bp from 1953
	*------------------------------------------------------------------------------	
	local cases "only_bp_from_1953"
	
	global id   = 1953    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 10{  // *--only rz
	*------------------------------------------------------------------------------	
	local cases "only_rz"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 11{  // *--only rz from 1953
	*------------------------------------------------------------------------------	
	local cases "only_rz_from_1953"
	
	global id   = 1953    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 12{  // *--lag 4
	*------------------------------------------------------------------------------	
	local cases "lag_4"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 4 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 13{  // *--lag 2
	*------------------------------------------------------------------------------	
	local cases "lag_2"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 2 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 14{  // *--wa 8
	*------------------------------------------------------------------------------	
	local cases "wa_8"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 8    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 15{  // *--wa 16
	*------------------------------------------------------------------------------	
	local cases "wa_16"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 16    // 12 quarters forward for state indicator
	global w_dw = 8     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 16{  // *--wb 4
	*------------------------------------------------------------------------------	
	local cases "wb_4"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 4     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

else if $case == 17{  // *--wb 12
	*------------------------------------------------------------------------------	
	local cases "wb_12"
	
	global id   = 1913    // final date to use in computations        
	global fd   = 2006.75 // final date to use in computations 
	
	global LAG  = 8 	// number of lags for the regression
	
	global w_up = 12    // 12 quarters forward for state indicator
	global w_dw = 12     // 8 quarters backward for state indicator
	
	local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_* 
	
	local ivar bp_P bp_R news_P news_R
	
	local varreg lngdp_L*_P  lngdp_L*_R lngov_L*_P lngov_L*_R mtr_L*_P mtr_L*_R c_P c_R trend_*
	*------------------------------------------------------------------------------
}

global idbp  = $id    // final date to use in bp shock computations        
global fdbp  = $fd // final date to use in bp shock computations    

****************************************************************************************
*---variables computations
****************************************************************************************
rename quarter yq  
gen year    = int(yq)
gen quarter = 1 + (yq - year)*4
gen date    = yq(year,quarter)
format date %tq
tsset date, quarterly 

order date

*---Time trend
gen trend = _n
forvalues tt = 1 (1) $TP {
	gen trend_`tt' = trend^`tt' 
}


*---logs and normalization
gen lngdp   = ln(gdp)
gen lngov   = ln(gov) 
gen lninv   = ln(rinvfx) 
gen lnwgenf = ln(wgenofarm)
gen lnhours = ln(hours)

gen rdefgdp = rdef/L.gdp
 
 
****************************************************************************************
*---progressivity index and by-prog variables
****************************************************************************************

 
*---Compute indicator variable
sum trend
local t1 = r(min)
local t2 = r(max)
disp "t1 = `t1' and t2 = `t2'"
gen ind = .

global tt1 = `t1' + 1
global tt2 = `t2' 


forvalues tt = $tt1 (1) $tt2 { 
	disp "tt = `tt'"
	egen m_up = mean(`pvar') if ( trend >= `tt' & trend <= `tt' + $w_up )
	egen m_dw = mean(`pvar') if ( trend <  `tt' & trend >= `tt' - $w_dw )
	
	qui sum m_up 
	local mn_up = r(mean)
	qui sum m_dw 
	local mn_dw = r(mean)
	
	replace ind = 1 if `mn_up' >= `mn_dw' in `tt' 
	replace ind = 0 if `mn_up' <  `mn_dw' in `tt' 
	
	drop m_up m_dw
	
}



*---variables x state
forvalues ll = 1 (1) $LAG {
	// ln GDP
	gen lngdp_L`ll'_P =    ind * L`ll'.lngdp
	gen lngdp_L`ll'_R = (1-ind)* L`ll'.lngdp
	
	// ln GOV
	gen lngov_L`ll'_P =    ind * L`ll'.lngov
	gen lngov_L`ll'_R = (1-ind)* L`ll'.lngov
	
	// ln INV
	gen lninv_L`ll'_P =    ind * L`ll'.lninv
	gen lninv_L`ll'_R = (1-ind)* L`ll'.lninv
	
	// ln wgenf
	gen lnwgenf_L`ll'_P =    ind * L`ll'.lnwgenf
	gen lnwgenf_L`ll'_R = (1-ind)* L`ll'.lnwgenf
	
	// lnhours
	gen lnhours_L`ll'_P =    ind * L`ll'.lnhours
	gen lnhours_L`ll'_R = (1-ind)* L`ll'.lnhours
	
	// GOV
	gen gov_L`ll'_P =    ind * L`ll'.gov
	gen gov_L`ll'_R = (1-ind)* L`ll'.gov
	
	// NEWS
	gen news_L`ll'_P =    ind * L`ll'.news
	gen news_L`ll'_R = (1-ind)* L`ll'.news
	
	// Progressivity
	gen p_L`ll'_P =    ind * L`ll'.`pvar'
	gen p_L`ll'_R = (1-ind)* L`ll'.`pvar'	
	
	// MTR
	gen mtr_L`ll'_P =    ind * L`ll'.mtr
	gen mtr_L`ll'_R = (1-ind)* L`ll'.mtr
	
	// ATR
	gen atr_L`ll'_P =    ind * L`ll'.atr
	gen atr_L`ll'_R = (1-ind)* L`ll'.atr
	
	// rdef
	gen rdef_L`ll'_P =    ind * L`ll'.rdef
	gen rdef_L`ll'_R = (1-ind)* L`ll'.rdef
	
	// rdefgdp
	gen rdefgdp_L`ll'_P =    ind * L`ll'.rdefgdp
	gen rdefgdp_L`ll'_R = (1-ind)* L`ll'.rdefgdp
	
	
	// tb3
	gen tb3_L`ll'_P =    ind * L`ll'.tb3
	gen tb3_L`ll'_R = (1-ind)* L`ll'.tb3
}

// GOV
gen lngov_P =    ind *lngov
gen lngov_R = (1-ind)*lngov
// GOV
gen gov_P =    ind *gov
gen gov_R = (1-ind)*gov
// NEWS
gen news_P =    ind *news
gen news_R = (1-ind)*news
// constant
gen c_P = ind
gen c_R = 1-ind



****************************************************************************************
*---BP shocks
****************************************************************************************
// local varregBP declared above 

reg lngov `varregBP'  if yq>=$idbp & yq <= $fdbp

predict bp, residuals

// bp
gen bp_P =    ind *bp
gen bp_R = (1-ind)*bp

****************************************************************************************
*---rhs variables for regression
****************************************************************************************
 // defined above alreadt


****************************************************************************************
*---regression
****************************************************************************************

*---Declare space for coeffients
preserve
clear 
global Lmax = $Hmax + 1
set obs $Lmax
gen m_P    = .
gen m_se_P = .
gen m_lb_P = .
gen m_ub_P = .

gen m_R    = .
gen m_se_R = .
gen m_lb_R = .
gen m_ub_R = .

gen pp_coef = .
gen pp_se   = .
gen pp_df   = .
gen pp_val  = . 


save "$auxdir${sep}coef_save.dta", replace 

restore

*---Loop for regressions
#delimit ;
gen dy = 0;
gen dg = 0; gen dg_P = 0; gen dg_R = 0; 
gen resid = 0;
#delimit cr

pause on

local residreg resid

gen cumuly = 0
gen cumulg = 0

 
forvalues ih = 0 (1) $Hmax {					

	local h = `ih' 
	
	disp "h = "`h'
	
	gen f`h'cumuly = ((F`h'.`yvar'-L.`yvar')/L.gdp) + cumuly
	gen f`h'cumulg = ((F`h'.gov-L.gov)/L.gdp)       + cumulg
	
	replace cumuly = f`h'cumuly
	replace cumulg = f`h'cumulg
	
	replace dy = f`h'cumuly
	replace dg = f`h'cumulg
	
	
	replace dy = . if yq + (`h' / 4) > $fd
	replace dg = . if yq + (`h' / 4) > $fd
	
	replace dg_P =    ind *dg
	replace dg_R = (1-ind)*dg					
	
	ivreg2 dy `varreg' (dg_P dg_R = `ivar') if yq>=$id & yq <= $fd, robust bw(auto) 
	
	lincom dg_P - dg_R
	local cc = r(estimate) // _b[UUF`ff'_f] // r(estimate)
	local ss = r(se)	   // _se[UUF`ff'_f] // r(se)				
	local dg = r(df)	
	
	test dg_P - dg_R = 0
	local pp = r(p)
	
	// store results //	
	preserve
	
	use "$auxdir${sep}coef_save.dta", clear
	
	local index = `h' + 1
	
	replace m_P     = _b[dg_P]  in `index'
	replace m_se_P  = _se[dg_P] in `index' 
	
	replace m_R     = _b[dg_R]  in `index'
	replace m_se_R  = _se[dg_R] in `index' 
	
	replace pp_coef = `cc'      in `index' 
	replace pp_se   = `ss'      in `index' 
	replace pp_df   = `dg'      in `index' 
	replace pp_val  = `pp'      in `index'
	
	save "$auxdir${sep}coef_save.dta", replace 
	
	restore	
}

****************************************************************************************
*---plot and save output
****************************************************************************************


use "$auxdir${sep}coef_save.dta", clear


*--plot
replace m_lb_P = m_P - m_se_P
replace m_ub_P = m_P + m_se_P
replace m_lb_R = m_R - m_se_R
replace m_ub_R = m_R + m_se_R

gen zero_line = 0
gen H_horiz   = _n - 1 
label variable H_horiz "Quarters"

#delimit ;
twoway (rarea m_ub_P m_lb_P H_horiz) (rarea m_ub_R m_lb_R H_horiz)
	   (line m_P m_R zero_line H_horiz,
	    lc(dknavy cranberry black)
		lwidth(thick thick thin)
		), // ylabel(-0.5(0.5)1) ymtick(-0.5(0.5)1) xlabel(0 (1) $Hmax ) ), // yaxis(-0.75(0.25)0.3) ),		  
		legend(order(1 "Progressive" 2 "Regressive") ring(0) position(8) rows(2) region(color(none))) 
		name("MACRO", replace) ;		    	
#delimit cr

*--save

keep m_P m_se_P m_R m_se_R pp_val H_horiz
keep if H_horiz == 4 | H_horiz == 8 | H_horiz == 12

outsheet using "$outputdir${sep}TABLES_11_12_`cases'.csv", comma noquote replace
