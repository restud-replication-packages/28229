
clear matrix
clear
set more off

****************************************************************************************
*---set up and directories
****************************************************************************************

global sep = "/"  // sep = "/" for mac, and sep = "\" for windows

*--set baseline directory
global baseline_path "/Users/gastonnavarro/Dropbox/Axelle-Gaston/GovernmentSpending/REPLICATION/EMPIRICAL"

global outputdir "$baseline_path${sep}output"
global auxdir    "$baseline_path${sep}auxiliar"


****************************************************************************************
*---load data
****************************************************************************************

import delimited "$baseline_path${sep}DATA_MACRO_FN.csv", delimiter(comma) varnames(1) 


****************************************************************************************
*---LPM set-up
****************************************************************************************
global Hmax = 16 	// Maximum horizon for the LPM  
global LAG  = 8 	// number of lags for the regression
global TP   = 4     // trend = t^[1, 2 , ... , TP]
global QLAG = 4     // number of lags for the Newey-West correction matrix

global w_up = 12    // 12 quarters forward for state indicator
global w_dw = 8     // 8 quarters backward for state indicator


global id  = 1913    // final date to use in computations        
global fd  = 2006.75 // final date to use in computations        

global idbp  = 1913    // final date to use in bp shock computations        
global fdbp  = 2006.75 // final date to use in bp shock computations    


*---unemp threshold for dummy
global Utresh = 6.5  // unemp state = 1 if UNEMP >= Utresh

*--lhs variable
local yvar "gdp"

*---can switch pvar: pfed (federal taxes only), and psoc (federal + social security taxes)
local pvar "pfed"
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
gen lngdp = ln(gdp)
gen lngov = ln(gov)
 

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

****************************************************************************************
*---BP shocks
****************************************************************************************
local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_*  // L(1/$LAG).rdefgdp

reg lngov `varregBP'  if yq>=$idbp & yq <= $fdbp

predict bp, residuals

gen isgn = 1
replace isgn = 0 if bp <0

****************************************************************************************
*---variables x state
****************************************************************************************

gen indPP = ind     * isgn
gen indPN = ind     * (1-isgn)
gen indRP = (1-ind) * isgn
gen indRN = (1-ind) * (1-isgn)


forvalues ll = 1 (1) $LAG {
	// ln GDP
	gen lngdp_L`ll'_PP =    indPP * L`ll'.lngdp
	gen lngdp_L`ll'_PN =    indPN * L`ll'.lngdp
	gen lngdp_L`ll'_RP =    indRP * L`ll'.lngdp
	gen lngdp_L`ll'_RN =    indRN * L`ll'.lngdp		
	
	// ln GOV
	gen lngov_L`ll'_PP =    indPP * L`ll'.lngov
	gen lngov_L`ll'_PN =    indPN * L`ll'.lngov
	gen lngov_L`ll'_RP =    indRP * L`ll'.lngov
	gen lngov_L`ll'_RN =    indRN * L`ll'.lngov		
	
	
	// NEWS
	gen news_L`ll'_PP =    indPP * L`ll'.news
	gen news_L`ll'_PN =    indPN * L`ll'.news
	gen news_L`ll'_RP =    indRP * L`ll'.news
	gen news_L`ll'_RN =    indRN * L`ll'.news	
	
	
	// MTR
	gen mtr_L`ll'_PP =    indPP * L`ll'.mtr
	gen mtr_L`ll'_PN =    indPN * L`ll'.mtr
	gen mtr_L`ll'_RP =    indRP * L`ll'.mtr
	gen mtr_L`ll'_RN =    indRN * L`ll'.mtr		
	
}

// ln GOV
gen lngov_PP =  indPP *lngov
gen lngov_PN =  indPN *lngov
gen lngov_RP =  indRP *lngov
gen lngov_RN =  indRN *lngov


// GOV
gen gov_PP =  indPP *gov
gen gov_PN =  indPN *gov
gen gov_RP =  indRP *gov
gen gov_RN =  indRN *gov


// NEWS
gen news_PP = indPP *news
gen news_PN = indPN *news
gen news_RP = indRP *news
gen news_RN = indRN *news

// bp
gen bp_PP = indPP *bp
gen bp_PN = indPN *bp
gen bp_RP = indRP *bp
gen bp_RN = indRN *bp



// constant
gen c_PP = indPP
gen c_PN = indPN
gen c_RP = indRP
gen c_RN = indRN

****************************************************************************************
*---rhs variables for regression
****************************************************************************************
#delimit ;
; local varreg lngdp_L*_P*  lngdp_L*_R* lngov_L*_P* lngov_L*_R* mtr_L*_P* mtr_L*_R* c_P* c_R* trend_*; //local control "controlbench"; 
#delimit cr


****************************************************************************************
*---regression
****************************************************************************************

*---Declare space for coeffients
preserve
clear 
global Lmax = $Hmax + 1
set obs $Lmax
gen m_PP    = .
gen m_se_PP = .
gen m_lb_PP = .
gen m_ub_PP = .

gen m_PN    = .
gen m_se_PN = .
gen m_lb_PN = .
gen m_ub_PN = .

gen m_RP    = .
gen m_se_RP = .
gen m_lb_RP = .
gen m_ub_RP = .

gen m_RN    = .
gen m_se_RN = .
gen m_lb_RN = .
gen m_ub_RN = .

gen pp_val_P  = . 
gen pp_val_N  = . 


save "$auxdir${sep}coef_save.dta", replace 

restore

*---Loop for regressions
#delimit ;
gen dy = 0;
gen dg = 0; gen dg_PP = 0; gen dg_PN = 0; gen dg_RP = 0;  gen dg_RN = 0;
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
	
	replace dg_PP =  indPP *dg
	replace dg_PN =  indPN *dg
	replace dg_RP =  indRP *dg
	replace dg_RN =  indRN *dg	
	
	ivreg2 dy `varreg' (dg_PP dg_PN dg_RP dg_RN = bp_PP bp_PN bp_RP bp_RN news_PP news_PN news_RP news_RN) if yq>=$id & yq <= $fd, robust bw(auto) 
	
	test dg_PP - dg_RP = 0
	local ppP = r(p)
	
	test dg_PN - dg_RN = 0
	local ppN = r(p)	
	
	// store results //	
	
	preserve
	
	use "$auxdir${sep}coef_save.dta", clear
	
	local index = `h' + 1
	
	replace m_PP     = _b[dg_PP]  in `index'
	replace m_se_PP  = _se[dg_PP] in `index' 	
	replace m_PN     = _b[dg_PN]  in `index'
	replace m_se_PN  = _se[dg_PN] in `index' 
	
	replace m_RP     = _b[dg_RP]  in `index'
	replace m_se_RP  = _se[dg_RP] in `index' 
	replace m_RN     = _b[dg_RN]  in `index'
	replace m_se_RN  = _se[dg_RN] in `index' 
	
	replace pp_val_P = `ppP'      in `index' 
	replace pp_val_N = `ppN'      in `index' 
	
	save "$auxdir${sep}coef_save.dta", replace 
	
	restore	
}

****************************************************************************************
*---plot and save output
****************************************************************************************


use "$auxdir${sep}coef_save.dta", clear

*--plot

#delimit ;
replace m_lb_PP = m_PP - m_se_PP; replace m_ub_PP = m_PP + m_se_PP;
replace m_lb_PN = m_PN - m_se_PN; replace m_ub_PN = m_PN + m_se_PN;
replace m_lb_RP = m_RP - m_se_RP; replace m_ub_RP = m_RP + m_se_RP;
replace m_lb_RN = m_RN - m_se_RN; replace m_ub_RN = m_RN + m_se_RN;
#delimit cr

gen zero_line = 0
gen H_horiz   = _n - 1 
label variable H_horiz "Quarters"

#delimit ;
twoway (rarea m_ub_PP m_lb_PP H_horiz) (rarea m_ub_PN m_lb_PN H_horiz) (rarea m_ub_RP m_lb_RP H_horiz)  (rarea m_ub_RN m_lb_RN H_horiz)
	   (line m_PP m_PN m_RP m_RN zero_line H_horiz,
	    lc(dknavy cranberry black)
		lwidth(thick thick thin)
		), // ylabel(-0.5(0.5)1) ymtick(-0.5(0.5)1) xlabel(0 (1) $Hmax ) ), // yaxis(-0.75(0.25)0.3) ),		  
		legend(order(1 "Progressive+Positive" 2 "Progressive+Negative" 3 "Regressive+Positive" 4 "Regressive+Negative") ring(0) position(8) rows(2) region(color(none))) 
		name("MACRO_sign", replace) ;		    	
#delimit cr

*--save
keep m_P* m_se_P* m_R* m_se_R* pp_val_* H_horiz
keep if H_horiz == 4 | H_horiz == 8 | H_horiz == 12


outsheet using "$outputdir${sep}LPM_TABLE_11_sign.csv", comma noquote replace
