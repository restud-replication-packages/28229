
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

*---unemp dummy
gen     indU = 0
replace indU = 1 if L.unemp >= $Utresh

gen indPS = ind     * indU
gen indPE = ind     * (1-indU)
gen indRS = (1-ind) * indU
gen indRE = (1-ind) * (1-indU)


forvalues ll = 1 (1) $LAG {
	// ln GDP
	gen lngdp_L`ll'_PS =    indPS * L`ll'.lngdp
	gen lngdp_L`ll'_PE =    indPE * L`ll'.lngdp
	gen lngdp_L`ll'_RS =    indRS * L`ll'.lngdp
	gen lngdp_L`ll'_RE =    indRE * L`ll'.lngdp		
	
	// ln GOV
	gen lngov_L`ll'_PS =    indPS * L`ll'.lngov
	gen lngov_L`ll'_PE =    indPE * L`ll'.lngov
	gen lngov_L`ll'_RS =    indRS * L`ll'.lngov
	gen lngov_L`ll'_RE =    indRE * L`ll'.lngov	
	
	// GOV
	gen gov_L`ll'_PS =    indPS * L`ll'.gov
	gen gov_L`ll'_PE =    indPE * L`ll'.gov
	gen gov_L`ll'_RS =    indRS * L`ll'.gov
	gen gov_L`ll'_RE =    indRE * L`ll'.gov	
	
	
	// NEWS
	gen news_L`ll'_PS =    indPS * L`ll'.news
	gen news_L`ll'_PE =    indPE * L`ll'.news
	gen news_L`ll'_RS =    indRS * L`ll'.news
	gen news_L`ll'_RE =    indRE * L`ll'.news	
	
	
	// MTR
	gen mtr_L`ll'_PS =    indPS * L`ll'.mtr
	gen mtr_L`ll'_PE =    indPE * L`ll'.mtr
	gen mtr_L`ll'_RS =    indRS * L`ll'.mtr
	gen mtr_L`ll'_RE =    indRE * L`ll'.mtr
	
	
	// ATR
	gen atr_L`ll'_PS =    indPS * L`ll'.atr
	gen atr_L`ll'_PE =    indPE * L`ll'.atr
	gen atr_L`ll'_RS =    indRS * L`ll'.atr
	gen atr_L`ll'_RE =    indRE * L`ll'.atr	
	
	// rdef
	gen rdef_L`ll'_PS =    indPS * L`ll'.rdef
	gen rdef_L`ll'_PE =    indPE * L`ll'.rdef
	gen rdef_L`ll'_RS =    indRS * L`ll'.rdef
	gen rdef_L`ll'_RE =    indRE * L`ll'.rdef	
	
	// rdefgdp
	gen rdefgdp_L`ll'_PS =    indPS * L`ll'.rdefgdp
	gen rdefgdp_L`ll'_PE =    indPE * L`ll'.rdefgdp
	gen rdefgdp_L`ll'_RS =    indRS * L`ll'.rdefgdp
	gen rdefgdp_L`ll'_RE =    indRE * L`ll'.rdefgdp
	
	
	// tb3
	gen tb3_L`ll'_PS =    indPS * L`ll'.tb3
	gen tb3_L`ll'_PE =    indPE * L`ll'.tb3
	gen tb3_L`ll'_RS =    indRS * L`ll'.tb3
	gen tb3_L`ll'_RE =    indRE * L`ll'.tb3
	
}

// ln GOV
gen lngov_PS =  indPS *lngov
gen lngov_PE =  indPE *lngov
gen lngov_RS =  indRS *lngov
gen lngov_RE =  indRE *lngov


// GOV
gen gov_PS =  indPS *gov
gen gov_PE =  indPE *gov
gen gov_RS =  indRS *gov
gen gov_RE =  indRE *gov


// NEWS
gen news_PS = indPS *news
gen news_PE = indPE *news
gen news_RS = indRS *news
gen news_RE = indRE *news


// constant
gen c_PS = indPS
gen c_PE = indPE
gen c_RS = indRS
gen c_RE = indRE


****************************************************************************************
*---BP shocks
****************************************************************************************
local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr trend_*  // L(1/$LAG).rdefgdp

reg lngov `varregBP'  if yq>=$idbp & yq <= $fdbp

predict bp, residuals

// bp
gen bp_PS = indPS *bp
gen bp_PE = indPE *bp
gen bp_RS = indRS *bp
gen bp_RE = indRE *bp


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
gen m_PS    = .
gen m_se_PS = .
gen m_lb_PS = .
gen m_ub_PS = .

gen m_PE    = .
gen m_se_PE = .
gen m_lb_PE = .
gen m_ub_PE = .

gen m_RS    = .
gen m_se_RS = .
gen m_lb_RS = .
gen m_ub_RS = .

gen m_RE    = .
gen m_se_RE = .
gen m_lb_RE = .
gen m_ub_RE = .

gen pp_val_S  = . 
gen pp_val_E  = . 

save "$auxdir${sep}coef_save.dta", replace 

restore

*---Loop for regressions
#delimit ;
gen dy = 0;
gen dg = 0; gen dg_PS = 0; gen dg_PE = 0; gen dg_RS = 0;  gen dg_RE = 0;
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
	
	replace dg_PS =  indPS *dg
	replace dg_PE =  indPE *dg
	replace dg_RS =  indRS *dg
	replace dg_RE =  indRE *dg	
	
	ivreg2 dy `varreg' (dg_PS dg_PE dg_RS dg_RE = bp_PS bp_PE bp_RS bp_RE news_PS news_PE news_RS news_RE) if yq>=$id & yq <= $fd, robust bw(auto) 
	
	test dg_PS - dg_RS = 0
	local ppS = r(p)
	
	test dg_PE - dg_RE = 0
	local ppE = r(p)	
	
	// store results //	
	
	preserve
	
	use "$auxdir${sep}coef_save.dta", clear
	
	local index = `h' + 1
	
	replace m_PS     = _b[dg_PS]  in `index'
	replace m_se_PS  = _se[dg_PS] in `index' 	
	replace m_PE     = _b[dg_PE]  in `index'
	replace m_se_PE  = _se[dg_PE] in `index' 
	
	replace m_RS     = _b[dg_RS]  in `index'
	replace m_se_RS  = _se[dg_RS] in `index' 
	replace m_RE     = _b[dg_RE]  in `index'
	replace m_se_RE  = _se[dg_RE] in `index' 
	
	replace pp_val_S = `ppS'      in `index' 
	replace pp_val_E = `ppE'      in `index' 
	
	save "$auxdir${sep}coef_save.dta", replace 
	
	restore	
}

****************************************************************************************
*---plot and save output
****************************************************************************************


use "$auxdir${sep}coef_save.dta", clear

*--plot

#delimit ;
replace m_lb_PS = m_PS - m_se_PS; replace m_ub_PS = m_PS + m_se_PS;
replace m_lb_PE = m_PE - m_se_PE; replace m_ub_PE = m_PE + m_se_PE;
replace m_lb_RS = m_RS - m_se_RS; replace m_ub_RS = m_RS + m_se_RS;
replace m_lb_RE = m_RE - m_se_RE; replace m_ub_RE = m_RE + m_se_RE;
#delimit cr

gen zero_line = 0
gen H_horiz   = _n - 1 
label variable H_horiz "Quarters"

#delimit ;
twoway (rarea m_ub_PS m_lb_PS H_horiz) (rarea m_ub_PE m_lb_PE H_horiz) (rarea m_ub_RS m_lb_RS H_horiz)  (rarea m_ub_RE m_lb_RE H_horiz)
	   (line m_PS m_PE m_RS  m_RE zero_line H_horiz,
	    lc(dknavy cranberry black)
		lwidth(thick thick thin)
		), // ylabel(-0.5(0.5)1) ymtick(-0.5(0.5)1) xlabel(0 (1) $Hmax ) ), // yaxis(-0.75(0.25)0.3) ),		  
		legend(order(1 "Progressive+Slack" 2 "Progressive+Expansion" 3 "Regressive+Slack" 4 "Regressive+Expansion") ring(0) position(8) rows(2) region(color(none))) 
		name("MACRO_unemp", replace) ;		    	
#delimit cr

*--save
keep m_P* m_se_P* m_R* m_se_R* pp_val_* H_horiz
keep if H_horiz == 4 | H_horiz == 8 | H_horiz == 12


outsheet using "$outputdir${sep}LPM_TABLE_11_unemp.csv", comma noquote replace
