
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

global id  = 1913    // final date to use in computations        
global fd  = 2006.75 // final date to use in computations        


*--lhs variable
// yvar should be: atr_psz (average rate), atr_t50_psz (top-50 average rate), atr_b50_psz (bottom-50 average rate), or datr_t50b50 (top-50 minus bottom-50 rate)
local yvar "atr_psz"  // atr_psz  atr_t50_psz  atr_b50_psz  datr_t50b50 

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


*---tax rates distribution from Piketty, Saez, Zucman (2019)
gen atr_t50_psz = ((1/0.5)*atr_psz) - atr_b50_psz
gen datr_t50b50 =  (atr_t50_psz) - (atr_b50_psz)
 

****************************************************************************************
*---BP shocks
****************************************************************************************
local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov L(1/$LAG).mtr trend_* 

reg lngov `varregBP'  if yq>=$id & yq <= $fd

predict bp, residuals

****************************************************************************************
*---rhs variables for regression
****************************************************************************************
#delimit ;
local varreg L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr  trend_*;  //local control "controlbench"; 
#delimit cr


****************************************************************************************
*---regression
****************************************************************************************

*---Declare space for coeffients
preserve
clear 
global Lmax = $Hmax + 1
set obs $Lmax
gen m    = .
gen m_se = .
gen m_lb = .
gen m_ub = .



save "$auxdir${sep}coef_save.dta", replace 

restore

*---Loop for regressions
#delimit ;
gen dy = 0;
gen dg = 0; 
gen resid = 0;
#delimit cr

pause on

local residreg resid
 
forvalues ih = 0 (1) $Hmax {					

	local h = `ih' 
	
	disp "h = "`h'
	
	replace dy = F`h'.`yvar' - L.`yvar'	
	replace dg = F`h'.lngov  - L.lngov			
	
	
	replace dy = . if yq + (`h' / 4) > $fd
	replace dg = . if yq + (`h' / 4) > $fd
	
	
	ivreg2  dy `varreg' (dg = bp news)  if yq>=$id & yq <= $fd,   robust bw(auto)	
		
	// store results //	
	preserve
	
	use "$auxdir${sep}coef_save.dta", clear
	
	local index = `h' + 1
	
	replace m     = _b[dg]  in `index'
	replace m_se  = _se[dg] in `index' 
	
	
	save "$auxdir${sep}coef_save.dta", replace 
	
	restore	
}

****************************************************************************************
*---plot and save output
****************************************************************************************


use "$auxdir${sep}coef_save.dta", clear

*--plot
replace m_lb = m - m_se
replace m_ub = m + m_se

gen zero_line = 0
gen H_horiz   = _n - 1 
label variable H_horiz "Quarters"


*--impute response to AR(1) shock as in model
gen rhog = 0.9^H_horiz
replace m    = rhog*m
replace m_se = rhog*m_se
replace m_lb = rhog*m_lb
replace m_ub = rhog*m_ub


#delimit ;
twoway (rarea m_ub m_lb H_horiz) 
	   (line m zero_line H_horiz,
	    lc(dknavy black)
		lwidth(thick thin)
		ylabel(-0.01(0.01)0.05) ymtick(-0.01(0.01)0.05) xlabel( 0 (1) $Hmax ) ), // yaxis(-0.75(0.25)0.3) ),		  
		legend(off) // legend(order(1 "Progressive" 2 "Regressive") ring(0) position(8) rows(2) region(color(none))) 
		name("LINEAR", replace) ;		    	
#delimit cr


*--save
outsheet using "$outputdir${sep}LPM_FIGURE_10_`yvar'.csv", comma noquote replace
