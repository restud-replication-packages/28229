
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
local yvar "gdp"

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
*---BP shocks
****************************************************************************************
local varregBP L(1/$LAG).lngdp  L(1/$LAG).lngov L(1/$LAG).mtr trend_* 

reg lngov `varregBP'  if yq>=$id & yq <= $fd

predict bp, residuals

preserve
keep date bp
outsheet using "$outputdir${sep}BP_shock.csv", comma noquote replace
restore


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

#delimit ;
twoway (rarea m_ub m_lb H_horiz) 
	   (line m zero_line H_horiz,
	    lc(dknavy black)
		lwidth(thick thin)
		ylabel(-0.5(0.5)1) ymtick(-0.5(0.5)1) xlabel( 0 (1) $Hmax ) ), // yaxis(-0.75(0.25)0.3) ),		  
		legend(off) // legend(order(1 "Progressive" 2 "Regressive") ring(0) position(8) rows(2) region(color(none))) 
		name("LINEAR", replace) ;		    	
#delimit cr


*--save
outsheet using "$outputdir${sep}LPM_FIGURE_9.csv", comma noquote replace
