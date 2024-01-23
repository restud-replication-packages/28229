
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
global QLAG = 4     // number of lags for the Newey-West correction matrix

// full sample: id=1913 & fd=2006.75 // post accord: id=1951 & fd=2006.75 // post 80s: id=1980 & fd=2006.75 
global case = 3 // case = 1 for full sample 1913Q1-2006Q4, case = 2 for post-accord sample 1951Q1-2006Q4, case = 3 for post-80s sample 1980Q1-2006Q4
if $case == 1{  // *--full sample 
	*------------------------------------------------------------------------------	
	global id  = 1913    // final date to use in computations        
	global fd  = 2006.75 // final date to use in computations        
	local ids = "1913"    
	local fds = "2006"	
	
	global TP   = 4      // trend = t^[1, 2 , ... , TP]
	*------------------------------------------------------------------------------	
}
if $case == 2{  // *--post-accord sample 
	*------------------------------------------------------------------------------	
	global id  = 1951    // final date to use in computations        
	global fd  = 2006.75 // final date to use in computations        
	local ids = "1951"    
	local fds = "2006"
	
	global TP   = 4      // trend = t^[1, 2 , ... , TP]
	*------------------------------------------------------------------------------	
}

if $case == 3{  // *--post-80s sample 
	*------------------------------------------------------------------------------	
	global id  = 1980    // final date to use in computations        
	global fd  = 2006.75 // final date to use in computations        
	local ids = "1980"    
	local fds = "2006"
	
	global TP   = 2      // trend = t^[1, 2 , ... , TP]
	*------------------------------------------------------------------------------	
}

global idbp  = $id    // initial date to use for bp computations        
global fdbp  = $fd    // final date to use for bp computations    

global idgg  = 1913    // initial date to use dg normalization
global fdgg  = 2006.75 // final date to use dg normalization
	
*--lhs variable
local yvar "tb3"


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

reg lngov `varregBP'  if yq>=$idbp & yq <= $fdbp

predict bp, residuals


****************************************************************************************
*---rhs variables for regression
****************************************************************************************
#delimit ;
local varreg L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr L(1/$LAG).tb3 trend_*;  //local control "controlbench"; 
#delimit cr


****************************************************************************************
*---normalize so that dg = 1% in full selcted case
****************************************************************************************

gen dyx = lngov - L.lngov
gen dgx = (gov-L.gov)/L.gdp

replace dyx = . if yq > $fdgg
replace dgx = . if yq > $fdgg	
	
ivreg2 dyx L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr  trend_* (dgx = bp news) if yq>=$idgg & yq <= $fdgg, robust bw(auto)   // robust bw(auto)

global g_adj = 0.01/_b[dgx]
//-------------------------------------------------------------------------------------------------------------------------//

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
	replace dg = (gov-L.gov)/L.gdp 
	
	
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

*--normalize st dg = 1% in full selcted sample case
#delimit ;
replace m = m*$g_adj ; replace m_se = m_se*$g_adj ;
#delimit cr


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
		ylabel(-0.01(0.01)0.05) ymtick(-0.01(0.01)0.05) xlabel( 0 (1) $Hmax ) ), // yaxis(-0.75(0.25)0.3) ),		  
		legend(off) // legend(order(1 "Progressive" 2 "Regressive") ring(0) position(8) rows(2) region(color(none))) 
		name("LINEAR", replace) ;		    	
#delimit cr


*--save
outsheet using "$outputdir${sep}LPM_FIGURE_23_`ids'_`fds'.csv", comma noquote replace
