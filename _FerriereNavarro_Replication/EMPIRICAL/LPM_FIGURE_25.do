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

*--State data
import delimited "$baseline_path${sep}DATA_STATE_FN.csv", delimiter(comma) varnames(1) 

****************************************************************************************
*---LPM set-up
****************************************************************************************
global Hmax = 16 	// Maximum horizon for the LPM  
global LAG  = 4 	// number of lags for the regression
global TP   = 2     // trend = t^[1, 2 , ... , TP]
global QLAG = 4     // number of lags for the Newey-West correction matrix

global id  = 1948    // final date to use in computations        
global fd  = 2006.75 // final date to use in computations        
// global fd  = 2012.00 // final date to use in computations        

global idbp  = $id    // final date to use in computations        
global fdbp  = $fd // final date to use in computations    

global LAGbp = $LAG
global TPbp  = $TP 

*--lhs variable
local yvar "state_taxrate"

****************************************************************************************
*---convert data to quarterly
****************************************************************************************

save "$auxdir${sep}DATA_STATE.dta", replace

*--add quarters
gen qq = 1
save "$auxdir${sep}temp.dta", replace

forvalues qqx = 2 (1) 4 {	
	use "$auxdir${sep}DATA_STATE.dta", clear
	gen qq = `qqx'
	
	merge 1:1 year id qq using "$auxdir${sep}temp.dta"
	
	drop _merge	
	save "$auxdir${sep}temp.dta", replace
}

sort id year qq
order id year qq name

gen quarter = year + (qq-1)/4

save "$auxdir${sep}temp.dta", replace

clear

****************************************************************************************
*---BP shocks
****************************************************************************************
import delimited "$baseline_path${sep}DATA_MACRO_FN.csv", delimiter(comma) varnames(1) 

*--rename quarter yq  
gen year    = int(quarter)
gen qq      = 1 + (quarter - year)*4
gen date    = yq(year,qq)
format date %tq
tsset date, quarterly 


gen yq = quarter
*** Time trend
gen trend = _n
forvalues tt = 1 (1) $TPbp {
	gen trend_`tt' = trend^`tt' 
}

gen lngdp = ln(gdp)
gen lngov = ln(gov)

local varregBP L(1/$LAGbp).lngdp  L(1/$LAGbp).lngov L(1/$LAGbp).mtr trend_* 

reg lngov `varregBP'  if yq>=$idbp & yq <= $fdbp

predict bp, residuals

order date

drop if quarter == .

drop yq

save "$auxdir${sep}DATA_MACRO_temp.dta", replace


*--merge macor data with panel data
use "$auxdir${sep}temp.dta", clear

merge m:1 quarter using "$auxdir${sep}DATA_MACRO_temp.dta" 

drop _merge

save "$auxdir${sep}DATA_MACRO_STATE.dta", replace


****************************************************************************************
*---set-up data
****************************************************************************************
drop date
gen date    = yq(year,qq)
format date %tq
xtset id date, quarterly 

sort id date
order date id name  income  state_inc_tax  local_inc_tax


****************************************************************************************
*---Variables computations
****************************************************************************************

***state tax rate
gen state_taxrate = state_inc_tax/income

***local tax rate
gen local_taxrate = local_inc_tax/income

***local tax rate
gen stateandlocal_taxrate = state_taxrate+local_taxrate

***state income
gen lnincome_state = ln(income)

***state income per capita
gen incpc    = income/pop
gen lnincpop = ln(incpc)


gen yq = quarter

gen inreg = 1
replace inreg = 0 if id == 0     // US
replace inreg = 0 if id == 2000  // Alaska
replace inreg = 0 if id == 15000 // Hawaii
// replace inreg = 0 if state_taxrate == 0

****************************************************************************************
*---rhs variables for regression
****************************************************************************************
local varreg L(1/$LAG).lngdp  L(1/$LAG).lngov  L(1/$LAG).mtr   L(1/$LAG).lnincpop  L(1/$LAG).`yvar' 

**************************
***** 	Regression	 *****
**************************

*** Declare space for coeffients

preserve
clear 
global Lmax = $Hmax + 1   // global Lmax = $Hmax + 1
set obs $Lmax
gen m    = 0
gen m_se = 0
gen m_lb = 0
gen m_ub = 0

gen pp_coef = .
gen pp_se   = .
gen pp_df   = .
gen pp_pv   = .
gen H_horiz   = _n - 1 

save "$auxdir${sep}coef_save.dta", replace 

restore

****************************************************************************************
*---regression
****************************************************************************************

#delimit ;
gen dy = 0;
gen dg = 0; 
gen resid = 0;
gen residhm1 = 0;
gen residhm2 = 0;
#delimit cr

pause on

local residreg resid

 
// foreach ih of local nlist {					 // forvalues ih = -$placebo (1) $Hmax {					

forvalues h = 0 (1) $Hmax {					
	disp "***** h = " `h' " ******"
	
	
	replace dy = F`h'.`yvar' - L.`yvar'	
	replace dg = F`h'.lngov  - L.lngov
	
	ivreg2 dy `varreg' i.id c.trend_*#id  (dg = bp news) if (yq>=$id & yq<= $fd) & inreg == 1,  dkraay(4)
		
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

replace m_lb = m - m_se
replace m_ub = m + m_se

gen rhog = 0.9^H_horiz //1  0.9^H_horiz  gen gmodel   = 0.90^H_horiz

replace m    = rhog*m
replace m_se = rhog*m_se
replace m_lb = rhog*m_lb
replace m_ub = rhog*m_ub


gen zero_line = 0
// gen H_horiz   = _n - 1 
label variable H_horiz "Quarters"


#delimit ;
twoway (rarea m_ub m_lb H_horiz) //  (line gmodel H_horiz)
	   (line m  zero_line H_horiz, // gmodel
	    lc(dknavy black)
		lwidth(thick thin)
		xlabel(0 (1) $Hmax ) ), // ylabel(-0.0005(0.0002)0.0007) ymtick(-0.0005(0.0002)0.0007) xlabel(0 (1) $Hmax ) ), // yaxis(-0.75(0.25)0.3) ),		  		
		title(`yvar')
		legend(off) // legend(order(1 "Progressive" 2 "Regressive") ring(0) position(8) rows(2) region(color(none))) 
		name(`yvar', replace) ;		    	 // name(`tvar', replace) ;		    	 // name("LINEAR `tvar'", replace) ;		    	
#delimit cr

*--save
outsheet using "$outputdir${sep}LPM_FIGURE_25.csv", comma noquote replace
