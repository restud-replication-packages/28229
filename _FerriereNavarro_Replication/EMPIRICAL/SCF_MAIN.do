********************
*Housekeeping
********************

clear matrix
clear
set more off
set scheme s2color
set more off


//----------------------------------------------------
//---settings
//----------------------------------------------------

global sep = "/"  // sep = "/" for mac, and sep = "\" for windows

*--set baseline directory
global baseline_path "/Users/gastonnavarro/Dropbox/Axelle-Gaston/GovernmentSpending/REPLICATION/EMPIRICAL"

global localdir  "${baseline_path}"
global outputdir "${baseline_path}${sep}output"

cd "$localdir"

use "scf83b.dta", clear

//----------------------------------------------------
//---variables definition
//----------------------------------------------------
// see here: https://www.federalreserve.gov/econres/files/1983_codebk83.txt

rename b4503 age_head

//---weights
gen wgt_all = b3005
gen wgt_FRB = b3016

//---variables
*--assets
// b3303 = paper assets
// b3801 = gross value of other properies
// b3902 = value of vehicles

*--liabilities
// b3319 = total consumer debt

gen assets = b3303 + b3801 + b3902
gen debt   = b3319
gen nw     = assets - debt

gen nwt = b3324  // SCF net worh variable, which further includes housing and SS payments

drop if age_head < 18
drop if age_head > 65

//----------------------------------------------------
//---computation
//----------------------------------------------------

//---weight 
gen ww = wgt_FRB

drop if ww == 0

gen nwvar = nw

//---empirical cdf
global NQ = 5


sort nwvar
egen wsum = sum(ww)
gen  wpdf = ww/wsum
gen  wcdf = sum(wpdf) 

replace wcdf = 1 if wcdf > 1

drop wsum

//---assign quintile
gen quint = 189			// decile id

forvalues qq = 1(1)$NQ {
	di "qq" `qq'
	replace quint = `qq' if ( wcdf > (`qq'-1)/$NQ ) & ( wcdf <= `qq'/$NQ )
}
gen ones = 1


gen nwxww = nwvar*ww

collapse (sum) nwxww,by(quint)

egen nwsum = sum(nwxww)
gen  nwsh  = nwxww/nwsum
 
