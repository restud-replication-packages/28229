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

global localdir "$baseline_path"
global outdir   "$baseline_path${sep}output"
global auxdir   "$baseline_path${sep}auxiliar"

cd "$localdir"


//----------------------------------------------------
//---merge 2010 and 2016
//----------------------------------------------------

use "rfam16.dta", clear

keep nquest y yl yl1 yl2 yc yca yt yta ym ycf

rename y   yx
rename yl  ylx
rename yl1 yl1x
rename yl2 yl2x

gen year = 2016


save "$auxdir${sep}datax_2016.dta", replace

clear

use "rfam10.dta", clear

keep nquest y yl yl1 yl2 yc yca yt yta ym ycf

rename y   yx
rename yl  ylx
rename yl1 yl1x
rename yl2 yl2x

gen year = 2010

save "$auxdir${sep}datax_2010.dta", replace

clear

use "mpc_2010_2016.dta", clear

// keep if year == 2016

merge 1:1 nquest year using "$auxdir${sep}datax_2016.dta"
keep if _merge == 1 | _merge == 3
drop _merge

merge 1:1 nquest year using "$auxdir${sep}datax_2010.dta"
drop if _merge == 2
drop _merge

replace yx  = yx/1000
replace ylx = ylx/1000
replace yc  = yc/1000
replace yca = yca/1000
replace yt  = yt/1000
replace yta = yta/1000
replace ym  = ym/1000
replace ycf = ycf/1000

//----------------------------------------------------
//---computations
//----------------------------------------------------

gen ylabor = ylx + ym + yt 
gen ytotal = ylabor + ycf + yca


di "corr mpc and labor income"
cor ylabor mpc 


di "corr mpc and total income"
cor ytotal mpc 

