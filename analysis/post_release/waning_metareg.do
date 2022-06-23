*cd "/Users/eh1415/Documents/covid-change-ve-over-time/release20220505"

*ssc install metareg

log using waning_metareg, replace

set more off
clear


import delimited data_metareg.csv

replace estimate="" if estimate=="NA"
destring estimate, replace

replace confhigh="" if confhigh=="NA"
destring confhigh, replace

replace conflow="" if conflow=="NA"
destring conflow, replace

gen loghr=estimate
gen lnu=confhigh
gen lnl=conflow

gen se1=(lnu-loghr)/1.96
gen se2=(loghr-lnl)/1.96

corr se1 se2

gen seloghr=(se1+se2)/2
// lines 26, 27, and 31 are usually simply coded as: gen seloghr = (lnu - lnl)/(2 * 1.96)
drop se1 se2 lnu lnl

* create strata from subgroup
rename subgroup stratum 

* create outcome
rename outcome temp

gen outcome=1 if temp=="covidadmitted"
replace outcome=2 if temp=="coviddeath"
replace outcome=3 if temp=="postest"
replace outcome=4 if temp=="noncoviddeath"
replace outcome=5 if temp=="anytest"

label define outcome 1 "COVID-19 hospitalisation" 2 "COVID-19 death" 3 "Positive test" ///
 4 "Non-COVID death" 5  "Any test"
label values outcome outcome
drop temp

* create vaccine from comparison
encode comparison, gen(vaccine)
tab vaccine
label list vaccine

* create model variable
rename model temp
gen model=1 if temp=="unadjusted"
replace model=2 if temp=="part_adjusted"
replace model=3 if temp=="max_adjusted"

label define model 1 "unadjusted" 2 "part_adjusted" 3 "max_adjusted" 
label values model model
drop temp

* create prior infection variable
rename prior temp
gen prior=1 if temp=="FALSE"
replace prior=2 if temp=="TRUE"

label define prior 1 "FALSE" 2 "TRUE" 
label values prior prior
drop temp

replace k=k-1

sort model outcome stratum vaccine
save waning_metareg.dta, replace

* what is the point of lines 63-72?
// metareg loghr k if outcome==2 & stratum==2 & vaccine==1, wsse(seloghr)
// local a=_b[k]
// local b=_se[k]
// local c=_b[_cons]
// local d=_se[_cons]
//
// di `a'
// di `b'
// di `c'
// di `d'

tempname memhold
postfile `memhold' model outcome stratum vaccine prior logrhr selogrhr loghr1 seloghr1 using results, replace

forvalues m=1/3 {
	forvalues i=1/5 {
		forvalues v=1/3 {
			forvalues s=1/4 {
				forvalues p=1/2 {
					di "A: " `m' `i' `s' `v' `p'
					count if model==`m' & outcome==`i' & stratum==`s' & vaccine==`v' & prior==`p' & loghr<.
					if r(N)>2 {
					di "B: " `a' /*`m' `i' `s' `a' `v' `p'*/
					metareg loghr k if model==`m' & outcome==`i' & stratum==`s' & vaccine==`v' & prior==`p', wsse(seloghr)
					local a=_b[k]
					local b=_se[k]
					local c=_b[_cons]
					local d=_se[_cons]
					post `memhold' (`m') (`i') (`s') (`v') (`p') (`a') (`b') (`c') (`d')
					}	
				}		
			}				
		}
	}	
}

postclose `memhold'
di "`memhold'"

use results, clear
sort model outcome stratum vaccine
save results, replace

log close
