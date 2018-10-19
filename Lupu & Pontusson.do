* Noam Lupu and Jonas Pontusson
* The Structure of Inequality and the Politics of Redistribution
* Last revised March 2011

* Replication by Filippo, Julian, and Zara
clear
use "/Users/zarariaz/Dropbox/Fall 2018/PED/Replication/LupPon_APSR.dta"
set more 1
set varabbrev off
version 10
*ssc install std01


/*********************************************************************
1. Generate variables
**********************************************************************/

* Sort dataset
sort id year
tsset id year


* Invert disproportionality measure
replace disp_gall = -1 * disp_gall


* Rescale variables
replace fempar = fempar*100
replace union = union*100


* Standardize variables to [0,1]
egen stdpjoint = std01(pjoint)
egen stddisp_gall = std01(disp_gall)


* Interpolate missing values
by id: ipolate ratio9050 year, gen(pratio9050)
by id: ipolate ratio5010 year, gen(pratio5010)
by id: ipolate ratio9050s year, gen(pratio9050s)
by id: ipolate ratio5010s year, gen(pratio5010s)
by id: ipolate foreign year, gen(pforeign)
by id: ipolate voc year, gen(pvoc)


* Generate immigration measure (percent foreign)
replace pforeign = pforeign*1000
gen fpop = (pforeign/pop)*100


* Generate inequality measures
gen ratio9010 = pratio9050 * pratio5010
gen ratio9010s = pratio9050s * pratio5010s
gen skew = pratio9050/pratio5010
gen skews = pratio9050s/pratio5010s


* Generate averages for redistribution models
gen since=.
foreach n of numlist 1(1)10 {
	by id: replace since = year - year[_n-`n'] if redist!=. & redist[_n-`n']!=. & since==.
}
foreach x of varlist ratio9010 pratio9050 pratio5010 stdpjoint skew stddisp_gall pvoc union fpop fempar unempl turnout {
	foreach y of numlist 1(1)10{
		tssmooth ma m`y'`x' = `x', window (`y' 0)
	}
	gen dv`x'=.
	foreach y of numlist 2(1)9{
		replace dv`x' = l1.`x' if since==1
		replace dv`x' = m`y'`x' if since==`y' & m`y'`x'!=.
		replace dv`x' = m10`x' if (since==10 & m10`x'!=.) | (since==. & redist!=. & m10`x'!=.)
	}
	foreach y of numlist 1(1)10{
		drop m`y'`x'
	}
}


* Generate five-year moving averages for social-spending models
foreach x of varlist socspend pratio9050s pratio5010s ratio9010s skews dreher pop65 stdpjoint stddisp_gall fempar unempl union turnout pvoc fpop {
	gen ma`x' = (l1.`x' + l2.`x' + l3.`x' + l4.`x' + l5.`x')/5
}



/*********************************************************************
2. Analysis
**********************************************************************/

** Redistribution models (Table 2)

preserve
keep if redist!=.
sort id year
by id: egen order = seq()
tsset id order

xtpcse redist l1.redist dvpratio9050 dvpratio5010 dvturnout dvfempar dvstddisp_gall dvpvoc dvunion dvunempl, pairwise cor(ar1)
predict pred if e(sample), xb
gen resid = redist-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse redist l1.redist dvpratio9050 dvpratio5010 dvturnout dvfempar dvstddisp_gall dvpvoc dvunion dvunempl if outlier!=1, pairwise cor(ar1) hetonly
drop pred resid stresid outlier

xi: xtpcse redist dvpratio9050 dvpratio5010 i.id, pairwise cor(ar1)
predict pred if e(sample), xb
gen resid = redist -pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse redist dvpratio9050 dvpratio5010 i.id if outlier!=1, pairwise cor(ar1)
drop pred resid stresid outlier

xtpcse redist l1.redist dvratio9010 dvskew dvturnout dvfempar dvstddisp_gall dvpvoc dvunion dvunempl, pairwise cor(ar1)
predict pred if e(sample), xb
gen resid = redist -pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse redist l1.redist dvratio9010 dvskew dvturnout dvfempar dvstddisp_gall dvpvoc dvunion dvunempl if outlier!=1, pairwise cor(ar1)
drop pred resid stresid outlier

xi: xtpcse redist dvratio9010 dvskew i.id, pairwise cor(ar1)
predict pred if e(sample), xb
gen resid = redist -pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse redist dvratio9010 dvskew i.id if outlier!=1, pairwise cor(ar1)
drop pred resid stresid outlier

restore


** Social spending models (Table 3)

xtpcse socspend l1.socspend mapratio9050s mapratio5010s mapop65 mafempar maturnout mastddisp_gall mapvoc maunion maunempl madreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = socspend-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse socspend l1.socspend mapratio9050s mapratio5010s mapop65 mafempar maturnout mastddisp_gall mapvoc maunion maunempl madreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xi: xtpcse socspend mapratio9050s mapratio5010s gdpgrowth i.id, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = socspend-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse socspend mapratio9050s mapratio5010s gdpgrowth i.id if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xtpcse socspend l1.socspend maratio9010s maskews mapop65 mafempar maturnout mastddisp_gall mapvoc maunion maunempl madreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = socspend-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse socspend l1.socspend maratio9010s maskews mapop65 mafempar maturnout mastddisp_gall mapvoc maunion maunempl madreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xi: xtpcse socspend maratio9010s maskews gdpgrowth i.id, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = socspend-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse socspend maratio9010s maskews gdpgrowth i.id if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier


** Immigration models (Table 4)

preserve
keep if redist!=.
sort id year
by id: egen order = seq()
tsset id order

xtpcse redist l1.redist dvratio9010 dvskew dvfpop dvfempar dvturnout dvstddisp_gall dvpvoc dvunion dvunempl, pairwise cor(ar1)
predict pred if e(sample), xb
gen resid = redist -pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse redist l1.redist dvratio9010 dvskew dvfpop dvfempar dvturnout dvstddisp_gall dvpvoc dvunion dvunempl if outlier!=1, pairwise cor(ar1)
drop pred resid stresid outlier

restore

xtpcse socspend l1.socspend maratio9010s maskews mafpop mapop65 mafempar maturnout mastddisp_gall mapvoc maunion maunempl madreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = socspend-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse socspend l1.socspend maratio9010s maskews mafpop mapop65 mafempar maturnout mastddisp_gall mapvoc maunion maunempl madreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier


** Partisanship (Table 5)

xtpcse stdpjoint maskews mastddisp_gall maturnout if year>1979, pairwise

xtpcse stdpjoint maskews mastddisp_gall maturnout madreher if year>1979, pairwise

xtpcse stdpjoint maskews mastddisp_gall maturnout madreher mafpop if year>1979, pairwise




preserve
drop if year<1980
collapse (mean) stdpjoint skew stddisp_gall dreher fpop pop65 turnout fempar, by(id)

reg stdpjoint skew stddisp_gall turnout, robust

reg stdpjoint skew stddisp_gall turnout dreher, robust

reg stdpjoint skew stddisp_gall turnout fpop dreher, robust

restore


** Redistribution and social spending with partisanship (Table 6)

preserve
keep if redist!=.
sort id year
by id: egen order = seq()
tsset id order

xtpcse redist l1.redist dvstdpjoint dvratio9010 dvskew dvfpop dvfempar dvturnout dvstddisp_gall dvpvoc dvunion dvunempl, pairwise cor(ar1)
predict pred if e(sample), xb
gen resid = redist -pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse redist l1.redist dvstdpjoint dvratio9010 dvskew dvfpop dvfempar dvturnout dvstddisp_gall dvpvoc dvunion dvunempl if outlier!=1, pairwise cor(ar1)
drop pred resid stresid outlier

restore

xtpcse socspend l1.socspend mastdpjoint maratio9010s maskews mafpop mapop65 mafempar maturnout mastddisp_gall mapvoc maunion maunempl madreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = socspend-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse socspend l1.socspend mastdpjoint maratio9010s maskews mafpop mapop65 mafempar maturnout mastddisp_gall mapvoc maunion maunempl madreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

***EXTENSION***

clear all
import delimited /Users/zarariaz/Documents/GitHub/Lupu---Pontusson/merged.csv


destring gdpgrowth socspend ma_socspend ma_pratio9050s ma_pratio9010s ma_skews ma_dreher ma_pop65 ma_stdpjoint ma_stddisp_gall ma_fempar ma_unempl ma_union ma_turnout ma_pvoc ma_fpop, replace force
destring ma_pratio5010s  fambenefits_lag incapacity_lag pubspendinglbr_lag pubunemplag, replace force


*we recreate table 3 results using 4 disaggregated spending DVs


** Social spending models (Table 3)

sort id year
by id: egen order = seq()
tsset id order

**Family Benefits**

xtpcse fambenefits fambenefits_lag ma_pratio9050s ma_pratio5010s ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = fambenefits-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse fambenefits fambenefits_lag ma_pratio9050s ma_pratio5010s ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xi: xtpcse fambenefits ma_pratio9050s ma_pratio5010s gdpgrowth i.id, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = fambenefits-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse fambenefits ma_pratio9050s ma_pratio5010s gdpgrowth i.id if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xtpcse fambenefits fambenefits_lag ma_ratio9010s ma_skews ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = fambenefits-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse fambenefits fambenefits_lag ma_ratio9010s ma_skews ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xi: xtpcse fambenefits fambenefits_lag ma_ratio9010s ma_skews gdpgrowth i.id, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = fambenefits-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse socspend ma_ratio9010s ma_skews gdpgrowth i.id if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

**Incapacity 

xtpcse incapacity incapacity_lag ma_pratio9050s ma_pratio5010s ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = incapacity-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse incapacity incapacity_lag mapratio9050s mapratio5010s mapop65 mafempar maturnout mastddisp_gall mapvoc maunion maunempl madreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xi: xtpcse incapacity ma_pratio9050s ma_pratio5010s gdpgrowth i.id, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = incapacity-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse incapacity ma_pratio9050s ma_pratio5010s gdpgrowth i.id if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xtpcse incapacity incapacity_lag  ma_ratio9010s ma_skews ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = incapacity-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse incapacity incapacity_lag ma_ratio9010s ma_skews ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xi: xtpcse incapacity ma_ratio9010s ma_skews gdpgrowth i.id, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = incapacity-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse incapacity ma_ratio9010s ma_skews gdpgrowth i.id if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

**Public Spending on Labor

xtpcse pubspending_labor pubspendinglbr_lag ma_pratio9050s ma_pratio5010s ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = pubspending_labor-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse pubspending_labor pubspendinglbr_lag ma_pratio9050s ma_pratio5010s ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xi: xtpcse pubspending_labor ma_pratio9050s ma_pratio5010s gdpgrowth i.id, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = pubspending_labor-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse pubspending_labor ma_pratio9050s ma_pratio5010s gdpgrowth i.id if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xtpcse pubspending_labor pubspendinglbr_lag  ma_ratio9010s ma_skews ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = pubspending_labor-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse pubspending_labor pubspendinglbr_lag maratio9010s maskews ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xi: xtpcse pubspending_labor ma_ratio9010s ma_skews gdpgrowth i.id, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = pubspending_labor-pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse pubspending_labor ma_ratio9010s ma_skews gdpgrowth i.id if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

**Unemployment


xtpcse public_unemp pubunemplag ma_pratio9050s ma_pratio5010s ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = public_unemp -pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse public_unemp  pubunemplag mapratio9050s mapratio5010s ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xi: xtpcse public_unemp  ma_pratio9050s ma_pratio5010s gdpgrowth i.id, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = public_unemp -pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse public_unemp  mapratio9050s mapratio5010s gdpgrowth i.id if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xtpcse public_unemp pubunemplag ma_ratio9010s ma_skews ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = public_unemp -pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xtpcse public_unemp pubunemplag ma_ratio9010s ma_skews ma_pop65 ma_fempar ma_turnout ma_stddisp_gall ma_pvoc ma_union ma_unempl ma_dreher gdpgrowth if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier

xi: xtpcse public_unemp  ma_ratio9010s ma_skews gdpgrowth i.id, pairwise cor(psar1)
predict pred if e(sample), xb
gen resid = public_unemp -pred
egen stresid=std(resid)
gen outlier = 0 if e(sample)
replace outlier = 1 if abs(stresid)>1.5
xi: xtpcse public_unemp maratio9010s maskews gdpgrowth i.id if outlier!=1, pairwise cor(psar1)
drop pred resid stresid outlier


