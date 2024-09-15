clear
set maxvar 10000
use "C:\Users\USER\Desktop\NDHS DATA 2018\NG_2018_DHS_02172020_1656_137126\NGIR7ADT\NGIR7AFL.DTA"

keep v013 v106 v716 v190 v025 v024 v201 v151 v225 v131 v021 v364 v130 v190 v384a v384b v384c v384d v005 v732 v312 v218 v228 v730 v171a sstate

save "data1.dta", replace

use "data1.dta"

recode v364 (4=0 "No") (3=1 "Yes"), gen(intention)
ta intention, mis
drop if intention==2

gen weight=v005/1000000

recode v130 (1 2=1 "Christian") (3=2 "Islam") (4 96=3 "Others"), gen(religion)

recode v190 (1 2=1 "pro-poor") (3/5=2 "pro-rich"), gen(w_index)

gen exp= v384a+v384b+v384c+v384d

recode exp (0=0 "No") (1/max=1 "Yes"), gen(fp_msg)

recode v312 (0=0 "Not Currently Using") (1/18 =1 "Currently Using"), gen(contrap_use)

recode v218 (0/5=1 "0-5") (6/10=2 "6-10") (11/15=2 "11-15"), gen(num_chdrnlivng)

recode v730 (14/35=1 "14-35") (36/55=2 "36-55") (56/75=3 "56-75") (76/95=4 "76-95"), gen(partner_age)

recode v201 (0/5=1 "0-5")(6/10=2 "6-10")(10/17=3 "10-17"), gen(childborn)

rename (v013 v106 v732 v025 v024 v151 v225 v384a v384b v384c v384d v228 v171a) (age edu employ resid region hsex cpregwant fpmsg_radio fpmsg_tv fpmsg_paper fpmsg_sms abort_preg use_web)

order intention age partner_age edu employ childborn cpregwant religion num_chdrnlivng use_web contrap_use fpmsg_radio fpmsg_tv fpmsg_paper fpmsg_sms abort_preg w_index resid region sstate



save "data1.dta", replace

********************************************************************************************************************************************************************************************
*  The decomposition code starts here
********************************************************************************************************************************************************************************************

clear
set maxvar 10000
use "data1.dta"


// Generates dummy variables
global D age partner_age edu employ childborn cpregwant religion num_chdrnlivng use_web contrap_use fpmsg_radio fpmsg_tv fpmsg_paper fpmsg_sms abort_preg w_index resid region
foreach d of global D {
tabulate `d', generate(`d')
}

global X "age2 age3 age4 age5 age6 age7 partner_age1 partner_age2 partner_age3 edu1 edu2 edu3 employ1 employ2 childborn1 childborn2 cpregwant1 cpregwant2 religion2 religion3 num_chdrnlivng1 use_web1 use_web2 fpmsg_radio1 fpmsg_tv2 fpmsg_paper2 fpmsg_sms1 abort_preg1 resid2 region1 region2 region3 region4 region6"

// Not thate collinear dummy variables were already dropped. These collinear variables were detected when the marginal probit regression was exploratorily done.

//To generate fractional rank variable
egen raw_rank=rank(v190), unique
sort raw_rank
quietly sum weight
gen wi=weight/r(sum)
gen cusum=sum(wi)
gen wj=cusum[_n-1]
replace wj=0 if wj==.
gen rank=wj+0.5*wi

qui sum intention
scalar mean=r(mean)
cor intention rank, cov
sca wi= (1-0)/((1-mean)*(mean-0))*2*r(cov_12)  //Wagstaff Index
sca list wi

dprobit intention $X [pw=weight]
matrix dfdx=e(dfdx) 

sca f=0
foreach x of global X {
qui {
mat b_`x' = dfdx[1,"`x'"]
sca b_`x' = _b[`x']
corr rank `x', cov
sca cov_`x' = r(cov_12)
sum `x' [aw=weight]
sca m_`x' = r(mean)
sca elas_`x' = (b_`x'*m_`x')/m_intention
sca CI_`x' = (1-0)/((1-m_`x')*(m_`x'-0))*2*cov_`x'
sca con_`x' = elas_`x'*CI_`x'
sca prcnt_`x' = con_`x'/wi
sca f=f+con_`x'
}

di "`x' elasticity:", elas_`x'
di "`x' concentration index:", CI_`x'
di "`x' contribution:", con_`x'
di "`x' percentage contribution:", prcnt_`x'
}

sca HI = wi - f
di "Horizontal Inequity Index:", HI