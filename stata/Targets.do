
***************************************************************************
/* -------------------------------------------------------------------------- */
/* ------------------- PSID - Targets --------------------------------------- */
/* -------------------------------------------------------------------------- */

clear
clear matrix
set more off, perm
version 13
//cap log close
//log using results, replace text
* ssc install ginidesc
* ssc install sumdist

***************************************************************************
/*--- Set globals ---*/

global DataDir "C:\Users\Alessandro Di Nola\Dropbox\1 - RESEARCH PROJECTS\PROJECTMARIA\SHARED\codes and data\M files\rev_v16_v2_replication_package\stata\data" //*need to customize this directory
global OutputDir "C:\Users\Alessandro Di Nola\Dropbox\1 - RESEARCH PROJECTS\PROJECTMARIA\SHARED\codes and data\M files\rev_v16_v2_replication_package\stata"    // *need to customize this directory


global userestr = "1"         //1-use data cleaning (age and female HH), any other num-dont use
global hourswork_res = "260"  //restriction for annual hours worked
global incometrim = "0"       // restriction for min income
global includeunempl = "0"    // 1-include unemployed to workers, 0-exclude unemployed from the sample

/*---Choose definition of entrepreneurs----*/

global DefEntre = "2" //1-Business Ownes, 2-Self-Employed, 3-Both (Active business owners), 4-BO or SE

qui{
use "$DataDir\PSID_Amerged_new.dta", clear
gen income2 = redpregovinc 

/*-------------Restriction; Data cleaning -------------------*/

if $userestr==1 { 
*Use data cleaning 
drop if hsex==2               // keep only families where the HEAD is male 
keep if hage > 24 & hage < 65 // keep individuals of working age
}


/*---Choose definition of Entrepreneurs----*/


if $DefEntre==1 & $includeunempl==1 {

/*----First definition--Business Owners----*/
gen occhead=.
replace occhead=0 if entre1 == 5 
replace occhead=1 if entre1 == 1
drop if occhead==.
}


if $DefEntre==1 & $includeunempl==0 {
drop if hannhrs < $hourswork_res //first drop unemployed based on hours worked
drop if entre2==0     //drop those 'not working for money' - unemployed

/*----First definition--Business Owners----*/
gen occhead=.
replace occhead=0 if entre1 == 5 
replace occhead=1 if entre1 == 1
drop if occhead==.
}



/*----Second Definition--Self-Employed---*/

if $DefEntre==2 & $includeunempl==1{
gen occhead=.
replace occhead=0 if entre2 == 1 | entre2 == 2 | entre2 == 0  
replace occhead=1 if entre2 == 3
drop if occhead==.
}


/*----Second Definition--Self-Employed---*/

if $DefEntre==2 & $includeunempl==0{
drop if hannhrs < $hourswork_res //first drop unemployed
drop if entre2==0     //drop those 'not working for money' - unemployed

gen occhead=.
replace occhead=0 if entre2 ==1 | entre2 ==2  
replace occhead=1 if entre2 ==3 
drop if occhead==.
}



/*----Third Definition--Both-------*/

if $DefEntre==3 & $includeunempl==1{
gen occhead=.
replace occhead=1 if entre1==1 & entre2==3 
replace occhead=0 if entre1==5 | entre2 ==1 | entre2 ==2 | entre2 ==0  
drop if occhead==.
}

if $DefEntre==3 & $includeunempl==0{
drop if hannhrs < $hourswork_res //first drop unemployed
drop if entre2==0     //drop those 'not working for money' - unemployed

gen occhead=.
replace occhead=1 if entre1 ==1 & entre2==3  // 
replace occhead=0 if entre1==5 | entre2 ==1 | entre2 ==2
drop if occhead==.
}


/*-------Forth Definition--Union (OR)------*/
if $DefEntre==4 & $includeunempl==1{
gen occhead=0
replace occhead=1 if entre1==1 | entre2==3 

}

if $DefEntre==4 & $includeunempl==0{
drop if hannhrs < $hourswork_res //first drop unemployed
drop if entre2==0     //drop those 'not working for money' - unemployed

gen occhead=0
replace occhead=1 if entre1 ==1 | entre2==3  // 

}

gen s006 = 1        //PSID is a representative sample (no weights)
gen taxes = fiitax  //federal taxes
drop if taxes < 0



/*--------------------TARGETS-----------------------*/

/*-----Fraction of entrepreneurs ----*/
 
 tab occhead, matcell(x)
 matrix list x
 gen numwork=x[1,1]
 gen numentr=x[2,1]
 gen fracentr=numentr/(numwork+numentr)
 
 display "fraction of self-employed = "  fracentr
 gen frac_of_entre = fracentr*100
 gen T1=frac_of_entre //Target 1


summarize income2 if occhead==1
return list
gen sumincse = r(sum)
gen meanincse = r(mean)

summarize income2 if occhead==0
return list
gen sumincw = r(sum)
gen meanincw = r(mean)



/*------- Share of entrepreneurs' income----*/
gen shareentrinc = sumincse/ (sumincse+sumincw) 
display "share of entrepreneurial income = "  shareentrinc
gen share_entr_inc=shareentrinc*100
gen T2=share_entr_inc


/*------Ratio of mean income (entr to worker)-----*/
gen ratiomeaninc = meanincse/meanincw
display "ratio of mean income = "  ratiomeaninc

/*------Mean Income for the whole population----*/
sum income2 
gen meaninc=r(mean)

/*----Summary statistics of income----*/
sum income2, detail
sum income2 if occhead==0, detail
sum income2 if occhead==1, detail
gen income_w=income2 if occhead==0


/*------Income dist entre-----*/
gen income_e=income2 if occhead==1	

sum income_e, detail
return list
gen meaninc_e   = r(mean)
gen medianinc_e = r(p50)   
gen meanmedinc_e = meaninc_e/medianinc_e
}


/*---Distribution of Self-Employed Income, Table 5---*/
display "mean to median income entre = "  meanmedinc_e
	
/*-----GINI coefficient for gross INCOME, for workers and self-employed---*/
ginidesc income2, by (occhead)  m(a1) gk(a2)


qui{
/*----INCOME ENTRE-----*/

sumdist income_e, ngp(100)
return list
matrix BB = r(cush40)*100 \ (r(cush100)-r(cush80))*100 \ (r(cush100)-r(cush90))*100 // Bottom 40, 20, 10
matrix list BB

sumdist income_e, ngp(100)
return list
matrix BC = r(sh100)*100  // Top 1
matrix list BC

*----income entre dist-----*/
matrix BCC = BB \ BC
matrix rownames BCC =  Bottom_40\%  Top_20\%  Top_10\%  Top_1\%
matrix colnames BCC =  income_e
}
matrix list BCC

qui{
/*--------WEALTH-------*/

sum wealth
gen totalwealth = r(sum)
gen numwlth = r(N) 

/*----zero or negative wealth---*/
sum wealth if wealth<=0
gen numnegwlth = r(N)
gen frac_neg_wlth = (numnegwlth/numwlth)*100

/*recode negative wealth to zero wealth (if in the model no neg wealth)*/
replace wealth=0 if wealth < 0

sum wealth if occhead==1
gen totalwealthe=r(sum)

sum wealth if occhead==0
gen totalwealthw=r(sum)

/*---Summary statistics of wealh----*/
sum wealth
sum wealth, detail

gen wealth_w=wealth if occhead==0
gen wealth_e=wealth if occhead==1

/*---summary statistics of wealth for workers---*/
sum wealth if occhead==0, detail
return list
gen sumwealth_w = r(sum)
gen meanwlth_w  = r(mean)
gen medianwlth_w = r(p50)

/*---summary statistics of wealth for entrepreneurs---*/
sum wealth if occhead==1, detail
return list
gen sumwealth_e = r(sum)
gen meanwlth_e = r(mean)
gen medianwlth_e = r(p50)   //median wealth of entrepr.
gen meanmed_e = meanwlth_e/medianwlth_e
display "mean to median assets entre = "  meanmed_e
gen medwealthratio = medianwlth_e/medianwlth_w
display "median wealth ratio = "   medwealthratio

/*--the ratio of mean/median wealth of entrepreneurs to workers---*/
gen ratiomeanwlth = meanwlth_e/meanwlth_w
display "ratio of mean wealth entrepreneurs to workers = "  ratiomeanwlth

gen ratiomedianwlth=medianwlth_e/medianwlth_w
gen T4 = ratiomedianwlth
display "ratio of median wealth entrepreneurs to workers = "  ratiomedianwlth


/*---- assets owned by entrepreneurs----*/
gen assetsownedentr = totalwealthe/totalwealth
display "Assets owned by entrepreneurs = "  assetsownedentr
gen assets_owned_entr = assetsownedentr*100
gen T3 =assets_owned_entr
}
/*----- Wealth Distribution, Table 6 ----*/
/*----- Gini coefficient for wealth, for workers and self-employed------*/
ginidesc wealth, by (occhead) m(a3) gk(a4)

/*----Aggregate distribution statistics---*/
qui{
sum wealth
gen wlthmean = r(mean)
sum wealth, detail
gen wlthmed = r(p50)
gen meanmedwlth = wlthmean/wlthmed
gen T7 = meanmedwlth
}
disp "mean to median ratio, wealth = "  meanmedwlth

qui{
sum income2
gen incmean = r(mean)
sum income2, detail
gen incmed = r(p50)
gen meanmedinc = incmean/incmed
gen T11 = meanmedinc
disp "mean to median ratio, income = "  meanmedinc

sumdist wealth, 
return list
matrix A = (r(sh1)+r(sh2)+r(sh3)+r(sh4))*100 \ (r(sh9)+r(sh10))*100 \ r(sh10)*100 //bottom 40, top 20, top 10
matrix list A

sumdist income2, ngp(100)
return list
matrix B = r(cush40)*100 \ (r(cush100)-r(cush80))*100 \ (r(cush100)-r(cush90))*100 
matrix list B

sumdist wealth_e,
return list
matrix D =(r(sh1)+r(sh2)+r(sh3)+r(sh4))*100 \ (r(sh9)+r(sh10))*100 \ r(sh10)*100 //bottom 40, top 20, top 10
matrix list D


/*--- wealth and income distribution---*/
matrix C = A,B,D
matrix rownames C =  Bottom_40\%  Top_20\%  Top_10\% 
matrix colnames C = wealth income wealth_e 
matrix list C

sumdist wealth, ngp(100)   //top 1 percent
return list
matrix A1 = r(sh100)*100
matrix list A1

sumdist income2, ngp(100)
return list
matrix B1 = r(sh100)*100
matrix list B1

sumdist wealth_e, ngp(100)
return list
matrix D1 = r(sh100)*100
matrix list D1

/*---top 1 percent ----*/

matrix C1 = A1,B1,D1
matrix rownames C1 = Top_1\% 
matrix colnames C1 = wealth income wealth_e 
matrix list C1


/*--- wealth distribution ----*/
matrix CC1 = A \ A1
matrix rownames CC1 = Bottom_40\%  Top_20\%  Top_10\% Top_1\%
matrix colnames CC1 = wealth 
}
matrix list CC1


qui{
/*----wealth-income ratios----*/
gen wlth_inc_se = sumwealth_e/sumincse 
gen wlth_inc_w  = sumwealth_w/sumincw    

gen wlth_inc_se1 = wealth_e/income_e
gen wlth_inc_w1 = wealth_w/income_w

sum wlth_inc_se1, detail
gen meanwise= r(mean)
gen medianwise=r(p50)

sum wlth_inc_w1, detail
gen meanwiw=r(mean)
gen medianwiw=r(p50)
}



/*------ Summary Statistics: Table 12 -----*/
disp "fraction of entrepreneurs = "  frac_of_entre
disp "share of entre income     = "  share_entr_inc
disp "assets owned by entre     = "  assets_owned_entr
disp "ratio of median wealth    = "  ratiomedianwlth


/*------ Exit Rates, Table 13 -----*/
do exitrate.do 



//log close



