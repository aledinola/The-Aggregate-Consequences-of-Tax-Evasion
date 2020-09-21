qui{
keep	idperson year occhead						



/* ------------- transition rates ------------------------------------------- */


levelsof year, local(levels)
mat 	trans = J(10,6,.)							/* matrix for results */
mat	rownames trans = `levels'
mat 	colnames trans = EM nEM ERW SE nSE ERE


local j = 1989
while `j' <= 2000 {
	qui{
	preserve
	
	if	`j' == 1989 | `j' == 1990 | `j' == 1993| `j' == 1994| `j' == 1995{					
		drop if year != `j' & year != (`j' + 1)	
		}
	else {	
	if	`j' == 1991 | `j' == 1996  | `j' == 1998 | `j' == 2000  {
		drop if year != `j' & year != (`j' + 2)	
		}
	}
	
	by	idperson, sort:	gen nr = _n				
	by	idperson:		egen mnr = max(nr)
	drop if mnr == 1							/* delete if only 1 appearance */

	by	idperson:		   egen labstat = sum(occhead)			/* 0 = 2 * employee, 1 = employee and self-employed, 2 = 2 * self-employed */
	gen	EM = 1		if labstat == 0 & year == `j'			/* stayed employed */
	gen	SE = 1 		if labstat == 2 & year == `j'			/* stayed self-employed */
		
	count if EM == 1
	local	EMstay = r(N)							/* # stayed employed */
	count if occhead == 0 & year == `j'				
	local	EMall = r(N)							/* # employed */
		
	count if SE == 1							/* # stayed self-employed */
	local	SEstay = r(N)	
	count if occhead == 1 & year ==`j'					/* # self-employed */
	local	SEall = r(N)


	
	mat	trans[rownumb(trans,"`j'"),colnumb(trans,"EM")] = `EMstay' / `EMall'
	mat	trans[rownumb(trans,"`j'"),colnumb(trans,"nEM")] =  `EMall'
	mat	trans[rownumb(trans,"`j'"),colnumb(trans,"ERW")] =  1-`EMstay' / `EMall'
	mat	trans[rownumb(trans,"`j'"),colnumb(trans,"SE")] = `SEstay' / `SEall'
	mat	trans[rownumb(trans,"`j'"),colnumb(trans,"nSE")] = `SEall'
    mat	trans[rownumb(trans,"`j'"),colnumb(trans,"ERE")] = 1-`SEstay' / `SEall'

	
		
	if	`j' == 1989 | `j' == 1990 | `j' == 1993| `j' == 1994| `j' == 1995{					/* go to next year */
		local	j = `j' + 1	
		}
	else {
	if	`j' == 1991 | `j' == 1996  | `j' == 1998 | `j' == 2000 {
		local	j = `j' + 2	
		}
	
		}
	restore	
	}
	}
	
/* -------------------------------------------------------------------------- */


}
/* ------------ table ------------------------------------------------------- */
*frmttable using "$OutputDir/Exit_rates", statmat(trans) /*
**/		varlabels replace findcons coljust(lc) /*
**/		nocenter sfmt(fc) statfont(normalsize) tex sdec($decimals1) /*
**/		ctitles("Year","\% stayed employed","N Workers", /*
**/		"Exit Rate W","\% stayed Entre","N Entre" "Exit Rate E") 		
/* -------------------------------------------------------------------------- */



