##  rrcov : Scalable Robust Estimators with High Breakdown Point
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##

## "FIXME": If you would change this, you must "sync" with
##          1) covMcd()'s         default in ./covMcd.R
##          2) ltsReg.default()'s default in ./ltsReg.R
rrcov.control <-
    function(alpha = 1/2, nsamp = 500, nmini = 300,
             seed = NULL, tolSolve = 1e-14,
	     trace = FALSE, wgtFUN = "01.original",
             use.correction = identical(wgtFUN, "01.original"),
             adjust = FALSE)
{
    list(alpha=alpha, nsamp=nsamp, nmini=nmini, seed = as.integer(seed),
	 tolSolve=tolSolve, trace=trace, wgtFUN=wgtFUN,
	 use.correction=use.correction, adjust=adjust)
}

## Only for back compatibility, as some new args did not exist pre 2013-04,
## and callers of covMcd() may use a "too small"  'control' list:
getDefCtrl <- function(nm) {
    callerEnv <- parent.frame()
    if(is.null(get(nm, envir = callerEnv)))
	assign(nm, rrcov.control()[[nm]], envir=callerEnv)
}
