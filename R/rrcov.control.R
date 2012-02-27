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

rrcov.control <-
    function(alpha = 1/2, nsamp = 500, nmini = 300,
             seed = NULL, tolSolve = 1e-14,
	     trace = FALSE, use.correction = TRUE, adjust = FALSE)
{
    list(alpha = alpha, nsamp = nsamp, nmini = nmini, seed = as.integer(seed),
	 tolSolve = tolSolve,
	 trace = trace, use.correction = use.correction, adjust = adjust)
}
