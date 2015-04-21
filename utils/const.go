/* Main code for NEGF based electron transport simulator for magnetic tunnel junction */

package utils

import (
    "math"
)

var (
    Pi		    = float64(4.0 * math.Atan(1.0));
    K_B		    = float64(1.3806488e-23);
    Hplanck	    = float64(6.62606957e-34);
    Hbar	    = Hplanck/2.0/Pi;
    Echarge	    = float64(1.60217657e-19);
    Mu0		    = Pi*4.0e-7;
    MuB		    = float64(9.27400968e-24);
    Zplus	    = 1e-12;
    M0		    = float64(9.10938291e-31);
    K_q		    = K_B / Echarge;
    ECurrFactor = Echarge * Echarge / Hplanck;
    SCurrFactor = Echarge / (4*Pi);
)
