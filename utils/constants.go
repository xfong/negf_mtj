/* Main code for NEGF based electron transport simulator for magnetic tunnel junction */

package utils

import (
    "math"
)

var (
    Pi		= float64(4.0)*math.Atan(1.0);
    k_B		= float64(1.3806488e-23);
    hplanck	= float64(6.62606957e-34);
    hbar	= hplanck/2.0/Pi;
    echarge	= float64(1.60217657e-19);
    mu0		= Pi*4.0e-7;
    muB		= float64(9.27400968e-24);
    zplus	= 1e-9;
    m0		= float64(9.10938291e-31);
    k_q		= k_B / echarge;
);
