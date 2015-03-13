/* Main code for NEGF based electron transport simulator for magnetic tunnel junction */

package main

import (
    "math"
    "math/cmplx"
    "fmt"
)

var (
    Pi		= float64(4.0)*math.Atan(1.0);
    hplanck	= float64(6.62606957e-34);
    hbar	= hplanck/2.0/Pi;
    echarge	= float64(1.60217657e-19);
    mu0		= Pi*4.0e-7;
    muB		= float64(9.27400968e-24);
    zplus	= complex128(cmplx.Sqrt(-1.0)*1e-9);
    m0		= float64(9.10938291e-31);
)

var (
    m_ox, m_fm_L, m_fm_R	float64;
    aSpace			float64;
    d_ox, d_fm_L, d_fm_R	float64;
);

func main() {
    fmt.Println(Pi)
    fmt.Println(hplanck)
    fmt.Println(hbar)
    fmt.Println(echarge)
    fmt.Println(mu0)
    fmt.Println(muB)
    fmt.Println(zplus)
    fmt.Println(m0)

    m_ox := float64(0.315);
    m_fm_L := float64(0.72);
    m_fm_R := m_fm_L;

    aSpace := float64(1.0e-11);
    d_ox := float64(1.15e-9);
    d_fm_L := float64(2.0e-11);
    d_fm_R := float64(2.0e-11);

    N_fm_L := int(d_fm_L/aSpace) - 1;
    N_fm_R := int(d_fm_R/aSpace) - 1;
    N_ox := int(d_ox/aSpace) - 1;

    t_base := hbar*hbar/2/echarge/m0/aSpace/aSpace;
    t_fm_L := t_base/m_fm_L;
    t_fm_R := t_base/m_fm_R;
    t_ox := t_base/m_ox;
}
