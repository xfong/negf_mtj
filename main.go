/* Main code for NEGF based electron transport simulator for magnetic tunnel junction */

package main

import (
    "math"
    "math/cmplx"
    "fmt"
    "github.com/negf_mtj/negf_mtj/cmplxSparse"
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
    Ub				float64;
    E_split_L, E_split_R	float64;
    Vmtj			float64;
    th_F,ph_F			float64;
    th_H,ph_H			float64;
);

func main() {
    fmt.Println("Pi =", Pi)
    fmt.Println("Planck constant =", hplanck)
    fmt.Println("Reduced Planck constant =", hbar)
    fmt.Println("Elementary charge =", echarge)
    fmt.Println("Permeability of free space =", mu0)
    fmt.Println("Bohr magneton =", muB)
    fmt.Println("Small imaginary number =", zplus)
    fmt.Println("Free electron mass =", m0)

    // Material parameters
    m_ox = float64(0.315);
    m_fm_L = float64(0.72);
    m_fm_R = m_fm_L;

    Ub = 0.79;
    E_split_L = 0.930;
    E_split_R = 0.930;

    fmt.Println("----------------------------------------");
    fmt.Printf("Ub = %g, E_split_L = %g, E_split_R = %g\n", Ub, E_split_L, E_split_R);

    // Configuration parameters
    th_F, ph_F = Pi/2.0, 0.0;    
    th_H, ph_H = 0.0, 0.0;    

    fmt.Println("----------------------------------------");
    fmt.Println("Left FM:");
    fmt.Printf("th = %g, ph = %g\n", th_F, ph_F);

    fmt.Println("----------------------------------------");
    fmt.Println("Right FM:")
    fmt.Printf("th = %g, ph = %g\n", th_H, ph_H);
    fmt.Println("----------------------------------------");

/*
    aSpace := float64(1.0e-11);
    d_ox := float64(1.15e-9);
    d_fm_L := float64(2.0e-11);
    d_fm_R := float64(2.0e-11);
*/
    aSpace := float64(0.25e-9);
    d_ox := float64(1.0e-9);
    d_fm_L := float64(0.5e-9);
    d_fm_R := float64(0.5e-9);

    N_fm_L := int(d_fm_L/aSpace);
    N_fm_R := int(d_fm_R/aSpace);
    N_ox := int(d_ox/aSpace) - 1;
    fmt.Printf("N_ox = %d, N_fm_L = %d, N_fm_R = %d\n",N_ox, N_fm_L, N_fm_R);
    fmt.Println("----------------------------------------");

    t_base := hbar*hbar/2/echarge/m0/aSpace/aSpace;
    t_fm_L := t_base/m_fm_L;
    t_fm_R := t_base/m_fm_R;
    t_ox := t_base/m_ox;


    fmt.Printf("t_ox = %f, t_fm_L = %f, t_fm_R = %f\n",t_ox, t_fm_L, t_fm_R);
    fmt.Printf("2t_ox = %f, 2t_fm_L = %f, 2t_fm_R = %f\n",2*t_ox, 2*t_fm_L, 2*t_fm_R);
    fmt.Println("----------------------------------------");

    // Initialize a base matrix template
    Hamiltonian := cmplxSparse.New();
    cmplxSparse.MakeHamTriDiag(N_fm_L+N_fm_R+N_ox+2,Hamiltonian);

    // Build the main diagonal first
    cmplxSparse.ScaleRangeSparseMatrixIP(0, 2*N_fm_L-1, 0, complex(t_fm_L, 0.0), Hamiltonian); // Left FM contact
    cmplxSparse.ScaleRangeSparseMatrixIP(2*N_fm_L, 2*N_fm_L+1, 0, complex(0.5*(t_fm_L+t_ox), 0.0), Hamiltonian); // Interface of left FM contact with barrier
    cmplxSparse.ScaleRangeSparseMatrixIP(2*(N_fm_L+1), 2*(N_fm_L+N_ox)+1, 0, complex(t_ox, 0.0), Hamiltonian); // Barrier
    cmplxSparse.ScaleRangeSparseMatrixIP(2*(N_fm_L+N_ox+1), 2*(N_fm_L+N_ox)+3, 0, complex(0.5*(t_ox+t_fm_R), 0.0), Hamiltonian); // Interface of right FM contact with barrier
    cmplxSparse.ScaleRangeSparseMatrixIP(2*(N_fm_L+N_ox)+4, 2*(N_fm_L+N_ox+N_fm_R)+3, 0, complex(t_fm_R, 0.0), Hamiltonian); // Right FM contact

    // Build the upper diagonal
    cmplxSparse.ScaleRangeSparseMatrixIP(0, 2*N_fm_L-1, 2, complex(t_fm_L, 0.0), Hamiltonian);
    cmplxSparse.ScaleRangeSparseMatrixIP(2*N_fm_L, 2*(N_fm_L+N_ox)+1, 2, complex(t_ox, 0.0), Hamiltonian);
    cmplxSparse.ScaleRangeSparseMatrixIP(2*(N_fm_L+N_ox+1), 2*(N_fm_L+N_ox+N_fm_R)+1, 2, complex(t_fm_R, 0.0), Hamiltonian);

    // Build the lower diagonal
    cmplxSparse.ScaleRangeSparseMatrixIP(2, 2*N_fm_L+1, -2, complex(t_fm_L, 0.0), Hamiltonian);
    cmplxSparse.ScaleRangeSparseMatrixIP(2*(N_fm_L+1), 2*(N_fm_L+N_ox)+3, -2, complex(t_ox, 0.0), Hamiltonian);
    cmplxSparse.ScaleRangeSparseMatrixIP(2*(N_fm_L+N_ox+2), 2*(N_fm_L+N_ox+N_fm_R)+3, -2, complex(t_fm_R, 0.0), Hamiltonian);

    // Include barrier
    cmplxSparse.AddBarrierProfile(N_fm_L, N_ox, Ub, Hamiltonian);

    // Include band splitting to left FM contact
    mx_, my_, mz_ := math.Sin(th_F)*math.Cos(ph_F), math.Sin(th_F)*math.Sin(ph_F), math.Cos(th_F);
    Hamiltonian = cmplxSparse.AddBandSplitLeftFM(mx_, my_, mz_, E_split_L, N_fm_L, Hamiltonian);

    // Include band splitting to right FM contact
    mx_, my_, mz_ = math.Sin(th_H)*math.Cos(ph_H), math.Sin(th_H)*math.Sin(ph_H), math.Cos(th_H);
    Hamiltonian = cmplxSparse.AddBandSplitRightFM(mx_, my_, mz_, E_split_R, N_fm_L, Hamiltonian);

    cmplxSparse.PrintSparseMatrix(Hamiltonian);

}
