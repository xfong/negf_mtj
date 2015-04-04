/* Main code for NEGF based electron transport simulator for magnetic tunnel junction */

package main

import (
    "math"
    "fmt"
    "github.com/negf_mtj/negf_mtj/cmplxSparse"
    "github.com/negf_mtj/negf_mtj/utils"
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
    fmt.Println("Pi =", utils.Pi)
    fmt.Println("Planck constant =", utils.Hplanck)
    fmt.Println("Reduced Planck constant =", utils.Hbar)
    fmt.Println("Elementary charge =", utils.Echarge)
    fmt.Println("Permeability of free space =", utils.Mu0)
    fmt.Println("Bohr magneton =", utils.MuB)
    fmt.Println("Small imaginary number =", utils.Zplus)
    fmt.Println("Free electron mass =", utils.M0)

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
    th_F, ph_F = utils.Pi/2.0, 0.0;    
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

    t_base := utils.Hbar*utils.Hbar/2/utils.Echarge/utils.M0/aSpace/aSpace;
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

    // Calculate matrices for basis transformation. Done here since it is only needed once.
//    BT_Mat_L := cmplxSparse.BasisTransform(th_F, ph_F);
//    BT_Mat_R := cmplxSparse.BasisTransform(th_H, ph_H);

    cmplxSparse.PrintSparseMatrix(Hamiltonian);

    // TODO: Add potential profile due to applied voltage

    // TODO: Begin integration over mode space
    GreensFunc := cmplxSparse.CalcGreensFunc(0, Hamiltonian);

    GreensSize := len(GreensFunc);
    fmt.Println("Printing Green's function matrix");
    for idx0 := 0; idx0 < GreensSize; idx0++ {
        VecSize := len(GreensFunc[idx0]);
        for idx1 := 0; idx1 < VecSize; idx1++ {
            fmt.Printf("%f  ", GreensFunc[idx0][idx1]);
        }
        fmt.Printf("\n");
    }
    fmt.Printf("\n");
}

