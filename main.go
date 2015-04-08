/* Main code for NEGF based electron transport simulator for magnetic tunnel junction */

package main

import (
    "math"
    "fmt"
    "github.com/negf_mtj/negf_mtj/cmplxSparse"
    "github.com/negf_mtj/negf_mtj/vecQuad"
    "github.com/negf_mtj/negf_mtj/utils"
)

var (
    m_ox, m_fm_L, m_fm_R	float64;
    aSpace              	float64;
    d_ox, d_fm_L, d_fm_R	float64;
    Ub, E_Fermi				float64;
    E_split_L, E_split_R	float64;
    Vmtj                	float64;
    th_F,ph_F               float64;
    th_H,ph_H       		float64;
    Temperature             float64;
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

    // Set material parameters
    m_ox = float64(0.315);
    m_fm_L = float64(0.72);
    m_fm_R = m_fm_L;

    Ub = 0.79;
    E_split_L = 0.930;
    E_split_R = 0.930;
    E_Fermi = 2.25;

    fmt.Println("----------------------------------------");
    fmt.Printf("Ub = %g, E_split_L = %g, E_split_R = %g\n", Ub, E_split_L, E_split_R);

    // Configuration parameters
    th_F, ph_F = utils.Pi/1.0, 0.0;    
    th_H, ph_H = 0.0, 0.0;
    Vmtj, Temperature = 0.0, 300.0;

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

    fmt.Printf("t_fm_L = %f, t_ox = %f, t_fm_R = %f\n",t_fm_L, t_ox, t_fm_R);
    fmt.Printf("2t_fm_L = %f, 2t_ox = %f, 2t_fm_R = %f\n",2*t_fm_L, 2*t_ox, 2*t_fm_R);
    fmt.Println("----------------------------------------");

    // Create data structure for use with integration functions.
    ProblemSet := vecQuad.CreateIntegStruct();

    // Initialize a base matrix template
    HamBuffer := cmplxSparse.New();
    cmplxSparse.MakeHamTriDiag( N_fm_L+N_fm_R+N_ox+2, HamBuffer);

    // Build the main diagonal first
    cmplxSparse.ScaleRangeSparseMatrix(0, 2*N_fm_L-1, 0, complex(t_fm_L, 0.0), HamBuffer); // Left FM contact
    cmplxSparse.ScaleRangeSparseMatrix(2*N_fm_L, 2*N_fm_L+1, 0, complex(0.5*(t_fm_L+t_ox), 0.0), HamBuffer); // Interface of left FM contact with barrier
    cmplxSparse.ScaleRangeSparseMatrix(2*(N_fm_L+1), 2*(N_fm_L+N_ox)+1, 0, complex(t_ox, 0.0), HamBuffer); // Barrier
    cmplxSparse.ScaleRangeSparseMatrix(2*(N_fm_L+N_ox+1), 2*(N_fm_L+N_ox)+3, 0, complex(0.5*(t_ox+t_fm_R), 0.0), HamBuffer); // Interface of right FM contact with barrier
    cmplxSparse.ScaleRangeSparseMatrix(2*(N_fm_L+N_ox)+4, 2*(N_fm_L+N_ox+N_fm_R)+3, 0, complex(t_fm_R, 0.0), HamBuffer); // Right FM contact

    // Build the upper diagonal
    cmplxSparse.ScaleRangeSparseMatrix(0, 2*N_fm_L-1, 2, complex(t_fm_L, 0.0), HamBuffer);
    cmplxSparse.ScaleRangeSparseMatrix(2*N_fm_L, 2*(N_fm_L+N_ox)+1, 2, complex(t_ox, 0.0), HamBuffer);
    cmplxSparse.ScaleRangeSparseMatrix(2*(N_fm_L+N_ox+1), 2*(N_fm_L+N_ox+N_fm_R)+1, 2, complex(t_fm_R, 0.0), HamBuffer);

    // Build the lower diagonal
    cmplxSparse.ScaleRangeSparseMatrix(2, 2*N_fm_L+1, -2, complex(t_fm_L, 0.0), HamBuffer);
    cmplxSparse.ScaleRangeSparseMatrix(2*(N_fm_L+1), 2*(N_fm_L+N_ox)+3, -2, complex(t_ox, 0.0), HamBuffer);
    cmplxSparse.ScaleRangeSparseMatrix(2*(N_fm_L+N_ox+2), 2*(N_fm_L+N_ox+N_fm_R)+3, -2, complex(t_fm_R, 0.0), HamBuffer);

    // Include barrier
    cmplxSparse.AddBarrierProfile(N_fm_L, N_ox, -1.0*Ub, HamBuffer);

    // Include band splitting to left FM contact
    mx_, my_, mz_ := math.Sin(th_F)*math.Cos(ph_F), math.Sin(th_F)*math.Sin(ph_F), math.Cos(th_F);
    HamBuffer = cmplxSparse.AddBandSplitLeftFM(mx_, my_, mz_, -1.0*E_split_L, N_fm_L, HamBuffer);

    // Include band splitting to right FM contact
    mx_, my_, mz_ = math.Sin(th_H)*math.Cos(ph_H), math.Sin(th_H)*math.Sin(ph_H), math.Cos(th_H);
    HamBuffer = cmplxSparse.AddBandSplitRightFM(mx_, my_, mz_, -1.0*E_split_R, N_fm_L, HamBuffer);

    // Calculate matrices for basis transformation. Done here since it is only needed once.
    ProblemSet.SetHamiltonian(HamBuffer);
    ProblemSet.SetParams(Vmtj, E_Fermi, Temperature, E_split_L, E_split_R, m_fm_L, m_ox, m_fm_R, N_fm_L, N_ox, N_fm_R, cmplxSparse.BasisTransform(th_F, ph_F), cmplxSparse.BasisTransform(th_H, ph_H));
    ProblemSet.SetMu();

    cmplxSparse.PrintSparseMatrix(ProblemSet.ReturnHamiltonianPtr());

    // TODO: integrate over mode energies
    // i.e. We should integrate over the region E_mode = [0, +inf)
    E_mode := 0.0;
    ProbDup := vecQuad.CreateIntegStruct();
    ProbDup.CopyIntegStruct(ProblemSet);
    currents := ProbDup.NEGF_ModeIntegFunc(E_mode);

    for idx0 := 0; idx0 < len(currents); idx0++ {
        fmt.Println("currents[", idx0, "] =", currents[idx0]);
    }
}

