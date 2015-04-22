
package main

import (
	"fmt"
	"github.com/negf_mtj/negf_mtj/cmplxSparse"
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
	fmt.Println("Testing self-energy calculation in cmplxSparse package...\n\n");

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
    fmt.Printf("Ub = %.15g, E_Fermi = %.15g, E_split_L = %.15g, E_split_R = %.15g\n", Ub, E_Fermi, E_split_L, E_split_R);

    // Configuration parameters
    th_F, ph_F = utils.Pi/2.0, 0.0;    
    th_H, ph_H = 0.0, 0.0;
    Vmtj, Temperature = 0.20, 300.0;

    fmt.Println("----------------------------------------");
    fmt.Println("Left FM:");
    fmt.Printf("th = %.15g, ph = %.15g\n", th_F, ph_F);

    fmt.Println("----------------------------------------");
    fmt.Println("Right FM:")
    fmt.Printf("th = %.15g, ph = %.15g\n", th_H, ph_H);
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
    fmt.Printf("Grid spacing = %.15g, N_ox = %d, N_fm_L = %d, N_fm_R = %d\n", aSpace, N_ox, N_fm_L, N_fm_R);
    fmt.Println("----------------------------------------");

    t_base := utils.Hbar*utils.Hbar/2/utils.Echarge/utils.M0/aSpace/aSpace;
    t_fm_L := t_base/m_fm_L;
    t_fm_R := t_base/m_fm_R;
    t_ox := t_base/m_ox;

    fmt.Printf("m_fm_L = %.15g, m_ox = %.15g, m_fm_R = %.15g\n", m_fm_L, m_ox, m_fm_R);
    fmt.Println("----------------------------------------");
    fmt.Printf("t_fm_L = %.15g, t_ox = %.15g, t_fm_R = %.15g\n", t_fm_L, t_ox, t_fm_R);
    fmt.Printf("2t_fm_L = %.15g, 2t_ox = %.15g, 2t_fm_R = %.15g\n", 2*t_fm_L, 2*t_ox, 2*t_fm_R);
    fmt.Println("----------------------------------------");

    BT_Ptr := cmplxSparse.BasisTransform(th_F, ph_F)
    BT_L := *BT_Ptr;
    fmt.Printf("BT Matrix (Left):\n %.15g  %.15g\n %.15g  %.15g\n\n", BT_L[0][0], BT_L[0][1], BT_L[1][0], BT_L[1][1]);

    SelfEnergyLeft := cmplxSparse.SelfEnergyEntries(1.8915, 2.1897e-07, 0.0, E_split_L, Vmtj, complex(t_fm_L,0.0), BT_Ptr);
    fmt.Printf("Self-energy (Left):\n %.15g  %.15g\n %.15g  %.15g\n\n", SelfEnergyLeft[0][0], SelfEnergyLeft[0][1], SelfEnergyLeft[1][0], SelfEnergyLeft[1][1]);

    BT_Ptr = cmplxSparse.BasisTransform(th_H, ph_H)
    BT_L = *BT_Ptr;
    fmt.Printf("BT Matrix (Left):\n %.15g  %.15g\n %.15g  %.15g\n\n", BT_L[0][0], BT_L[0][1], BT_L[1][0], BT_L[1][1]);

    SelfEnergyLeft = cmplxSparse.SelfEnergyEntries(1.8915, 2.1897e-07, 0.0, E_split_R, -1.0*Vmtj, complex(t_fm_R,0.0), BT_Ptr);
    fmt.Printf("Self-energy (Right):\n %.15g  %.15g\n %.15g  %.15g\n\n", SelfEnergyLeft[0][0], SelfEnergyLeft[0][1], SelfEnergyLeft[1][0], SelfEnergyLeft[1][1]);

}