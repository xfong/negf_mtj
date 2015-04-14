// Package for sparse complex matrix operations for NEGF calculations

package vecQuad

import (
    "math"
	"math/cmplx"
    "errors"
	"fmt"
	"github.com/negf_mtj/negf_mtj/cmplxSparse"
)

type IntegStruct struct {
    E_mode, V_MTJ, E_Fermi      float64;
    Temperature, deltaL     	float64;
    deltaR, m_fmL, m_ox     	float64;
    m_fmR, mu1, mu2, f1, f2 	float64;
    f1_prime, f2_prime      	float64;
    N_fmL, N_ox, N_fmR      	int;
    BT_Mat_L, BT_Mat_R      	*[2][2]complex128;
    Hamiltonian             	*cmplxSparse.SparseMat
}

var (
    InitialIntervalCount, MaxIntervalCount      int;
    AbsTol, RelTol                              float64;
    Nodes, Wt15, Wt7, EWts                      [15]float64;
);

func init() {
    InitialIntervalCount = 10;
    MaxIntervalCount = 16384;
    AbsTol = 1e-10;
    RelTol = 1e-8;
    // rule.nodes
    Nodes = [...]float64{ -0.9914553711208126, -0.9491079123427585, -0.8648644233597691, -0.7415311855993944, -0.5860872354676911, -0.4058451513773972, -0.2077849550078985, 0, 0.2077849550078985, 0.4058451513773972, 0.5860872354676911, 0.7415311855993944, 0.8648644233597691, 0.9491079123427585, 0.9914553711208126};
    // rule.HighWeights
    Wt15 = [...]float64{ 0.02293532201052922, 0.06309209262997855, 0.1047900103222502, 0.1406532597155259, 0.1690047266392679, 0.1903505780647854, 0.2044329400752989, 0.2094821410847278, 0.2044329400752989, 0.1903505780647854, 0.1690047266392679, 0.1406532597155259, 0.1047900103222502, 0.06309209262997855, 0.02293532201052922};
    // rule.LowWeights
    Wt7 = [...]float64{ 0, 0.1294849661688697, 0, 0.2797053914892767, 0, 0.3818300505051189, 0, 0.4179591836734694, 0, 0.3818300505051189, 0, 0.2797053914892767, 0, 0.1294849661688697, 0};
    for idx0 := 0; idx0 < 15; idx0++ {
        EWts[idx0] = Wt15[idx0] - Wt7[idx0];
    }
}

// TODO: Since the integration over energies use the same function parameters,
// (i.e., for each function call, we only need to change the mode energy and 
// value of energy at which we are evaluating the function) and we know how to
// implement quadrature methods for functions taking scalar inputs and returns
// a scalar, we can consider using a data structure to hold the function
// parameters, and use a function call to set the energies. These functions
// can be called SetAndIntegMode(float64) and SetAndIntegEnergy(float64),
// respectively. To support parallization, these functions may need a local
// copy of the data structure holding the input parameters.

// Function that allows main package to initialize IntegStruct data structure
func CreateIntegStruct() *IntegStruct {
    s := new(IntegStruct);
    return s;
}

// Function needed to perform quadrature. The idea is that numerically, we want
// to perform the integration over the interval [a, inf). However, we map the
// entire function down to the inteval [0, 1]. In the first step, we map
// [a, inf) -> [0, inf). Finally, we map [0, inf) -> [0, 1]
func IntervalA2InfTransform(A, t float64) (x, w float64) {
    tt := t / (1 - t);
    x = A + tt*tt;
    w = 2 * tt / ((1 - t) * (1 - t));
    return;
}

func (t *IntegStruct) CopyIntegStruct(s *IntegStruct) {
    t.E_mode = s.E_mode;
    t.mu1 = s.mu1;
    t.mu2 = s.mu2;
    t.f1 = s.f1;
    t.f2 = s.f2;
    t.f1_prime = s.f1_prime;
    t.f2_prime = s.f2_prime;
    t.V_MTJ = s.V_MTJ;
    t.Temperature = s.Temperature;
    t.deltaL = s.deltaL;
    t.deltaR = s.deltaR;
    t.m_fmL = s.m_fmL;
    t.m_ox = s.m_ox;
    t.m_fmR = s.m_fmR;
    t.N_fmL = s.N_fmL;
    t.N_ox = s.N_ox;
    t.N_fmR = s.N_fmR;
    t.BT_Mat_L = s.BT_Mat_L;
    t.BT_Mat_R = s.BT_Mat_R;
    t.Hamiltonian = s.Hamiltonian;
    return;
}

func (s *IntegStruct) SetMode( ModeEnergy float64 ) {
    s.E_mode = ModeEnergy;
}

func (s *IntegStruct) SetHamiltonian( Hamiltonian *cmplxSparse.SparseMat ) {
    s.Hamiltonian = Hamiltonian;
}

func (s *IntegStruct) ReturnHamiltonianPtr( ) *cmplxSparse.SparseMat {
    return s.Hamiltonian;
}

func (s *IntegStruct) SetParams( V_MTJ, E_Fermi, Temperature, deltaL, deltaR, m_fmL, m_ox, m_fmR float64, N_fmL, N_ox, N_fmR int, BT_Mat_L, BT_Mat_R *[2][2]complex128 ) {
    s.V_MTJ, s.E_Fermi, s.Temperature, s.deltaL, s.deltaR = V_MTJ, E_Fermi, Temperature, deltaL, deltaR;
    s.m_fmL, s.m_ox, s.m_fmR = m_fmL, m_ox, m_fmR;
    s.N_fmL, s.N_ox, s.N_fmR = N_fmL, N_ox, N_fmR;
    s.BT_Mat_L, s.BT_Mat_R = BT_Mat_L, BT_Mat_R;
}

func (s *IntegStruct) SetMu( ) {
    s.mu1, s.mu2 = s.E_Fermi + 0.5*s.V_MTJ, s.E_Fermi - 0.5*s.V_MTJ;
    s.f1, s.f2 = cmplxSparse.FermiEnergy(s.E_Fermi, s.mu1, s.Temperature), cmplxSparse.FermiEnergy(s.E_Fermi, s.mu2, s.Temperature);
    s.f1_prime, s.f2_prime = 1.0 - s.f1, 1.0 - s.f2;
}

// This function is called to perform integration over mode energies.
func (s *IntegStruct) NEGF_AutoModeInteg() *[4]float64 {
    ProbDup := CreateIntegStruct();
    ProbDup.CopyIntegStruct(s);
    E_mode := 0.0;
    fmt.Println("Checking internal structure: Nodes[", len(Nodes), "], Wt15[", len(Wt15), "], Wt7[", len(Wt7), "]");
    t := ProbDup.NEGF_ModeIntegFunc(E_mode);
    return t;
}

// TODO: inside this function, we need to integrate over energy. This integration should call the
// function NEGF_EnergyIntegFunc
func (s *IntegStruct) NEGF_ModeIntegFunc( E_mode float64 ) *[4]float64 {
    // Initialize return value
    t := new([4]float64);
    s.SetMode(E_mode);

    // TODO: Begin integration over energy space

    fmt.Println("Inside NEGF_ModeIntegFunc. Calling NEGF_EnergyIntegFunc...")
    // TODO: the integration over energy should occur from the highest minimum conduction band to infinity
    // i.e. the integral should be over the region Eval in [max(mu1 - E_Fermi, mu2 - E_Fermi), +inf)
    t = s.NEGF_EnergyIntegFunc(0.0);

    // Return result
    return t;
}

// This function is called during integration over energies.
func (s *IntegStruct) NEGF_EnergyIntegFunc( EnergyValue float64 ) *[4]float64 {
    // Initialize the return value
    t := new([4]float64);

    // Get t0 on left and right of Hamiltonian matrix
    MatrixSize := len(s.Hamiltonian.Data);
    MaxDiagIdx := len(s.Hamiltonian.Data[0]);
    MainDiagIdx := (MaxDiagIdx - 1)/2;

    // Store f1*gam1 and f2*gam2. Note that these matrices are non-zero at the top
    // left and bottom right, respectively. Hence, it is sufficient to store them
    // as just 2 x 2 matrices.
    f1gam1, f2gam2 := new([2][2]complex128), new([2][2]complex128);

    // Calculate f1*gam1 first and store in a buffer
    SelfEnergyBuffer := cmplxSparse.SelfEnergyEntries(EnergyValue, s.E_mode, 0.0, s.deltaL, -0.5*s.V_MTJ, s.Hamiltonian.Data[0][MaxDiagIdx-1], s.BT_Mat_L);
    f1gam1[0][0] = complex(0.0, s.f1) * (SelfEnergyBuffer[0][0] - cmplx.Conj(SelfEnergyBuffer[0][0]));
    f1gam1[0][1] = complex(0.0, s.f1) * (SelfEnergyBuffer[0][1] - cmplx.Conj(SelfEnergyBuffer[1][0]));
    f1gam1[1][0] = complex(0.0, s.f1) * (SelfEnergyBuffer[1][0] - cmplx.Conj(SelfEnergyBuffer[0][1]));
    f1gam1[1][1] = complex(0.0, s.f1) * (SelfEnergyBuffer[1][1] - cmplx.Conj(SelfEnergyBuffer[1][1]));

    // Adjust Hamiltonian with f1gam1 to obtain the inverse of Green's function matrix
    InvGMatrix := cmplxSparse.SparseCopy(s.Hamiltonian);
    InvGMatrix.Data[0][MainDiagIdx] -= SelfEnergyBuffer[0][0]
    InvGMatrix.Data[0][MainDiagIdx+1] -= SelfEnergyBuffer[0][1]
    InvGMatrix.Data[1][MainDiagIdx-1] -= SelfEnergyBuffer[1][0]
    InvGMatrix.Data[1][MainDiagIdx] -= SelfEnergyBuffer[1][1]

    // Calculate f2*gam2 next and store in a buffer
    SelfEnergyBuffer = cmplxSparse.SelfEnergyEntries(EnergyValue, s.E_mode, 0.0, s.deltaR, 0.5*s.V_MTJ, s.Hamiltonian.Data[MatrixSize-MaxDiagIdx+MainDiagIdx][MaxDiagIdx-1], s.BT_Mat_R);
    f2gam2[0][0] = complex(0.0, s.f2) * (SelfEnergyBuffer[0][0] - cmplx.Conj(SelfEnergyBuffer[0][0]));
    f2gam2[0][1] = complex(0.0, s.f2) * (SelfEnergyBuffer[0][1] - cmplx.Conj(SelfEnergyBuffer[1][0]));
    f2gam2[1][0] = complex(0.0, s.f2) * (SelfEnergyBuffer[1][0] - cmplx.Conj(SelfEnergyBuffer[0][1]));
    f2gam2[1][1] = complex(0.0, s.f2) * (SelfEnergyBuffer[1][1] - cmplx.Conj(SelfEnergyBuffer[1][1]));

    // Final adjustment of Hamiltonian with f2gam2 to obtain the inverse of Green's function matrix
    InvGMatrix.Data[MatrixSize-2][MainDiagIdx] -= SelfEnergyBuffer[0][0]
    InvGMatrix.Data[MatrixSize-2][MainDiagIdx+1] -= SelfEnergyBuffer[0][1]
    InvGMatrix.Data[MatrixSize-1][MainDiagIdx-1] -= SelfEnergyBuffer[1][0]
    InvGMatrix.Data[MatrixSize-1][MainDiagIdx] -= SelfEnergyBuffer[1][1]

    // Calculate Green's Function matrix
    GMatrix := cmplxSparse.CalcGreensFunc(EnergyValue, InvGMatrix);

    // Calculate Gn = G * (f1*gam1 + f2*gam2) * (G^+)...
    // Since f1*gam1 and f2*gam2 are sparse matrices, we do not need to do a full matrix multiplication
    // Calculate the non-zero portion of (f1*gam1 + f2*gam2) * (G^+) and store in RearMultBuffer:
    // i.e. First two rows of the product are obtained from f1*gam1*(G^+) and are stored as the first
    // two rows of RearMultBuffer. The last two rows of the product are obtained from f1*gam2*(G^+) and
    // are stored as the last two rows of RearMultBuffer.
    RearMultBuffer := make([][]complex128,4);
    RearIdx0, RearIdx1 := MatrixSize - 2, MatrixSize - 1;
    // Initialize buffer storage
    for idx0 := 0; idx0 < 4; idx0++ {
        RearMultBuffer[idx0] = make([]complex128, MatrixSize);
    }
    // Calculate each row of RearMultBuffer
    for idx0 := 0; idx0 < MatrixSize; idx0++ {
        RearMultBuffer[0][idx0] = f1gam1[0][0] * GMatrix[idx0][0] + f1gam1[0][1] * GMatrix[idx0][1];
        RearMultBuffer[1][idx0] = f1gam1[1][0] * GMatrix[idx0][0] + f1gam1[1][1] * GMatrix[idx0][1];
        if ((idx0 > 1) && (idx0 < RearIdx0)) {
            RearMultBuffer[2][idx0] = f2gam2[0][0] * GMatrix[idx0][2] + f2gam2[0][1] * GMatrix[idx0][3];
            RearMultBuffer[3][idx0] = f2gam2[1][0] * GMatrix[idx0][2] + f2gam2[1][1] * GMatrix[idx0][3];
        } else {
            RearMultBuffer[2][idx0] = f2gam2[0][0] * GMatrix[idx0][RearIdx0] + f2gam2[0][1] * GMatrix[idx0][RearIdx1];
            RearMultBuffer[3][idx0] = f2gam2[1][0] * GMatrix[idx0][RearIdx0] + f2gam2[1][1] * GMatrix[idx0][RearIdx1];
        }
    }

    // Finally, compute the required entries of Gn.
    GnF, GnR := new([2][2]complex128), new([2][2]complex128);
    GnF[0][0], GnF[0][1], GnF[1][0], GnF[1][1] = 0.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i;
    GnR[0][0], GnR[0][1], GnR[1][0], GnR[1][1] = 0.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i;
    PtIdx0 := 2*s.N_fmL;
    PtIdx1, PtIdx2, PtIdx3 := PtIdx0 - 1, PtIdx0 + 1, PtIdx0 + 2;

    // GMatrix is actually the transpose of G. Due to the nature of (f1*gam1 + f2*gam2)
    // the entries of G that we need are:
    // 1. The first two rows.
    // 2. The first two columns.
    // 3. The last two rows.
    // 4. The last two columns.
    // Scanning the first index of GMatrix scans through the columns of G.
    // The second index of GMatrix selects the row.
    for idx0 := 0; idx0 < 4; idx0++ {
        if (idx0 < 2) {
            // When multiplication involves first two rows of RearMultBuffer,
            // we need to access first two columns of G. 
            GnF[0][0] += GMatrix[idx0][PtIdx1] * RearMultBuffer[idx0][PtIdx2];
            GnF[0][1] += GMatrix[idx0][PtIdx1] * RearMultBuffer[idx0][PtIdx3];
            GnF[1][0] += GMatrix[idx0][PtIdx0] * RearMultBuffer[idx0][PtIdx2];
            GnF[1][1] += GMatrix[idx0][PtIdx0] * RearMultBuffer[idx0][PtIdx3];
            GnR[0][0] += GMatrix[idx0][PtIdx2] * RearMultBuffer[idx0][PtIdx1];
            GnR[0][1] += GMatrix[idx0][PtIdx2] * RearMultBuffer[idx0][PtIdx0];
            GnR[1][0] += GMatrix[idx0][PtIdx3] * RearMultBuffer[idx0][PtIdx1];
            GnR[1][1] += GMatrix[idx0][PtIdx3] * RearMultBuffer[idx0][PtIdx0];
        } else {
            // When multiplication involves last two rows of RearMultBuffer,
            // we need to access last two columns of G. 
            targIdx := MatrixSize - 4;
            GnF[0][0] += GMatrix[targIdx + idx0][PtIdx1] * RearMultBuffer[idx0][PtIdx2];
            GnF[0][1] += GMatrix[targIdx + idx0][PtIdx1] * RearMultBuffer[idx0][PtIdx3];
            GnF[1][0] += GMatrix[targIdx + idx0][PtIdx0] * RearMultBuffer[idx0][PtIdx2];
            GnF[1][1] += GMatrix[targIdx + idx0][PtIdx0] * RearMultBuffer[idx0][PtIdx3];
            GnR[0][0] += GMatrix[targIdx + idx0][PtIdx2] * RearMultBuffer[idx0][PtIdx1];
            GnR[0][1] += GMatrix[targIdx + idx0][PtIdx2] * RearMultBuffer[idx0][PtIdx0];
            GnR[1][0] += GMatrix[targIdx + idx0][PtIdx3] * RearMultBuffer[idx0][PtIdx1];
            GnR[1][1] += GMatrix[targIdx + idx0][PtIdx3] * RearMultBuffer[idx0][PtIdx0];
        }
    }

    // We know the Hamiltonian is stored in diagonal format. Extract the wanted values of Hamiltonian
    // and multiply accordingly to GnF and GnR.
    HamF, HamR, MatrixBuffer := new([2][2]complex128), new([2][2]complex128), new([2][2]complex128);
    HamR[0][1] = s.Hamiltonian.Data[PtIdx2][MainDiagIdx-1];
    HamR[0][0] = s.Hamiltonian.Data[PtIdx2][MainDiagIdx-2];
    HamR[1][1] = s.Hamiltonian.Data[PtIdx3][MainDiagIdx-2];
    HamF[1][0] = s.Hamiltonian.Data[PtIdx0][MainDiagIdx+1];
    HamF[0][0] = s.Hamiltonian.Data[PtIdx1][MainDiagIdx+2];
    HamF[1][1] = s.Hamiltonian.Data[PtIdx0][MainDiagIdx+2];
    if (MaxDiagIdx > 6) {
        HamR[1][0] = s.Hamiltonian.Data[PtIdx3][MainDiagIdx-3];
        HamF[0][1] = s.Hamiltonian.Data[PtIdx1][MainDiagIdx+3];
    } else {
        HamR[1][0] = 0.0 + 0.0i;
        HamF[0][1] = 0.0 + 0.0i;
    }

    // Calculate term from H*Gn - Gn*H;
    MatrixBuffer[0][0] = HamR[0][0] * GnF[0][0] + HamR[0][1] * GnF[1][0] - GnR[0][0] * HamF[0][0] - GnR[0][1] * HamF[1][0];
    MatrixBuffer[0][1] = HamR[0][0] * GnF[0][1] + HamR[0][1] * GnF[1][1] - GnR[0][0] * HamF[0][1] - GnR[0][1] * HamF[1][1];
    MatrixBuffer[1][0] = HamR[1][0] * GnF[0][0] + HamR[1][1] * GnF[1][0] - GnR[1][0] * HamF[0][0] - GnR[1][1] * HamF[1][0];
    MatrixBuffer[1][1] = HamR[1][0] * GnF[0][1] + HamR[1][1] * GnF[1][1] - GnR[1][0] * HamF[0][1] - GnR[1][1] * HamF[1][1];

    // Calculate current density, and components of spin current
    t[0] = real(MatrixBuffer[0][0] + MatrixBuffer[1][1]); // Charge current density
    t[1] = real(MatrixBuffer[1][0] + MatrixBuffer[0][1]);
    t[2] = real(1.0i * MatrixBuffer[0][1] - 1.0i * MatrixBuffer[1][0]);
    t[3] = real(MatrixBuffer[0][0] - MatrixBuffer[1][1]); // z-component of spin current density

    // Return result
    return t;
}

/*
func IntegralCalc2B(f func(float64) []float64, a, b float64) {
}
*/
// Function for performing integration over the interval [a, inf)
// using Gauss-Konrod quadrature (7-point estimator,
// 15 point corrector). This is used in MATLAB as well.
func IntegralCalc2Inf(f func(float64) *[]float64, a float64, expectSize int) (q, errbnd []float64) {
    // Generate 10 subintervals first
    nsubs := 10;
    var (
        subs, qsubs, errsubs, xx                             	[][]float64;
        qsubsk, errsubsk, t, w, fxj, q_ok, err_ok, err_not_ok	[]float64;
        fTmp                                                    *[]float64;
        NNodes, nleft                                           int;
        too_close                                               bool;
    );
    fxj = make([]float64, expectSize);
    // Set up buffer for accumulating results over one subinterval
    qsubsk, errsubsk = make([]float64, expectSize), make([]float64, expectSize);

    // subs[0][nn] and subs[1][nn] stores the start and end points of the
    // nn-th subinterval, respectively. In the first step, we generate 10
    // subintervals in [0, 1].
    subs = make([][]float64, 2);
    subs[0] = make([]float64, nsubs);
    subs[1] = make([]float64, nsubs);

    // midpts[nn] and halfh[nn] stores the midpts and length of the nn-th
    // subinterval, respectively.
    midpts, halfh := make([]float64, nsubs), make([]float64, nsubs);

    // Set up arrays for the first subinterval
    subs[0][0] = 0.0;
    subs[1][0] = 0.1;

    // Finish setting up for the rest of the subintervals
    for idx0 := 1; idx0 < nsubs; idx0++ {
        subs[1][idx0] = subs[1][idx0-1] + 0.1;
        subs[0][idx0] = subs[1][idx0-1];
    }

    // Initialize more buffers
    q, q_ok, err_ok, err_not_ok, errbnd = make([]float64, expectSize), make([]float64, expectSize), make([]float64, expectSize), make([]float64, expectSize), make([]float64, expectSize);
    for idx0 := range q {
        q[idx0], q_ok[idx0], err_ok[idx0], err_not_ok[idx0], errbnd[idx0] = 0.0, 0.0, 0.0, 0.0, 0.0;
    }

    // Begin "infinite" loop. Loop breaks out when error tolerances are met
    // or when the interation process meets/fails certain conditions.
    for {

        // Update q with the previous OK value before going through loop
        for idx0 := range q_ok {
            q[idx0] = q_ok[idx0];
        }
        // Set up arrays defining midpoints and half path lengths of every
        // subinterval. subs and nsubs are updated at the end of every
        // iteration. Hence, we need to compute the midpoints and lengths
        // of every subinterval at the beginning of the iteration.
        for idx0 := 0; idx0 < nsubs; idx0++ {
            midpts[idx0], halfh[idx0] = 0.5*(subs[0][idx0] + subs[1][idx0]), 0.5*(subs[1][idx0] - subs[0][idx0]);
        }

        // Set up arrays for storing results over each subinterval
        qsubs, errsubs = make([][]float64, nsubs), make([][]float64, nsubs);
        for idx0 := 0; idx0 < nsubs; idx0++ {
            qsubs[idx0], errsubs[idx0] = make([]float64, expectSize), make([]float64, expectSize);
        }

        // Set up arrays storing the actual nodes at which we are
        // evaluating the results
        NNodes = len(Nodes);
        xx = make([][]float64, NNodes);
        t, w = make([]float64, NNodes), make([]float64, NNodes);

        // Begin going through every subinterval to calculate integral over
        // each of them
        for idx0 := 0; idx0 < nsubs; idx0++ {
            // Zero the buffer prior to scanning through the nodes
            for idx1 := range qsubsk {
                qsubsk[idx1] = 0.0;
                errsubsk[idx1] = 0.0;
            }

            // For each subinterval, scan through the nodes to compute values
            for idx1 := range Nodes {
                xx[idx1] = make([]float64, nsubs);
                xx[idx1][idx0] = Nodes[idx1]*halfh[idx0] + midpts[idx0];
                t[idx1], w[idx1] = IntervalA2InfTransform(a, xx[idx1][idx0]);
                fTmp = f(t[idx1]);
                fxj = *fTmp;
                for idx2 := range qsubsk {
                    qsubsk[idx2] += fxj[idx2] * w[idx1] * Wt15[idx1] * halfh[idx0];
                    errsubsk[idx2] += fxj[idx2] * w[idx1] * EWts[idx1] * halfh[idx0];
                }
            }
            // At this point, we have the estimated integral and the error
            // for the subnterval. So store into the bigger buffer
            for idx1 := range qsubsk {
                qsubs[idx0][idx1] = qsubsk[idx1];
                errsubs[idx0][idx1] = errsubsk[idx1];
                q[idx1] += qsubsk[idx1];
            }

            // Terminate and exit if the spacing is too close
            too_close = checkSpacing(&t);
            if too_close {
                break;
            }
        }
        // Terminate and exit if the spacing is too close. The previous
        // too_close check breaks out of the scan through subintervals.
        // This check here breaks out of the entire "infinite" for loop
        if too_close {
            break;
        }

        // Scan through the subintervals and reiterate for those
        // subintervals that are insufficiently accurate
        nleft = 0;
        tmpSub := make([][]float64, 2);
        tmpMids := make([]float64, 1);
        tmpSub[0] = make([]float64, 1);
        tmpSub[1] = make([]float64, 1);

        tol, tolr, tola := make([]float64, expectSize), make([]float64, expectSize), 2.0*AbsTol;
        for idx0 := range q {
            tol[idx0] = RelTol * math.Abs(q[idx0]);
            tolr[idx0] = 2*tol[idx0];
        }
        for idx0 := 0; idx0 < nsubs; idx0++ {
            abserrsubsk := make([]float64, expectSize);
            flagSum0 := int(0);
            for idx1 := range q {
                abserrsubsk[idx1] = math.Abs(errsubs[idx0][idx1]);
                if ((abserrsubsk[idx1] > tolr[idx1] * halfh[idx0]) && (abserrsubsk[idx1] > tola * halfh[idx0])) {
                    flagSum0++;
                }
            }
            if (flagSum0 == 0) {
                for idx1 := 0; idx1 < expectSize; idx1++ {
                    q_ok[idx1] += qsubs[idx0][idx1];
                    err_ok[idx1] += errsubs[idx0][idx1];
                }
            } else {
                if (nleft == 0) {
                    tmpSub[0][0] = subs[0][idx0];
                    tmpSub[1][0] = subs[1][idx0];
                    tmpMids[0] = midpts[idx0];
                    nleft = 1;
                } else {
                    tmpMids = append(tmpMids, midpts[idx0]);
                    tmpSub[0] = append(tmpSub[0], subs[0][idx0]);
                    tmpSub[1] = append(tmpSub[1], subs[1][idx0]);
                    nleft++;
                }
                for idx1 := range err_not_ok {
                    err_not_ok[idx1] += abserrsubsk[idx1];
                }
            }
        }
        // By this point, we have figured out the subintervals that failed
        // tolerance checks. We will divide the subintervals into two and
        // reiterate the "infinite" for loop.
        flagSum1 := int(0);
        for idx0 := range errbnd {
            errbnd[idx0] = math.Abs(err_ok[idx0]) + err_not_ok[idx0];
            if ((errbnd[idx0] > tol[idx0]) && (errbnd[idx0] > AbsTol)) {
                flagSum1++;
            }
        }

        // Break out of infinite loop if we are within error bounds.
        if ((nleft < 1) || (flagSum1 == 0)) {
            break;
        }

        // Dividing subintervals before reiterating
        nsubs = 2 * nleft;
        if (nsubs > MaxIntervalCount ) {
            errors.New("ERROR: MaxIntervalCount reached!");
            break;
        }
        subs_div := make([][]float64, 2);
        subs_div[0] = make([]float64, nsubs);
        subs_div[1] = make([]float64, nsubs);
        for idx0 := 0; idx0 < nleft; idx0++ {
            targIdxL := idx0*2;
            targIdxH := targIdxL+1;
            subs_div[0][targIdxL] = tmpSub[0][idx0];
            subs_div[1][targIdxL] = tmpMids[idx0];
            subs_div[0][targIdxH] = tmpMids[idx0];
            subs_div[1][targIdxH] = tmpSub[1][idx0];
        }
        subs = subs_div;
    }
    return;
}


// Check to ensure spacing between integration nodes is sufficiently large
func checkSpacing(t *[]float64) bool {
    tmpT := *t;
    vecLength := len(tmpT);
    ax := make([]float64, vecLength);
    for idx0 := range tmpT {
        ax[idx0] = math.Abs(tmpT[idx0]);
    }
    adx := make([]float64, (vecLength - 1));
    for idx0 := range adx {
        adx[idx0] = math.Abs(ax[idx0+1] - ax[idx0]);
        if (adx[idx0] < (2.2204460492503131e-14)*math.Max(ax[idx0+1],ax[idx0])) {
            return true;
        }
    }
    return false;
}
