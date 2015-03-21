// Package for sparse complex matrix operations for NEGF calculations

package cmplxSparse

import (
	"errors"
	"math"
)

// Sparse matrix data structure stores a matrix in Diagonal form, which is convenient
// since the Hamiltonian is built from a tridiagonal matrix. Each column of the 2D-array
// given by Data is a diagonal of the actual matrix. Hence, the number of columns must be
// odd. The middle column is the main diagonal. Every column right of the middle is the
// off-diagonal to the right, while every column left of the middle is the off-diagonal
// to the left. The row index of an entry corresponds to the matrix row the entry belongs.
// i.e. the identity matrix is given as:
//      A[i][0] = 1.0
// i.e. the standard 4x4 tridiagonal matrix is given as:
//      A =  0.0 2.0 -1.0
//          -1.0 2.0 -1.0
//          -1.0 2.0 -1.0
//          -1.0 2.0  0.0

type sparseMat struct {
	Data			[][]complex128
}

func init() {
}

func New() *sparseMat {
	tmp := new(sparseMat);
	return tmp;
}

// Function to access matrix elements
func AccessMatrix( m, n int, s *sparseMat ) complex128 {
	if ((m < 0) || (n < 0)) {
		errors.New("Access indices cannot be negative!");
	} else {
		matDim := len(s.Data)
		if ((m >= matDim) || (n >= matDim)) {
			errors.New("Access index out of range!");
		}
	}
	totalDiags := len(s.Data[0]);
	maxOffset := (totalDiags-1)/2;
	if ((n > m+maxOffset) || (n < m-maxOffset)) {
		return 0.0
	} else {
		idx0 := n-m+maxOffset
		return s.Data[m][idx0]
	}
}

// Function to initialize an identity matrix in Diagonal format
func MakeIdentity( matSize int, s *sparseMat ) {
	// TODO: error if matsize is zero and less
	if (matSize <= 0) {
		errors.New("ERROR: Invalid size for matrix (<=0)");
	}
	s.Data = make([][]complex128, matSize);
	for idx0 := 0; idx0 < matSize; idx0++ {
		s.Data[idx0] = make([]complex128,1);
		s.Data[idx0][0] = 1.0;
	}
}

// Function to initialize a standard tridiagonal matrix
func MakeTriDiag ( matSize int, s *sparseMat ) {
	// TODO: error if matSize is zero and less
	if (matSize <= 0) {
		errors.New("ERROR: matrix size must be 1 or larger!");
	}
	s.Data = make([][]complex128, matSize);
	for idx0 := 0; idx0 < matSize; idx0++ {
		if (matSize == 1) {
			s.Data[idx0] = make([]complex128,1);
			s.Data[idx0][0] = 2.0;
		} else {
			s.Data[idx0] = make([]complex128,3);
			if (idx0 == 0) {
				s.Data[idx0][0],s.Data[idx0][1],s.Data[idx0][2] =  0.0, 2.0, -1.0;
			} else if (idx0 == matSize - 1) {
				s.Data[idx0][0],s.Data[idx0][1],s.Data[idx0][2] = -1.0, 2.0,  0.0;
			} else {
				s.Data[idx0][0],s.Data[idx0][1],s.Data[idx0][2] = -1.0, 2.0, -1.0;
			}
		}
	}
}

// Function for creating basic tridiagonal Hamiltonian for MTJ
func MakeHamTriDiag( grdSize int, s *sparseMat ) {
	if (grdSize < 5) {
		errors.New("ERROR: Creating the Hamiltonian requires at least 5 points!");
	}
	s.Data = make([][]complex128, 2*grdSize);
	for idx0 := 0; idx0 < 2*grdSize; idx0++ {
		s.Data[idx0] = make([]complex128,5);
		s.Data[idx0][0], s.Data[idx0][1], s.Data[idx0][2], s.Data[idx0][3], s.Data[idx0][4] = -1.0, 0.0, 2.0, 0.0, -1.0;
	}
	for idx0 := 0; idx0 < 2; idx0++ {
		offSet := 2*grdSize-1-idx0;
		s.Data[idx0][0], s.Data[offSet][4] = 0.0, 0.0;
	}
}

// Function for scalar multiplication of sparse matrix
func ScaleSparseMatrix(A complex128, B *sparseMat) *sparseMat {
	s := B;

	for idx0 := 0 ; idx0 < len(s.Data); idx0++ {
		for idx1 := 0; idx1 < len(s.Data[idx0]); idx1++ {
			s.Data[idx0][idx1] *= A;
		}
	}

	return s
}

// Function for in-place scalar multiplication of sparse matrix
func ScaleSparseMatrixIP(A complex128, s *sparseMat) {
	for idx0 := 0 ; idx0 < len(s.Data); idx0++ {
		for idx1 := 0; idx1 < len(s.Data[idx0]); idx1++ {
			s.Data[idx0][idx1] *= A;
		}
	}
}

// Function for adding applied voltage profile to Hamiltonian
func AddVoltagePotential( N_fm, N_ox int, voltage float64, s *sparseMat ) *sparseMat {
	tmpLength := len(s.Data);
	totalPts := tmpLength/2;
	voltageProfile := make([]float64,totalPts)
	voltageDelta := voltage/float64(N_ox+1);
	t := s;
	for idx0 := 0; idx0 < totalPts; idx0++ {
		if (idx0 <= N_fm) {
			voltageProfile[idx0] = 0.5*voltage;
		} else if (idx0 < N_fm+1+N_ox) {
			voltageProfile[idx0] = voltageProfile[idx0-1]-voltageDelta;
		} else {
			voltageProfile[idx0] = -0.5*voltage;
		}
		t.Data[2*idx0][2] += complex(voltageProfile[idx0], 0.0);
		t.Data[2*idx0+1][2] += complex(voltageProfile[idx0], 0.0);
	}
	return t
}

// Function for adding barrier potential profile to Hamiltonian
func AddBarrierProfile( N_fm, N_ox int, Eb float64, s *sparseMat ) {
	tmpLength := len(s.Data);
	totalPts := tmpLength/2;
	for idx0 := 0; idx0 < totalPts; idx0++ {
		if ((idx0 == N_fm) || (idx0 == N_fm+1+N_ox)) {
			s.Data[2*idx0][2] += complex(0.5*Eb,0);
			s.Data[2*idx0+1][2] += complex(0.5*Eb,0);
		} else if ((idx0 > N_fm) && (idx0 < N_fm+1+N_ox)) {
			s.Data[2*idx0][2] += complex(Eb,0);
			s.Data[2*idx0+1][2] += complex(Eb,0);
		}
	}
}

// Function for adding band splitting for up-spin and down-spin conduction bands on left contact
func AddBandSplitLeftFM(mx, my, mz, deltE float64, N_fm int, s *sparseMat) *sparseMat {
	Mnorm := mx*mx + my*my + mz*mz;
	if (Mnorm == 0) {
		errors.New("ERROR: invalid mx, my and mz combination!");
	} else if (N_fm < 0) {
		errors.New("ERROR: invalid value of N_fm!");
	}
	Mnorm = math.Sqrt(Mnorm);
	m_x := mx/Mnorm;
	m_y := my/Mnorm;
	m_z := mz/Mnorm;
	BT := make([][]complex128,2);
	for idx0 := 0; idx0 < 2; idx0++ {
		BT[idx0] = make([]complex128,2);
	}
	scaleFac := 0.5*deltE;
	BT[0][0] = complex(scaleFac*(1.0 - m_z),0.0);
	BT[1][1] = complex(scaleFac*(1.0 + m_z),0.0);
	BT[0][1] = complex(m_x*scaleFac,-m_y*scaleFac);
	BT[1][0] = complex(m_x*scaleFac,m_y*scaleFac);

	t := s;

	for idx0 := 0; idx0 < N_fm; idx0++ {
		currIdx := 2*idx0;
		t.Data[currIdx][2] += BT[0][0];
		t.Data[currIdx+1][2] += BT[1][1];
		t.Data[currIdx][3] += BT[0][1];
		t.Data[currIdx+1][1] += BT[1][0];
	}
	t.Data[2*N_fm][2] += 0.5*BT[0][0];
	t.Data[2*N_fm+1][2] += 0.5*BT[1][1];
	t.Data[2*N_fm][3] += 0.5*BT[0][1];
	t.Data[2*N_fm+1][1] += 0.5*BT[1][0];

	return t;
}

// Function for adding band splitting for up-spin and down-spin conduction bands on right contact
func AddBandSplitRightFM(mx, my, mz, deltE float64, N_fm int, s *sparseMat) *sparseMat {
	Mnorm := mx*mx + my*my + mz*mz;
	if (Mnorm == 0) {
		errors.New("ERROR: invalid mx, my and mz combination!");
	} else if (N_fm < 0) {
		errors.New("ERROR: invalid value of N_fm!");
	}
	Mnorm = math.Sqrt(Mnorm);
	m_x := mx/Mnorm;
	m_y := my/Mnorm;
	m_z := mz/Mnorm;
	BT := make([][]complex128,2);
	for idx0 := 0; idx0 < 2; idx0++ {
		BT[idx0] = make([]complex128,2);
	}
	scaleFac := 0.5*deltE;
	BT[0][0] = complex(scaleFac*(1.0 - m_z),0.0);
	BT[1][1] = complex(scaleFac*(1.0 + m_z),0.0);
	BT[0][1] = complex(m_x*scaleFac,-m_y*scaleFac);
	BT[1][0] = complex(m_x*scaleFac,m_y*scaleFac);

	t := s;

	grdSz := len(s.Data)/2;
	startIdx := grdSz-N_fm;
	for idx0 := startIdx; idx0 < grdSz; idx0++ {
		currIdx := 2*idx0;
		t.Data[currIdx][2] += BT[0][0];
		t.Data[currIdx+1][2] += BT[1][1];
		t.Data[currIdx][3] += BT[0][1];
		t.Data[currIdx+1][1] += BT[1][0];
	}
	t.Data[2*(startIdx-1)][2] += 0.5*BT[0][0];
	t.Data[2*startIdx-1][2] += 0.5*BT[1][1];
	t.Data[2*(startIdx-1)][3] += 0.5*BT[0][1];
	t.Data[2*startIdx-1][1] += 0.5*BT[1][0];

	return t;
}
