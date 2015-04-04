
package main

import (
	"fmt"
	"github.com/negf_mtj/negf_mtj/cmplxSparse"
)

func main() {
	fmt.Println("Testing cmplxSparse package...\n\n");

	// Test LU factorization...
	testLU(6);

	// Test diagonal matrix inversion...
	testDiagInversion(3);

}


func testLU( n int ) {
        fmt.Println("Initializing Sparse matrix structure\n");
        tmp := cmplxSparse.New();

        matrixSize := int(2*n);
        fmt.Println("Making sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix\n");
        cmplxSparse.MakeHamTriDiag(n, tmp);
	for idx0 := 0; idx0 < matrixSize; idx0++ {
		tmp.Data[idx0][2] -= complex(0.0,5.0);
	}

	LU_mat := cmplxSparse.SparseDiagLU(tmp);
        fmt.Println("Accessing matrix elements (m,n):");

        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,LU_mat);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }
	fmt.Println("\n");
}

func testDiagInversion( n int ) {
    fmt.Println("Initializing Sparse matrix structure\n");
    tmp := cmplxSparse.New();

    matrixSize := int(2*n);
    fmt.Println("Making sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix\n");
    cmplxSparse.MakeHamTriDiag(n, tmp);
	for idx0 := 0; idx0 < matrixSize; idx0++ {
		tmp.Data[idx0][2] -= complex(0.0,5.0);
	}

    InvBuffer := make([][]complex128, matrixSize);
    BVector := make([]complex128, matrixSize);
    for idx0 := 0; idx0 < matrixSize; idx0++ {
        BVector[idx0] = 0;
    }
    for idx0 := 0; idx0 < matrixSize; idx0++ {
        if (idx0 > 0) {
            BVector[idx0-1] = 0.0 + 0.0i;
        }
        BVector[idx0] = 1;
        InvBuffer[idx0] = make([]complex128, matrixSize);
        InvBuffer[idx0] = SparseDiagLinearSolver(tmp, BVector);
    }
    fmt.Println("Accessing matrix elements (m,n):");

    for idx0 := 0; idx0 < matrixSize; idx0++ {
        for idx1 := 0; idx1 < matrixSize; idx1++ {
            fmt.Printf("%f  ", InvBuffer[idx1][idx0]);
        }
        fmt.Printf("\n");
    }
	fmt.Println("\n");
}
