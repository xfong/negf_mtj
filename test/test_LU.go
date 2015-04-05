
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

	LU_mat := cmplxSparse.SparseCopy(tmp);
	cmplxSparse.SparseDiagLU(LU_mat);
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
    fmt.Println("Accessing matrix elements of Hamiltonian (m,n):");

    for idx0 := 0; idx0 < matrixSize; idx0++ {
        for idx1 := 0; idx1 < matrixSize; idx1++ {
            test := cmplxSparse.AccessMatrix(idx0,idx1,tmp);
            fmt.Printf("%f  ", test);
        }
        fmt.Printf("\n");
    }
    fmt.Printf("\n");

    tmp0 := cmplxSparse.SparseCopy(tmp);
    cmplxSparse.SparseDiagLU(tmp0);
    fmt.Println("Accessing matrix elements of LU (m,n):");
    for idx0 := 0; idx0 < matrixSize; idx0++ {
        for idx1 := 0; idx1 < matrixSize; idx1++ {
            test := cmplxSparse.AccessMatrix(idx0,idx1,tmp0);
            fmt.Printf("%f  ", test);
        }
        fmt.Printf("\n");
    }
    fmt.Printf("\n");

    InvBuffer := make([][]complex128, matrixSize);
    for idx0 := 0; idx0 < matrixSize; idx0++ {
        InvBuffer[idx0] = make([]complex128, matrixSize);
        for idx1 := 0; idx1 < matrixSize; idx1++ {
            InvBuffer[idx0][idx1] = 0.0 + 0.0i;
        }
    }

    for idx0 := 0; idx0 < matrixSize; idx0++ {
        InvBuffer[idx0][idx0] = 1.0 + 0.0i;
        BufferMatrix := cmplxSparse.SparseDiagLinearSolver(tmp0, InvBuffer[idx0]);
        for idx1 := 0; idx1 < matrixSize; idx1++ {
            InvBuffer[idx0][idx1] = BufferMatrix[idx1];
        }
        fmt.Printf("\n");
    }

    fmt.Println("Accessing matrix elements if inverse (m,n):");

    for idx0 := 0; idx0 < matrixSize; idx0++ {
        for idx1 := 0; idx1 < matrixSize; idx1++ {
            fmt.Printf("%f  ", InvBuffer[idx1][idx0]);
        }
        fmt.Printf("\n");
    }
	fmt.Println("\n");
}
