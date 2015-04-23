
package main

import (
	"fmt"
	"github.com/negf_mtj/negf_mtj/cmplxSparse"
)

func main() {
	fmt.Println("Testing cmplxSparse package...\n\n");

	// Test identity matrix...
	testIdentity();

	// Test Tridiagonal matrix...
	testTriDiag(2);

	// Test Tridiagonal matrix...
	testTriDiag(4);

	// Test Tridiagonal matrix...
	testTriDiag(6);

	// Test Tridiagonal matrix...
	testHamTri(5);

	// Test matrix scaling...
	testScaleMatrix(3);

	// Test adding potential profile...
	testVoltage(6);

	// Test adding potential profile...
	testBarrier(6);

	// Test adding band split on left contact...
	testLeftSplit(3);

	// Test adding band split on right contact...
	testRightSplit(3);

	// Test scaling portions of matrix...
	testRangeScale(complex(1.5,0.0), 0, 6, 0, 10 );

	testRangeScale(complex(0.5,0.0), 13, 19, -2, 10 );

	testRangeScale(complex(4.5,0.0), 17, 17, 2, 10 );

	testAddModeProfile(3.0, 2, 1.0, 2, 2.0, 2, 3.0);

	testAddModeProfile(0.10, 2, 0.75, 4, 0.25, 2, 0.75);
}

func testIdentity() {
	fmt.Println("Initializing Sparse matrix structure\n");
	tmp := cmplxSparse.New();

	fmt.Println("Making sparse 5x5 identity matrix\n");
	matrixSize := int(5);
	cmplxSparse.MakeIdentity(matrixSize, tmp);

	fmt.Println("Accessing matrix elements (m,n):");

	for idx0 := 0; idx0 < matrixSize; idx0++ {
		for idx1 := 0; idx1 < matrixSize; idx1++ {
			test := cmplxSparse.AccessMatrix(idx0,idx1,tmp);
			fmt.Printf("%f  ", test);
		}
		fmt.Printf("\n");
	}
	fmt.Println("\n");
}

func testTriDiag( n int ) {
	fmt.Println("Initializing Sparse matrix structure\n");
	tmp := cmplxSparse.New();

	fmt.Println("Making sparse ",n,"x",n," tridiagonal matrix\n");
	matrixSize := int(n);
	cmplxSparse.MakeTriDiag(matrixSize, tmp);

	fmt.Println("Accessing matrix elements (m,n):");

	for idx0 := 0; idx0 < matrixSize; idx0++ {
		for idx1 := 0; idx1 < matrixSize; idx1++ {
			test := cmplxSparse.AccessMatrix(idx0,idx1,tmp);
			fmt.Printf("%f  ", test);
		}
		fmt.Printf("\n");
	}
	fmt.Println("\n");
}

func testHamTri( n int ) {
        fmt.Println("Initializing Sparse matrix structure\n");
        tmp := cmplxSparse.New();

        matrixSize := int(2*n);
        fmt.Println("Making sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix\n");
        cmplxSparse.MakeHamTriDiag(n, tmp);

        fmt.Println("Accessing matrix elements (m,n):");

        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,tmp);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }
	fmt.Println("\n");
}

func testScaleMatrix( n int ) {
        fmt.Println("Initializing Sparse matrix structure\n");
        tmp := cmplxSparse.New();

        matrixSize := int(2*n);
        fmt.Println("Making sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix\n");
        cmplxSparse.MakeHamTriDiag(n, tmp);

        fmt.Println("Scaling matrix by factor of 5+1*im\n");
	cmplxSparse.ScaleSparseMatrix(complex(5,1), tmp);
        fmt.Println("Accessing matrix elements (m,n):");

        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,tmp);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }

	fmt.Printf("\n");
        fmt.Println("Scaling matrix by factor of 1/5\n");
	cmplxSparse.ScaleSparseMatrix(complex128(0.2),tmp);
        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,tmp);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }
	fmt.Println("\n");

}

func testVoltage( n int ) {
        fmt.Println("Initializing Sparse matrix structure\n");
        tmp := cmplxSparse.New();

        matrixSize := int(2*n);
        fmt.Println("Making sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix for voltage profile test\n");
        cmplxSparse.MakeHamTriDiag(n, tmp);

        fmt.Println("Adding applied potential profile of 0.3 \n");
	tmp2 := cmplxSparse.AddVoltagePotential(1, n-4, 0.30, tmp);
        fmt.Println("Accessing matrix elements (m,n):");

        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,tmp2);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }
	fmt.Println("\n");

}

func testBarrier( n int ) {
        fmt.Println("Initializing Sparse matrix structure\n");
        tmp := cmplxSparse.New();

        matrixSize := int(2*n);
        fmt.Println("Making sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix for barrier test\n");
        cmplxSparse.MakeHamTriDiag(n, tmp);

        fmt.Println("Adding barrier profile of 5 \n");
	cmplxSparse.AddBarrierProfile(1, n-4, 5.0, tmp);
        fmt.Println("Accessing matrix elements (m,n):");

        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,tmp);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }
	fmt.Println("\n");

}

func testLeftSplit( n int ) {
        fmt.Println("Initializing Sparse matrix structure\n");
        tmp := cmplxSparse.New();

        matrixSize := int(2*(n+n+4));
        fmt.Println("Making sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix for band split test (left)\n");
        cmplxSparse.MakeHamTriDiag(n+n+4, tmp);

        fmt.Println("Adding splitting of 1.25 along +z direction\n");
	tmp2 := cmplxSparse.AddBandSplitLeftFM(0, 0, 0.1, 1.25, n, tmp);
        fmt.Println("Accessing matrix elements (m,n):");

        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,tmp2);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }
	fmt.Println("\n");

        fmt.Println("Adding splitting of 0.75 along (1,1,1) direction\n");
	tmp2 = cmplxSparse.AddBandSplitLeftFM(1.0, 1.0, 1.0, 0.75, n, tmp);
        fmt.Println("Accessing matrix elements (m,n):");

        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,tmp2);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }
	fmt.Println("\n");
}

func testRightSplit( n int ) {
        fmt.Println("Initializing Sparse matrix structure\n");
        tmp := cmplxSparse.New();

        matrixSize := int(2*(n+n+4));
        fmt.Println("Making sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix for band split test (right)\n");
        cmplxSparse.MakeHamTriDiag(n+n+4, tmp);

        fmt.Println("Adding splitting of 1.25 along +z direction\n");
	tmp2 := cmplxSparse.AddBandSplitRightFM(0, 0, 0.1, 1.25, n, tmp);
        fmt.Println("Accessing matrix elements (m,n):");

        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,tmp2);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }
	fmt.Println("\n");

        fmt.Println("Adding splitting of 0.75 along (1,1,1) direction\n");
	tmp2 = cmplxSparse.AddBandSplitRightFM(1.0, 1.0, 1.0, 0.75, n, tmp);
        fmt.Println("Accessing matrix elements (m,n):");

        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,tmp2);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }
	fmt.Println("\n");
}

func testRangeScale(A complex128, startIdx, endIdx, diagIdx, n int ) {
        fmt.Println("Initializing Sparse matrix structure\n");
        tmp := cmplxSparse.New();

        matrixSize := int(2*n);
        fmt.Println("Making sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix for range scale test\n");
        cmplxSparse.MakeHamTriDiag(n, tmp);

        fmt.Println("Scaling index ", startIdx, " to ", endIdx, " on the ", diagIdx," diagonal of matrix by ", A);
	cmplxSparse.ScaleRangeSparseMatrix(startIdx, endIdx, diagIdx, A, tmp);
        fmt.Println("Accessing matrix elements (m,n):");

        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,tmp);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }
	fmt.Println("\n");
}

func testAddModeProfile(E_mode float64, N_fmL int, m_fmL float64, N_ox int, m_ox float64, N_fmR int, m_fmR float64) {
        fmt.Println("Initializing Sparse matrix structure\n");
        tmp := cmplxSparse.New();

        matrixSize := int(2*(N_fmL+N_ox+N_fmR+2));
        fmt.Println("Making sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix for mode energy test\n");
        cmplxSparse.MakeHamTriDiag(matrixSize/2, tmp);

        fmt.Println("Printing sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix before action\n");
	cmplxSparse.PrintSparseMatrix(tmp);

	cmplxSparse.AddModeEnergy(E_mode, N_fmL, m_fmL, N_ox, m_ox, N_fmR, m_fmR, tmp);

        fmt.Println("Printing sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix after action\n");
	cmplxSparse.PrintSparseMatrix(tmp);
}
