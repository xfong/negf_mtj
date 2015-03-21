
package main

import (
	"fmt"
	"cmplxSparse"
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
}

func testScaleMatrix( n int ) {
        fmt.Println("Initializing Sparse matrix structure\n");
        tmp := cmplxSparse.New();

        matrixSize := int(2*n);
        fmt.Println("Making sparse ",matrixSize,"x",matrixSize," Hamiltonian tridiagonal matrix\n");
        cmplxSparse.MakeHamTriDiag(n, tmp);

        fmt.Println("Scaling matrix by factor of 5+1*im\n");
	cmplxSparse.ScaleSparseMatrixIP(complex(5,1), tmp);
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
	tmp2 := cmplxSparse.ScaleSparseMatrix(complex128(0.2),tmp);
        for idx0 := 0; idx0 < matrixSize; idx0++ {
                for idx1 := 0; idx1 < matrixSize; idx1++ {
                        test := cmplxSparse.AccessMatrix(idx0,idx1,tmp2);
                        fmt.Printf("%f  ", test);
                }
                fmt.Printf("\n");
        }

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
}