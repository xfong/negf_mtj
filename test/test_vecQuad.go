
package main

import (
	"math"
	"fmt"
	"github.com/negf_mtj/negf_mtj/vecQuad"
)

func main() {
	fmt.Println("Testing vecQuad package...\n\n");

	// Test integrate exp(-x) over [0, inf)
	testExp0(0.0)
	fmt.Printf("\n");

	// Test integrate exp(-x) over [1, inf)
	testExp0(1.0)
	fmt.Println("Actual result should be =", math.Exp(-1.0));
	fmt.Printf("\n");

	// Test vectorized integration over [1, inf)
	testExp1(1.0)
	fmt.Println("Actual result should be =", math.Exp(-1.0), " and ", float64(0.5)*math.Exp(-1.0) );
	fmt.Printf("\n");
}

func SimpleExp( n float64 ) *[]float64 {
	test := make([]float64, 1);
	test[0] = math.Exp(-1.0*n);
	return &test;
}

func VectorizedFunc( n float64 ) *[]float64 {
	test := make([]float64, 2);
	test[0] = math.Exp(-1.0*n);
	test[1] = math.Exp(-2.0*n);
	return &test;
}

func testExp0(n float64) {
	answer, errbnd := vecQuad.IntegralCalc2Inf(SimpleExp, n, 1);
	fmt.Println("Integral of 1/Exp(x) from 0 to inf =", answer);
	fmt.Println("Error of integral =", errbnd);
}

func testExp1(n float64) {
	answer, errbnd := vecQuad.IntegralCalc2Inf(VectorizedFunc, n, 2);
	fmt.Println("Integral of 1/Exp(x) from 0 to inf =", answer[0]);
	fmt.Println("Error of first integral =", errbnd[0]);
	fmt.Println("Integral of 1/Exp(2*x) from 0 to inf =", answer[1]);
	fmt.Println("Error of second integral =", errbnd[1]);
}
