
package main

import (
	"math"
	"fmt"
	"github.com/negf_mtj/negf_mtj/vecQuad"
)

func main() {
	fmt.Println("Testing vecQuad package...\n\n");

	n := float64(0.0);

	// Test integrate exp(-x) over [0, inf)
	testExp0(n)
	fmt.Printf("\n");

	// Test integrate exp(-x) over [1, inf)
	n = 1.0;
	testExp0(n)
	fmt.Println("Actual result should be =", math.Exp(-1.0));
	fmt.Printf("\n");

	// Test vectorized integration over [1, inf)
	testExp1(n)
	fmt.Println("Actual result should be =", math.Exp(-n), " and ", float64(0.5)*math.Exp(-2.0*n) );
	fmt.Printf("\n");

	// Test vectorized integration over [1, inf)
	n = 0.0;
	testExp2(n)
	fmt.Println("Actual result should be =", 1.0, ", ", float64(0.5), ", ", float64(0.5), " and ", float64(1.0/3.0));
	fmt.Printf("\n");
}

func SimpleExp( n float64 ) *[]float64 {
	test := make([]float64, 1);
	test[0] = math.Exp(-1.0*n);
	return &test;
}

func VectorizedFunc0( n float64 ) *[]float64 {
	test := make([]float64, 2);
	test[0] = math.Exp(-1.0*n);
	test[1] = math.Exp(-2.0*n);
	return &test;
}

func VectorizedFunc1( n float64 ) *[]float64 {
	test := make([]float64, 4);
	test[0] = math.Exp(-1.0*n);
	test[1] = math.Exp(-2.0*n);
	test[2] = math.Exp(-0.5*n*n)/math.Sqrt(2.0*math.Pi);
	test[3] = math.Exp(-3.0*n);
	return &test;
}

func testExp0(n float64) {
	answer, errbnd := vecQuad.IntegralCalc2Inf(SimpleExp, n, 1);
	fmt.Println("Integral of 1/Exp(x) from 0 to inf =", answer);
	fmt.Println("Error of integral =", errbnd);
}

func testExp1(n float64) {
	answer, errbnd := vecQuad.IntegralCalc2Inf(VectorizedFunc0, n, 2);
	fmt.Println("Integral of 1/Exp(x) from 0 to inf =", answer[0]);
	fmt.Println("Error of first integral =", errbnd[0]);
	fmt.Println("Integral of 1/Exp(2*x) from 0 to inf =", answer[1]);
	fmt.Println("Error of second integral =", errbnd[1]);
}

func testExp2(n float64) {
	answer, errbnd := vecQuad.IntegralCalc2Inf(VectorizedFunc1, n, 4);
	fmt.Println("Integral of 1/Exp(x) from 0 to inf =", answer[0]);
	fmt.Println("Error of first integral =", errbnd[0]);
	fmt.Println("Integral of 1/Exp(2*x) from 0 to inf =", answer[1]);
	fmt.Println("Error of second integral =", errbnd[1]);
	fmt.Println("Integral of 1/(Exp(0.5*x*x) * Sqrt(2.0*Pi)) from 0 to inf =", answer[2]);
	fmt.Println("Error of third integral =", errbnd[2]);
	fmt.Println("Integral of 1/Exp(3*x) from 0 to inf =", answer[3]);
	fmt.Println("Error of fourth integral =", errbnd[3]);
}
