
package main

import (
	"math"
	"fmt"
	"github.com/negf_mtj/negf_mtj/vecQuad"
)

func main() {
	fmt.Println("Testing vecQuad package...\n\n");

	n, m := float64(0.0), float64(0.0);

	// Test integrate exp(-x) over [0, inf)
	testExp0(n)
	fmt.Printf("\n");

	// Test integrate exp(-x) over [1, inf)
	n = 1.0;
	testExp0(n)
	fmt.Printf("\n");

	// Test vectorized integration over [1, inf)
	testExp1(n)
	fmt.Printf("\n");

	// Test vectorized integration over [1, inf)
	n = 0.0;
	testExp2(n)
	fmt.Printf("\n");

	n, m = 0.0, 1.0;
	testExp3(n, m);
	fmt.Printf("\n");

	testExp4(n, m);
	fmt.Printf("\n");

	n, m = 2.0, 3.0;
	testExp5(n, m);
	fmt.Printf("\n");

	n, m = 1.5, 3.0;
	testExp6(n, m);
	fmt.Printf("\n");

}

func SimpleExp( n float64 ) *[]float64 {
	test := make([]float64, 1);
	test[0] = math.Exp(-1.0*n);
	return &test;
}

func SimpleQuadr( n float64 ) *[]float64 {
	test := make([]float64, 1);
	test[0] = n*n;
	return &test
}

func SimpleLinear( n float64 ) *[]float64 {
	test := make([]float64, 1);
	test[0] = 2.0*n;
	return &test
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

func VectorizedFunc2( n float64 ) *[]float64 {
	test := make([]float64, 2);
	test[0] = n*n;
	test[1] = 3.0*n;
	return &test;
}

func VectorizedFunc3( n float64 ) *[]float64 {
	test := make([]float64, 4);
	test[0] = math.Exp(-1.0*n);
	test[1] = 2.0*n;
	test[2] = n*n*n;
	test[3] = math.Exp(-2.0*n);
	return &test;
}

func testExp0(n float64) {
	nn := make([]float64, 1);
	nn[0] = n;
	answer, errbnd := vecQuad.IntegralCalc(SimpleExp, &nn, 1);
	fmt.Println("Integral of 1/Exp(x) from 0 to inf =", answer);
	fmt.Println("Actual result should be =", math.Exp(-n));
	fmt.Println("Error of integral =", errbnd);
}

func testExp1(n float64) {
	nn := make([]float64, 1);
	nn[0] = n;
	answer, errbnd := vecQuad.IntegralCalc(VectorizedFunc0, &nn, 2);
	fmt.Println("Integral of 1/Exp(x) from 0 to inf =", answer[0]);
	fmt.Println("Actual result should be =", math.Exp(-n));
	fmt.Println("Error of first integral =", errbnd[0]);
	fmt.Println("Integral of 1/Exp(2*x) from 0 to inf =", answer[1]);
	fmt.Println("Actual result should be =", float64(0.5)*math.Exp(-2.0*n) );
	fmt.Println("Error of second integral =", errbnd[1]);
}

func testExp2(n float64) {
	nn := make([]float64, 1);
	nn[0] = n;
	answer, errbnd := vecQuad.IntegralCalc(VectorizedFunc1, &nn, 4);
	fmt.Println("Integral of 1/Exp(x) from 0 to inf =", answer[0]);
	fmt.Println("Actual result should be =", math.Exp(-1.0*n));
	fmt.Println("Error of first integral =", errbnd[0]);
	fmt.Println("Integral of 1/Exp(2*x) from 0 to inf =", answer[1]);
	fmt.Println("Actual result should be =", 0.5*math.Exp(-2.0*n));
	fmt.Println("Error of second integral =", errbnd[1]);
	fmt.Println("Integral of 1/(Exp(0.5*x*x) * Sqrt(2.0*Pi)) from 0 to inf =", answer[2]);
	fmt.Println("Actual result should be =", float64(0.5));
	fmt.Println("Error of third integral =", errbnd[2]);
	fmt.Println("Integral of 1/Exp(3*x) from 0 to inf =", answer[3]);
	fmt.Println("Actual result should be =", float64(1.0/3.0)*math.Exp(-3.0*n));
	fmt.Println("Error of fourth integral =", errbnd[3]);
}

func testExp3(n, m float64) {
	LimitsList := make([]float64, 2);
	LimitsList[0], LimitsList[1] = n, m;
	answer, errbnd := vecQuad.IntegralCalc(SimpleLinear, &LimitsList, 1);
	fmt.Println("Integral of x from ", n, " to ", m, " =", answer[0]);
	fmt.Println("Actual result should be =", (m*m - n*n));
	fmt.Println("Error of integral =", errbnd);
}

func testExp4(n, m float64) {
	LimitsList := make([]float64, 2);
	LimitsList[0], LimitsList[1] = n, m;
	answer, errbnd := vecQuad.IntegralCalc(SimpleQuadr, &LimitsList, 1);
	fmt.Println("Integral of x^2 from ", n, " to ", m, " =", answer[0]);
	fmt.Println("Actual result should be =", float64(1.0/3.0) * (m*m*m - n*n*n));
	fmt.Println("Error of integral =", errbnd);
}

func testExp5(n, m float64) {
	LimitsList := make([]float64, 2);
	LimitsList[0], LimitsList[1] = n, m;
	answer, errbnd := vecQuad.IntegralCalc(VectorizedFunc2, &LimitsList, 2);
	fmt.Println("Integral of x^2 from ", n, " to ", m, " =", answer[0]);
	fmt.Println("Actual result should be =", float64(1.0/3.0)*(m*m*m - n*n*n));
	fmt.Println("Error of first integral =", errbnd[0]);
	fmt.Println("Integral of 3*x from ", n, " to ", m, " =", answer[1]);
	fmt.Println("Actual result should be =", float64(1.5) * (m*m - n*n));
	fmt.Println("Error of second integral =", errbnd[1]);
}

func testExp6(n, m float64) {
	LimitsList := make([]float64, 2);
	LimitsList[0], LimitsList[1] = n, m;
	answer, errbnd := vecQuad.IntegralCalc(VectorizedFunc3, &LimitsList, 4);
	fmt.Println("Integral of 1/Exp(x) from ", n, " to ", m, " =", answer[0]);
	fmt.Println("Actual result should be =", math.Exp(-1.0*n) - math.Exp(-1.0*m));
	fmt.Println("Error of first integral =", errbnd[0]);
	fmt.Println("Integral of 2*x from ", n, " to ", m, " =", answer[1]);
	fmt.Println("Actual result should be =", (m*m - n*n));
	fmt.Println("Error of second integral =", errbnd[1]);
	fmt.Println("Integral of x^3 from ", n, " to ", m, " =", answer[2]);
	fmt.Println("Actual result should be =", float64(0.25)*(m*m*m*m - n*n*n*n));
	fmt.Println("Error of third integral =", errbnd[2]);
	fmt.Println("Integral of 1/Exp(2*x) from ", n, " to ", m, " =", answer[3]);
	fmt.Println("Actual result should be =", float64(0.5)*(math.Exp(-2.0*n) - math.Exp(-2.0*m)));
	fmt.Println("Error of fourth integral =", errbnd[3]);
}

