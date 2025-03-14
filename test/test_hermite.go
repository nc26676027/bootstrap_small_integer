package main

import (
	"fmt"
	"math"
	"math/cmplx"
)

func main() {
	t := 256.0
	nodesCount := int(t)
	testPoints := []float64{0, 1.001, 2.01, 3.02, 4., 6.001}

	coeffs := hermiteInterpolation(nodesCount)
	fmt.Println("Interpolation polynomial degree: ", )
	testInterpolation(coeffs, testPoints, t)
}

func generateRoots(nodesCount int) []complex128 {
	roots := make([]complex128, nodesCount)
	for k := 0; k < nodesCount; k++ {
		angle := 2 * math.Pi * float64(k) / float64(nodesCount)
		roots[k] = cmplx.Rect(1, angle)
	}
	return roots
}

func buildHermiteMatrix(nodesCount int, roots []complex128) ([][]complex128, []complex128) {
	matrixSize := 2 * nodesCount
	matrix := make([][]complex128, matrixSize)
	for i := range matrix {
		matrix[i] = make([]complex128, matrixSize)
	}

	rhs := make([]complex128, matrixSize)

	// Generate fValues (identity function)
	fValues := make([]complex128, nodesCount)
	for i := range fValues {
		fValues[i] = complex(float64(i), 0)
	}

	// All derivatives are zero
	fDerivatives := make([]complex128, nodesCount)

	// Populate matrix
	for i := 0; i < nodesCount; i++ {
		for j := 0; j < matrixSize; j++ {
			matrix[i][j] = complexIntPower(roots[i], j)
			if j > 0 {
				matrix[nodesCount+i][j] = complex(float64(j), 0) * complexIntPower(roots[i], j-1)
			} else {
				matrix[nodesCount+i][j] = complex(0, 0)
			}
		}
	}

	// Populate RHS vector
	for i := 0; i < nodesCount; i++ {
		rhs[i] = fValues[i]
		rhs[nodesCount+i] = fDerivatives[i]
	}
	return matrix, rhs
}

func hermiteInterpolation(nodesCount int) []complex128 {
	roots := generateRoots(nodesCount)
	matrix, rhs := buildHermiteMatrix(nodesCount, roots)
	coeffs := solveMatrix(matrix, rhs)
	return coeffs
}

func testInterpolation(coeffs []complex128, testPoints []float64, t float64) {
	for _, x := range testPoints {
		z := expComplex(x, t)
		interpolated := evalPolynomial(coeffs, z)
		expected := complex(x, 0)
		errorVal := cmplx.Abs(interpolated - expected)
		fmt.Printf("Test point x=%.4f and Interpolation result: %.4f, error: %.2e\n", x, interpolated, errorVal)
	}
}

func complexIntPower(base complex128, exp int) complex128 {
	if exp == 0 {
		return complex(1, 0)
	}
	result := base
	for i := 1; i < exp; i++ {
		result *= base
	}
	return result
}

func solveMatrix(A [][]complex128, b []complex128) []complex128 {
	n := len(A)
	x := make([]complex128, n)
	fmt.Printf("Matrix dimension: %d x %d\n", n, n)
	for i := 0; i < n; i++ {
		pivot := A[i][i]
		if pivot == 0 {
			fmt.Printf("Warning: Line %d main element is zero, invalid solution\n", i)
			continue
		}
		for j := i + 1; j < n; j++ {
			mul := A[j][i] / pivot
			for k := i; k < n; k++ {
				A[j][k] -= mul * A[i][k]
			}
			b[j] -= mul * b[i]
		}
	}
	for i := n - 1; i >= 0; i-- {
		sum := complex(0, 0)
		for j := i + 1; j < n; j++ {
			sum += x[j] * A[i][j]
		}
		if A[i][i] == 0 {
			fmt.Printf("Fail: Line %d main element is zero, cannot solve\n", i)
			x[i] = complex(0, 0)
			continue
		}
		x[i] = (b[i] - sum) / A[i][i]
	}
	return x
}

func expComplex(x float64, t float64) complex128 {
	angle := 2 * math.Pi * x / t
	return cmplx.Exp(complex(0, angle))
}

func evalPolynomial(coeffs []complex128, z complex128) complex128 {
	result := complex(0, 0)
	zPower := complex(1, 0)
	for _, coeff := range coeffs {
		result += coeff * zPower
		zPower *= z
	}
	return result
}
