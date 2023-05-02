package main

import (
	"fmt"
	"math"
	"time"
)

func gaussSolver(n int, A [][]float64, b []float64) []float64 {
	var i, j, k, l, m int
	x := make([]float64, n)
	//ETAPA DE ESCALONAMENTO
	for k = 0; k < n-1; k++ {
		max := math.Abs(A[k][k])
		indiceMax := k
		for i = k + 1; i < n; i++ {
			if max < math.Abs(A[i][k]) {
				max = math.Abs(A[i][k])
				indiceMax = i
			}
		}
		if indiceMax != k {

			for j = 0; j < n; j++ {
				temp := A[k][j]
				A[k][j] = A[indiceMax][j]
				A[indiceMax][j] = temp
			}
			temp := b[k]
			b[k] = b[indiceMax]
			b[indiceMax] = temp
		}
		if A[k][k] == 0 {
			fmt.Println("A matriz dos coeficientes é singular")
			return x
		} else {
			for m = k + 1; m < n; m++ {
				F := -A[m][k] / A[k][k]
				A[m][k] = 0
				b[m] = b[m] + F*b[k]

				for l = k + 1; l < n; l++ {
					A[m][l] = A[m][l] + F*A[k][l]
				}
			}
		}
	}
	for i = n - 1; i >= 0; i-- {
		x[i] = b[i]
		for j = i + 1; j < n; j++ {
			x[i] = x[i] - x[j]*A[i][j]
		}
		x[i] = x[i] / A[i][i]
	}
	return x
}
func gauss_solver_executionTime(n int, A [][]float64, b []float64) (tempoDecorrido int64) {
	inicio := time.Now().UnixNano() / int64(time.Millisecond)

	for i := 0; i < 100000; i++ {
		gaussSolver(n, A, b)
	}
	fim := time.Now().UnixNano() / int64(time.Millisecond)
	tempoDecorrido = fim - inicio
	fmt.Printf("Time elapsed in gauss_solver() is %v us\n\n", tempoDecorrido)
	return tempoDecorrido
}
func main() {
	var tempo_decorrido int64 = 0
	n := 3
	A := [][]float64{{2, 1, -1}, {1, 2, 1}, {1, 1, 1}}
	b := []float64{-3, 3, 2}
	tempo_decorrido = tempo_decorrido + gauss_solver_executionTime(n, A, b)

	var A2 = [][]float64{{1, 1, 1}, {-2, 1, 1}, {1, 3, 1}}
	var b2 = []float64{2, 5, 4}
	tempo_decorrido = tempo_decorrido + gauss_solver_executionTime(n, A2, b2)

	A3 := [][]float64{{3, 2, -1}, {2, -2, 4}, {-1, 0.5, -1}}
	b3 := []float64{1, -2, 0}
	tempo_decorrido = tempo_decorrido + gauss_solver_executionTime(n, A3, b3)

	n = 2
	A4 := [][]float64{{2, 3}, {4, 9}}
	b4 := []float64{6, 15}
	tempo_decorrido = tempo_decorrido + gauss_solver_executionTime(n, A4, b4)

	n = 4
	A5 := [][]float64{{4, 1, 2, -3}, {-3, 3, -1, 4}, {-1, 2, 5, 1}, {5, 4, 3, -1}}
	b5 := []float64{-16, 20, -4, -10}
	tempo_decorrido = tempo_decorrido + gauss_solver_executionTime(n, A5, b5)
	n = 5
	A6 := [][]float64{{4, 1, 2, -3, 5}, {-3, 3, -1, 4, -2}, {-1, 2, 5, 1, 3}, {5, 4, 3, -1, 2}, {1, -2, 3, -4, 5}}
	b6 := []float64{-16, 20, -4, -10, 3}
	tempo_decorrido = tempo_decorrido + gauss_solver_executionTime(n, A6, b6)

	fmt.Printf("Total time elapsed in program is %v us", tempo_decorrido)
}
