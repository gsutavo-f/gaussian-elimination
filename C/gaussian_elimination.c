//GitHub: HenriqueIni
//www.blogcyberini.com

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * Algoritmo para resolução de sistemas lineares via eliminação de Gauss
 * Complexidade no tempo: O(n^3)
 * A é a matriz dos coeficientes do sistema
 * b é a matriz dos coeficientes dos termos independentes
 * n é a ordem do sistema
 * X é o vetor onde a solução será armazenada
 * Forma do sistema (matricial): Ax = b  
 */
void gaussSolver(int n, double A[n][n], double b[n], double X[n]) {
    int i, j, k, l, m;
    //ETAPA DE ESCALONAMENTO
    for (k = 0; k < n - 1; k++) {
        double max = fabs(A[k][k]);
        int maxIndex = k;
        //procura o maior k-ésimo coeficiente em módulo
        for (i = k + 1; i < n; i++) {
            if (max < fabs(A[i][k])) {
                max = fabs(A[i][k]);
                maxIndex = i;
            }
        }
        if (maxIndex != k) {
            /*
             troca a equação k pela equação com o
             maior k-ésimo coeficiente em módulo
             */
            for (j = 0; j < n; j++) {
                double temp = A[k][j];
                A[k][j] = A[maxIndex][j];
                A[maxIndex][j] = temp;
            }
            double temp = b[k];
            b[k] = b[maxIndex];
            b[maxIndex] = temp;
        }
        //Se A[k][k] é zero, então a matriz dos coeficiente é singular
        //det A = 0
        if (A[k][k] == 0) {
            printf("A matriz dos coeficientes é singular\n");
            return;
        } else {
            //realiza o escalonamento
            for (m = k + 1; m < n; m++) {
                double F = -A[m][k] / A[k][k];
                A[m][k] = 0; //evita uma iteração
                b[m] = b[m] + F * b[k];
                for (l = k + 1; l < n; l++) {
                    A[m][l] = A[m][l] + F * A[k][l];
                }
            }
        }
    }
    //ETAPA DE RESOLUÇÃO DO SISTEMA
    for (i = n - 1; i >= 0; i--) {
        X[i] = b[i];
        for (j = i + 1; j < n; j++) {
            X[i] = X[i] - X[j] * A[i][j];
        }
        X[i] = X[i] / A[i][i];
    }
}
//Código de testes
int main(){
    int n = 3;
    double A[3][3] = {{2, 1, -1}, {1, 2, 1}, {1, 1, 1}};
    double b[3] = {-3, 3, 2};
    double x[3];
    gaussSolver(n, A, b, x);
    printf("x1 = %f\nx2 = %f\nx3 = %f\n", x[0], x[1], x[2]);
    
    double A2[3][3] = {{1, 1, 1}, {-2, 1, 1}, {1, 3, 1}};
    double b2[3] = {2, 5, 4};
    gaussSolver(n, A2, b2, x);
    printf("\nx1 = %f\nx2 = %f\nx3 = %f\n", x[0], x[1], x[2]);
    
    double A3[3][3] = {{3, 2, -1}, {2, -2, 4}, {-1, 0.5, -1}};
    double b3[3] = {1, -2, 0};
    gaussSolver(n, A3, b3, x);
    printf("\nx1 = %f\nx2 = %f\nx3 = %f\n", x[0], x[1], x[2]);
    
    n=2;
    double A4[2][2] = {{2,3}, {4,9}};
    double b4[2] = {6,15};
    gaussSolver(n, A4, b4, x);
    printf("\nx1 = %f\nx2 = %f\n", x[0], x[1]);
    
    n = 4;
    double A5[4][4] = {{4, 1, 2, -3}, {-3, 3, -1, 4}, {-1, 2, 5, 1}, {5, 4, 3, -1}};
    double b5[4] = {-16, 20, -4, -10};
    gaussSolver(n, A5, b5, x);
    printf("\nx1 = %f\nx2 = %f\nx3 = %f\nx4 = %f\n", x[0], x[1], x[2], x[3]);

    n = 5;
    double A6[5][5] = {{4, 1, 2, -3, 5}, {-3, 3, -1, 4, -2}, {-1, 2, 5, 1, 3}, {5, 4, 3, -1, 2}, {1, -2, 3, -4, 5}};
    double b6[5] = {-16, 20, -4, -10, 3};
    gaussSolver(n, A6, b6, x);
    printf("\nx1 = %f\nx2 = %f\nx3 = %f\nx4 = %f\nx5 = %f\n", x[0], x[1], x[2], x[3], x[4]);
    return 0;
}