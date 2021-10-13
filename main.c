#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linear_func.h"

int main() {
    int n, is_solved,option = 0;
    double **A, *b, *x,*r, tol = 0;
    while(option != 3){
        printf("---------Matrix Solver---------\n");
        printf("[1]. Matriz triangular superior\n");
        printf("[2]. Matriz triangular inferior\n");
        printf("[3]. Salir\n");
        printf("-------------------------------\n");
        scanf("%d",&option);
        switch (option) {
            case 1:
                printf("Ingrese la dimensión de la matriz:\n");
                scanf("%d", &n);
                A = (double **) malloc(sizeof(double) * n);
                b = (double *) malloc(sizeof(double) * n);
                x = (double *) malloc(sizeof(double) * n);
                r = (double *) malloc(sizeof(double) * n);
                if (A == NULL || b == NULL || x == NULL) {
                    printf("No hay suficiente memoria\n");
                    exit(2);
                }
                for (int i = 0; i < n; i++) {
                    A[i] = (double *) malloc(n * sizeof(double));
                    if (A[i] == NULL) {
                        printf("No hay suficiente memoria\n");
                        exit(2);
                    }
                }
                printf("Ingrese los elementos de la matriz:\n");
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        scanf("%lf", &A[i][j]);
                    }
                }
                if (checkupptri(A, n)) {
                    printf("La matriz es triangular superior\n");
                } else {
                    printf("\033[0;31m");
                    printf("La matriz no es triangular superior\n");
                    exit(2);
                }
                printf("Ingrese los elementos del vector B:\n");
                for (int i = 0; i < n; ++i) {
                    scanf("%lf", &b[i]);
                }
                printf("Ingrese la tolerancia aceptada:\n");
                scanf("%lf",&tol);
                printf("Matriz ingresada:\n");
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        printf("%lf\t", A[i][j]);
                    }
                    printf("|\t %lf",b[i]);
                    printf("\n");
                }
                is_solved = resoltrisup(n, A, b, x, tol);
                if(is_solved == 0){
                    printf("El sistema se ha resuelto: Status %d\n",is_solved);
                    printf("\n");
                    printf("-----------Soluciones----------\n");
                    for (int i = 0; i < n; ++i) {
                        printf("X%d: %lf\n",i+1, x[i]);
                    }
                    printf("-------------------------------\n");
                    double *Ax = prodMatVect(A,x,n);
                    printf("Comprobación:\n");
                    for (int i = 0; i < n; ++i) {
                        printf("%lf(Ax) = %lf(b)\n",Ax[i],b[i]);
                    }
                    for (int i = 0; i < n; ++i) {
                        for (int j = 0; j < n; ++j) {
                            r[i] += sqrt(pow(((A[i][j]*x[i])-b[i]),2));
                        }
                    }
                    printf("\n");
                    printf("Vector Residuo:\n");
                    for (int i = 0; i < n; ++i) {
                        printf("%lf\n",r[i]);
                    }
                    printf("\n");
                }else{
                    printf("\033[0;31m");
                    printf("El sistema no se puede resolver\n");
                    printf("\033[0;37m");
                }
                break;
            case 2:
                printf("Ingrese la dimensión de la matriz:\n");
                scanf("%d", &n);
                A = (double **) malloc(sizeof(double) * n);
                b = (double *) malloc(sizeof(double) * n);
                x = (double *) malloc(sizeof(double) * n);
                r = (double *) malloc(sizeof(double) * n);
                if (A == NULL || b == NULL || x == NULL) {
                    printf("No hay suficiente memoria\n");
                    exit(2);
                }
                for (int i = 0; i < n; i++) {
                    A[i] = (double *) malloc(n * sizeof(double));
                    if (A[i] == NULL) {
                        printf("No hay suficiente memoria\n");
                        exit(2);
                    }
                }
                printf("Ingrese los elementos de la matriz:\n");
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        scanf("%lf", &A[i][j]);
                    }
                }
                if (checklowertri(A, n)) {
                    printf("La matriz es triangular inferior\n");
                } else {
                    printf("\033[0;31m");
                    printf("La matriz no es triangular inferior\n");
                    exit(2);
                }
                printf("Ingrese los elementos del vector B:\n");
                for (int i = 0; i < n; ++i) {
                    scanf("%lf", &b[i]);
                }
                printf("Ingrese la tolerancia aceptada:\n");
                scanf("%lf",&tol);
                printf("Matriz ingresada:\n");
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        printf("%lf\t", A[i][j]);
                    }
                    printf("|\t %lf",b[i]);
                    printf("\n");
                }
                is_solved = resoltriinf(n, A, b, x, tol);
                if(is_solved == 0){
                    printf("El sistema se ha resuelto: Status %d\n",is_solved);
                    printf("\n");
                    printf("-----------Soluciones----------\n");
                    for (int i = 0; i < n; ++i) {
                        printf("X%d: %lf\n",i+1, x[i]);
                    }
                    printf("-------------------------------\n");
                    double *Ax = prodMatVect(A,x,n);
                    printf("Comprobación:\n");
                    for (int i = 0; i < n; ++i) {
                        printf("%lf(Ax) = %lf(b)\n",Ax[i],b[i]);
                    }
                    for (int i = 0; i < n; ++i) {
                        for (int j = 0; j < n; ++j) {
                            r[i] += sqrt(pow(((A[i][j]*x[i])-b[i]),2));
                        }
                    }
                    printf("\n");
                    printf("Vector Residuo:\n");
                    for (int i = 0; i < n; ++i) {
                        printf("%lf\n",r[i]);
                    }
                    printf("\n");
                }else{
                    printf("\033[0;31m");
                    printf("El sistema no se puede resolver\n");
                    printf("\033[0;37m");
                }
                break;

        }
    }

    free(A);
    free(b);
    free(x);
    free(r);
    return 0;
}

