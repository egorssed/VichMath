#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void Gauss(double *Matr, int Nx, double *free, double *sol) {
	//Метод Гаусса ----------------------------------------------------
	double *BufMatr = (double *)calloc(Nx*Nx, sizeof(double));
	double *MatrRev = (double *)calloc(Nx*Nx, sizeof(double));
	for (int i = 0; i < Nx; i++) {
		MatrRev[i*Nx + i] = 1;
	}
	int N = Nx;
	//проверка матрицы
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%lf ", Matr[i*N + j]);
		}
		printf("%lf", free[i]);
		printf("\n");
	}
	printf("\n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%lf ", MatrRev[i*N + j]);
		}
		printf("%lf", free[i]);
		printf("\n");
	}
	//прямой ход Гаусса
	for (int q = 0; q < N; q++) {
		for (int i = q; i < N - 1; i++) {
			double a = Matr[(q)*N + q];
			double b = Matr[(i + 1)*N + q];
			for (int j = q; j < N; j++) {
				Matr[(i + 1)*N + j] += Matr[(q)*N + j] * (-b / a);
				MatrRev[(i + 1)*N + j] += MatrRev[(q)*N + j] * (-b / a);
			}
			free[i + 1] += free[q] * (-b / a);
		}
	}
	//проверка матрицы
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%lf ", Matr[i*N + j]);
		}
		printf("%lf", free[i]);
		printf("\n");
	}
	printf("\n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%lf ", MatrRev[i*N + j]);
		}
		printf("%lf", free[i]);
		printf("\n");
	}
	printf("\n");

	//обратная матрица-----------------------------
	//буферная матрица
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			BufMatr[i*N + j] = Matr[(i)*N + j];
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%lf ", BufMatr[i*N + j]);
		}
		printf("\n");
	}
	printf("\n");
	//обратный ход в приведениее А к единичному
	for (int q = N - 1; q >= 0; q--) {
		for (int i = q; i >= 1; i--) {
			double a = BufMatr[(q)*N + q];
			double b = BufMatr[(i - 1)*N + q];
			for (int j = q; j < N; j++) {
				BufMatr[(i - 1)*N + j] += BufMatr[(q)*N + j] * (-b / a);
				MatrRev[(i - 1)*N + j] += MatrRev[(q)*N + j] * (-b / a);
			}
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%lf ", BufMatr[i*N + j]);
		}
		printf("\n");
	}
	printf("\n");
	//нормировка А на 1
	for (int q = 0; q < N; q++) {
		for (int i = q; i < N; i++) {
			for (int j = q; j < N; j++) {
				BufMatr[i*N + j] = BufMatr[i*N + j] / BufMatr[i*N + i];
				MatrRev[i*N + j] = MatrRev[i*N + j] / BufMatr[i*N + i];
			}
		}
	}
	printf("\n");
	//-------------------------------------------

	printf("\n");
	//обратный ход Гаусса
	for (int q = N - 1; q >= 0; q--) {
		sol[q] = free[q];
		for (int d = q + 1; d < N; d++) {
			sol[q] += -sol[d] * Matr[q*N + d];
		}
		sol[q] = sol[q] / Matr[q*N + q];
	}
	//решения Гаусс
	for (int q = 0; q < N; q++) {
		printf("solution %lf  \n", sol[q]);
	}
}
int main() {
	int Nx = 3;
	//double *A = (double *)calloc(Nx*Nx, sizeof(double));
	//double *B = (double *)calloc(Nx, sizeof(double));
	double *S = (double *)calloc(Nx, sizeof(double));
	double A[9] = { 1,2,4,4,7,6,7,8,9 };
	double B[3] = { 1,2,23 };
	printf("\n");
	Gauss(A, Nx, B, S);
	for (int i = 0; i < Nx; i++) printf("%lf ", S[i]);
	scanf_s("%d", &Nx);
	return 0;
}