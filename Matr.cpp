#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
//матрицу создаем так что второй индекс - столбцы , первый - строки

int main() {
	int N = 10;
	double Matr[10][10];
	double BufMatr[10][10];
	double MatrRev[10][10];
	double free[10];
	double sol[10];
	double solZ[10];
	double solZ1[10];
	//заполнение матрицы
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			MatrRev[i][j] = 0;
			if ((i != j) && (abs(i - j) < 5)) {
				Matr[i][j] = 1;
			}
			else {
				Matr[i][j] = 0;
			}
			if (i == j) {
				Matr[i][j] = 10;
				MatrRev[i][j] = 1;
			}
		}
	}
	for (int q = 0; q < N; q++) {
		free[q] = q + 1;
	}
	 //проверка матрицы
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%lf ", Matr[i][j]);
		}
		printf("%lf",free[i]);
		printf("\n");
	}
	printf("\n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%lf ", MatrRev[i][j]);
		}
		printf("%lf", free[i]);
		printf("\n");
	}
	printf("\n");
	//Метод Гаусса ----------------------------------------------------
	//прямой ход Гаусса
	for (int q = 0; q < N; q++) {
		for (int i = q; i < N-1; i++) {
			double a = Matr[q][q]; 
			double b = Matr[i + 1][q]; 
			for (int j = q; j < N; j++) {
				Matr[i+1][j] += Matr[q][j] * (-b / a);
				MatrRev[i+1][j]+= MatrRev[q][j] * (-b / a);
			}
			free[i+1] += free[q] * (-b / a);
		}
	}
	//проверка матрицы
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%lf ", Matr[i][j]);
		}
		printf("%lf", free[i]);
		printf("\n");
	}
	printf("\n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%lf ", MatrRev[i][j]);
		}
		printf("%lf", free[i]);
		printf("\n");
	}

	//обратная матрица-----------------------------
	//буферная матрица
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			BufMatr[i][j] = Matr[i][j];
		}
	}
	//обратный ход в приведениее А к единичному
	for (int q = N-1; q>=0; q--) {
		for (int i = q; i >= 1; i--) {
			double a = BufMatr[q][q];
			double b = BufMatr[i - 1][q];
			for (int j = q; j < N; j++) {
				BufMatr[i - 1][j] += BufMatr[q][j] * (-b / a);
				MatrRev[i - 1][j] += MatrRev[q][j] * (-b / a);
			}
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			//printf("%lf ", BufMatr[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	//нормировка А на 1
	for (int q = 0; q < N; q++) {
		for (int i = q; i < N; i++) {
			for (int j = q; j < N; j++) {
				BufMatr[i][j] = BufMatr[i][j] / BufMatr[i][i];
				MatrRev[i][j] = MatrRev[i][j] /BufMatr[i][i];
			}
		}
	}
	printf("\n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			//printf("%lf ", BufMatr[i][j]);
		}
		printf("\n");
	}
	//-------------------------------------------

	printf("\n");
	//обратный ход Гаусса
	for (int q =N-1; q >= 0; q--) {
		sol[q] = free[q];
		for (int d = q + 1; d < N; d++) {
			sol[q] += -sol[d] * Matr[q][d];
		}
		sol[q] = sol[q] / Matr[q][q];
	}
	//решения Гаусс
	for (int q = 0; q < N; q++) {
		printf("solution %lf  \n", sol[q]);
	}

	//заполнение матрицы
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if ((i != j) && (abs(i - j) < 5)) {
				Matr[i][j] = 1;
			}
			else {
				Matr[i][j] = 0;
			}
			if (i == j) Matr[i][j] = 10;
		}
	}
	for (int q = 0; q < N; q++) {
		free[q] = q + 1;
	}
	//подсчет невязок
	double r[10];
	for (int q = 0; q < N; q++) {
		r[q] = q + 1;
		//printf("q %d r %lf \n",q, r[q]);
		for (int j = 0; j < N; j++) {
			r[q] += -Matr[q][j] * sol[j];
			//printf("q %d r %lf \n",q,r[q]);
		}
	}
	//невязки
	for (int q = 0; q < N; q++) {
		printf("r %lf  \n", r[q]);
	}
	//--------------------------------------------------------

	//Метод Зейделя ------------------------------------------
	//начальное приближение Зейделя
	for (int q = 0; q < N; q++) {
		solZ[q] = 0;
	}
	//Зейдель
	for (int q = 0; q < 10; q++) {
		for (int i = 0; i < N; i++) {
			solZ[i] = free[i];
			for (int j = 0; j < N; j++) {
				if (i != j) {
					solZ[i] += -Matr[i][j] * solZ[j];
				}
			}
			solZ[i] = solZ[i] / Matr[i][i];
		}
	}
	for (int q = 0; q < N; q++) {
		printf("solutionZ %lf  \n", solZ[q]);
	}
	//подсчет невязок
	for (int q = 0; q < N; q++) {
		r[q] = q + 1;
		for (int j = 0; j < N; j++) {
			r[q] += -Matr[q][j] * solZ[j];
		}
	}
	//невязки Зейделя
	for (int q = 0; q < N; q++) {
		printf("r %lf  \n", r[q]);
	}
	//----------------------------------------------------------

	//максимальное собственное значение-------------------------
	for (int q = 0; q < N; q++) {
		solZ1[q] = 1;
	}
	for (int d = 0; d < 10; d++) {
		for (int q = 0; q < N; q++) {
			solZ[q] = solZ1[q];
			solZ1[q] = 0;
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				solZ1[i] += Matr[i][j] * solZ[j];
			}
		}
	}
	//Lmax
	double sc1 = 0;
	for (int i = 0; i < N; i++) {
		sc1 += solZ1[i] * solZ[i];
	}
	double sc2 = 0;
	for (int i = 0; i < N; i++) {
		sc2 += solZ[i] * solZ[i];
	}
	double Lmax = sc1 / sc2;
	//------------------------------------------------------------

	//минимальное собственное значение------------------------- minL=1/maxL(A^(-1))
	for (int q = 0; q < N; q++) {
		solZ1[q] = 1;
	}
	for (int d = 0; d < 10; d++) {
		for (int q = 0; q < N; q++) {
			solZ[q] = solZ1[q];
			solZ1[q] = 0;
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				solZ1[i] += MatrRev[i][j] * solZ[j];
			}
		}
	}
	//Lmin
	sc1 = 0;
	for (int i = 0; i < N; i++) {
		sc1 += solZ1[i] * solZ[i];
	}
	sc2 = 0;
	for (int i = 0; i < N; i++) {
		sc2 += solZ[i] * solZ[i];
	}
	double Lmin = sc1 / sc2;
	//------------------------------------------------------------
	printf("Lmax %lf Lmin %lf\n", Lmax , Lmin);
	printf("m %lf", Lmax / Lmin);
	scanf_s("%lf", sol[1]);
		return 0;
}