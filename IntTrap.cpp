#include <stdio.h>
#include <math.h>
double htrap(double e) {
	//достаточный шаг для формулы трапеций
	double x1 = sqrt(log(2));
	double maxddF = 1 / (4 * x1*pow((1 + x1), 2)) + 1 / (8 * pow(x1, 3)*pow((1 + x1), 2)) + 1 / (4 * pow(x1, 2)*pow((1 + x1), 3));
	double h = sqrt(12 * e / maxddF);
	return h;
}
double hSim(double e) {
	//достаточный шаг для формулы трапеций
	double x1 = sqrt(log(2));
	double maxddddF = (48 * pow(x1, 9) + 276 * pow(x1, 7) + 400 * pow(x1, 5) + 285 * pow(x1, 3) + 144 * pow(x1, 8) + 356 * pow(x1, 6) + 368 * pow(x1, 4) +177*pow(x1,2)+ 75 * x1 + 15 )/ (16 * pow(2, 4)*pow((x1 + 1), 5)*pow(x1, 7));
	double h = pow((2880 * e / maxddddF),1/4);
	return h;
}
double f(double x) {
	return 1 / (1 + sqrt(log(x)));
}

int main() {
	double e = 0.0001;
	//трапеции
	double h = htrap(e);
	printf("h %lf \n", h);
	double x = 2;
	double Itrap=0;
	int i = 0;
	while (x < 3) {
		double y = x + h;
		if (y < 3) {
			Itrap += h * (f(y) + f(x)) / 2;
			x = y;
		}
		else {
			Itrap += (3-x) * (f(3) + f(x)) / 2;
			x = 3;
		}
		i++;
		printf("i - %d x - %lf Itrap  %lf \n", i,x , Itrap);
	}
	printf("i - %d Itrap  %lf \n", i, Itrap);
	//симпсон
	h=hSim(e);
	printf("h %lf \n", h);
	x = 2;
	double ISim = 0;
	i = 0;
	while (x < 3) {
		if (x+h < 3) {
			ISim += h * (f(x+h) + 4*f(x+h/2)+f(x)) / 6;
			x+=2*h;
		}
		else {
			ISim += (3 - x) * (f(3) +4*f((3+x)/2)+ f(x)) / 6;
			x = 3;
		}
		i++;
		printf("i - %d x - %lf ISim  %lf \n", i, x, ISim);
	}
	scanf_s("%lf", &e);
	return 0;
}