#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*
	ALGUMAS CARACTERÍSTICAS PARA COMPILAÇÃO

	- gcc rtn_simple.c -o test -lm -std=c99
	dessa forma, o compilador consegue entender o que significa a função round()

*/

double* rtn_test(double *t, double dt, const double *tc, int impact, int coeff_1){
	double tau_local;
	double *trap_state, *x_return;
	double avg = 0.0;
	srand(time(0));
	tau_local = pow(((1/ tc[1]) + (1/ tc[0])), -1);
	trap_state = (double*)malloc((coeff_1/dt)*sizeof(double));
	x_return = (double*)malloc((coeff_1/dt)*sizeof(double));

	for (int i = 0; i <= (coeff_1/dt); i++){
		trap_state[i] = round((float)rand() /RAND_MAX);
	}

	for (int i = 1; i <= (coeff_1/dt); i++){
		if(trap_state[i-1]){
			if(((float)rand() /RAND_MAX >= ((1-exp(-dt/tau_local))*(tc[1]/tc[1]+tc[0])))){
				trap_state[i] = 1;
			} else
				trap_state[i] = 0;
		} else if(!trap_state[i-1]){
			if(((float)rand() /RAND_MAX <= ((1-exp(-dt/tau_local))*(tc[0]/tc[1]+tc[0])))){
				trap_state[i] = 1;	
			} else
				trap_state[i] = 0;
		}

	}
	
	for (int i = 0; i <= (coeff_1/dt); i++){
		if((int)trap_state[i] == 1){
			x_return[i] = (double)impact/2;
		}
		else if((int)trap_state[i] == 0){
			x_return[i] = -(double)impact/2;
		}
		avg = avg + x_return[i];
	}

	for (int i = 0; i <= (coeff_1/dt); i++){
		x_return[i] = x_return[i] - (avg/(coeff_1/dt));
	}

	for (int i = 0; i <= (coeff_1/dt); i++){
		printf("%lf\n", x_return[i]);
	}

	printf("\n\n%lf", avg);

	return x_return;
}

int main (){
	double dt, t_aux, tau_aux;
	double *t, *tau, *Rx;
	int coeff_1 = 1000;
	int coeff_2 = 400;
	int MC = 100;
	int impact = 1;
	double tc[2] = {0.1, 1};
	dt = 0.001;

	t = (double*)malloc((coeff_1/dt) * sizeof(double));
	tau = (double*)malloc((coeff_2/dt)*sizeof(double));
	Rx = (double*)malloc((coeff_2/dt)*sizeof(double));
	t_aux = -(coeff_1/2);
	tau_aux = -(coeff_2/2);

	for (int i = 0; i <= (coeff_1/dt); i++){
		t[i] = t_aux;
		t_aux = t_aux + dt;
	}

	for (int i = 0; i <= (coeff_2/dt); i++){
		tau[i] = tau_aux;
		tau_aux += dt;
	}

	/*
	for (int i = 0; i <= (coeff_2/dt); i++){
		printf("%lf\n", tau[i]);
	}*/


	//printf("%lf\n", t[10]);
	double *save_variable;
	save_variable = rtn_test(t, dt, tc, impact, coeff_1);
	free(t);
	free(tau);


	return 0;

}

