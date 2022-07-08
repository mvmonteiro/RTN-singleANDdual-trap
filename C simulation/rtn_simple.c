#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>

/*
	ALGUMAS CARACTERÍSTICAS PARA COMPILAÇÃO

	- gcc rtn_simple.c -o test -lm -std=c99
	dessa forma, o compilador consegue entender o que significa a função round()

*/

// define de alguns parâmetros padrões para quando o usuário não os define
#define COEFF_1 1000
#define COEFF_2 400
#define MC 100
#define IMPACT 1
#define DT 0.001

// função para receber os valores a partir do usuário
static void start_rtn (int argc, char **argv, int *coeff_1, int *coeff_2, int *mc, int *impact, int *dt, char **path) {

	const char *opt = "b:d:p:f:c:"; // definição das flag
	int c; // ver o caso da flag que foi modificada e entrar no switch

// Inicialização dos valores padrões quando não modificado os parâmetros
	*coeff_1 = COEFF_1;
	*coeff_2 = COEFF_2;
	*mc = MC;
	*impact = IMPACT;
	*dt = DT;

// controle dos parâmetros advindo do terminal
	while ((c = getopt (argc, argv, opt)) != EOF) {
		switch (c) {
			case 'b':
				printf(".");
		}
	}
}


double* rtn_calc(double *t, double dt, const double *tc, int impact, int coeff_1){
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

	return x_return;
}

int main (){
	double dt, t_aux, tau_aux;
	double *t, *tau, *Rx, *x_new;
	int coeff_1 = COEFF_1;
	int coeff_2 = COEFF_2;
	int mc = MC;
	int impact = IMPACT;
	double tc[2] = {0.1, 1};
	double avg, pin, P = 0.0;
	dt = DT;

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

	
	for (int i = 0; i <= MC; i++){
		x_new = rtn_calc(t, dt, tc, impact, coeff_1);
		avg = 0.0;
		for (int j = 0; j <= (coeff_1/dt); j++){
			avg += (dt*x_new[j]);
		}

		for (int j = 0; j <= (coeff_1/dt); j++){
			x_new[j] = x_new[j] - (avg/1000); // mudar para length t
		}
		pin = 0.0;

		for (int j = 0; j <= (coeff_1/dt); j++){
			pin += dt*(pow(x_new[j], 2));
		}
		P += pin/1000; // aqui também

		for (int j = 0; j <= (coeff_2/dt); j++){
			Rx[j] = x_new[(int)((coeff_1/dt)+1)/2] * x_new[(int)((coeff_1/dt)+1)/2 - (int)(((coeff_2/dt)+1)/2)+j] + Rx[j];
		}
	}


	free(t);
	free(tau);
	free(Rx);
	free(x_new);


	return 0;

}

