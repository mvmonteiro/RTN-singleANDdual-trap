// author: Marcus V. Monteiro

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <complex.h>

/*
	COMPILATION METHOD

	- gcc rtn_dual.c -o rtn -lm -std=c99 -lfftw3l
	dessa forma, o compilador consegue entender o que significa a função round()
	e ter acesso a library de fft fftw3

*/

// define de alguns parâmetros padrões para quando o usuário não os define
#define COEFF_1 1000	// time lenght 1 by 1 (case t = -500 to 500 we have coeff_1 = 1000)
#define COEFF_2 400		// tau lenght 1 by 1 (case t = -200 to 200 we have coeff_2 = 400)
#define MC 100			// Monte Carlo
#define IMPACT 0.005	// impact value
#define DT 0.001		// step

// complex number definition
#define REAL 0
#define IMAG 1


// função para receber os valores a partir do usuário
static void start_rtn (int argc, char **argv, int *coeff_1, int *coeff_2, int *mc, double *dt) {

	const char *opt = "c:f:m:i:d:"; // definição das flag
	int c; // ver o caso da flag que foi modificada e entrar no switch

	// Inicialização dos valores padrões quando não modificado os parâmetros
	*coeff_1 = COEFF_1;
	*coeff_2 = COEFF_2;
	*mc = MC;
	*dt = DT;

	// controle dos parâmetros advindo do terminal
	while ((c = getopt (argc, argv, opt)) != EOF) {
		switch (c) {
			case 'c':
				sscanf(optarg, "%d", coeff_1);
				break;
			case 'f':
				sscanf(optarg, "%d", coeff_2);
				break;
			case 'm':
				sscanf(optarg, "%d", mc);
				break;
			case 'd':
				sscanf(optarg, "%lf", dt);
				break;
		};
	}
}

double* rtn_calc(double *t, double dt, const double *tc_T1, const double *tc_T2, double *impact, int coeff_1){
	// The time constants (tc) are [avg_time_in_high avg_time_in_low]
	
	double tau_local_T1, tau_local_T2;
	double *trap_state_1, *trap_state_2, *x_return, *z, *k;
	double avg = 0.0, aux;
	tau_local_T1 = pow(((1/ tc_T1[1]) + (1/ tc_T1[0])), -1);				// tau definition
	tau_local_T2 = pow(((1/ tc_T2[1]) + (1/ tc_T2[0])), -1);

	// seed of random numbers based on a file of unix (using time was not so efficient because they repeat in each one second)
	uint32_t seed=0;
  	FILE *devrnd = fopen("/dev/random","r");
  	fread(&seed, 4, 1, devrnd);
  	fclose(devrnd);
  	srand(seed);

	// memmory allocation
	trap_state_1 = (double*)malloc((coeff_1/dt)*sizeof(double));
    trap_state_2 = (double*)malloc((coeff_1/dt)*sizeof(double));
	x_return = (double*)malloc((coeff_1/dt)*sizeof(double));
	z = (double*)malloc((coeff_1/dt)*sizeof(double));
	k = (double*)malloc((coeff_1/dt)*sizeof(double));

	// inialize the randomic variation of the trap state between 0 or 1

	for (int i=0; i <= (coeff_1/dt); i++){		// x = zeros
		x_return[i] = 0.0;
	}

    // zeros traps vector

    for (int i=0; i<= (coeff_1/dt); i++){
        trap_state_1[i] = 0.0;
    }

     for (int i=0; i<= (coeff_1/dt); i++){
        trap_state_2[i] = 0.0;
    }

	trap_state_1[0] = (rand() > RAND_MAX / 2);   // 0 false and 1 true
    trap_state_2[0] = (rand() > RAND_MAX / 2);   // 0 false and 1 true

    // RTN chance calculation
	for (int i = 1; i <= (coeff_1/dt); i++){
		if(trap_state_1[i-1] == 1.0){
			z[i] = ((double)rand() / (double)RAND_MAX);
			trap_state_1[i] = ( z[i] >= ((1-exp(-dt/tau_local_T1))*(tc_T1[1]/(tc_T1[1]+tc_T1[0]))));
		}  
		else{
			k[i] = ((double)rand() / (double)RAND_MAX);
			trap_state_1[i] = ( k[i] <= ((1-exp(-dt/tau_local_T1))*(tc_T1[0]/(tc_T1[1]+tc_T1[0]))));	
		}
	}

    	for (int i = 1; i <= (coeff_1/dt); i++){
		if(trap_state_2[i-1] == 1.0){
			z[i] = ((double)rand() / (double)RAND_MAX);
			trap_state_2[i] = ( z[i] >= ((1-exp(-dt/tau_local_T2))*(tc_T2[1]/(tc_T2[1]+tc_T2[0]))));
		}  
		else{
			k[i] = ((double)rand() / (double)RAND_MAX);
			trap_state_2[i] = ( k[i] <= ((1-exp(-dt/tau_local_T2))*(tc_T2[0]/(tc_T2[1]+tc_T2[0]))));	
		}
	}

	avg = 0.0;

    // impact calculation
	for (int i = 0; i <= (coeff_1/dt); i++){
        if (trap_state_1[i]==1 && trap_state_2[i]==1)
            x_return[i] = impact[3];
        else if(trap_state_1[i])
            x_return[i] = impact[1];
        else if(trap_state_2[i])
            x_return[i] = impact[2];
        else
			x_return[i] = impact[0];

        avg = x_return[i] + avg;
	}

	for (int i = 0; i <= (coeff_1/dt); i++){
		x_return[i] = x_return[i] - (avg / ((coeff_1/dt)+1));
	}

	free(z); free(k); free(trap_state_1); free(trap_state_2);
	return x_return;

}

int main (){
	double dt = DT, impact[4]={0.0, 1.0, 2.0, 4.0}, t_aux, tau_aux;
	double *t, *tau, *Rx, *x_new;
	int coeff_1 = COEFF_1;
	int coeff_2 = COEFF_2;
	int mc = MC;
	double tc_T1[2] = {0.01, 0.01}, tc_T2[2] = {10, 10};
	double avg, pin, P = 0.0, temp=1.0;
	FILE *file_data;
	double lengthT = (coeff_1/dt), lengthTau = (coeff_2/dt);
	long double *freq, *abs;
	int aux;

	// fft var
    fftwl_complex *in_fft;      // input array - equivalent to double x[n][2]
    fftwl_complex *out_fft;     // output array - equivalent to double x[n][2]
    fftwl_plan plan;        	// fft plant declaration
	fftwl_complex *ssd;			// spectral density var

	// memmory allocation
	t = (double* )malloc(((lengthT)+1) *sizeof(double));
	tau = (double* )malloc(((lengthTau)+1) *sizeof(double));
	Rx = (double* )malloc(((lengthTau)+1) *sizeof(double));
	freq = (long double* )malloc(((lengthTau)+1) *sizeof(long double));
	abs = (long double* )malloc(((lengthTau/2)+1) *sizeof(long double));
	x_new = (double* )malloc(((lengthT)+1) *sizeof(double));
	t_aux = -(coeff_1/2);
	tau_aux = -(coeff_2/2);

	// time definition
	for (int i = 0; i <= (lengthT); i++){
		t[i] = t_aux;
		t_aux = t_aux + dt;
	}

	// tau definition
	for (int i = 0; i <= (lengthTau); i++){
		tau[i] = tau_aux;
		tau_aux += dt;
	}

	// Rx = zeros
	for (int i = 0; i <= (lengthTau); i++){
		Rx[i] = 0.0;
	}

	// impact calculation
	for (int i=0; i < 21; i++){
		temp = temp + 0.1;
	}

	printf("\n %f", temp);

	impact[0] = 0; impact[1] = 1; impact[2] = 2; impact[3] = temp;

	/*for (int i = 0; i <= lengthT; i++){
		x_new[i] = i*0.37;
	}*/ // TEST

	// power and Rx function
	for (int i = 0; i < MC; i++){
		x_new = rtn_calc(t, dt, tc_T1, tc_T2, impact, coeff_1);

		avg = 0.0;
		for (int j = 0; j <= (lengthT); j++){
			avg += (dt*x_new[j]);
		}

		for (int j = 0; j <= (lengthT); j++){
			x_new[j] = x_new[j] - (avg/(lengthT+1));
		}

		pin = 0.0;

		for (int j = 0; j <= (lengthT); j++){
			pin += dt*(pow(x_new[j], 2));
		}
		P += pin/(coeff_1);

		for (int j = 0; j <= (lengthTau); j++){
			Rx[j] = x_new[(int)(((lengthT)+1)/2)+1] * x_new[(int)(((((lengthT)+1)/2)+1) - (int)(((((lengthTau)+1)/2))+1))+j+1] + Rx[j];
		}
			
	}

	// P calc and print
	printf("P = %le", P/mc);

    // memmory allocation for fft
    in_fft = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * (lengthTau+1));
    out_fft = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * (lengthTau+1));
	ssd = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * (lengthTau+1));

    // plan definition
    plan = fftwl_plan_dft_1d((lengthTau)+1, in_fft, out_fft, FFTW_FORWARD, FFTW_ESTIMATE);

    // fill the in array with data
    for(int i = 0; i <= (lengthTau); i++){
        in_fft[i][REAL] = (Rx[i]/mc);
        in_fft[i][IMAG] = 0;
    }

    // initial array print
    /*printf("\n");
    printf("\n   Input FFT coefficients: \n");
    for(int i=0; i <= lengthTau; i++){
        printf("%2d - %10Le    + %10Le\n", i, in_fft[i][REAL], in_fft[i][IMAG]);
    }*/

    // fft with one dimesion execution
    fftwl_execute(plan);

	// SSD definition
	out_fft[0][IMAG] = 0.0; 				// for some reason the fft of the first position of the vector doesn't return a 0 in the imaginary part
	for(int i=0; i<=(lengthTau); i++){
		ssd[i][REAL] = out_fft[i][REAL]/(lengthTau+1);
		ssd[i][IMAG] = out_fft[i][IMAG]/(lengthTau+1);
	}

	// absolut calc from ssd
	for(int i = 0; i<=lengthTau/2; i++){
		abs[i] = sqrt(pow(ssd[i][REAL], 2) + pow(ssd[i][IMAG], 2));
		//printf("\n %d = %Le + %Le i", i, creall(absolut[i]), cimagl(absolut[i]));		// in case the user want to print all the array
	}

	// frequency definition
	for (int i=0; i <= (lengthTau/2); i++) { 
		freq[i]=i/(2*dt*(lengthTau)/2); 
	}

	// the plot of this abs(signal) is made usying python tools, so the code writes the abs data and its frequency in a txt file
	// printing the data into a file
  	file_data = fopen("psd_dual.txt", "w");	// create a file with write permission

	if(file_data == NULL)					// is the file working?
    	printf("\n\nERRO! O arquivo não foi aberto!\n" );
	else
     	printf("\n\nO arquivo foi aberto com sucesso!\n");

	fprintf(file_data, "frequency,power\n");		// first line of the two coloms: freq x power
	for(int i = 0; i <= lengthTau/2; i++){
		fprintf(file_data, "%Le,%Le\n", freq[i], abs[i]);
	}

	printf("Dados gravados com sucesso!\n");

	fclose(file_data);						// close file

	// by this point the code create and wrote a txt file with the abs((fft(Rx/mc)) / (length(tau)+1))

	// memmory clean
    fftwl_destroy_plan(plan);
    fftwl_free(in_fft); fftwl_free(out_fft); fftwl_free(ssd);
	free(t); free(tau);
	free(Rx); free(x_new);
	free(freq); free(abs);

	return 0;

}