#include <stdio.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <complex.h>

/*
	ALGUMAS CARACTERÍSTICAS PARA COMPILAÇÃO

	- gcc rtn_simple.c -o test -lm -std=c99 -lfftw3l
	dessa forma, o compilador consegue entender o que significa a função round()

*/

// define de alguns parâmetros padrões para quando o usuário não os define
#define COEFF_1 1000
#define COEFF_2 400
#define MC 5
#define IMPACT 1
#define DT 0.001

// complex number definition
#define REAL 0
#define IMAG 1


// função para receber os valores a partir do usuário
static void start_rtn (int argc, char **argv, int *coeff_1, int *coeff_2, int *mc, int *impact, double *dt) {

	const char *opt = "c:f:m:i:d:"; // definição das flag
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
			case 'c':
				sscanf(optarg, "%d", coeff_1);
				break;
			case 'f':
				sscanf(optarg, "%d", coeff_2);
				break;
			case 'm':
				sscanf(optarg, "%d", mc);
				break;
			case 'i':
				sscanf(optarg, "%d", impact);
				break;
			case 'd':
				sscanf(optarg, "%lf", dt);
				break;
		};
	}
}

double* rtn_calc(double *t, double dt, const double *tc, int impact, int coeff_1){
	// The time constants (tc) are [avg_time_in_high avg_time_in_low]
	
	double tau_local;
	double *trap_state, *x_return;
	double avg = 0.0;
	srand(time(0));
	tau_local = pow(((1/ tc[1]) + (1/ tc[0])), -1);				// tau definition - ok
	trap_state = (double*)malloc((coeff_1/dt)*sizeof(double));
	x_return = (double*)malloc((coeff_1/dt)*sizeof(double));

	// inialize the randomic variation of the trap state between 0 or 1

	for (int i=0; i <= (coeff_1/dt); i++){		// x = zeros
		x_return[i] = 0.0;
	}

	trap_state[0] = round((double)rand() /RAND_MAX);

	for (int i = 1; i <= (coeff_1/dt); i++){
		if(trap_state[i-1] == 1.0){
			trap_state[i] = ((double)rand() /RAND_MAX >= ((1-exp(-dt/tau_local))*(tc[1]/(tc[1]+tc[0]))));
		}  
		else{
			trap_state[i] = ((double)rand() /RAND_MAX <= ((1-exp(-dt/tau_local))*(tc[0]/(tc[1]+tc[0]))));	
		}
	}

	avg = 0.0;

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
		x_return[i] = x_return[i] - (avg / ((coeff_1/dt)+1));
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
	FILE *file_data;
	double lengthT = (coeff_1/dt), lengthTau = (coeff_2/dt);
	long double _Complex *test;
	long double *freq;

	// fft var
    fftwl_complex *in_fft;       // input array - equivalent to double x[n][2]
    fftwl_complex *out_fft;       // output array - equivalent to double x[n][2]
    fftwl_plan plan;        // fft plant declaration
    double fs = 0;          // sampling frequency
    double t_aux2 = 0;      // time aux
    double t_1 = -200;      // first step from time
    double t_2 = 200;       // second step from time
	double pi = 3.14159265359;
	fftwl_complex *ssd;
	fftwl_complex *FTsiga;

	t = (double*)malloc(((lengthT)+1) * sizeof(double));
	tau = (double*)malloc(((lengthTau)+1)*sizeof(double));
	Rx = (double*)malloc(((lengthTau)+1)*sizeof(double));
	freq = (long double*)malloc(((lengthTau)+1)*sizeof(long double));
	test = (long double _Complex*)malloc(((lengthTau)+1)*sizeof(long double _Complex));
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

	for (int i = 0; i <= (lengthTau); i++){		// Rx = zeros
		Rx[i] = 0.0;
	}


	for (int i = 0; i < MC; i++){
		x_new = rtn_calc(t, dt, tc, impact, coeff_1);

		avg = 0.0;
		for (int j = 0; j <= (lengthT); j++){
			avg += (dt*x_new[j]);
		}

		for (int j = 0; j <= (lengthT); j++){
			x_new[j] = x_new[j] - (avg/(lengthT+1)); // mudar para length t
		}

		pin = 0.0;

		for (int j = 0; j <= (lengthT); j++){
			pin += dt*(pow(x_new[j], 2));
		}
		P += pin/(lengthT); // aqui também

		for (int j = 0; j <= (lengthTau); j++){
			Rx[j] = x_new[(int)(((lengthT)+1)/2)+1] * x_new[(int)(((lengthT)+1)/2 - (int)(((lengthTau)+1)/2))+j+1] + Rx[j];
		}
	}

	// P calc and print
	printf("P = %le", P/MC);

	//printf("\n01 = %d\n", (int)(((lengthT)+1)/2 - (int)(((lengthTau)+1)/2)));

	/*printf("\n   RX: \n");
    for(int i=0; i<=(lengthTau); i++){
        printf("%2d - %10f\n", i, Rx[i]);
    }*/

    // memmory allocation for fft
    in_fft = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * (lengthTau+1));
    out_fft = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * (lengthTau+1));
	ssd = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * (lengthTau+1));
	FTsiga = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * (lengthTau+1));

    // plan definition
    plan = fftwl_plan_dft_1d((lengthTau)+1, in_fft, out_fft, FFTW_FORWARD, FFTW_ESTIMATE);

    // fill the array with data

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
	for(int i=0; i<=(lengthTau); i++){
		ssd[i][REAL] = out_fft[i][REAL]/(lengthTau+1);
		ssd[i][IMAG] = out_fft[i][IMAG]/(lengthTau+1);
	}

	// abs function for FTsiga
	/*for(int i=0; i<=(lengthTau); i++){
		FTsiga[i] = sqrt(pow(ssd[i][REAL], 2) + (pow(ssd[i][IMAG], 2)));
	}*/

	for(int i = 0; i<=lengthTau; i++){
		test[i] = sqrt(pow(ssd[i][REAL], 2) + (pow((ssd[i][IMAG])*I, 2)));
		//printf("\n %d = %Le + %Le i", i, creall(test[i]), cimagl(test[i]));
	}

	// C language cant apply abs function in float, double or long double numbers, so we need to write the data into a text file so we can
	// work with abs function in python. Thats fine because the plot are made in python.

    // fft print
    /*printf("\n   FFT array coefficients:\n");
    for(int i=0; i <= lengthTau; i++) {
        printf("%2d - %10Le    + %10Le \n", i, FTsiga[i][REAL], FTsiga[i][IMAG]);
    }*/

	// frequency definition
	for(int i=1; i<=lengthTau; i++){		// i começar no 1
		freq[i-1] = 1/-(coeff_2/2) + (i/(dt)); // 1 / 2*i/dt
		//printf("\n f = %Le\n", freq[0]);
	}

	//printf("\n f = %Le\n", freq[0]);

	// print the ssd data in a txt file
  	file_data = fopen("ssd_out.txt", "w");	// create a file with write permission

	if(file_data == NULL)					// is the file working?
    	printf("\nERRO! O arquivo não foi aberto!\n");
	else
     	printf("\nO arquivo foi aberto com sucesso!\n");

	//fprintf(file_data, "	Frequency  -  abs(SSD)\n");		// ssd print
	for(int i = 0; i <= lengthTau; i++){
		//fprintf(file_data, "ssd(%d)		-		%Le  +  %Le\n", i, ssd[i][REAL], ssd[i][IMAG]);
		fprintf(file_data, "%Le,%Le\n", freq[i], creall(test[i]));
	}

	printf("Dados gravados com sucesso!\n");

	fclose(file_data);						// close file

	// memmory clean
    fftwl_destroy_plan(plan);
    fftwl_free(in_fft); fftwl_free(out_fft); fftwl_free(ssd);
	free(t);
	free(tau);
	free(Rx);
	free(x_new);

	return 0;

}

