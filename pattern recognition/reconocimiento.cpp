//PROGRAMA QUE DETERMINA LA CAPACIDAD DE RECONOCER CADA LETRA
//INPUT: Matrices originales "a.txt", ..., "a_dist.txt", ...
//OUTPUT: Fichero de solapamiento para cada patrón, en cada iteración y con su parámetro lambda "reconocimiento.txt"

#include<cmath>
#include<iostream>
#include<ctime>
#include"gsl_rng.h"

#define SEED 986433257998741 //Semilla

#define N 64 //Tamaño matriz
#define N_PATRONES 3 //Número de patrones almacenados
#define ITER 10 //Número de pasos Monte Carlo

float T = 10e-4; //Temperatura
float lambda = 0; //Parámetro de deformación (se tomará aleatorio posteriormente)

float w[N][N][N][N] = {0}, theta[N][N] = {0}; 
float a[N_PATRONES] = {0};
float H = 0;
float sol[N_PATRONES] = {0};
int d[N_PATRONES][N][N], s[N][N]; 

using namespace std;

gsl_rng *tau;

//Guarda un patrón de un fichero (fname) en la matriz (mat)
void load_pattern(const char* fname, int mat[N][N]) {
    FILE *file;
    file = fopen(fname, "r");
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if ( fscanf(file, "%d", &mat[i][j]) );
    fclose(file);
    return;
}

//Mínimo entre dos float
float minimo(float a, float b) {
    if (a < b) return a;
    if (a > b) return b;
    if (a == b) return a;
    return 0;
}

//Escribe matriz (NxN) en un fichero (fname)
void outmat(FILE *file, int x[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(file, "%d ", x[i][j]);
        }
        fprintf(file, "\n");
        if (i == N-1) fprintf(file, "\n");
    }
    return;
}

//Aplicamos condiciones periódicas a la matriz (NxN)
void period(int x[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            x[N-1][N-1] = x[0][0];
            x[i][N-1] = x[i][0];
            x[N-1][j] = x[0][j];
        }
    }
    return;
}

//Solapamiento matriz con patrones
void solapamiento(FILE *file, float letra, float lambda) {

    for (int mu = 0; mu < N_PATRONES; mu++) {
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            sol[mu] += (d[mu][i][j]-a[mu])*(s[i][j]-a[mu]);
        sol[mu] *= 1.0/(N*N*a[mu]*(1-a[mu]));
    }
    fprintf(file, "%.10f\t", letra+1);
    for (int mu = 0; mu < N_PATRONES; mu++) {
        fprintf(file, "%.15lf\t", sol[mu]);
    }
    fprintf(file, "%.5f\n", lambda);

    return;

}

//Generador de enteros aleatorios {0,1}
int rnd_int() {
    return gsl_rng_uniform_int(tau, 2);
}

//Genera un patrón aleatorio en la matriz (mat)
void random_pattern(int mat[N][N]) {
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
        mat[i][j] = rnd_int();
    period(mat);
    return;
}

int main() {

    FILE *fsol;
    fsol = fopen("reconocimiento.txt", "w");

    extern gsl_rng *tau;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, SEED);

    //Cargamos los patrones en la matriz
    load_pattern("a.txt", d[0]);
    load_pattern("b.txt", d[1]);
    load_pattern("c.txt", d[2]);

    //Calculamos parámetro a^mu
    for (int mu = 0; mu < N_PATRONES; mu++) {
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            a[mu] += d[mu][i][j];
        a[mu] *= 1.0/(N*N);
    }

    //Calculamos pesos sinápticos
    for (int i = 0; i < N; i++) 
    for (int j = 0; j < N; j++) 
    for (int k = 0; k < N; k++) 
    for (int l = 0; l < N; l++)
    for (int mu = 0; mu < N_PATRONES; mu++) 
        if ((i == k) && (j == l)) w[i][j][k][l] = 0.0;
        else w[i][j][k][l] += (1.0/(N*N))*(d[mu][i][j]-a[mu])*(d[mu][k][l]-a[mu]);
    
    //Calculamos umbral de disparo
    for (int i = 0; i < N; i++) 
    for (int j = 0; j < N; j++) 
    for (int k = 0; k < N; k++) 
    for (int l = 0; l < N; l++) 
        theta[i][j] += 0.5*w[i][j][k][l];

    //Bucle simulaciones realizadas
    for (int z = 0; z <= 100; z++) {

        cout << z << endl; //Escribimos el número de la simulación para comprobar que el programa funciona

        //Escogemos aleatoriamente una de las tres letras
        int letra = gsl_rng_uniform_int(tau, 3);
        if (letra == 0) load_pattern("a_dist.txt", s);
        if (letra == 1) load_pattern("b_dist.txt", s);
        if (letra == 2) load_pattern("c_dist.txt", s);

        //Parámetro lambda aleatorio
        lambda = gsl_rng_uniform(tau);

        //Deformamos el patrón
        for (int m = 0; m < int(lambda*N*N); m++) {
            int i = gsl_rng_uniform_int(tau, N);
            int j = gsl_rng_uniform_int(tau, N);
            s[i][j] = (1-s[i][j]);
        }
        period(s);

        //Bucle principal (Sobre los pasos Monte-Carlo)
        for (int b = 0; b < ITER; b++) {

            //Bucle secundario (N*N para cada paso Monte-Carlo)
            for (int t = 0; t < N*N; t++) {

                //Escogemos punto (i,j) aleatorio
                int i = gsl_rng_uniform_int(tau, N);
                int j = gsl_rng_uniform_int(tau, N);

                //Calculamos el Hamiltoniano (Delta H)
                H = theta[i][j]*(1-2*s[i][j]);
                for (int k = 0; k < N; k++) for (int l = 0; l < N; l++) H += 0.5*w[i][j][k][l]*s[k][l]*(2*s[i][j]-1);

                float p = minimo(1.0, exp(-H/T)); //Probabilidad de cambio de spin

                float r = gsl_rng_uniform(tau);

                //Cambiamos el spin ( en el caso r < p)
                if (r < p) {
                    s[i][j] = (1 - s[i][j]);
                    period(s);
                }

            }

        }    

        //Calculamos solapamiento con los patrones y lo escribimos en fichero
        solapamiento(fsol, letra, lambda);

    }

    fcloseall(); //Cerramos todos los ficheros

    return 0;

}