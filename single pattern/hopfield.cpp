//PROGRAMA PRINCIPAL DE LA RED NEURONAL DE HOPFIELD
//INPUT: Matriz binaria "moon.txt"
//OUTPUT: Matriz salida para cada paso "output.txt"
//OUTPUT: Salida solapamiento para cada paso "solapamiento.txt"

#include<cmath>
#include<iostream>
#include<ctime>
#include"gsl_rng.h"

#define SEED 84982490479247 //Semilla de números aleatorios

#define FIN "moon.txt" //Fichero con matriz de entrada

#define N 128 //Tamaño imagen
#define ITER 20 //Número iteraciones(Pasos Monte-Carlo)

float T = 10e-3; //Temperatura del sistema

float w[N][N][N][N] = {0}, theta[N][N] = {0}; //Pesos sinápticos y umbral de disparo
float a = 0, H = 0, slpm = 0; //Parámetro a^mu, hamiltoniano, y solapamiento
int d[N][N], s[N][N]; //Matriz patrón y sistema

bool deformado = 1; //ALEATORIO == 0, DEFORMADO == 1
float lambda = 0.9; //Deformación (Entre 0 y 1)

using namespace std;

gsl_rng *tau;

//Mínimo entre dos float
float minimo(float a, float b) {

    if (a < b) return a;
    if (a > b) return b;
    if (a == b) return a;

    return 0;

}

//Escribe matriz (NxN) en el fichero (file)
void outmat(FILE *file, int x[N][N]) {

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(file, "%d\t", x[i][j]);
        }
        fprintf(file, "\n");
        if (i == N-1) fprintf(file, "\n");
    }

    return;
}

//Condiciones periódicas para una matriz (x) de dimensiones (NxN)
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

//Generador int aleatorio {0,1}
int rnd_int() {

    return gsl_rng_uniform_int(tau, 2);

}

int main() {

    extern gsl_rng *tau;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, SEED);

    //Abrimos ficheros
    FILE *in;
    in = fopen(FIN, "r"); //Fichero entrada (Con patrón generado)
    FILE *out;
    out = fopen("output.txt", "w"); //Fichero salida (Con iteraciones)
    FILE *sola;
    sola = fopen("solapamiento.txt", "w"); //Fichero salida (Solapamiento con patrón)

    int d[N][N]; //Matriz entrada
    int s[N][N]; //Matriz sistema

    //Leemos fichero entrada en matriz
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if ( fscanf(in, "%d", &d[i][j]) );

    //Configuración inicial
    if (deformado == 0) { //Aleatorio

        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) s[i][j] = rnd_int();
        period(s); //Aplicamos condiciones periódicas

    } else { //Deformado

        //Igualamos matriz patrón y matriz sistema
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) s[i][j] = d[i][j];
        for (int m = 0; m < int(lambda*N*N); m++) { //Deformamos valores aleatorios
            int i = gsl_rng_uniform_int(tau, N);
            int j = gsl_rng_uniform_int(tau, N);
            s[i][j] = (1-s[i][j]);
        }
        period(s); //Aplicamos condiciones periódicas

    }

    //Escribimos matriz sistema inicial en fichero
    outmat(out, s);

    //Calculamos parámetro a^mu
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) a += d[i][j];
    a *= 1.0/(N*N);

    //Calculamos pesos sinápticos w_ij,kl
    for (int i = 0; i < N; i++) 
    for (int j = 0; j < N; j++) 
    for (int k = 0; k < N; k++) 
    for (int l = 0; l < N; l++) 
        if ((i == k) && (j == l)) w[i][j][k][l] = 0.0;
        else w[i][j][k][l] = (1.0/(N*N))*(d[i][j]-a)*(d[k][l]-a);

    //Calculamos umbral de disparo theta_ij
    for (int i = 0; i < N; i++) 
    for (int j = 0; j < N; j++) 
    for (int k = 0; k < N; k++) 
    for (int l = 0; l < N; l++) 
        theta[i][j] += 0.5*w[i][j][k][l];

    //Calculamos solapamiento entre matriz sistema y matriz patrón
    slpm = 0.0;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) slpm += (1.0/(N*N*a*(1.0-a)))*(d[i][j]-a)*(s[i][j]-a);
    fprintf(sola, "%f\n", slpm);

    //Bucle principal (Sobre los pasos Monte-Carlo)
    for (int b = 0; b < ITER; b++) {

        //Bucle secundario (N*N para cada paso Monte-Carlo)
        for (int t = 0; t < N*N; t++) {

            //Punto aleatorio (i,j) de matriz sistema
            int i = gsl_rng_uniform_int(tau, N);
            int j = gsl_rng_uniform_int(tau, N);

            //Calculamos el hamiltoniano (Delta H)
            H = theta[i][j]*(1-2*s[i][j]);
            for (int k = 0; k < N; k++) for (int l = 0; l < N; l++) H += 0.5*w[i][j][k][l]*s[k][l]*(2*s[i][j]-1);

            float p = minimo(1.0, exp(-H/T)); //Probabilidad cambio de spin

            float r = gsl_rng_uniform(tau); 

            //Cambiamos el spin si se cumple (r < p)
            if (r < p) {
                s[i][j] = (1 - s[i][j]); //Cambio de valor
                period(s); //Condiciones periódicas
            }

        }

        //Calculamos solapamiento entre matriz sistema y matriz patrón
        slpm = 0.0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) slpm += (1.0/(N*N*a*(1.0-a)))*(d[i][j]-a)*(s[i][j]-a);
        fprintf(sola, "%f\n", slpm);

        //Escribimos matriz sistema modificada en el fichero de salida
        outmat(out, s);

    }

    fclose(out);
    fclose(sola);
    fclose(in);

    return 0;

}