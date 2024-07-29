//PROGRAMA PRINCIPAL DE LA RED NEURONAL DE HOPFIELD PARA VARIOS PATRONES
//INPUT: Matrices patrón y deformadas "a.txt", ..., "a_dist.txt", ...
//OUTPUT: Matriz salida para cada paso "output.txt"
//OUTPUT: Salida solapamiento para cada paso y patrón "solapamiento.txt"

#include<cmath>
#include<iostream>
#include<ctime>
#include"gsl_rng.h"

#define SEED 986453218741 //Semilla de números aleatorios

#define N 64 //Tamaño imagen
#define N_PATRONES 3 //Número de patrones almacenados
#define ITER 20 //Número iteraciones(Pasos Monte-Carlo)

float T = 0.0001; //Temperatura

float w[N][N][N][N] = {0}, theta[N][N] = {0}; //Pesos sinápticos y umbral de disparo
float a[N_PATRONES] = {0}; //Parámetro a^mu
float H = 0; //Hamiltoniano
float sol[N_PATRONES] = {0}; //Solapamiento m^mu
int d[N_PATRONES][N][N], s[N][N]; //Matriz patrones y sistema

bool deformado = false; //ALEATORIO == 0, DEFORMADO == 1
bool aleatorio = false; //INICIAL DEFORMADO == 0, INICIAL ALEATORIO == 1
float lambda = 0.01; //Deformación (Entre 0 y 1)

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
void solapamiento(FILE *file) {

    for (int mu = 0; mu < N_PATRONES; mu++) {
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            sol[mu] += (d[mu][i][j]-a[mu])*(s[i][j]-a[mu]);
        sol[mu] *= 1.0/(N*N*a[mu]*(1-a[mu]));
    }
    for (int mu = 0; mu < N_PATRONES; mu++) {
        fprintf(file, "%.15lf\t", sol[mu]);
    }
    fprintf(file, "\n");

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

    //Abrimos ficheros
    FILE *out;
    out = fopen("output.txt", "w"); //Fichero salida (Con iteraciones)
    FILE *fsol;
    fsol = fopen("solapamiento.txt", "w"); //Fichero salida (Solapamiento entre patrones)

    extern gsl_rng *tau;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, SEED);

    //Cargamos los patrones en la matriz
    load_pattern("a.txt", d[0]);
    load_pattern("b.txt", d[1]);
    load_pattern("c.txt", d[2]);

    if (aleatorio == true) { //Configuración inicial aleatoria
        
        for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++) 
            s[i][j] = rnd_int();
        period(s);

    } else { //Configuración patrón inicial (letra B deformada)
        
        load_pattern("b_dist.txt", s);
        if (deformado = true) { //Deformamos aleatoriamente

            for (int m = 0; m < int(lambda*N*N); m++) { //Deformamos valores aleatorios
                int i = gsl_rng_uniform_int(tau, N);
                int j = gsl_rng_uniform_int(tau, N);
                s[i][j] = (1-s[i][j]);
            }
            period(s); //Aplicamos condiciones periódicas
        }

    }

    outmat(out, s);

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

    //Escribimos solapamiento inicial
    solapamiento(fsol);

    //Bucle principal (Sobre los pasos Monte-Carlo)
    for (int b = 0; b < ITER; b++) {

        //Bucle secundario (N*N para cada paso Monte-Carlo)
        for (int t = 0; t < N*N; t++) {

            //Punto aleatorio (i,j) de matriz sistema
            int i = gsl_rng_uniform_int(tau, N);
            int j = gsl_rng_uniform_int(tau, N);

            //Calculamos el Hamiltoniano (Delta H)
            H = theta[i][j]*(1-2*s[i][j]);
            for (int k = 0; k < N; k++) for (int l = 0; l < N; l++) H += 0.5*w[i][j][k][l]*s[k][l]*(2*s[i][j]-1);

            float p = minimo(1.0, exp(-H/T)); //Probabilidad cambio de spin

            float r = gsl_rng_uniform(tau); 

            //Cambiamos el spin ( en el caso r < p)
            if (r < p) {
                s[i][j] = (1 - s[i][j]);
                period(s);
            }

        }

        //Escribimos matriz modificada en el fichero
        outmat(out, s);

        // Calculamos solapamiento con los patrones y lo escribimos en fichero
        solapamiento(fsol);

    }    

    fclose(out);
    fclose(fsol);

    return 0;

}