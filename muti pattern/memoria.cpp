//PROGRAMA PARA DETERMINAR LA CAPACIDAD DE RECORDAR DE LA RED EN FUNCIÓN DEL NÚMERO DE PATRONES
//OUTPUT: Fichero con solapamiento en función del número de patrones "memoria.txt"

#include<cmath>
#include<iostream>
#include<ctime>
#include"gsl_rng.h"

#define SEED 980460521793 //Semilla

#define N 20 //Tamaño de las matrices
#define N_PATRONES 50 //Número de patrones máximo

float T = 10e-4; //Temperatura del sistema

float w[N][N][N][N] = {0}, theta[N][N] = {0}; //Matriz pesos y umbral de disparo
float a[N_PATRONES] = {0}; //Parámetro a^mu
float H = 0; //Hamiltoniano
float m[N_PATRONES] = {0}; //Solapamiento para cada patrón
int d[N_PATRONES][N][N], s[N][N]; //Matrices de patrones y sistema

using namespace std;

gsl_rng *tau;

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
    fsol = fopen("memoria.txt", "w");

    extern gsl_rng *tau;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, SEED);

    //Bucle patrones
    for (int n_pat = 1; n_pat <= N_PATRONES; n_pat++) {

        //Generamos patrones aleatorios
        for (int i = 0; i < n_pat; i++)
            random_pattern(d[i]);

        //Configuración inicial aleatoria
        for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++) 
            s[i][j] = rnd_int();
        period(s);

        //Calculamos parámetro a^mu
        for (int mu = 0; mu < n_pat; mu++) a[mu] = 0.0;
        for (int mu = 0; mu < n_pat; mu++) {
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
            w[i][j][k][l] = 0.0;
        for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++) 
        for (int k = 0; k < N; k++) 
        for (int l = 0; l < N; l++)
        for (int mu = 0; mu < n_pat; mu++) 
            if ((i == k) && (j == l)) w[i][j][k][l] = 0.0;
            else w[i][j][k][l] += (1.0/(N*N))*(d[mu][i][j]-a[mu])*(d[mu][k][l]-a[mu]);

        //Calculamos umbral de disparo
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            theta[i][j] = 0.0;
        for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++) 
        for (int k = 0; k < N; k++) 
        for (int l = 0; l < N; l++) 
            theta[i][j] += 0.5*w[i][j][k][l];

        //Calculamos solapamiento inicial
        for (int mu = 0; mu < n_pat; mu++) m[mu] = 0;
        for (int mu = 0; mu < n_pat; mu++) {
            for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                m[mu] += (d[mu][i][j]-a[mu])*(s[i][j]-a[mu]);
            m[mu] *= 1.0/(N*N*a[mu]*(1-a[mu]));
        }

        int cont = 0;

        //Bucle principal (pasos Monte Carlo)
        for (int i = 0; i < 100; i++) {

            //Bucle sobre todos los puntos
            for (int t = 0; t < N*N; t++) {

                int i = gsl_rng_uniform_int(tau, N);
                int j = gsl_rng_uniform_int(tau, N);

                //Calculamos el Hamiltoniano (Delta H)
                H = theta[i][j]*(1-2*s[i][j]);
                for (int k = 0; k < N; k++) for (int l = 0; l < N; l++) H += 0.5*w[i][j][k][l]*s[k][l]*(2*s[i][j]-1);

                float p = minimo(1.0, exp(-H/T)); 
                float r = gsl_rng_uniform(tau); 

                //Cambiamos el spin (en el caso r < p)
                if (r < p) {
                    s[i][j] = (1 - s[i][j]);
                    period(s);
                }

            }

        }

        //Calculamos solapamientos finales
        for (int mu = 0; mu < n_pat; mu++) {
            for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                m[mu] += (d[mu][i][j]-a[mu])*(s[i][j]-a[mu]);
            m[mu] *= 1.0/(N*N*a[mu]*(1-a[mu]));
        }
        //Solapamientos con valor > 0.75
        for (int mu = 0; mu < n_pat; mu++) if (m[mu] > 0.75 || m[mu] < -0.75) cont++;

        cout << n_pat << "\t" << cont << endl;
        fprintf(fsol, "%d\t%d\n", n_pat, cont);

    }

    fcloseall(); //Cerramos todos los ficheros

    return 0;

}