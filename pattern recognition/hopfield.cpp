//  Main Hopfield Neural Network for multiple patterns
//  Input: Pattern matrices "a.txt",... and distorted (handwritten) "a_dist.txt",...
//  Output: System matrix "output.txt" for each step, overlap for each step and pattern "overlap.txt"

#include<cmath>
#include<iostream>
#include<ctime>
#include"gsl_rng.h"

#define SEED 986453218741 // Random number generator seed

#define N 64            // Image size
#define N_PATRONES 3    // Number of patterns saved
#define ITER 20         // Number of iterations (Monte-Carlo steps)

float T = 0.0001; // Temperature

float w[N][N][N][N] = {0}, theta[N][N] = {0}; // Synaptic weights and trigger threshold
float a[N_PATRONES] = {0};                    // Parameter a^mu
float H = 0;                                  // Hamiltonian
float sol[N_PATRONES] = {0};                  // Overlap m^mu
int d[N_PATRONES][N][N], s[N][N];             // System and pattern matrices

bool deformado = false; // Random == 0, Deformed == 1
bool aleatorio = false; // Initially deformed == 0, Initially random == 1
float lambda = 0.01;    // Deformation (0,1)

using namespace std;

gsl_rng *tau;

//  Save pattern from file (fname) in matrix (mat)
void load_pattern(const char* fname, int mat[N][N]) {
    FILE *file;
    file = fopen(fname, "r");
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if ( fscanf(file, "%d", &mat[i][j]) );
    fclose(file);
    return;
}

//  Minimum value between floats
float minimo(float a, float b) {
    if (a < b) return a;
    if (a > b) return b;
    if (a == b) return a;
    return 0;
}

//  Write (NxN) matrix onto file (fname)
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

//  Apply periodic boundary conditions to the (NxN) matrix (x)
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

//  Overlap of matrix with patterns
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

//  Random number generator (0,1)
int rnd_int() {
    return gsl_rng_uniform_int(tau, 2);
}

//  Generate random patterin in matrix (mat)
void random_pattern(int mat[N][N]) {
    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
        mat[i][j] = rnd_int();
    period(mat);
    return;
}

int main() {

    //  Open files
    FILE *out;
    out = fopen("output.txt", "w");         // Output file (iterations)
    FILE *fsol;
    fsol = fopen("solapamiento.txt", "w");  // Output file (overlap with patterns)

    extern gsl_rng *tau;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, SEED);

    //  Load patterns onto matrices
    load_pattern("a.txt", d[0]);
    load_pattern("b.txt", d[1]);
    load_pattern("c.txt", d[2]);

    if (aleatorio == true) {        // Random initial configuration
        
        for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++) 
            s[i][j] = rnd_int();
        period(s);

    } else {                        // Deformed pattern initial configuration
        
        load_pattern("b_dist.txt", s);
        if (deformado = true) { // Randomly add noise

            for (int m = 0; m < int(lambda*N*N); m++) {
                int i = gsl_rng_uniform_int(tau, N);
                int j = gsl_rng_uniform_int(tau, N);
                s[i][j] = (1-s[i][j]);
            }
            period(s); // Apply periodic boundary conditions
        }

    }

    outmat(out, s);

    //  Calculate parameter a^mu
    for (int mu = 0; mu < N_PATRONES; mu++) {
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            a[mu] += d[mu][i][j];
        a[mu] *= 1.0/(N*N);
    }

    //  Calculate synaptic weights
    for (int i = 0; i < N; i++) 
    for (int j = 0; j < N; j++) 
    for (int k = 0; k < N; k++) 
    for (int l = 0; l < N; l++)
    for (int mu = 0; mu < N_PATRONES; mu++) 
        if ((i == k) && (j == l)) w[i][j][k][l] = 0.0;
        else w[i][j][k][l] += (1.0/(N*N))*(d[mu][i][j]-a[mu])*(d[mu][k][l]-a[mu]);
    
    //  Calculate trigger threshold
    for (int i = 0; i < N; i++) 
    for (int j = 0; j < N; j++) 
    for (int k = 0; k < N; k++) 
    for (int l = 0; l < N; l++) 
        theta[i][j] += 0.5*w[i][j][k][l];

    //  Write initial overlap
    solapamiento(fsol);

    //  Main loop (over Monte-Carlo steps)
    for (int b = 0; b < ITER; b++) {

        //  Secondary loop (N*N for each Monte-Carlo step)
        for (int t = 0; t < N*N; t++) {

            //  Random point (i,j) of system matrix
            int i = gsl_rng_uniform_int(tau, N);
            int j = gsl_rng_uniform_int(tau, N);

            //  Obtain Hamiltonian (Delta H)
            H = theta[i][j]*(1-2*s[i][j]);
            for (int k = 0; k < N; k++) for (int l = 0; l < N; l++) H += 0.5*w[i][j][k][l]*s[k][l]*(2*s[i][j]-1);

            float p = minimo(1.0, exp(-H/T)); // Probability change of spin

            float r = gsl_rng_uniform(tau); 

            //  Change spin if (r < p)
            if (r < p) {
                s[i][j] = (1 - s[i][j]);
                period(s);
            }

        }

        //  Write modified matrix onto file
        outmat(out, s);

        //  Calculate overlap and write onto file
        solapamiento(fsol);

    }    

    fclose(out);
    fclose(fsol);

    return 0;

}