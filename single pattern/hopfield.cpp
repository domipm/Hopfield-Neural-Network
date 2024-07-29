//  HOPFIELD'S NEURAL NETWORK MODEL
//  Input: Binary matrix "moon.txt"
//  Output: Output binary matrix "output.txt", overlap for each step "overlap.txt"

#include<cmath>
#include<iostream>
#include<ctime>
#include"gsl_rng.h"

#define SEED 84982490479247 // Random number seed

#define FIN "moon.txt" // Input matrix filename

#define N 128   // Image size
#define ITER 20 // Number of Monte-Carlo steps to perform (iterations)

float T = 10e-3; // Temperature of the system

float w[N][N][N][N] = {0}, theta[N][N] = {0};   // Synaptic steps and trigger threshold
float a = 0, H = 0, overlap = 0;                   // Parameter a^mu, hamiltonian, and overlap
int d[N][N], s[N][N];                           // Pattern matrix and system matrix

bool noise_deform = 1; // Whether to deform (lambda) initial matrix
float lambda = 0.9; // Random noise deformation

using namespace std;

gsl_rng *tau;

//  Minimum between floats
float minimo(float a, float b) {

    if (a < b) return a;
    if (a > b) return b;
    if (a == b) return a;

    return 0;

}

//  Write (NxN) matrix onto (file)
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

//  Periodic conditions for an (NxN) matrix (x)
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

//  Random integer generator in range (0,1)
int rnd_int() {

    return gsl_rng_uniform_int(tau, 2);

}

int main() {

    extern gsl_rng *tau;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, SEED);

    //  Open necessary files
    FILE *in;
    in = fopen(FIN, "r");                   // Input file (with generated pattern)
    FILE *out;
    out = fopen("output.txt", "w");         // Output file (with iterations)
    FILE *foverlap;
    foverlap = fopen("solapamiento.txt", "w");  // Output file (overlap with pattern)

    int d[N][N]; // Input matrix
    int s[N][N]; // System matrix

    //  Read input file onto matrix
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if ( fscanf(in, "%d", &d[i][j]) );

    //  Initial configuration
    if (noise_deform == 0) { // Random

        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) s[i][j] = rnd_int();
        period(s); // Apply periodic conditions

    } else { // Deformed

        //  Set patterns and system matrix equal
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) s[i][j] = d[i][j];
        for (int m = 0; m < int(lambda*N*N); m++) { // Trigger random points
            int i = gsl_rng_uniform_int(tau, N);
            int j = gsl_rng_uniform_int(tau, N);
            s[i][j] = (1-s[i][j]);
        }
        period(s); // Apply periodic conditions

    }

    //  Write initial matrix state onto file
    outmat(out, s);

    //  Calculate parameter a^mu
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) a += d[i][j];
    a *= 1.0/(N*N);

    //  Calculate synaptic weghts w_ij,kl
    for (int i = 0; i < N; i++) 
    for (int j = 0; j < N; j++) 
    for (int k = 0; k < N; k++) 
    for (int l = 0; l < N; l++) 
        if ((i == k) && (j == l)) w[i][j][k][l] = 0.0;
        else w[i][j][k][l] = (1.0/(N*N))*(d[i][j]-a)*(d[k][l]-a);

    //  Calculate trigger threshold theta_ij
    for (int i = 0; i < N; i++) 
    for (int j = 0; j < N; j++) 
    for (int k = 0; k < N; k++) 
    for (int l = 0; l < N; l++) 
        theta[i][j] += 0.5*w[i][j][k][l];

    //  Calculate overlap between system and pattern matrices
    overlap = 0.0;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) overlap += (1.0/(N*N*a*(1.0-a)))*(d[i][j]-a)*(s[i][j]-a);
    fprintf(foverlap, "%f\n", overlap);

    //  Main loop (Over the Monte-Carlo steps)
    for (int b = 0; b < ITER; b++) {

        //  Secondary loop (N*N for each Monte-Carlo step)
        for (int t = 0; t < N*N; t++) {

            //  Random point (i,j) of the system matrix
            int i = gsl_rng_uniform_int(tau, N);
            int j = gsl_rng_uniform_int(tau, N);

            //  Calculate Hamiltonian (Delta H)
            H = theta[i][j]*(1-2*s[i][j]);
            for (int k = 0; k < N; k++) for (int l = 0; l < N; l++) H += 0.5*w[i][j][k][l]*s[k][l]*(2*s[i][j]-1);

            float p = minimo(1.0, exp(-H/T)); // Change of spin probability

            float r = gsl_rng_uniform(tau); 

            //  Change the spin if (r < p)
            if (r < p) {
                s[i][j] = (1 - s[i][j]);    // Change value
                period(s);                  // Periodic conditions
            }

        }

        //  Calculate overlap between system and pattern matrices
        overlap = 0.0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) overlap += (1.0/(N*N*a*(1.0-a)))*(d[i][j]-a)*(s[i][j]-a);
        fprintf(foverlap, "%f\n", overlap);

        //  Write modified system matrix onto output file
        outmat(out, s);

    }

    fclose(out);
    fclose(foverlap);
    fclose(in);

    return 0;

}