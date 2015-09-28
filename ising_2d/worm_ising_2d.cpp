/*
program: worm_ising_2d.cpp
created: 2015-09-28 -- 23 CEST
author: tc
notes: worm-algorithm simulation of the Ising model on the two-dimensional
       square lattice
*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "mt19937ar.c"

int bond_number(int i1, int i2, int *x, int *y, int L);


int main(){

    // seed random number
    unsigned long int seed;
    seed=time(NULL);
    init_genrand(seed);

    //read parameters from the input file
    int L;
    unsigned long int nsteps;
    double T;
    FILE *in_file = fopen("input.dat","r");
    if (in_file == NULL ){
        printf("Error! Could not open input file (input.dat: L, T, nsteps)\n"); 
        return 2;
    } 
    fscanf(in_file,"%i %lf %li", &L, &T, &nsteps);
    fclose(in_file);


    double K = 1.0 / T, inv_K = 1.0 / K, P_acc;
    int N = L * L;
    int ira = 0, masha = 0, new_masha;
    int nb, ibond, delta_nb, Nb = 0, Nb_tot = 0;
    unsigned long int step, Z = 0;

    // build site->x and site->y tables
    int x[N], y[N];
    for (int site=0; site<N; site++){
        x[site] = site - (site / L) * L;
        y[site] = site / L;
    }
    // build neighbors table
    int nbr[N][4];
    for (int site=0; site<N; site++){
        nbr[site][0] = (site / L) * L + (site + 1 + N) % L;
        nbr[site][1] = (site + L) % N;
        nbr[site][2] = (site / L) * L + (site - 1 + N) % L;
        nbr[site][3] = (site - L + N) % N;
    }
    // initialize bond weights
    int bonds[N * 2];
    for (int b=0; b<N * 2; b++){
        bonds[b] = 0;
    }

    // main Monte Carlo loop
    for (step=0; step<nsteps; step++){
        if (ira == masha){
            ira = (int)floor(genrand() * N);
            masha = ira;
            Z += 1;
            Nb_tot -= Nb;
        }
        // shift move -- start
        new_masha = nbr[masha][(int)floor(genrand() * 4)];
        ibond = bond_number(masha, new_masha, x, y, L);
        nb = bonds[ibond];
        if (genrand() < 0.5){
            delta_nb = 1;
            P_acc = K / (nb + 1.0);
        }
        else {
            delta_nb = - 1;
            P_acc = nb * inv_K;
        }
        if (genrand() < P_acc){
            bonds[ibond] += delta_nb;
            Nb += delta_nb;
            masha = new_masha;
        }
        // shift move -- end
    } // end MC loop

    // print output
    FILE *out_file = fopen("output.dat","w");
    fprintf(out_file, "%4i %.12f %.12f %.12lf %.12lf\n", L, T, K, Z / (nsteps * 1.0), Nb_tot * T / (Z * N * 1.0));
    fclose(out_file);

} //end main


// takes two sites i1,i2 and returns the index of the i1-i2 bond
int bond_number(int i1, int i2, int *x, int *y, int L){
    int x1 = x[i1], x2 = x[i2], y1 = y[i1], y2 = y[i2];
    if (y1 == y2){
        if (x2 == x1 + 1){
            return 2 * i1;
        }
        else if (x1 == x2 + 1){
            return 2 * i2;
        }
        else if (x1 == L - 1){
            return 2 * i1;
        }
        else if (x2 == L - 1){
            return 2 * i2;
        }
    }
    else if (x1 == x2){
        if (y2 == y1 + 1){
            return 2 * i1 + 1;
        }
        else if (y1 == y2 + 1){
            return 2 * i2 + 1;
        }
        else if (y1 == L - 1){
            return 2 * i1 +1;
        }
        else if (y2 == L - 1){
            return 2 * i2 + 1;
        }
    }
    // if you got here, something went wrong
    printf("ERROR!\n");
    return 2 * L * L + 1;
}
