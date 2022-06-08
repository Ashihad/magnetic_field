#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define m 1.0
#define pi atan(1)*4
#define q 1.0
#define B 1.0
#define MODE 0

void pochodne(double t, double *s, double *k) {
    k[0] = s[3]/m;
    k[1] = s[4]/(m*pow(s[0], 2)) - q*B/(2*m);
    k[2] = s[5]/m;
    k[3] = pow(s[4], 2)/(m*pow(s[0], 3)) - pow(q, 2)*pow(B, 2)*s[0]/(4*m);
    k[4] = 0;
    k[5] = 0;
}

void rk4_vec(double t, double dt, int n, double *s, void (*f)(double , double *, double  *)) {
    #define M 1000
    static double k1[M], k2[M], k3[M], k4[M], w[M];
    int i;
    // Kopia tablicy
    for(i = 0; i < n; i++) {
        w[i] = s[i];
    }
    // get rk1
    f(t, w, k1);
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt/2*k1[i];
    }
    // get rk2
    f(t+dt/2, w, k2);
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt/2*k2[i];
    }
    // get rk3
    f(t+dt/2, w, k3);
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt*k3[i];
    }
    // get rk4
    f(t+dt, w, k4);
    // do s[]
    for(i = 0; i < n; i++) {
        s[i]= s[i] + dt/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
}

void changeSystem(double* s, double* s2) {
    s2[0] = s[0]*cos(s[1]);
    s2[1] = s[0]*sin(s[1]);
    s2[2] = s[2];
}

int main() {
    // inicjalizacja  parametrów:
    int n = 6;              // ilość zmiennych w układzie RRZ1
    int N = 5000;           // ilość kroków  czasowych
    double wc = q*B/m;
    double T = 2*pi/wc;
    double dt = 5*T/N;
    double tmax = N*dt;
    double t = 0;

    void (*f)(double, double *, double *);
    f = pochodne;

    double *s;
    s = (double *)malloc(n*sizeof(double));    // tablica  rozwiązań

    double *s2;
    s2 = (double *)malloc(3*sizeof(double));

    FILE *fp;
    FILE *fp2;

    switch(MODE){
        case 0:
            // plik
            
            fp = fopen("wyniki0.dat", "w");
            if(fp == NULL) {
                printf("Error: No file found\n");
                return EXIT_FAILURE;
            }

            // transformacja do układu laboratoryjnego
            fp2 = fopen("laboratoryjny0.dat", "w");
            if(fp == NULL) {
                printf("Error: No file found\n");
                return EXIT_FAILURE;
            }

            // warunki początkowe
            s[0] = 1.5; // r
            s[1] = 1.25*pi; // fi
            s[2] = 0; // z
            s[3] = 0;
            s[4] = q*B*pow(s[0], 2)/2;
            s[5] = 0;
            break;
        case 1:
            // plik
            
            fp = fopen("wyniki1.dat", "w");
            if(fp == NULL) {
                printf("Error: No file found\n");
                return EXIT_FAILURE;
            }

            // transformacja do układu laboratoryjnego
            fp2 = fopen("laboratoryjny1.dat", "w");
            if(fp == NULL) {
                printf("Error: No file found\n");
                return EXIT_FAILURE;
            }

            s[0] = 1.0; // r
            s[1] = 0; // fi
            s[2] = 0; // z
            s[3] = 0;
            s[4] = -q*B*pow(s[0], 2)/2;
            s[5] = 0;
            break;
        case 2:
            // plik
            
            fp = fopen("wyniki2.dat", "w");
            if(fp == NULL) {
                printf("Error: No file found\n");
                return EXIT_FAILURE;
            }

            // transformacja do układu laboratoryjnego
            fp2 = fopen("laboratoryjny2.dat", "w");
            if(fp == NULL) {
                printf("Error: No file found\n");
                return EXIT_FAILURE;
            }
            s[0] = 2.0; // r
            s[1] = 0; // fi
            s[2] = 0; // z
            s[3] = 0;
            s[4] = -q*B*pow(s[0], 2)/2;
            s[5] = 0;
            break;
        case 3:
            // plik
            
            fp = fopen("wyniki3.dat", "w");
            if(fp == NULL) {
                printf("Error: No file found\n");
                return EXIT_FAILURE;
            }

            // transformacja do układu laboratoryjnego
            fp2 = fopen("laboratoryjny3.dat", "w");
            if(fp == NULL) {
                printf("Error: No file found\n");
                return EXIT_FAILURE;
            }
            s[0] = 2.0; // r
            s[1] = 0; // fi
            s[2] = 0; // z
            s[3] = 2;
            s[4] = -q*B*pow(s[0], 2)/2;
            s[5] = 0;
            break;
    }
    
    // headers
    fprintf(fp, "%s %s %s %s %s %s %s\n", "czas", "r", "fi", "z", "pr", "pfi", "pz");
    fprintf(fp2, "%s %s %s\n", "x", "y", "z");
    // simulation loop
    for(size_t i = 1; i <= N; i++) {
        rk4_vec(t, dt, n, s, f);
        // double first = (pow(s[3], 2) + pow(s[4], 2)/pow(s[0], 2) + pow(s[5], 2))/(2.0*m);
        // double second = -q*B*s[4]/(2.0*m);
        // double third = 
        double e = 1/(2.0*m)*(pow(s[3], 2) + pow(s[4], 2)/pow(s[0], 2) + pow(s[5], 2)) - q*B*s[4]/(2.0*m) + pow(q, 2)*pow(B, 2)*pow(s[0], 2)/(8.0*m);
        t=t+dt;
        fprintf(fp, "%f %f %f %f %f %f %f %f\n", t, s[0], s[1], s[2], s[3], s[4], s[5], e);
        changeSystem(s, s2);
        fprintf(fp2, "%f %f %f\n", s2[0], s2[1], s2[2]);
    }
    fclose(fp);
    fclose(fp2);
    return EXIT_SUCCESS;
}