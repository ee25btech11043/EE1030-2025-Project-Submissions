#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double **new(int r, int c) {
    
    double **m = malloc(r * sizeof(double *));
    for (int i = 0; i < r; i++)
        
        m[i] = calloc(c, sizeof(double));
    
    return m;
}

void f(double **m, int r) {
    for (int i = 0; i < r; i++)
        
        free(m[i]);
    free(m);
}

// transpose
double **transpose(double **a, int r, int c) {
    double **t = new(c, r);
    
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            
            t[j][i] = a[i][j];
    return t;
}

//product of two matrix
double **product(double **A, double **B, int r, int p, int c) {

    
    double **R = new(r, c);
    
    for (int i = 0; i < r; i++) {
        
        for (int j = 0; j < c; j++) {
            double s = 0;
            
            for (int k = 0; k < p; k++)
                s += A[i][k] * B[k][j];
            R[i][j] = s;
        }
    }
    return R;
}

void jacobi(double **M, int n, double *eval, double **V) {
    for (int i = 0; i < n; i++)
        
        for (int j = 0; j < n; j++)
            
            V[i][j] = (i == j);

    for (int it = 0; it < 100; it++) {
        
        double mx = 0;
        
        int p = -2;
        int q = -2;
        for (int i = 0; i < n; i++)
            
            for (int j = i + 1; j < n; j++)
                
                if (fabs(M[i][j]) > mx) {
                    mx = fabs(M[i][j]);
                    p = i;
                    q = j;
                }

        if (mx < 1e-9)
         break;

        double bpp = M[p][p];
        double bqq = M[q][q];
        
        double bpq = M[p][q];
        
        //finding angle
        double ang = 0.5 * atan2(2 * bpq, (bqq - bpp));
            double c = cos(ang);
        double s = sin(ang);

        for (int i = 0; i < n; i++) {
             if (i == p || i == q) 
            continue;
            double mip = M[i][p];
            double miq = M[i][q];
            
              M[i][p] = c * mip + s * miq;
            
              M[p][i] = M[i][p];
            M[i][q] = -s * mip + c * miq;
              M[q][i] = M[i][q];
        }

        M[p][p] = c * c * bpp + 2 * s * c * bpq + s * s * bqq;
           M[q][q] = s * s * bpp - 2 * s * c * bpq + c * c * bqq;
        M[p][q] = M[q][p] = 0;

          for (int i = 0; i < n; i++) {
            double vip = V[i][p], viq = V[i][q];
              
            
            V[i][p] = c * vip + s * viq;
            
            V[i][q] = -s * vip + c * viq;
        }
    }

    for (int i = 0; i < n; i++)
        eval[i] = M[i][i];
}

void sort(double *eval, double **V, int n) {
    for (int i = 0; i < n - 1; i++)
  for (int j = i + 1; j < n; j++)
            if (eval[j] > eval[i]) {
                
                double t = eval[i];
                eval[i] = eval[j];
                
                
                eval[j] = t;
                for (int k = 0; k < n; k++) {
                    double tmp = V[k][i];
                    
                    V[k][i] = V[k][j];
                    
                    V[k][j] = tmp;
                }
            }
}

//for genrating U matrix
double **Umat(double **A, int r, int c, double **V, double *s) {
    
    double **U = new(r, c);
    for (int j = 0; j < c; j++) {
        for (int i = 0; i < r; i++) {
            
            double sm = 0;
            for (int t = 0; t < c; t++)
                
                sm = sm + A[i][t] * V[t][j];
            
            U[i][j] = (s[j] > 1e-10) ? sm / s[j] : 0;
        }
    }
    return U;
}

//making Ak matrix
double **form(double **U, double *s, double **V, int r, int c, int k) {
    
    double **R = new(r, c);
    
    for (int t = 0; t < k; t++) {
        if (s[t] == 0) 
        continue;
        
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                R[i][j] += s[t] * U[i][t] * V[j][t];
    }
    return R;
}


double **read(char *name, int *r, int *c, int *maxv) {
    
    FILE *fp = fopen(name, "r");
    if (!fp) {
        printf("Cannot open file %s\n", name);
        
        return 0;
    }
    
    char ty[3];
    fscanf(fp, "%2s", ty);

    int ch = fgetc(fp);
    while (ch == '#') {
        
        while (fgetc(fp) != '\n');
        ch = fgetc(fp);
    }

    ungetc(ch, fp);

    fscanf(fp, "%d %d %d", c, r, maxv);

    
    double **A = new(*r, *c);
    
    for (int i = 0; i < *r; i++)
        
        for (int j = 0; j < *c; j++) {
            int v;
            fscanf(fp, "%d", &v);
            A[i][j] = v;
        }
    fclose(fp);
    
    return A;
}

//for U matrix
void generate(char *name, double **A, int r, int c, int maxv) {
    
    FILE *fp = fopen(name, "w");
    
    fprintf(fp, "P2\n%d %d\n%d\n", c, r, maxv);
    
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            int v = (int)A[i][j];
            if (v < 0) 
            v = 0;
            if (v > maxv) 
                
            v = maxv;
            fprintf(fp, "%d ", v);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

int main() {
    int r;
    int c;
    int maxp;

    char *name = "einstein.pgm";

    double **A = read(name, &r, &c, &maxp);


    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            
            A[i][j] /= maxp;

    double **At = transpose(A, r, c);


    double **AtA = product(At, A, c, r, c);

    double *eval = malloc(c * sizeof(double));


    double **V = new(c, c);

    jacobi(AtA, c, eval, V);
    sort(eval, V, c);

    double *s = malloc(c * sizeof(double));

        
    for (int i = 0; i < c; i++) {
        
    if (eval[i] > 0)
        s[i] = sqrt(eval[i]);
    else
        s[i] = 0;
}

    double **U = Umat(A, r, c, V, s);

    int k;
    printf("Enter k: ");
        
    scanf("%d", &k);
    if (k < 1) 
        
    k = 1;
    if (k > c) 
    k = c;

    double **A2 = form(U, s, V, r, c, k);

        
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++) {
            A2[i][j] *= maxp;
            if (A2[i][j] < 0)
             A2[i][j] = 0;
            if (A2[i][j] > maxp) 
            A2[i][j] = maxp;
        }


        
    double diff = 0;
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++) {
            
            double d = A[i][j] * maxp - A2[i][j];
            
            diff += d * d;
        }
    printf("Frobenius error = %.3f\n", sqrt(diff));

    generate("out.pgm", A2, r, c, maxp);
        
    f(A, r);
    f(U, r);
        
    f(At, c);
        
    f(AtA, c);
    f(V, c);
        
    f(A2, r);
    free(eval);
        
    free(s);
    return 0;
}
