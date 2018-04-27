// srad.cpp : Defines the entry point for the console application.
//

//#define OUTPUT

//#define	ITERATION
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void random_matrix(float *I, int rows, int cols);

void usage(int argc, char **argv)
{
	fprintf(stderr, "Usage: %s <rows> <cols> <y1> <y2> <x1> <x2> <no. of threads><lamda> <no. of iter>\n", argv[0]);
	fprintf(stderr, "\t<rows>\t- number of rows\n");
	fprintf(stderr, "\t<cols>\t- number of cols\n");
	fprintf(stderr, "\t<y1>\t- y1 value of the speckle\n");
	fprintf(stderr, "\t<y2>\t- y2 value of the speckle\n");
	fprintf(stderr, "\t<x1>\t- x1 value of the speckle\n");
	fprintf(stderr, "\t<x2>\t- x2 value of the speckle\n");
	fprintf(stderr, "\t<no. of threads>  - no. of threads\n");
	fprintf(stderr, "\t<lamda>   - lambda (0,1)\n");
	fprintf(stderr, "\t<no. of iter>   - number of iterations\n");
}

int main(int argc, char* argv[])
{   
	int rows, cols, size_I, size_R, niter = 10, iter, k;
    float *J, q0sqr, sum, sum2, tmp, meanROI,varROI ;
	float Jc, G2, L, num, den, qsqr;
	float *dN,*dS,*dW,*dE;
	int r1, r2, c1, c2;
	float cN,cS,cW,cE;
	float *c, D;
	float lambda;
	int i, j, iN_value, iS_value;
    int nthreads;

	if (argc == 10)
	{
		rows = atoi(argv[1]); //number of rows in the domain
		cols = atoi(argv[2]); //number of cols in the domain
		if ((rows%16!=0) || (cols%16!=0)){
			fprintf(stderr, "rows and cols must be multiples of 16\n");
			exit(1);
		}
		r1   = atoi(argv[3]); //y1 position of the speckle
		r2   = atoi(argv[4]); //y2 position of the speckle
		c1   = atoi(argv[5]); //x1 position of the speckle
		c2   = atoi(argv[6]); //x2 position of the speckle
		nthreads = atoi(argv[7]); // number of threads
		lambda = atof(argv[8]); //Lambda value
#ifdef ITERATION
		niter = atoi(argv[9]); //number of iterations
#endif //ITERATION
	}
    else{
		usage(argc, argv);
		exit(1);
    }

	size_I = cols * rows;
    size_R = (r2-r1+1)*(c2-c1+1);   

    J = (float *)malloc( size_I * sizeof(float) );
    if (J == NULL) {
		printf("malloc failed.\n");
		exit(1);
	}
	c  = (float *)malloc(sizeof(float)* size_I) ;
  	if (c == NULL) {
		printf("malloc failed.\n");
		exit(1);
	}

    int iN[rows];
    int iS[rows];
    int jW[rows];
    int jE[rows];
    
	dN = (float *)malloc(sizeof(float)* size_I) ;
  	if (dN == NULL) {
		printf("malloc failed.\n");
		exit(1);
	}
    dS = (float *)malloc(sizeof(float)* size_I) ;
  	if (dS == NULL) {
		printf("malloc failed.\n");
		exit(1);
	}
    dW = (float *)malloc(sizeof(float)* size_I) ;
  	if (dW == NULL) {
		printf("malloc failed.\n");
		exit(1);
	}
    dE = (float *)malloc(sizeof(float)* size_I) ;    
  	if (dE == NULL) {
		printf("malloc failed.\n");
		exit(1);
	}

    for (i=0; i< rows; i++) {
        iN[i] = i-1;
        iS[i] = i+1;
    }    
    for (int j=0; j< cols; j++) {
        jW[j] = j-1;
        jE[j] = j+1;
    }
    iN[0]    = 0;
    iS[rows-1] = rows-1;
    jW[0]    = 0;
    jE[cols-1] = cols-1;
	
	printf("Randomizing the input matrix\n");

    random_matrix(J, rows, cols);

    for (k = 0;  k < size_I; k++ ) {
     	J[k] = (float)exp(J[k]) ;
    }
   
	printf("Start the SRAD main loop\n");

#ifdef ITERATION
	for (iter=0; iter< niter; iter++) {
#endif        
		sum=0; sum2=0;     
		for (i=r1; i<=r2; i++) {
            for (j=c1; j<=c2; j++) {
            	tmp   = J[i * cols + j];
                sum  += tmp;
                sum2 += tmp*tmp;
            }
        }
        meanROI = sum / size_R;
        varROI  = (sum2 / size_R) - meanROI*meanROI;
        q0sqr   = varROI / (meanROI*meanROI);

		for (i = 0 ; i < rows ; i++) {
			for (iN_value = iN[i], iS_value = iS[i], j = 0; j < cols; j++) { 
		
				k = i * cols + j;
				Jc = J[k];
 				
				// directional derivates
                dN[k] = J[iN_value * cols + j] - Jc;
                dS[k] = J[iS_value * cols + j] - Jc;
                dW[k] = J[i * cols + jW[j]] - Jc;
                dE[k] = J[i * cols + jE[j]] - Jc;
			
                G2 = (dN[k]*dN[k] + dS[k]*dS[k] + dW[k]*dW[k] + dE[k]*dE[k]) / (Jc*Jc);

   		        L = (dN[k] + dS[k] + dW[k] + dE[k]) / Jc;

				num  = (0.5*G2) - ((1.0/16.0)*(L*L)) ;
                den  = 1 + (.25*L);
                qsqr = num/(den*den);
 				
                // diffusion coefficent (equ 33)
                den = (qsqr-q0sqr) / (q0sqr * (1+q0sqr)) ;

                c[k] = 1.0 / (1.0+den) ;
                  
                // saturate diffusion coefficent
                if (c[k] < 0)
                	c[k] = 0;
                else if (c[k] > 1) 
                	c[k] = 1;
   			}
      	} 

		for (i = 0; i < rows; i++) {
            for (iS_value = iS[i], j = 0; j < cols; j++) {        

                // current index
                k = i * cols + j;
                
                // diffusion coefficent
				cN = c[k];
				cS = c[iS_value * cols + j];
				cW = c[k];
				cE = c[i * cols + jE[j]];

                // divergence (equ 58)
                D = cN * dN[k] + cS * dS[k] + cW * dW[k] + cE * dE[k];
                
                // image update (equ 61)
                J[k] = J[k] + 0.25*lambda*D;
#ifdef OUTPUT
                printf("%.5f ", J[k]); 
#endif //output
            }
#ifdef OUTPUT
                printf("\n"); 
#endif //output
	    }

#ifdef ITERATION
	}
#endif

#ifdef OUTPUT
	for(i = 0 ; i < rows ; i++) {
		for (j = 0 ; j < cols ; j++) {
			printf("%.5f ", J[i * cols + j]); 
    	}
        printf("\n"); 
    }
#endif 
	printf("Computation Done\n");
	free(J);
    free(dN); free(dS); free(dW); free(dE);
	free(c);
	return 0;
}

void random_matrix(float *I, int rows, int cols) {
	int i, j;
	srand(7);
	for(i = 0 ; i < rows ; i++) {
		for (j = 0 ; j < cols ; j++) {
			I[i * cols + j] = rand()/(float)RAND_MAX ;
#ifdef OUTPUT
			printf("%.2f ", I[i * cols + j]); 
#endif 
		}
#ifdef OUTPUT
		printf("\n"); 
#endif 
	}
}