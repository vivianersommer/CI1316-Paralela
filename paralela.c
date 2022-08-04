#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
#include <math.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#define NUM_THREADS 2

int SUBY_SIZE = 500000;
typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

char** read_seq(char *fname, int type, int *size_temp) {

	FILE *fseq = NULL;
	char **seq = NULL, *line = NULL;
	int i = 0, size = 0;
    size_t lnlen;

	fseq = fopen(fname, "rt");
	if (fseq == NULL ) {
		printf("Error reading file %s\n", fname);
		exit(1);
	}

    size = getline(&line,&lnlen,fseq);

    if (line[size-1] == '\n')
		size--;


	int SUBMAT_SIZE = ceil((double)size / NUM_THREADS); 
	int num_blocks, type_size;

    if (type == 0) {
	    num_blocks = ceil((double)size / SUBMAT_SIZE); 
        type_size = SUBMAT_SIZE;
    } else {
	    num_blocks = ceil((double)size / SUBY_SIZE); 
        type_size = SUBY_SIZE;
    }

	seq = (char **)malloc(sizeof(char *) * num_blocks);

	for (int b = 0; b < num_blocks; b++) {
		seq[b] = (char *)malloc(sizeof(char) * (type_size+1));
		if (((b+1) * type_size) <= size) {
			strncpy(seq[b],line+(b*type_size),type_size);
			seq[b][type_size] = '\0';
		}
		else {
			strncpy(seq[b],line+(b*type_size),
					size - b*type_size);
			seq[b][size - b*type_size] = '\0';
		}
	}

	fclose(fseq);

    *size_temp = size;

	return seq;
}

mtype ** allocateScoreMatrix(int sizeA, int sizeB) {
	int i;
	//Allocate memory for LCS score matrix
	mtype ** scoreMatrix = (mtype **) malloc((sizeB + 1) * sizeof(mtype *));
	for (i = 0; i <= sizeB; i++)
		scoreMatrix[i] = (mtype *) malloc((sizeA + 1) * sizeof(mtype));
	return scoreMatrix;
}

void initScoreMatrix(mtype ** scoreMatrix, int sizeA, int sizeB) {
	int i, j;

	//Fill first line of LCS score matrix with zeroes    
    #pragma omp parallel for
	for (j = 0; j <= sizeA; j++){
		scoreMatrix[0][j] = 0;
    }

	//Do the same for the first collumn
    #pragma omp parallel for
	for (i = 1; i <= sizeB; i++){
		scoreMatrix[i][0] = 0;
    }

}

void printMatrix(char * seqA, char * seqB, mtype ** scoreMatrix, int sizeA,
		int sizeB) {
	int i, j;

	printf("Score Matrix:\n");
	printf("========================================\n");

	printf("    ");
	printf("%5c   ", ' ');

	for (j = 0; j < sizeA; j++)
		printf("%5c   ", seqA[j]);
	printf("\n");
	for (i = 0; i < sizeB + 1; i++) {
		if (i == 0)
			printf("    ");
		else
			printf("%c   ", seqB[i - 1]);
		for (j = 0; j < sizeA + 1; j++) {
			printf("%5d   ", scoreMatrix[i][j]);
		}
		printf("\n");
	}
	printf("========================================\n");
}

int LCS(mtype ** scoreMatrix, int sizeA, int sizeB, char **seqA, char **seqB) {

	int d,r,c,rSize,cSize,addR,addC, i,j;
	int SUBMAT_SIZE = ceil((double)sizeA / NUM_THREADS); 
	int num_blocksY = ceil((double)sizeB / SUBY_SIZE); 
	int num_blocksX = ceil((double)sizeA / SUBMAT_SIZE);  
 
	#pragma omp parallel private(i,d,j,r,c,rSize,cSize,addR,addC) num_threads(NUM_THREADS)
	{
		d = 0;

		#pragma omp for schedule(dynamic)
    	for (r = 0 ; r < num_blocksY; r++) {
            addR = r * SUBY_SIZE;
            rSize = strlen(seqB[r]);

            for (c = 0; c <= num_blocksX && d < (num_blocksY+num_blocksX-1) ; c++) { 
        
                if (c == num_blocksX && d < (num_blocksY-1)) {
                    break;
                } else if (c == num_blocksX && d >= (num_blocksY-1)) {
                    #pragma omp critical 
                    d++;
                    c--;
                    continue;
                }

                while (c > (d-r) ) {
                    #pragma omp critical
                    d++;
                }

                cSize = strlen(seqA[c]);
                addC = c * SUBMAT_SIZE;

                for (i = 0; i < rSize; i++) {
                    for (j = 0; j < cSize ; j++) {

                        if (seqA[c][j] == seqB[r][i]) {
                            scoreMatrix[i+addR+1][j+addC+1] =
                                    scoreMatrix[i+addR][j+addC] + 1;
                        }	
                        else {
                            scoreMatrix[i+addR+1][j+addC+1] =
                                max(scoreMatrix[i+addR][j+addC+1],
                                    scoreMatrix[i+addR+1][j+addC]);
                        }
                    }
                }	
            }
        }
	}

	// for ( i = 0; i <= sizeB; i++) {
	// 	for ( j = 0; j <= sizeA; j++) {
	// 		printf("%d\t",scoreMatrix[i][j]);
	// 	}
	// 	printf("\n");
	// }


	return scoreMatrix[sizeB][sizeA];
}

void freeScoreMatrix(mtype **scoreMatrix, int sizeB) {
	int i;
	#pragma omp parallel for
	for (i = 0; i < (sizeB + 1); i++)
		free(scoreMatrix[i]);
	free(scoreMatrix);
}

int main(int argc, char ** argv) {

	char **seqA, **seqB;
	int sizeA, sizeB;

	double start, start_read_seq, start_allocateScoreMatrix; 
	double end , end_read_seq, end_allocateScoreMatrix; 
	start = omp_get_wtime(); 
	
	omp_set_num_threads(NUM_THREADS);

	start_read_seq = omp_get_wtime();
	seqA = read_seq("sequenciaA.in", 0, &sizeA);
	seqB = read_seq("sequenciaB.in", 1, &sizeB);
	end_read_seq = omp_get_wtime();

	start_allocateScoreMatrix = omp_get_wtime();
	mtype ** scoreMatrix = allocateScoreMatrix(sizeA, sizeB);
	end_allocateScoreMatrix = omp_get_wtime();

	initScoreMatrix(scoreMatrix, sizeA, sizeB);

	mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);

	// printf("\nScore: %d\n", score);

	freeScoreMatrix(scoreMatrix, sizeB);

	end = omp_get_wtime(); 

	printf("%f\n", end - start);
	// printf("read_seq took %f seconds\n", end_read_seq - start_read_seq);
	// printf("allocateScoreMatrix took %f seconds\n", end_allocateScoreMatrix - start_allocateScoreMatrix);

	return EXIT_SUCCESS;
}
