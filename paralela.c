#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#define NUM_THREADS 8

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

char* read_seq(char *fname) {
	//file pointer
	FILE *fseq = NULL;
	//sequence size
	long size = 0;
	//sequence pointer
	char *seq = NULL;
	//sequence index
	int i = 0;

	//open file
	fseq = fopen(fname, "rt");
	if (fseq == NULL ) {
		printf("Error reading file %s\n", fname);
		exit(1);
	}

	//find out sequence size to allocate memory afterwards
	fseek(fseq, 0L, SEEK_END);
	size = ftell(fseq);
	rewind(fseq);

	//allocate memory (sequence)
	seq = (char *) calloc(size + 1, sizeof(char));
	if (seq == NULL ) {
		printf("Erro allocating memory for sequence %s.\n", fname);
		exit(1);
	}

	//read sequence from file
	while (!feof(fseq)) {
		seq[i] = fgetc(fseq);
		if ((seq[i] != '\n') && (seq[i] != EOF))
			i++;
	}
	//insert string terminator
	seq[i] = '\0';

	//close file
	fclose(fseq);

	//return sequence pointer
	return seq;
}

mtype ** allocateScoreMatrix(int sizeA, int sizeB) {
	int i;
	//Allocate memory for LCS score matrix
	mtype ** scoreMatrix = (mtype **) malloc((sizeB + 1) * sizeof(mtype *));
	for (i = 0; i < (sizeB + 1); i++)
		scoreMatrix[i] = (mtype *) malloc((sizeA + 1) * sizeof(mtype));
	return scoreMatrix;
}

void initScoreMatrix(mtype ** scoreMatrix, int sizeA, int sizeB) {
	int i, j;

	//Fill first line of LCS score matrix with zeroes
    omp_set_num_threads(NUM_THREADS);
    
    #pragma omp parallel for
	for (j = 0; j < (sizeA + 1); j++){
		scoreMatrix[0][j] = 0;
    }

	//Do the same for the first collumn
    #pragma omp parallel for
	for (i = 1; i < (sizeB + 1); i++){
		scoreMatrix[i][0] = 0;
    }

}

int LCS(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB) {
	int i, j;

	#pragma omp parallel for collapse(2) private(i,j) shared(scoreMatrix)
	for (i = 1; i < sizeB + 1; i++) {
		for (j = 1; j < sizeA + 1; j++) {
			if (seqA[j - 1] == seqB[i - 1]) {
				#pragma omp critical
				{
					scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
				}
			}
		}
	}

	#pragma omp parallel for collapse(2) private(i,j) shared(scoreMatrix)
	for (i = 1; i < sizeB + 1; i++) {
		for (j = 1; j < sizeA + 1; j++) {
			if (seqA[j - 1] != seqB[i - 1]) {
				scoreMatrix[i][j] = max(scoreMatrix[i-1][j], scoreMatrix[i][j-1]);
			}
		}
	}

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

	char *seqA, *seqB;
	int sizeA, sizeB;

	double start; 
	double end; 
	start = omp_get_wtime(); 

	seqA = read_seq("fileA.in");
	seqB = read_seq("fileB.in");

	sizeA = strlen(seqA);
	sizeB = strlen(seqB);

	mtype ** scoreMatrix = allocateScoreMatrix(sizeA, sizeB);

	initScoreMatrix(scoreMatrix, sizeA, sizeB);

	mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);

	printf("\nScore: %d\n", score);

	freeScoreMatrix(scoreMatrix, sizeB);

	end = omp_get_wtime(); 
	printf("Work took %f seconds\n", end - start);

	return EXIT_SUCCESS;
}
