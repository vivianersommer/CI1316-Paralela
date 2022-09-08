#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef unsigned short mtype;


char* read_seq(char *fname) {
	
	FILE *fseq = NULL;
	long size = 0;
	char *seq = NULL;
	int i = 0;

	fseq = fopen(fname, "rt");
	if (fseq == NULL ) {
		printf("Error reading file %s\n", fname);
		exit(1);
	}

	fseek(fseq, 0L, SEEK_END);
	size = ftell(fseq);
	rewind(fseq);
	seq = (char *) calloc(size + 1, sizeof(char));
    
	if (seq == NULL ) {
		printf("Erro allocating memory for sequence %s.\n", fname);
		exit(1);
	}

	while (!feof(fseq)) {
		seq[i] = fgetc(fseq);
		if ((seq[i] != '\n') && (seq[i] != EOF)){
			i++;
        }
	}

	seq[i] = '\0';
	fclose(fseq);

	return seq;
}

mtype ** allocateScoreMatrix(int sizeA, int sizeB) {

	int i;
	mtype ** scoreMatrix = (mtype **) malloc((sizeB + 1) * sizeof(mtype *));

	for (i = 0; i < (sizeB + 1); i++){
		scoreMatrix[i] = (mtype *) malloc((sizeA + 1) * sizeof(mtype));
    }

	return scoreMatrix;
}

void initScoreMatrix(mtype ** scoreMatrix, int sizeA, int sizeB) {

	int i, j;

	for (j = 0; j < (sizeA + 1); j++){
		scoreMatrix[0][j] = 0;
    }

	for (i = 1; i < (sizeB + 1); i++){
		scoreMatrix[i][0] = 0;
    }
}

int LCS(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB) {

	int i, j;

	for (i = 1; i < sizeB + 1; i++) {
		for (j = 1; j < sizeA + 1; j++) {
			if (seqA[j - 1] == seqB[i - 1]) {
				scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
			} else {
				scoreMatrix[i][j] =max(scoreMatrix[i-1][j], scoreMatrix[i][j-1]);
			}
		}
	}

	return scoreMatrix[sizeB][sizeA];
}

void printMatrix(char * seqA, char * seqB, mtype ** scoreMatrix, int sizeA,int sizeB) {
	
    int i, j;
	printf("Score Matrix:\n");
	printf("========================================\n");
	printf("    ");
	printf("%5c   ", ' ');

	for (j = 0; j < sizeA; j++){
		printf("%5c   ", seqA[j]);
    }

	printf("\n");

	for (i = 0; i < sizeB + 1; i++) {
		if (i == 0){
			printf("    ");
        }
		else{
			printf("%c   ", seqB[i - 1]);
        }
		for (j = 0; j < sizeA + 1; j++) {
			printf("%5d   ", scoreMatrix[i][j]);
		}
		printf("\n");
	}

	printf("========================================\n");
}

void freeScoreMatrix(mtype **scoreMatrix, int sizeB) {

	int i;

	for (i = 0; i < (sizeB + 1); i++){
		free(scoreMatrix[i]);
    }

	free(scoreMatrix);
}

int main(int argc, char ** argv) {

	char *seqA, *seqB;
	int sizeA, sizeB;
	double start, start_read_seq, start_allocateScoreMatrix, end , end_read_seq, end_allocateScoreMatrix; 

	start = omp_get_wtime(); 
	omp_set_num_threads(1);

	seqA = read_seq("sequenciaA.in");
	seqB = read_seq("sequenciaB.in");

	sizeA = strlen(seqA);
	sizeB = strlen(seqB);

	mtype ** scoreMatrix = allocateScoreMatrix(sizeA, sizeB);

	initScoreMatrix(scoreMatrix, sizeA, sizeB);

	mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);

	// printMatrix(seqA, seqB, scoreMatrix, sizeA, sizeB);

	// printf("\nScore: %d\n", score);

	freeScoreMatrix(scoreMatrix, sizeB);

	end = omp_get_wtime(); 
	
	return EXIT_SUCCESS;
}
