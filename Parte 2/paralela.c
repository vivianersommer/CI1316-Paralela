/*
	Viviane da Rosa Sommer - GRR20182564
	
    Como rodar:
        mpicc paralela.c -o paralela  
        mpirun -np 2 ./paralela
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

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
	// mtype ** scoreMatrix = (mtype **) malloc((sizeB + 1) * sizeof(mtype *));

	// for (i = 0; i < (sizeB + 1); i++){
	// 	scoreMatrix[i] = (mtype *) malloc((sizeA + 1) * sizeof(mtype));
    // }

	mtype ** scoreMatrix = malloc((sizeB + 1) * sizeof(mtype *));        /*allocating pointers */
	scoreMatrix[0] = malloc((sizeB + 1) * (sizeA + 1) * sizeof(mtype *));  /* allocating data */
	for(i=1; i<(sizeA + 1); i++){
		scoreMatrix[i]=&(scoreMatrix[0][i*(sizeA + 1)]);
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

int LCS(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB, int numero_processos, int rank) {

	int coluna1 = 0, coluna2 = 1, l = 0, destino = 0;

	mtype *A = malloc(sizeof(unsigned char) * (sizeB + 1));
	mtype *B = malloc(sizeof(unsigned char) * (sizeB + 1));

    MPI_Datatype col_matrix;
	MPI_Status status;
    MPI_Type_vector((sizeB + 1), 1, (sizeA + 1), MPI_UNSIGNED_SHORT, &col_matrix);
    MPI_Type_commit(&col_matrix);

	if(rank == 0){

		for(l=0; l< sizeA; l++){
            if(l != rank){ 

				destino = l%(numero_processos - 1) + 1;
				printf("destino = %d\n", destino);

				MPI_Send(&scoreMatrix[0][coluna1], 1, col_matrix, destino, 0, MPI_COMM_WORLD);
				MPI_Send(&scoreMatrix[0][coluna2], 1, col_matrix, destino, 0, MPI_COMM_WORLD);
				MPI_Send(&coluna1, 1, MPI_INT, destino, 0, MPI_COMM_WORLD);
				MPI_Send(&coluna2, 1, MPI_INT, destino, 0, MPI_COMM_WORLD);

				MPI_Recv(&scoreMatrix[0][coluna2], 1, col_matrix, destino, 0, MPI_COMM_WORLD, &status);
				puts("RECEBI!!");
				coluna1++;
				coluna2++;
            }
		}

		printMatrix(seqA, seqB, scoreMatrix, sizeA, sizeB);


	} else {

		int i, j;
        MPI_Recv(&A[0], (sizeB + 1), MPI_UNSIGNED_SHORT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&B[0], (sizeB + 1), MPI_UNSIGNED_SHORT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&coluna1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&coluna2, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

		for (i = 1; i < sizeB + 1; i++) { // percorre pawheae
			if (seqA[coluna2 -1] == seqB[i - 1]) {
				B[i] = scoreMatrix[i - 1][coluna2] + 1;
			} else {
				B[i] = max(B[i-1], A[i]);
			}
		}

		printf("\n --------------------------------- \n");
		for (i = 0; i < (sizeB + 1); i++){
			printf(" %d ", B[i]);
		}
		printf("\n --------------------------------- \n");

		MPI_Send(&B[0], (sizeB + 1), MPI_UNSIGNED_SHORT, 0, 0, MPI_COMM_WORLD);
		puts("ENVIEI!!!");

	}	

	return scoreMatrix[sizeB][sizeA];
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
	int sizeA, sizeB, rank, numero_processos;
	double start, start_read_seq, start_allocateScoreMatrix, end , end_read_seq, end_allocateScoreMatrix; 
    double tempo_total, tempo_sequencial, inicio, fim;	

	seqA = read_seq("sequenciaA.in");
	seqB = read_seq("sequenciaB.in");

	sizeA = strlen(seqA);
	sizeB = strlen(seqB);

	mtype ** scoreMatrix = allocateScoreMatrix(sizeA, sizeB);

	initScoreMatrix(scoreMatrix, sizeA, sizeB);

	// MPI ----------------------------------------------------------------
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numero_processos);
    if(numero_processos != 2)
    {
        printf("Não foi possível abrir 2 processos!\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

	inicio = MPI_Wtime();	

	mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB, numero_processos, rank);

	// printMatrix(seqA, seqB, scoreMatrix, sizeA, sizeB);

	printf("\nScore: %d\n", score);

	fim = MPI_Wtime();
	tempo_total = inicio - fim;
    printf("Tempo total: %f\n", tempo_total);
    
    MPI_Finalize();
	// --------------------------------------------------------------------

	// freeScoreMatrix(scoreMatrix, sizeB);
	
	return EXIT_SUCCESS;
}
