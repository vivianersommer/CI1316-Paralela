/*
	Viviane da Rosa Sommer - GRR20182564
	
    Como rodar:
        mpicc paralela.c -o paralela -O3
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

	mtype *mem = malloc((sizeB + 1) * (sizeA + 1) * sizeof(mtype));
	mtype **scoreMatrix = malloc((sizeB + 1) * sizeof(mtype*));
	scoreMatrix[0] = mem;

	for(int i = 1; i < (sizeB + 1); i++){
		scoreMatrix[i] = scoreMatrix[i-1] + (sizeA + 1);
	}

	return scoreMatrix;
}

void initScoreMatrix(mtype ** scoreMatrix, int sizeA, int sizeB) {

	for (int j = 0; j < (sizeA + 1); j++){
		scoreMatrix[0][j] = 0;
    }

	for (int i = 1; i < (sizeB + 1); i++){
		scoreMatrix[i][0] = 0;
    }
}

int LCS(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB, int numero_processos, int rank) {

	int coluna1 = 0, coluna2 = 1, l = 0, destino = 0;

	mtype result;

    MPI_Datatype col_matrix;
	MPI_Status status;
    MPI_Type_vector((sizeB + 1), 1, (sizeA + 1), MPI_UNSIGNED_SHORT, &col_matrix);
    MPI_Type_commit(&col_matrix);

	if(rank == 0){

		for(l=0; l< (sizeA + 1); l++){
            if(l != rank){ 

				destino = (l%(numero_processos - 1) == 0) ? (numero_processos - 1): (l%(numero_processos - 1));

				MPI_Send(&scoreMatrix[0][coluna1], 1, col_matrix, destino, 0, MPI_COMM_WORLD);
				MPI_Send(&scoreMatrix[0][coluna2], 1, col_matrix, destino, 0, MPI_COMM_WORLD);
				MPI_Send(&coluna1, 1, MPI_INT, destino, 0, MPI_COMM_WORLD);
				MPI_Send(&coluna2, 1, MPI_INT, destino, 0, MPI_COMM_WORLD);

				MPI_Recv(&scoreMatrix[0][coluna2], 1, col_matrix, destino, 0, MPI_COMM_WORLD, &status);

				coluna1++;
				coluna2++;
            }
		}

	} else {

		mtype *A = (mtype *) malloc(sizeof(mtype) * (sizeB + 1));
		mtype *B = (mtype *) malloc(sizeof(mtype) * (sizeB + 1));

		int termina_mod = sizeA % (numero_processos - 1);
		int termina_div = sizeA / (numero_processos - 1);
		int fim = 0;

		if(termina_mod == 0){
			fim = termina_div;
		} else {
			if (termina_mod >= rank){
				fim = termina_div + termina_mod;
			} else {
				fim = termina_div;
			}
		}

		do{

			fim--;

			int i, j;
			MPI_Recv(&A[0], (sizeB + 1), MPI_UNSIGNED_SHORT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&B[0], (sizeB + 1), MPI_UNSIGNED_SHORT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&coluna1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&coluna2, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

			for (i = 1; i < sizeB + 1; i++) { // LCS
				if (seqA[coluna2 -1] == seqB[i - 1]) {
					B[i] = A[i - 1] + 1;
				} else {
					B[i] = max(B[i-1], A[i]);
				}
			}

			MPI_Ssend(&B[0], (sizeB + 1), MPI_UNSIGNED_SHORT, 0, 0, MPI_COMM_WORLD);
			
		} while (fim != 0);
	}	

	return scoreMatrix[sizeB][sizeA];
}

int main(int argc, char ** argv) {
	

	char *seqA, *seqB;
	int sizeA, sizeB, rank, numero_processos;
	double inicio_serial, fim_serial, inicio_paralela, fim_paralela, inicio_total, fim_total;

	MPI_Init(&argc, &argv);

	inicio_total = MPI_Wtime();

	// SERIAL --------------------------------------------------------------------------------
	inicio_serial = MPI_Wtime();	

		seqA = read_seq("sequenciaA.in");
		seqB = read_seq("sequenciaB.in");
		sizeA = strlen(seqA);
		sizeB = strlen(seqB);
		mtype ** scoreMatrix = allocateScoreMatrix(sizeA, sizeB);
		initScoreMatrix(scoreMatrix, sizeA, sizeB);

	fim_serial = MPI_Wtime() - inicio_serial;
	//----------------------------------------------------------------------------------------

	// MPI -----------------------------------------------------------------------------------
	inicio_paralela = MPI_Wtime();	

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &numero_processos);

		mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB, numero_processos, rank);

	fim_paralela = MPI_Wtime() - inicio_paralela;
	//----------------------------------------------------------------------------------------

	fim_total= MPI_Wtime() - inicio_total;

	if (rank == 0){
		printf("%f\n", fim_total);
	}

	MPI_Finalize();
	
	return EXIT_SUCCESS;
}
