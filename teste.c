int d,r,c,rSize,cSize,addR,addC;
#pragma omp parallel private(i,d,j,r,c,rSize,cSize,addR,addC) num_threads(NUM_THREADS)
{	
	d = 0;
	#pragma omp for schedule(dynamic)
    for (r = 0 ; r < sizeB; r++) {
		addR = r * SUBY_SIZE;
		rSize = strlen(subY[r]);

		for (c = 0; c <= sizeA && d < (sizeB+sizeA-1) ; c++) { 
	
			if (c == sizeA && d < (sizeB-1)) {
				break;
			}

			else if (c == sizeA && d >= (sizeB-1)) {
				#pragma omp barrier
				d++;
				c--;
				continue;
			}

			while (c > (d-r) ) {
				#pragma omp barrier
				d++;
			}

			cSize = strlen(subX[c]);
			addC = c * SUBMAT_SIZE;

			for(i = 0; i < rSize; i++){
				for(j = 0; j < cSize ; j++){
					if (subX[c][j] == subY[r][i]) {
						fTab[i+addR+1][j+addC+1] =
				       			fTab[i+addR][j+addC] + 1;
					}	
					else {
						fTab[i+addR+1][j+addC+1] =
							util_max(fTab[i+addR][j+addC+1],
								fTab[i+addR+1][j+addC]);
					}
				}
			}	
   		}
	}
}