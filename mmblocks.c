#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

#define MSIZE 256
#define NPROC 4
#define GSIZE 2
 
int malloc2dfloat(float ***array, int n, int m) {

    /* allocate the n*m contiguous items */
    float *p = (float *)malloc(n*m*sizeof(float));
    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*array) = (float **)malloc(n*sizeof(float*));
    if (!(*array)) {
       free(p);
       return -1;
    }

    /* set up the pointers into the contiguous memory */
    for (int i=0; i<n; i++)
       (*array)[i] = &(p[i*m]);

    return 0;
}

int free2dfloat(float ***array) {
    /* free the memory - the first element of the array is at the start */
    free(&((*array)[0][0]));

    /* free the pointers into the memory */
    free(*array);

    return 0;
}

int main(int argc, char **argv) {
    float **ma, **mb, **mc, **la, **lb, **lc;

    int rank, size;        // rank of current process and no. of processes
		double t1, t2;  //timing
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (size != NPROC) {
        printf("not enough processors\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    
    if (rank == 0) {
        /* fill in the array, and print it */
	malloc2dfloat(&ma, MSIZE, MSIZE);
	malloc2dfloat(&mb, MSIZE, MSIZE);
	malloc2dfloat(&mc, MSIZE, MSIZE);
	int counter = 1;
        for (int i=0; i<MSIZE; i++) {
            for (int j=0; j<MSIZE; j++){
                ma[i][j] = 0.0 + counter;
								mb[i][j] = 0.0 + counter;
								mc[i][j] = 0.0;
								counter++;
	    }
        }

/*
        printf("ma is:\n");
        for (int i=0; i<MSIZE; i++) {
            for (int j=0; j<MSIZE; j++)
                printf("%.1f ", ma[i][j]);

            printf("\n");
        }
	printf("mb is:\n");
        for (int i=0; i<MSIZE; i++) {
            for (int j=0; j<MSIZE; j++)
                printf("%.1f ", mb[i][j]);

            printf("\n");
        }
	printf("mc is:\n");
        for (int i=0; i<MSIZE; i++) {
            for (int j=0; j<MSIZE; j++)
                printf("%.1f ", mc[i][j]);

            printf("\n");
        }
		*/
    }

    /* create the local array which we'll process */
    malloc2dfloat(&la, MSIZE/GSIZE, MSIZE);
		malloc2dfloat(&lb, MSIZE, MSIZE/GSIZE);
    malloc2dfloat(&lc, MSIZE/GSIZE, MSIZE/GSIZE);
    for(int i = 0; i < MSIZE/GSIZE;i++){
			for(int j = 0;j < MSIZE/GSIZE;j++){
					lc[i][j] = 0;
			}
    }
		for(int i = 0; i < MSIZE/GSIZE;i++){
			for(int j = 0;j < MSIZE;j++){
					la[i][j] = 0;
					lb[j][i] = 0;
			}
    }
		/* now all processors print their local data: */
		/*
    for (int p=0; p<size; p++) {
				//MPI_Barrier(MPI_COMM_WORLD);
        if (rank == p) {
            printf("Local a on rank %d is:\n", rank);
            for (int i=0; i<MSIZE/GSIZE; i++) {
                putchar('|');
                for (int j=0; j<MSIZE; j++) {
									printf("%.1f", la[i][j]);
			//printf("%d", rank);
                    //putchar(local[i][j]);
                }
                printf("| %d\n", rank);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
		
		for (int p=0; p<size; p++) {
				//MPI_Barrier(MPI_COMM_WORLD);
        if (rank == p) {
            printf("Local b on rank %d is:\n", rank);
            for (int i=0; i<MSIZE; i++) {
                putchar('|');
                for (int j=0; j<MSIZE/GSIZE; j++) {
									printf("%.1f", lb[i][j]);
			//printf("%d", rank);
                    //putchar(local[i][j]);
                }
                printf("| %d\n", rank);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
		for (int p=0; p<size; p++) {
				//MPI_Barrier(MPI_COMM_WORLD);
        if (rank == p) {
            printf("Local c on rank %d is:\n", rank);
            for (int i=0; i<MSIZE/GSIZE; i++) {
                putchar('|');
                for (int j=0; j<MSIZE/GSIZE; j++) {
									printf("%.1f", lc[i][j]);
			//printf("%d", rank);
                    //putchar(local[i][j]);
                }
                printf("| %d\n", rank);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
		*/

    /* create a datatype to describe the subarrays of the global array */

    int sizes[2]    = {MSIZE, MSIZE};         /* global size */
    int subsizesla[2] = {MSIZE/GSIZE, MSIZE};     /* local size */
		int subsizeslb[2] = {MSIZE, MSIZE/GSIZE};     /* local size */
		int subsizeslc[2] = {MSIZE/GSIZE, MSIZE/GSIZE};     /* local size */
    int starts[2]   = {0,0};                        /* where this one starts */
    MPI_Datatype typela, typelb, typelc, subarrtypela, subarrtypelb, subarrtypelc;
    
		//Create la type		
		MPI_Type_create_subarray(2, sizes, subsizesla, starts, MPI_ORDER_C, MPI_FLOAT, &typela);
    MPI_Type_create_resized(typela, 0, MSIZE*MSIZE/GSIZE*sizeof(float), &subarrtypela);
    MPI_Type_commit(&subarrtypela);

		//Create lb type		
		MPI_Type_create_subarray(2, sizes, subsizeslb, starts, MPI_ORDER_C, MPI_FLOAT, &typelb);
    MPI_Type_create_resized(typelb, 0, MSIZE/GSIZE*sizeof(float), &subarrtypelb);
    MPI_Type_commit(&subarrtypelb);
	
		//Create lc type		
		MPI_Type_create_subarray(2, sizes, subsizeslc, starts, MPI_ORDER_C, MPI_FLOAT, &typelc);
    MPI_Type_create_resized(typelc, 0, sizeof(float), &subarrtypelc);
    MPI_Type_commit(&subarrtypelc);

		float *aptr=NULL;
    float *bptr=NULL;
    float *cptr=NULL;
    if (rank == 0){
			aptr = &(ma[0][0]);
			bptr = &(mb[0][0]);
			cptr = &(mc[0][0]);
    }

		//get time
		t1 = MPI_Wtime(); 
    /* scatter the array to all processors */
    int sendcounts[NPROC];
    int displsa[NPROC];
		int displsb[NPROC];
		int displsc[NPROC];

    if (rank == 0) {
        for (int i=0; i<NPROC; i++) sendcounts[i] = 1;
        //int disp = 0;
        for (int i=0; i<GSIZE; i++){            
					for(int j=0; j<GSIZE; j++){
            displsa[i*GSIZE+j] = i;
						displsb[i*GSIZE+j] = j;
						displsc[i*GSIZE+j] = i*MSIZE/GSIZE+j*MSIZE*MSIZE/GSIZE;
            //disp += 1;
        	}
    		}
			}

		MPI_Scatterv(aptr, sendcounts, displsa, subarrtypela, &(la[0][0]),
                 MSIZE*MSIZE/GSIZE, MPI_FLOAT,
                 0, MPI_COMM_WORLD);

    MPI_Scatterv(bptr, sendcounts, displsb, subarrtypelb, &(lb[0][0]),
                 MSIZE*MSIZE/GSIZE, MPI_FLOAT,
                 0, MPI_COMM_WORLD);

	
    /* now each processor has its local array, and can process it */
    float tempc = 0.0;
		int rowstart = rank/GSIZE*(MSIZE/GSIZE);
		int colstart = rank%GSIZE*(MSIZE/GSIZE);
    for (int i=0; i<MSIZE/GSIZE; i++) {
        for (int j=0; j<MSIZE/GSIZE; j++) {
					tempc = 0.0;
					for(int k = 0; k < MSIZE; k ++){
						//tempc = tempc + la[i][k+rowstart] * lb[k+colstart][j]; 
						tempc = tempc + la[i][k] * lb[k][j]; 
					}
					lc[i][j] = tempc;
        }
    }

    /* it all goes back to process 0 */
    MPI_Gatherv(&(lc[0][0]),  MSIZE/GSIZE*MSIZE/GSIZE,  MPI_FLOAT,
                 cptr, sendcounts, displsc, subarrtypelc,
                 0, MPI_COMM_WORLD);


		//get time
		t2 = MPI_Wtime(); 

    /* don't need the local data anymore */
		free2dfloat(&la);    
		free2dfloat(&lb);
    free2dfloat(&lc);

    /* or the MPI data type */
    MPI_Type_free(&subarrtypela);
		MPI_Type_free(&subarrtypelb);
		MPI_Type_free(&subarrtypelc);

    if (rank == 0) {
				/*
        printf("Processed grid:\n");
        for (int i=0; i<MSIZE; i++) {
            for (int j=0; j<MSIZE; j++) {
		printf("%.1f", mc[i][j]);
                //putchar(global[i][j]);
            }
            printf("\n");
        }
				*/
				free2dfloat(&ma);
        free2dfloat(&mb);
				free2dfloat(&mc);
    }
    
		if (rank == 0) {
			printf( "Elapsed time is %f\n", t2 - t1 ); 
		}

    MPI_Finalize();

    return 0;
}
