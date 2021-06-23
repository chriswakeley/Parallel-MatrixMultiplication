#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

#define MSIZE 144
#define NPROC 1
 
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
		double t1, t2;  //timing
    float **ma, **mb, **mc, **lb, **lc;

    int rank, size;        // rank of current process and no. of processes

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (size != NPROC) {
        printf("not enough processors\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    malloc2dfloat(&ma, MSIZE, MSIZE);
    if (rank == 0) {
        /* fill in the array, and print it */
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
        }*/
    }

    /* create the local array which we'll process */
    malloc2dfloat(&lb, MSIZE, MSIZE/NPROC);
    malloc2dfloat(&lc, MSIZE, MSIZE/NPROC);
    for(int i = 0; i < MSIZE;i++){
	for(int j = 0;j < MSIZE/NPROC;j++){
	    lc[i][j] = 0;
	}
    }
/*
    for (int p=0; p<size; p++) {
        if (rank == p) {
            printf("Local process on rank %d is:\n", rank);
            for (int i=0; i<MSIZE; i++) {
                putchar('|');
                for (int j=0; j<MSIZE/NPROC; j++) {
			printf("%.1f", lc[i][j]);
                    //putchar(local[i][j]);
                }
                printf("|\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
*/
    /* create a datatype to describe the subarrays of the global array */

    int sizes[2]    = {MSIZE, MSIZE};         /* global size */
    int subsizes[2] = {MSIZE, MSIZE/NPROC};     /* local size */
    int starts[2]   = {0,0};                        /* where this one starts */
    MPI_Datatype type, subarrtype;
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_FLOAT, &type);
    MPI_Type_create_resized(type, 0, MSIZE/NPROC*sizeof(float), &subarrtype);
    MPI_Type_commit(&subarrtype);

    float *bptr=NULL;
    float *cptr=NULL;
    if (rank == 0){
	bptr = &(mb[0][0]);
	cptr = &(mc[0][0]);
    }

    /* scatter the array to all processors */
		//get time
		t1 = MPI_Wtime(); 
    int sendcounts[NPROC];
    int displs[NPROC];

    if (rank == 0) {
        for (int i=0; i<NPROC; i++) sendcounts[i] = 1;
        int disp = 0;
        for (int i=0; i<NPROC; i++) {            
            displs[i] = disp;
            disp += 1;
        }
    }


    MPI_Scatterv(bptr, sendcounts, displs, subarrtype, &(lb[0][0]),
                 MSIZE*MSIZE/NPROC, MPI_FLOAT,
                 0, MPI_COMM_WORLD);

    /* broadcast a*/
    MPI_Bcast (&(ma[0][0]), MSIZE*MSIZE, MPI_FLOAT, 0, MPI_COMM_WORLD);

    /* now all processors print their local data: */
	/*
    for (int p=0; p<size; p++) {
        if (rank == p) {
            printf("Local process on rank %d is:\n", rank);
            for (int i=0; i<MSIZE; i++) {
                putchar('|');
                for (int j=0; j<MSIZE/NPROC; j++) {
			printf("%.1f", lb[i][j]);
			//printf("%d", rank);
                    //putchar(local[i][j]);
                }
                printf("|\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
*/
	
    /* now each processor has its local array, and can process it */
    float tempc = 0.0;
    for (int i=0; i<MSIZE; i++) {
        for (int j=0; j<MSIZE/NPROC; j++) {
	    tempc = 0.0;
	    for(int k = 0; k < MSIZE; k ++){
		tempc = tempc + ma[i][k] * lb[k][j]; 
	    }
	    lc[i][j] = tempc;
        }
    }
    
    /* now all processors print their local data: */
	/*
    for (int p=0; p<size; p++) {
        if (rank == p) {
            printf("Local process on rank %d is:\n", rank);
            for (int i=0; i<MSIZE; i++) {
                putchar('|');
                for (int j=0; j<MSIZE/NPROC; j++) {
			printf("%.1f", lc[i][j]);
                    //putchar(local[i][j]);
                }
                printf("|\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
	*/
    /* it all goes back to process 0 */
    MPI_Gatherv(&(lc[0][0]),  MSIZE*MSIZE/NPROC,  MPI_FLOAT,
                 cptr, sendcounts, displs, subarrtype,
                 0, MPI_COMM_WORLD);

		//get time
		t2 = MPI_Wtime(); 
    /* don't need the local data anymore */
    free2dfloat(&lb);
    free2dfloat(&lc);

    /* or the MPI data type */
    MPI_Type_free(&subarrtype);
		
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


        free2dfloat(&mb);
	free2dfloat(&mc);
    }
    free2dfloat(&ma);

if (rank == 0) {
			printf( "Elapsed time is %f\n", t2 - t1 ); 
		}
    MPI_Finalize();

    return 0;
}
