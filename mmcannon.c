#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

#define MSIZE 144
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
		MPI_Status status; 
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
    malloc2dfloat(&la, MSIZE/GSIZE, MSIZE/GSIZE);
		malloc2dfloat(&lb, MSIZE/GSIZE, MSIZE/GSIZE);
    malloc2dfloat(&lc, MSIZE/GSIZE, MSIZE/GSIZE);
		for(int i = 0; i < GSIZE;i++){
			for(int j = 0; j < GSIZE; j++){
				lc[i][j] = 0;
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
		int subsizeslchunk[2] = {MSIZE/GSIZE, MSIZE/GSIZE};     /* local size */
    int starts[2]   = {0,0};                        /* where this one starts */
    MPI_Datatype typelchunk, subarrtypelchunk;
	
		//Create lc type		
		MPI_Type_create_subarray(2, sizes, subsizeslchunk, starts, MPI_ORDER_C, MPI_FLOAT, &typelchunk);
    MPI_Type_create_resized(typelchunk, 0, sizeof(float), &subarrtypelchunk);
    MPI_Type_commit(&subarrtypelchunk);

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
		int displschunk[NPROC];

    if (rank == 0) {
        for (int i=0; i<NPROC; i++) sendcounts[i] = 1;
        //int disp = 0;
        for (int i=0; i<GSIZE; i++){            
					for(int j=0; j<GSIZE; j++){
						displschunk[i*GSIZE+j] = i*MSIZE/GSIZE+j*MSIZE*MSIZE/GSIZE;
            //disp += 1;
        	}
    		}
			}

		MPI_Scatterv(aptr, sendcounts, displschunk, subarrtypelchunk, &(la[0][0]),
                 MSIZE/GSIZE*MSIZE/GSIZE, MPI_FLOAT,
                 0, MPI_COMM_WORLD);

		MPI_Scatterv(bptr, sendcounts, displschunk, subarrtypelchunk, &(lb[0][0]),
                 MSIZE/GSIZE*MSIZE/GSIZE, MPI_FLOAT,
                 0, MPI_COMM_WORLD);

	
    /* now each processor has its local array, and can process it */
		for(int c = 0;c < GSIZE;c++){
		  float tempc = 0.0;
		  for (int i=0; i<MSIZE/GSIZE; i++) {
		      for (int j=0; j<MSIZE/GSIZE; j++) {
						tempc = 0.0;
						for(int k = 0; k < MSIZE/GSIZE; k ++){
							tempc = tempc + la[i][k] * lb[k][j]; 
						}
						lc[i][j] = lc[i][j]+tempc;
		      }
		  }
			
			//rotate rows and columns
			int send_tag = 0;
			int recv_tag = 0;
			int source = rank/GSIZE;
			source = source*GSIZE + (rank+1)%GSIZE;
			int destination = rank/GSIZE;
			destination = destination*GSIZE + abs((rank - 1)%GSIZE);
/*
			source = (rank + 1)%NPROC;
			destination = (rank - 1);
			if(rank == 0){
				destination = NPROC - 1;
			}
*/
			MPI_Sendrecv_replace(&(la[0][0]), MSIZE/GSIZE*MSIZE/GSIZE, MPI_FLOAT,
			destination,send_tag, source, recv_tag, MPI_COMM_WORLD, &status);
 			

			/*MPI_Send(&(la[0][0]), MSIZE/GSIZE*MSIZE/GSIZE, MPI_FLOAT,
			destination,send_tag,MPI_COMM_WORLD);
			MPI_Recv(&(la[0][0]), MSIZE/GSIZE*MSIZE/GSIZE, MPI_FLOAT,
			source,send_tag,MPI_COMM_WORLD, &status);
	*/
			
			source = rank/GSIZE;
			source = (source+1)%GSIZE*GSIZE + (rank)%GSIZE;
			destination = rank/GSIZE;
			destination = abs((destination-1)%GSIZE)*GSIZE + (rank)%GSIZE;
			
			/*source = (rank + 1)%NPROC;
			destination = (rank - 1);
			if(rank == 0){
				destination = NPROC - 1;
			}
*/
			MPI_Sendrecv_replace(&(lb[0][0]), MSIZE/GSIZE*MSIZE/GSIZE, MPI_FLOAT,
			destination,send_tag, source, recv_tag, MPI_COMM_WORLD, &status);			
			
			/*MPI_Send(&(lb[0][0]), MSIZE/GSIZE*MSIZE/GSIZE, MPI_FLOAT,
			destination,send_tag,MPI_COMM_WORLD);
			MPI_Recv(&(lb[0][0]), MSIZE/GSIZE*MSIZE/GSIZE, MPI_FLOAT,
			source,send_tag,MPI_COMM_WORLD, &status);
			*/
		}

    /* it all goes back to process 0 */
    MPI_Gatherv(&(lc[0][0]),  MSIZE/GSIZE*MSIZE/GSIZE,  MPI_FLOAT,
                 cptr, sendcounts, displschunk, subarrtypelchunk,
                 0, MPI_COMM_WORLD);


		//get time
		t2 = MPI_Wtime(); 

    /* don't need the local data anymore */
		free2dfloat(&la);    
		free2dfloat(&lb);
    free2dfloat(&lc);

    /* or the MPI data type */
		MPI_Type_free(&subarrtypelchunk);

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
