#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

    int
    main(int argc, char **argv)
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    int **buf = NULL;
    int root = 0;
    int *global_buf = NULL;
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    global_buf = malloc(4 * sizeof(int));
    buf = malloc(2 * sizeof(int));

    // if (world_rank == 0)
    // {
    //     global_buf = malloc(world_size * sizeof(int));
    // }
    // if (world_rank == root){
    // 	buf = 777;
    // }

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;

    MPI_Get_processor_name(processor_name, &name_len);

    // MPI Broadcast Example
    // MPI_Bcast(&buf, 1, MPI_INT, root, MPI_COMM_WORLD);

    

    if (world_rank == root){
        for (int i = 0; i < 4; i++)
        {
            global_buf[i] = i;
        }
    }else{
        for (int i = 0; i < 4; i++)
        {
            buf[i] = malloc(4 * sizeof(int));
        }
    }

    MPI_Request recv_request;
    MPI_Request send_request;
    MPI_Status send_status;

    // MPI All Reduce is just like reduce but each process does the reduction
    // Locally
    // MPI_Reduce(&buf, &global_buf, 1, MPI_INT, MPI_SUM, 0,
    //        MPI_COMM_WORLD);
    // MPI_Gather(&buf, 4, MPI_INT, global_buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Print off a hello world message
    if (world_rank == root){
        MPI_Isend(global_buf, 4, MPI_DOUBLE, 1,
                  1, MPI_COMM_WORLD, &send_request);
        MPI_Wait(&send_request, &send_status);
        // printf("Processor %d has data: ", world_rank);
        // for (int i = 0; i < world_size; i++)
        //     printf("%d ", global_buf[i]);
        // printf("\n");
    }else{
        MPI_Irecv(buf[0], 4, MPI_DOUBLE, 0,
                  1, MPI_COMM_WORLD, &recv_request);
        MPI_Wait(&recv_request, &send_status);

        printf("Hello world from rank %d out of %d processors sum = ", world_rank, world_size);
        for (int i = 0; i < 4; i++)
        {
            printf("%d ", buf[0][i]);
        }
        printf("\n");
    }

    // MPI_Barrier(MPI_COMM_WORLD);
    // Finalize the MPI environment.
    MPI_Finalize();
}
