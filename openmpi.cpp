#include <mpi.h>
#include <omp.h>
#include<iostream>
#include <stdio.h>
#include<bits/stdc++.h>
#include <sys/time.h>   
#include <string>
using namespace std; 
// Return current wallclock time, for performance measurement
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}
int **new2d (int width, int height)
{
	int **dp = new int *[width];
	size_t size = width;
	size *= height;
	int *dp0 = new int [size];
	if (!dp || !dp0)
	{
	    std::cerr << "getMinimumPenalty: new failed" << std::endl;
	    exit(1);
	}
	dp[0] = dp0;
	for (int i = 1; i < width; i++)
	    dp[i] = dp[i-1] + height;

	return dp;
}
int max(int a, int b);
int min(int a, int b); 
/* Returns length of LCS for X[0..m-1], Y[0..n-1] */
 
/* Utility function to get max of 2 integers */
int max(int a, int b)
{
    return (a > b)? a : b;
}
int min(int a, int b)
{
    return (a >= b)? b : a;
}
int main(int argc, char** argv) {
    // std::string X = "AGGTAB";
    // std::string Y = "GXTXAYB";

    uint64_t start = GetTimeStamp();
    std::string X, Y;
    std::cin >> X;
    std::cin >> Y;
    int m = X.length(); 
    int n = Y.length();
    int provided;
    MPI_Init_thread(NULL,NULL,MPI_THREAD_FUNNELED,&provided);
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int root = 0;

    // printf("The char at process ranked %d is : %c\n",rank, X[rank]);
    MPI_Bcast(&m,1,MPI_INT,root,MPI_COMM_WORLD);
    MPI_Bcast(&n,1,MPI_INT,root,MPI_COMM_WORLD);

    int matrix_length = m+n+1;
    int matrix_width = m+1;
    int subsection_length = 3;
    //int **L = new2d (subsection_length, matrix_width); 
    int L[subsection_length][matrix_width];
    int i, j;

    if (rank==root)
    { 
        for (i = 1; i < world_size; i++)
        {            
            MPI_Send(X.c_str(),m+1,MPI_CHAR,i,0,MPI_COMM_WORLD);
            MPI_Send(Y.c_str(),n+1,MPI_CHAR,i,0,MPI_COMM_WORLD);      
        }
        for (i = 0; i < subsection_length; i++)
        {           
            for (j = 0; j < matrix_width; j++)
            {
                L[i][j] = 0;
            }
        }
    }
    
    if (rank != root)
    {
        char *x_buf = (char*) malloc((m+1)*sizeof(char));
        char *y_buf = (char*) malloc((n+1)*sizeof(char));
        MPI_Recv(x_buf,m+1,MPI_CHAR,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(y_buf,n+1,MPI_CHAR,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        X = x_buf;
        Y = y_buf;
    }

    MPI_Bcast(L,subsection_length*matrix_width,MPI_INT,root,MPI_COMM_WORLD);
    
    
//    for (int i = 0; i < subsection_length; i++)
//                     {
//                         for (int j = 0; j < matrix_width; j++)
//                         {
//                             cout << L[i][j] << " ";
//                         }
                
//                     // Newline for new row
//                     cout << endl;
//                     }
    MPI_Barrier(MPI_COMM_WORLD);
    int total_load;
    int quotient = matrix_length/world_size;
    int remain = matrix_length%world_size;
    int local_load =rank<remain? quotient + 1:quotient;
    int start_point;
    int end_point;
    int local_start_point;
    int local_end_point;
    int k = subsection_length-1;
    int discard;
    int buffer_size = 2*(sizeof(int) + MPI_BSEND_OVERHEAD);
    int *buffer =  (int *) malloc(buffer_size);
    MPI_Status status;
    int position;
    int source;
    for (i=2; i<=m+n; i++)
        {
            start_point= j=max(1,i-n);
            end_point= min(i-1,m);
            total_load = end_point - start_point + 1;
            quotient =total_load/world_size;
            remain = total_load%world_size;
            local_load =rank<remain? quotient + 1:quotient;
            if (rank <remain)
            {
                local_start_point = start_point + (quotient+1)*rank;
                local_end_point = local_start_point +local_load;
            }
            else  
            {
                local_start_point = start_point + quotient*rank + remain;
                local_end_point = local_start_point +local_load;
            }
            #pragma omp parallel for schedule(static) shared(X,Y,L) private(j) num_threads(8) 
            for (j=local_start_point; j<local_end_point; j++)
            {
                if (Y[i-j-1] == X[j-1])
                    L[k][j] = L[k-2][j-1] + 1;
            
                else
                    L[k][j] = max(L[k-1][j], L[k-1][j-1]);
            }

            //printf("i is %d\n", i);
            //Send results to all other processes
            MPI_Buffer_attach(buffer, buffer_size);
            int tag = local_load == 0? 1:0;
            if (rank > 0 && rank < world_size -1)
            {
                MPI_Bsend(&L[k][local_end_point -1],1,MPI_INT,rank+1,tag,MPI_COMM_WORLD); 
                MPI_Bsend(&L[k][local_start_point],1,MPI_INT,rank-1,tag,MPI_COMM_WORLD);   
            } else if(rank == 0 && world_size >1)
            {
                //printf("source: %d, dest: %d\n",rank, rank+1);
                MPI_Bsend(&L[k][local_end_point-1],1,MPI_INT,rank+1,tag,MPI_COMM_WORLD);
            } else if(rank ==  world_size -1 && world_size >1)
            {
                MPI_Bsend(&L[k][local_start_point],1,MPI_INT,rank-1,tag,MPI_COMM_WORLD); 
            }
            
            if (rank > 0 && rank < world_size -1)
            {   
                for (j = 0; j < 2; j++)
                {
                    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    source = status.MPI_SOURCE;
                    if (status.MPI_TAG == 1)
                    {
                        MPI_Recv(&discard,1,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    } else if (rank >source)
                    {
                        position = local_start_point-1;
                        MPI_Recv(&L[k][position],1,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);                  
                    } else if (rank < source)
                    {
                        position = local_end_point;
                        MPI_Recv(&L[k][position],1,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);   
                    }
                    
                    
                }        
            } else if((rank == 0 || rank==world_size-1) && world_size >1)
            {            
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                source = status.MPI_SOURCE;
                if (status.MPI_TAG == 1)
                {
                    MPI_Recv(&discard,1,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                } else if (rank >source)
                {
                    position = local_start_point-1;
                    MPI_Recv(&L[k][position],1,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);                  
                } else if (rank < source)
                {
                    position = local_end_point;
                    MPI_Recv(&L[k][position],1,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);   
                }  
            }
            // if (rank == 0)
            // {
            //     for (int i = k; i < subsection_length; i++)
            //     {
            //         for (int j = 0; j < matrix_width; j++)
            //         {
            //             cout << L[i][j] << " ";
            //         }
            
            //     // Newline for new row
            //     cout << endl;
            //     }
            // }
            
            
                
              
            //MPI_Bcast(L[k],matrix_width,MPI_INT,root,MPI_COMM_WORLD);   
            for(j=0; j<matrix_width; ++j)
                L[0][j] = L[1][j];
            for(j=0; j<matrix_width; ++j)
                L[1][j] = L[2][j];
            for(j=0; j<matrix_width; ++j)
                L[2][j] = 0;

            MPI_Barrier(MPI_COMM_WORLD);

            // if (i==4042)
            // {
            //     printf("barrier reached at %d\n", rank);
            // }
        }
    if (rank == root)
    {
        /* Following steps build L[m+1][n+1] in bottom up fashion. Note
        that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1] */
        

        int result = L[k-1][m];

        printf("Length of LCS is %d\n", result);

        // print the time taken to do the computation
	    printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
    }
    
    
    MPI_Finalize();

    return 0;
  
}