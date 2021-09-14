#include"exchange.h"

//1D exchange along the rows
void exchange(int g[m][n], int s, int e, int nbrtop, int nbrbot, MPI_Comm comm){
        MPI_Request reqs[4];

        MPI_Irecv(&g[modulo(e+1,m)][0], n, MPI_INT, nbrbot, 0, comm, &reqs[0]);
        MPI_Irecv(&g[modulo(s-1,m)][0], n, MPI_INT, nbrtop, 1, comm, &reqs[1]);
        MPI_Isend(&g[modulo(s,m)][0], n, MPI_INT, nbrtop, 0, comm, &reqs[2]);
        MPI_Isend(&g[modulo(e,m)][0], n, MPI_INT, nbrbot, 1, comm, &reqs[3]);

//	Alternative blocking exchange function - used for testing
//	MPI_Sendrecv(&g[modulo(s,m)][0], n, MPI_INT, nbrtop, 0, &g[modulo(e+1,m)][0], n, MPI_INT, nbrbot, 0, comm, MPI_STATUS_IGNORE);
//	MPI_Sendrecv(&g[modulo(e,m)][0], n, MPI_INT, nbrbot, 1, &g[modulo(s-1,m)][0], n, MPI_INT, nbrtop, 1, comm, MPI_STATUS_IGNORE);

        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
}

//2D exchange along the rows and cols
void exchange2d(int g[m][n], int s, int e, int s2, int e2, int nbrtop, int nbrbot, int nbrleft, int nbrright, MPI_Comm comm){
        MPI_Request reqs[8];

        //Send columns
        MPI_Datatype col;
        MPI_Type_vector(e-s, 1, n, MPI_INT, &col);//e-s
        MPI_Type_commit(&col);
        MPI_Irecv(&g[s][modulo(e2+1,n)], 1, col, nbrright, 0, comm, &reqs[0]);
        MPI_Irecv(&g[s][modulo(s2-1,n)], 1, col, nbrleft, 1, comm, &reqs[1]);
        MPI_Isend(&g[s][s2], 1, col, nbrleft, 0, comm, &reqs[2]);
        MPI_Isend(&g[s][e2], 1, col, nbrright, 1, comm, &reqs[3]);

        //Send rows
        MPI_Irecv(&g[modulo(e+1,m)][s2], e2-s2+1, MPI_INT, nbrbot, 2, comm, &reqs[4]);
        MPI_Irecv(&g[modulo(s-1,m)][s2], e2-s2+1, MPI_INT, nbrtop, 3, comm, &reqs[5]);
        MPI_Isend(&g[modulo(s,m)][s2], e2-s2+1, MPI_INT, nbrtop, 2, comm, &reqs[6]);
        MPI_Isend(&g[modulo(e,m)][s2], e2-s2+1, MPI_INT, nbrbot, 3, comm, &reqs[7]);

        MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
        MPI_Type_free(&col);
}

//2D exchange along the rows and cols as well as corner integers
void exchange2d_tri(int g[m][n], int s, int e, int s2, int e2, int nbrtop, int nbrbot, int nbrleft, int nbrright, int nbrd1, int nbrd2, int nbrd3, int nbrd4, MPI_Comm comm){
        MPI_Request reqs[16];

        //Send columns
        MPI_Datatype col;
        MPI_Type_vector(e-s, 1, n, MPI_INT, &col);//e-s
        MPI_Type_commit(&col);
        MPI_Irecv(&g[s][modulo(e2+1,n)], 1, col, nbrright, 0, comm, &reqs[0]);
        MPI_Irecv(&g[s][modulo(s2-1,n)], 1, col, nbrleft, 1, comm, &reqs[1]);
        MPI_Isend(&g[s][s2], 1, col, nbrleft, 0, comm, &reqs[2]);
        MPI_Isend(&g[s][e2], 1, col, nbrright, 1, comm, &reqs[3]);

        //Send rows
        MPI_Irecv(&g[modulo(e+1,m)][s2], e2-s2+1, MPI_INT, nbrbot, 2, comm, &reqs[4]);
        MPI_Irecv(&g[modulo(s-1,m)][s2], e2-s2+1, MPI_INT, nbrtop, 3, comm, &reqs[5]);
        MPI_Isend(&g[modulo(s,m)][s2], e2-s2+1, MPI_INT, nbrtop, 2, comm, &reqs[6]);
        MPI_Isend(&g[modulo(e,m)][s2], e2-s2+1, MPI_INT, nbrbot, 3, comm, &reqs[7]);

	//Send corncer integers to diagonal neighbours

	//Receive top corners
	MPI_Irecv(&g[modulo(s-1,m)][modulo(s2-1,n)], 1, MPI_INT, nbrd1, 4, comm, &reqs[8]);
	MPI_Irecv(&g[modulo(s-1,m)][modulo(e2+1,n)], 1, MPI_INT, nbrd2, 5, comm, &reqs[9]);
	MPI_Isend(&g[e][e2], 1, MPI_INT, nbrd4, 4, comm, &reqs[10]);
	MPI_Isend(&g[e][s2], 1, MPI_INT, nbrd3, 5, comm, &reqs[11]);

	//Receive bottom corners
        MPI_Irecv(&g[modulo(e+1,m)][modulo(s2-1,n)], 1, MPI_INT, nbrd3, 6, comm, &reqs[12]);
        MPI_Irecv(&g[modulo(e+1,m)][modulo(e2+1,n)], 1, MPI_INT, nbrd4, 7, comm, &reqs[13]);
        MPI_Isend(&g[s][e2], 1, MPI_INT, nbrd2, 6, comm, &reqs[14]);
        MPI_Isend(&g[s][s2], 1, MPI_INT, nbrd1, 7, comm, &reqs[15]);


        MPI_Waitall(16, reqs, MPI_STATUSES_IGNORE);
        MPI_Type_free(&col);
}


//Another 1D exchange function which exchanges odd and even processes separately (Obsolete)
void exchange2(int g[m][n], int s, int e, int nbrtop, int nbrbot, int rank, int turn, MPI_Comm comm){
        if(turn%2==0){
                if(rank%2==0){
                        MPI_Send(&g[modulo(e,m)][0], n, MPI_INT, nbrbot, 0, comm);
                        MPI_Send(&g[modulo(s,m)][0], n, MPI_INT, nbrtop, 0, comm);
                }
                if(rank%2==1){
                        MPI_Recv(&g[modulo(e+1,m)][0], n, MPI_INT, nbrbot, 0, comm, MPI_STATUS_IGNORE);
                        MPI_Recv(&g[modulo(s-1,m)][0], n, MPI_INT, nbrtop, 0, comm, MPI_STATUS_IGNORE);
                }
        }
        if(turn%2==1){
                if(rank%2==1){
                        MPI_Send(&g[modulo(s,m)][0], n, MPI_INT, nbrtop, 0, comm);
                        MPI_Send(&g[modulo(e,m)][0], n, MPI_INT, nbrbot, 0, comm);
                }
                if(rank%2==0){
                        MPI_Recv(&g[modulo(s-1,m)][0], n, MPI_INT, nbrtop, 0, comm, MPI_STATUS_IGNORE);
                        MPI_Recv(&g[modulo(e+1,m)][0], n, MPI_INT, nbrbot, 0, comm, MPI_STATUS_IGNORE);
                }
        }
}

