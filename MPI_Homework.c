#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"

// global variables
#define N  (16384+2)		// number of rows/cols in matrix..

// Prototypes
void Initialize(float (*x)[n], unsigned int n);
void Smooth(float (*x)[n], float (*y)[n], unsigned int n, float a, float b, float c);
void Count(float (*Mat)[n],unsigned int n, float t, int *num_below_t);
void Print_Matrix(float (*x)[n], float (*y)[n], unsigned int n);

int main(int argc, char **argv) {
	// MPI variables
	MPI_Comm comm = MPI_COMM_WORLD;
	int rank, nprocs, root = 0;

	// Initialize MPI, get rank and size
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &nprocs);

	/////////////////////////////////////////////////////////////////////////////
	// Set up a 2d grid toplolgy

	// Set up topology variables
	int dims[2];						// Two components are {rows, cols}
	int wrap_around[2] = {0,0};
	char reorder = 1;
	MPI_Comm Grid;
	signed char coords[2];

	// Set up dims array depending on number of procs
	switch nprocs {
		case 2:
			dims = {1,2};
			break;
		case 4:
			dims = {2,2};
			break;
		case 8: 
			dims = {2,4};
			break;
		case 16: 
			dims = {4,4};
			break;
		case 32:
			dims = {4,8};
			break;
		case 64:
			dims = {8,8};
			break;
		default:
			printf("Pick a different number of dimensions!");
			return 0;
	} // switch nprocs {

	// Now create the carteesian topology
	MPI_Cart_create(comm, 2, dims, wrap_around, reorder, &Grid);

	// Now get coordinates
	MPI_Cart_coords(Grid, rank, 2, coords);

	/////////////////////////////////////////////////////////////////////////////
	// Get random seeds on each proc

	// If we are the root, set a random seed using time, use this to populate
	// an array with random seeds
	int Seeds[nprocs];
	int Seed;
	if(rank == root) {
		// First, start the runtime timer. 
		runtime_timer = clock();

		// Next, set a random seed based on the current time
		srand(time(0));

		// Populate our random array.
		for(int i = 0; i < nprocs; i++) {
			Seeds[i] = rand();
		}
	} // if(rank == root) {

	// Now scatter the results to each proc
	MPI_Scatter(Seeds, 1, MPI_INT, &Seed, 1, MPI_INT, root, comm); 

	// Now set seed based on received seed
	srand(Seed);

	printf("rank: %d \tseed: %d ",rank,Seed);			// For testing purposes

	////////////////////////////////////////////////////////////////////////////////
	// Allocate and Iniitalize each task's arrays.

	// declare smoothing variables
	float (*x)[n], (*y)[n];
	float Ar_Size;
	const float a = .05, b = .1, c = .4, t = .1;
	unsigned int count_x = 0, count_y = 0;
	time_t timer, runtime_timer, t_alloc_x, t_alloc_y, t_init_x, t_smooth, t_count_x, t_count_y, t_runtime;

	// Allocate matricies and time it
	if(rank == root) { timer = clock(); }
	x = malloc(sizeof(float[N][N])); 
	if(rank == root) { t_alloc_x = clock() - timer; }

	if(rank == root) { timer = clock(); }
	y = malloc(sizeof(float[N][N]));
	if(rank == root) { t_alloc_y = clock() - timer; }

	// Calculate size of arrays in GB (only on root)
	if(rank == root) { Ar_Size = ((float)sizeof(float[N][N]))/((float)1024*1024*1024); }

	// Initialize matrix and time it.
	if(rank == root) { timer = clock(); }
	Initialize(x,n);
	if(rank == root) { t_init_x = clock()-timer; }

	////////////////////////////////////////////////////////////////////////////////
	// Now pass border elements between arrays.

	// Remember that we are really solving this problem on a big grid. In order to 
	// smooth this grid, we need the grid pieces (the matrix on each proc to agree.
	// If any of the borders of our grid touch another grid piece, we need to set the 
	// Boundary elements along that grid to the first interior row/column of the 
	// Adjacent grid piece. To do this, we must first store the border elements in 
	// arrays. We use one array for each side of the border.

	//------------------------------------------------------------------
	// Set up send, receive variables
	float L_Col_send[N], R_Col_send[N], T_Row_send[N], B_Row_send[N];
	float L_Col_recv[N], R_Col_recv[N], T_Row_recv[N], B_Row_recv[N];
	int dest_rank;
	MPI_Request send_req[4], recv_req[4];
	MPI_Status status;

	for(int i = 0; i < N; i++) {
		L_Col_send[i] = x[i][1];				// First interior column
		R_col_send[i] = x[i][(N-2)]			// Last interior column
		T_row_send[i] = x[1][i]				// First interior row
		B_Row_send[i] = x[(N-2)][i]			// LAst interior row
	} 

	//------------------------------------------------------------------
	// Send/recv borders 

	// If not leftmost column in grid: send left border, receive right border
	if(coords[1] != 0) {			// if true, we have a neighbor to the left
		// First, get coordiantes of proc to the left of us (in grid) 
		MPI_Car_shift(comm, 1, -1, &rank, &dest_rank);

		// Send our left interior column to proc to our left (in grid)
		MPI_Isend(L_Col_send, N, MPI_INT, dest_rank, 1, comm, &send_req[0]);

		// Receive right interior column of the proc to our left (in grid)
		MPI_Isend(R_Col_recv, N, MPI_INT, dest_rank, 1, comm, &recv_req[0]);
	} // if(coords[1] != 0) {

	// If not rightmost column in grid: send right border, receive left border
	if(coords[1] != dims[1]) {
		// First, get coordiantes of the proc to the right of us (in grid) 
		MPI_Car_shift(comm, 1, 1, &rank, &dest_rank);

		// Send our right interior column to the proc to our right (in grid)
		MPI_Isend(R_Col_send, N, MPI_INT, dest_rank, 1, comm, &send_req[0]);

		// Receive left interior column of the proc to our right (in grid)
		MPI_Isend(L_Col_recv, N, MPI_INT, dest_rank, 1, comm, &recv_req[0]);
	} // if(coords[1] != dims[1]) {

	// If not 1st row in grid: send top border, receive bottom border.
	if(coords[0] != 0) {
		// First, get coordiantes of proc in previous row (in grid)
		MPI_Car_shift(comm, 0, -1, &rank, &dest_rank);

		// Send our top interior row to the proc in previous row (in grid)
		MPI_Isend(T_Row_send, N, MPI_INT, dest_rank, 1, comm, &send_req[0]);

		// Receive bottom interior row of the proc in previous row (in grid)
		MPI_Isend(B_Row_recv, N, MPI_INT, dest_rank, 1, comm, &recv_req[0]);
	} // if(coords[0] != 0) {

	// If not last (bottom) row in grid: send bottom border, receive top border.
	if(coords[0] != dims[0]) {
		// First, get coordiantes of the proc in next row (in grid)
		MPI_Car_shift(comm, 0, 1, &rank, &dest_rank);

		// Send our bottom interior row to the proc in the next row (in grid)
		MPI_Isend(B_Row_send, N, MPI_INT, dest_rank, 1, comm, &send_req[0]);

		// Receive top interior rof of the proc in the next row (in grid)
		MPI_Isend(T_Row_recv, N, MPI_INT, dest_rank, 1, comm, &recv_req[0]);
	} // if(coords[0] != 0) {

	//------------------------------------------------------------------
	// Now, wait until all sending/receiving has finished
	for(i = 0; i < 4; i++) {
		// Wait until send is finished
		MPI_Wait(&send_req[i], &status);

		// Wait until recv is finished
		MPI_Wait(&recv_req[i], &status);
	} // for(i = 0; i < 4; i++) {

	//------------------------------------------------------------------
	// Now, integrate received borders into x array

	if(coords[1] != 0) {			// Case for not leftmost column
		// In this case, we received the right column of the proc to our left
		// We want to store this in our left row
		for(int i = 0; i < N; i++) {
			x[i][0] = R_Col_recv[i];
		} // for(int i = 0; i < N; i++) {
	} // if(coords[1] != 0) {

	if(coords[1] != dims[1]) {			// Case for not rightmost column
		// In this case, we received the left column of the proc to our right
		// We want to store this in our right row
		for(int i = 0; i < N; i++) {
			x[i][(N-1)] = L_Col_recv[i];
		} // for(int i = 0; i < N; i++) {
	} // if(coords[1] != dims[1]) {

	if(coords[0] != 0) {			// Case for not 1st (top) row
		// In this case, we received the bottom row of the proc in the previous row.
		// We want to store this in our bottom row
		for(int i = 0; i < N; i++) {
			x[0][i] = B_Row_recv[i];
		} // for(int i = 0; i < N; i++) {
	} // if(coords[0] != 0) {

	if(coords[0] != dims[0]) {			// Case for not last (bottom) row
		// In this case, we received the top row of the proc in the next row.
		// We want to store this in our top row
		for(int i = 0; i < N; i++) {
			x[(N-1)][i] = T_Row_recv[i];
		} // for(int i = 0; i < N; i++) {
	} // if(coords[0] != dims[0]) {

	////////////////////////////////////////////////////////////////////////////////
	// Smooth the matrix and run the count.

	// Smooth matrix and time it
	timer = clock();
	Smooth(x,y,n,a,b,c);
	t_smooth = clock() - timer;

	// Count number of elements in x,y below t.
	timer = clock();
	Count(x,n,t,&count_x);
	t_count_x = clock() - timer;

	timer = clock();
	Count(y,n,t,&count_y);
	t_count_y = clock() - timer;

	// Finish runtime timer
	t_runtime = clock() - runtime_timer;

	////////////////////////////////////////////////////////////////////////////////
	// Print results

	// Print how long these various operations took.
	printf(" --= Summary =-- \n");
	printf("Elements in a row/column \t::\t%d\n",n);
	printf("Inner elements in a row/col \t::\t%d\n",(n-2));
	printf("Total elements \t\t\t::\t%d\n",(n*n));
	printf("Memory used per node\t::\t%f (Gb)\n",Ar_Size);
	printf("threshold \t\t\t::\t%f\n",t);
	printf("Smoothing constants (a,b,c) \t::\t(%f,\t%f,\t%f)\n",a,b,c);
	printf("Number of elements below threshold (X)\t\t::\t%d\n",count_x);
	printf("Proportion of elements below threshold (X)\t::\t%f\n",(((float)count_x)/((float)n*n)));
	printf("Number of elements below threshold (Y)\t\t::\t%d\n",count_y);
	printf("Proportion of elements below threshold (Y)\t::\t%f\n",(((float)count_y)/((float)(n-2)*(n-2))));


	printf("\n--= Actions =-- \n");
	printf("Allocating x took\t::\t%f (s).\n"
		,(float)(t_alloc_x)/((float)CLOCKS_PER_SEC));
	printf("Acllocating y took\t::\t%f (s).\n"
		,(float)(t_alloc_y)/((float)CLOCKS_PER_SEC));
	printf("Initializing x took\t::\t%f (s).\n"
		,(float)(t_init_x)/(CLOCKS_PER_SEC));
	printf("Smoothing x into y took\t::\t%f (s).\n"
		,(float)(t_smooth)/(CLOCKS_PER_SEC));
	printf("Counting x took \t::\t%f (s).\n"
		,(float)(t_count_x)/(CLOCKS_PER_SEC));
	printf("Counting y took \t::\t%f (s).\n"
		,(float)(t_count_y)/(CLOCKS_PER_SEC));
	printf("Total runtime \t::\t%f (s).\n"
		,(float)(t_runtime)/(CLOCKS_PER_SEC));

	/*
	// Print out resulting matrix (for debugging purposes)
	printf("\n\n");
	Print_Matrix(x,y,n);
	*/

	return 0;
} // int main(void) {

void Initialize(float (*x)[n],unsigned int n) {
	int i,j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			x[i][j] = rand()/(float)RAND_MAX;
		} // for(j = 0; j < n; j++) {
	} // for(i = 0; i < n; i++) {
 } // void Initialize(float (*x)[n],unsigned int n) {

 void Smooth(float (*x)[n], float (*y)[n], unsigned int n, float a, float b, float c) {
 	int i, j;
 	for(i = 1; i < (n-1); i++) { 			// Note: only update interior elements (2nd to n-1st row) 
 		for(j = 1; j < (n-1); j++) {		// Only update interior elements
 			y[i][j] = 	a*(x[i-1][j-1] + x[i-1][j+1] + x[i+1][j-1] + x[i+1][j+1]) + 
 						b*(x[i-1][j] + x[i+1][j] + x[i][j-1] + x[i][j+1]) +
 					  	c*(x[i][j]);
 		} // for(j = 1; j < (n-1); j++) {
 	} // for(i = 1; i < (n-1); i++) { 
 } //  void Smooth(float (*x)[n], float (*y)[n], unsigned int n, float a, float b, float c) {

 void Count(float (*Mat)[n], unsigned int n, float t, int *num_below_t) {
 	int i,j;
 	for(i = 1; i < (n-1); i++) {
 		for(j = 1; j < (n-1); j++) {
 			if(Mat[i][j] < t) {
 				(*num_below_t)++;
 			} // if(Mat[i][j] < t) {
 		} // for(j = 1; j < (n-1); j++) {
 	} // for(i = 1; i < (n-1); i++) {
 } //  void Count(float (*Mat)[n], unsigned int n, float t, int *num_below_t) {

 void Print_Matrix(float (*x)[n], float (*y)[n],unsigned int n) {
 	int i,j;
 	for(i = 0; i < n; i++) {
 		for(j = 0; j < n; j++) {
 			printf(" %.3f ",x[i][j]);
 		} // for(j = 0; j < n; j++) {
 		printf("\n");
 	} // for(i = 0; i < n; i++) {

 	printf("\n \n");			// create some space between the two matricies

  	for(i = 1; i < (n-1); i++) {	// only print interior elements
  		 printf("\n\t");
 		for(j = 1; j < (n-1); j++) {
 			printf(" %.3f ",y[i][j]);
 		} // for(j = 1; j < (n-1); j++) {
 	} // for(i = 1; i < (n-1); i++) {
 	printf("\n");
 } // void Print_Matrix(float (*x)[n], float (*y)[n],unsigned int n) {