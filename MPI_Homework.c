#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
//#include "unistd.h"			// For testing

// global variables
#define N (16384+2)		// number of rows/cols in matrix..

// Prototypes
void Initialize(float *x, unsigned int n);
void Smooth(float *x, float *y, unsigned int n, float a, float b, float c);
void Count(float *Mat,unsigned int n, float t, int *num_below_t);
void Print_Matrix(float *x, float  *y, unsigned int n);
void Edge_Splice(float *x, MPI_Comm Grid, int *dims);
void Corner_Splice(float *x, MPI_Comm Grid, int *dims);

int main(int argc, char **argv) {
	// MPI variables
	MPI_Comm comm = MPI_COMM_WORLD;
	int rank, n_procs, root = 0;

	// Timing variables
	time_t timer, runtime_timer, t_alloc_x, t_alloc_y, t_init_x, t_smooth, t_count_x, t_count_y, t_runtime;

	// Loop variables
	int i;

	// Initialize MPI, get rank and size
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &n_procs);

	/////////////////////////////////////////////////////////////////////////////
	// Set up a 2d grid toplolgy

	// Set up topology variables
	int dims[2];						// Two components are {rows, cols}
	int wrap_around[2] = {0,0};
	char reorder = 1;
	MPI_Comm Grid;
	int coords[2];
	int grid_rank;

	// Set up dims array depending on number of procs
	switch(n_procs) {
		case 1: 
			dims[0] = 1;	// rows
			dims[1] = 1;	// cols
			break;
		case 2:
			dims[0] = 1;	// rows
			dims[1] = 2;	// cols
			break;
		case 4:
			dims[0] = 2;	// rows
			dims[1] = 2;	// cols
			break;
		case 8: 
			dims[0] = 2;	// rows
			dims[1] = 4;	// cols
			break;
		case 16: 
			dims[0] = 4;	// rows
			dims[1] = 4;	// cols
			break;
		case 32:
			dims[0] = 4;	// rows
			dims[1] = 8;	// cols
			break;
		case 64:
			dims[0] = 8;	// rows
			dims[1] = 8;	// cols
			break;
		default:
			printf("Pick a different number of dimensions!");
			return 0;
	} // switch(n_procs) {

	// Now create the carteesian topology
	MPI_Cart_create(comm, 2, dims, wrap_around, reorder, &Grid);

	// Now get coordinates, rank in new topology
	MPI_Cart_coords(Grid, rank, 2, coords);
	MPI_Comm_rank(Grid, &grid_rank);

	/////////////////////////////////////////////////////////////////////////////
	// Get random seeds on each proc

	// If we are the root, set a random seed using time, use this to populate
	// an array with random seeds
	int Seeds[n_procs];
	unsigned int Seed;
	if(rank == root) {
		// First, start the runtime timer. 
		runtime_timer = clock();

		// Next, set a random seed based on the current time
		srand(time(0));

		// Populate our random array.
		for(i = 0; i < n_procs; i++) {
			Seeds[i] = rand();
		}
	} // if(rank == root) {

	// Now scatter the results to each proc
	MPI_Scatter(Seeds, 1, MPI_INT, &Seed, 1, MPI_INT, root, comm); 

	// Now set seed based on received seed
	srand(Seed);

	printf("rank: %d \tseed: %d \n",rank,Seed);			// For testing purposes

	////////////////////////////////////////////////////////////////////////////////
	// Allocate and Iniitalize each task's arrays.

	// declare smoothing variables
	float *x, *y;
	float Ar_Size;
	const float a = .05, b = .1, c = .4, t = .1;
	int count_x_local = 0, count_y_local = 0, count_x, count_y;

	// Allocate matricies and time it
	if(rank == root) { timer = clock(); }
	x = (float *)malloc(sizeof(float)*N*N); 
	if(rank == root) { t_alloc_x = clock() - timer; }

	if(rank == root) { timer = clock(); }
	y = (float *)malloc(sizeof(float)*N*N);
	if(rank == root) { t_alloc_y = clock() - timer; }

	// Calculate size of arrays in GB (only on root)
	if(rank == root) { Ar_Size = ((float)sizeof(float)*N*N)/((float)1024*1024*1024); }

	// Initialize matrix and time it.
	if(rank == root) { timer = clock(); }
	Initialize(x,N);
	if(rank == root) { t_init_x = clock()-timer; }

	////////////////////////////////////////////////////////////////////////////////
	// Splice together grids
	Edge_Splice(x, Grid, dims);
	Corner_Splice(x, Grid, dims);

	////////////////////////////////////////////////////////////////////////////////
	// Smooth the matrix and get local count. 

	//sleep(1*rank+1);				// For testing
	//Print_Matrix(x,y,N);

	// Smooth matrix and time it
	if(rank == root) { timer = clock(); }
	Smooth(x,y,N,a,b,c);
	if(rank == root) { t_smooth = clock() - timer; }

	// Count number of elements in x,y below t.
	if(rank == root) { timer = clock(); }
	Count(x,N,t,&count_x_local);
	if(rank == root) { t_count_x = clock() - timer; }

	if(rank == root) { timer = clock(); }
	Count(y,N,t,&count_y_local);
	if(rank == root) { t_count_y = clock() - timer; }

	// Finish runtime timer
	if(rank == root) { t_runtime = clock() - runtime_timer; }

	////////////////////////////////////////////////////////////////////////////////
	// Reduce the local counts, get global counts

	MPI_Reduce(&count_x_local, &count_x, 1, MPI_INT, MPI_SUM, root, comm);
	MPI_Reduce(&count_y_local, &count_y, 1, MPI_INT, MPI_SUM, root, comm);

	////////////////////////////////////////////////////////////////////////////////
	// Print results

	// Print how long these various operations took.
	if(rank == root) {
		printf(" --= Summary =-- \n");
		printf("Elements in a row/column \t::\t%d\n",N);
		printf("Inner elements in a row/col \t::\t%d\n",(N-2));
		printf("Total elements \t\t\t::\t%d\n",(N*N));
		printf("Memory used per node\t::\t%f (Gb)\n",Ar_Size);
		printf("threshold \t\t\t::\t%f\n",t);
		printf("Smoothing constants (a,b,c) \t::\t(%f,\t%f,\t%f)\n",a,b,c);
		printf("Number of elements below threshold (X)\t\t::\t%d\n",count_x);
		printf("Proportion of elements below threshold (X)\t::\t%f\n",(((float)count_x)/((float)n_procs*N*N)));
		printf("Number of elements below threshold (Y)\t\t::\t%d\n",count_y);
		printf("Proportion of elements below threshold (Y)\t::\t%f\n",(((float)count_y)/((float)n_procs*(N-2)*(N-2))));


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
	} // if(rank == root) {

	/*
	// Print out resulting matrix (for debugging purposes)
	printf("\n\n");
	Print_Matrix(x,y,n);
	*/

	// Finish up MPI
	MPI_Finalize();

	return 0;
} // int main(void) {

void Initialize(float *x,unsigned int n) {
	int i,j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			x[i*n+j] = rand()/(float)RAND_MAX;
		} // for(j = 0; j < n; j++) {
	} // for(i = 0; i < n; i++) {
 } // void Initialize(float *x,unsigned int n) {


 void Smooth(float *x, float *y, unsigned int n, float a, float b, float c) {
 	int i, j;
 	for(i = 1; i < (n-1); i++) { 			// Note: only update interior elements (2nd to n-1st row) 
 		for(j = 1; j < (n-1); j++) {		// Only update interior elements
 			y[i*n+j] = 	a*(x[(i-1)*n+(j-1)] + x[(i-1)*n+(j+1)] + x[(i+1)*n+(j-1)] + x[(i+1)*n+(j+1)]) + 
 						b*(x[(i-1)*n+j] + x[(i+1)+j] + x[i*N+(j-1)] + x[i*N+(j+1)]) +
 					  	c*(x[i*N+j]);
 		} // for(j = 1; j < (n-1); j++) {
 	} // for(i = 1; i < (n-1); i++) { 
 } //  void Smooth(float *x, float *y, unsigned int n, float a, float b, float c) {

 void Count(float *Mat, unsigned int n, float t, int *num_below_t) {
 	int i,j;
 	for(i = 1; i < (n-1); i++) {
 		for(j = 1; j < (n-1); j++) {
 			if(Mat[i*n+j] < t) {
 				(*num_below_t)++;
 			} // if(Mat[i][j] < t) {
 		} // for(j = 1; j < (n-1); j++) {
 	} // for(i = 1; i < (n-1); i++) {
 } //  void Count(float *Mat, unsigned int n, float t, int *num_below_t) {

 void Print_Matrix(float *x, float *y,unsigned int n) {
 	int i,j;
 	for(i = 0; i < n; i++) {
 		for(j = 0; j < n; j++) {
 			printf(" %.3f ",x[i*n+j]);
 		} // for(j = 0; j < n; j++) {
 		printf("\n");
 	} // for(i = 0; i < n; i++) {

 	/*printf("\n \n");			// create some space between the two matricies

  	for(i = 1; i < (n-1); i++) {	// only print interior elements
  		 printf("\n\t");
 		for(j = 1; j < (n-1); j++) {
 			printf(" %.3f ",y[i*n+j]);
 		} // for(j = 1; j < (n-1); j++) {
 	} // for(i = 1; i < (n-1); i++) {
 	*/
 	printf("\n");

 } // void Print_Matrix(float *x, float *y,unsigned int n) {

 void Edge_Splice(float *x, MPI_Comm Grid, int *dims) {
 	int coords[2], grid_rank, i;

 	// get coordinates, rank in grid topology
 	MPI_Comm_rank(Grid, &grid_rank);
	MPI_Cart_coords(Grid, grid_rank, 2, coords);

	////////////////////////////////////////////////////////////////////////////////
	// Exchange border elements between procs.

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
	MPI_Status status[8];

	// First, set up send/recv requests. if these are not properly
	// initialized, bad things happen
	for(i = 0; i < 4; i++) {
		send_req[i] = MPI_REQUEST_NULL;
		recv_req[i] = MPI_REQUEST_NULL;
	} // for(i = 0; i < 4; i++) {

	// Populate send arrays with border elements of our block
	for(i = 0; i < N; i++) {
		L_Col_send[i] = x[i*N+1];				// First interior column
		R_Col_send[i] = x[i*N+(N-2)];			// Last interior column
		T_Row_send[i] = x[1*N+i];				// First interior row
		B_Row_send[i] = x[(N-2)*N + i];			// Last interior row
	} 

	//------------------------------------------------------------------
	// Send/recv borders 

	// If not leftmost column in grid: send left border, receive right border
	if(coords[1] > 0) { // if we have a neighbor to the left
		// First, get coordiantes of proc to the left of us (in grid) 
		MPI_Cart_shift(Grid, 1, -1, &grid_rank, &dest_rank);

		// Send our left interior column to proc to our left (in grid)
		MPI_Isend(L_Col_send, N, MPI_FLOAT, dest_rank, 1, Grid, &send_req[0]);

		// Receive right interior column of the proc to our left (in grid)
		MPI_Irecv(R_Col_recv, N, MPI_FLOAT, dest_rank, 1, Grid, &recv_req[0]);
	} // if(coords[1] != 0) {

	// If not rightmost column in grid: send right border, receive left border
	if(coords[1] < (dims[1]-1)) {	// if we have a neighbor to the right 
		// First, get coordiantes of the proc to the right of us (in grid) 
		MPI_Cart_shift(Grid, 1, 1, &grid_rank, &dest_rank);

		// Send our right interior column to the proc to our right (in grid)
		MPI_Isend(R_Col_send, N, MPI_FLOAT, dest_rank, 1, Grid, &send_req[1]);

		// Receive left interior column of the proc to our right (in grid)
		MPI_Irecv(L_Col_recv, N, MPI_FLOAT, dest_rank, 1, Grid, &recv_req[1]);
	} // if(coords[1] != dims[1]) {

	// If not 1st row in grid: send top border, receive bottom border.
	if(coords[0] > 0) {	// If we have a neighbor above
		// First, get coordiantes of proc in previous row (in grid)
		MPI_Cart_shift(Grid, 0, -1, &grid_rank, &dest_rank);

		// Send our top interior row to the proc in previous row (in grid)
		MPI_Isend(T_Row_send, N, MPI_FLOAT, dest_rank, 1, Grid, &send_req[2]);

		// Receive bottom interior row of the proc in previous row (in grid)
		MPI_Irecv(B_Row_recv, N, MPI_FLOAT, dest_rank, 1, Grid, &recv_req[2]);
	} // if(coords[0] != 0) {

	// If not last (bottom) row in grid: send bottom border, receive top border.
	if(coords[0] < (dims[0]-1)) { // if we have a neighbor below
		// First, get coordiantes of the proc in next row (in grid)
		MPI_Cart_shift(Grid, 0, 1, &grid_rank, &dest_rank);

		// Send our bottom interior row to the proc in the next row (in grid)
		MPI_Isend(B_Row_send, N, MPI_FLOAT, dest_rank, 1, Grid, &send_req[3]);

		// Receive top interior rof of the proc in the next row (in grid)
		MPI_Irecv(T_Row_recv, N, MPI_FLOAT, dest_rank, 1, Grid, &recv_req[3]);
	} // if(coords[0] != 0) {

	//------------------------------------------------------------------
	// Wait until all sending/receiving has finished
	for(i = 0; i < 4; i++) {
		// Wait until send is finished
		MPI_Wait(&send_req[i], &status[2*i]);

		// Wait until recv is finished
		MPI_Wait(&recv_req[i], &status[2*i+1]);
	} // for(i = 0; i < 4; i++) {

	//------------------------------------------------------------------
	// splice received borders into x array

	if(coords[1] > 0) {			// if not leftmost column
		// In this case, we received the right column of the proc to our left
		// We want to store this in our left row
		for(i = 0; i < N; i++) {
			x[i*N+0] = R_Col_recv[i];
		} // for(i = 0; i < N; i++) {
	} // if(coords[1] > 0) {

	if(coords[1] < (dims[1]-1)) { // if not rightmost column
		// In this case, we received the left column of the proc to our right
		// We want to store this in our right row
		for(i = 0; i < N; i++) {
			x[i*N+(N-1)] = L_Col_recv[i];
		} // for(i = 0; i < N; i++) {
	} // if(coords[1] < (dims[1]-1)) {

	if(coords[0] > 0) {			// if not 1st (top) row
		// In this case, we received the bottom row of the proc in the previous row.
		// We want to store this in our bottom row
		for(i = 0; i < N; i++) {
			x[0*N+i] = B_Row_recv[i];
		} // for(i = 0; i < N; i++) {
	} // if(coords[0] > 0) {

	if(coords[0] < (dims[0]-1)) {			// if not last (bottom) row
		// In this case, we received the top row of the proc in the next row.
		// We want to store this in our top row
		for(i = 0; i < N; i++) {
			x[(N-1)*N+i] = T_Row_recv[i];
		} // for(i = 0; i < N; i++) {
	} // if(coords[0] < (dims[0]-1)) {
}

void Corner_Splice(float *x, MPI_Comm Grid, int *dims) {
	int coords[2], grid_rank, i;

 	// get coordinates, rank in grid topology
 	MPI_Comm_rank(Grid, &grid_rank);
	MPI_Cart_coords(Grid, grid_rank, 2, coords);

	////////////////////////////////////////////////////////////////////////////////
	// Exchange corner elements between procs.

	float BR_Corner_Send, BL_Corner_Send, TR_Corner_Send, TL_Corner_Send,
	BR_Corner_Recv, BL_Corner_Recv, TR_Corner_Recv, TL_Corner_Recv;
	int dest_rank;
	int dest_coords[2];
	MPI_Status Corner_Status[8];
	MPI_Request send_req[4], recv_req[4];

	// First, set up send/recv requests. if these are not properly
	// initialized, bad things happen
	for(i = 0; i < 4; i++) {
		send_req[i] = MPI_REQUEST_NULL;
		recv_req[i] = MPI_REQUEST_NULL;
	} // for(i = 0; i < 4; i++) {	

	BR_Corner_Send = x[(N-2)*N+(N-2)];
	BL_Corner_Send = x[(N-2)*N+1];
	TR_Corner_Send = x[N+(N-2)];
	TL_Corner_Send = x[N+1];

	//------------------------------------------------------------------------------
	// Send/Recv corner elements

	// If there is a proc to the bottom right (if not last row or col)
	if(coords[0] < (dims[0]-1) && coords[1] < (dims[1]-1)) {
		// Get coords of proc to bottom right 
		dest_coords[0] = coords[0] + 1;
		dest_coords[1] = coords[1] + 1;
		MPI_Cart_rank(Grid, dest_coords, &dest_rank);

		// Send bottom right corner element to proc to our bottom right
		MPI_Isend(&BR_Corner_Send, 1, MPI_FLOAT, dest_rank, 1, Grid, &send_req[0]);

		// Recv top left corner element of proc to our bottom right
		MPI_Irecv(&TL_Corner_Recv, 1, MPI_FLOAT, dest_rank, 1, Grid, &recv_req[0]);
	} // if(coords[0] < (dims[0]-1) && coords[1] < (dims[1]-1)) {

	// If there is a proc to the bottom left (if not last row or 1st col)
	if(coords[0] < (dims[0]-1) && coords[1] > 0) {
		// Get coords of proc to bottom left 
		dest_coords[0] = coords[0] + 1;
		dest_coords[1] = coords[1] - 1;
		MPI_Cart_rank(Grid, dest_coords, &dest_rank);

		// Send bottom left corner element to proc to our bottom left
		MPI_Isend(&BL_Corner_Send, 1, MPI_FLOAT, dest_rank, 1, Grid, &send_req[1]);

		// Recv top right corner element of proc to our bottom left
		MPI_Irecv(&TR_Corner_Recv, 1, MPI_FLOAT, dest_rank, 1, Grid, &recv_req[1]);
	} // if(coords[0] < (dims[0]-1) && coords[1] > 0) {

	// If there is a proc to the top left (if not first row or last col)
	if(coords[0] > 0 && coords[1] < (dims[1]-1)) {
		// Get coords of proc to top right 
		dest_coords[0] = coords[0] - 1;
		dest_coords[1] = coords[1] + 1;
		MPI_Cart_rank(Grid, dest_coords, &dest_rank);

		// Send top right corner element to proc to our top right 
		MPI_Isend(&TR_Corner_Send, 1, MPI_FLOAT, dest_rank, 1, Grid, &send_req[2]);

		// Recv bottom left corner element of proc to our top right 
		MPI_Irecv(&BL_Corner_Recv, 1, MPI_FLOAT, dest_rank, 1, Grid, &recv_req[2]);
	} // if(coords[0] > 0 && coords[1] < (dims[1]-1)) {

	// If there is a proc to the top left (if not first row or col) 
	if(coords[0] > 0 && coords[1] > 0) {
		// Get coords of proc to top left 
		dest_coords[0] = coords[0] - 1;
		dest_coords[1] = coords[1] - 1;
		MPI_Cart_rank(Grid, dest_coords, &dest_rank);

		// Send top left corner element to proc to our top left
		MPI_Isend(&TL_Corner_Send, 1, MPI_FLOAT, dest_rank, 1, Grid, &send_req[3]);

		// Recv bottom right corner element of proc to our top left
		MPI_Irecv(&BR_Corner_Recv, 1, MPI_FLOAT, dest_rank, 1, Grid, &recv_req[3]);
	} // if(coords[0] > 0 && coords[1] > 0) {

	//------------------------------------------------------------------
	// Wait until all sending/receiving has finished
	for(i = 0; i < 4; i++) {
		// Wait until send is finished
		MPI_Wait(&send_req[i], &Corner_Status[2*i]);

		// Wait until recv is finished
		MPI_Wait(&recv_req[i], &Corner_Status[2*i+1]);
	} // for(i = 0; i < 4; i++) {

	//------------------------------------------------------------------
	// splice received borders into x array

	if(coords[0] < (dims[0]-1) && coords[1] < (dims[1]-1)) { // if not bottom row/col
		// In this case, we received the top left corner of the proc to our bottom right
		// We want to store this in our bottm right corner element
		x[N*N-1] = TL_Corner_Recv;
	} // if(coords[0] < (dims[0]-1) && coords[1] < (dims[1]-1)) {

	if(coords[0] < (dims[0]-1) && coords[1] > 0)  { // if not bottom row/first col
		// In this case, we received the top right corner element from the proc
		// to our bottom left. We want to store this in our bottom left corner element
		x[(N-1)*N] = TR_Corner_Recv;
	} // if(coords[0] < (dims[0]-1) && coords[1] > 0) { 

	if(coords[0] > 0 && coords[1] < (dims[1]-1)) {	// if not 1st row of last col 
		// In this casse, we received the bottom left corner element from the proc to our 
		// upper right. We want to store this element in our upper right corener element.
		x[N-1] = BL_Corner_Recv;
	} // if(coords[0] > 0 && coords[1] < (dims[1]-1)) {	

	if(coords[0] > 0 && coords[1] > 0) {			// if not 1st row/col
		// In this casse, we received the bottom right corner element from the proc to our 
		// upper left. We want to store this element in our upper left corener element.
		x[0] = BR_Corner_Recv;
	} // if(coords[0] > 0 && coords[1] > 0) {
}