#include "stdio.h"
#include "stdlib.h"
#include "time.h"

// global variables
const unsigned int n = 16386;		// number of rows/cols in matrix..

// Prototypes
void Initialize(float (*x)[n], unsigned int n);
void Smooth(float (*x)[n], float (*y)[n], unsigned int n, float a, float b, float c);
void Count(float (*Mat)[n],unsigned int n, float t, int *num_below_t);
void Print_Matrix(float (*x)[n], float (*y)[n], unsigned int n);

int main(void) {
	////////////////////////////////////////////////////////////////////////////////
	// declare local variables
	float (*x)[n], (*y)[n];
	float Ar_Size;
	const float a = .05, b = .1, c = .4, t = .1;
	int count_x = 0, count_y = 0;
	time_t timer,runtime_timer, time_alloc_x, time_alloc_y, time_init_x, time_smooth, time_count_x, time_count_y, time_runtime;

	// Set random seed
	srand(time(0));

	// Start runtime timer
	runtime_timr = clock();

	////////////////////////////////////////////////////////////////////////////////
	// Run calculations (initialize, smooth, count)

	// Allocate matricies and time it
	timer = clock();
	x = malloc(sizeof(float[n][n])); 
	time_alloc_x = clock() - timer;

	timer = clock();
	y = malloc(sizeof(float[n][n]));
	time_alloc_y = clock() - timer;

	// Calculate size of arrays in GB
	Ar_Size = ((float)sizeof(float[n][n]))/((float)1024*1024*1024);

	// Initialize matrix and time it.
	timer = clock();
	Initialize(x,n);
	time_init_x = clock()-timer;

	// Smooth matrix and time it
	timer = clock();
	Smooth(x,y,n,a,b,c);
	time_smooth = clock() - timer;

	// Count number of elements in x,y below t.
	timer = clock();
	Count(x,n,t,&count_x);
	time_count_x = clock() - timer;

	timer = clock();
	Count(y,n,t,&count_y);
	time_count_y = clock() - timer;

	// Finish runtime timer
	time_runtime = clock() - runtime_timer;

	////////////////////////////////////////////////////////////////////////////////
	// Print results

	// Print how long these various operations took.
	printf(" --= Summary =-- \n");
	printf("Elements in a row/column \t::\t%d\n",n);
	printf("Inner elements in a row/col \t::\t%d\n",(n-2));
	printf("Total elements \t\t\t::\t%d\n",(n*n));
	printf("Total memory used per array\t::\t%f (Gb)\n",Ar_Size);
	printf("threshold \t\t\t::\t%f\n",t);
	printf("Smoothing constants (a,b,c) \t::\t(%f,\t%f,\t%f)\n",a,b,c);
	printf("Number of elements below threshold (X)\t\t::\t%d\n",count_x);
	printf("Proportion of elements below threshold (X)\t::\t%f\n",(((float)count_x)/((float)n*n)));
	printf("Number of elements below threshold (Y)\t\t::\t%d\n",count_y);
	printf("Proportion of elements below threshold (Y)\t::\t%f\n",(((float)count_y)/((float)(n-2)*(n-2))));


	printf("\n--= Actions =-- \n");
	printf("Allocating x took\t::\t%f (s).\n"
		,(float)(time_alloc_x)/((float)CLOCKS_PER_SEC));
	printf("Acllocating y took\t::\t%f (s).\n"
		,(float)(time_alloc_y)/((float)CLOCKS_PER_SEC));
	printf("Initializing x took\t::\t%f (s).\n"
		,(float)(time_init_x)/(CLOCKS_PER_SEC));
	printf("Smoothing x into y took\t::\t%f (s).\n"
		,(float)(time_smooth)/(CLOCKS_PER_SEC));
	printf("Counting x took \t::\t%f (s).\n"
		,(float)(time_count_x)/(CLOCKS_PER_SEC));
	printf("Counting y took \t::\t%f (s).\n"
		,(float)(time_count_y)/(CLOCKS_PER_SEC));
	printf("Total runtime \t::\t%f (s).\n"
		,(float)(time_runtime)/(CLOCKS_PER_SEC));

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