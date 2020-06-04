#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

float gravitational_force(int x1, int y1, int z1, int mass1, int x2, int y2, int z2, int mass2); // function to calculate gravitational forces

int main(int argc, char *argv[]) {

	int my_rank , size;
	MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

	int i, j, number_of_stars;

	float global_grav_force = 0.0;
	float two_norm_sum = 0.0;
	float output; // to achieve the partial summations which calculated in each processes, in MPI_Reduce

	/* Some initializations */
	number_of_stars = atoi(argv[2]); 

	int force_info[number_of_stars][4];
	int force_of_each_process[size][number_of_stars/size];
	
	clock_t begin,end; // to determine time consuming 

	float gForce[number_of_stars/size]; // to hold forces of stars

	// FILE READÝNG

	FILE *ptr = fopen(argv[1],"r");

    if( ptr == NULL )
        printf("file could not be opened\n");

    else{
			int i = 0;
			fscanf(ptr,"%d %d %d %d",&force_info[i][0],&force_info[i][1],&force_info[i][2],&force_info[i][3]);
			while(i < number_of_stars)
			{
				i++;    
				fscanf(ptr,"%d %d %d %d",&force_info[i][0],&force_info[i][1],&force_info[i][2],&force_info[i][3]);
			}
			fclose(ptr);
        }
    int n = 0;
	// timer is started
	begin = clock(); 
	for (i = my_rank*(number_of_stars/size) ; i < (my_rank + 1)*number_of_stars/size; i++) {
		// Calculations of gravitational forces
		for (j = i+1; j < number_of_stars; j++) {

			global_grav_force += gravitational_force(force_info[i][0],force_info[i][1],force_info[i][2],force_info[i][3],
				 									 force_info[j][0],force_info[j][1],force_info[j][2],force_info[j][3]);
		}

		for (j = i - 1; j >= 0; j--) {

			global_grav_force += gravitational_force(force_info[i][0],force_info[i][1],force_info[i][2],force_info[i][3],
				 								  	 force_info[j][0],force_info[j][1],force_info[j][2],force_info[j][3]);
		}

		gForce[n] = global_grav_force;
		
		two_norm_sum += pow(gForce[n], 2);
		n++;
		global_grav_force = 0.0;

	}

	float *grav_forces = (float*)malloc(sizeof(float)*number_of_stars);

	MPI_Reduce(&two_norm_sum, &output, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
	// partial summations sended to master process
	
	MPI_Gather(gForce, number_of_stars/size, MPI_FLOAT, grav_forces, number_of_stars/size, MPI_FLOAT, 0, MPI_COMM_WORLD);
	//gravitational forces sended to master process
	
	if(my_rank == 0)
	{		
		printf("Two norm : %.10E\n", sqrtf(output));
		end = clock(); // timer stopped
		double time_taken = ((double)(end-begin))/CLOCKS_PER_SEC;
		printf("taken time :  %f\n",time_taken);

		// Writing gravitational forces to the forces.txt
		FILE *ptr2 = fopen("forces.txt","w");
		i = 0;
		while(i < number_of_stars)
		{
			fprintf(ptr,"%f\n",grav_forces[i]);
			i++;
		}
		fclose(ptr2);


	}
	
	MPI_Finalize();
	return 0;
}

float gravitational_force(int x1, int y1, int z1, int mass1, int x2, int y2,
		int z2, int mass2) {
	float force = 0.0;

	force = (0.0006674 * mass1 * mass2)
			/ (pow(
					(sqrtf(
							pow((x1 - x2), (2)) + pow((y1 - y2), (2))
									+ pow((z1 - z2), (2)))), 2));

	return force;
}

