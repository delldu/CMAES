/************************************************************************************
***
***	Copyright 2021 Dell(18588220928g@163.com), All Rights Reserved.
***
***	File Author: Dell, 2021-03-06 15:46:55
***
************************************************************************************/

// Gradient Free and Black Box Search Demo

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

#include "cmaes.h"
#include "boundary.h"

#define Z_SPACE_DIM 512

double fitfun(double const *x, unsigned long N)
{								/* function "cigtab" */
	unsigned i;
	double sum = 0;
	for (i = 0; i < N; i++) {
		if (i % 2 == 0) {
			sum += (x[i] - 0.8) * (x[i] - 0.8);		// +0.8
		} else {
			sum += (x[i] + 0.8) * (x[i] + 0.8);		// -0.8
		}
	}

	return sum;
}

double *search(int maxtimes)
{
	cmaes_t evo;				/* an CMA-ES type struct or "object" */
	cmaes_boundary_t boundaries;
	double *cost_values, *x_in_bounds, *const *pop;
	double lowerBounds[] = { -1.0, -100.0 };	/* last number is recycled for all remaining coordinates */
	double upperBounds[] = { +1.0, 100.0 };
	int nb_bounds = 1;			/* numbers used from lower and upperBounds */

	double xstart[Z_SPACE_DIM];
	double stddev[Z_SPACE_DIM];
	int lambda;

	unsigned long dimension;
	int i;

	/* initialize boundaries ... */
	cmaes_boundary_init(&boundaries, lowerBounds, upperBounds, nb_bounds);

	/* Initialize everything into the struct evo, 0 means default */
	for (i = 0; i < Z_SPACE_DIM; i++) {
		xstart[i] = 0.0;
		stddev[i] = 1.0;
	}
	lambda = 4 + 3*9;	// 4 + 3*log(Z_SPACE_DIM)
	cost_values = cmaes_init(&evo, Z_SPACE_DIM /*dimmesion*/, xstart, stddev, 0, lambda, "none");


	dimension = (unsigned long) cmaes_Get(&evo, "dimension");


	printf("%s\n", cmaes_SayHello(&evo));
	x_in_bounds = cmaes_NewDouble(dimension);	/* calloc another vector */
	// cmaes_ReadSignals(&evo, "cmaes_signals.par");	/* write header and initial values */

	/* Iterate until stop criterion holds */
	// stopMaxFunEvals  4.5e6  # stop after given number of function evaluations
	// stopMaxIter 3e3         # stop after given number of iterations (generations)
	// stopFitness 1e-6     # stop if function value is smaller than stopFitness
	// stopTolFun 1e-12     # stop if function value differences are small

	evo.sp.stopMaxFunEvals = 32*maxtimes;
	evo.sp.stopMaxIter = maxtimes;
	// evo.sp.stStopFitness.val = 1e-3;
	evo.sp.stopTolFun = 1e-6;

	while (!cmaes_TestForTermination(&evo)) {
		/* generate lambda new search points, sample population */
		pop = cmaes_SamplePopulation(&evo);	/* do not change content of pop */

		/* transform into bounds and evaluate the new search points */
		for (i = 0; i < cmaes_Get(&evo, "lambda"); ++i) {
			cmaes_boundary_transformation(&boundaries, pop[i], x_in_bounds, dimension);
			cost_values[i] = fitfun(x_in_bounds, dimension);	/* evaluate */
		}

		/* update the search distribution used for cmaes_SampleDistribution() */
		cmaes_UpdateDistribution(&evo, cost_values);	/* assumes that pop[i] has not been modified */

		/* read instructions for printing output or changing termination conditions */
		// cmaes_ReadSignals(&evo, "cmaes_signals.par");
		fflush(stdout);			/* useful in MinGW */
	}
	printf("Stop:\n%s\n", cmaes_TestForTermination(&evo));	/* print termination reason */
	cmaes_WriteToFile(&evo, "all", "allcmaes.dat");	/* write final results */

	/* get best estimator for the optimum, xmean */
	cmaes_boundary_transformation(&boundaries, (double const *) cmaes_GetPtr(&evo, "xmean"),	/* "xbestever" might be used as well */
								  x_in_bounds, dimension);

	cmaes_exit(&evo);			/* release memory */
	cmaes_boundary_exit(&boundaries);	/* release memory */

	return x_in_bounds;
}

void help(char *cmd)
{
	printf("Usage: %s [option]\n", cmd);
	printf("    -h, --help                   Display this help.\n");

	exit(1);
}

int main(int argc, char **argv)
{
	int optc;
	int option_index = 0;

	struct option long_opts[] = {
		{ "help", 0, 0, 'h'},
		{ 0, 0, 0, 0}

	};

	if (argc <= 1)
		help(argv[0]);
	
	while ((optc = getopt_long(argc, argv, "h", long_opts, &option_index)) != EOF) {
		switch (optc) {
		case 'h':	// help
		default:
			help(argv[0]);
			break;
	    	}
	}

	int i;
	double *result = search(3000);

	for (i = 0; i < 512; i++)
		printf("%.2lf ", result[i]);
	printf("\n");

	free(result);

	return 0;
}
