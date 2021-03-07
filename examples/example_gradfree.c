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

double *search(int epochs)
{
	int i, dimension;
	cmaes_t evo;
	double *cost_values, *xfinal, *const *pop;

	// cost_values = cmaes_init(&evo, Z_SPACE_DIM dimmesion, xstart, stddev, 0, / * lambda */, "none");
	cmaes_init_para(&evo, Z_SPACE_DIM /*dimmesion*/, NULL, NULL, 0, 0 /* lambda */, "none");

	evo.sp.stopMaxIter = epochs;	// stop after given number of iterations (generations)
	evo.sp.stopFitness.flg = 1;
	evo.sp.stopFitness.val = 1e-3;	// stop if function value is smaller than stopFitness
	evo.sp.stopTolFun = 1e-4;		// stop if function value differences are small
	cost_values =  cmaes_init_final(&evo);


	dimension = (unsigned int) cmaes_Get(&evo, "dimension");
	printf("%s\n", cmaes_SayHello(&evo));

	while (!cmaes_TestForTermination(&evo)) {
		/* generate lambda new search points, sample population */
		pop = cmaes_SamplePopulation(&evo);	/* do not change content of pop */

		/* transform into bounds and evaluate the new search points */
		for (i = 0; i < cmaes_Get(&evo, "lambda"); ++i) {
			cost_values[i] = fitfun(pop[i], dimension);	/* evaluate */
		}

		cmaes_UpdateDistribution(&evo, cost_values);	/* assumes that pop[i] has not been modified */

		if ((int)evo.gen % 10 == 0) {
			printf("Progress %6.2f %% ...\n", (float)(100.0 * evo.gen/epochs));
			fflush(stdout);
		}
	}
	printf("Stop Condition:\n%s\n", cmaes_TestForTermination(&evo));	/* print termination reason */

	cmaes_WriteToFile(&evo, "all", "/tmp/cmaes_results.txt");	/* write final results */

	/* get best estimator for the optimum, xmean */
	xfinal = cmaes_GetNew(&evo, "xmean");

	cmaes_exit(&evo);			/* release memory */

	return xfinal;
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
	double *result = search(1024);

	for (i = 0; i < 512; i++)
		printf("%.2lf ", result[i]);
	printf("\n");

	free(result);

	return 0;
}
