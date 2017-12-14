#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <conio.h>

#define PI 3.14159265358979

#define MAX 10
#define MAXNODES 10
#define MAXINPUTS 10
#define MAXOUTPUTS 10
#define MAXPROBES 10
#define MAXCONSTS 10
#define MAXELEMENTS 200
#define MAXSTEPBUFFER 30
#define TAMBUFFER 15
#define TYPE_INPUT 1
#define TYPE_OUTPUT 2
#define TYPE_CONST 3
#define TYPE_PROBE 4
#define TYPE_SUM 5
#define TYPE_GAIN 6
#define TYPE_INTEGRATOR 7
#define TYPE_DELAY 8

FILE *arq;

FILE *resfile;

void main (void) {

	int nodes[MAXNODES];
	long double inputs[MAXNODES][2*MAXSTEPBUFFER];
	long double inputparms[MAXNODES][10];
	long double outputs[MAXNODES][2*MAXSTEPBUFFER];
	long double probes[MAXNODES];
	long double consts[MAXNODES];
	long double nodevalue[MAXNODES][2*MAXSTEPBUFFER];
	int typeelemmtx[MAXNODES][MAXNODES];
	long double valuemtx[MAXNODES][MAXNODES];
	int nlinks[MAXNODES];
	int nsums[MAXNODES];
	int tobeprocessed[MAXNODES][MAXNODES];
	int nodeinput[MAXNODES];
	int nodeoutput[MAXNODES];
	int nodeprobe[MAXNODES];
	int nodeconst[MAXNODES];

	long double t, tmax, tstep;
	int stepbuffer, actualnode, resultnode;

	int i,j,k;

	int processflag;

	char parms1[5],parms2[5],parms3[10],parms4[16];
	//char parms5[16];
	int parmsi1, parmsi2;
	int nnodes, nsteps;
	long double parmsf4;
	//long double parmsf5;

	int ninputs, noutputs, nprobes, nconsts;

	// Limpeza das matrizes e vetores
	for (i=0; i<MAXNODES; i++) {
		nlinks[i] = 0;
		nsums[i] = -1; // -1 serve como flag para tobeprocessed
		nodeinput[i] = 0;
		nodeoutput[i] = 0;
		nodeprobe[i] = 0;
		nodes[i] = 0;
		for (j=0; j<MAXNODES; j++) {
			typeelemmtx[i][j] = -1;
			valuemtx[i][j] = 0;
			tobeprocessed[i][j] = -1;
		}
		for (k=0; k<2*MAXSTEPBUFFER; k++) {
			inputs[i][k] = 0;
			outputs[i][k] = 0;
			nodevalue[i][k] = 0;
		}
	}

	arq = fopen("sist04.bds","rt");
	resfile = fopen("result04.m","wt");
	fprintf(resfile, "saida = [");

    // Open for read (will fail if file "data" does not exist)
	if(arq == NULL ) {
      printf( "The file was not opened\n" );
	  //fprintf(resfile, "The file was not opened\n" );
	}
	else {
      printf( "The file was opened\n" );
	  //fprintf(resfile, "The file was opened\n" );
	}
	
	fscanf(arq,"%s",parms1);
	tmax = atof(parms1);
	fscanf(arq,"%s",parms2);
	tstep = atof(parms2);

	ninputs = 0;
	noutputs = 0;
	nprobes = 0;
	nconsts = 0;
	nnodes = 0;

	while (!feof(arq)) {

		// leitura do arquivo
		fscanf(arq,"%s",parms1);
		parmsi1 = atoi(parms1);

		fscanf(arq,"%s",parms2);
		parmsi2 = atoi(parms2);

		fscanf(arq,"%s",parms3);

		fscanf(arq,"%s",parms4);
		parmsf4 = atof(parms4);

		// contagem dos nohs existentes
		if (nodes[parmsi1] == 0) {
			nodes[parmsi1]++;
			nnodes++;
		}
		if (nodes[parmsi2] == 0) {
			nodes[parmsi2]++;
			nnodes++;
		}

		// se a função é de SOMATÓRIA (SUM)
		if (strcmp(parms3,"SUM")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_SUM;
			valuemtx[parmsi1][parmsi2] = 0;
			tobeprocessed[parmsi1][parmsi2] = 1;
			tobeprocessed[parmsi2][parmsi2]--; 
			nlinks[parmsi1]++;
			nsums[parmsi2]--;
		}

		// se é função de GANHO (GAIN)
		if (strcmp(parms3,"GAIN")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_GAIN;
			valuemtx[parmsi1][parmsi2] = parmsf4;
			tobeprocessed[parmsi1][parmsi2] = 1;
			nlinks[parmsi1]++;
		}

		// se é variável de OUTPUT
		if (strcmp(parms3,"OUTPUT")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_OUTPUT;
			nodeoutput[noutputs] = parmsi1;
			noutputs++;
		}

		// se é variável de INPUT
		if (strcmp(parms3,"INPUT")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_INPUT;
			valuemtx[parmsi1][parmsi2] = parmsf4;
			inputs[parmsi1][0] = parmsf4;
			inputparms[parmsi1][0] = parmsf4;
			nodeinput[ninputs] = parmsi1;
			ninputs++;
		}

		// se é variável de PROBE
		if (strcmp(parms3,"PROBE")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_PROBE;
			nodeprobe[nprobes] = parmsi1;
			nprobes++;
		}

		// se é variável de CONST
		if (strcmp(parms3,"CONST")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_CONST;
			nodeconst[nconsts] = parmsi1;
			consts[parmsi1] = parmsf4;
			nconsts++;
		}
		// se é variável de INTEGRATOR
		if (strcmp(parms3,"INTEGRATOR")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_INTEGRATOR;
			valuemtx[parmsi1][parmsi2] = parmsf4;
			tobeprocessed[parmsi1][parmsi2] = 1;
			nlinks[parmsi1]++;
		}
		// se é variável de DELAY
		if (strcmp(parms3,"DELAY")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_DELAY;
			valuemtx[parmsi1][parmsi2] = parmsf4;
			tobeprocessed[parmsi1][parmsi2] = 1;
			nlinks[parmsi1]++;
		}
	}
	
	// Aqui começa a análise dos blocos
	// Tentativa 3: análise por nós / matriz de relacionamentos
	
	t = 0;
	nsteps = 0;
	stepbuffer = 1; // começa de 1, pois existe a condição inicial 0
	while (t<tmax) {
		processflag = 0;
		// Transferir INPUTS e CONSTS se existirem
		for (i=0; i<ninputs; i++) {
			//inputs[nodeinput[i]][stepbuffer] = (i+1)*sin(2*PI*60*t);
			inputs[nodeinput[i]][stepbuffer] = inputparms[nodeinput[i]][0]*sin(2*PI*60*t);
		}

		for (actualnode = 0; actualnode<nnodes; actualnode++) {
			// Atualizar tobeprocessed
			if (typeelemmtx[actualnode][actualnode] == TYPE_INPUT) {
				valuemtx[actualnode][actualnode] = inputs[actualnode][stepbuffer];
				tobeprocessed[actualnode][actualnode] = nlinks[actualnode];
				processflag = 1;
			}
			if (typeelemmtx[actualnode][actualnode] == TYPE_CONST) {
				valuemtx[actualnode][actualnode] = consts[actualnode];
				tobeprocessed[actualnode][actualnode] = nlinks[actualnode];
				processflag = 1;
			}
			for (resultnode = 0; resultnode<nnodes; resultnode++) {
				if (typeelemmtx[actualnode][resultnode] == TYPE_SUM) {
					tobeprocessed[actualnode][resultnode] = 1;
					tobeprocessed[resultnode][resultnode] = nsums[resultnode];
					processflag = 1;
				}
				if (typeelemmtx[actualnode][resultnode] == TYPE_DELAY) {
					tobeprocessed[actualnode][resultnode] = 1;
					//tobeprocessed[actualnode][actualnode] = nlinks[actualnode];
					processflag = 1;
				}
			}
		}
		while (processflag) {
			processflag = 0;
			for (actualnode = 0; actualnode<nnodes; actualnode++) {
				for (resultnode = 0; resultnode<nnodes; resultnode++) {
					// GAIN
					if (typeelemmtx[actualnode][resultnode] == TYPE_GAIN) {
						if (tobeprocessed[actualnode][resultnode] !=0 && tobeprocessed[actualnode][actualnode] > 0) {
							tobeprocessed[resultnode][resultnode] = nlinks[resultnode];
							tobeprocessed[actualnode][actualnode]--;
							tobeprocessed[actualnode][resultnode] = 0;
							valuemtx[resultnode][resultnode] = valuemtx[actualnode][resultnode]*valuemtx[actualnode][actualnode];
							processflag = 1;
						}
					} // GAIN
					
					// SUM
					if (typeelemmtx[actualnode][resultnode] == TYPE_SUM) {
						if (tobeprocessed[actualnode][resultnode] !=0 && tobeprocessed[actualnode][actualnode] > 0) {
							tobeprocessed[actualnode][actualnode]--;
							tobeprocessed[actualnode][resultnode] = 0;
							tobeprocessed[resultnode][resultnode]++;
							if (tobeprocessed[resultnode][resultnode]==-1) {
								tobeprocessed[resultnode][resultnode] = nlinks[resultnode];
							}
							valuemtx[resultnode][resultnode] = valuemtx[resultnode][resultnode]+valuemtx[actualnode][actualnode];
							processflag = 1;
						}
					} // SUM

					// DELAY
					if (typeelemmtx[actualnode][resultnode] == TYPE_DELAY) {
						if (tobeprocessed[actualnode][resultnode] > 0) {
							//tobeprocessed[actualnode][actualnode]--;
							tobeprocessed[actualnode][resultnode] = 0;
							tobeprocessed[resultnode][resultnode] = nlinks[resultnode];
							valuemtx[resultnode][resultnode] = nodevalue[actualnode][TAMBUFFER+stepbuffer-1];
							processflag = 1;
						}
					} // DELAY
					
					// INTEGRATOR
					if (typeelemmtx[actualnode][resultnode] == TYPE_INTEGRATOR) {
						if (tobeprocessed[actualnode][resultnode] !=0 && tobeprocessed[actualnode][actualnode] > 0) {
							tobeprocessed[actualnode][actualnode]--;
							tobeprocessed[actualnode][resultnode] = 0;
							tobeprocessed[resultnode][resultnode]++;
							// valuemtx[resultnode][resultnode] = ;
							processflag = 1;
						}
					} // INTEGRATOR

				}
			}
		}

		// IMPLEMENTAR BUFFER CIRCULAR!!!

		for (i=0; i<nnodes; i++) {
			nodevalue[i][stepbuffer] = valuemtx[i][i];
			nodevalue[i][TAMBUFFER+stepbuffer] = valuemtx[i][i];
		}

		fprintf(resfile, "%lf ", t);
		for (i=0; i<nnodes; i++) {
			//for (j=0; j<ninputs; j++) {
			//	fprintf(resfile, "%lf ", inputs[nodeinput[i]][stepbuffer]);
			//}
			fprintf(resfile,"%.15lf ", nodevalue[i][stepbuffer]);
		}
		fprintf(resfile, "\n");


		t = t + tstep;
		stepbuffer++;
		if (stepbuffer >= TAMBUFFER) {
			stepbuffer = 0;
		}

		// limpa tobeprocessed e valores nos nós (valuemtx) e atualiza a parte de SUMS
		for (i=0; i<MAXNODES; i++) {
			for (j=0; j<MAXNODES; j++) {
				tobeprocessed[i][j] = -1;
				if (i==j) {
					valuemtx[i][j] = 0;
					if (nsums[j] < -1) {
						tobeprocessed[i][j] = nsums[j];
					}
				}
				if (typeelemmtx[i][j] == TYPE_DELAY) {
					tobeprocessed[i][j] = 1;
				}
			}
		}

	}
	fprintf(resfile,"];");
	if (arq != NULL) fclose(arq);
	if (resfile!=NULL) fclose(resfile);
	getch();
}	

