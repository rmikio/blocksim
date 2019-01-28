/********************************************************************************

	BLOCKSIM eh um programa limitado que realiza algumas simulacoees de 
	sistemas dinamicos que podem ser representados por diagramas de blocos.
    Copyright (C) 2004, Renato Mikio Nakagomi

    Este arquivo eh parte do programa BLOCKSIM

    Este programa eh software livre; voce pode redistribui-lo e/ou
    modifica-lo sob os termos da Licença Pública Geral GNU, conforme
    publicada pela Free Software Foundation; tanto a versão 2 da
    Licença como (a seu criterio) qualquer versao mais nova.

    Este programa eh distribuido na expectativa de ser util, mas SEM
    QUALQUER GARANTIA; sem mesmo a garantia implicita de
    COMERCIALIZACAO ou de ADEQUACAO A QUALQUER PROPOSITO EM
    PARTICULAR. Consulte a Licenca Publica Geral GNU para obter mais
    detalhes.
 
    Voce deve ter recebido uma copia da Licenca Publica Geral GNU
    junto com este programa; se nao, escreva para a Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
    02111-1307, USA.

********************************************************************************
	BLOCKSIM is a limited software that runs some simulations of dinamic systems
	that can be represented by block diagrams.
    Copyright (C) 2004, Renato Mikio Nakagomi

    This file is part of BLOCKSIM software. 

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

********************************************************************************
	Renato Mikio Nakagomi
	rmikio@gmail.com
********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "conio.h"
// Declaracao das Bibliotecas de Funcoes de Integracao
// Declaration of Integration Functions Libraries
#include "integ_trapz.h"

// Variavel de DEBUG. Soh para imprimir mensagens de DEBUG
// DEBUG variable
#define DEBUG 1
#define DEBUGCYCLES 5

// Definicao de Constantes
// Constants definition
#define PI 3.14159265358979

// Definicao de zero
// ZERO definition
#define ZERR 0.000001

// Definicao temporaria de tamanhos aceitaveis de matrizes/vetores
// Todas as definicoes de tamanho podem (e devem) ser trocadas
// por alocacao dinamica de memoria (malloc).

// Temporary definition for acceptable sizes of matrices/vectors
// All size definitions can (and must) be replaced 
// by dynamic memory allocation (malloc).

#define MAX 10
#define MAXNODES 40
#define MAXINPUTS 10
#define MAXOUTPUTS 10
#define MAXPROBES 10
#define MAXCONSTS 10
#define MAXELEMENTS 200
#define MAXSTEPBUFFER 30
#define TAMBUFFER 15

// Definicao dos tipos de blocos existentes.
// Tais definicoes servem para o algoritmo saber
// que tipo de bloco estah presente em determinado noh.

// Definition of existing blocks
// Such definitions are useful so that the algorithm knows
// which kind of block is present in a certain node.

#define TYPE_INPUT 1
#define TYPE_OUTPUT 2
#define TYPE_CONST 3
#define TYPE_PROBE 4
#define TYPE_SUM 5
#define TYPE_GAIN 6
#define TYPE_INTEGRATOR 7
#define TYPE_DELAY 8
#define TYPE_INV 9
#define TYPE_MULTIPLY 10

// Definicao do tipo de integracao a ser realizada pelo
// bloco integrador.

// Type of integration to be executed by
// integrator block
#define INTEGTYPE_RECT 1
#define INTEGTYPE_TRAPZ 2

// Abertura do arquivo de entrada de dados
// Open input file
FILE *arq;

// Abertura do arquivo de resultados
// Open output files
FILE *resfile;
FILE *dbg;

// main
void main (void) {

	// vetor que armazena os nohs existentes no sistema
	// Este vetor nao eh ainda utilizado nesta versao.
	int nodes[MAXNODES]; 
	// matriz de valores de entrada. Utiliza o conceito de buffer circular.
	// Este vetor atualiza a matriz de valores nodais a cada iteração.
	double inputs[MAXNODES][2*MAXSTEPBUFFER];
	// matriz que armazena até 10 parâmetros para variáveis de entrada
	// Esta matriz não é ainda utilizada nesta versão.
	double inputparms[MAXNODES][10];
	// matriz de valores de saída. Utiliza o conceito de buffer circular.
	// Esta matriz não é ainda utilizada nesta versão.
	double outputs[MAXNODES][2*MAXSTEPBUFFER];
	// matriz de valores de saída dos medidores.
	// Esta matriz não é ainda utilizada nesta versão.
	double probes[MAXNODES];
	// Vetor que armazena valores das constantes presentes no sistema.
	// Este vetor atualiza a matriz de valores nodais a cada iteração.
	double consts[MAXNODES];
	// Matriz que armazena os valores de cada nó ao final de cada iteração.
	double nodevalue[MAXNODES][2*MAXSTEPBUFFER];
	// Matriz que armazena os tipos de blocos que compõem o sistema a ser simulado
	int typeelemmtx[MAXNODES][MAXNODES];
	// Matriz que armazena os valores associados a blocos e nós
	// valuemtx[i][i]: valores atuais em cada nó
	// valuemtx[i][j]: (i!=j): valores de cálculos auxiliares para cada bloco que
	// liga o nó i com o j.
	// P.ex.: se um bloco de ganho de valor 3 vai do nó 'i' ao 'j', então 
	// valuemtx[i][i] = 3. Os blocos de somatória recebem o valor '0' e os blocos 
	// de multiplicação recebem o valor '1', pois estes são os operadores neutros
	// respectivos destas funções.
	double valuemtx[MAXNODES][MAXNODES];
	// Vetor fixo que armazena a quantidade de derivações de um nó, ou seja, 
	// o número de blocos alimentados por um determinado nó. P.ex.: um INPUT no nó 'i'
	// é utilizado por 2 somadores diferentes e um ganho. Portanto nlinks[i] = 3.
	int nlinks[MAXNODES];
	// Vetor fixo que guarda a quantidade de fatores (argumentos) que deverão
	// ser processados para determinado nó. P.ex.: Um nó 'i' que resulta na somatória de
	// 3 números posui nlinks[i] = 3.
	int nfactors[MAXNODES];
	// Matriz que é atualizada a cada operação realizada e que armazena a situação de
	// cada nó processado ou a ser processado. Essa é a principal matriz do sistema que
	// define prioridades de execução de operadores.
	int tobeprocessed[MAXNODES][MAXNODES];
	// Vetor fixo que armazena onde (em quais nós) se situam todos os INPUTS
	// Este vetor ainda não é utilizado nesta versão
	int nodeinput[MAXNODES];
	// Vetor fixo que armazena onde (em quais nós) se situam todos os OUTPUTS
	// Este vetor ainda não é utilizado nesta versão
	int nodeoutput[MAXNODES];
	// Vetor fixo que armazena onde (em quais nós) se situam todos os PROBES
	// Este vetor ainda não é utilizado nesta versão
	int nodeprobe[MAXNODES];
	// Vetor fixo que armazena onde (em quais nós) se situam todos as CONSTS
	// Este vetor ainda não é utilizado nesta versão
	int nodeconst[MAXNODES];
	// Matriz que armazena qual o tipo de integrador está associado ao bloco
	// de integração. Os tipos de integradores diferem-se pelo método de integração
	// implementado para cada um.
	int integtype[MAXNODES][MAXNODES];
	// Vetor que armazena a princípio valores utilizados para a integração.
	// Este vetor poderá ser atualizado em cada passo de integração e provavelmente
	// deve ser implementado utilizando-se a teoria de buffer circular.
	double f[MAXNODES];
	// Variáveis relacionadas ao tempo. t = tempo atual em [s]. tmax = tempo total de 
	// simulação em [s]. tstep = passo de integração em [s].
	double t, tmax, tstep;
	// Passo de integração numérico relacionado ao tamanho do buffer circular implementado.
	int stepbuffer;
	// Contadores relacionados aos nós dos blocos.
	int actualnode, resultnode;
	// Contadores
	int i,j,k;
	// Flag que garante a continuidade das operações em um mesmo passo de integração
	// até que todas as operações sejam realizadas e um novo passo de integração seja
	// contabilizado.
	int processflag;
	// Variável que define qual caso será simulado.
	// Esta versão contempla a simulação de até 6 casos pré-definidos.
	// Se o usuário quiser gerar novos casos, é só respeitar a nomenclatura
	// dos arquivos já utilizados.
	int caso;
	// Strings lidas do arquivo ASCII que descreve o sistema.
	char parms1[5],parms2[5],parms3[10],parms4[16];
	// Possível expansão de número de parâmetros, dependendo das implementações
	// futuras.
	//char parms5[16];
	// O primeiro e o segundo parâmetros são referentes ao número de identificação
	// do nó utilizado (origem e destino). Ambos os valores são do tipo int.
	int parmsi1, parmsi2;
	// Variáveis que armazenam a quantidade de nós existentes no sistema.
	int nnodes;
	// O quarto parâmetro do arquivo ASCII que descreve o sistema é um valor
	// numérico.
	double parmsf4;
	// Possível expansão de número de parâmetros, dependendo das implementações
	// futuras.
	//double parmsf5;
	// Quantidade de inputs, outputs, probes e constantes.
	int ninputs, noutputs, nprobes, nconsts;


	int arqlinha, debugging, idebug, jdebug;

	dbg = fopen("debug.txt","wt");

	if (DEBUG) {
		debugging = 1;
	}

	if (debugging) {
		printf(" \n");
	}

	// Limpeza das matrizes e vetores
	if (debugging) {
		fprintf(dbg," Zerando matrizes e vetores... \n");
		fprintf(dbg," Vetor nfactors[] recebe -1.\n");
		fprintf(dbg," Matrizes typeelemmtx[][] e tobeprocessed[][] recebem -1.\n");
		fprintf(dbg," Demais vetores e matrizes recebem 0.\n");
	}
	for (i=0; i<MAXNODES; i++) {
		nlinks[i] = 0;
		nfactors[i] = -1; 
		nodeinput[i] = 0;
		nodeoutput[i] = 0;
		nodeprobe[i] = 0;
		nodes[i] = 0;
		f[i] = 0;
		for (j=0; j<MAXNODES; j++) {
			integtype[i][j] = 0;
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

	// Escolha de caso a ser simulado.
	// Se o usuário quiser criar seu próprio caso, ele pode
	// utulizar um arquivo de mesmo nome de qquer um dos
	// arquivos listados abaixo.
	printf(" Caso a ser simulado (2,3,4,5,6,7): \n");
	scanf("%d", &caso);

	if (caso == 2) {
		arq = fopen("sist02.bds","rt");
		resfile = fopen("result02.m","wt");
	}
	if (caso == 3) {
		arq = fopen("sist03.bds","rt");
		resfile = fopen("result03.m","wt");
	}
	if (caso == 4) {
		arq = fopen("sist04.bds","rt");
		resfile = fopen("result04.m","wt");
	}
	if (caso == 5) {
		arq = fopen("sist05.bds","rt");
		resfile = fopen("result05.m","wt");
	}
	if (caso == 6) {
		arq = fopen("sist06.bds","rt");
		resfile = fopen("result06.m","wt");
	}
	if (caso == 7) {
		arq = fopen("sist07.bds","rt");
		resfile = fopen("result07.m","wt");
	}
	fprintf(resfile, "saida = [");

    // Gerenciador de erro para bertura de arquivo para leitura. 
	if(arq == NULL ) {
      printf( "O arquivo não foi aberto.\n" );
	  exit(0);
	}
	else {
      printf( "O arquivo foi aberto.\n" );
	}
	
	// Leitura da primeira linha: duração da simulação e passo de integração
	// em segundos.
	fscanf(arq,"%s",parms1);
	tmax = atof(parms1);
	fscanf(arq,"%s",parms2);
	tstep = atof(parms2);

	// Zera contadores
	if (debugging) {
		fprintf(dbg," Contadores auxiliares zerados.\n");
	}
	ninputs = 0;
	noutputs = 0;
	nprobes = 0;
	nconsts = 0;
	nnodes = 0;

	arqlinha = 1;
	while (!feof(arq)) {

		if (debugging) {
			fprintf(dbg," Lendo linha %d do arquivo.\n",arqlinha);
		}
		// leitura do arquivo
		if (debugging) {
			fprintf(dbg," Lendo parametro 1 da linha %d do arquivo.\n",arqlinha);
		}
		fscanf(arq,"%s",parms1);
		parmsi1 = atoi(parms1);

		if (debugging) {
			fprintf(dbg," Lendo parametro 2 da linha %d do arquivo.\n",arqlinha);
		}
		fscanf(arq,"%s",parms2);
		parmsi2 = atoi(parms2);

		if (debugging) {
			fprintf(dbg," Lendo parametro 3 da linha %d do arquivo.\n",arqlinha);
		}
		fscanf(arq,"%s",parms3);

		if (debugging) {
			fprintf(dbg," Lendo parametro 4 da linha %d do arquivo.\n",arqlinha);
		}
		fscanf(arq,"%s",parms4);
		parmsf4 = atof(parms4);

		// contagem dos nohs existentes
		if (nodes[parmsi1] == 0) {
			nodes[parmsi1]++;
			nnodes++;
			if (debugging) {
				fprintf(dbg," Marcando a existencia do noh %d.\n", parmsi1);
				fprintf(dbg," Somando mais um noh ao total: TOTAL %d nohs \n",nnodes);
			}
		}
		if (nodes[parmsi2] == 0) {
			nodes[parmsi2]++;
			nnodes++;
			if (debugging) {
				fprintf(dbg," Marcando a existencia do noh %d.\n", parmsi1);
				fprintf(dbg," Somando mais um noh ao total: TOTAL %d nohs \n",nnodes);
			}
		}

		// Gerando as condições iniciais do sistema.
		if (debugging) {
			fprintf(dbg,"-> Gerando condicoes iniciais do sistema.\n");
		}
		// se a função é de SOMATÓRIA (SUM)
		if (strcmp(parms3,"SUM")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_SUM;
			valuemtx[parmsi1][parmsi2] = 0;
			// o somador pode ser processado
			tobeprocessed[parmsi1][parmsi2] = 1; 
			// acumula ptos negativos para liberar o processamento do resultado do somador
			// a cada fator existente, subtrai-se um ponto .
			tobeprocessed[parmsi2][parmsi2]--; 
			nlinks[parmsi1]++;
			// acumulador (negativo) de fatores do somador
			nfactors[parmsi2]--;
			if (debugging) {
				fprintf(dbg,"...Achada funcao de SOMATORIA.\n");
				fprintf(dbg,"...Valor typeelemmtx[%d][%d] = TYPE_SUM <- Marcando tipo de elemento. \n",parmsi1,parmsi2);
				fprintf(dbg,"...Valor valuemtx[%d][%d] = 0 <- Condicao inicial de somatoria\n", parmsi1,parmsi2);
				fprintf(dbg,"...Valor tobeprocessed[%d][%d] = 1 <- SOMADOR pode ser processado\n",parmsi1,parmsi2);
				fprintf(dbg,"...Valor tobeprocessed[%d][%d] = %d <- Subtraindo um ponto para cada fator (input para somador) existente\n",parmsi2,parmsi2,tobeprocessed[parmsi2][parmsi2]);
				fprintf(dbg,"...nlinks[%d] = %d <- contabiliza numero de links no noh %d \n",parmsi1,nlinks[parmsi1],parmsi1);
				fprintf(dbg,"...nfactors[%d] = %d <- numero (negativado) de fatores do somador (permissao da saida do somador)\n",parmsi2,nfactors[parmsi2]);
			}
		}

		// se a função é de MULTIPLICAÇÃO (MULTIPLY)
		if (strcmp(parms3,"MULTIPLY")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_MULTIPLY;
			valuemtx[parmsi1][parmsi2] = 1;
			// o multiplicador pode ser processado
			tobeprocessed[parmsi1][parmsi2] = 1;
			// acumula ptos negativos para liberar o processamento do resultado do multiplicador
			// a cada fator existente, subtrai-se um ponto .
			tobeprocessed[parmsi2][parmsi2]--; 
			nlinks[parmsi1]++;
			// acumulador (negativo) de fatores do multiplicador
			nfactors[parmsi2]--;
			if (debugging) {
				fprintf(dbg,"...Achada funcao de MULTIPLICACAO.\n");
				fprintf(dbg,"...Valor typeelemmtx[%d][%d] = TYPE_MULTIPLY <- Marcando tipo de elemento. \n",parmsi1,parmsi2);
				fprintf(dbg,"...Valor valuemtx[%d][%d] = 1 <- Condicao inicial de produtorio\n", parmsi1,parmsi2);
				fprintf(dbg,"...Valor tobeprocessed[%d][%d] = 1 <- MULTIPLICADOR pode ser processado\n",parmsi1,parmsi2);
				fprintf(dbg,"...Valor tobeprocessed[%d][%d] = %d <- Subtraindo um ponto para cada fator (input para produtorio) existente\n",parmsi2,parmsi2,tobeprocessed[parmsi2][parmsi2]);
				fprintf(dbg,"...nlinks[%d] = %d <- contabiliza numero de links no noh %d \n",parmsi1,nlinks[parmsi1],parmsi1);
				fprintf(dbg,"...nfactors[%d] = %d <- numero (negativado) de fatores do produtorio (permissao da saida do produtorio)\n",parmsi2,nfactors[parmsi2]);
			}
		}

		
		// se é função de GANHO (GAIN)
		if (strcmp(parms3,"GAIN")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_GAIN;
			// armazena o valor do ganho
			valuemtx[parmsi1][parmsi2] = parmsf4;
			// o bloco de ganho pode ser processado
			tobeprocessed[parmsi1][parmsi2] = 1;
			nlinks[parmsi1]++;
			if (debugging) {
				fprintf(dbg,"...Achada funcao de GANHO.\n");
				fprintf(dbg,"...Valor typeelemmtx[%d][%d] = TYPE_GAIN <- Marcando tipo de elemento. \n",parmsi1,parmsi2);
				fprintf(dbg,"...Valor valuemtx[%d][%d] = %f <- Valor do GANHO\n", parmsi1,parmsi2,parmsf4);
				fprintf(dbg,"...Valor tobeprocessed[%d][%d] = 1 <- GANHO pode ser processado\n",parmsi1,parmsi2);
				fprintf(dbg,"...nlinks[%d]++ = %d <- contabiliza numero de links no noh %d \n",parmsi1,nlinks[parmsi1],parmsi1);
			}
		}

		// se é função de INV (INV)
		if (strcmp(parms3,"INV")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_INV;
			// O quarto parâmetro de INV não é utilizado
			//valuemtx[parmsi1][parmsi2] = parmsf4;
			// O bloco de inversão pode ser processado
			tobeprocessed[parmsi1][parmsi2] = 1;
			nlinks[parmsi1]++;
			if (debugging) {
				fprintf(dbg,"...Achada funcao de INVERSA.\n");
				fprintf(dbg,"...Valor typeelemmtx[%d][%d] = TYPE_INV <- Marcando tipo de elemento. \n",parmsi1,parmsi2);
				fprintf(dbg,"...Valor tobeprocessed[%d][%d] = 1 <- INV pode ser processado\n",parmsi1,parmsi2);
				fprintf(dbg,"...nlinks[%d]++ = %d <- contabiliza numero de links no noh %d \n",parmsi1,nlinks[parmsi1],parmsi1);
			}
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
			if (debugging) {
				fprintf(dbg,"...Achada variavel de INPUT.\n");
				fprintf(dbg,"...Valor typeelemmtx[%d][%d] = TYPE_INPUT <- Marcando tipo de elemento. \n",parmsi1,parmsi2);
				fprintf(dbg,"...Valor valuemtx[%d][%d] = %f <- Condicao inicial de INPUT\n", parmsi1,parmsi2,parmsf4);
				fprintf(dbg,"...Valor inputs[%d][0] = %f <- condicao inicial de INPUT no noh %d \n",parmsi1,parmsf4,parmsi1);
				fprintf(dbg,"...Valor inputparms[%d][0] = %f <- guarda parametro de INPUT do noh %d \n",parmsi1,parmsf4,parmsi1);
				fprintf(dbg,"...nodeinput[%d] = %d <- guarda noh de input \n",ninputs, parmsi1);
				fprintf(dbg,"...ninputs++ = %d <- contabilizando inputs \n",ninputs);
			}
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
			if (debugging) {
				fprintf(dbg,"...Achada CONSTANTE.\n");
				fprintf(dbg,"...Valor typeelemmtx[%d][%d] = TYPE_CONST <- Marcando tipo de elemento. \n",parmsi1,parmsi2);
				fprintf(dbg,"...nodeconsts[%d] = %d <- guarda noh de constante \n",nconsts, parmsi1);
				fprintf(dbg,"...nconsts++ = %d <- contabilizando constantes \n",nconsts);
			}
		}
		// se é variável de INTEGRATOR
		if (strcmp(parms3,"INTEGRATOR")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_INTEGRATOR;
			integtype[parmsi1][parmsi2] = (int)(parmsf4);
			// O integrador pode ser processado
			tobeprocessed[parmsi1][parmsi2] = 1;
			nlinks[parmsi1]++;
			if (debugging) {
				fprintf(dbg,"...Achada funcao de INTEGRADOR.\n");
				fprintf(dbg,"...Valor typeelemmtx[%d][%d] = TYPE_INTEGRATOR <- Marcando tipo de elemento. \n",parmsi1,parmsi2);
				fprintf(dbg,"...Valor integtype[%d][%d] = %d <- Tipo do metodo de integracao a ser utilizado\n", parmsi1,parmsi2,(int)(parmsf4));
				fprintf(dbg,"...Valor tobeprocessed[%d][%d] = 1 <- INTEGRADOR pode ser processado\n",parmsi1,parmsi2);
				fprintf(dbg,"...nlinks[%d] = %d <- contabiliza numero de links no noh %d \n",parmsi1,nlinks[parmsi1],parmsi1);
			}
		}
		// se é variável de DELAY
		if (strcmp(parms3,"DELAY")==0) {
			typeelemmtx[parmsi1][parmsi2] = TYPE_DELAY;
			valuemtx[parmsi1][parmsi2] = parmsf4;
			// O atraso pode ser processado
			tobeprocessed[parmsi1][parmsi2] = 1;
			nlinks[parmsi1]++;
		}
		arqlinha++;
	}
	
	

	// Aqui começa a análise dos blocos
	// Tentativa 3: análise por nós / matriz de relacionamentos
	
	t = 0;
	stepbuffer = 1; // começa de 1, pois existe a condição inicial 0
	if (debugging) {
		fprintf(dbg," ********** Inicializando t... t = 0 **********\n");
		fprintf(dbg," Buffer -> stepbuffer = 1 <- a condicao iniciada foi gravada em stepbuffer = ZERO\n");
		fprintf(dbg,">>> valuemtx:\n");
		for (idebug=0; idebug<nnodes; idebug++) {
			for (jdebug=0; jdebug<nnodes; jdebug++) {
				fprintf(dbg," %09.4lf ", valuemtx[idebug][jdebug]);
				if (jdebug == nnodes - 1) {
					fprintf(dbg,"\n");
				}
			}
		}
	}
	if (debugging) {
		fprintf(dbg,">>> tobeprocessed:\n");
		for (idebug=0; idebug<nnodes; idebug++) {
			for (jdebug=0; jdebug<nnodes; jdebug++) {
				fprintf(dbg," %03d ", tobeprocessed[idebug][jdebug]);
				if (jdebug == nnodes - 1) {
					fprintf(dbg,"\n");
				}
			}
		}
	}
	while (t<tmax) {
		processflag = 0;
		if (debugging) {
			fprintf(dbg," Passo t = %f [s] \n", t);
			fprintf(dbg," Flag de processo = 0 <- inicio do passo de simulacao. Nada pode ser processado por enquanto.\n");
		}
		// Transferir INPUTS e CONSTS se existirem
		if (debugging) {
			fprintf(dbg," Calculando INPUTS baseados no parametro de amplitude fornecido  \n");
		}
		for (i=0; i<ninputs; i++) {
			inputs[nodeinput[i]][stepbuffer] = inputparms[nodeinput[i]][0]*sin(2*PI*60*t);
			if (debugging) {
				fprintf(dbg,"-> INPUT %d @ t = %f : amplitude*sen(2*pi*f*t) = %f \n",i,t,inputs[nodeinput[i]][stepbuffer]);
				fprintf(dbg,"   Este valor eh inserido na variavel inputs[nodeinput[i]][stepbuffer]. nodeinput[i] fornece o noh onde INPUT eh aplicado. \n");
			}
		}
		if (debugging) {
			fprintf(dbg,"Inicia laco FOR para insercao de INPUTS e CONSTS. Sao estes valores que permitem o inicio da execucao dos demais processamentos. \n");
		}
		for (actualnode = 0; actualnode<nnodes; actualnode++) {
			// Atualiza INPUT
			if (typeelemmtx[actualnode][actualnode] == TYPE_INPUT) {
				valuemtx[actualnode][actualnode] = inputs[actualnode][stepbuffer];
				tobeprocessed[actualnode][actualnode] = nlinks[actualnode];
				if (debugging) {
					fprintf(dbg,">>> valuemtx:\n");
					for (idebug=0; idebug<nnodes; idebug++) {
						for (jdebug=0; jdebug<nnodes; jdebug++) {
							fprintf(dbg," %09.4lf ", valuemtx[idebug][jdebug]);
							if (jdebug == nnodes - 1) {
								fprintf(dbg,"\n");
							}
						}
					}
				}

				if (debugging) {
					fprintf(dbg,">>> tobeprocessed:\n");
					for (idebug=0; idebug<nnodes; idebug++) {
						for (jdebug=0; jdebug<nnodes; jdebug++) {
							fprintf(dbg," %+d ", tobeprocessed[idebug][jdebug]);
							if (jdebug == nnodes - 1) {
								fprintf(dbg,"\n");
							}
						}
					}
				}
				processflag = 1;
				if (debugging) {
					fprintf(dbg," Atualiza matriz valuemtx[%d][%d] = valor de input no noh %d = %f \n", actualnode, actualnode, actualnode,inputs[actualnode][stepbuffer]);
					fprintf(dbg," Atualiza matriz tobeprocessed[%d][%d] = nlinks[%d] = %d <- o noh jah possui um valor. A matriz tobeprocessed recebe o numero de links a serem processados.\n",actualnode,actualnode,actualnode,nlinks[actualnode]);
					fprintf(dbg," Flag de processo = 1 <- pelo menos um noh tem valor. Talvez seja possivel executar uma operacao com este valor.\n");
				}
			}
			// Atualiza CONST
			if (typeelemmtx[actualnode][actualnode] == TYPE_CONST) {
				valuemtx[actualnode][actualnode] = consts[actualnode];
				tobeprocessed[actualnode][actualnode] = nlinks[actualnode];
				if (debugging) {
					fprintf(dbg,">>> valuemtx:\n");
					for (idebug=0; idebug<nnodes; idebug++) {
						for (jdebug=0; jdebug<nnodes; jdebug++) {
							fprintf(dbg," %09.4lf ", valuemtx[idebug][jdebug]);
							if (jdebug == nnodes - 1) {
								fprintf(dbg,"\n");
							}
						}
					}
				}

				if (debugging) {
					fprintf(dbg,">>> tobeprocessed:\n");
					for (idebug=0; idebug<nnodes; idebug++) {
						for (jdebug=0; jdebug<nnodes; jdebug++) {
							fprintf(dbg," %+d ", tobeprocessed[idebug][jdebug]);
							if (jdebug == nnodes - 1) {
								fprintf(dbg,"\n");
							}
						}
					}
				}

				processflag = 1;
				if (debugging) {
					fprintf(dbg," Atualiza matriz valuemtx[%d][%d] = valor de constante no noh %d = %f \n", actualnode, actualnode, actualnode,inputs[actualnode][stepbuffer]);
					fprintf(dbg," Atualiza matriz tobeprocessed[%d][%d] = nlinks[%d] = %d <- o noh jah possui um valor. A matriz tobeprocessed recebe o numero de links a serem processados.\n",actualnode,actualnode,actualnode,nlinks[actualnode]);
					fprintf(dbg," Flag de processo = 1 <- pelo menos um noh tem valor. Talvez seja possivel executar uma operacao com este valor.\n");
				}
			}
			for (resultnode = 0; resultnode<nnodes; resultnode++) {
				// Atualiza somatória
				if (typeelemmtx[actualnode][resultnode] == TYPE_SUM) {
					tobeprocessed[actualnode][resultnode] = 1;
					tobeprocessed[resultnode][resultnode] = nfactors[resultnode];
					if (debugging) {
						fprintf(dbg,">>> tobeprocessed:\n");
						for (idebug=0; idebug<nnodes; idebug++) {
							for (jdebug=0; jdebug<nnodes; jdebug++) {
								fprintf(dbg," %+d ", tobeprocessed[idebug][jdebug]);
								if (jdebug == nnodes - 1) {
									fprintf(dbg,"\n");
								}
							}
						}
					}

					//processflag = 1;
					if (debugging) {
						fprintf(dbg," SOMADOR: \n");
						fprintf(dbg," Atualiza matriz tobeprocessed[%d][%d] = 1 \n",actualnode,resultnode);
						fprintf(dbg," Atualiza matriz tobeprocessed[%d][%d] = nfactors[resultnode] = %d <- noh resultante do somador possui n-fatores a serem utilizados\n",resultnode,resultnode,nfactors[resultnode]);
					}
				}
				// Atualiza multiplicador
				if (typeelemmtx[actualnode][resultnode] == TYPE_MULTIPLY) {
					tobeprocessed[actualnode][resultnode] = 1;
					tobeprocessed[resultnode][resultnode] = nfactors[resultnode];
					valuemtx[resultnode][resultnode] = 1;
					if (debugging) {
						fprintf(dbg,">>> valuemtx:\n");
						for (idebug=0; idebug<nnodes; idebug++) {
							for (jdebug=0; jdebug<nnodes; jdebug++) {
								fprintf(dbg," %09.4lf ", valuemtx[idebug][jdebug]);
								if (jdebug == nnodes - 1) {
									fprintf(dbg,"\n");
								}
							}
						}
					}

					if (debugging) {
						fprintf(dbg,">>> tobeprocessed:\n");
						for (idebug=0; idebug<nnodes; idebug++) {
							for (jdebug=0; jdebug<nnodes; jdebug++) {
								fprintf(dbg," %+d ", tobeprocessed[idebug][jdebug]);
								if (jdebug == nnodes - 1) {
									fprintf(dbg,"\n");
								}
							}
						}
					}

					//processflag = 1;
					if (debugging) {
						fprintf(dbg," MULTIPLICADOR: \n");
						fprintf(dbg," Atualiza matriz tobeprocessed[%d][%d] = 1 \n",actualnode,resultnode);
						fprintf(dbg," Atualiza matriz tobeprocessed[%d][%d] = nfactors[resultnode] = %d <- noh resultante do multiplicador possui n-fatores a serem utilizados\n",resultnode,resultnode,nfactors[resultnode]);
					}
				}
				// Atualiza atraso
				if (typeelemmtx[actualnode][resultnode] == TYPE_DELAY) {
					tobeprocessed[actualnode][resultnode] = 1;

					processflag = 1;
					if (debugging) {
						fprintf(dbg," DELAY: \n");
						fprintf(dbg," Atualiza matriz tobeprocessed[%d][%d] = 1 \n",actualnode,resultnode);
						fprintf(dbg," Flag de processo = 1 <- DELAY gera valor usando condicao inicial ou valores historicos (buffer). Talvez seja possivel executar uma operacao com este valor.\n");
					}
				}
				// Atualiza integrador
				if (typeelemmtx[actualnode][resultnode] == TYPE_INTEGRATOR) {
					tobeprocessed[actualnode][resultnode] = 1;
					processflag = 1;
					if (debugging) {
						fprintf(dbg," INTEGRADOR: \n");
						fprintf(dbg," Atualiza matriz tobeprocessed[%d][%d] = 1 \n",actualnode,resultnode);
						fprintf(dbg," Flag de processo = 1 <- INTEGRADOR gera valor usando condicao inicial ou valores historicos (buffer). Talvez seja possivel executar uma operacao com este valor.\n");
					}
				}
			}
		}
		if (debugging) {
			fprintf(dbg," Inicio do processamento para o passo de integracao atual. \n");
			fprintf(dbg," Flag de processamento deve ser 1 para iniciar o processo. \n");
		}
		// Início do processamento para o passo de integração atual
		while (processflag) {
			// Flag de processamento é zerado, pois espera-se que
			// o mesmo seja reativado enquanto existirem blocos a 
			// serem processados.
			processflag = 0;
			if (debugging) {
				fprintf(dbg," Dentro do WHILE. Executar enquanto houver necessidade de processamento.\n");
				fprintf(dbg," Flag de processamento eh ZERADO.\n");
			}

			for (actualnode = 0; actualnode<nnodes; actualnode++) {
				for (resultnode = 0; resultnode<nnodes; resultnode++) {
					// GAIN
					if (typeelemmtx[actualnode][resultnode] == TYPE_GAIN) {
						// Se o bloco pode ser processado e o nó de entrada do bloco está
						// liberado para processamento
						if (tobeprocessed[actualnode][resultnode] !=0 && tobeprocessed[actualnode][actualnode] > 0) {
							// O nó resultado é liberado para processamento.
							// Este nó resultado possui X processamentos a dever
							tobeprocessed[resultnode][resultnode] = nlinks[resultnode];
							// Subtrai-se um processo a dever do nó de entrada 
							tobeprocessed[actualnode][actualnode]--;
							// O bloco foi processado.
							tobeprocessed[actualnode][resultnode] = 0;
							// Calcula-se o valor do resultado
							valuemtx[resultnode][resultnode] = valuemtx[actualnode][resultnode]*valuemtx[actualnode][actualnode];
							// Reativa-se o flag de processamento
							processflag = 1;

							if (debugging) {
								fprintf(dbg,"Processando elemento GANHO em [%d][%d]\n", actualnode, resultnode);
								fprintf(dbg,"// O nó resultado é liberado para processamento. \n");
								fprintf(dbg,"// Este nó resultado possui X processamentos a dever \n");
								fprintf(dbg,"tobeprocessed[resultnode][resultnode] = nlinks[resultnode]; \n");
								fprintf(dbg,"// Subtrai-se um processo a dever do nó de entrada  \n");
								fprintf(dbg,"tobeprocessed[actualnode][actualnode]--; \n");
								fprintf(dbg,"// O bloco foi processado. \n");
								fprintf(dbg,"tobeprocessed[actualnode][resultnode] = 0; \n");
								fprintf(dbg,"// Calcula-se o valor do resultado \n");
								fprintf(dbg,"valuemtx[resultnode][resultnode] = valuemtx[actualnode][resultnode]*valuemtx[actualnode][actualnode]; \n");
								fprintf(dbg,"// Reativa-se o flag de processamento \n");
								fprintf(dbg,"processflag = 1; \n");
							}
							if (debugging) {
								fprintf(dbg,">>> valuemtx:\n");
								for (idebug=0; idebug<nnodes; idebug++) {
									for (jdebug=0; jdebug<nnodes; jdebug++) {
										fprintf(dbg," %09.4lf ", valuemtx[idebug][jdebug]);
										if (jdebug == nnodes - 1) {
											fprintf(dbg,"\n");
										}
									}
								}
							}

							if (debugging) {
								fprintf(dbg,">>> tobeprocessed:\n");
								for (idebug=0; idebug<nnodes; idebug++) {
									for (jdebug=0; jdebug<nnodes; jdebug++) {
										fprintf(dbg," %+d ", tobeprocessed[idebug][jdebug]);
										if (jdebug == nnodes - 1) {
											fprintf(dbg,"\n");
										}
									}
								}
							}
						}
					} // GAIN
					
					// INV
					if (typeelemmtx[actualnode][resultnode] == TYPE_INV) {
						// Se o bloco pode ser processado e o nó de entrada do bloco está
						// liberado para processamento
						if (tobeprocessed[actualnode][resultnode] !=0 && tobeprocessed[actualnode][actualnode] > 0) {
							// O nó resultado é liberado para processamento.
							// Este nó resultado possui X processamentos a dever
							tobeprocessed[resultnode][resultnode] = nlinks[resultnode];
							// Subtrai-se um processo a dever do nó de entrada 
							tobeprocessed[actualnode][actualnode]--;
							// O bloco foi processado.
							tobeprocessed[actualnode][resultnode] = 0;
							// Atenção para divisão por zero
							if (valuemtx[actualnode][actualnode] == 0) {
								valuemtx[actualnode][actualnode] = ZERR;
							}
							// Calcula-se o valor do resultado
							valuemtx[resultnode][resultnode] = 1/valuemtx[actualnode][actualnode];
							// Reativa-se o flag de processamento
							processflag = 1;
							if (debugging) {
								fprintf(dbg,"Processando elemento INV em [%d][%d]\n", actualnode, resultnode);
							}
						}
					} // GAIN

					// SUM
					if (typeelemmtx[actualnode][resultnode] == TYPE_SUM) {
						// Se o bloco pode ser processado e o nó de entrada do bloco está
						// liberado para processamento
						if (tobeprocessed[actualnode][resultnode] !=0 && tobeprocessed[actualnode][actualnode] > 0) {
							// Calcula-se o valor do resultado parcial da somatória
							valuemtx[resultnode][resultnode] = valuemtx[resultnode][resultnode]+valuemtx[actualnode][actualnode];
							// Subtrai-se um processo a dever do nó de entrada 
							tobeprocessed[actualnode][actualnode]--;
							// O bloco foi processado.
							tobeprocessed[actualnode][resultnode] = 0;
							// Soma-se um processo realizado para o nó resultado.
							tobeprocessed[resultnode][resultnode]++;
							// Reativa-se o flag de processamento
							processflag = 1;
							// O nó resultado é liberado para processamento somente se
							// o valor de processos atingir -1.
							if (tobeprocessed[resultnode][resultnode]==-1) {
								// Este nó resultado possui X processamentos a dever
								tobeprocessed[resultnode][resultnode] = nlinks[resultnode];
							}
							if (debugging) {
								fprintf(dbg,"Processando elemento SOMATORIA em [%d][%d]\n", actualnode, resultnode);
								fprintf(dbg,"// Calcula-se o valor do resultado parcial da somatória \n");
								fprintf(dbg,"valuemtx[resultnode][resultnode] = valuemtx[resultnode][resultnode]+valuemtx[actualnode][actualnode]; \n");
								fprintf(dbg,"// Subtrai-se um processo a dever do nó de entrada  \n");
								fprintf(dbg,"tobeprocessed[actualnode][actualnode]--; \n");
								fprintf(dbg,"// O bloco foi processado. \n");
								fprintf(dbg,"tobeprocessed[actualnode][resultnode] = 0; \n");
								fprintf(dbg,"// Soma-se um processo realizado para o nó resultado. \n");
								fprintf(dbg,"tobeprocessed[resultnode][resultnode]++; \n");
								fprintf(dbg,"// Reativa-se o flag de processamento \n");
								fprintf(dbg,"processflag = 1; \n");
								fprintf(dbg,"// O nó resultado é liberado para processamento somente se \n");
								fprintf(dbg,"// o valor de processos atingir -1. \n");
								fprintf(dbg,"if (tobeprocessed[resultnode][resultnode]==-1) { \n");
								fprintf(dbg,"	// Este nó resultado possui X processamentos a dever \n");
								fprintf(dbg,"	tobeprocessed[resultnode][resultnode] = nlinks[resultnode]; \n");
								fprintf(dbg,"} \n");
							}
							if (debugging) {
								fprintf(dbg,">>> valuemtx:\n");
								for (idebug=0; idebug<nnodes; idebug++) {
									for (jdebug=0; jdebug<nnodes; jdebug++) {
										fprintf(dbg," %09.4lf ", valuemtx[idebug][jdebug]);
										if (jdebug == nnodes - 1) {
											fprintf(dbg,"\n");
										}
									}
								}
							}

							if (debugging) {
								fprintf(dbg,">>> tobeprocessed:\n");
								for (idebug=0; idebug<nnodes; idebug++) {
									for (jdebug=0; jdebug<nnodes; jdebug++) {
										fprintf(dbg," %+d ", tobeprocessed[idebug][jdebug]);
										if (jdebug == nnodes - 1) {
											fprintf(dbg,"\n");
										}
									}
								}
							}

						}
					} // SUM

					// MULTIPLY
					if (typeelemmtx[actualnode][resultnode] == TYPE_MULTIPLY) {
						// Se o bloco pode ser processado e o nó de entrada do bloco está
						// liberado para processamento
						if (tobeprocessed[actualnode][resultnode] !=0 && tobeprocessed[actualnode][actualnode] > 0) {
							// Calcula-se o valor do resultado parcial da multiplicação
							valuemtx[resultnode][resultnode] = valuemtx[resultnode][resultnode]*valuemtx[actualnode][actualnode];
							// Subtrai-se um processo a dever do nó de entrada 
							tobeprocessed[actualnode][actualnode]--;
							// O bloco foi processado.
							tobeprocessed[actualnode][resultnode] = 0;
							// Soma-se um processo realizado para o nó resultado.
							tobeprocessed[resultnode][resultnode]++;
							// Reativa-se o flag de processamento
							processflag = 1;
							// O nó resultado é liberado para processamento somente se
							// atingir o valor -1.
							if (tobeprocessed[resultnode][resultnode]==-1) {
								// Este nó resultado possui X processamentos a dever
								tobeprocessed[resultnode][resultnode] = nlinks[resultnode];
							}
							if (debugging) {
								fprintf(dbg,"Processando elemento PRODUTORIO em [%d][%d]\n", actualnode, resultnode);
								fprintf(dbg,"// Calcula-se o valor do resultado parcial da multiplicação \n");
								fprintf(dbg,"valuemtx[resultnode][resultnode] = valuemtx[resultnode][resultnode]*valuemtx[actualnode][actualnode]; \n");
								fprintf(dbg,"// Subtrai-se um processo a dever do nó de entrada  \n");
								fprintf(dbg,"tobeprocessed[actualnode][actualnode]--; \n");
								fprintf(dbg,"// O bloco foi processado. \n");
								fprintf(dbg,"tobeprocessed[actualnode][resultnode] = 0; \n");
								fprintf(dbg,"// Soma-se um processo realizado para o nó resultado. \n");
								fprintf(dbg,"tobeprocessed[resultnode][resultnode]++; \n");
								fprintf(dbg,"// Reativa-se o flag de processamento \n");
								fprintf(dbg,"processflag = 1; \n");
								fprintf(dbg,"// O nó resultado é liberado para processamento somente se \n");
								fprintf(dbg,"// atingir o valor -1. \n");
								fprintf(dbg,"if (tobeprocessed[resultnode][resultnode]==-1) { \n");
								fprintf(dbg,"	// Este nó resultado possui X processamentos a dever \n");
								fprintf(dbg,"	tobeprocessed[resultnode][resultnode] = nlinks[resultnode]; \n");
								fprintf(dbg,"} \n");
							}
														if (debugging) {
								fprintf(dbg,">>> valuemtx:\n");
								for (idebug=0; idebug<nnodes; idebug++) {
									for (jdebug=0; jdebug<nnodes; jdebug++) {
										fprintf(dbg," %09.4lf ", valuemtx[idebug][jdebug]);
										if (jdebug == nnodes - 1) {
											fprintf(dbg,"\n");
										}
									}
								}
							}

							if (debugging) {
								fprintf(dbg,">>> tobeprocessed:\n");
								for (idebug=0; idebug<nnodes; idebug++) {
									for (jdebug=0; jdebug<nnodes; jdebug++) {
										fprintf(dbg," %+d ", tobeprocessed[idebug][jdebug]);
										if (jdebug == nnodes - 1) {
											fprintf(dbg,"\n");
										}
									}
								}
							}

						}
					} // MULTIPLY

					
					// DELAY
					if (typeelemmtx[actualnode][resultnode] == TYPE_DELAY) {
						// Se o bloco pode ser processado 
						// (DELAY não precisa ter o nó de entrada liberado
						// uma vez que lida com valores históricos)
						if (tobeprocessed[actualnode][resultnode] > 0) {
							// O bloco foi processado.
							tobeprocessed[actualnode][resultnode] = 0;
							// O nó resultado é liberado para processamento.
							// Este nó resultado possui X processamentos a dever
							tobeprocessed[resultnode][resultnode] = nlinks[resultnode];
							// Calcula-se o valor do resultado
							valuemtx[resultnode][resultnode] = nodevalue[actualnode][TAMBUFFER+stepbuffer-1];
							// Reativa-se o flag de processamento
							processflag = 1;
						}
							if (debugging) {
								fprintf(dbg,"Processando elemento ATRASO em [%d][%d]\n", actualnode, resultnode);
							}
					} // DELAY
					
					// INTEGRATOR
					if (typeelemmtx[actualnode][resultnode] == TYPE_INTEGRATOR) {
						// Se o bloco pode ser processado 
						// (INTEGRATOR não precisa ter o nó de entrada liberado
						// uma vez que lida com valores históricos)
						if (tobeprocessed[actualnode][resultnode] > 0) {
							// O nó resultado é liberado para processamento.
							// Este nó resultado possui X processamentos a dever
							tobeprocessed[resultnode][resultnode] = nlinks[resultnode];
							// O bloco foi processado.
							tobeprocessed[actualnode][resultnode] = 0;
							if (tobeprocessed[actualnode][resultnode] <= 0) { 
							// não conheço o valor atual a ser integrado
							// tentativa de estimativa
							// nodevalue[actualnode][TAMBUFFER+stepbuffer] é usado como
							// variável temporária
								if (t >= 2*tstep) {
									nodevalue[actualnode][TAMBUFFER+stepbuffer] = ((nodevalue[actualnode][TAMBUFFER+stepbuffer-1]-nodevalue[actualnode][TAMBUFFER+stepbuffer-2])/((t-tstep)-(t-2*tstep)))*(t-(t-tstep)) + nodevalue[actualnode][TAMBUFFER+stepbuffer-1];
									//nodevalue[actualnode][TAMBUFFER+stepbuffer] = nodevalue[actualnode][TAMBUFFER+stepbuffer-1];
								}
								else {
									nodevalue[actualnode][TAMBUFFER+stepbuffer] = 0;
								}
							}
							// Cálculo da integral (experimental) como média de
							// valores de derivadas
							//if (integtype[actualnode][resultnode] == INTEGTYPE_TRAPZ) {
							//	f[0] = nodevalue[actualnode][TAMBUFFER+stepbuffer-1];
							//	f[1] = nodevalue[actualnode][TAMBUFFER+stepbuffer];
							//	valuemtx[resultnode][resultnode] = nodevalue[resultnode][TAMBUFFER+stepbuffer-1]+Trapezoidal_Rule_Tab_Sum_RL(tstep,1,f);
							//	f[0] = nodevalue[actualnode][TAMBUFFER+stepbuffer-2];
							//	f[1] = nodevalue[actualnode][TAMBUFFER+stepbuffer-1];
							//	valuemtx[resultnode][resultnode] = (valuemtx[resultnode][resultnode] + nodevalue[resultnode][TAMBUFFER+stepbuffer-1]+Trapezoidal_Rule_Tab_Sum_RL(tstep,1,f))/2;
							//}
							
							// Cálculo da integral trapezoidal normal via função
							//if (integtype[actualnode][resultnode] == INTEGTYPE_TRAPZ) {
							//	f[0] = nodevalue[actualnode][TAMBUFFER+stepbuffer-1];
							//	f[1] = nodevalue[actualnode][TAMBUFFER+stepbuffer];
							//	valuemtx[resultnode][resultnode] = nodevalue[resultnode][TAMBUFFER+stepbuffer-1]+Trapezoidal_Rule_Tab_Sum_RL(tstep,1,f);
							//}

							// Cálculo da integral trapezoidal normal via implementação na mão
							if (integtype[actualnode][resultnode] == INTEGTYPE_TRAPZ) {
								f[0] = nodevalue[actualnode][TAMBUFFER+stepbuffer-1];
								f[1] = nodevalue[actualnode][TAMBUFFER+stepbuffer];
								if (t<=tstep) {
									valuemtx[resultnode][resultnode] = nodevalue[resultnode][TAMBUFFER+stepbuffer-1]+(f[0]+f[1])*tstep;
								}
								else {
									valuemtx[resultnode][resultnode] = nodevalue[resultnode][TAMBUFFER+stepbuffer-1]+(f[0]+f[1])*tstep/2;
								}

							}
							
							// Reativa-se o flag de processamento
							processflag = 1;

							if (debugging) {
								fprintf(dbg,"Processando elemento INTEGRADOR em [%d][%d]\n", actualnode, resultnode);
								fprintf(dbg,"// Se o bloco pode ser processado \n");
								fprintf(dbg,"// (INTEGRATOR não precisa ter o nó de entrada liberado \n");
								fprintf(dbg,"// uma vez que lida com valores históricos) \n");
								fprintf(dbg,"if (tobeprocessed[actualnode][resultnode] > 0) { \n");
								fprintf(dbg,"// O nó resultado é liberado para processamento. \n");
								fprintf(dbg,"// Este nó resultado possui X processamentos a dever \n");
								fprintf(dbg,"tobeprocessed[resultnode][resultnode] = nlinks[resultnode]; \n");
								fprintf(dbg,"// O bloco foi processado. \n");
								fprintf(dbg,"tobeprocessed[actualnode][resultnode] = 0; \n");
							}
							if (debugging) {
								fprintf(dbg,">>> valuemtx:\n");
								for (idebug=0; idebug<nnodes; idebug++) {
									for (jdebug=0; jdebug<nnodes; jdebug++) {
										fprintf(dbg," %09.4lf ", valuemtx[idebug][jdebug]);
										if (jdebug == nnodes - 1) {
											fprintf(dbg,"\n");
										}
									}
								}
							}

							if (debugging) {
								fprintf(dbg,">>> tobeprocessed:\n");
								for (idebug=0; idebug<nnodes; idebug++) {
									for (jdebug=0; jdebug<nnodes; jdebug++) {
										fprintf(dbg," %+d ", tobeprocessed[idebug][jdebug]);
										if (jdebug == nnodes - 1) {
											fprintf(dbg,"\n");
										}
									}
								}
							}

						}
					} // INTEGRATOR
				}
			}
		}

		// Buffer circular
		// Atualiza valores históricos
		for (i=0; i<nnodes; i++) {
			nodevalue[i][stepbuffer] = valuemtx[i][i];
			nodevalue[i][TAMBUFFER+stepbuffer] = valuemtx[i][i];
		}

		// Escrevendo o arquivo de resultado (.M)
		fprintf(resfile, "%f ", t);
		for (i=0; i<nnodes; i++) {
			fprintf(resfile,"%.15lf ", nodevalue[i][stepbuffer]);
		}
		fprintf(resfile, "\n");

		// Incrementando passo de integração
		t = t + tstep;
		stepbuffer++;

		// Reset do valor de stepbuffer
		// Tamanho do buffer é limitado
		if (stepbuffer >= TAMBUFFER) {
			stepbuffer = 0;
		}

		// limpa tobeprocessed e valores nos nós (valuemtx) e atualiza
		// a parte de somadores e multiplicadores quanto ao número de
		// fatores (processos) a serem calculados;
		for (i=0; i<MAXNODES; i++) {
			for (j=0; j<MAXNODES; j++) {
				tobeprocessed[i][j] = -1;
				if (i==j) {
					valuemtx[i][j] = 0;
					if (nfactors[j] < -1) {
						tobeprocessed[i][j] = nfactors[j];
					}
				}
				if (typeelemmtx[i][j] == TYPE_DELAY) {
					tobeprocessed[i][j] = 1;
				}
			}
		}
		if (debugging) {
			debugging++;
		}
		if (debugging == DEBUGCYCLES+1) {
			debugging = 0;
		}

	}

	fprintf(resfile,"];");

	// Fecha arquivos abertos
	if (arq != NULL) fclose(arq);
	if (resfile!=NULL) fclose(resfile);
	if (dbg!=NULL) fclose(dbg);

	
	//getch();
}	

