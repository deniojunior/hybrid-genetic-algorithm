#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>

#define n 42
#define tamPopulacao 2000 // ou 100
#define nomeInstancia "/TSP/swiss42.tsp"
#define reconhecimentoGrafo "/VMP/swiss42_reconhecimento.tsp"
#define nomeInstanciaResult "/RESULTS/swiss42_result.tsp"

#define quantSelecaoMelhores (int)(tamPopulacao*0.05)/100 //0.05% da população
#define quantSelecaoPais (int)(((tamPopulacao*40)/100) - quantSelecaoMelhores) //39.95% da população
#define quantSelecaoMutacao (int)(tamPopulacao*10)/100 //10% da população selecionado aleatoriamente
#define quantNovosIndividuosAleatorios (int)(tamPopulacao*10)/100 //10% de novos individuos aleatórios serão criados
#define quantOffsprings (int)(tamPopulacao*40)/100 //40% da população
#define fatorMutacao (int)trunc(((n+1)*5)/100)+1 //10% dos genes de um indivíduo sofrerao mutação
#define fatorConvergencia 500 //Gerações sem avanço na evolução da elite da população
#define locus 38 //Locus do cromossomo
#define quantExtincaoEmMassa 50

enum boolean{
    true = 1, false = 0
};

typedef enum boolean bool;

struct Define_Individuo
  {
    int individuo[n+1];
  };

void getch(void)
{
   system("read key");
}

void ordenaQuickSort(int vetor[tamPopulacao][2], int inicio, int fim){

    int pivo, aux, i, j, meio;

    i = inicio;
    j = fim;

    meio = (int) ((i + j) / 2);
    pivo = vetor[meio][0];

    do{
        while(vetor[i][0] < pivo){
            i = i + 1;
        }

        while(vetor[j][0] > pivo){
            j = j - 1;
        }

        if(i <= j){
            aux = vetor[i][0];
            vetor[i][0] = vetor[j][0];
            vetor[j][0] = aux;

            aux = vetor[i][1];
            vetor[i][1] = vetor[j][1];
            vetor[j][1] = aux;

            i = i + 1;
            j = j - 1;
        }

    }while(j > i);

    if(inicio < j){
        ordenaQuickSort(vetor, inicio, j);
    }

    if(i < fim){
        ordenaQuickSort(vetor, i, fim);
    }
}

void embaralhar(int *vet, int k){
    srand(k+(time(NULL)));//time(NULL)); //Gerar novos numeros aleatorios a cada instancia

    int i = 0, r = 0, temp = 0;

    for(i = 0; i < n; i++){

        r = rand() % n;

        temp = vet[i];
        vet[i] = vet[r];
        vet[r] = temp;
    }
}

void embaralharSemLocus(int *vet, int k, int controlCidadesPossiveis){
    srand(k+(time(NULL)));//time(NULL)); //Gerar novos numeros aleatorios a cada instancia

    int i = 0, r = 0, temp = 0;

    for(i = 0; i < controlCidadesPossiveis; i++){

        r = rand() % controlCidadesPossiveis;

        temp = vet[i];
        vet[i] = vet[r];
        vet[r] = temp;
    }
}

void embaralharGenes(int *genes, int k){
    srand(k+(time(NULL)));//time(NULL)); //Gerar novos numeros aleatorios a cada instancia

    int i = 0, r = 0, temp = 0;

    for(i = 0; i < (n-1); i++){

        r = rand() % (n-1);

        temp = genes[i];
        genes[i] = genes[r];
        genes[r] = temp;
    }
}

void trocaGenes(int *vet, int *random){

    int i = 0, r = 0, temp = 0;

    for(i = 0; i < fatorMutacao; i++){

        temp = vet[random[r]];
        vet[random[r]] = vet[random[r+1]];
        vet[random[r+1]] = temp;
        r++;

    }
}

void leMatriz(int matAdjacencia[n][n]){

    FILE *arquivo;
    char cwd[1024];
    int i, j;

    getcwd(cwd, sizeof(cwd));
    strcat(cwd, nomeInstancia);

    arquivo = fopen(cwd, "r");
    //printf("\n%s", cwd);

    if(arquivo == NULL){
        printf("\nErro ao abrir o arquivo!\n");
    }

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            fscanf(arquivo, "%d", &matAdjacencia[i][j]);
        }
    }

    fclose(arquivo);
}

void leReconhecimento(struct Define_Individuo * populacaoVMP){

    FILE *arquivo;
    char cwd[1024];
    int i, j;

    getcwd(cwd, sizeof(cwd));
    strcat(cwd, reconhecimentoGrafo);

    arquivo = fopen(cwd, "r");
    //printf("\n%s", cwd);

    if(arquivo == NULL){
        printf("\nErro ao abrir o arquivo!\n");
    }

    for(i = 0; i < n; i++){
        //printf("\nIndividuo %d: ", i);
        for(j = 0; j < (n+1); j++){
            fscanf(arquivo, "%d", &populacaoVMP[i].individuo[j]);
          //  printf("%d ", populacao[i].individuo[j]);
        }
    }
   // printf("\n");

    /*for(i = n; i < tamPopulacao; i++){
        for(j = 0; j < (n+1); j++){
            populacao[i].individuo[j] = populacao[i-n].individuo[j];
        }
    }*/

    fclose(arquivo);
}

int insereVirusVMP(struct Define_Individuo * populacao, struct Define_Individuo * populacaoVMP, int i){
    int j, indiceVirus = 0, controlReconhecimento = 0;

    srand(((i+11*time(NULL))*(j+19*time(NULL)))+(time(NULL)));
    indiceVirus = rand() % ((n-2)-locus)+1; //n-2, pois a ultima não entra; e -locus para evitar que o indiceVirus ultrapasse n
                                            //+1 porque a primeira posição não entra

    srand(((i+7*time(NULL))*(j+13*time(NULL)))+(time(NULL)));
    controlReconhecimento = rand() % n;

    for(j = 0; j < locus; j++){
        populacao[i].individuo[indiceVirus+j] = populacaoVMP[controlReconhecimento].individuo[indiceVirus+j];
    }

    return indiceVirus;
}

void gerarPopulacaoInicial(struct Define_Individuo * populacao, struct Define_Individuo * populacaoVMP){

    int i, j, k, indiceVirus = 0;
    int * cidadesPossiveis = calloc((n-locus), sizeof(int));
    bool contemCidade = false;
    int controlCidadesPossiveis = 0;

    for(i = 0; i < tamPopulacao; i++){//tamPopulacao; i++){

        indiceVirus = insereVirusVMP(populacao, populacaoVMP, i);
        controlCidadesPossiveis = 0;
        /*printf("\nindiceVirus = %d\n", indiceVirus);
        printf("\nk = %d\n", indiceVirus);
        printf("\nk+indiceVirus = %d\n",(indiceVirus+(locus-1)));*/

        for(j = 0; j < n; j++){

            contemCidade = false;
            //Pequeno for

            for(k = indiceVirus; k < indiceVirus+locus; k++){
                if(populacao[i].individuo[k] == j){
                    contemCidade = true;
                }
            }

            if(!contemCidade){
                cidadesPossiveis[controlCidadesPossiveis] = j;
                controlCidadesPossiveis++;
            }
        }

        /*printf("\nCidades possíveis: ");
        for(k = 0; k < controlCidadesPossiveis; k++){
            printf("%d ", cidadesPossiveis[k]);
        }
        printf("\n");*/

        embaralharSemLocus(cidadesPossiveis, i*tamPopulacao, controlCidadesPossiveis); //2º parametro é para aumentar o fator de aleatoriedade

        /*printf("\nCidades possíveis: ");
        for(k = 0; k < controlCidadesPossiveis; k++){
            printf("%d ", cidadesPossiveis[k]);
        }
        printf("\n");*/

        controlCidadesPossiveis = 0;

        for(j = 0; j < n; j++){

            if((j < indiceVirus)||(j > indiceVirus+(locus-1))){

                populacao[i].individuo[j] = cidadesPossiveis[controlCidadesPossiveis];
                controlCidadesPossiveis++;
            }
        }
    }


    //A ultima posição recebe a primeira(Rota retorna a cidade origem)
    for(i = 0; i < tamPopulacao; i++){
        populacao[i].individuo[n] = populacao[i].individuo[0];
    }
}

void calculaFitness(struct Define_Individuo * populacao, int matAdjacencia[n][n], int fitness[tamPopulacao][2]){

    int somaFitness;
    int linha, coluna;
    int i, j;

    somaFitness = 0;

    for(i = 0; i < tamPopulacao; i++){
        somaFitness = 0;

        for(j = 0; j < n; j++){
            linha = populacao[i].individuo[j+1];
            coluna = populacao[i].individuo[j];

            somaFitness = somaFitness + matAdjacencia[linha][coluna];
        }

        fitness[i][0] = somaFitness;
        fitness[i][1] = i;
    }
}

void selecao(int melhores[quantSelecaoMelhores], int resto[tamPopulacao - quantSelecaoMelhores], int pais[quantSelecaoPais][2],
        int popselecaoMutacao[tamPopulacao], int selecaoMutacao[quantSelecaoMutacao], int fitness[tamPopulacao][2], int g){

    int i, j;

    ordenaQuickSort(fitness, 0, tamPopulacao-1);

    //2% melhores da população
    for(i = 0; i < quantSelecaoMelhores; i++){
        melhores[i] = fitness[i][1];
    }

    //Separa o resto (98% da população) para que depois seja possível pegar os 38% para o pais
    for(i = quantSelecaoMelhores; i < tamPopulacao; i++){
        resto[i-quantSelecaoMelhores] = fitness[i][1];
    }

    embaralhar(resto, quantSelecaoMelhores*g); //2º parametro e para aumentar o fator de aleatoriedade

    //Seleciona 38% aleatoriamente para o pais
    for(i = 0; i < quantSelecaoPais; i++){
        for(j = 0; j < tamPopulacao; j++){

            if(resto[i] == fitness[j][1]){
                pais[i][0] = fitness[j][0];
                pais[i][1] = resto[i];
                break;
            }

        }
    }

    //ordenaQuickSort(pais, 0, quantSelecaoPais-1);

    //auxilia a seleção aleatória da população
    for(i = 0; i < tamPopulacao; i++){
        popselecaoMutacao[i] = i;
    }

    embaralhar(popselecaoMutacao, quantSelecaoMutacao*g); //2º parametro é para aumentar o fator de aleatoriedade

    //Seleciona os 10% aleatórios do restante para aplicar a mutação
    for(i = 0; i < quantSelecaoMutacao; i++){
        selecaoMutacao[i] = popselecaoMutacao[i];
    }

    //Totaliza acima 50%
    // Gerar 40% por crossover
    // Gerar 10% novos individuos aleatorios
}

void opCrossover(struct Define_Individuo * populacao, int melhores[quantSelecaoMelhores], int pais[quantSelecaoPais][2],
    struct Define_Individuo * crossover, struct Define_Individuo * offsprings, int g){

         int pos, cidadeAtual, controlCyrcle, i, j, k, l;
         int cyrcle[n];
         bool contain;

    for(i = 0; i < quantSelecaoMelhores; i++){
        for(j = 0; j < (n+1); j++){
            crossover[i].individuo[j] = populacao[melhores[i]].individuo[j];
        }
    }

    for(i = quantSelecaoMelhores; i < (quantSelecaoPais+quantSelecaoMelhores); i++){
        for(j = 0; j < (n+1); j++){
            crossover[i].individuo[j] = populacao[pais[i-quantSelecaoMelhores][1]].individuo[j];
        }
    }

    //inicializa variáveis
    pos = 0;
    cidadeAtual = -1;
    controlCyrcle = 0;

    k = 0;
    l = 0;

    i = 0;
    j = 1;

    while(j < quantOffsprings){
        controlCyrcle = 0;
        cidadeAtual = -1;

        //Gene n+1 ficam fora do crossover
        while(cidadeAtual != crossover[i].individuo[0]){

            cyrcle[controlCyrcle] = pos;
            cidadeAtual = crossover[j].individuo[pos];

            for(k = 0; k < n; k++){
                if(crossover[i].individuo[k] == cidadeAtual){
                    pos = k;
                    break;
                }
            }
            controlCyrcle++;
        }

        for(k = 0; k < controlCyrcle; k++){
            offsprings[i].individuo[cyrcle[k]] = crossover[i].individuo[cyrcle[k]];
            offsprings[j].individuo[cyrcle[k]] = crossover[j].individuo[cyrcle[k]];
        }

        if((controlCyrcle == 1)||(controlCyrcle == n)){

            //Se não houver ciclo -> Gerar dois novos indivíduos para evitar que a população se estaguine
            for(k = 0; k < n; k++){
                offsprings[i].individuo[k] = k;
                offsprings[j].individuo[k] = k;
            }

            embaralhar(offsprings[i].individuo, (i+7)*g); //2º parametro é para aumentar o fator de aleatoriedade
            embaralhar(offsprings[j].individuo, (j+7)*g); //2º parametro é para aumentar o fator de aleatoriedade

        }else{

            for(k = 0; k < n; k++){
                contain = false;
                for(l = 0; l < controlCyrcle; l++){
                    if(cyrcle[l] == k){
                        contain = true;
                        break;
                    }
                }

                if(!contain){
                    offsprings[i].individuo[k] = crossover[j].individuo[k];
                    offsprings[j].individuo[k] = crossover[i].individuo[k];
                }
            }
        }

        //Retorna à cidade inicial
        offsprings[i].individuo[n] = offsprings[i].individuo[0];
        offsprings[j].individuo[n] = offsprings[j].individuo[0];

        i += 2;
        j += 2;
    }
}

void opMutacao(struct Define_Individuo * populacao, int selecaoMutacao[quantSelecaoMutacao], struct Define_Individuo * mutacao, int g){

    int i;
    int j;
    int genes[n-1];

    for(i = 0; i < quantSelecaoMutacao; i++){
        for(j = 0; j < (n+1); j++){
            mutacao[i].individuo[j] = populacao[selecaoMutacao[i]].individuo[j];
        }
    }

    for(i = 1; i < n; i++){
        genes[i-1] = i;
    }

    for(i = 0; i < quantSelecaoMutacao; i++){
        embaralharGenes(genes, g*fatorMutacao);
        trocaGenes(mutacao[i].individuo, genes);
    }
}

void opNovosIndividuos(struct Define_Individuo * novosIndividuosAleatorios, int g){

    int i;
    int j;

    for(i = 0; i < quantNovosIndividuosAleatorios; i++){
       for(j = 0; j < n; j++){
            novosIndividuosAleatorios[i].individuo[j] = j;
        }

        embaralhar(novosIndividuosAleatorios[i].individuo, (i+7)*g); //2º parametro é para aumentar o fator de aleatoriedade
    }

   //A ultima posição recebe a primeira(Rota retorna a cidade origem)
    for(i = 0; i < quantNovosIndividuosAleatorios; i++){
        novosIndividuosAleatorios[i].individuo[n] = novosIndividuosAleatorios[i].individuo[0];
    }
}

void novaGeracao(struct Define_Individuo * populacao, struct Define_Individuo * crossover, struct Define_Individuo * offsprings,
    struct Define_Individuo * mutacao, struct Define_Individuo * novosIndividuosAleatorios){

    int i;
    int j;
    int k;

    k = 0;

    for(i = 0; i < quantOffsprings; i++){
        for(j = 0; j < n+1; j++){
            populacao[k].individuo[j] = crossover[i].individuo[j];
        }
        k++;
    }

    for(i = 0; i < quantOffsprings; i++){
        for(j = 0; j < n+1; j++){
            populacao[k].individuo[j] = offsprings[i].individuo[j];
        }
        k++;
    }

    for(i = 0; i < quantSelecaoMutacao; i++){
        for(j = 0; j < n+1; j++){
            populacao[k].individuo[j] = mutacao[i].individuo[j];
        }
        k++;
    }

    for(i = 0; i < quantNovosIndividuosAleatorios; i++){
        for(j = 0; j < n+1; j++){
            populacao[k].individuo[j] = novosIndividuosAleatorios[i].individuo[j];
        }
        k++;
    }
}

void imprimeResultado(struct Define_Individuo * populacao, int fitness[tamPopulacao][2], int tempo, int g){

        int i;
        FILE *arquivo;
        char cwd[1024];

        getcwd(cwd, sizeof(cwd));
        strcat(cwd, nomeInstanciaResult);

        //printf("\n %s", cwd);

        arquivo = fopen(cwd, "a");

        if(arquivo == NULL){
            arquivo = fopen(cwd, "w");

            if(arquivo == NULL){
                printf("\nErro ao abrir o arquivo para gravar os resultados!\n");
            }
        }

        //printf("\nExecutado em %d milisegundos!\n", tempo);

        printf("\nFitness: %d", fitness[0][0]);
        fprintf(arquivo, "Fitness: %d km", fitness[0][0]);

        float minutos = (float)((float)((float)tempo / 1000 ) / 60);
        printf("\nExecutado em %f minutos", minutos);
        fprintf(arquivo, "\nTempo: %f minutos", minutos);

        printf("\nRota %d: ", g);
        fprintf(arquivo, "\nRota: ");
        for(i = 0; i < (n+1); i++){
           printf("%d ", populacao[fitness[0][1]].individuo[i]);
           fprintf(arquivo, "%d ", populacao[fitness[0][1]].individuo[i]);
        }

        fprintf(arquivo, "\n------------------------------------------------------------------------------------------------------------------------");
        printf("\n");

        fclose(arquivo);
}

bool analiseConvergencia(int * fitnessMelhorIndividuo, int fitness[tamPopulacao][2], int * quantFitnessSemMelhora){
    //printf("\n");
    //printf("fitnessMelhorIndividuo: %i\n", *fitnessMelhorIndividuo);
    //printf("fitness[0][0]: %i\n", fitness[0][0]);
    //printf("quantFitnessSemMelhora: %i\n", *quantFitnessSemMelhora);

    if(*fitnessMelhorIndividuo == fitness[0][0]){
        *quantFitnessSemMelhora += 1;
    }else{
        *fitnessMelhorIndividuo = fitness[0][0];
        *quantFitnessSemMelhora = 0;
    }

    if(*quantFitnessSemMelhora == fatorConvergencia){
        return true;
    }else{
        return false;
    }

    getch();
}


void extincaoEmMassa(struct Define_Individuo * populacao, int melhorIndividuo){
    struct Define_Individuo aux;
    int i, j;

    for(i = 0; i < (n+1); i++){
        aux.individuo[i] = populacao[melhorIndividuo].individuo[i];
    }

    for(i = 0; i < tamPopulacao; i++){
          for(j = 0; j < n; j++){
             populacao[i].individuo[j] = j;
          }
          embaralhar(populacao[i].individuo, i+tamPopulacao); //2º parametro é para aumentar o fator de aleatoriedade
      }

     //A ultima posição recebe a primeira(Rota retorna a cidade origem)
     for(i = 0; i < tamPopulacao; i++){
         populacao[i].individuo[n] = populacao[i].individuo[0];
     }

     for(i = 0; i < (n+1); i++){
        populacao[0].individuo[i] = aux.individuo[i];
    }
}

int main(){

    ///VARIÁVEIS PEGAR O TEMPO
    struct timeval inicio, final;
    int tempo;

    ///VARIÁVEIS GERAIS
    int g; //GERAÇÕES
    struct Define_Individuo ind;

    ///VARIÁVEIS LER MATRIZ DE ADJACÊNCIA
    //FILE *arquivo;
    int matAdjacencia[n][n];

    ///VARIÁVEIS GERAR POPULAÇÃO INICIAL
    struct Define_Individuo * populacao = calloc(tamPopulacao, sizeof(ind));
    struct Define_Individuo * populacaoVMP = calloc(n, sizeof(ind));

    ///VARIÁVEIS FITNESS
    int fitness[tamPopulacao][2];

    ///VARIÁVEIS SELEÇÂO
    int melhores[quantSelecaoMelhores];
    int resto[tamPopulacao - quantSelecaoMelhores];
    int pais[quantSelecaoPais][2];
    int popselecaoMutacao[tamPopulacao];
    int selecaoMutacao[quantSelecaoMutacao];

    ///VARIÁVEIS CROSSOVER
    struct Define_Individuo * crossover = calloc ((quantSelecaoPais+quantSelecaoMelhores), sizeof(ind));
    struct Define_Individuo * offsprings = calloc ((quantOffsprings), sizeof(ind));

    ///VARIÁVEIS MUTAÇÃO
    struct Define_Individuo * mutacao = calloc(quantSelecaoMutacao, sizeof(ind));

    ///VARIÁVEIS GERAR 10% DE NOVOS INDIVÍDUOS PARA A NOVA GERAÇÃO
    struct Define_Individuo * novosIndividuosAleatorios = calloc(quantNovosIndividuosAleatorios, sizeof(ind));

    ///CRITÉRIO DE PARADA
    int quantFitnessSemMelhora;
    int fitnessMelhorIndividuo;
    bool criterioParada;

/// ------------------------------------------ ALGORITMO GENÉTICO -------------------------------------------------
    leMatriz(matAdjacencia);

    gettimeofday(&inicio, NULL);

//    printf("\nquantSelecaoMelhores: %d", quantSelecaoMelhores);
//    printf("\nquantSelecaoPais: %d", quantSelecaoPais);
//    printf("\nquantSelecaoMutacao: %d", quantSelecaoMutacao);
//    printf("\nquantNovosIndividuosAleatorios: %d", quantNovosIndividuosAleatorios);
//    printf("\nquantOffsprings: %d", quantOffsprings);
//    printf("\nquantSelecaoMelhores: %d", quantSelecaoMelhores);

    leReconhecimento(populacaoVMP);
    gerarPopulacaoInicial(populacao, populacaoVMP);

    int i, j;

    /*for(i = 0; i < tamPopulacao; i++){
        printf("\nIndivíduo %d: ", i);
        for(j = 0; j < n+1; j++){
            printf("%d ", populacao[i].individuo[j]);
        }
    }

    printf("\n");*/



    ///LOOP
    //g = 0;
    //for(g = 0; g < quantGeracoes; g++){ ///g Gerações

    fitnessMelhorIndividuo = fitness[0][0];
    quantFitnessSemMelhora = 0;

    for(int controlExtincao = 0; controlExtincao < quantExtincaoEmMassa; controlExtincao++){
        while(!criterioParada){

            criterioParada = analiseConvergencia(&fitnessMelhorIndividuo, fitness, &quantFitnessSemMelhora);

            calculaFitness(populacao, matAdjacencia, fitness);

            selecao(melhores, resto, pais, popselecaoMutacao, selecaoMutacao, fitness, g);

            ///OPERADORES
            opCrossover(populacao, melhores, pais, crossover, offsprings, g);
            opMutacao(populacao, selecaoMutacao, mutacao, g);
            opNovosIndividuos(novosIndividuosAleatorios, g);

            novaGeracao(populacao, crossover, offsprings, mutacao, novosIndividuosAleatorios);

           /* printf("\nGerações sem melhora: %d\n", quantFitnessSemMelhora);
            printf("Melhor fitness: %d\n", fitness[0][0]);*/
        }

        criterioParada = false;
        extincaoEmMassa(populacao, fitness[0][1]);
        quantFitnessSemMelhora = 0;
        //printf("\n\nEXTINÇÃO EM MASSA\n\n");
    }

    calculaFitness(populacao, matAdjacencia, fitness);

    ordenaQuickSort(fitness, 0, tamPopulacao-1);

    /*for(i = 0; i < 2000; i++){
        printf("\nIndividuo %d: %d ", i, fitness[i][0]);
    }
    printf("\nMelhor Solução: %d ", fitness[0][0]);*/

    gettimeofday(&final, NULL);
    tempo = (int) (1000 * (final.tv_sec - inicio.tv_sec) + (final.tv_usec - inicio.tv_usec) / 1000);

    imprimeResultado(populacao, fitness, tempo, g);

    return(0);
}
