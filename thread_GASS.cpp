// Genetic Active Site Search - GASS
// Author: Sandro Carvalho Izidoro
// Last update: 28/10/2020 (Vinicius Paiva)

/*

Funcionamento dessa versão:
São executados 4 GAs (que podem ser paralelizados), onde cada execução corresponde a 1 dos quadrantes do espaço 3D da proteína.
Em cada execução, a cada intervalo de 'ngeracoes/10', os 10 primeiros individuos ranqueados vao para a população final.
Assim, ao final de cada execução de quadrantes, 1/4 da população final é construída.
Por fim, um quinto GA é executado sobre a população final.

A versão 3 realiza a operação de mutação somente dentro do quadrante em que está executando a atual chamada do GA.

*/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <thread>
#include <vector>
#include <mutex>

using std::ifstream;
using std::ofstream;
using namespace std;

//---> STRUCTURE ATOMO - GENE OF THE INDIVIDUAL <---//
struct Atomo
{
	char amino[5];
	char atomo[4];
	char cadeia;
	int atomo_ID;
	float x, y, z;
};

//---> STRUCTURE SUBSTITUICAO - CONSERVATIVE MUTATIONS <---//
struct Substituicao
{
	char lista[3][10];
};

//---> GLOBAL VARIABLES <---//

char teste[20][4] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

// Variables - AG
int tpopulacao, torneio, ngeracoes, ranking, templates; // size of population, , tournament, number of generations, size of ranking, number of templates.
float txmutacao, txcruzamento;							// rates of crossover and mutation.

// Variables - Active Site reference
char proteina_sitiof[10], cadeia_sitiof[1]; // Name of protein, chain.
int tam_sitiof;								// Size of active site.
char *testecadeia;							// chain.
int *id_residuo_sitiof;						// Size of Active Site, ID of atoms
int tam_dist;
int *posicaolista;
char **nome_aminoacido;
Atomo *sitio_ativo_fitness;
char nome_proteina[14];

Substituicao Matriz[567];
int linhaarquivoref;

char ecnumber[100], uniprot[100], resolution[100], nome_enzima[100];

float centroide[3];

float centroide_aux[3] = {-9.980, -15.233, 17.398}; //1euu
//float centroide_aux[3] = {16.065,10.299,42.423}; //3nos
//float centroide_aux[3] = {39.504,45.474,29.456}; //1hpl
//float centroide_aux[3] = {22.677,-13.246,15.607}; //5cnv
//float centroide_aux[3] = {27.358,0.424,15.948}; //1mhl
//float centroide_aux[3] = {41.959,24.060,12.098}; //1ah7
//float centroide_aux[3] = {47.225,35.245,25.108}; //1n20
//float centroide_aux[3] = {49.699,36.726,23.317}; //1su4
//float centroide_aux[3] = {47.495, 21.292, 40.862}; //1a8h
//float centroide_aux[3] = {30.364 46.828 56.363}; //2frv
//float centroide_aux[3] = {2.151, 12.796, 10.278}; //1lba

//---> FUNCTIONS <---//
#include "fontes/openfile.hpp"
#include "fontes/setup.hpp"
#include "fontes/escreve.hpp"
#include "fontes/iniciaRef.hpp"
#include "fontes/repositorio.hpp"
#include "fontes/poslista.hpp"
#include "fontes/IniciaPopulacao.hpp"
#include "fontes/torneio.hpp"
#include "fontes/crossover.hpp"
#include "fontes/mutation.hpp"
#include "fontes/substituicao.hpp"
#include "fontes/fitref.hpp"
#include "fontes/fitness.hpp"
#include "fontes/best.hpp"
#include "fontes/quick.hpp"
#include "fontes/iniciaTop.hpp"
#include "fontes/crialistasub.hpp"
#include "atomic"
atomic<int> conta;
mutex mut;
void RunGASS(int quadrante, Atomo **pop_final, ofstream &ofs, Atomo *refe, long nrec, int *numero_aminos, Atomo **repositorio, int tam, char** substitutos, int inicio = 0);

int main(int argc, char *argv[])
{
	
	if (argc != 6 && argc != 3)
	{
		cout << "Error! Parameters missing!" << endl;
		cout << argc << endl;
		exit(1);
	}

	if (argc == 3)
	{
		for (int i = 0; i < 3; i++)
		{
			centroide[i] = centroide_aux[i];
		}
	}
	else
	{
		char coord[100];
		for (int i = 0; i < 3; i++)
		{
			strcpy(coord, argv[i + 3]);
			centroide[i] = stof(coord);
			//cout << centroide[i] << " ";
		}
	}
	srand(time(NULL));
	ofstream ofs;

	char caminho[200];
	strcpy(caminho, argv[1]);
	strcat(caminho, "/ActiveSitesFound.txt");
	ofs.open(caminho); // Open file to save the active sites found
	//ofs.open("ActiveSitesFound.txt"); // Open file to save the active sites found
	if (ofs.is_open() == false)
	{
		cout << "File of active sites found can not be open.\n";
		exit(1);
	}
	linhaarquivoref = 0;
	ConfiguracaoGA(argv[2]);
	//cout << tpopulacao << " " << tam_sitiof << " " << ranking << endl;

	for (int ntemplates = 0; ntemplates < templates; ntemplates++)
	{
		conta = 0;
		Atomo **pop, **pop_final;
		ConfiguracaoSitioAtivo(argv[2]); //Active Site Reference - Configuration
		Atomo *refe = new Atomo[tam_sitiof];
		IniciaReferencia(refe, argv[2]); //
										 //cout << "Sitio de referencia: " << endl;

		pop_final = new Atomo *[tpopulacao];
		for (int i = 0; i < tpopulacao; i++)
			pop_final[i] = new Atomo[tam_sitiof];

		//top = new Atomo*[ranking];
		//for(int i = 0; i < ranking; i++)
		//top[i] = new Atomo[tam_sitiof];

		strcpy(nome_proteina, argv[1]);
		strcat(nome_proteina, "/protein.dat\0");

		ifstream linhasFile(nome_proteina, ios::in); //Open file from filtroC_Dados.cpp
		if (!linhasFile)
		{
			cout << "End of Program or Fault on File of Proteins!" << endl;
			exit(1);
		}
		linhasFile.seekg(0, ios::end);					  //put pointer in the end of file
		long nrec = (linhasFile.tellg()) / sizeof(Atomo); //return number of register
		linhasFile.close();
		Atomo ***repositorio = new Atomo **[5];
		for (int j = 0; j < 5; j++)
		{
			repositorio[j] = new Atomo *[20];
			for (int i = 0; i < 20; i++)
				repositorio[j][i] = new Atomo[nrec];
		}

		int **numero_aminos = new int *[5];
		for (int j = 0; j < 5; j++)
		{
			numero_aminos[j] = new int[20];
			for (int i = 0; i < 20; i++)
				numero_aminos[j][i] = 0;
		}
		GeraRepositorio(repositorio, numero_aminos);
		// Substitution matrix
		CriaListaSubstituicao(Matriz, argv[2]);
		char** substitutos = new char*[30];
		for(int i=0;i<30;i++){
			substitutos[i] = new char[4];
		}
		int tam = ListaSubstituicao(substitutos);
		vector<thread> threads;
		for (int i = 1; i < 5; i++)
		{
			threads.push_back(thread([](int i,int tam,char** substitutos,long nrec,  Atomo **pop_final, ofstream &ofs, Atomo *refe, int *numero_aminos, Atomo **repositorio) {
				RunGASS(i, pop_final, ofs, refe, nrec, numero_aminos, repositorio, tam,substitutos, 0);
			},
									 i, tam,substitutos,nrec, pop_final, ref(ofs), refe, numero_aminos[i], repositorio[i]));

			//cout << "------------------------------------" << endl;
		}
		for (auto &x : threads)
		{
			x.join();
		}

		RunGASS(0, pop_final, ref(ofs), refe, nrec, numero_aminos[0], repositorio[0],tam,substitutos, conta);
		// Deleting pointers
		delete[] testecadeia;
		delete[] id_residuo_sitiof;
		delete[] refe;
		for (int i = 0; i < tam_sitiof; i++)
			delete nome_aminoacido[i];
		for (int i = 0; i < tpopulacao; i++)
			delete pop_final[i];
	}
	cout << "GASS finished! Please, verify the file ActivesSitesFound.txt \n";
}

void RunGASS(int quadrante, Atomo **pop_final, ofstream &ofs, Atomo *refe, long nrec, int *numero_aminos, Atomo **repositorio, int tam, char** substitutos, int inicio)
{
	int contador = 0, pos;
	char buff[20]; // line from file
	Atomo **pop;
	if (quadrante != 0)
	{
		pop = new Atomo *[tpopulacao];
		for (int i = 0; i < tpopulacao; i++)
			pop[i] = new Atomo[tam_sitiof];
	}
	else
	{
		pop = pop_final;
	}

	tam_dist = ((tam_sitiof * tam_sitiof) - tam_sitiof) / 2;
	float fitr, fitnessvet[tpopulacao], distanciasref[tam_dist], distanciasitio[tam_dist];
	float fitnessvetfinal[ngeracoes], fitnesspior[ngeracoes], fitnessmedio[ngeracoes];

	// Starting repository enzyme

	// Starting ranking
	float *fittop = new float[ranking];
	IniciaTop(fittop);

	// Number of each residues int the protein

	// Starting top
	Atomo **top = new Atomo *[ranking];
	for (int i = 0; i < ranking; i++)
		top[i] = new Atomo[tam_sitiof];

	contador++;

	// Starting population
	/*
	   if(pop == NULL){
	   		pop = new Atomo*[tpopulacao];
	   		for(int i = 0; i < tpopulacao; i++)
		  		pop[i] = new Atomo[tam_sitiof];
		}
	   */
	Atomo **ms = new Atomo *[ngeracoes];
	for (int i = 0; i < ngeracoes; i++)
		ms[i] = new Atomo[tam_sitiof];

	//cout << nome_enzima << endl <<  ecnumber << endl << uniprot << endl <<  resolution << endl;

	int ResiduosNoRepositorioMenorQueNoSitio = 1;
	int ResiduosNoRepositorioZerado = 1;

	// Checks if there are sufficient number of residues in the repository
	for (int i = 0; i < tam_sitiof; i++)
	{
		int k = 0;
		for (int j = 0; j < tam_sitiof; j++)
			if (!strcmp(nome_aminoacido[i], nome_aminoacido[j]))
				k++;
		int posicao = posicaonalista(teste, nome_aminoacido[i]);
		if (numero_aminos[posicao] < k)
		{
			ResiduosNoRepositorioMenorQueNoSitio = 0;
			i = tam_sitiof;
		}
	}

	// If there is a substitution matrix, check if there are residues in repository
	if (tam > 0)
	{
		for (int i = 0; i < tam; i++)
		{
			int posicao = posicaonalista(teste, substitutos[i]);
			if (numero_aminos[posicao] == 0)
			{
				ResiduosNoRepositorioZerado = 0;
				i = tam;
			}
		}
		if (ResiduosNoRepositorioZerado)
			for (int i = 0; i < tam_sitiof; i++)
			{
				int posicao = posicaonalista(teste, nome_aminoacido[i]);
				if (numero_aminos[posicao] == 0)
				{
					ResiduosNoRepositorioZerado = 0;
					i = tam_sitiof;
				}
			}
	}

	// Sem lista de substituicao e Residuos em quantidade igual a referencia OU
	// Lista de substituicao e ao menos um residuo de cada tipo no repositorio
	if (((tam == 0) && (ResiduosNoRepositorioMenorQueNoSitio == 1)) || ((tam > 0) && (ResiduosNoRepositorioZerado == 1)))
	{
		//cout << "vez +++ "  << endl;
		// Starting population

		IniciaPopulacao(pop, repositorio, numero_aminos, refe, inicio);

		// Fitness of reference
		fitr = fitnessRef(refe, distanciasref);
		float fitnessmelhores[ngeracoes][2];
		int contadorfitness = 0;
		float fitnessparamudar = -100;
		int contador10 = 0;
		int pedaco = 0;

		// ----------------------------------------------------------------------------------------------------------------------------------------------------

		//GA - Fitness, crossover, mutation

		for (int tempo = 0; tempo < ngeracoes; tempo++)
		{
			// Fitness function
			fitness(pop, fitnessvet, distanciasitio, distanciasref, tam_dist);

			//Quicksort
			quickSort(pop, fitnessvet, 0, tpopulacao - 1);

			// Top ranking
			int containdividuos = 0;
			int conttop = 0;
			while ((containdividuos < tpopulacao) && (fitnessvet[containdividuos] < fittop[ranking - 1]))
			{
				if (fitnessvet[containdividuos] < fittop[conttop])
				{
					for (int i = ranking - 1; i > conttop; i--)
					{
						fittop[i] = fittop[i - 1];
						for (int j = 0; j < tam_sitiof; j++)
							top[i][j] = top[i - 1][j];
					}
					fittop[conttop] = fitnessvet[containdividuos];
					for (int j = 0; j < tam_sitiof; j++)
						top[conttop][j] = pop[containdividuos][j];
					containdividuos++;
				}
				else if (fitnessvet[containdividuos] == fittop[conttop])
					containdividuos++;
				else if (fitnessvet[containdividuos] > fittop[conttop])
					conttop++;
			}

			// Best of the generation
			pos = melhordageracao(fitnessvet);
			for (int j = 0; j < tam_sitiof; j++)
				ms[tempo][j] = pop[pos][j];

			// Final fitness
			fitnessvetfinal[tempo] = fitnessvet[pos];

			// Worst fitness
			pos = piordageracao(fitnessvet);
			fitnesspior[tempo] = fitnessvet[pos];

			// Average fitness
			float g = 0;
			for (int j = 0; j < tpopulacao; j++)
				g = g + fitnessvet[j];
			fitnessmedio[tempo] = g / tpopulacao;

			if (fitnessparamudar == fittop[0])
				contadorfitness++;
			else
			{
				fitnessparamudar = fittop[0];
				contadorfitness = 0;
			}

			// Swap the better if it does not change for 5 generations

			if (contadorfitness == 5)
			{
				char find[5];
				int posicao, j, i = 0, conta = 0;
				int tentativas = 0;
				int quadrante_AA;
				while (fittop[0] == fitnessvet[i])
				{
					for (j = 0; j < tam_sitiof; j++)
					{
						strcpy(find, pop[i][j].amino);
						pos = posicaonalista(teste, find);
						while (1)
						{
							pop[i][j] = repositorio[pos][rand() % numero_aminos[pos]];
							quadrante_AA = MostraQuadrante(pop[i][j], centroide);
							if (quadrante_AA == quadrante || quadrante == -1)
								break;
							else
							{
								tentativas++;
								if (tentativas > 1000)
									break;
							}
						}
					}
					i++;
				}
				contadorfitness = 0;
			}

			// Crossover Operator
			crossover(pop, fitnessvet);

			// Mutation Operator
			mutation(pop, numero_aminos, repositorio);

			// Runs substitution if there are residues in substitution matrix
			if ((tam > 0) && (ResiduosNoRepositorioZerado == 1))
				substituicao(pop, numero_aminos, repositorio, substitutos, tam, refe);

			// Fitness of population
			fitness(pop, fitnessvet, distanciasitio, distanciasref, tam_dist);

			//Quicksort
			quickSort(pop, fitnessvet, 0, tpopulacao - 1);

			pos = melhordageracao(fitnessvet);

			//Elitism
			for (int j = 0; j < tam_sitiof; j++)
				pop[pos][j] = ms[tempo][j];

			if (quadrante != 0 && (tempo + 1) % 20 == 0)
			{

				//cout << "Vez " << tempo << endl;
				/*for(int i=0; i<10; i++){
					for(int j=0; j<tam_sitiof; j++){
						cout <<top[i][j].amino<<" "<< top[i][j].atomo_ID <<" "<<top[i][j].cadeia << " ";
					}
					cout << endl;
				}*/
				for (int i = 0; i < 10; i++)
				{
					mut.lock();
					for (int j = 0; j < tam_sitiof; j++)
					{
						pop_final[conta][j] = pop[i][j];
					}
					conta++;
					pedaco++;
					mut.unlock();
					//cout << pedaco << endl;
				}
			}

			/*
			if(quadrante == 1 && tempo == 199){
				cout << "pop_final###########: " << endl;
				for(int i=0; i<100;i++){
					for(int j = 0; j < tam_sitiof; j++)
						cout <<pop_final[i][j].amino<<" "<< pop_final[i][j].atomo_ID <<" "<<pop_final[i][j].cadeia << " ";
					cout << endl;
				}
				
				cout << "top#################: " << endl;
				for(int i=0; i<100;i++){
					for(int j = 0; j < tam_sitiof; j++)
						cout <<top[i][j].amino<<" "<< top[i][j].atomo_ID <<" "<<top[i][j].cadeia << " ";
					cout << endl;
				}
				
			}
			cout << "top#################: " << endl;
			if(quadrante == 4){
				for(int i=0; i<400;i++){
					for(int j = 0; j < tam_sitiof; j++)
						
						cout <<pop[i][j].amino<<" "<< pop[i][j].atomo_ID <<" "<<pop[i][j].cadeia << " Quadrante: " << MostraQuadrante(pop[i][j], centroide) << " | ";
					cout << endl;
				}
			}
			*/
		}

		// ----------------------------------------------------------------------------------------------------------------------------------------------------

		pos = melhordageracaof(fitnessvetfinal);

		//ofs << "Template: " << proteina_sitiof << endl;
		if (quadrante == 0)
		{
			for (int i = 0; i < ranking; i++)
			{
				ofs << contador << "\t" << i << "\t";
				ofs << setprecision(3) << fittop[i] << "\t";
				// ofs << nome_proteina << " " << fittop[i] << " ";
				if (fittop[i] != 1000)
				{
					ofs << top[i][0].amino << " " << top[i][0].atomo_ID << " " << top[i][0].cadeia;
					for (int j = 1; j < tam_sitiof; j++)
						ofs << ";" << top[i][j].amino << " " << top[i][j].atomo_ID << " " << top[i][j].cadeia;
					ofs << "\t" << nome_enzima << "\t";
					ofs << refe[0].amino << " " << refe[0].atomo_ID << " " << refe[0].cadeia;
					for (int j = 1; j < tam_sitiof; j++)
						ofs << ";" << refe[j].amino << " " << refe[j].atomo_ID << " " << refe[j].cadeia;
				}
				else
				{
					ofs << " - "
						<< " "
						<< " - "
						<< " "
						<< " - ";
					for (int j = 1; j < tam_sitiof; j++)
						ofs << ";"
							<< " - "
							<< " "
							<< " - "
							<< " "
							<< " - ";
					ofs << "\t" << nome_enzima << "\n";
					ofs << refe[0].amino << " " << refe[0].atomo_ID << " " << refe[0].cadeia;
					for (int j = 1; j < tam_sitiof; j++)
						ofs << ";" << refe[j].amino << " " << refe[j].atomo_ID << " " << refe[j].cadeia;
				}
				ofs << "\t" << ecnumber << "\t" << uniprot << "\t" << resolution << "\t" << tam_sitiof;
				ofs << endl;
			}
		}

		//if (quadrante == 0)
			//ofs << " -- TEMPLATE -- " << endl;
		//ofs << endl;
		// Deleting pointers
		delete[] fittop;

		for (int i = 0; i < ranking; i++)
			delete top[i];
		for (int i = 0; i < ngeracoes; i++)
			delete ms[i];
	}
}

//----> FUNCTIONS - IMPLEMENTATION <---- //
