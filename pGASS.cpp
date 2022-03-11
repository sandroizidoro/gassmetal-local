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

using std::ifstream;
using std::ofstream;
using namespace std;


//---> STRUCTURE ATOMO - GENE OF THE INDIVIDUAL <---//
struct Atomo{
	char amino[5];
	char atomo[4];
	char cadeia;
	int atomo_ID;
	float x,y,z;
};

//---> STRUCTURE SUBSTITUICAO - CONSERVATIVE MUTATIONS <---//
struct Substituicao{
	char lista[3][10];
};

//---> GLOBAL VARIABLES <---//

char teste[20][4] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };

// Variables - AG
int tpopulacao, torneio, ngeracoes, ranking, templates; // size of population, , tournament, number of generations, size of ranking, number of templates.
float txmutacao, txcruzamento; // rates of crossover and mutation.

// Variables - Active Site reference
char proteina_sitiof[10], cadeia_sitiof [1]; // Name of protein, chain.
int tam_sitiof; // Size of active site.
char *testecadeia; // chain.
int *id_residuo_sitiof; // Size of Active Site, ID of atoms
int tam_dist;
int *posicaolista;
char **nome_aminoacido;
Atomo *sitio_ativo_fitness;
char nome_proteina[14];

Substituicao Matriz [567];
int linhaarquivoref;

char ecnumber[100], uniprot[100], resolution[100], nome_enzima[100];

float centroide[3];

//float centroide_aux[3] = {-9.980,-15.233,17.398}; //1euu
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
//float centroide_aux[3] = {132.709,14.901,22.418}; //1f1d
float centroide_aux[3] = {159.937,92.626,51.358}; //2gvb


//---> FUNCTIONS <---//
void ConfiguracaoGA(char *argv); // GA Parameters
void ConfiguracaoSitioAtivo(char *argv); // Active Site Reference - Configuration
void IniciaReferencia(Atomo *refe, char *argv); // Setup reference active site
int  ListaSubstituicao(char substitutos [30][4]); // Builds substitution matrix
void CriaListaSubstituicao(Substituicao Matriz[567], char *argv);
void GeraRepositorio(Atomo **repositorio, int *numero_aminos); // Builds repository
void IniciaPopulacao(Atomo **pop,  Atomo **repositorio, int *numero_aminos, Atomo *refe, int quadrante); // Starting population
void substituicao (Atomo **pop, int *numero_aminos, Atomo **repositorio, char substitutos [30][4], int tam, Atomo *refe); // Substitution
void crossover(Atomo **pop, float fitnessvet[]); //Crossover operator
int funcaotorneio(float fitnessvet[], int escolhidos_torneio[]);
void mutation (Atomo **pop, int *numero_aminos, Atomo **repositorio, int quadrante); //Mutation operator
float fitnessRef (Atomo *refe, float distanciasref[]);
void fitness(Atomo **p, float fitnessv[], float distanciasitio[], float distanciasref[]);
int posicaonalista (char teste[20][4], char find[5]);
int melhordageracao(float fitnessvet[]);
int melhordageracaof(float fitnessvet[]);
int piordageracao(float fitnessvet[]);
void EscreveAtomo(Atomo atomo);
void quickSort(Atomo **tempop, float [], int, int);
int partition(Atomo **tempop, float [], int, int);
void IniciaTop(float fittop[]); //Starting top
bool AbrirArquivo(ifstream &ifs, char nomeArquivo[]); //Open file

void ImprimeRepositorio(Atomo **repositorio, int *numero_aminos);
int MostraQuadrante(Atomo amino, float *centroide);
void IniciaPopulacao2(Atomo **pop,  Atomo **repositorio, int *numero_aminos, Atomo *refe, int quadrante);
void RunGASS(char *argv[], Atomo **pop, int quadrante, Atomo **pop_final, ofstream &ofs, Atomo *refe);

int main(int argc, char *argv[]){
	if(argc != 6 && argc != 3){
		cout << "Error! Parameters missing!" << endl;
		cout << argc << endl;
		exit(1);
   	}
   	
   	if(argc == 3){
   		for(int i = 0; i < 3; i++){
   			centroide[i] = centroide_aux[i];
   		}
   	}
   	else{
   		char coord[100];
   		for(int i = 0; i < 3; i++){
   			strcpy(coord, argv[i+3]);
   			centroide[i] = stof(coord);
   			//cout << centroide[i] << " ";
   		}
   	}
   	srand (time (NULL));
   	ofstream ofs;

   	char caminho[200];
   	strcpy(caminho, argv[1]);
   	strcat (caminho, "/ActiveSitesFound.txt");
   	ofs.open(caminho); // Open file to save the active sites found
   	//ofs.open("ActiveSitesFound.txt"); // Open file to save the active sites found
   	if(ofs.is_open() == false){
	 	cout << "File of active sites found can not be open.\n";
	 	exit(1);
   	}
   	linhaarquivoref = 0;
	ConfiguracaoGA(argv[2]);
	//cout << tpopulacao << " " << tam_sitiof << " " << ranking << endl;
	
	for(int ntemplates = 0; ntemplates < templates; ntemplates++){
		Atomo **pop, **pop_final;
		ConfiguracaoSitioAtivo(argv[2]); //Active Site Reference - Configuration
		Atomo *refe = new Atomo[tam_sitiof];
	  	IniciaReferencia(refe, argv[2]); //
	  	//cout << "Sitio de referencia: " << endl;
		for(int i = 0; i < tam_sitiof; i++){
			int quadrante_AA = MostraQuadrante(refe[i], centroide);
			//cout << "Quadrante: " << quadrante_AA << " ";
			//EscreveAtomo(refe[i]);
			//cout << endl;
		}
		
		pop_final = new Atomo*[tpopulacao];
		for(int i = 0; i < tpopulacao; i++)
			pop_final[i] = new Atomo[tam_sitiof];
		
		//top = new Atomo*[ranking];
		//for(int i = 0; i < ranking; i++)
			//top[i] = new Atomo[tam_sitiof];
		pop = new Atomo*[tpopulacao];
			for(int i = 0; i < tpopulacao; i++)
				pop[i] = new Atomo[tam_sitiof];

		for(int i = 0; i < 4; i++){
			RunGASS(argv, pop, i+1, pop_final, ofs, refe);
			//cout << "------------------------------------" << endl;
		}

		//cout << " ULTIMA " << endl;
		RunGASS(argv, pop_final, -1, pop_final, ofs, refe);
		// Deleting pointers
	  	delete[] testecadeia;
	  	delete[]id_residuo_sitiof;
	  	delete[]refe;
	  	for(int i = 0; i < tam_sitiof; i++)
		 	delete nome_aminoacido[i];
		for(int i = 0; i < tpopulacao; i++)
		   delete pop[i];
		for(int i = 0; i < tpopulacao; i++)
		   delete pop_final[i];
		

	}
	cout << "GASS finished! Please, verify the file ActivesSitesFound.txt \n";
}


void RunGASS(char *argv[], Atomo **pop, int quadrante, Atomo **pop_final, ofstream &ofs, Atomo *refe){
  	int contador = 0, pos;
  	char buff[20];  // line from file
  	

  	tam_dist = ((tam_sitiof * tam_sitiof) - tam_sitiof)/2;
  	float fitr, fitnessvet[tpopulacao], distanciasref [tam_dist], distanciasitio[tam_dist];
  	float fitnessvetfinal[ngeracoes], fitnesspior[ngeracoes], fitnessmedio[ngeracoes];

  	// Substitution matrix
  	CriaListaSubstituicao(Matriz, argv[2]);
  	char substitutos [30][4];
  	int tam = ListaSubstituicao(substitutos);
	strcpy(nome_proteina,argv[1]);
	strcat(nome_proteina,"/protein.dat\0");

	   // Starting repository enzyme
	ifstream linhasFile(nome_proteina,ios::in); //Open file from filtroC_Dados.cpp
	if(!linhasFile){
 		cout << "End of Program or Fault on File of Proteins!" << endl;
 		exit(1);
	}
	linhasFile.seekg(0,ios::end); //put pointer in the end of file
	long nrec = (linhasFile.tellg( ))/sizeof(Atomo); //return number of register
	linhasFile.close();
	Atomo **repositorio;
	repositorio = new Atomo*[20];
	for(int i = 0; i < 20; i++)
  		repositorio[i] = new Atomo[nrec];

	// Starting ranking
	float *fittop = new float[ranking];
	IniciaTop(fittop);

	// Number of each residues int the protein
	int *numero_aminos = new int [20];
	for(int i = 0; i < 20; i++)
  		numero_aminos[i] = 0;

   	// Starting top
   	Atomo **top = new Atomo*[ranking];
   	for(int i = 0; i < ranking; i++)
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
   	Atomo **ms = new Atomo*[ngeracoes];
   	for(int i = 0; i < ngeracoes; i++)
		ms[i] = new Atomo[tam_sitiof];

   	GeraRepositorio(repositorio, numero_aminos); // Starting repository
   	
   	if(quadrante == 100)
   		ImprimeRepositorio(repositorio, numero_aminos);

   	//cout << nome_enzima << endl <<  ecnumber << endl << uniprot << endl <<  resolution << endl;

   	int ResiduosNoRepositorioMenorQueNoSitio = 1;
   	int ResiduosNoRepositorioZerado = 1;

   	// Checks if there are sufficient number of residues in the repository
   	for(int i = 0; i < tam_sitiof; i++){
		int k = 0;
	  	for(int j = 0; j < tam_sitiof; j++)
		 	if(!strcmp(nome_aminoacido[i], nome_aminoacido[j]))
		   	k++;
	  	int posicao = posicaonalista(teste, nome_aminoacido[i]);
	  	if(numero_aminos[posicao] < k){
		  	cout << "uno" << endl;
		  	ResiduosNoRepositorioMenorQueNoSitio = 0;
		  	i = tam_sitiof;
	  	}
	}

	// If there is a substitution matrix, check if there are residues in repository
	//cout << substitutos[1] << endl;
	//2jcw,HIS,CYS,THR,GLU,VAL,ASP
	if(tam > 0){
		//HIS CYS LEU MET THR GLU VAL ASP
		int verifica = 1;
		while(verifica == 1){
			verifica = 0;
			/*cout << "Vetor inicial: ";
			for(int i = 0; i < tam; i++)
				cout << substitutos[i] << " ";
			cout << tam << endl;*/
			for(int i = 0; i < tam; i++){
				int posicao = posicaonalista(teste, substitutos[i]);
			 	if(numero_aminos[posicao] == 0){
			 		verifica = 1;
			 		//cout << "Excluindo " << substitutos[i] << endl;
			   		if(i >= tam-1){
						tam = tam - 2;
					}
					else{
						if(i % 2 == 0){
							for(int j = i; j < tam; j=j+2){
								strcpy(substitutos[j],substitutos[j+2]);
								strcpy(substitutos[j+1],substitutos[j+3]);
							}
						}
						else{
							for(int j = i; j < tam; j=j+2){
								strcpy(substitutos[j],substitutos[j+2]);
								strcpy(substitutos[j-1],substitutos[j+1]);
							}
						}
						tam = tam - 2;
					}
				/*cout << "Novo vetor: ";
				for(int i = 0; i < tam; i++)
					cout << substitutos[i] << " ";
				cout << tam << endl;*/
			 	}
		  	}	
		}
	  	/*for(int i = 0; i < tam; i++)
			cout << substitutos[i] << " ";
		cout << tam << endl;
		cout << endl;*/
	  	if(ResiduosNoRepositorioZerado)
			for(int i = 0; i < tam_sitiof; i++){
		   		int posicao = posicaonalista(teste, nome_aminoacido[i]);
		   		if(numero_aminos[posicao] == 0){
		   			cout << "tres" << endl;
			 		ResiduosNoRepositorioZerado = 0;
			 	i = tam_sitiof;
		   		}
			}
	}

 	// Sem lista de substituicao e Residuos em quantidade igual a referencia OU
 	// Lista de substituicao e ao menos um residuo de cada tipo no repositorio
 	//cout << "resi: " << ResiduosNoRepositorioZerado << " tam: " << tam << endl;
 	if(((tam == 0) && (ResiduosNoRepositorioMenorQueNoSitio == 1)) || ((tam > 0) && (ResiduosNoRepositorioZerado == 1))){
 		//cout << "vez +++ "  << endl;
		// Starting population
		if(quadrante > 0 && quadrante < 5){
			IniciaPopulacao(pop, repositorio, numero_aminos, refe, quadrante);
		}
		

		// Fitness of reference
		fitr = fitnessRef (refe, distanciasref);
		float fitnessmelhores[ngeracoes][2];
		int contadorfitness = 0;
		float fitnessparamudar = -100;
		int contador10 = 0;
		int pedaco = 0;

		// ----------------------------------------------------------------------------------------------------------------------------------------------------

		//GA - Fitness, crossover, mutation
		
		for(int tempo=0; tempo < ngeracoes; tempo++){
			// Fitness function
		   	fitness(pop, fitnessvet, distanciasitio, distanciasref);

		   	//Quicksort
		   	quickSort(pop, fitnessvet, 0, tpopulacao-1);

		   	// Top ranking
		   	int containdividuos = 0;
		   	int conttop  = 0;
		   	while((containdividuos<tpopulacao) && (fitnessvet[containdividuos] < fittop[ranking-1])){
				if(fitnessvet[containdividuos] < fittop[conttop]){
				  	for(int i = ranking-1; i > conttop; i--){
					 	fittop[i] = fittop[i-1];
					 	for(int j = 0;j<tam_sitiof;j++)
							top[i][j] = top[i-1][j];
				  	}
				  	fittop[conttop] = fitnessvet[containdividuos];
				  	for(int j = 0;j<tam_sitiof;j++)
					 	top[conttop][j] = pop[containdividuos][j];
				  	containdividuos++;
				}
				else if(fitnessvet[containdividuos]==fittop[conttop])
						containdividuos++;
					 else if(fitnessvet[containdividuos]>fittop[conttop])
							conttop++;
		   	}

		   	// Best of the generation
		   	pos = melhordageracao(fitnessvet);
		   	for(int j=0 ; j < tam_sitiof; j++)
			  	ms[tempo][j] = pop[pos][j];

		   	// Final fitness
		   	fitnessvetfinal[tempo] = fitnessvet[pos];

		   	// Worst fitness
		   	pos = piordageracao(fitnessvet);
		   	fitnesspior[tempo] = fitnessvet[pos];

		   	// Average fitness
		   	float g = 0;
		   	for(int j=0 ; j < tpopulacao ; j++)
			  	g = g + fitnessvet[j];
		   	fitnessmedio[tempo] = g / tpopulacao;

		   	if(fitnessparamudar == fittop[0])
			 	contadorfitness++;
		   	else{
			   	fitnessparamudar = fittop[0];
			   	contadorfitness = 0;
		   	}

		   	// Swap the better if it does not change for 5 generations
		   	
		   	if(contadorfitness == 5) {
			 	char find[5];
			 	int posicao, j, i = 0, conta = 0;
			 	int tentativas = 0;
			 	int quadrante_AA;
			 	while(fittop[0] == fitnessvet[i]){
				  	for(j=0; j < tam_sitiof; j++){
					 	strcpy(find, pop[i][j].amino);
					 	pos = posicaonalista(teste, find);
					 	while(1){
							pop[i][j] = repositorio[pos][rand() % numero_aminos[pos]];
							quadrante_AA = MostraQuadrante(pop[i][j], centroide);
							if(quadrante_AA == quadrante || quadrante == -1)
								break;
							else{
			  					tentativas++;
			  					if(tentativas > 1000)
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
			mutation(pop, numero_aminos,repositorio, quadrante);

			// Runs substitution if there are residues in substitution matrix
			if((tam > 0) && (ResiduosNoRepositorioZerado == 1)){
			  	substituicao(pop, numero_aminos, repositorio, substitutos, tam, refe);
			}

			// Fitness of population
			fitness(pop, fitnessvet, distanciasitio, distanciasref);

			//Quicksort
			quickSort(pop, fitnessvet, 0, tpopulacao-1);

			pos = piordageracao(fitnessvet);

			//Elitism
			for(int j=0 ; j<tam_sitiof; j++)
			   pop[pos][j] = ms [tempo][j];

			if(quadrante != -1 && (tempo+1)%20 == 0){
				
				//cout << "Vez " << tempo << endl;
				/*for(int i=0; i<10; i++){
					for(int j=0; j<tam_sitiof; j++){
						cout <<top[i][j].amino<<" "<< top[i][j].atomo_ID <<" "<<top[i][j].cadeia << " ";
					}
					cout << endl;
				}*/
				for(int i = pedaco; i < (pedaco+10); i++){
					for(int j = 0; j < tam_sitiof; j++){
						if(fittop[(pedaco+i)%10] == 1000){
							pop_final[((quadrante-1)*(tpopulacao/4))+i][j] = pop[rand()%tpopulacao][j];
						}
						else{
							pop_final[((quadrante-1)*(tpopulacao/4))+i][j] = top[(pedaco+i)%10][j];
						}
						//if(j == 0)
							//cout << "pop[" << ((quadrante-1)*(tpopulacao/4))+i << "][" << j << "] top[" << (pedaco+i)%10 << "][" << j << "]" << endl;
					}
				}
				//cout << pedaco << endl;
				pedaco = pedaco + 10;

			}/*
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
		if(quadrante == -1){
			for(int i=0;i<ranking;i++){
				ofs << contador << "\t"<< i << "\t";
				ofs << setprecision(3) << fittop[i] << "\t";
				// ofs << nome_proteina << " " << fittop[i] << " ";
				if(fittop[i] != 1000){
				    ofs <<top[i][0].amino<<" "<< top[i][0].atomo_ID <<" "<<top[i][0].cadeia;
				    for(int j=1; j<tam_sitiof; j++)
					   ofs << ";" << top[i][j].amino<<" "<< top[i][j].atomo_ID <<" "<<top[i][j].cadeia;
				   	ofs << "\t" << nome_enzima << "\t";
					ofs << refe[0].amino<<" "<< refe[0].atomo_ID <<" "<<refe[0].cadeia;
				   	for(int j=1; j<tam_sitiof; j++)
					   	ofs << ";" <<refe[j].amino<<" "<< refe[j].atomo_ID <<" "<<refe[j].cadeia;
				}
				else{
					  	ofs << " - " << " " << " - " <<" "<<" - ";
					  	for(int j=1; j<tam_sitiof; j++)
						  	ofs << ";" << " - " << " " << " - " <<" "<<" - ";
					  	ofs << "\t" << nome_enzima << "\n";
			 			ofs <<refe[0].amino<<" "<< refe[0].atomo_ID <<" "<<refe[0].cadeia;
					  	for(int j=1; j<tam_sitiof; j++)
						  	ofs << ";" <<refe[j].amino<<" "<< refe[j].atomo_ID <<" "<<refe[j].cadeia;
				}
			   ofs << "\t" << ecnumber <<"\t"<<uniprot<<"\t"<< resolution << "\t" << tam_sitiof;
			   ofs << endl;
			}
		}
		//ofs << endl;
		// Deleting pointers
		delete[]fittop;
		for(int i = 0; i < 20; i++)
		   delete repositorio[i];
		//for(int i = 0; i < tpopulacao; i++)
		   //delete pop[i];
		delete[] numero_aminos;
		for(int i = 0; i < ranking; i++)
		   delete top[i];
		for(int i = 0; i < ngeracoes; i++)
		   delete ms[i];

    }
}




//----> FUNCTIONS - IMPLEMENTATION <---- //

// GA PARAMETERS
void ConfiguracaoGA(char *argv){
  ifstream fin;
  char buff[255];
  int nlinha = 0;
  char localtemplate[200];
  strcpy(localtemplate,argv);
  strcat(localtemplate,"/GA_Conf.txt");
  fin.open(localtemplate);
  //fin.open("GA_Conf.txt");
  if(fin.is_open()){
	 while(!fin.eof()){
		   fin.getline(buff, 255);
		   switch (nlinha) {
			   case 1: ngeracoes = atoi(buff); break;
			   case 3: tpopulacao = atoi(buff); break;
			   case 5: torneio = atoi(buff); break;
			   case 7: txcruzamento = atof(buff); break;
			   case 9: txmutacao = atof(buff); break;
			   case 11: ranking = atoi(buff); break;
			   case 13: templates = atoi(buff); break;
		   }
		   nlinha++;
	  }
	  fin.close();
  }
  else{
	   cout << "File GA_Conf.txt can not be open.\n";
	   exit(1);
  }
}


// CONFIGURATION ACTIVE SITE REFERENCE
void ConfiguracaoSitioAtivo(char *argv){
   ifstream fin;
   char buff[255];
   int nlinha = 0;
   char localtemplate[200];
   strcpy(localtemplate,argv);
   strcat(localtemplate,"/Templates.txt");
   // fin.open("Templates.txt");
   fin.open(localtemplate);
   if(fin.is_open()){
	 for(int i=0; i<linhaarquivoref;i++)
		fin.getline(buff, 25);
	 while(nlinha <= 4){
		  fin.getline(buff, 25);
		  switch (nlinha) {
			  case 1: strcpy(proteina_sitiof,buff); break;
			  case 3: tam_sitiof = atoi(buff); break;
		  }
		  nlinha++;
	 }
	 id_residuo_sitiof = new int[tam_sitiof];
	 nome_aminoacido = new char*[tam_sitiof];
	 testecadeia = new char[tam_sitiof];

	 for(int i = 0; i < tam_sitiof; i++)
		nome_aminoacido[i] = new char[4];

	 int i = 0, j = 0, k = 0, l = 0, linha = 1;

	 while((!fin.eof())&&(k<tam_sitiof*3)){
		  fin.getline(buff, 25);
		  if(linha == 1){
			 strcpy(nome_aminoacido[j], buff);
			 nome_aminoacido[j][3] = '\0';
			 j++;
			 linha++;
		  }
		  else if(linha == 2){
				 id_residuo_sitiof[i++] = atoi (buff);
				 linha++;
			   }
			   else{
				   strcpy(cadeia_sitiof,buff);
				   testecadeia[l++] = buff[0];
				   linha = 1;
			   }
		  nlinha++;
		  k++;
	 }
	 fin.close();
	 linhaarquivoref = linhaarquivoref + 5 + (tam_sitiof * 3);
   }
   else{
	   cout << "File templates.txt can not be open.\n";
	   exit(1);
   }
}

// SETUP REFERENCE
void IniciaReferencia(Atomo *refe, char *argv){
   Atomo aux;
   ifstream ifs;
   char localtemplate[50];
   char linha[100];
   strcpy(localtemplate,argv);
   strcat(localtemplate,"/");
   strcat(localtemplate, proteina_sitiof);
   //if(AbrirArquivo(ifs,proteina_sitiof) == true){
   if(AbrirArquivo(ifs,localtemplate) == true){
	 int i=0;
	 int contador = 0;
	 while((!ifs.eof()) && (i<tam_sitiof)){
		  if(contador<=3){
			ifs.read((char*)&linha, sizeof(linha));
			if(contador == 0) {strcpy(nome_enzima, linha);cout << nome_enzima << " : ";}
			   else if(contador == 1) strcpy(ecnumber, linha);
					else if(contador == 2) strcpy(uniprot, linha);
						 else if(contador == 3) strcpy(resolution, linha);
			contador++;
		  }
		  else{
		  ifs.read((char*)&aux, sizeof(Atomo));
		  for(int j=0;j<tam_sitiof;j++)
			 if((aux.atomo_ID == id_residuo_sitiof[j])&&(aux.cadeia==testecadeia[j])){
				refe[i++] = aux;
		        //EscreveAtomo(aux);
			 }
		  }

	 }
   }
   else{ cout << "File can not be open: " << proteina_sitiof << endl;}
   ifs.close();
   cout << endl;
   
}

// SUBSTITUTION MATRIX
int ListaSubstituicao(char substitutos [30][4]){
	int i = 0, ok = 0;
   	char proteina [5];
   	for(int j=0;j<4;j++){
	  	proteina[j] = (proteina_sitiof[j]);
	}
   	proteina[4] = '\0';
   	//cout << proteina << endl;
   	for(int k= 0 ; k < 567; k++){	
	   	//cout << Matriz[k].lista[0] << endl;
	   	if(!strcmp(proteina, Matriz[k].lista[0])){
		 	//cout << "aqui";
		 	strcpy(substitutos[i], Matriz[k].lista[1]);
		 	strcpy(substitutos[i+1], Matriz[k].lista[2]);
		 	//cout << substitutos[i] << " " << substitutos[i+1] << endl;
		 	i += 2;
		 	ok = 1;
	   	}
	   	else {
			if(ok == 1)
			   k = 567;
	   	}
   	}
   	return i;
}

// BUILDING REPOSITORY OF PROTEINS
void GeraRepositorio(Atomo **repositorio, int *numero_aminos){
	Atomo aux;
	ifstream inFile(nome_proteina,ios::binary);
	if(!inFile){
	  cout<<"File can not be open: " << nome_proteina << endl;
	   exit(1);
	}
	int contador = 0;
	char linha[100];
	while(inFile){
		if(contador<=3){
		   inFile.read((char*)&linha, sizeof(linha));
		//   if(contador == 0) strcpy(nome_enzima, linha);
		  //   else if(contador == 1) strcpy(ecnumber, linha);
			//        else if(contador == 2) strcpy(uniprot, linha);
			  //             else if(contador == 3) strcpy(resolution, linha);
		   contador++;
	  }
		else{
			inFile.read((char *)&aux, sizeof(Atomo));
			//EscreveAtomo(aux); int bunda; cin >> bunda;
		  for(int i = 0; i < 20; i++)
			   if(!strcmp(aux.amino,teste[i])){
				 repositorio[i][numero_aminos[i]] = aux;
				 numero_aminos[i]++;
				 /*while(1){
				 	repositorio[i][numero_aminos[i]] = aux;
				 	if(MostraQuadrante(repositorio[i][numero_aminos[i]], centroide) == quadrante)
				 		numero_aminos[i]++;
				 		break;
				 }*/
				 //EscreveAtomo(aux);
				 //cout << endl;
			   }
		}
   }
   inFile.close();

}

void ImprimeRepositorio(Atomo **repositorio, int *numero_aminos){
	//float centroide[3] = {-9.980,-15.233,17.398};
	for(int j = 0; j < 3; j++)
		cout << centroide[j] << endl;
	
	for(int i = 0; i < 20; i++){
		cout << repositorio[i][0].amino << ": " << numero_aminos[i] << endl;
		for(int j = 0; j < numero_aminos[i]; j++){
			cout << "quadrante: " << MostraQuadrante(repositorio[i][j], centroide) << " ";
			EscreveAtomo(repositorio[i][j]);
			cout << endl;
		}
	}   

}

int MostraQuadrante(Atomo amino, float *centroide){
	if(amino.x > centroide[0] && amino.y > centroide[1]){ //primeiro quadrante
		return 1;
	}
	else if(amino.x > centroide[0] && amino.y < centroide[1]){
		return 2;
	}
	else if(amino.x < centroide[0] && amino.y < centroide[1]){
		return 3;
	}
	else{
		return 4;
	}   
}

// STARTING POPULATION
void IniciaPopulacao(Atomo **pop,  Atomo **repositorio, int *numero_aminos, Atomo *refe, int quadrante){
	// IMPLEMENTACAO B
	int pos;
	int tentativas = 0;
	int total_quad[4][20] = {0};
	int quadrante_AA;
	int vez = 0;
	for(int i = 0; i < tpopulacao; i++){
	   for(int j=0; j<tam_sitiof; j++){
			if(j == 0){
				//cout << "Individuo " << i+j << endl;
			}
			// Calcula o quadrante ideal de cada individuo e todos os AA desse individuo devem pertencer a esse quadrante
			// Caso nao tenha uma AA nesse quadrante, escolhe outro aleatoriamente.
			pos = posicaonalista(teste, refe[j].amino);
			while(1){
				pop[i][j] = repositorio[pos][rand() % numero_aminos[pos]];
				quadrante_AA = MostraQuadrante(pop[i][j], centroide);
				if(quadrante_AA == quadrante){
					total_quad[quadrante_AA-1][pos]++;
					break;
				}
				else{
			  		tentativas++;
			  		if(tentativas > 1000){
			  			total_quad[quadrante_AA-1][pos]++;
			  			break;
			  		}
			  	}
			}
			//EscreveAtomo(pop[i][j]);
			//cout << "quadrante " << MostraQuadrante(pop[i][j], centroide) << endl;
			tentativas = 0;
	   }
	   //cout << "\n" << endl;
	}/*
	for(int i = 0; i < 20; i++)
		cout << teste[i] << " ";
	cout << endl;
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 20; j++){
			cout << total_quad[i][j] << "\t";
		}
		cout << endl;
	}*/
	
}

// STARTING POPULATION
void IniciaPopulacao2(Atomo **pop,  Atomo **repositorio, int *numero_aminos, Atomo *refe, int quadrante){
	int pos;
	int quadrante_AA;
	int total_quad[4][20] = {0};
	for(int i = 0; i < tpopulacao; i++){
	   for(int j=0; j<tam_sitiof; j++){
		  pos = posicaonalista(teste, refe[j].amino);
		  pop[i][j] = repositorio[pos][rand() % numero_aminos[pos]];
		  quadrante_AA = MostraQuadrante(pop[i][j], centroide);
		  EscreveAtomo(pop[i][j]);
		  cout << "quadrante " << quadrante_AA << endl;
		  total_quad[quadrante_AA-1][pos]++;
		  cout << endl;
	   }
	   cout << endl;
	}
	for(int i = 0; i < 20; i++)
		cout << teste[i] << " ";
	cout << endl;
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 20; j++){
			cout << total_quad[i][j] << "\t";
		}
		cout << endl;
	}
}

// CROSSOVER OPERATOR
void crossover(Atomo **pop, float fitnessvet[]){
	int i2;
	Atomo pn[tpopulacao][tam_sitiof];
	Atomo pais[2][tam_sitiof];
	int escolhidos_torneio [torneio];
	int escolhido;
	for(int k=0; k<tpopulacao; k+=2){
	   for(int i =0 ; i < 2; i++){
		   for(int i1=0;i1<torneio;i1++)
			  escolhidos_torneio[i1]= rand()%tpopulacao;
		   i2 = funcaotorneio(fitnessvet, escolhidos_torneio);
		   for(int j = 0;j<tam_sitiof;j++)
			  pais[i][j] = pop[i2][j];
	   }
	   if(rand()%100 <= txcruzamento){
		   int pointcross = rand() % tam_sitiof;
		   Atomo aux;
		   for(int i=0; i<=pointcross;i++){
			  aux = pais[0][i];
			  pais[0][i] = pais[1][i];
			  pais[1][i] = aux;
		   }
		}
	   for(int i = 0; i<tam_sitiof;i++){
		   pn[k][i] = pais[0][i];
		   pn[k+1][i] = pais[1][i];
	   }
	 }
	 for (int i = 0; i< tpopulacao; i++)
		 for (int j = 0 ; j < tam_sitiof; j++)
			  pop [i][j] = pn[i][j];
}

// TOURNEMENT
int funcaotorneio(float fitnessvet[], int escolhidos_torneio[]){
	 int escolhido = escolhidos_torneio[0];
	 float menor = fitnessvet[escolhido];
	 for(int i=1;i<torneio;i++)
		 if(fitnessvet[escolhidos_torneio[i]]< menor){
			 menor = fitnessvet[escolhidos_torneio[i]];
			 escolhido = escolhidos_torneio[i];
		 }
return escolhido;
}


// MUTATION OPERATOR
void mutation (Atomo **pop, int *numero_aminos, Atomo **repositorio, int quadrante){
   char find[5];
   int pos, j;
   int quadrante_AA;
   int tentativas = 0;
   for(int i=0; i < tpopulacao; i++){
	  	if(rand()%100<=txmutacao){
			for(j=0; j<tam_sitiof; j++){
		   		if(rand()%100<=50){
			 		strcpy(find, pop[i][j].amino);
			 		pos = posicaonalista(teste, find);
			 		while(1){
						pop[i][j] = repositorio[pos][rand() % numero_aminos[pos]];
						quadrante_AA = MostraQuadrante(pop[i][j], centroide);
						if(quadrante_AA == quadrante || quadrante == -1)
							break;
						else{
			  				tentativas++;
			  				if(tentativas > 1000)
			  					break;
			  			}
					}
		   		}
			}
	  	}	
   	}
}

// SUBSTITUTION
void substituicao (Atomo **pop, int *numero_aminos, Atomo **repositorio, char substitutos [30][4], int tam, Atomo *refe){
   int posicao_da_troca, pos = 0, sorteio[tam], w;

   for(int i = 0; i < tpopulacao; i++){
	  if(rand()%100<=5){   // era 10
		posicao_da_troca = rand()%tam_sitiof; // escolhe aleatoriamente a posicao no individuo para trocar
		w = 0;
		for(int j=0; j<tam;j+=2)
		   if(!strcmp(pop[i][posicao_da_troca].amino,substitutos[j]) || !strcmp(pop[i][posicao_da_troca].amino,substitutos[j+1]))
			   sorteio[w++]=j;
		if(w!=0){
			int posicao = rand() % w;
			if(!strcmp(refe[posicao_da_troca].amino,substitutos[sorteio[posicao]])){
			   pos = posicaonalista (teste, substitutos[sorteio[posicao]+1]);
			   pop[i][posicao_da_troca] = repositorio[pos][rand()%numero_aminos[pos]];
			}
			else if(!strcmp(refe[posicao_da_troca].amino,substitutos[sorteio[posicao]+1])){
				   pos = posicaonalista (teste, substitutos[sorteio[posicao]]);
				   pop[i][posicao_da_troca] = repositorio[pos][rand()%numero_aminos[pos]];
				 }
		}
	  }
   }
}

// FITNESS OF REFERENCE
float fitnessRef (Atomo *refe, float distanciasref[]){
   float fit = 0.0, aux;
   int pos = 0;
   for(int i=0;i<tam_sitiof;i++)
	  for(int j=i+1;j<tam_sitiof;j++){
		 aux = sqrt(pow(refe[i].x - refe[j].x, 2) + pow(refe[i].y - refe[j].y,2)  +  pow(refe[i].z - refe[j].z,2));
		 distanciasref[pos++] = aux;
		 fit = fit + aux;
	 }
   return fit;
}

// FITNESS OF POPULATION
void fitness(Atomo **p, float fitnessvet[], float distanciasitio[], float distanciasref[]){
   float fit = 0.0, aux;
   int pos;
   for(int k=0; k < tpopulacao; k++){
	  pos = 0;
	  for(int i=0; i<tam_sitiof; i++)
		 for(int j=i+1; j<tam_sitiof; j++){
			aux = sqrt(pow(p[k][i].x - p[k][j].x, 2)  + pow(p[k][i].y - p[k][j].y, 2)  + pow(p[k][i].z - p[k][j].z, 2));
			distanciasitio[pos++] = aux;
		 }
	  for(int i=0;i<tam_dist;i++)
		 fit = fit + fabs(distanciasitio[i]-distanciasref[i]);
	  fitnessvet [k]= fit;
	  fit = 0.0;

	  for(int j=0; j<tam_sitiof; j++)
		  for(int i=j+1; i < tam_sitiof; i++)
			 //if (p[k][j].atomo_ID == p[k][i].atomo_ID){
		 if ((p[k][j].atomo_ID == p[k][i].atomo_ID) && (p[k][j].cadeia == p[k][i].cadeia)){
				  fitnessvet[k] = fitnessvet[k] + 100.0;
				  i = tam_sitiof;
				  j = tam_sitiof;
			  }
   }
}

int posicaonalista (char teste[20][4], char find[5]){
	  for(int i=0; i<20; i++)
		 if (!strcmp(find, teste[i]))
			return i;
}

// BEST OF GENERATION
int melhordageracao(float fitnessvet[]){
	float melhor = fitnessvet[0];
	int posicao = 0;
	for (int i = 1; i<tpopulacao; i++)
		 if (melhor > fitnessvet[i]){
			melhor = fitnessvet[i];
			posicao = i;
		 }
	return posicao;
}

// FOUND THE BEST OF GENERATION
int melhordageracaof(float fitnessvet[]){
	  float melhor = fitnessvet[0];
	  int posicao = 0;
	  for (int i = 1; i<ngeracoes; i++)
		   if (melhor > fitnessvet[i]){
			  melhor = fitnessvet[i];
			  posicao = i;
		   }
	   return posicao;
}

// WORST OF GENERATION
int piordageracao(float fitnessvet[]){
	float melhor = fitnessvet[0];
	int posicao = 0;
	for (int i = 1; i<tpopulacao; i++)
		 if (melhor < fitnessvet[i]){
			melhor = fitnessvet[i];
			posicao = i;
		 }
	return posicao;
}

// WRITE ATOMS ON SCREEN
void EscreveAtomo2(Atomo atomo){
cout << atomo.amino << "  " << atomo.atomo_ID << "  " << atomo.x << "  " << atomo.y << "  " << atomo.z << "  " << atomo.cadeia << "  -  ";
}
void EscreveAtomo(Atomo atomo){
cout << atomo.amino << "  " << atomo.atomo_ID << "  ";
}

// QUICKSORT
void quickSort(Atomo **tempop, float a[], int l, int r){
   int j;
   if( l < r ){
	 j = partition( tempop, a, l, r);
	 quickSort(tempop, a, l, j-1);
	 quickSort(tempop, a, j+1, r);
   }
}

int partition(Atomo **tempop, float a[], int l, int r) {
   float pivot, t;
   int i, j;
   Atomo z;
   pivot = a[l];
   i = l; j = r+1;
   while(1) {
	do ++i; while( a[i] <= pivot && i <= r );
	do --j; while( a[j] > pivot );
	if( i >= j ) break;
	t = a[i]; a[i] = a[j]; a[j] = t;
	for(int k = 0; k<tam_sitiof; k++){
	   z = tempop[i][k];
	   tempop[i][k] = tempop[j][k];
	   tempop[j][k] = z;
	}
   }
   t = a[l]; a[l] = a[j]; a[j] = t;
   for(int k = 0; k<tam_sitiof; k++){
	   z = tempop[l][k];
	   tempop[l][k] = tempop[j][k];
	   tempop[j][k] = z;
   }
   return j;
}


// STARTING TOP
void IniciaTop(float *fittop){
	for(int j=0;j<ranking;j++)
	   fittop[j] = 1000;
}

// OPEN FILE
bool AbrirArquivo(ifstream &ifs, char nomeArquivo[]){
	ifs.open(nomeArquivo, ios::binary);
	return ifs.is_open();
}

void CriaListaSubstituicao(Substituicao Matriz[], char *argv){
	ifstream arquivotexto;  // variável para controlar o fluxo de entrada

	char localtemplate[200];
	strcpy(localtemplate,argv);
	strcat(localtemplate,"/SubstitutionMatrix.txt");
   // fin.open("Templates.txt");

	arquivotexto.open(localtemplate);
	int l1 = 0;
	if(!arquivotexto.is_open()){
	   cout<<"File can not be open: SubstitutionMatrix.txt.\n";

	   exit(1);
	}
	else{
		char linha[30];
		//int l1 = 0;
		while(!arquivotexto.eof()){
			 arquivotexto.getline(linha,30);
			 char temp[10];
			 int c1 = 0;
			 int c2 = 0;

			 while(linha[c2]!='\0'){
				  while(linha[c2]!=','){
					   temp[c1] = (linha[c2]);
					   c1++; c2++;
				  }
				  temp[c1]='\0';
				  strcpy(Matriz[l1].lista[0], temp);
				  //cout << temp << endl;

				  c1 = 0; c2++;
				  while(linha[c2]!=','){
					   temp[c1] = linha[c2];
					   c1++; c2++;
				  }
				  temp[c1]='\0';
				  strcpy(Matriz[l1].lista[1],temp);
				  //cout << temp << endl;

				  c1 = 0; c2++;
				  while(linha[c2]!='\0'){
					   temp[c1] = linha[c2];
					   c1++; c2++;
				  }
				  temp[c1]='\0';
				  strcpy(Matriz[l1].lista[2],temp);
				  //cout << temp << endl;
			 }
			 l1++;
		}
	}
}