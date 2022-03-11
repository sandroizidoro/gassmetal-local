// filtroC01.cpp
// Importa os dados gerados pelo filtroperl01.pl para um arquivo
// que possa ser manipulado na linguagem C
// Ultima atualização: 10/08/11

#include <iostream>
#include <fstream>
#include <string.h>
#include <cstdlib>

using std::ifstream;
using std::ofstream;
using namespace std;

struct Atomo{
    char amino[5];
    char atomo[4];
    char cadeia;
    int atomo_ID;
    float x,y,z;
};

void transfer(char *linha, Atomo &aux);

int main(){
    char nome_proteina[10], nome_arquivo[10];
	char nometemp [6];

    ifstream arquivostexto;  // variável para controlar o fluxo de entrada
    char buff[10];  // variável para receber cada linha do arquivo
    int nlinha = 0; // conta o numero de linhas
    Atomo aux;

	arquivostexto.open("LISTAENZIMAS.TXT");  // associa um arquivo com o fluxo de entrada
	if(!arquivostexto.is_open()){    // verifica se o arquivo pode ser aberto para leitura
      // caso haja algum erro para abrir para leitura
      cout << "Arquivo LISTAENZIMAS.TXT não pode ser aberto.\n";
      exit(1);
    }
    else{
        while(!arquivostexto.eof()){
             arquivostexto.getline(buff,10);  // le uma linha do arquivo, no máximo 255 caracteres
 			 strncpy(nome_proteina, buff, 5);
			 strncpy(nome_arquivo, buff,10);

			 nome_proteina[5]='\0';
			 strcat(nome_proteina, ".dat" );
			 nome_proteina[9]='\0';

			 ifstream inFileProteina;
             cout <<"escrita: "<<nome_proteina<< "   -   " << nome_arquivo << endl;

			 inFileProteina.open(nome_arquivo);
   			 if(!inFileProteina){
                cout<<"Erro ao abrir o arquivo"<<endl;
                exit(1);
             }
			  else{
			       char linha[100];
                   ofstream outFile(nome_proteina,ios::binary);
                   int countline = 0;
                   while(!inFileProteina.eof()){
                        inFileProteina.getline(linha,100);
                        if(countline<=3){
                           outFile.write(reinterpret_cast<const char *>(&linha), sizeof(linha));
                           countline++;

                         } else{
                               transfer(linha,aux);
                               countline++;
                               outFile.write(reinterpret_cast<const char *>(&aux), sizeof(Atomo));
                         }
                    }
                    outFile.write(reinterpret_cast<const char *>(&countline), sizeof(int));
                    outFile.close();
                    inFileProteina.close();
				}
          }
        }
	arquivostexto.close(); // fecha o fluxo de entrada
    return 0;
}

void transfer(char *linha,Atomo &aux){
    aux.cadeia=linha[21]; //Adicionando a cadeia
    //Pegando o nome do atomo;
    char temp5[4];
    int j=0;
    for(int i=12;i<16;i++){
        if(linha[i]!=' '){
            temp5[j]=linha[i];
            j++;
        }
    }
    temp5[j]='\0';
    strcpy(aux.atomo,temp5);

    // Pegando o nome do aminoacido;
    char temp[4];
    j=0;
    for(int i=17;i<20;i++){
        if(linha[i]!=' '){
            temp[j]=linha[i];
            j++;
        }
    }
    temp[j]='\0';
    strcpy(aux.amino,temp);

    // Coordenada X;
    char temp2[8];
    j=0;
    for(int i=30;i<38;i++){
        if(linha[i]!=' '){
            temp2[j]=linha[i];
            j++;
        }
    }
    temp2[j]='\0';
    float X = atof(temp2);
    aux.x=X;

    //Coordenada Y
    char temp3[8];
    j=0;
    for(int i=38;i<46;i++){
        if(linha[i]!=' '){
            temp3[j]=linha[i];
            j++;
        }
    }
    temp3[j]='\0';
    float Y = atof(temp3);
    aux.y=Y;

    // Coordenada Z
    char temp4[8];
    j=0;
    for(int i=46;i<54;i++){
        if(linha[i]!=' '){
            temp4[j]=linha[i];
            j++;
        }
    }
    temp4[j]='\0';
    float Z = atof(temp4);
    aux.z=Z;

    // ID do átomo
    char temp6[4];
    j=0;
    for(int i=22;i<26;i++){
        if(linha[i]!=' '){
            temp6[j]=linha[i];
            j++;
        }
    }
    temp6[j]='\0';
    int ID = atoi(temp6);
    aux.atomo_ID=ID;
};

