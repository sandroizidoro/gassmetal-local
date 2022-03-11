#pragma once 

void substituicao (Atomo **pop, int *numero_aminos, Atomo **repositorio, char** substitutos, int tam, Atomo *refe){
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
