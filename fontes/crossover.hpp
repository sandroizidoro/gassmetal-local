#pragma once

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
