#pragma once

// STARTING POPULATION
void IniciaPopulacao(Atomo **pop,  Atomo **repositorio, int *numero_aminos, Atomo *refe, int inicio=0){
	// IMPLEMENTACAO B
	int pos;
	int total_quad[4][20] = {0};
	for(int i = inicio; i < tpopulacao; i++){
	   for(int j=0; j<tam_sitiof; j++){
			if(j == 0){
				//cout << "Individuo " << i+j << endl;
			}
			// Calcula o quadrante ideal de cada individuo e todos os AA desse individuo devem pertencer a esse quadrante
			// Caso nao tenha uma AA nesse quadrante, escolhe outro aleatoriamente.
			pos = posicaonalista(teste, refe[j].amino);
			pop[i][j] = repositorio[pos][rand() % numero_aminos[pos]];	
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