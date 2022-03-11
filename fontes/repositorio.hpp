#pragma once

void GeraRepositorio(Atomo ***repositorio, int **numero_aminos){
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
				 repositorio[0][i][numero_aminos[0][i]] = aux;
				 repositorio[MostraQuadrante(aux, centroide)][i][numero_aminos[MostraQuadrante(aux, centroide)][i]] = aux;
				 numero_aminos[0][i]++;
				 numero_aminos[MostraQuadrante(aux, centroide)][i]++;
				 
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
