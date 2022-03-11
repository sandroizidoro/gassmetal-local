#pragma once

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
			if(contador == 0) strcpy(nome_enzima, linha);
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
		//        EscreveAtomo(aux); cout << endl;
			 }
		  }
	 }
   }
   else{ cout << "File can not be open: " << proteina_sitiof << endl;}
   ifs.close();
}
