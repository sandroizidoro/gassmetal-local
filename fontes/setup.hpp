#pragma once

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
