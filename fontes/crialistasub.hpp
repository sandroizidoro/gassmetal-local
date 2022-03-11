#pragma once

int ListaSubstituicao(char **substitutos ){
   int i = 0, ok = 0;
   char proteina [5];
   for(int j=0;j<4;j++)
	  proteina[j] = proteina_sitiof[j];
   proteina[4] = '\0';
   for(int k= 0 ; k < 567; k++){
	   if(!strcmp(proteina, Matriz[k].lista[0])){
		 strcpy(substitutos[i], Matriz[k].lista[1]);
		 strcpy(substitutos[i+1], Matriz[k].lista[2]);
		 i+=2;
		 ok = 1;
	   }
	   else {
			if(ok == 1)
			   k = 567;
	   }
   }
   return i;
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
					   temp[c1] = toupper(linha[c2]);
					   c1++; c2++;
				  }
				  temp[c1]='\0';
				  strcpy(Matriz[l1].lista[0], temp);

				  c1 = 0; c2++;
				  while(linha[c2]!=','){
					   temp[c1] = linha[c2];
					   c1++; c2++;
				  }
				  temp[c1]='\0';
				  strcpy(Matriz[l1].lista[1],temp);

				  c1 = 0; c2++;
				  while(linha[c2]!='\0'){
					   temp[c1] = linha[c2];
					   c1++; c2++;
				  }
				  temp[c1]='\0';
				  strcpy(Matriz[l1].lista[2],temp);
			 }
			 l1++;
		}
	}
}