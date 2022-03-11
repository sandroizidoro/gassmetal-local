#########################################################
#Nome: pdb_to_txt.                                      #
#Descrição: Gera a partir de arquivos PDB, arquivos     #
#texto com nome da enzima, EC, e posição do último      #
#átomo mais pesado da cadeia lateral.                   #
#Data: 27/09/2013                                       #
#########################################################

#Inicializando vetor de comparação de aminoacidos contendo o atomo mais pesado de cada aminoácido
open (Amino, "AtomoMaisPesado.txt") or die("Arquivo AtomoMaisPesado.txt não encontrado.\n");
@aminoacidos;
while(<Amino>){
     chomp $_;  
     push(@aminoacidos,$_);
};
close(Amino);
pop(@aminoacidos); 
$tam = @aminoacidos;

# Abrindo diretório
use Cwd qw(getcwd);
my $diretorio = getcwd($0);
opendir(MEUDIR,$diretorio);
@pegoodir = readdir(MEUDIR);
closedir (MEUDIR);

foreach (@pegoodir){
    $dados = $_;
    if($dados =~ /.pdb$/){	
       $ArquivoE=$dados;
       $ArquivoS=substr($dados,0,4).".txt";
       chomp $ArquivoE;
       chomp $ArquivoS;
       #print $ArquivoS;
       open (File, "$ArquivoE") or die ("Arquivo de entrada não encontrado.\n");
       open (my $Out, ">$ArquivoS") or die ("Impossível abrir arquivo de saída.\n");
       $countline = 0;
       $nome = substr($dados,0,4);
       $contador = 0;
       $ec = 0; $Resolucao = 0; $uni = 0;
       while(<File>){
          if($_=~/^HEADER/){
             @matriz = split(" ",$_); 
             $tamanho = $#matriz;
             $NomeProteina = $matriz[$tamanho];
             #print $Out $teste;
	     $countline++;
	  }  elsif(($_ =~ /EC:/)and($ec==0)){
                 @matriz = split(" ",$_); 
                 $ECNumber = $matriz[3];
                 $temp = length($ECNumber);  
                 #print substr($ECNumber,0,$merda-1);
                 #print " $merda \n";
                 if(substr($ECNumber,$temp-1,1) == ";"){
                    $ECNumber = substr($ECNumber,0,$temp-1);
                 }
                 #print $Out $teste;
   	         $countline++;
	         $ec = 1;
             } elsif(($_ =~ /REMARK   2 /) and ($_ =~ /RESOLUTION./)and($Resolucao==0)){
                 @matriz = split(" ",$_); 
                 $RESOLUTION = $matriz[3];
                 #print $Out $teste;
   	         $countline++;
	         $Resolucao = 1;
               } elsif(($_ =~ /DBREF/)and($uniprot==0)and($_ =~ /UNP/)){
                   @matriz = split(" ",$_); 
                   $UNIPROT = $matriz[6];
                   # print $Out $teste;
   	           $countline++;
	           $uni = 1;
                } else {$countline++;}
       }
       close (File);

       open (File, "$ArquivoE") or die ("Arquivo de entrada não encontrado.\n");
       print $Out "$NomeProteina\n";
       if($ec == 1){print $Out "$ECNumber\n";} else { print $Out "NULL\n";}
       if($uni== 1){print $Out "$UNIPROT\n";} else { print $Out "NULL\n";}
       if($Resolucao == 1){print $Out "$RESOLUTION\n";} else { print $Out "NULL\n";}
       while(<File>){
           if($_=~/^ATOM/){
   	      $cadeia = substr($_,21,1);
   	      for($i=0;$i<$tam;$i=$i+2){
                 $teste = $aminoacidos[$i];
		 $teste2 = $aminoacidos[$i+1];
		 chomp $teste; chomp $teste2;
		 if(($_=~/$teste/)and($_=~/$teste2/)){
                    if($countline < 2){  
		       print $Out "\n$_";
                    } else{
		          print $Out $_;
		    }
 		    $countline++;
  		 }
	     }
         }
       }
       close (File);
       close $Out;
    }
}
