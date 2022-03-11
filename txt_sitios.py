'''
Script que gera um arquitvo .txt de todos os .pdb do diretório para os templates do GASS.
O arquivo .txt contem apenas os atomos de CA dos residuos dos templates dos sitios metalicos.
Para gerar o .txt da proteina inteira, utilizar o script 'completo_dat.py'
'''

import os
arq_templates = open("sortSitiosLIT.txt","r")

templates = arq_templates.readlines()
templates = [t.replace('\n','').split() for t in templates]
cont = 1
for t in templates:
    cont = cont + 1
    nome_txt = t[0] + '_.txt'
    nome_pdb = t[0] + '.pdb'
    sitio = [t[s] for s in range(3,len(t))] # arquivo de sitios onde o primeiro residuo comeca no indice 3
    try:
        arq_pdb = open(nome_pdb, 'r')
        linhas_pdb = arq_pdb.readlines()
        linhas_pdb = [linha.replace('\n','') for linha in linhas_pdb]
        arq_pdb.close()
    except:
        continue
    if os.path.exists(nome_txt) == False:
        arq_txt = open(nome_txt,'a')
        ec_proteina = 'NULL'
        unp_proteina = 'NULL'
        for l_pdb in linhas_pdb:
            if l_pdb.find('HEADER') != -1:
                nome_proteina = l_pdb.strip().split()[-1]
            if l_pdb.find('EC:') != -1 and l_pdb.find('COMPND') != -1:
                ec_proteina = l_pdb.strip().split()[-1].replace(';','')
            if l_pdb.find('UNP') != -1 and l_pdb.find('DBREF') != -1:
                unp_proteina = l_pdb.strip().split()[6]
            if l_pdb.find('RESOLUTION.') != -1 and l_pdb.find('REMARK') != -1:
                reso_proteina = l_pdb.strip().split()[3]
        arq_txt.write(nome_proteina)
        arq_txt.write('\n')
        arq_txt.write(ec_proteina)
        arq_txt.write('\n')
        arq_txt.write(unp_proteina)
        arq_txt.write('\n')
        arq_txt.write(reso_proteina)
        arq_txt.write('\n')
        
    arq_txt = open(nome_txt,'a')
    # loops que percorrem todos os atomos de CA e verifica se existem dados duplicados    
    for r in range(0, len(sitio),3):
        for i in range(0,len(linhas_pdb)):
            encontrou = False
            aux_i = linhas_pdb[i].strip().split()
            if aux_i[0] == 'ATOM' and aux_i[2] == 'CA':
                if len(aux_i[4]) > 1:
                    tmp = aux_i[4]
                    del aux_i[4]
                    aux_i.insert(4,tmp[0])
                    aux_i.insert(5,tmp[1:]) 
                if aux_i[3] == sitio[r] and aux_i[5] == sitio[r+1] and aux_i[4] == sitio[r+2]:
                    for j in range(i+1,len(linhas_pdb)):
                        aux_j = linhas_pdb[j].strip().split() 
                        if aux_j[0] == 'ATOM' and aux_j[2] == 'CA':
                            if len(aux_j[4]) > 1:
                                tmp = aux_j[4]
                                del aux_j[4]
                                aux_j.insert(4,tmp[0])
                                aux_j.insert(5,tmp[1:])
                            if aux_j[3] == sitio[r] and aux_j[5] == sitio[r+1] and aux_j[4] == sitio[r+2]:
                                if aux_i[4] == aux_j[4] and aux_i[5] == aux_j[5]:
                                    encontrou = True 
                    if encontrou == False:
                        if(len(aux_i[3]) > 3):
                            linhas_pdb[i] = linhas_pdb[i][:16] + " " + linhas_pdb[i][17:]
                        arq_txt.write("%s\n" % linhas_pdb[i])
    arq_txt.close()

