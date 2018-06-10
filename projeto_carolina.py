# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 15:50:01 2018

@author: carol
"""

# Anotacao do genoma
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Entrez
import pprint


#Genome annotation

def parsing_fasta_files(file):
    with open(file,"r") as ficheiro:
        for seq_record in SeqIO.parse(ficheiro,"fasta"):
            print(seq_record.id)
            print(repr(seq_record.seq))
            print(len(seq_record))
            

def blast_sequences(file):
    fasta_string = open(file).read()
    result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)


def genbank_parsing(file):
    for seq_record in SeqIO.parse(file, "genbank"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))
        

def genbank_extraction(file):
    record_iterator = SeqIO.parse(file, "genbank")
    first_record = next(record_iterator)
    print(first_record)
    print(first_record.annotations)

     
     
#Phylogenetic analysis 

#download sequences of rRNA 16S and danJ, gyrB, hsp60, recA and rpoB genes
def download_fasta_files(id_list,gene):
    Entrez.email = 'carolina.dias.silva@gmail.com'
    for id_gene in id_list:
        with open("C:/Users/carol/Desktop/mestrado/2o semestre/Projeto/files_scripts/fasta_files/"+gene+id_gene+".fasta","w") as f:
            handle = Entrez.efetch(db="nucleotide",id=id_gene,rettype="fasta",retmode="text")
            f.write(handle.read())

    


if __name__ == '__main__':
    #parsing_fasta_files("Siavash_KHB_HGAP2_rm.fasta")
    #parsing_fasta_files("Siavash_X5_HGAP2.fasta")
    #blast_sequences("Siavash_KHB_HGAP2_rm.fasta")
    #blast_sequences("Siavash_X5_HGAP2.fasta")
    #genbank_parsing("annotation_khb.gbk")
    #genbank_parsing("annotation_X5.gbk")
    #genbank_extraction("annotation_khb.gbk") 
    #genbank_extraction("annotation_X5.gbk")
    L_dnaJ = ['AB619862','AB619875','AB619877','AB619878','AB619886','AB619893','AB619894','AB619895','AB619896']
    L_gyrB = ['AB619900','AB619913','AB619915','AB619916','AB619924','AB619931','AB619932','AB619933','AB619934']
    L_hsp60=['AB510670','AB510677','AB510678','AB547556','AB510686','AB510692','AB510693','AB510694','AB510695']
    L_recA = ['AB619938','AB619951','AB619953','AB619954','AB619962','AB619969','AB619970','AB619971','AB619972']
    L_rpoB = ['AB619976','AB619989','AB619991','AB619992','AB620000','AB620007','AB620008','AB620009','AB620010']
    L_RNA16s = ['AB510698','AB510701','AB253732','AB547643','AB200229','AB510710','AB510711','AB510712','AB510713']
    #download_fasta_files(L_dnaJ,'dnaJ')
    #download_fasta_files(L_gyrB,'gyrB')
    #download_fasta_files(L_hsp60,'hsp60')
    #download_fasta_files(L_recA,'recA')
    #download_fasta_files(L_rpoB,'rpoB')
    #download_fasta_files(L_RNA16s,'RNA16s')
    
    
    
    
    
    
    
    
    