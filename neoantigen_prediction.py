# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 13:43:01 2022

@author: zhanglu
"""

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from collections import OrderedDict
import re
from Aligner import Aligner


def read_fasta(path):
    seq_dic = OrderedDict()
    for line in open(path).xreadlines():
        if line.startswith('>'):
            nm = line.strip().split('|')[1].split('.')[0]
            seq = ''
            seq_dic[nm] = seq
            
        else:
            seq_dic[nm] += line.strip()
    return seq_dic

def chopchop(aaSeq, peptide_length):
    peptides = []    
    for i in range(len(aaSeq)):
        pep = aaSeq[i:i + peptide_length]
        if len(pep) < peptide_length:           
            break
        peptides.append(pep)
    return peptides

def get_missense():
    b16_f = open(r'C:/D/pMHC/mutation/mutation/common_between_DNA_RNA/B16-DNARNAcommon.snp.annot.csv')
    mouse_f = open(r'C:/D/pMHC/mutation/mutation/common_between_DNA_RNA/Mouse_DNARNAcommon.snp.annot.csv')
    fo = open(r'C:/D/pMHC/mutation/mutation/common_between_DNA_RNA/missense_snp.txt','w')
    fo.write(b16_f.readline()+'\n')
    for line in b16_f.readlines()[1:]:
        line = line.strip().split()
        if line[2] == 'MISSENSE':
            transcript = line[3].split(':')[0]
            aa_change = line[3].split(':')[1].split('/')[0][2:]
            nt_change = line[3].split(':')[1].split('/')[1][2:]
            fo.write('b16'+'\t'+line[1]+'\t'+transcript+'\t'+aa_change+'\t'+nt_change+'\t'+'\t'.join(line[4:7])+'\t'+'\t'.join(line[7:12])+'\n')
    b16_f.close()
    
    for line in mouse_f.readlines()[1:]:
        line = line.strip().split()
        if line[2] == 'MISSENSE':
            transcript = line[3].split(':')[0]
            aa_change = line[3].split(':')[1].split('/')[0][2:]
            nt_change = line[3].split(':')[1].split('/')[1][2:]
            fo.write('mouse'+'\t'+line[1]+'\t'+transcript+'\t'+aa_change+'\t'+nt_change+'\t'+'\t'.join(line[4:7])+'\t'+'\t'.join(line[7:12])+'\n')
    mouse_f.close()
    fo.close()
def get_frameshift():
    b16_f = open(r'C:/D/pMHC/mutation/mutation/common_between_DNA_RNA/B16-DNARNAcommon.indel.annot.csv')
    mouse_f = open(r'C:/D/pMHC/mutation/mutation/common_between_DNA_RNA/Mouse_DNARNAcommon.indel.annot.csv')
    fo = open(r'C:/D/pMHC/mutation/mutation/common_between_DNA_RNA/frameshift_indel.txt','w')
    fo.write(b16_f.readline())
    for line in b16_f.readlines()[1:]:
        line = line.strip().split()
        if line[0] == 'frameshift_variant':
            transcript = line[3].split(':')[0]
            aa_change = line[3].split(':')[1].split('/')[0][2:]
            
            nt_change = line[3].split(':')[1].split('/')[1][2:]
            fo.write('b16'+'\t'+line[1]+'\t'+transcript+'\t'+aa_change+'\t'+nt_change+'\t'+'\t'.join(line[4:7])+'\t'+'\t'.join(line[7:12])+'\n')
    b16_f.close()
    
    for line in mouse_f.readlines()[1:]:
        line = line.strip().split()
        if line[0] == 'frameshift_variant':
            transcript = line[3].split(':')[0]
            
            aa_change = line[3].split(':')[1].split('/')[0][2:]
            nt_change = line[3].split(':')[1].split('/')[1][2:]
            fo.write('mouse'+'\t'+line[1]+'\t'+transcript+'\t'+aa_change+'\t'+nt_change+'\t'+'\t'.join(line[4:7])+'\t'+'\t'.join(line[7:12])+'\n')
    mouse_f.close()
    fo.close()        

def get_RawMutSeq():
    protein_seq_path = r'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/B16_mouse_missense_snp_protein.txt'
    nm_seq_dic = read_fasta(protein_seq_path)
    snp_record_f = open(r'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/B16_mouse_missense_snp.txt')
    snp_seq_f = open(r'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/snp_seq_29mer.txt','w')
 
    for line in snp_record_f.readlines()[1:]:
        line = line.strip().split('\t')
        gene_name = line[1]
        nm = line[3]
        subinfo = line[4][2:]
        aa_dict = {
        'GLY': 'G', 'SER': 'S', 'ALA': 'A', 'THR': 'T', 'VAL': 'V',\
        'ILE': 'I', 'LEU': 'L', 'TYR': 'Y', 'PHE': 'F', 'HIS': 'H',\
        'PRO': 'P', 'ASP': 'D', 'MET': 'M', 'GLU': 'E', 'TRP': 'W',\
        'LYS': 'K', 'CYS': 'C', 'ARG': 'R', 'ASN': 'N', 'GLN': 'Q'}
        raw_aa = aa_dict[subinfo[:3].upper()]
        mut_aa = aa_dict[subinfo[-3:].upper()]
        site = int(subinfo[3:-3])
        protein_seq = nm_seq_dic[nm]
        if protein_seq[site-1] == raw_aa:
            if site-15<0:
                raw_seq = protein_seq[:site+14]
                mut_seq = raw_seq[:site-1]+mut_aa+raw_seq[site:]
            else:                
                raw_seq = protein_seq[site-15:site+14]            
                mut_seq = raw_seq[:14]+mut_aa+raw_seq[15:]
            if site-28<0:
                mut_seq_55 = protein_seq[:site-1]+ mut_aa+ protein_seq[site:site+27]
            else:               
                mut_seq_55 = protein_seq[site-28:site-1]+ mut_aa+ protein_seq[site:site+27] 
        snp_seq_f.write(gene_name+'\t'+nm+'\t'+subinfo+'\t'+raw_seq+'\t'+mut_seq+'\t'+mut_seq_55+'\n')
        #raw_seq,mut_seq,mut_seq_45 = get_seq(subinfo,nm,nm_seq_dic)
    snp_record_f.close()    
    snp_seq_f.close()    


########### step 1 get raw/mut 29mer peptide #################
snp_seq_f = open(r'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/snp_seq_29mer.txt')
#raw_seq_f = open(r'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/raw_seq_29mer.txt','w')
mut_seq_f = open(r'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/mut_seq_29mer.txt','w')
for line in snp_seq_f.readlines():
    line = line.strip().split('\t')
    info = '_'.join(line[:2])
    raw_seq = line[3]
    mut_seq = line[4]
    #raw_seq_f.write('>'+info+'\n'+raw_seq+'\n')
    mut_seq_f.write('>'+info+'\n'+mut_seq+'\n')
snp_seq_f.close()
#raw_seq_f.close()
mut_seq_f.close()    


############# step 2 select peptides that with affinity<=500 #################### 
allele1_li = ['H2-Dd','H2-Kd','H2-Ld','H2-Db','H2-Kb']
allele2_li = ['H2-IAb','H2-IAd']
snp_seq_f = open(r'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/snp_seq_29mer.txt')
mhcpan_path1 = 'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/netMHCout/'
mhcpan_path2 = 'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/netMHCIIout/'
gene_aasub_dic = OrderedDict()
gene_pep_dic = OrderedDict()
for line in snp_seq_f.readlines():
    line = line.strip().split('\t')
    gene = line[0]
    aasub = line[2]
    pep = line[-1]
    gene_aasub_dic[gene] = aasub
    gene_pep_dic[gene] = pep
snp_seq_f.close()

for allele in allele1_li:
    mut_f = open(mhcpan_path1+allele+'_mut.xls')
    raw_f = open(mhcpan_path1+allele+'_raw.xls')
    out_f = open(mhcpan_path1+allele+'_selected.txt','w')
    raw_line = raw_f.readline()
    raw_line = raw_f.readline()
    for line in mut_f.readlines()[2:]:
        line = line.strip().split('\t')
        pep = line[1]
        gene = line[2].split('_')[0]
        mut_aff = float(line[6])
        raw_aff = float(raw_f.readline().strip().split('\t')[6])
        if mut_aff<=500:            
            minus = raw_aff-mut_aff
            divide = raw_aff/mut_aff
            out_f.write(gene+'\t'+gene_aasub_dic[gene]+'\t'+pep+'\t'+gene_pep_dic[gene]+'\t'+str(mut_aff)+'\t'+str(raw_aff)+'\t'+str(minus)+'\t'+str(divide)+'\n')
    mut_f.close()
    raw_f.close()
    out_f.close()
     
for allele in allele2_li:
    mut_f = open(mhcpan_path2+allele+'_mut.xls')
    raw_f = open(mhcpan_path2+allele+'_raw.xls')
    out_f = open(mhcpan_path2+allele+'_selected.txt','w')
    raw_dic = OrderedDict()
    for line in raw_f.readlines()[2:]:
        line = line.strip().split('\t')
        info = line[2]+'_'+line[0]
        raw_dic[info] = float(line[4])
               
    for line in mut_f.readlines()[2:]:
        line = line.strip().split('\t')
        if len(line[1]) == 15:
            info = line[2]+'_'+line[0]
            pep = line[1]
            gene = line[2].split('_')[0]
            mut_aff = float(line[4])
            raw_aff = raw_dic[info]
            if mut_aff<=500:            
                minus = raw_aff-mut_aff
                divide = raw_aff/mut_aff
                out_f.write(gene+'\t'+gene_aasub_dic[gene]+'\t'+pep+'\t'+gene_pep_dic[gene]+'\t'+str(mut_aff)+'\t'+str(raw_aff)+'\t'+str(minus)+'\t'+str(divide)+'\n')
    mut_f.close()
    raw_f.close()
    out_f.close()            
    

############# step3: get all allele info and seq####################
allele1_li = ['H2-Dd','H2-Kd','H2-Ld','H2-Db','H2-Kb']
allele2_li = ['H2-IAb','H2-IAd']
mhcpan_path1 = 'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/netMHCout/'
mhcpan_path2 = 'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/netMHCIIout/'

allele1_fo1 = open(mhcpan_path1+'all_selected.txt','w')
allele1_fo2 = open(mhcpan_path1+'all_selected.fasta','w')
allele2_fo1 = open(mhcpan_path2+'all_selected.txt','w')
allele2_fo2 = open(mhcpan_path2+'all_selected.fasta','w')

for allele in allele1_li:
    f = open(mhcpan_path1+allele+'_selected.txt')
    tmp_info = []
    for line in f.readlines():
        line = line.strip().split('\t')
        info = allele+'_'+line[0]
        tmp_info.append(info)
        count = tmp_info.count(info)
        if count>1:
            info = info+'_'+str(count)
        allele1_fo1.write(info+'\t'+'\t'.join(line[1:])+'\n')
        allele1_fo2.write('>'+info+'\n'+line[2]+'\n')
    f.close()

for allele in allele2_li:
    f = open(mhcpan_path2+allele+'_selected.txt')
    tmp_info = []
    for line in f.readlines():
        line = line.strip().split('\t')
        info = allele+'_'+line[0]
        tmp_info.append(info)
        count = tmp_info.count(info)
        if count>1:
            info = info+'_'+str(count)
        allele2_fo1.write(info+'\t'+'\t'.join(line[1:])+'\n')
        allele2_fo2.write('>'+info+'\n'+line[2]+'\n')
    f.close()        
allele1_fo1.close()
allele1_fo2.close()
allele2_fo1.close()
allele2_fo2.close()


################# step4: Compute TCR-recognition probabilities for all neoantigens ############
def calu_tcr_pro(fpath,sample):
    fo = open(fpath+sample+'out.txt','w')
    a = 26
    k = 1
    aligner=Aligner()    
    aligner.readAllBlastAlignments(fpath+sample+'.xml')    
    aligner.computeR(a, k)    
    
    Ai=OrderedDict()
    data = OrderedDict()
    selected_f = open(fpath+sample+'.txt')
    for line in selected_f.readlines():
        line = line.strip().split('\t')
        Ai[line[0]] = line[7]
        data[line[0]] = line[2]
    
    #Compute neoantigen quality
    nids=Ai.keys()
    header=["NeoantigenID","MT.Peptide.Form","NeoantigenQuality",
            "NeoantigenAlignment","IEDB_EpitopeAlignment","AlignmentScore"]
    header="\t".join(header)
    fo.write(header+'\n')
    for i in nids:
        A=float(Ai[i])
        [R,species,alignment]=aligner.getR(i)
        
        neoAlignment=alignment[0]
        epitopeAlignment=alignment[1]
        score=alignment[2]
        #print A,R
        l=[i, data[i], A*R, neoAlignment, epitopeAlignment, score]
        l="\t".join(map(lambda s: str(s),l))
        #print l
        fo.write(l+'\n')        
    fo.close()

fpath = 'C:/D/pMHC/mutation/common_between_DNA_RNA/selected_mut/netMHCIIout/'
calu_tcr_pro(fpath,'mhc2_all_selected')





