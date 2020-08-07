import subprocess
import logging
import pandas as pd
import argparse
import pickle
from pathlib import Path

def input_miRanda(interactionpair,save_dir,mrna_interact=False):
    '''
    input pandas dataframe interaction pair
    generate input file for miRanda 
    '''
    mir=interactionpair.loc[:,['miR_ID','miRNA_seq']].drop_duplicates()
    site=interactionpair.loc[:,['mRNA_ID','site_seq']].drop_duplicates()
    mrna=interactionpair.loc[:,['mRNA_ID','mRNA Seq']].drop_duplicates()
    
    print(interactionpair)
    with open(save_dir / 'all_site.fa','w') as w:
        word=''
        for _,e in site.iterrows():
            utrid=e['mRNA_ID']
            utrseq=e["site_seq"][::-1]
            word+='>'+utrid+'\n'
            word+=utrseq+'\n'
        w.write(word)
            
    with open(save_dir / 'all_mir.fa','w') as w:
        word=''
        for _,e in mir.iterrows():
            mirid=e['miR_ID']
            mirseq=e["miRNA_seq"]
            word+='>'+mirid+'\n'
            word+=mirseq+'\n'
        w.write(word)
    with open(save_dir / 'interaction_pair.txt','w') as w:
        for _,e in interactionpair.iterrows():
            utrid=e['mRNA_ID']
            siteid=e['mRNA_ID']
            mirid=e['miR_ID']
            
            w.write(mirid+'\t'+siteid+'\n')
            
    #test
    if mrna_interact:
        with open(save_dir / 'all_mrna.fa','w') as w:
            word=''
            for _,e in mrna.iterrows():
                utrid=e['mRNA_ID']
                utrseq=e["mRNA Seq"][::-1]
                word+='>'+utrid+'\n'
                word+=utrseq+'\n'
            w.write(word)



def interaction2fasta(interactionpair,save_dir):
    '''
    input pandas dataframe interaction pair
    generate input data for PITA and RNAhybrid
    '''
    prev_mireid=''
    for idx,e in interactionpair.iterrows():
        mirid=e['miR_ID']
        if(prev_mireid!=mirid):
            word=''
            prev_mireid=mirid
            mirseq=e["miRNA_seq"]
            word+='>'+mirid+'\n'
            word+=mirseq
            save_str=prev_mireid+'.fa'
            save_str1=prev_mireid+'_site.fa'

            
            with open(save_dir / save_str,'w') as w1:
                w1.write(word)
            try:
                f.close()
            except: 
                pass
            f = open(save_dir / save_str1, 'w')
        word=''
        siteid=e['mRNA_ID']
        siteseq=e["site_seq"][::-1]
        word+='>'+siteid+'_site\n'
        word+=siteseq+'\n'
        f.write(word) 
        
    f.close()

def p_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("all_data",help="pickle dataframe file for all data")
    parser.add_argument('pita_out',type=Path,help='pita input file dir')
    parser.add_argument('miranda_out',type=Path,help='miranda input file dir')
    parser.add_argument('rnahybrid_out',type=Path,help='rnahybrid input file dir')
    args=parser.parse_args()

    return args

def main():
    args=p_args()
    with open(args.all_data,'rb') as r:
        data=pickle.load(r)

    interaction2fasta(data,args.pita_out)
    interaction2fasta(data,args.rnahybrid_out)
    input_miRanda(data,args.miranda_out)


if __name__=='__main__':
    main()