import pandas as pd
import argparse
from utils_annotate import add_mock_id,add_site_id,rename_utrseq,rename_mirseq,seq2index,load_miRNA2seq
from generate_mock import gen_mock_mir,gen_mock_data
from comparetools import input_miRanda
import subprocess
import logging
from pathlib import Path
import random

def load_deepmirtar(pos_file,neg_file):
    pos=pd.read_excel(pos_file)
    # rename miRNA ID
    pos['miRNA ID']=pos['miRNA ID'].map(lambda x:x.replace('U','T'))
    neg=pd.read_excel(neg_file)

    neg=neg.rename(columns={"miRNA ID":'miR_ID',"miRNA Seq (mock, 3'-->5')":'miRNA_seq',"mRNA Accession Number":"mRNA_ID","Target Site":"site_seq",\
            "Start":"Start_position",'End':'End_position'})

    neg=neg.drop_duplicates(['miRNA_seq','site_seq'])
    neg['miR_ID']=neg['miR_ID'].map(lambda x:x+'_mock')
    neg['label']=len(neg['miR_ID'])*[0]
    #sort dataframe to add mock easily
    neg=neg.sort_values('miR_ID')
    neg['miR_ID']=add_mock_id(neg['miR_ID'])

    pos=pos.rename(columns={"miRNA ID":'miR_ID',"miRNA Seq (3'-->5')":'miRNA_seq',"mRNA Accession Number":"mRNA_ID","Target Site":"site_seq",\
            "Start":"Start_position",'End':'End_position'})
    pos=pos.drop_duplicates(['miRNA_seq','site_seq'])
    pos['label']=len(pos)*[1]

    all_deepmirtar=pd.concat((pos,neg)).drop_duplicates(['miRNA_seq','site_seq']).reset_index(drop=True)

    all_deepmirtar['miRNA_seq']=all_deepmirtar['miRNA_seq'].map(rename_mirseq)
    all_deepmirtar['site_seq']=all_deepmirtar['site_seq'].map(rename_utrseq)
    all_deepmirtar['mRNA Seq']=all_deepmirtar['mRNA Seq'].map(rename_utrseq)


    #in oreder to match seed sequence pairing
    # get miRNA 5'->3' sequence
    all_deepmirtar['miRNA_seq']=all_deepmirtar['miRNA_seq'].map(lambda x:x[::-1])
    # get mRNA 3' -> 5' sequence
    all_deepmirtar['mRNA Seq']=all_deepmirtar['mRNA Seq'].map(rename_utrseq).map(lambda x:x[::-1])
    # get site 3' -> 5' sequence
    all_deepmirtar['site_seq']=all_deepmirtar['site_seq'].map(rename_utrseq).map(lambda x:x[::-1])
    # sort dataframe to add site id easily
    all_deepmirtar=all_deepmirtar.sort_values('mRNA_ID')
    all_deepmirtar['mRNA_ID']=add_site_id(all_deepmirtar['mRNA_ID'])

    all_deepmirtar=all_deepmirtar.sort_values(by=['miR_ID']).reset_index(drop=True)
    all_deepmirtar['miRNA_seq_len']=all_deepmirtar['miRNA_seq'].map(len)
    all_deepmirtar['site_seq_len']=all_deepmirtar['site_seq'].map(len)
    all_deepmirtar['mir_idx']=all_deepmirtar['miRNA_seq'].map(seq2index)
    all_deepmirtar['site_idx']=all_deepmirtar['site_seq'].map(seq2index)

    return all_deepmirtar

def load_mirtarbase(pos_file,v22_miRNA2seq):
    target_site=pd.read_excel(pos_file)
    # get homo sapiens data
    mi2site=target_site.loc[target_site['Species (miRNA)']=='Homo sapiens']
    mi2site=mi2site.rename(columns={'Target Site':'site_seq','miRNA':'miR_ID'})


    def annotate_mir(mir):
        '''
        check if mirtarbase mir id exist in v22 mirBase
        '''
        if mir in v22_miRNA2seq:
            return v22_miRNA2seq[mir]
        else:
            return False

    #to match seed pairing
    #get mir sequence 5'->3'
    # site sequence 3'-> 5'
    mi2site['site_seq']=mi2site['site_seq'][::-1]
    mi2site['miRNA_seq']=mi2site['miR_ID'].map(annotate_mir)
    mi2site=mi2site.loc[mi2site['miRNA_seq']!=False]

    mi2site['miRNA_seq_len']=mi2site['miRNA_seq'].map(len)
    mi2site['site_seq_len']=mi2site['site_seq'].map(len)


    mi2site=mi2site.loc[(mi2site['site_seq_len']<=50) & (mi2site['miRNA_seq_len']<=30)]

    mi2site['mir_idx']=mi2site['miRNA_seq'].map(seq2index)
    mi2site['site_idx']=mi2site['site_seq'].map(seq2index)
    mi2site['label']=[1]*len(mi2site['site_idx'])


    mi2site=mi2site.sort_values('Target Gene')
    mi2site['mRNA_ID']=add_site_id(mi2site['Target Gene'])


    mi2site=mi2site.drop_duplicates(['site_seq','miRNA_seq'])
    mi2site=mi2site.reset_index(drop=True)

    return mi2site

def seed_27_38(v22_miRNA2seq):
    seed_27=[]
    seed_38=[]
    for i in v22_miRNA2seq.values():
        seed_27.append(i[1:7])
        seed_38.append(i[2:8])
    return seed_27,seed_38


def filter_pos(all_pos):
    '''
    filter dataset by miRanda
    '''
    def check_filter_score(file,original_df):
        '''
        filter out pairs that miranda interaction score < 100  
        '''
        no_hits={}
        with open(file,'r') as r:
            line=r.readline()
            while line:
                if(line.startswith('Performing Scan:')):
                    mirna=line.split(' ')[2]
                    site=line.split(' ')[4].rstrip()
                elif(line.startswith('No Hits')):
                    if mirna in no_hits:
                        no_hits[mirna].append(site)
                    else:
                        no_hits[mirna]=[site]
                else:
                    pass
                line=r.readline()
        arr=[]
        for i,e in original_df.iterrows():
            if e['miR_ID'] in no_hits:
                if e['mRNA_ID'] in no_hits[e['miR_ID']]:
                    arr.append(False)
                else:
                    arr.append(True)
            else:
                arr.append(True)
        return arr

    input_miRanda(all_pos,Path('./tmp/'))
    subprocess.call(["miranda", "./tmp/all_mir.fa", "./tmp/all_site.fa","-sc","100","-restrict","./tmp/interaction_pair.txt","-out","./tmp/miranda_result.txt"])
    filter_idx=check_filter_score('./tmp/miranda_result.txt',all_pos)
    subprocess.call(["rm","-r",'./tmp'])
    all_pos_filter=all_pos.loc[filter_idx]

    return all_pos_filter



def p_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("deepmirtar_pos",help="deepmirtar positive interaction pairs data  (source: https://academic.oup.com/bioinformatics/article/34/22/3781/5026656#supplementary-data)")
    parser.add_argument("deepmirtar_neg",help="deepmirtar negative interaction mock pairs data  (source: https://academic.oup.com/bioinformatics/article/34/22/3781/5026656#supplementary-data)")
    parser.add_argument("mirtarbase_pos",help="mirtarbase positive interaction pair data  (source: http://mirtarbase.cuhk.edu.cn/php/download.php)")
    parser.add_argument("v22_mirna",help="mirBase v22 mature microRNA file (source: http://www.mirbase.org/ftp.shtml)")
    parser.add_argument("output_file",help="output whole dataset in pickle format after adding mock data")
    parser.add_argument("--mock_threshold",default=2000, help="maximun shuffle times for mock data generation")
    parser.add_argument("--filter",default=False,type=bool,help="determine if the positive data filted by miranda score 100")
    parser.add_argument("--seed",default=9487,type=bool,help="seed for fisher shuffle in order to reproduce same mock data")
    parser.add_argument("--only_deepmirtar",default=False,type=bool,help="generate only deepmirtar data")

    args = parser.parse_args()
    return args


def main():
    args=p_args()
    logging.basicConfig(format='%(asctime)s | %(levelname)s | %(message)s',
                        level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

    v22_miRNA2seq=load_miRNA2seq(args.v22_mirna)

    logging.info('loading DeepMirTar data')
    all_deepmirtar=load_deepmirtar(args.deepmirtar_pos,args.deepmirtar_neg)
    if args.only_deepmirtar:
        all_data=all_deepmirtar
    else:
        logging.info('loading mirRTarBase data')
        pos_mirtarbase=load_mirtarbase(args.mirtarbase_pos,v22_miRNA2seq)

        seed_27,seed_38=seed_27_38(v22_miRNA2seq)
        all_pos=pd.concat((pos_mirtarbase,all_deepmirtar.loc[all_deepmirtar['label']==1]))
        logging.info('total positive data before filter: %s from deepmirtar %d from miRTarBase %d'%(str(len(all_pos)),len(all_deepmirtar.loc[all_deepmirtar['label']==1]),len(pos_mirtarbase)))

        logging.info('filtering data with miRanda')
        if args.filter==True:
            all_pos=filter_pos(all_pos)
            logging.info('total positive data after filtering: '+str(len(all_pos)))


        random.seed(args.seed)
        logging.info('generating mock data')
        mock_mir=gen_mock_mir(all_pos,seed_27,seed_38,int(args.mock_threshold))
        all_neg=gen_mock_data(all_pos,mock_mir)
        

        all_data=pd.concat([all_pos,all_neg]).reset_index(drop=True)
        all_data=all_data.loc[(all_data['site_idx']+all_data['mir_idx']).apply(tuple, 1).drop_duplicates().index].reset_index(drop=True)
    
        
    all_data.to_pickle(args.output_file)
    logging.info('*** total data: %d  ***'%(len(all_data)))

if __name__ == '__main__':
    main()





