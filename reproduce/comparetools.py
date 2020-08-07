import subprocess
import glob
import numpy as np
from sklearn.metrics import f1_score
import logging
import pandas as pd
import argparse
from pathlib import Path

def input_miRanda(interactionpair,save_dir,mrna_interact=False):
    '''
    input pandas dataframe interaction pair
    generate input file for miRanda 
    '''
    subprocess.call(["mkdir",save_dir])
    mir=interactionpair.loc[:,['miR_ID','miRNA_seq']].drop_duplicates()
    site=interactionpair.loc[:,['mRNA_ID','site_seq']].drop_duplicates()
    mrna=interactionpair.loc[:,['mRNA_ID','mRNA Seq']].drop_duplicates()
    
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
            with open(save_dir / save_str,'w') as w1:
                w1.write(word)
            try:
                f.close()
                f2.close()
            except: 
                pass
            f = open(save_dir+prev_mireid+'_site.fa', 'w')
            f2 =open(save_dir+prev_mireid+'_utr.fa','w')
        word=''
        siteid=e['mRNA_ID']
        siteseq=e["site_seq"][::-1]
        word+='>'+siteid+'_site\n'
        word+=siteseq+'\n'
        f.write(word) 
        
        word=''
        utrid=e['mRNA_ID']
        utrseq=e["mRNA Seq"][::-1]
        word+='>'+utrid+'_utr\n'
        word+=utrseq+'\n'
        f2.write(word) 
    f2.close()
    f.close()




# calculate performance

def cal_lines(file):
    p = subprocess.Popen(["wc","-l",file], stdout=subprocess.PIPE)
    lines=p.communicate()[0].decode('ascii').split(' ')[0]
    return int(lines)

#RNAhybrid
def calculate_rnahybrid_result(rnahybrid_result_dir,rnahybrid_site_input_dir):
    def parse_result(result_file):
        count=0
        with open(result_file,'r') as r:
            prev_site_id=''
            line=r.readline()
            while line:
                cur_site_id=line.split(':')[0]
                if(prev_site_id==cur_site_id):
                    pass
                else:
                    count+=1
                prev_site_id=cur_site_id
                line=r.readline()
        return count


    negative_total=0
    negative_correct=0
    positive_total=0
    positive_correct=0

    def sort_value(x):
        # print(x.split('/')[-1].split('_')[0])
        try:
            s2=x.split('/')[-1].split('_')[2]
        except:
            s2='a'
        return (x.split('/')[-1].split('_')[0],s2)

    interaction_files=sorted(glob.glob(str(rnahybrid_site_input_dir / '*_site.fa')))
    result_files=sorted(glob.glob(str(rnahybrid_result_dir / '*')))
    assert len(interaction_files)==len(result_files)
    for result_file,interact_file in zip(result_files,interaction_files):
        if result_file.split('/')[-1].split('.')[0] in interact_file:
            # lines=!wc -l $interact_file
            lines=cal_lines(interact_file)
            site_lines=lines//2
            
            
            #negative data
            if 'mock' in interact_file:
                
                negative_total+=site_lines
                negative_correct+=site_lines-parse_result(result_file)
            #positive data
            else:
                positive_total+=site_lines
                positive_correct+=parse_result(result_file)
        else:
            pass
            # print('err',result_file,interact_file)
    return positive_total,positive_correct,negative_total,negative_correct


import glob
def calculate_pita_result(pita_result_dir,pita_site_input_dir):
    negative_total=0
    negative_correct=0
    positive_total=0
    positive_correct=0
    

    interaction_files=sorted(glob.glob(str(pita_site_input_dir / '*_site.fa')))
    result_files=sorted(glob.glob(str(pita_result_dir / '*_targets.tab')))
    assert len(interaction_files)==len(result_files)

    for result_file,interact_file in zip(result_files,interaction_files):
        if result_file.split('/')[-1].split('_pita_results_targets.tab')[0] in interact_file:
            lines=cal_lines(interact_file)            
            site_lines=lines//2
            

            lines=cal_lines(result_file)
            result_lines=lines-1
            
            #negative data
            if 'mock' in interact_file:
                negative_total+=site_lines
                negative_correct+=site_lines-result_lines
            #positive data
            else:
                positive_total+=site_lines
                positive_correct+=result_lines
    return positive_total,positive_correct,negative_total,negative_correct



def calculate_miranda_performance(miranda_result):
    def parse_predict(line):
        if 'No Hits Found above Threshold' in line:
            return 0
        elif 'Forward:' in line:
            return 1
        else:
            raise Exception('miranda input file error')

    label=[]
    predict=[]
    with open(miranda_result,'r') as r:
        line=r.readline()
        while line:
            if(line.startswith('Performing Scan:')):
                #negative data
                if 'mock' in line:
                    label.append(0)
                    r.readline()
                    r.readline()
                    line=r.readline()
                    predict.append(parse_predict(line))
                #positive data
                else:
                    label.append(1)
                    r.readline()
                    r.readline()
                    line=r.readline()
                    predict.append(parse_predict(line))

            line=r.readline()
    return np.array(label),np.array(predict)

def score_matrix(ptotal,pcorrect,ntotal,ncorrect):
    '''
    reutrn (acc,tpr==precision,tnr==Specificity,prec==ppr)
    '''
    return (pcorrect+ncorrect)/(ptotal+ntotal)\
,pcorrect/ptotal,ncorrect/ntotal,pcorrect/(pcorrect+(ntotal-ncorrect))
    
def p_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("miranda_result",help="miranda result to calculate performance")
    parser.add_argument('pita_result',type=Path,help='pita result dir to calculate performance')
    parser.add_argument('pita_input',type=Path,help='pita original input dir')
    parser.add_argument('rnahybrid_result',type=Path,help='rnahybrid result dir to calculate performance')
    parser.add_argument('rnahybrid_input',type=Path,help='rnahybrid original input dir')
    args=parser.parse_args()

    return args

def main():

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(level=logging.INFO,format=FORMAT)
    args=p_args()

    logging.info('calculating miRanda result')
    label,predict=calculate_miranda_performance(args.miranda_result)
    miranda_acc=np.sum(label==predict)/len(predict)
    miranda_tpr=np.sum(label[label==1]==predict[label==1])/len(label[label==1])
    miranda_tnr=np.sum(label[label==0]==predict[label==0])/len(label[label==0])
    miranda_f1=f1_score(label,predict)



    logging.info('calculating RNAhybrid result')
    rnahybrid_positive_total,rnahybrid_positive_correct,rnahybrid_negative_total,rnahybrid_negative_correct=calculate_rnahybrid_result(args.rnahybrid_result,args.rnahybrid_input)
    print(rnahybrid_positive_total,rnahybrid_positive_correct,rnahybrid_negative_total,rnahybrid_negative_correct)
    rnahybrid_acc,rnahybrid_tpr,rnahybrid_tnr,rnahybrid_prec=score_matrix(rnahybrid_positive_total,rnahybrid_positive_correct,rnahybrid_negative_total,rnahybrid_negative_correct)
    rnahybrid_f1=2*(rnahybrid_prec*rnahybrid_tpr)/(rnahybrid_prec+rnahybrid_tpr)



    logging.info('calculating PITA result')
    pita_positive_total,pita_positive_correct,pita_negative_total,pita_negative_correct=calculate_pita_result(args.pita_result,args.pita_input)
    print(pita_positive_total,pita_positive_correct,pita_negative_total,pita_negative_correct)
    pita_acc,pita_tpr,pita_tnr,pita_prec=score_matrix(pita_positive_total,pita_positive_correct,pita_negative_total,pita_negative_correct)
    pita_f1=2*(pita_prec*pita_tpr)/(pita_prec+pita_tpr)



    data=[[miranda_acc,rnahybrid_acc,pita_acc],\
    [miranda_tpr,rnahybrid_tpr,pita_tpr],\
    [miranda_tnr,rnahybrid_tnr,pita_tnr],\
        [miranda_f1,rnahybrid_f1,pita_f1]]

    columns=['miranda','rnahybrid','pita']
    print(pd.DataFrame(data,columns=columns,index=['acc','tpr','tnr','f1']))

if __name__=='__main__':
    main()