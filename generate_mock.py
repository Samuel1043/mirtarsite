import random
import subprocess
from tqdm import tqdm
import pandas as pd
import argparse
from utils_annotate import load_miRNA2seq,add_mock_id,seq2index


def knuth_shuffle(mir_seq):
    # no extra space
    for i in range(len(mir_seq) - 1, 0, -1):
        p = random.randrange(0, i + 1)
        mir_seq[i], mir_seq[p] = mir_seq[p], mir_seq[i]
    return mir_seq



def miranda_filter(mir_seq,site_seq):
    mock_able=0
    with open('./tmp_mir.fa','w') as w:
        w.write('>tmp_mir\n'+mir_seq)
    with open('./tmp_site.fa','w') as w:
        w.write('>tmp_site\n'+site_seq)
    subprocess.call(["miranda", "./tmp_mir.fa", "./tmp_site.fa","-sc","100","-out","./tmp_out.txt"])
#     ! ../comparetools/miRanda-3.3a/bin/miranda ./tmp_mir.fa ./tmp_site.fa -sc 100 -out ./tmp_out.txt
    with open('./tmp_out.txt','r') as r:
        line=r.readline()
        while line:
            if line.startswith("No Hits Found"):
                mock_able=0
                break
            if line.startswith("   Forward:"):
                mock_able=1
                break
            line=r.readline()
    subprocess.call(["rm",'./tmp_mir.fa','tmp_site.fa','./tmp_out.txt'])
    return mock_able
       
def gen_mock_mir(mi2site_arr,seed_27,seed_38,threshold):
    mock_mir=[]
    mock_fail=0
    for i,e in tqdm(mi2site_arr.iterrows(),total=mi2site_arr.shape[0]):
        mir_seq=e['miRNA_seq']
        r1=mir_seq[1:7]
        r2=mir_seq[2:8]
        mock_able=0
        cnt=0
        while(True):
            cnt+=1
            if cnt>threshold:
                mir=[]
                mock_fail+=1
                break
            mir=knuth_shuffle(list(mir_seq))
            r1=mir[1:7]
            if r1 in seed_27:
                continue
            r2=mir[2:8]
            if r2 in seed_38:
                continue
            mock_able=miranda_filter(''.join(mir),e['site_seq'])
            if mock_able==0:
                continue
            break
        mock_mir.append(''.join(mir))
    print('***generate mock fail*** counts: %d'%(mock_fail))
    return mock_mir

def gen_mock_data(all_pos,mock_mir):
    all_neg=pd.DataFrame({'miRTarBase ID':all_pos['miRTarBase ID'],
                            'miR_ID':all_pos['miR_ID'],'mRNA_ID':all_pos['mRNA_ID']
                            ,'Species (miRNA)':all_pos['Species (miRNA)']
                            ,'site_seq':all_pos['site_seq'],'miRNA_seq':mock_mir})
    all_neg=all_neg.sort_values('miR_ID')
    all_neg['miR_ID']=all_neg['miR_ID'].map(lambda x:x+'_mock')
    all_neg['miR_ID']=add_mock_id(all_neg['miR_ID'])

    all_neg['miRNA_seq_len']=all_neg['miRNA_seq'].map(len)
    all_neg['site_seq_len']=all_neg['site_seq'].map(len)
    all_neg['mir_idx']=all_neg['miRNA_seq'].map(seq2index)
    all_neg['site_idx']=all_neg['site_seq'].map(seq2index)
    all_neg['label']=[0]*len(all_neg['site_idx'])
    return all_neg
       

def p_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("v22_mirna",help="mirBase v22 mature microRNA file (http://www.mirbase.org/ftp.shtml)")
    parser.add_argument("positive_data_pickle",help="interaction pair of positive data in pickle file format")
    parser.add_argument("output_file",help="output whole dataset in pickle format after adding mock data")
    parser.add_argument("--mock_threshold",default=1000, help="maximun shuffle times for mock data generation")
    
    args = parser.parse_args()
    return args


def main():
    args=p_args()
    
    v22_miRNA2seq=load_miRNA2seq(args.v22_mirna)

    seed_27=[]
    seed_38=[]
    for i in v22_miRNA2seq.values():
        seed_27.append(i[1:7])
        seed_38.append(i[2:8])

    all_pos=pd.read_pickle(args.positive_data_pickle)

    mock_mir=gen_mock_mir(all_pos,seed_27,seed_38,int(args.mock_threshold))
    all_neg=gen_mock_data(all_pos,mock_mir)


    all_data=pd.concat([all_pos,all_neg]).reset_index(drop=True)
    all_data.to_pickle(args.output_file)

if __name__ =='__main__':
    main()
