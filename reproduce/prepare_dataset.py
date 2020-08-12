import numpy as np
import pandas as pd
from torch.utils.data import DataLoader,Dataset
from torch.utils.data.sampler import SubsetRandomSampler
from sklearn.model_selection import train_test_split
import argparse
import torch
import pickle
import os
from pathlib import Path



class ClassifyDataset(Dataset):
    def __init__(self,X,site_max_length,mir_max_length,label):
        # padded or slice to max sequence length
        self.site_max_length=site_max_length
        self.mir_max_length=mir_max_length
        
        site_padded=[]
        site_mask=[]
        mir_padded=[]
        mir_mask=[]
        # idx [0,1,2,3,4] -> [pad,A,T,C,G]
        for site,mir in X:
            site_padded.append(np.hstack((np.array(site[:site_max_length])+1,np.zeros(site_max_length-len(site[:site_max_length])))))
            if len(site)>=site_max_length:
                site_mask.append([1]*site_max_length)
            else:
                site_mask.append([1]*len(site)+[0]*(site_max_length-len(site)))
            mir_padded.append(np.hstack((np.array(mir[:mir_max_length])+1,np.zeros(mir_max_length-len(mir[:mir_max_length])))))
            if len(mir)>=mir_max_length:
                mir_mask.append([1]*mir_max_length)
            else:
                mir_mask.append([1]*len(mir)+[0]*(mir_max_length-len(mir)))
    
        self.label=np.array(label,dtype=np.float)
        self.site_seq=np.array(site_padded)
        self.mir_seq=np.array(mir_padded)
        self.mir_mask=np.array(mir_mask)
        self.site_mask=np.array(site_mask)

    def __len__(self):
        return len(self.site_seq)
    def __getitem__(self,idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        data=(self.site_seq[idx],self.mir_seq[idx],self.site_mask[idx],self.mir_mask[idx])
        label=self.label[idx]
        return data,label
    
    def get_batch(self,idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        data=(torch.tensor(self.site_seq[idx]),torch.tensor(self.mir_seq[idx]),torch.tensor(self.site_mask[idx]),torch.tensor(self.mir_mask[idx])) 
        label=torch.tensor(self.label[idx])
        return data,label


def create_dataset(X,y,max_len_site,max_len_mirna,out_file):
    dataset = ClassifyDataset(X,max_len_site,max_len_mirna,y)
    with open(out_file, 'wb') as f:
        pickle.dump(dataset, f)



def p_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("all_data",help="pickle file pandas dataframe which consist of all data(created by genrate_data.py)")
    parser.add_argument("output_dir",help="output dir for dataset",type=Path)
    parser.add_argument("out_comparetools",help="generate pickle dataframe data for comparetools")
    parser.add_argument("--max_len_site",default=50,help="max length for site sequence")
    parser.add_argument("--max_len_mirna",default=30,help="max length for mirna sequence")
    parser.add_argument("--valid_split",default=0.1)
    parser.add_argument("--test_split",default=0.1)
    args = parser.parse_args()

    return args






def main():
    args=p_args()
    all_data=pd.read_pickle(args.all_data)

    X=all_data.loc[:,['site_idx','mir_idx']].values
    y=all_data['label']
    idx=range(len(all_data))

    split_train=args.valid_split+args.test_split
    X_train, X_test, y_train, y_test,idx_train,idx_test = train_test_split(X, y,idx, test_size=split_train,random_state=9487,stratify=y)
    split_test=args.test_split/split_train
    X_test, X_val, y_test, y_val,idx_test,idx_val = train_test_split(X_test, y_test,idx_test ,test_size=split_test,random_state=9487,stratify=y_test)
    print('*** created dataset train: %d valid: %d test: %d ***'%(len(X_train),len(X_val),len(X_test)))
    create_dataset(X_train,y_train,args.max_len_site,args.max_len_mirna, args.output_dir / 'train.pkl')
    create_dataset(X_val,y_val,args.max_len_site,args.max_len_mirna, args.output_dir / 'valid.pkl')
    create_dataset(X_test,y_test,args.max_len_site,args.max_len_mirna, args.output_dir / 'test.pkl')

    with open(args.out_comparetools,'wb') as w:

        pickle.dump(all_data.loc[idx_test].sort_values('miR_ID').reset_index(drop=True),w)
if __name__ == '__main__':
    main()


    
