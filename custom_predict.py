import argparse
import torch
import numpy as np
import pickle
from torch.utils.data import DataLoader,TensorDataset
from pathlib import Path


def parse_fasta(fa_file):
    fa_dict={}
    with open(fa_file,'r') as r:
        line=r.readline()
        while line:
            if(line.startswith('>')):
                fa_seq=''
                fa_id=line[1:].split(' ')[0].strip()
                fa_dict[fa_id]=''
            else:
                fa_seq+=line.strip()
                fa_dict[fa_id]=fa_seq
            line=r.readline()
    return fa_dict

def parse2numpy(interaction_pair,mir_dict,site_dict):
    # idx [0,1,2,3,4] -> [pad,A,T,C,G]
    max_mir_len=30
    max_site_len=50

    one_hot_dict={'a':1,'t':2,'c':3,'g':4,'u':2}
    mir_input=np.zeros((len(interaction_pair),max_mir_len))
    site_input=np.zeros((len(interaction_pair),max_site_len))
    mir_mask=np.zeros((len(interaction_pair),max_mir_len))
    site_mask=np.zeros((len(interaction_pair),max_site_len))

    for pair_idx,line in enumerate(interaction_pair):
        line=line.strip()
        mirid,siteid=line.split('\t')
        mir_idx=[one_hot_dict[seq] for seq in mir_dict[mirid].lower()]
        site_idx=[one_hot_dict[seq] for seq in site_dict[siteid].lower()]
        # max length of input to model pad to max len
        mir_idx=mir_idx[:max_mir_len]
        site_idx=site_idx[:max_site_len]
        mir_mask[pair_idx]=[1]*len(mir_idx)+[0]*(max_mir_len-len(mir_idx))
        site_mask[pair_idx]=[1]*len(site_idx)+[0]*(max_site_len-len(site_idx))
        mir_idx=mir_idx+(max_mir_len-len(mir_idx))*[0]
        site_idx=site_idx+(max_site_len-len(site_idx))*[0]

        mir_input[pair_idx]=mir_idx
        site_input[pair_idx]=site_idx
    return mir_input,site_input,mir_mask,site_mask

def predict(model,testLoader,threshold):
    model.eval()
    with torch.no_grad():
        pred=torch.DoubleTensor([])
        
        for idx,seq in enumerate(testLoader):
            out=model(seq)
            prob=torch.sigmoid(out).detach().cpu()
            predict=prob
            predict[predict>threshold]=1
            predict[predict<=threshold]=0
            pred=torch.cat((pred,predict),0)

    return pred.numpy().astype(int)

def output_result(out_file,interaction_pair,preds):
    word=''
    for pair,pred in zip(interaction_pair,preds):
        word+=pair.rstrip()+'\t'+str(pred)+'\n'
    with open(out_file,'w') as w:
        w.write(word[:-1])

def p_args():
    # default is the best hyperparameter
    parser=argparse.ArgumentParser(description='''mirtarsite \nend-to-end microRNA target mrna site prediction tool \n\nInput: microRNA,mRNA target site fasta file,pretrained state dict file \n(Default hyperparameter seeting is the same as my pretrained state dict file)''',\
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("state_dict_file",help="state dict for model")
    parser.add_argument('miRNA_file',help="microRNA fasta file 5'->3' sequence")
    parser.add_argument('site_file',help="mRNA target site fasta file 3'->5' sequence")
    parser.add_argument('interaction_pair',help='interaction pairs for the given microRNA and mRNA target site fasta file')
    parser.add_argument('output_file',help='output predction file')
    parser.add_argument('--batch_size',type=int,default=100)
    parser.add_argument('--model_type',default='similarity matrix')
    parser.add_argument('--embedding_hidden',type=int,default=100,help='size of hidden in embedding')
    parser.add_argument('--rnn_hidden',type=int,default=100,help='size of hidden in RNN layer')
    parser.add_argument('--rnn_layer',type=int,default=1,help='num of layer for RNN')
    parser.add_argument('--class_dropout',type=float,default=0.3,help='dropout in classify layer')
    parser.add_argument('--threshold',type=float,default=0.3,help='threshold for binary classification')
    args=parser.parse_args()
    return args


def main():

    args=p_args()
    mirna_dict=parse_fasta(args.miRNA_file)
    site_dict=parse_fasta(args.site_file)

    with open(args.interaction_pair,'r') as r:
        interaction_pair=r.readlines()

    if args.model_type=='similarity matrix':
        from models import SimilarityMatrixMask
        modelclassify=SimilarityMatrixMask(5,args.rnn_hidden,args.embedding_hidden,args.rnn_layer,args.class_dropout).double().cuda()

    mir_input,site_input,mir_mask,site_mask=parse2numpy(interaction_pair,mirna_dict,site_dict)

    custom_dataset=TensorDataset(torch.tensor(site_input),torch.tensor(mir_input),torch.tensor(site_mask),torch.tensor(mir_mask))

    modelclassify.load_state_dict(torch.load(args.state_dict_file))
    testLoader=DataLoader(custom_dataset,batch_size=args.batch_size,num_workers=4)

    pred=predict(modelclassify,testLoader,args.threshold)
    output_result(args.output_file,interaction_pair,pred)

if __name__=='__main__':
    main()