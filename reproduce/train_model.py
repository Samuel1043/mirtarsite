import subprocess
import time
import torch
import numpy as np 
import argparse
from pathlib import Path
from prepare_dataset import ClassifyDataset
import pickle
from torch.utils.data import DataLoader
from models import SimilarityMatrixMask
import torch.optim as optim
import torch.nn as nn
import logging

def train_classify(model,optimizer,criterion,num_epoch,trainLoader,valLoader,size_train,size_val,state_dict_path=None,device='cpu'):
    '''
    train model
    '''
    if state_dict_path:
        subprocess.call(["mkdir",state_dict_path])
        f=open(state_dict_path / 'train_epoch.txt','w')
    model=model
    
    best_acc=0
    val_acc=0
    for epoch in range(1,num_epoch+1):
        epoch_start_time=time.time()
        train_loss=0
        right=0
        model.train()
        for idx,seq in enumerate(trainLoader):
            out=model(seq[0])
            pred=torch.sigmoid(out)
            pred[pred>=0.5]=1
            pred[pred<0.5]=0

            right+=np.sum(pred.detach().cpu().numpy()==seq[1].numpy())
            
            batch_loss=criterion(out,seq[1].to(device))
            optimizer.zero_grad()
            batch_loss.backward()
            optimizer.step()

            train_loss+=batch_loss.item()
            
            progress = ('#' * int(float(idx)/len(trainLoader)*40)).ljust(40)
            print('[%03d|%03d] %2.2f sec(s) | %s |' %(epoch,num_epoch,time.time()-epoch_start_time,progress),end='\r',flush=True)
        
        right_val=0
        val_loss=0
        model.eval()
        with torch.no_grad():
            for idx,seq in enumerate(valLoader):
                out=model(seq[0])
                
                pred=torch.sigmoid(out)
                pred[pred>=0.5]=1
                pred[pred<0.5]=0
                right_val+=np.sum(pred.detach().cpu().numpy()==seq[1].numpy())

                batch_loss=criterion(out,seq[1].to(device))

                val_loss+=batch_loss.item()
                progress = ('#' * int(float(idx)/len(valLoader)*40)).ljust(40)
                print('[%03d|%03d] %2.2f sec(s) | %s |' %(epoch,num_epoch,time.time()-epoch_start_time,progress),end='\r',flush=True)
        

        if state_dict_path:
            file='classify'+str(epoch)+'.pth'
            torch.save(model.state_dict(),state_dict_path / file)
            f.write('[%03d|%03d] %2.2f sec(s) | train loss: %2.5f  | train acc: %2.5f | val loss: %2.5f | val acc: %2.5f\n' \
              %(epoch,num_epoch,time.time()-epoch_start_time,train_loss/len(trainLoader),right/size_train,val_loss/len(valLoader),right_val/size_val))
        print('[%03d|%03d] %2.2f sec(s) | train loss: %2.5f  | train acc: %2.5f | val loss: %2.5f | val acc: %2.5f' \
              %(epoch,num_epoch,time.time()-epoch_start_time,train_loss/len(trainLoader),right/size_train,val_loss/len(valLoader),right_val/size_val))
    if state_dict_path:
        f.close()
    return model

def load_dataset(dataset_dir,batch_size):
    with open(dataset_dir / 'train.pkl','rb') as r:
        trainData=pickle.load(r)
    with open(dataset_dir / 'valid.pkl','rb') as r:
        valData=pickle.load(r)
    
    ClassifytrainData = DataLoader(trainData, batch_size=batch_size,num_workers=4)
    ClassifyvalData =DataLoader(valData ,batch_size=batch_size,num_workers=4)
    return ClassifytrainData,ClassifyvalData

def p_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset_dir",type=Path,help="dataset dir")
    parser.add_argument("state_dict_path",help="path for storing state dict",type=Path)
    parser.add_argument("--lr",type=float,default=0.001,help='learning rate for model')
    parser.add_argument('--embedding_hidden',type=int,default=100,help='size of hidden in embedding')
    parser.add_argument('--rnn_hidden',type=int,default=100,help='size of hidden in RNN layer')
    parser.add_argument('--rnn_layer',type=int,default=1,help='num of layer for RNN')
    parser.add_argument('--class_dropout',type=float,default=0,help='dropout in classify layer')
    parser.add_argument('--train_epochs',type=int,default=30,help="num of epoch to train")
    parser.add_argument('--batch_size',type=int,default=32,help="batch size for training model")
    parser.add_argument("--model_type",default='similarity matrix')
    parser.add_argument('--cuda',default=False)
    args=parser.parse_args()
    return args

def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
def main():
    set_seed(9487)
    args=p_args()
    logging.basicConfig(format='%(asctime)s | %(levelname)s | %(message)s',
                        level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
                        
    logging.info('loading dataset ...')
    trainLoader,valLoader=load_dataset(args.dataset_dir,args.batch_size)
    logging.info('train data size '+str(len(trainLoader.dataset)))
    logging.info('val data size '+str(len(valLoader.dataset)))

    if args.cuda:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    else:
        device='cpu'

    modelclassify=SimilarityMatrixMask(5,args.rnn_hidden,args.embedding_hidden,args.rnn_layer,args.class_dropout,device=device).double().to(device)
    optimizer=optim.Adam(params=modelclassify.parameters(),lr=args.lr)
    criterion=nn.BCEWithLogitsLoss()
    logging.info('training model ...')
    train_classify(modelclassify,optimizer,criterion,args.train_epochs,trainLoader,valLoader,\
        len(trainLoader.dataset),len(valLoader.dataset),args.state_dict_path,device=device)

    
if __name__=='__main__':
    main()

