import argparse
import torch
import numpy as np
from sklearn.metrics import f1_score
import pickle
from torch.utils.data import DataLoader
from prepare_dataset import ClassifyDataset
from pathlib import Path

def calculate_tpr_fpr(model,testLoader,thresholds):
    '''
    calculate tprs and fprs for different threshold to plot auc curve
    '''
    tprs=[]
    fprs=[]
    f1s=[]
    accs=[]
    for threshold in thresholds:
        pred,target=evaluate_model(model,testLoader,threshold)
        accs.append(np.sum(pred==target)/len(pred))
        #tpr tp/tp+fn
        tpr=np.sum(target[np.where(target==1)[0]]==pred[np.where(target==1)[0]])/len(np.where(target==1)[0])
        #fpr fp/fp+tn
        fpr= np.sum(target[np.where(target==0)[0]]!=pred[np.where(target==0)[0]])/len(np.where(target==0)[0])
        f1s.append(f1_score(target,pred))
        tprs.append(tpr)
        fprs.append(fpr)
    return tprs,fprs,f1s,accs


def evaluate_model(model,testLoader,threshold):
    model.eval()
    with torch.no_grad():
        pred=torch.DoubleTensor([])
        target=torch.FloatTensor([])
        for idx,seq in enumerate(testLoader):
            out=model(seq[0])
            prob=torch.sigmoid(out).detach().cpu()
            predict=prob
            predict[predict>threshold]=1
            predict[predict<=threshold]=0
            pred=torch.cat((pred,predict),0)
            tar=seq[1].detach().cpu()
            target=torch.cat((target,tar),0)

    return target.numpy(),pred.numpy()

def load_dataset(pred_file,batch_size):
    with open(pred_file,'rb') as r:
        testData=pickle.load(r)
    ClassifytestData =DataLoader(testData ,batch_size=batch_size,num_workers=4)
    return ClassifytestData

def calculate_performance(target,pred):
    model_acc=np.sum(target==pred)/len(pred)
    model_tpr=np.sum(target[target==1]==pred[target==1])/len(target[target==1])
    model_tnr=np.sum(target[target==0]==pred[target==0])/len(target[target==0])
    model_f1=f1_score(target,pred)

    print("total testing data %d | accuracy %.4f | TPR %.4f | TNR %.4f | F1 %.4f "%(len(target),model_acc,model_tpr,model_tnr,model_f1))

def p_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("state_dict_file",help="state dict for model")
    parser.add_argument('pred_file')
    parser.add_argument('--batch_size',type=int,default=100)
    parser.add_argument('--model_type',default='similarity matrix')
    parser.add_argument('--embedding_hidden',type=int,default=100,help='size of hidden in embedding')
    parser.add_argument('--rnn_hidden',type=int,default=100,help='size of hidden in RNN layer')
    parser.add_argument('--rnn_layer',type=int,default=1,help='num of layer for RNN')
    parser.add_argument('--class_dropout',type=float,default=0,help='dropout in classify layer')
    parser.add_argument('--threshold',type=float,default=0.5,help='threshold for binary classification')

    args=parser.parse_args()
    return args


def main():
    
    args=p_args()
    if args.model_type=='similarity matrix':
        from models import SimilarityMatrixMask
        modelclassify=SimilarityMatrixMask(5,args.rnn_hidden,args.embedding_hidden,args.rnn_layer,args.class_dropout).double().cuda()

    
    modelclassify.load_state_dict(torch.load(args.state_dict_file))
    testLoader=load_dataset(args.pred_file,args.batch_size)

    target,pred=evaluate_model(modelclassify,testLoader,args.threshold)


    calculate_performance(target,pred)

if __name__=='__main__':
    main()
