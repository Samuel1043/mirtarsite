# MirTarSite
end-to-end microRNA target mrna site prediction

A tool for predicting microRNA target site, given a pair of microRNA and mRNA target sites predict if it's a functional effective miRNA:mRNA interaction.



## Predict custom data

### Requirements
#### Dependencies
* Python 3.7+
* argparse
* torch
* numpy

#### Input 
* microRNA fasta file
```
>hsa-miR-183-5p
UAUGGCACUGGUAGAAUUCACU
>hsa-miR-33a-5p
GUGCAUUGUAGUUGCAUUGCA
...
...
```
* mRNA target site fasta file
```
>AANAT_site_0
TTACCGTGCTGATTACTGTGCTA
>ABCA1_site_9
AACGTAACTTAACGTAACGTAACTTATCATAGTCATGT
>ABCB1_site_1
GAAAAACGATGAAGATATTAGAA
```
* interaction pair file
```
hsa-miR-183-5p  AANAT_site_0    
hsa-miR-33a-5p  ABCA1_site_9    
hsa-miR-129-5p  ABCB1_site_1    
hsa-miR-451a    ABCB1_site_2 
hsa-miR-129-5p  ABCB1_site_4
...
...
```
* pretrain state dict file


#### Usage
```
usage: custom_predict.py [-h] [--batch_size BATCH_SIZE]
                         [--model_type MODEL_TYPE]
                         [--embedding_hidden EMBEDDING_HIDDEN]
                         [--rnn_hidden RNN_HIDDEN] [--rnn_layer RNN_LAYER]
                         [--class_dropout CLASS_DROPOUT]
                         [--threshold THRESHOLD]
                         state_dict_file miRNA_file site_file interaction_pair
                         output_file

mirtarsite
end-to-end microRNA target mrna site prediction tool

Input: microRNA,mRNA target site fasta file,pretrained state dict file
(Default hyperparameter seeting is the same as my pretrained state dict file)

positional arguments:
  state_dict_file       state dict for model
  miRNA_file            microRNA fasta file 5'->3' sequence
  site_file             mRNA target site fasta file 3'->5' sequence
  interaction_pair      interaction pairs for the given microRNA and mRNA
                        target site fasta file
  output_file           output predction file

optional arguments:
  -h, --help            show this help message and exit
  --batch_size BATCH_SIZE
  --model_type MODEL_TYPE
  --embedding_hidden EMBEDDING_HIDDEN
                        size of hidden in embedding
  --rnn_hidden RNN_HIDDEN
                        size of hidden in RNN layer
  --rnn_layer RNN_LAYER
                        num of layer for RNN
  --class_dropout CLASS_DROPOUT
                        dropout in classify layer
  --threshold THRESHOLD
                        threshold for binary classification

```

#### example command
```
python3 custom_predict.py ./state_dict/b16_lr0.001_embd100_rnnlayer1_rnnhidden100_drop0.3_ep47.pth ./example/all_mir.fa ./example/all_site.fa ./example/interaction_pair.txt ./out.txt
```


## Thesis Reproduce

```
cd reproduce
```
### download require data
* require gdrive
```
sh download.sh
```
* manually download
https://drive.google.com/drive/folders/1jOEVfwPUoNPLQ6V26bJLKNMVGq0i6HtR?usp=sharing
putting as the below file structure
```
download   
└───deepmirtar...
└───mirtarsite...
└───miRBase...
└───miRanda-aug2010.tar.gz
```

### set miRanda to PATH
```
mkdir comparetools
tar -zxvf download/miRanda-aug2010.tar.gz -C comparetools/
export PATH=$PATH:$(pwd)/comparetools/miRanda-3.3a/bin
```

### V1 dataset
```
python3 generate_data.py ./download/deepmirtar/Supplementary\ File\ S1.\ Positive\ dataset.xlsx ./download/deepmirtar/Supplementary\ File\ S2.\ Negative\ dataset.xlsx ./download/mirtarbase/MicroRNA_Target_Sites.xlsx ./download/miRBase/v22_mature.fa  ./deepmirtar_filter.pkl --filter True --only_deepmirtar True

mkdir dataset_deepmirtar

python3 prepare_dataset.py ./deepmirtar_filter.pkl ./dataset_deepmirtar/ ./dataset_deepmirtar/deepmirtar_filter_test.pkl

mkdir thesis_hyper_step2

# train for five kinds of dropout (step2 hyperparameter tuning)
declare -a dropout_arr=(0 0.1 0.2 0.3 0.4 0.5)
for i in "${dropout_arr[@]}"; do
    python3 train_model.py --lr 0.001 --embedding_hidden 100 --rnn_hidden 400 --rnn_layer 1 --train_epoch 50 --batch_size 64 --class_dropout $i ./dataset_deepmirtar './thesis_hyper_step2/deepmirtar_b64_lr0.001_embd100_rnnlayer1_rnnhidden400_drop'$i
done

# evaluate on the best validation accuracy state dict
python3 evaluate.py ./thesis_hyper_step2/deepmirtar_b64_lr0.001_embd100_rnnlayer1_rnnhidden400_drop0/classify9.pth ./dataset_deepmirtar/valid.pkl --embedding_hidden 100 --rnn_hidden 400 --rnn_layer 1 --class_dropout 0 --threshold 0.7
```


### V2 dataset

```
python3 generate_data.py ./download/deepmirtar/Supplementary\ File\ S1.\ Positive\ dataset.xlsx ./download/deepmirtar/Supplementary\ File\ S2.\ Negative\ dataset.xlsx ./download/mirtarbase/MicroRNA_Target_Sites.xlsx ./download/miRBase/v22_mature.fa  ./all_data_filter.pkl --filter True


mkdir ./all_data

python3 prepare_dataset.py ./all_data_filter.pkl ./all_data/ ./all_data/all_data_filter_test.pkl

mkdir thesis_hyper_step2


# train for five kinds of dropout 
declare -a dropout_arr=(0 0.1 0.2 0.3 0.4 0.5)
for i in "${dropout_arr[@]}"; do
    python3 train_model.py --lr 0.001 --embedding_hidden 100 --rnn_hidden 100 --rnn_layer 1 --train_epoch 50 --batch_size 16 --class_dropout $i ./all_data './thesis_hyper_step2/b16_lr0.001_embd100_rnnlayer1_rnnhidden100_drop'$i;
done

# evaluate on the best validation accuracy state dict
python3 evaluate.py ./thesis_hyper_step2/b16_lr0.001_embd100_rnnlayer1_rnnhidden100_drop0.3/classify47.pth ./all_data/valid.pkl --embedding_hidden 100 --rnn_hidden 100 --rnn_layer 1 --class_dropout 0.3 --threshold 0.3

```

### compare with other prediction tools
#### require tools
* PITA
* RNAhybrid
* miRanda
```
python3 gen_compare.py  ./dataset_deepmirtar/deepmirtar_filter_test.pkl ../comparetools/pita/thesis_v2_deepmirtar ../comparetools/miRanda-3.3a/thesis_v2_deepmirtar ../comparetools/RNAhybrid-2.1.2/thesis_v2_deepmirtar

miranda ./thesis_v2_deepmirtar/all_mir.fa ./thesis_v2_deepmirtar/all_site.fa -sc 100 -en -18 -restrict ./thesis_v2_deepmirtar/interaction_pair.txt > ./thesis_v2_deepmirtar_out/en-18.txt

bash run_rnahybrid.sh ./thesis_v2_deepmirtar/ ./thesis_v2_deepmirtar_out/

bash run_pita.sh  ./thesis_v2_deepmirtar/ ./thesis_v2_deepmirtar_out/

python3 comparetools.py ../comparetools/miRanda-3.3a/thesis_v2_deepmirtar_out/en-18.txt ../comparetools/pita/thesis_v2_deepmirtar_out/ ../comparetools/pita/thesis_v2_deepmirtar/ ../comparetools/RNAhybrid-2.1.2/thesis_v2_deepmirtar_out ../comparetools/RNAhybrid-2.1.2/thesis_v2_deepmirtar
```

