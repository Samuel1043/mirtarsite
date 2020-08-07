
def load_miRNA2seq(file_mature):
    '''
    input miRbase mature file (mature.fa)
    return a dict of miRNA(key) to sequence(value) which is 5'->3' sequence
    '''
    miRNA2seq={}
    miRNA=""
    with open(file_mature,'r') as r:
        line=r.readline()
        while line:
            line=line.rstrip().lstrip()
            if(line.startswith('>')):
                miRNA=line.split(' ')[0][1:]
                miRNA=miRNA.replace('-mir-','-miR-')
                miRNA2seq[miRNA]=""
            else:
                miRNA2seq[miRNA]+=line
            line=r.readline()
    
    return miRNA2seq

def add_mock_id(mirnas):
    '''
    add mock id for site level microRNA data 
    '''
    tmp=[]
    cnt=0
    prev_id=''
    for mirna in mirnas:
        if(mirna!=prev_id):
            cnt=0
        mirid=mirna+'_'+str(cnt)
        tmp.append(mirid)
        prev_id=mirna
        cnt+=1
    return tmp
def add_site_id(sites):
    '''
    add site id for site level site data
    '''
    tmp=[]
    cnt=0
    prev_id=''
    for site in sites:
        if(site!=prev_id):
            cnt=0
        siteid=site+'_site_'+str(cnt)
        tmp.append(siteid)
        prev_id=site
        cnt+=1
    return tmp
def rename_utrseq(x):
    '''
    make sequence identical (ATCG)
    '''
    tmp=''
    for seq in x:
        if(seq=='a'):
            tmp+='A'
        elif(seq=='t' or seq=='u' or seq=='U'):
            tmp+='T'
        elif(seq=='c'):
            tmp+='C'
        elif(seq=='g'):
            tmp+='G'
        elif(seq=='A' or seq=='T' or seq=='C' or seq=='G'):
            tmp+=seq
            continue
        else:
            print('err',seq)
    return tmp

def rename_mirseq(x):
    '''
    make sequence identical (AUCG)
    '''
    tmp=''
    for seq in x:
        if(seq=='a'):
            tmp+='A'
        elif(seq=='t' or seq=='u' or seq=='T'):
            tmp+='U'
        elif(seq=='c'):
            tmp+='C'
        elif(seq=='g'):
            tmp+='G'
        elif(seq=='A' or seq=='U' or seq=='C' or seq=='G'):
            tmp+=seq
            continue
        else:
            print('err',seq)
    return tmp

def seq2index(seq):
    '''
    convert sequence to index
    '''
    tmp=[]
    for i in seq:
        if(i=='A'):
            tmp.append(0)
        elif(i=='U' or i=='T'):
            tmp.append(1)
        elif(i=='C'):
            tmp.append(2)
        elif(i=='G'):
            tmp.append(3)
        else:
             raise Exception()
    return tmp