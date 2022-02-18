import sys
import numpy as np
import pandas as pd


def hmp2num(gt):
    nmaize=len(gt.columns)-12
    nmarker=len(gt)
    
    # keep only biallelic sites
    gt['alleleN']=[len(x) for x in gt.alleles]
    gt=gt.loc[gt.alleleN==3].copy()


    # recode homo-ref to 0 and homo-alt to 1 and others to np.nan 
    for col in range(11,len(gt.columns)-1):
        gt.iloc[:,col]=[-1 if x==a[0] else 1 if x==a[2] else 0 for (x,a) in zip(gt.iloc[:,col],gt.alleles)]
    
    # filter out irrelevant columns
    gt=gt.iloc[:,np.r_[0,11:len(gt.columns)-1]]
    gt=gt.set_index('rs#')
    
    # maize inbred lines as rows, genotypic calls as columns
    gt=gt.T
    
    # drop any SNP sites with nan
    gt=gt.dropna(axis='columns').copy()

    return gt


def convertname(gt, namedict):
    namedict = {SNP:field for SNP,field in zip(namedict.SNPname,namedict.fieldname)}
    for i in range(len(gt)):
        name = gt.iloc[i,0]
        if name not in namedict: continue
        gt.iloc[i,0] = namedict[gt.iloc[i,0]]
    gt = gt[-pd.isnull(gt.iloc[:,0])]
    return gt


if __name__ == "__main__":
    
    gt_in=sys.argv[1]
    gt_out=sys.argv[2]
    
    gt=pd.read_csv(gt_in,sep="\t")
    gt=hmp2num(gt)
    
#    if len(sys.argv) == 4:
    namedict = sys.argv[3]
    namedict = pd.read_csv(namedict,sep=",")
    gt = convertname(gt,namedict)
    gt.to_csv(gt_out)
