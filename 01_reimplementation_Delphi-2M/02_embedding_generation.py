'''
This script implements embedding generation on individual-wise dataset
'''
import os
os.chdir("<path_to_file>")
from Delphi.utils import get_batch, get_p2i
import numpy as np
import torch
from Delphi.model import DelphiConfig, Delphi
from tqdm import tqdm
import math
import pandas as pd
from scipy.stats import norm
from collections import Counter

# ================================
# Load Delphi-2M
# ================================
out_dir = 'Delphi/Delphi-2M-real'
device = 'cuda:0' # examples: 'cpu', 'cuda', 'cuda:0', 'mps', etc.
dtype ='float32' #'bfloat16' # 'float32' or 'bfloat16' or 'float16'
dtype = {'float32': torch.float32, 'float64': torch.float64, 'bfloat16': torch.bfloat16, 'float16': torch.float16}[dtype]
seed = 1337

torch.manual_seed(seed)
torch.cuda.manual_seed(seed)

ckpt_path = os.path.join(out_dir, 'ckpt.pt')
checkpoint = torch.load(ckpt_path, map_location=device)
conf = DelphiConfig(**checkpoint['model_args'])
model = Delphi(conf)
state_dict = checkpoint['model']
model.load_state_dict(state_dict)

model.eval()
model = model.to(device)

checkpoint['model_args']

# ================================
# Generate patients' embeddings
# ================================
# load data
val = np.fromfile('Delphi/data/ukb_real_data/all.bin', dtype=np.uint32).reshape(-1,3)
# data statistics
#token_list = val[:,2].tolist() # min: 1, max: 1268
#uniq_list = list(set(token_list)) # 1261
#full_list = [i for i in range(1, 1269)]
#list(set(full_list) - set(uniq_list)) # [935, 781, 1038, 340, 1141, 1046, 1047] 

val_p2i = get_p2i(val) # mapping trajectory id to its position in the dataset

dataset_subset_size = len(val_p2i)

d100k = get_batch(range(dataset_subset_size), val, val_p2i,  
              select='left', block_size=128, 
              device=device, padding='random')



PAD_ID = 0
B  = d100k[0].shape[0]
bs = 256

embs = []
model.to(device).eval()

with torch.no_grad():
    # split each tensor the same way
    splits = [x.split(bs) for x in d100k]
    for dd in tqdm(zip(*splits), total=math.ceil(B/bs)):
        idx = dd[0].to(device)   # tokens [b, T], int
        age = dd[1].to(device)   # ages   [b, T], float

        logits, loss, att, hidden = model(idx, age, return_hiddens=True)  # hidden = x = [b, T, D]

        # one vector per patient
        T = idx.size(1)
        mask = (idx != PAD_ID)                                        # [b, T]
        pos  = torch.arange(T, device=idx.device).unsqueeze(0).expand_as(idx)
        last_ix = (pos * mask).max(dim=1).values.long()               # [b]
        attn = mask.unsqueeze(-1)              # [b, T, 1]
        emb  = (hidden * attn).sum(1) / attn.sum(1).clamp(min=1)

        embs.append(emb.cpu().numpy())

embs = np.vstack(embs)  # [num_patients, D]

patient_indices = range(dataset_subset_size)
starts = val_p2i[np.array(patient_indices), 0].astype(np.int64)
i2p = val[starts, 0].astype(np.uint32)

np.savez("result/patient_embeddings.npz", eid=i2p, emb=embs)


# ================================
# Reformatted and processed patient embeddings
# ================================
z = np.load(f"result/patient_embeddings.npz")
eids, X = z["eid"], z["emb"]

# build column names: eid, emb1..emb120
emb_cols = [f"emb{i+1}" for i in range(X.shape[1])]
emb_df = pd.DataFrame(X, columns=emb_cols)
emb_df.insert(0, "FID", eids) 
emb_df.insert(1, "IID", eids)

# apply normalisation on embeddings then do GWAS, to make embeddings comparable
def inverse_rank_normal_transform(df, columns=None):
    """
    Apply inverse rank normal transformation to specified columns.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe
    columns : list or None
        Columns to transform. If None, all numeric columns are used.
    
    Returns
    -------
    pandas.DataFrame
        Transformed dataframe
    """
    df_out = df.copy()

    if columns is None:
        columns = df_out.select_dtypes(include=[np.number]).columns

    for col in columns:
        x = df_out[col]

        # Rank (average for ties), keep NaNs
        ranks = x.rank(method="average", na_option="keep")

        # Number of non-missing values
        n = ranks.notna().sum()

        # Inverse normal transformation
        df_out[col] = norm.ppf((ranks - 0.5) / n)

    return df_out

cols = emb_df.columns[2:122]
df_int = inverse_rank_normal_transform(emb_df, columns=cols)

df_int.iloc[:, 2:].describe()

df_int.to_csv(f'/gwas/data/norm_embeds_phenotype_441189.txt', sep='\t', index=False)


# load normalised embedding dataset and separate by sex
df_int = pd.read_csv(f'/gwas/data/norm_embeds_phenotype_441189.txt', sep='\t')

# individuals sex info
meta_df = pd.read_csv(f"/data/ukb45266_cov.tab", sep="\t", low_memory=False)
merge_df = pd.merge(df_int, meta_df[["f.eid", "f.31.0.0"]], left_on="FID", right_on="f.eid", how="left")
merge_df["f.31.0.0"].isna().sum()
Counter(merge_df["f.31.0.0"])

female_df = merge_df[merge_df["f.31.0.0"]=="Female"]
female_df = female_df.drop(columns=["f.eid", "f.31.0.0"])
female_df.reset_index(drop=True, inplace=True)
female_df.iloc[:, 2:].describe()
female_df.to_csv(f'/gwas/data/female_norm_embeds_phenotype.txt', sep='\t', index=False)

male_df = merge_df[merge_df["f.31.0.0"]=="Male"]
male_df = male_df.drop(columns=["f.eid", "f.31.0.0"])
male_df.reset_index(drop=True, inplace=True)
male_df.iloc[:, 2:].describe()
male_df.to_csv(f'/gwas/data/male_norm_embeds_phenotype.txt', sep='\t', index=False)
