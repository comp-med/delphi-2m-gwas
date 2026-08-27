import pandas as pd
import tqdm
import numpy as np

folder_dir = '<path_to_file>'

labels_file = f'{folder_dir}/Delphi/data/ukb_simulated_data/labels.csv'
ukb_field_to_icd10_map_file = f'{folder_dir}/Delphi/data/ukb_simulated_data/icd10_codes_mod.tsv'
ubk_basket_tab_file = f'{folder_dir}/data/ukb45266_cov.tab'
ubk_dictionary = f'{folder_dir}/data/Data.dictionary.UKBB.first.occurence.45266.txt'

train_proportion = 0.8 # proportion of full data set to use for training (the rest will be used for validation)
output_prefix = 'ukb_real'


# Read icd10 mapping file and defined index label link

icdict ={}
icdcodes = []
with open(ukb_field_to_icd10_map_file,'r') as f:
    for l in f:
        lvals=l.strip().split()
        icdict[lvals[0]]=lvals[5]
        icdcodes.append(lvals[5])

i = -1
label_dict = {}
with open(labels_file,'r') as f:
    for l in f:
        label_dict[l.strip().split(' ')[0]]=i
        i += 1

# prep monthe dict
month_dict = {
    "January": 1, "February": 2, "March": 3, "April": 4,
    "May": 5, "June": 6, "July": 7, "August": 8,
    "September": 9, "October": 10, "November": 11, "December": 12,
}

# hard coded sex and dob
icdict['f.31.0.0'] = "sex"
icdict['f.34.0.0'] = "YEAR"
icdict['f.52.0.0'] = "MONTH"
icdict['f.40000.0.0'] = "Death"

# cancer fields
for j in range(17):
    icdict['f.40005.'+str(j)+'.0'] = "cancer_date_"+str(j)
    icdict['f.40006.'+str(j)+'.0'] = "cancer_type_"+str(j)

# cancer hes fields 
#for j in range(213):
#    icdict['f.41270.0.'+str(j)] = "hicd_"+str(j)
#    icdict['f.41280.0.'+str(j)] = "hicd_date_"+str(j)

icdict['f.53.0.0'] = "assessment_date"
icdict['f.21001.0.0']="BMI"
icdict['f.1239.0.0']="smoking"
icdict['f.1558.0.0']="alcohol"

# remove the values which not in our dataset, F83 ('f.130968.0.0'), M90 ('f.131982.0.0')
icdcodes = [c for c in icdcodes if c not in {'F83', 'M09'}]

len_icd = len(icdcodes)
#icdcodes.extend(['Death','assessment_date']+['cancer_date_'+str(j) for j in range(17)]+['hicd_date_'+str(j) for j in range(213)])
icdcodes.extend(['Death','assessment_date']+['cancer_date_'+str(j) for j in range(17)])

smoking_to_label = {
    "No": "Smoking_low",
    "Only occasionally": "Smoking_mid",
    "Yes, on most or all days": "Smoking_high",
    "Prefer not to answer": np.nan,   # drop
}

alcohol_to_label = {
    "Never": "Alcohol_low",
    "Special occasions only": "Alcohol_low",

    "One to three times a month": "Alcohol_mid",
    "Once or twice a week": "Alcohol_mid",

    "Three or four times a week": "Alcohol_high",
    "Daily or almost daily": "Alcohol_high",

    "Prefer not to answer": np.nan,
}



# Read ukb basket file in chunks, select icd10 code occurance and dates, format for delphi
data_list = []
ukb_iterator = pd.read_csv(ubk_basket_tab_file, sep='\t',chunksize=1000,index_col=0,low_memory=False)
for _, dd in tqdm.tqdm(enumerate(ukb_iterator)):
    dd = dd.rename(columns=icdict)
    dd.dropna(subset=['sex'], inplace=True)
    dd['sex'] = dd['sex'].map(label_dict)
    #dd['sex'] += 1
    dd = dd[[col for col in dd.columns if not col.startswith('f.')]]
    dd['MONTH'] = dd['MONTH'].map(month_dict)
    dd['dob'] =  pd.to_datetime(dd[['YEAR', 'MONTH']].assign(DAY=1))
    dd[icdcodes] = dd[icdcodes].apply(pd.to_datetime, format="%Y-%m-%d")
    dd[icdcodes]=dd[icdcodes].sub(dd['dob'], axis=0)
    dd[icdcodes]=dd[icdcodes].apply(lambda x : x.dt.days)

    for col in icdcodes[:len_icd+1]:
        X = dd[col].dropna().reset_index().to_numpy().astype(int)
        data_list.append(np.hstack((X,label_dict[col]*np.ones([X.shape[0],1],X.dtype))))

    X = dd['sex'].reset_index().to_numpy().astype(int)
    data_list.append(np.c_[X[:,0],np.zeros(X.shape[0]),X[:,1]].astype(int))

    for j in range(17):
        dd_cancer = dd[['cancer_date_'+str(j),'cancer_type_'+str(j)]].dropna().reset_index()
        if not dd_cancer.empty:
            dd_cancer['cancer'] = dd_cancer['cancer_type_'+str(j)].str.slice(0,3)
            dd_cancer['cancer_label'] = dd_cancer["cancer"].map(label_dict)
            data_list.append(dd_cancer[['f.eid','cancer_date_'+str(j),'cancer_label']].dropna().astype(int).to_numpy())

    #for j in range(213):
    #    dd_hicd = dd[['hicd_date_'+str(j),'hicd_'+str(j)]].dropna().reset_index()
    #    if not dd_hicd.empty:
    #        dd_hicd['hicd'] = dd_hicd['hicd_'+str(j)].str.slice(0,3)
    #        dd_hicd['hicd_label'] = dd_hicd["hicd"].map(label_dict)
    #        data_list.append(dd_hicd[['f.eid','hicd_date_'+str(j),'hicd_label']].dropna().astype(int).to_numpy())

    dd_bmi = dd[['assessment_date','BMI']].dropna().reset_index()
    dd_bmi['bmi_status'] = np.where(dd_bmi['BMI']>28,5,np.where(dd_bmi.BMI>22,4,3))
    data_list.append(dd_bmi[['f.eid','assessment_date','bmi_status']].astype(int).to_numpy())

    dd['smoking'] = dd['smoking'].map(smoking_to_label) # convert to authors value
    dd['smoking'] = dd['smoking'].map(label_dict)
    dd_sm = dd[['assessment_date','smoking']].dropna().reset_index()
    dd_sm = dd_sm[dd_sm['smoking']!=-3]
    dd_sm['smoking_status'] = dd_sm['smoking']
    #dd_sm['smoking_status'] = np.where(dd_sm['smoking']==1,8,np.where(dd_sm.smoking==2,7,6))
    data_list.append(dd_sm[['f.eid','assessment_date','smoking_status']].astype(int).to_numpy())

    dd['alcohol'] = dd['alcohol'].map(alcohol_to_label)
    dd['alcohol'] = dd['alcohol'].map(label_dict)
    dd_al = dd[['assessment_date','alcohol']].dropna().reset_index()
    dd_al = dd_al[dd_al['alcohol']!=-3]
    dd_al['alcohol_status'] = dd_al['alcohol']
    #dd_al['alcohol_status'] = np.where(dd_al['alcohol']==1,11,np.where(dd_al.alcohol < 4,10,9))
    data_list.append(dd_al[['f.eid','assessment_date','alcohol_status']].astype(int).to_numpy())



# reformat, split train and val and output to delphi format
data= np.vstack(data_list)
data = data[np.lexsort((data[:,1], data[:,2]==data[:,2].max(), data[:,0]))]
data = data[data[:,1]>=0]
data = pd.DataFrame(data).drop_duplicates([0,2]).values
data = data.astype(np.uint32)
data.tofile(f'{folder_dir}/Delphi/data/{output_prefix}_data/all.bin')
ids = list(set(data[:,0]))
ids.sort()

train_val_split = data[:,0] <= ids[int(len(ids)*train_proportion)]
data[train_val_split].tofile(f'{folder_dir}/Delphi/data/{output_prefix}_data/train.bin')
data[~train_val_split].tofile(f'{folder_dir}/Delphi/data/{output_prefix}_data/val.bin')
