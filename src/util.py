import codecs
import os
import pickle
import warnings

import joblib
import numpy as np
import pandas as pd
from mordred import Calculator, descriptors
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, MACCSkeys, PandasTools
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from tqdm import tqdm

warnings.simplefilter('ignore')

def baseline_data(df, is_train=True):
    X=None
    y=None
    if is_train:
        X = df.drop(["SMILES","log P (octanol-water)"],axis=1)
        y = df["log P (octanol-water)"]
    else:
        X = df.drop(["SMILES"],axis=1)
    return X,y

def base44(df, is_train=True):
    X=None
    y=None
    select44=joblib.load(filename=os.path.dirname(__file__) + "/li_select44in200.joblib")
    if is_train:
        X = df[select44]
        y = df["log P (octanol-water)"]
    else:
        X = df[select44]
    return X,y

def load_data(df):
    y = None
    df = df.dropna()
    desc_df = create200(df)
    # model = pickle.load(open(os.path.dirname(__file__) + "/model_q2.pkl", "rb"))
    # desc_df["pred_q2"] = model.predict(desc_df)
    PandasTools.AddMoleculeColumnToFrame(df,'SMILES')
    df = df.dropna()

    # １列にfingerprintのリストを追加する場合
    df['FP'] = df.apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x.ROMol, 2, 1024), axis=1)

    # fingerprintの各値を各列に格納する場合
    # 個別に０１をデータフレームに格納する
    FP = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024) for mol in df.ROMol]
    df_FP = pd.DataFrame(np.array(FP))

    # フィンガープリントをもとのデータフレームに結合
    df_FP.index = df.index
    df = pd.concat([df, df_FP, desc_df], axis=1)
    df = df.loc[:,~df.columns.duplicated()]
    # df = add_count_col(df)
    if "IC50 (nM)" in df.columns:
        y = df["IC50 (nM)"]
        y = np.log1p(y)
        X = df.drop(["SMILES", "ROMol", "FP", "IC50 (nM)"], axis=1)
    else:
        X = df.drop(["SMILES", "ROMol", "FP"], axis=1)
    f = open("./src/imp0.txt","rb")
    imp0 = pickle.load(f)
    print(X.shape)
    X = X.drop(imp0, axis=1)
    print(X.shape)
    # X = X.drop(X.columns[np.isnan(X).any()], axis=1)
    print(X.shape)
    ss = StandardScaler()
    X = ss.fit_transform(X)
    # X = desc_df
    # X = create_mordred(df=df)
    # X = create_features(df["SMILES"])
    # X = pd.concat(objs=[X,df_des],axis=1)
    # X = X.drop(labels=["MOL"],axis=1)

    return X,y

def count_elements(df):
    df_ = None
    df_add = pd.DataFrame(index=df.index, columns=[])
    elements=['C', 'O', 'N', 'Cl', 'Na', 'P', 'S', 'Si', 'Br', 'F', 'I', 'As', 'B', 'Hg', 'Se', 'Li', 'Mg', 'K', 'H']
    l_l = []
    for i in elements:
        l=[]
        for smile in df["SMILES"]:
            if (i=='C'):
                l.append(smile.count('C')+smile.count('c')-smile.count('Cl'))
            elif (i=='N'):
                l.append(smile.count('N')+smile.count('n')-smile.count('Na'))
            elif (i=='S'):
                l.append(smile.count('S')+smile.count('s')-smile.count('Si'))
            elif (i=='H'):
                l.append(smile.count('H')-smile.count('Hg'))
            else:
                l.append(smile.count(i))
        l_l.append(l)
    df_ = pd.DataFrame(l_l,index=elements).T
    return df_

def count_ionized(df):
    df_add = pd.DataFrame(index=df.index, columns=[])
    lll=[]
    for smile in df["SMILES"]:
        lll.append(smile.count('+')+smile.count('-'))
    df_add['Ionized number']=lll
    feature_l = ['=O', '(', 'O-', 'N+', '1', '2', '3', '4', '5']
    for i in feature_l:
        l=[]
        for smile in df["SMILES"]:
            if (i=='1' or i=='2' or i=='3' or i=='4' or i=='5'):
                l.append(smile.count(i) / 2)
            else:
                l.append(smile.count(i))
        df_add['add_'+i] = l
    return df_add

def create200(df):
    smiles_all = df.iloc[:, 0].values
    mols = []
    for smi in smiles_all:
        mol = Chem.MolFromSmiles(smi)
        if mol != None:
            mols.append(mol)
        else:
            print(smi)
    descs = [desc_name[0] for desc_name in Descriptors._descList]
    desc_calc = MoleculeDescriptors.MolecularDescriptorCalculator(descs)
    result = []
    for i, mol in enumerate(mols):
        result.append(desc_calc.CalcDescriptors(mol))

    desc_df = pd.DataFrame(result)
    desc_df.columns = descs
    return desc_df

def count_num_SMILES(df):
    s = df["SMILES"]
    count_num = lambda x: len(x)
    return s.map(count_num)

def make_fingerprint(df,n=1):
    PandasTools.AddMoleculeColumnToFrame(df,'SMILES')
    # df = df.dropna()
    FP = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, n*1024) for mol in df["ROMol"]]
    df_FP = pd.DataFrame(np.array(FP))
    df_FP.index = df.index
    df_FP = df_FP.add_prefix('FP_')
    return df_FP

def create_mordred(df):
    df['MOL'] = df['SMILES'].apply(Chem.MolFromSmiles)
    calc = Calculator(descriptors, ignore_3D=False)
    df_descriptors_mordred = calc.pandas(df['MOL'])
    df_descriptors = df_descriptors_mordred.astype(str)
    masks = df_descriptors.apply(lambda d: d.str.contains('[a-zA-Z]' ,na=False))
    df_descriptors = df_descriptors[~masks]
    df_descriptors = df_descriptors.astype(float)
    return df_descriptors

def top200(X,num=200):
    pf = PolynomialFeatures(degree=2, interaction_only=True, include_bias=False)
    df_pf = pd.DataFrame(pf.fit_transform(X.fillna(-99999)))
    # cols = pickle.load(open(os.path.dirname(__file__) + "/sidfcol.pkl", "rb"))
    cols = joblib.load(filename=os.path.dirname(__file__) + "/list_sidfcol20100.joblib")
    X = df_pf[cols[:num]]
    X = X.add_prefix('top200_')
    return X

def morder(df):
    mols = df['SMILES'].apply(Chem.MolFromSmiles)
    df_descriptors_mordred = Calculator(descriptors, ignore_3D=False).pandas(mols)
    # df_descriptors = df_descriptors_mordred.astype(str)
    # masks = df_descriptors.apply(lambda d: d.str.contains('[a-zA-Z]' ,na=False))
    # df_descriptors = df_descriptors[~masks]
    # df_descriptors = df_descriptors.astype(float)
    return df_descriptors_mordred

def mord(df,num=250):
    sidfcolvalues_mord = joblib.load(filename=os.path.dirname(__file__) + "/sidfcolvalues_mord.joblib")
    mols = df['SMILES'].apply(Chem.MolFromSmiles)
    calc_dummy = Calculator(descriptors, ignore_3D=False)
    my_desc_names = sidfcolvalues_mord[:250]
    my_descs = []
    for i, desc in enumerate(calc_dummy.descriptors):
        if desc.__str__()  in my_desc_names:
            my_descs.append(desc)

    calc_real = Calculator(my_descs, ignore_3D=False)
    df_mord = calc_real.pandas(mols, nproc=-1)
    df_mord = df_mord.replace(['[a-z](.*)'], np.nan, regex=True).astype(float)
    df_mord = df_mord.fillna(-99999)
    df_mord = df_mord.add_prefix('mordred_')
    return df_mord

def mord66(df):
    # sidfcolvalues_mord = joblib.load(filename=os.path.dirname(__file__) + "/sidfcolvalues_mord.joblib")
    select_col_mord = joblib.load(filename=os.path.dirname(__file__) + '/li_select66in1826.joblib')
    mols = df['SMILES'].apply(Chem.MolFromSmiles)
    calc_dummy = Calculator(descriptors, ignore_3D=False)
    my_desc_names = select_col_mord
    my_descs = []
    for i, desc in enumerate(calc_dummy.descriptors):
        if desc.__str__()  in my_desc_names:
            my_descs.append(desc)

    calc_real = Calculator(my_descs, ignore_3D=False)
    df_mord = calc_real.pandas(mols, nproc=8)
    df_mord = df_mord.replace(['[a-z](.*)'], np.nan, regex=True).astype(float)
    df_mord = df_mord.fillna(0)
    df_mord = df_mord.add_prefix('mordred_')
    return df_mord

def top200_fp(df):
    X_top200 = top200(df,200,)
    X_fp = make_fingerprint(df,2,)
    X = pd.concat([X_top200,X_fp],axis=1)
    X.columns = list(range(len(X.columns)))
    return X

def MACCS(df):
    fpss = []
    for smi in tqdm(df.SMILES):
        mol = Chem.MolFromSmiles(smi)
        fps = MACCSkeys.GenMACCSKeys(mol)
        fp_bits = tuple(fps.GetOnBits())
        fp_arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fps,fp_arr)
        fpss.append(np.array(fp_arr))
    df_fpss = pd.DataFrame(fpss)
    df_fpss = df_fpss.add_prefix('MACCS_')
    return df_fpss

def del_VT(X,is_train=True):
    X_VT_df = None
    VT = VarianceThreshold()
    X_VT = VT.fit_transform(X.fillna(-99999))
    if is_train:
        X_VT_df = pd.DataFrame(X_VT, columns=X.columns[VT.get_support()==True])
        joblib.dump(value=X.columns[VT.get_support()==True],filename=open(file=os.path.dirname(__file__)+"/X_VT_col.joblib",mode="wb"),compress=3)
    else:
        X_VT_col = joblib.load(filename=open(file=os.path.dirname(__file__)+"/X_VT_col.joblib",mode="rb"))
        X_VT_df = X[X_VT_col]
    return X_VT_df

def convert_job(X,converter):
    X_conv = pd.DataFrame(converter.fit_transform(X))
    if X_conv.shape[1] == X.shape[1]:
        X_conv.columns = X.columns
    return X_conv

def prepare_data(df,is_train=True):
    df = df.dropna()
    df_,_=baseline_data(df,is_train)
    df_200 = top200(df_,num=50)
    df.index = list(range(df.shape[0]))
    X,y=base44(df,is_train)
    X["cout_smi"]=count_num_SMILES(df)
    # df_mord66 = mord66(df)
    # X,y=baseline_data(df,is_train)
    dfs = [
        X,
        # mord66(df),
        count_ionized(df),
        count_elements(df),
        # top200(X,num=100),
        df_200,
        # mord(df,200)
    ]
    X = pd.concat(dfs,axis=1)

    # # X = convert_job(X,PolynomialFeatures(degree=2, interaction_only=True, include_bias=False))
    # # X = convert_job(X,TruncatedSVD(n_components=1000, n_iter=7, random_state=42))
    X = convert_job(X,StandardScaler())
    fps = [
        X,
        MACCS(df),
        # make_fingerprint(df,2),
    ]
    X = pd.concat(fps,axis=1)
    # X = MACCS(df)
    # X=make_fingerprint(df)
    # X = del_VT(X,is_train)
    # col = joblib.load(os.path.dirname(__file__)+"/sidf_2875.joblib")
    # col.remove("SMILES")
    # X = X[col[:650]]
    # X=X.T.drop_duplicates().T
    # X = convert_job(X,TruncatedSVD(n_components=50, n_iter=7, random_state=42))
    # X = convert_job(X,StandardScaler())
    # X=X.loc[:,~X.columns.duplicated(keep='last')]
    return X,y
