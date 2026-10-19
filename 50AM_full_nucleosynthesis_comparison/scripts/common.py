"""Shared read-only audit utilities. No source dataset is opened for writing."""
from pathlib import Path
import hashlib, json, re, time
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SYMBOLS = 'n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og'.split()
ZS = {s:z for z,s in enumerate(SYMBOLS)}
MSUN = 1.9889e33

def log(s): print(time.strftime('%Y-%m-%d %H:%M:%S'),s,flush=True)
def stamp(p):
    s=Path(p).stat()
    return dict(path=str(Path(p).resolve()),bytes=s.st_size,mtime_ns=s.st_mtime_ns,ctime_ns=s.st_ctime_ns,inode=s.st_ino)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        while b:=f.read(8*1024*1024): h.update(b)
    return h.hexdigest()
def clean(v):
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,np.integer):return int(v)
    if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
    if isinstance(v,np.bool_):return bool(v)
    return v
def save(name,v):
    p=ROOT/name;p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(clean(v),indent=2,allow_nan=False)+'\n')
def csvsave(name,rows):
    p=ROOT/name;p.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(p,index=False)
def canon(s):
    s=s.strip()
    if s.startswith('X_{') and s.endswith('}'):s=s[3:-1]
    if s.lower() in ['id','density','peak temp','peak density']:s=s.lower()
    return s
def species(s):
    s=canon(s)
    if s=='n':return ('isotope',0,1)
    m=re.fullmatch(r'(\d+)([A-Z][a-z]?)',s)
    if m and m[2] in ZS:return ('isotope',ZS[m[2]],int(m[1]))
    if s in ZS:return ('element',ZS[s],None)
    return ('physical',None,None)
def percentiles(v):
    a=np.asarray(v);a=a[np.isfinite(a)]
    if not len(a):return {}
    return dict(zip(['min','p01','p05','p25','p50','p75','p95','p99','max'],np.quantile(a,[0,.01,.05,.25,.5,.75,.95,.99,1]).tolist()))
def table(label):
    meta=json.loads((ROOT/f'provenance/{label}_scan.json').read_text())
    a=np.load(ROOT/f'tables/{label}_values.npy',mmap_mode='r')
    ids=np.load(ROOT/f'tables/{label}_ids.npy',mmap_mode='r')
    return meta,a,ids
def col(label,name):
    m,a,ids=table(label)
    return a[:,m['columns'].index(canon(name))]
def unique_rows(ids):
    u,first,n=np.unique(ids,return_index=True,return_counts=True)
    ok=(u>=0)&(n==1)
    return u[ok],first[ok]
def matched():
    _,_,oi=table('old');_,_,ni=table('new')
    ou,orow=unique_rows(oi);nu,nrow=unique_rows(ni)
    ids,i,j=np.intersect1d(ou,nu,return_indices=True)
    return ids,orow[i],nrow[j]

