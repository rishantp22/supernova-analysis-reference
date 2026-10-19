#!/usr/bin/env python3
"""Census all legacy HDF5 particles/species, without executing a nuclear network."""
import concurrent.futures as futures
import sys,json,re,resource
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'software'))
import h5py
import numpy as np
from common import *
BASE=Path('/data/pyoung3/supernova_data/50Am/50Am.dir')

def flat(ds):
    a=ds[()]
    return np.asarray(a[a.dtype.names[0]] if a.dtype.names else a).ravel()
def work(path):
    p=Path(path);dest=ROOT/'tables/raw_old'/p.stem;dest.mkdir(parents=True,exist_ok=True)
    receipt=dest/'receipt.json'
    if receipt.exists():return json.loads(receipt.read_text())
    before=stamp(p);log('Reading '+p.name)
    with h5py.File(p,'r') as f:
        z=flat(f['nz']).astype('i4');nn=flat(f['nn']).astype('i4');A=z+nn
        keys=[k for k in f if k.startswith('cycle')];ids=np.array([int(k[5:]) for k in keys],dtype='i8');n=len(ids)
        values=np.lib.format.open_memmap(dest/'mass_fractions.npy',mode='w+',dtype='f8',shape=(n,len(z)),fortran_order=True)
        mass=np.zeros(n);stop=np.zeros(n);loop=np.zeros(n,dtype='i8');sums=np.zeros(n);ye=np.zeros(n);history=[];offset=[0];negative=np.zeros(n,dtype='i4');nonfinite=negative.copy()
        for i,key in enumerate(keys):
            g=f[key];x=flat(g['fmass']);assert len(x)==len(z);values[i,:]=x
            mass[i]=np.asarray(g.attrs['mass']).ravel()[0];stop[i]=np.asarray(g.attrs['tfinal']).ravel()[0];loop[i]=np.asarray(g.attrs.get('looperror',[-999])).ravel()[0]
            sums[i]=x.sum();ye[i]=np.sum(x*z/A);negative[i]=(x<0).sum();nonfinite[i]=(~np.isfinite(x)).sum()
            h=g['SE_DATASET'][()];history.append(np.column_stack([h['time'],h['temp'],h['rho']]));offset.append(offset[-1]+len(h))
            if i%10000==0:log(f'{p.name}: {i:,}/{n:,}')
        values.flush();np.savez_compressed(dest/'metadata.npz',id=ids,Z=z,A=A,mass_g=mass,tfinal_s=stop,looperror=loop,sumX=sums,Ye=ye,negative=negative,nonfinite=nonfinite)
        np.savez_compressed(dest/'histories.npz',id=ids,offset=np.array(offset,dtype='i8'),time_temp_rho=np.concatenate(history))
        attrs={k:np.asarray(v).astype(str).tolist() for k,v in f.attrs.items()}
    assert stamp(p)==before
    r=dict(source=before,rows=n,species=len(z),duplicate_species=len(z)-len(set(zip(z.tolist(),A.tolist()))),invalid_values=int(nonfinite.sum()),negative_values=int(negative.sum()),sumX=percentiles(sums),Ye=percentiles(ye),tfinal=percentiles(stop),looperror_counts={str(k):int(v) for k,v in zip(*np.unique(loop,return_counts=True))},attributes=attrs)
    receipt.write_text(json.dumps(clean(r),indent=2)+'\n');return r

def combine():
    dirs=sorted((ROOT/'tables/raw_old').glob('out*'));ids=[];zs=[];aa=[];metas=[]
    for p in dirs:
        d=np.load(p/'metadata.npz');ids.append(d['id']);zs.append(d['Z']);aa.append(d['A']);metas.append(d)
    assert all(np.array_equal(zs[0],z) and np.array_equal(aa[0],a) for z,a in zip(zs,aa))
    allids=np.concatenate(ids);unique,first,counts=np.unique(allids,return_index=True,return_counts=True)
    save('provenance/old_raw_population.json',dict(rows=len(allids),unique_ids=len(unique),duplicate_ids=int((counts>1).sum()),files=len(dirs),species=len(zs[0])))
    csvsave('tables/old_raw_duplicate_ids.csv',{'id':unique[counts>1],'count':counts[counts>1]})
    om,old,oi=table('old');nm,new,ni=table('new');ou,orp=unique_rows(oi);nu,nrp=unique_rows(ni)
    mapping={int(v):i for i,v in enumerate(allids)}
    # A duplicate raw particle cannot be given an unambiguous run association.
    valid=unique[counts==1];shared,oix=np.intersect1d(ou,valid,return_indices=True)[:2]
    values=[np.load(p/'mass_fractions.npy',mmap_mode='r') for p in dirs]
    rawmass=np.concatenate([d['mass_g'] for d in metas]);tfinal=np.concatenate([d['tfinal_s'] for d in metas]);rawsum=np.concatenate([d['sumX'] for d in metas]);rawye=np.concatenate([d['Ye'] for d in metas])
    inventory=[];checks=[];yields=[];raw_new=[]
    from population_audit import stats
    shared_new,ix,iy=np.intersect1d(allids,nu,return_indices=True)
    np.savez_compressed(ROOT/'tables/raw_old_population_metadata.npz',id=allids,mass_g=rawmass,tfinal_s=tfinal,sumX=rawsum,Ye=rawye)
    for j,(z,A) in enumerate(zip(zs[0],aa[0])):
        name='n' if z==0 and A==1 else f'{A}{SYMBOLS[z]}'
        x=np.concatenate([v[:,j] for v in values]);inventory.append(dict(index=j,Z=z,A=A,name=name,in_old_plot=name in om['columns'],in_new_plot=name in nm['columns']))
        yields.append(dict(species=name,Z=z,A=A,raw_particles=len(x),nonzero=int((x>0).sum()),yield_msun=float(np.sum(x*rawmass)/MSUN),**percentiles(x)))
        if name in om['columns']:
            shared_ids,ri,pi=np.intersect1d(allids,ou,return_indices=True);xp=old[orp[pi],om['columns'].index(name)];xr=x[ri]
            # Burn query writes six digits after the decimal in scientific notation.
            ok=np.abs(xp-xr)<=5.1e-7*np.abs(xr)
            checks.append(dict(species=name,pairs=len(xp),matches_query_precision=int(ok.sum()),raw_nonzero_plotted_zero=int(((xr>0)&(xp==0)).sum()),max_abs_delta=float(np.max(abs(xp-xr))),raw_max_omitted=float(np.max(np.where(xp==0,xr,0)))))
        if name in nm['columns']:
            xn=new[nrp[iy],nm['columns'].index(name)]
            row,_=stats(x[ix],xn,shared_new,rawmass[ix],new[nrp[iy],nm['columns'].index('mass')],1e-5);row.update(species=name,Z=z,A=A);raw_new.append(row)
    csvsave('tables/old_raw_species_inventory.csv',inventory);csvsave('tables/old_raw_species_yields.csv',yields);csvsave('tables/old_raw_vs_plotting.csv',checks);csvsave('tables/old_raw_vs_new_selected_species.csv',raw_new)
    save('provenance/raw_audit_summary.json',dict(raw_rows=len(allids),old_plot_rows=len(ou),old_plot_ids_without_raw=len(np.setdiff1d(ou,allids)),new_ids_with_raw=len(shared_new),raw_sumX=percentiles(rawsum),raw_Ye=percentiles(rawye),raw_tfinal_s=percentiles(tfinal),raw_over_unity_1e5=int((rawsum>1.00001).sum()),raw_under_unity_1e5=int((rawsum<.99999).sum()),max_raw_sum_id=int(allids[np.argmax(rawsum)])))
    log('Legacy raw output census complete')

def main():
    paths=sorted(BASE.glob('out*.h5'));(ROOT/'tables/raw_old').mkdir(exist_ok=True)
    with futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(work,map(str,paths)))
    save('provenance/old_hdf5_files.json',results);combine()
if __name__=='__main__':main()
