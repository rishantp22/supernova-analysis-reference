#!/usr/bin/env python3
"""Audit every supplied A/B chunk and every ID/category/species key on Slurm."""
import concurrent.futures as futures
import hashlib,io,json,re,traceback
from pathlib import Path
import numpy as np
import pandas as pd
from common import *

CATS=['closest_to_100','closest_to_1000','last']
BASE=Path('/scratch/rpraka15/final_chunks')

def process(number):
    out=ROOT/f'tables/chunks/chunk_{number}.npz';report=ROOT/f'provenance/chunks/chunk_{number}.json'
    previous=json.loads(report.read_text()) if report.exists() else None
    if previous and previous['status']=='complete':return previous
    a_path=BASE/f'chunk_{number}.csv';b_path=BASE/'final_chunks_B'/f'chunk_{number}_B.csv'
    record={'chunk':number};frames=[]
    if previous:record['previous_attempt_error']=previous.get('traceback')
    try:
        for label,p,value in [('A',a_path,'C_sum'),('B',b_path,'mass_times_Csum')]:
            before=stamp(p);data=p.read_bytes();f=pd.read_csv(io.BytesIO(data),dtype={'ID':'i8','category':pd.CategoricalDtype(CATS),'A':'i2','B':'i2'})
            after=stamp(p);digest=hashlib.sha256(data).hexdigest()
            if after!=before:
                record[label+'_metadata_changed_during_read']={'before':before,'after':after}
                assert sha(p)==digest, 'Source bytes changed during read'
                after=stamp(p)
            record[label+'_source']={**after,'sha256':digest};del data
            assert not f.isna().any().any();assert (f.ID>=0).all();assert np.isfinite(f[value]).all();assert (f[value]>=0).all()
            key=(f.ID.to_numpy()*3+f.category.cat.codes.to_numpy(dtype='i8'))*65536+f.A.to_numpy(dtype='i8')*512+f.B.to_numpy(dtype='i8')
            order=np.argsort(key);keys=key[order];assert not (np.diff(keys)==0).any()
            frames.append((f.iloc[order].reset_index(drop=True),keys))
        a,ka=frames[0];b,kb=frames[1];assert np.array_equal(ka,kb);del ka,kb,frames
        ids=np.unique(a.ID.to_numpy());code=a.A.to_numpy(dtype='i8')*512+a.B.to_numpy(dtype='i8');sp=np.unique(code)
        newmeta,newtable,newids=table('new');newcols=newmeta['columns']
        expected=sorted(species(c)[1]*512+species(c)[2] for c in newcols if species(c)[0]=='isotope' and c!='n')
        assert sp.tolist()==expected
        count=a.groupby(['ID','category'],observed=True).size();assert len(count)==len(ids)*3 and (count==len(sp)).all()
        n=len(ids);k=len(sp);values=a.C_sum.to_numpy().reshape(n,3,k);weighted=b.mass_times_Csum.to_numpy().reshape(n,3,k)
        times=a.actual_time.to_numpy().reshape(n,3,k);steps=a.actual_timestep.to_numpy().reshape(n,3,k)
        assert (times==times[:,:,:1]).all() and (steps==steps[:,:,:1]).all()
        assert (np.diff(times[:,:,0],axis=1)>=0).all()
        om,ot,oi=table('old');ou,op=unique_rows(oi);loc=np.searchsorted(ou,ids);assert (loc<len(ou)).all() and np.array_equal(ou[loc],ids)
        mass=ot[op[loc],om['columns'].index('mass')]/1.9889e27
        expected_mass=values*mass[:,None,None];b_bad=abs(weighted-expected_mass)>2e-5*abs(expected_mass)
        sums=values.sum(axis=2);ye=(values*(sp//512/(sp%512))[None,None,:]).sum(axis=2)
        nu,npidx=unique_rows(newids);loc=np.searchsorted(nu,ids);safe=np.minimum(loc,len(nu)-1);present=(loc<len(nu))&(nu[safe]==ids)
        plotcols=[newcols.index(f'{v%512}{SYMBOLS[v//512]}') for v in sp]
        actual=np.asarray(newtable[npidx[safe[present]],:])[:,plotcols]
        delta=actual-values[present,2,:];mismatches=int(np.sum(delta!=0));nearbad=int(np.sum(abs(delta)>1e-10*abs(values[present,2,:])))
        global_ids=ids[present];discrepant=np.any(delta!=0,axis=1)
        record.update(status='complete',ids=n,A_rows=len(a),B_rows=len(b),species=k,categories=3,A_B_key_mismatches=0,duplicate_keys=0,missing_grid_cells=0,negative_or_nonfinite_values=0,mass_weight_mismatches=int(b_bad.sum()),mass_weight_rtol=2e-5,plot_shared_ids=int(present.sum()),plot_missing_ids=ids[~present].tolist(),plot_exact_mismatch_cells=mismatches,plot_mismatches_rtol_1e10=nearbad,plot_mismatch_ids=global_ids[discrepant].tolist(),source_zeros=int((values==0).sum()),sum_min=float(sums.min()),sum_max=float(sums.max()),last_sum_max=float(sums[:,2].max()),last_sum_max_id=int(ids[np.argmax(sums[:,2])]),last_over_unity=int((sums[:,2]>1.00001).sum()),last_time=percentiles(times[:,2,0]))
        out.parent.mkdir(exist_ok=True)
        np.savez_compressed(out,id=ids,time=times[:,:,0],step=steps[:,:,0],isotope_sum=sums,partial_ye=ye,zero_count=(values==0).sum(axis=2),plot_present=present,plot_mismatch_count=np.pad((delta!=0).sum(axis=1),(0,n-int(present.sum()))) if present.all() else np.array([int(np.sum(delta[np.flatnonzero(np.flatnonzero(present)==j)[0]]!=0)) if present[j] else -1 for j in range(n)]))
    except Exception as e:record.update(status='error',error=repr(e),traceback=traceback.format_exc())
    save(str(report.relative_to(ROOT)),record);return record

def main():
    (ROOT/'provenance/chunks').mkdir(exist_ok=True);(ROOT/'tables/chunks').mkdir(exist_ok=True)
    aa={int(re.search(r'chunk_(\d+)\.csv',p.name)[1]) for p in BASE.glob('chunk_*.csv')};bb={int(re.search(r'chunk_(\d+)_B\.csv',p.name)[1]) for p in (BASE/'final_chunks_B').glob('chunk_*_B.csv')}
    save('provenance/chunk_sets.json',dict(A_chunks=sorted(aa),B_chunks=sorted(bb),A_only=sorted(aa-bb),B_only=sorted(bb-aa),missing_numbers=sorted(set(range(1,max(aa|bb)+1))-(aa|bb))))
    results=[]
    with futures.ProcessPoolExecutor(max_workers=4) as pool:
        for r in pool.map(process,sorted(aa|bb)):
            results.append(r)
            if len(results)%25==0 or r['status']=='error':log(f'Chunk {r["chunk"]}: {r["status"]}; {len(results)}/{len(aa|bb)}')
    csvsave('tables/chunk_integrity.csv',[{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in results])
    valid=[r for r in results if r['status']=='complete'];allids=[];allchunk=[];epochs=[]
    for r in valid:
        with np.load(ROOT/f'tables/chunks/chunk_{r["chunk"]}.npz') as d:
            allids.append(d['id']);allchunk.append(np.full(len(d['id']),r['chunk']));epochs.append({k:d[k].copy() for k in d.files})
    arrays={k:np.concatenate([d[k] for d in epochs]) for k in epochs[0]};arrays['chunk']=np.concatenate(allchunk)
    np.savez_compressed(ROOT/'tables/all_chunk_particles.npz',**arrays)
    u,counts=np.unique(arrays['id'],return_counts=True)
    csvsave('tables/duplicate_ids_across_chunks.csv',{'id':u[counts>1],'count':counts[counts>1]})
    ni=table('new')[2];oi=table('old')[2]
    save('provenance/chunk_summary.json',dict(chunks=len(results),complete_chunks=len(valid),failed_chunks=len(results)-len(valid),A_bytes=sum(r['A_source']['bytes'] for r in valid),B_bytes=sum(r['B_source']['bytes'] for r in valid),A_rows=sum(r['A_rows'] for r in valid),B_rows=sum(r['B_rows'] for r in valid),unique_ids=len(u),duplicate_ids_across_chunks=int((counts>1).sum()),A_B_mass_mismatches=sum(r['mass_weight_mismatches'] for r in valid),plot_exact_mismatch_cells=sum(r['plot_exact_mismatch_cells'] for r in valid),plot_approx_mismatch_cells=sum(r['plot_mismatches_rtol_1e10'] for r in valid),chunk_ids_absent_from_valid_new=np.setdiff1d(u,ni).tolist(),valid_new_ids_absent_from_chunks=np.setdiff1d(ni[ni>=0],u).tolist(),old_ids_absent_from_chunks_count=len(np.setdiff1d(oi[oi>=0],u)),last_time=percentiles(arrays['time'][:,2]),last_sum=percentiles(arrays['isotope_sum'][:,2]),last_sum_over_unity_ids=arrays['id'][arrays['isotope_sum'][:,2]>1.00001].tolist(),archive_run_identity='Not established by chunk-only keys; no model/run identifier in original CSV'))
    log('Complete chunk grid audit finished')
if __name__=='__main__':main()
