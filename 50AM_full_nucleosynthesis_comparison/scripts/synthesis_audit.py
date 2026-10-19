#!/usr/bin/env python3
"""Full-population attribution and independent raw-output crosschecks on Slurm."""
import sys,json,hashlib,collections
from pathlib import Path
import numpy as np
from common import *
from population_audit import stats
sys.path.insert(0,'/scratch/rpraka15/supernova-grain-analysis/code')
from readers import time_blocks

def coverage():
    om,ot,oi=table('old');nm,nt,ni=table('new');ids,orr,nrr=matched()
    rd=np.load(ROOT/'tables/raw_old_population_metadata.npz');covered=np.isin(ids,rd['id']);allcovered=np.isin(oi,rd['id'])
    mass=ot[orr,om['columns'].index('mass')]/MSUN;rows=[];nonzero=np.zeros(len(oi),bool)
    for j,s in enumerate(om['columns']):
        if species(s)[0]!='isotope':continue
        nonzero|=np.isfinite(ot[:,j])&(ot[:,j]!=0)
        if s not in nm['columns']:continue
        x=ot[orr,j];y=nt[nrr,nm['columns'].index(s)];r=dict(species=s)
        for label,mask in [('raw_covered',covered),('raw_absent',~covered)]:
            good=mask&np.isfinite(x)&np.isfinite(y)
            r.update({label+'_pairs':int(good.sum()),label+'_old_yield_msun':float(np.sum(mass[good]*x[good])) if good.any() else np.nan,label+'_new_yield_msun':float(np.sum(mass[good]*y[good])) if good.any() else np.nan,label+'_significant':int((good&(abs(y-x)>1e-5*abs(x))).sum())})
        rows.append(r)
    csvsave('tables/yield_coverage_decomposition.csv',rows)
    flags=collections.Counter();history=collections.Counter()
    for p in sorted((ROOT/'tables/raw_old').glob('out*')):
        d=np.load(p/'metadata.npz');flags.update({str(k):int(v) for k,v in zip(*np.unique(d['looperror'],return_counts=True))})
        h=np.load(p/'histories.npz');v=h['time_temp_rho'];history.update(points=len(v),nonfinite=int((~np.isfinite(v)).sum()),zero_temperature=int((v[:,1]==0).sum()),zero_density=int((v[:,2]==0).sum()))
    save('provenance/coverage_attribution.json',dict(old_ids_without_raw=int((~allcovered).sum()),old_ids_without_raw_all_exported_isotopes_zero=int((~allcovered&~nonzero).sum()),old_ids_without_raw_any_nonzero=int((~allcovered&nonzero).sum()),shared_raw_covered=int(covered.sum()),shared_raw_absent=int((~covered).sum()),mass_shared_raw_covered_msun=float(mass[covered].sum()),mass_shared_raw_absent_msun=float(mass[~covered].sum()),legacy_looperror_counts=dict(flags),legacy_history_validation=dict(history)))
    np.savez_compressed(ROOT/'tables/raw_coverage_alignment.npz',id=ids,legacy_raw_available=covered)

def eoc_database():
    inv=pd.read_csv(ROOT/'tables/eoc_species_inventory.csv');codes=inv.code.to_numpy();names=inv.species.tolist();e=np.load(ROOT/'tables/eoc_particle_validation.npz');ids=e['id'];order=np.argsort(ids);sortedids=ids[order]
    assert len(np.unique(ids))==len(ids)
    v=np.lib.format.open_memmap(ROOT/'tables/eoc_values.npy',mode='w+',dtype='f8',shape=(len(ids),len(codes)),fortran_order=True);v[:]=np.nan
    present=np.lib.format.open_memmap(ROOT/'tables/eoc_row_present.npy',mode='w+',dtype='bool',shape=v.shape,fortran_order=True);present[:]=False
    seen=[]
    for p in sorted((ROOT/'tables/eoc').glob('block_*.npz')):
        d=np.load(p);ii=d['id'];seen.append(ii);rr=order[np.searchsorted(sortedids,ii)];r=np.repeat(rr,np.diff(d['offset']));c=np.searchsorted(codes,d['Z'].astype('i8')*512+d['A']);v[r,c]=d['X'];present[r,c]=True
    seen=np.concatenate(seen);assert len(seen)==len(ids) and np.array_equal(np.sort(seen),sortedids),'Missing, duplicated, or stale EOC blocks'
    v.flush();present.flush();np.save(ROOT/'tables/eoc_ids.npy',ids);log('Memory-mapped EOC database created; all ragged blocks validated')
    dirs=sorted((ROOT/'tables/raw_old').glob('out*'));rawids=np.concatenate([np.load(p/'metadata.npz')['id'] for p in dirs]);mass=np.concatenate([np.load(p/'metadata.npz')['mass_g'] for p in dirs]);rawvalues=[np.load(p/'mass_fractions.npy',mmap_mode='r') for p in dirs];meta=np.load(dirs[0]/'metadata.npz');rcodes=meta['Z'].astype('i8')*512+meta['A'];shared,ri,ei=np.intersect1d(rawids,ids,return_indices=True)
    summary=[]
    for j,c in enumerate(rcodes):
        if c not in codes:continue
        k=int(np.searchsorted(codes,c));x=np.concatenate([a[:,j] for a in rawvalues])[ri];y=v[ei,k]
        s,_=stats(x,y,shared,mass[ri],mass[ri],1e-5);s.update(species=names[k],Z=int(c//512),A=int(c%512));summary.append(s)
    csvsave('tables/old_raw_vs_new_eoc_all_shared_species.csv',summary)
    nm,nt,ni=table('new');nu,nr=unique_rows(ni);common,ei,ni2=np.intersect1d(ids,nu,return_indices=True);newrows=nr[ni2]
    reported=nt[newrows,nm['columns'].index('Y_{e}')];delta=e['Ye'][ei]-reported
    save('provenance/eoc_database_summary.json',dict(rows=len(ids),species=len(codes),legacy_raw_common_ids=len(shared),legacy_raw_common_species=len(summary),new_plot_common_ids=len(common),composition_minus_reported_Ye=percentiles(delta),composition_vs_reported_abs_gt_1e5=int((abs(delta)>1e-5).sum()),missingness='eoc_row_present distinguishes absent isotope rows from explicit NaN; neither is zero'))
    rows=[];cells=[]
    for k in np.flatnonzero(e['plot_mismatches'][ei]>0):
        ident=common[k];er=ei[k];pr=newrows[k]
        for j,name in enumerate(names):
            if name not in nm['columns']:continue
            x=v[er,j];y=nt[pr,nm['columns'].index(name)]
            if np.isfinite(x) and np.isfinite(y) and abs(x-y)>5.1e-6*abs(x):cells.append(dict(id=int(ident),species=name,final_eoc_X=x,plot_X=y,ratio=y/x if x!=0 else np.nan))
        rows.append(dict(id=int(ident),mismatched_species=int(e['plot_mismatches'][er]),final_sumX=e['sumX'][er],final_Ye=e['Ye'][er]))
    csvsave('tables/eoc_plot_mismatch_cells.csv',cells);csvsave('tables/eoc_plot_mismatch_particles.csv',rows)
    special=pd.DataFrame(cells);special=special[special.id==7753]
    save('provenance/particle_7753.json',dict(compared_mismatching_species=len(special),plot_over_final_ratio=percentiles(special.ratio),Fe54=special[special.species=='54Fe'].to_dict('records')))

def runtime_provenance():
    nm,nt,ni=table('new');nu,nr=unique_rows(ni);hydro=np.load(ROOT/'tables/hydro_histories.npy',mmap_mode='r');hm=json.loads((ROOT/'provenance/hydro_history_sources.json').read_text());iteration=np.array([x['iteration'] for x in hm['snapshots']]);rows=[]
    for p in sorted((ROOT/'provenance/raw_archive_prefixes').glob('*/*/xtime.txt')):
        ident=int(p.parent.name);times=[];headers=[];last=None;first=None;grouped={}
        try:
            with p.open() as f:
                for h,x in time_blocks(f):
                    if first is None:first=(h,x.copy())
                    last=(h,x);times.append(h['time_s']);headers.append(h);key=format(h['time_s'],'.6g');grouped.setdefault(key,collections.defaultdict(float))
                    for pair,val in x.items():grouped[key][pair]+=val
            h,x=last;idslot=np.searchsorted(nu,ident);plot=nt[nr[idslot]] if idslot<len(nu) and nu[idslot]==ident else None
            bad=groupbad=compared=0
            for j,name in enumerate(nm['columns']):
                kind,Z,A=species(name)
                if kind!='isotope' or name=='n' or plot is None:continue
                if (Z,A) in x:
                    ref=x[(Z,A)];compared+=1;bad+=abs(plot[j]-ref)>5.1e-6*abs(ref)
                ref=grouped[format(h['time_s'],'.6g')].get((Z,A),0);groupbad+=abs(plot[j]-ref)>5.1e-6*abs(ref)
            firsth=first[0];t=np.array(times);T=np.array([q['temp_GK'] for q in headers]);rho=np.array([q['rho_g_cm3_as_written'] for q in headers]);valid=(t>=0)&(t<=iteration[-1]);native=np.interp(t[valid],iteration,hydro[ident,:,1]);nativeT=np.interp(t[valid],iteration,hydro[ident,:,0]/1e9)
            rows.append(dict(id=ident,blocks=len(t),first_time=firsth['time_s'],last_time=h['time_s'],first_T_GK=firsth['temp_GK'],first_rho_cgs_as_written=firsth['rho_g_cm3_as_written'],hydro_initial_T_GK=hydro[ident,0,0]/1e9,hydro_initial_rho_native=hydro[ident,0,1],hydro_initial_rho_cgs=hydro[ident,0,1]*(1.9889e27/6.955e10**3),duplicate_printed_times=len(t)-len(np.unique(t)),six_sig_time_key_collisions=len(t)-len(grouped),sum_available_last_X=sum(x.values()),finite_last_vs_plot_pairs=compared,last_vs_plot_mismatches=int(bad),timekey_aggregate_vs_plot_mismatches=int(groupbad),runtime_prefix_sha256=sha(p),candidate_interpolation_T_matches=int(np.isclose(T[valid],nativeT,rtol=2e-5,atol=0).sum()),candidate_interpolation_rho_matches=int(np.isclose(rho[valid],native,rtol=2e-5,atol=5.1e-6).sum()),candidate_interpolation_points=int(valid.sum())))
            pd.DataFrame(headers).to_csv(ROOT/f'tables/runtime_headers_{ident}.csv',index=False)
        except Exception as ex:rows.append(dict(id=ident,error=repr(ex)))
    csvsave('tables/archive_runtime_validation.csv',rows)
    controls=json.loads((ROOT/'provenance/archived_control_records.json').read_text());crecs=[];paths=set()
    for r in controls:
        c=r['control'];datasets=c['nuclear']['datasets'];crecs.append(dict(id=r['particle'],nuclear_configuration_sha256=hashlib.sha256(json.dumps(datasets,sort_keys=True).encode()).hexdigest(),active_datasets=sum(d.get('active',True) for d in datasets),stop=c['network']['stop'],conditions=c['conditions'],extent=c['network']['extent']))
        for d in datasets:
            if 'path' in d:paths.add(d['path'])
        paths.add(c['network']['extent']['path'])
    files=[]
    for name in sorted(paths):
        try:files.append(dict(**stamp(name),sha256=sha(name)))
        except OSError as ex:files.append(dict(path=name,error=str(ex)))
    save('provenance/current_nuclear_input_fingerprints.json',files);save('provenance/control_comparison.json',crecs)

def final_stability():
    checks=[]
    for p in sorted((ROOT/'provenance/chunks').glob('chunk_*.json')):
        d=json.loads(p.read_text())
        for label in ['A','B']:
            before=d[label+'_source'];now=stamp(before['path']);changed=any(now[k]!=before[k] for k in ['bytes','mtime_ns','ctime_ns','inode']);checks.append(dict(path=before['path'],metadata_changed=changed,content_rechecked=changed,content_matches=sha(before['path'])==before['sha256'] if changed else None))
    for label in ['old','new']:
        d=json.loads((ROOT/f'provenance/{label}_scan.json').read_text());checks.append(dict(path=d['source']['path'],content_rechecked=True,content_matches=sha(d['source']['path'])==d['sha256']))
    save('provenance/source_stability.json',checks)
    assert all(r.get('content_matches') is not False for r in checks),'Input content changed; results need review'
    log('Source stability verified against original receipts')

def reduce_inventory():
    rows=[];exact=[];dirs=[]
    for df in pd.read_csv(ROOT/'provenance/additional_50AM_candidates.csv',chunksize=100000):
        for row in df.itertuples(index=False):
            p=Path(row.path);name=p.name.lower()
            if p.suffix.lower() in ['.py','.sh','.ipynb','.f','.f90','.par','.cfg'] and row.type=='file':rows.append(dict(path=str(p),bytes=row.bytes))
            if name in ['50amplotting.out','50am_plotting.out','50am_plotting_all_isotopes.out','final_chunks.out']:exact.append(dict(path=str(p),bytes=row.bytes))
            if row.type=='directory' and any(s in name for s in ['50am','final_chunks']):dirs.append(dict(path=str(p)))
    csvsave('provenance/focused_source_candidates.csv',rows);csvsave('provenance/all_group_dataset_paths.csv',exact);csvsave('provenance/all_50am_directories.csv',dirs)

if __name__=='__main__':
    coverage();eoc_database();runtime_provenance();final_stability();reduce_inventory();log('Full-population attribution complete')
