#!/usr/bin/env python3
"""Audit all archived corrected AM trajectories against raw hydro and legacy histories.

Archive times are tested as printed and separately as snapshot identifiers. No
correction is applied to any source. A match to an archive is not proof of the
exact historical input revision used by every PRISM run.
"""
import sys,re,json,zipfile,hashlib
from pathlib import Path
import numpy as np
from common import *
sys.path.insert(0,'/scratch/rpraka15/supernova-grain-analysis/code')
from readers import open_sdf
HYDRO=Path('/data/pyoung3/supernova_data/50Am/run3g_50Am6')
RHO_CF=1.9889e27/6.955e10**3

def prepare_hydro():
    cache=ROOT/'tables/hydro_histories.npy';meta=ROOT/'provenance/hydro_history_sources.json'
    if cache.exists() and meta.exists():return json.loads(meta.read_text())
    paths=sorted(p for p in HYDRO.glob('run3g_50Am6_sph.*') if re.fullmatch(r'run3g_50Am6_sph\.\d+',p.name))
    header,raw=open_sdf(paths[-1]);limit=int(raw['ident'].max())+1;del raw
    history=np.lib.format.open_memmap(cache,mode='w+',dtype='f8',shape=(limit,len(paths),2));history[:]=np.nan;records=[]
    for i,p in enumerate(paths):
        before=stamp(p);h,raw=open_sdf(p);ids=np.asarray(raw['ident'],dtype='i8');assert ids.max()<limit and len(np.unique(ids))==len(ids)
        history[ids,i,0]=raw['temp'];history[ids,i,1]=raw['rho'];assert stamp(p)==before
        records.append(dict(source=before,iteration=int(h['meta']['iter']),time_code=float(h['meta']['tpos']),header_sha256=hashlib.sha256(h['header'].encode()).hexdigest(),rows=len(ids)));del raw;log(f'Hydro history {i+1}/{len(paths)}')
    history.flush();save('provenance/hydro_history_sources.json',dict(snapshots=records,shape=list(history.shape),temperature_units='K',density_units='code',time_conversion_s=100,conversion_provenance='/data/pyoung3/gsvance/data_mining/sn_utils.py; 50Am headers omit explicit factors'))
    return json.loads(meta.read_text())

def tau(t,v):
    if len(v)<2 or not np.isfinite(v).all():return np.nan
    k=int(np.argmax(v));later=np.flatnonzero((np.arange(len(v))>k)&(v<=v[k]/np.e))
    return float(t[later[0]]-t[k]) if len(later) else np.nan

def relative_max(a,b):
    valid=np.isfinite(a)&np.isfinite(b)&(b!=0)
    if not valid.any():return np.nan
    with np.errstate(over='ignore'):
        return float(np.max(abs(a[valid]-b[valid])/abs(b[valid])))

def main():
    hm=prepare_hydro();hydro=np.load(ROOT/'tables/hydro_histories.npy',mmap_mode='r');it=np.array([s['iteration'] for s in hm['snapshots']]);physical_time=np.array([s['time_code']*100 for s in hm['snapshots']])
    old={};old_arrays=[]
    for p in sorted((ROOT/'tables/raw_old').glob('out*/histories.npz')):
        with np.load(p) as d:
            h=d['time_temp_rho'];off=d['offset'];ids=d['id'];old_arrays.append(h)
            for k,ident in enumerate(ids):old[int(ident)]=(len(old_arrays)-1,int(off[k]),int(off[k+1]))
    source=Path('/scratch/svierlin/prism-v02/input/conditions/50Am_output_tempcorr.zip');before=stamp(source);records=[];errors=[];member_records=[]
    with zipfile.ZipFile(source) as z:
        for member in z.infolist():
            if '__MACOSX' in member.filename:continue
            m=re.fullmatch(r'(\d+)_trajectory\.txt',Path(member.filename).name)
            if not m:continue
            ident=int(m[1])
            try:
                data=z.read(member);v=np.fromstring(data.split(b'\n',1)[1].replace(b'D',b'E'),sep=' ').reshape(-1,3);t=v[:,0];T=v[:,1];rho=v[:,2]
                if not np.isfinite(v).all() or (T<0).any() or (rho<0).any():raise ValueError('invalid trajectory values')
                loc=np.searchsorted(it,t);safe=np.minimum(loc,len(it)-1);mapped=(loc<len(it))&(it[safe]==t)
                r=dict(id=ident,points=len(t),time_first=t[0],time_last=t[-1],time_nonincreasing=int((np.diff(t)<=0).sum()),peak_temp_GK=float(T.max()),peak_density_as_labelled=float(rho.max()),cooling_efold_sampled_s=tau(t,T),density_efold_sampled_s=tau(t,rho),snapshot_index_matches=int(mapped.sum()),hydro_temperature_matches=0,hydro_native_density_matches=0,hydro_CGS_density_matches=0,hydro_physical_time_matches=0,old_history_available=ident in old)
                if ident<len(hydro) and mapped.any():
                    raw=hydro[ident,safe[mapped],:];tv=T[mapped];rv=rho[mapped];finite=np.isfinite(raw).all(axis=1)
                    r['hydro_points_available']=int(finite.sum());r['hydro_temperature_matches']=int(np.isclose(tv,raw[:,0]/1e9,rtol=5.1e-6,atol=0).sum());r['hydro_native_density_matches']=int(np.isclose(rv,raw[:,1],rtol=5.1e-6,atol=0).sum());r['hydro_CGS_density_matches']=int(np.isclose(rv,raw[:,1]*RHO_CF,rtol=5.1e-6,atol=0).sum());r['hydro_physical_time_matches']=int(np.isclose(t[mapped],physical_time[safe[mapped]],rtol=5.1e-6,atol=0).sum())
                    # Source rho tokens have five decimal places after a prior six-significant-digit export.
                    r['hydro_native_density_matches_export_precision']=int(np.isclose(rv,raw[:,1],rtol=5.1e-6,atol=5.1e-6).sum())
                    r['hydro_physical_cooling_efold_sampled_s']=tau(physical_time[safe[mapped]],raw[:,0])
                if ident in old:
                    b,lo,hi=old[ident];h=old_arrays[b][lo:hi];r.update(old_points=len(h),old_time_last=h[-1,0],old_peak_temp_GK=h[:,1].max()/1e9,old_peak_density=h[:,2].max(),old_cooling_efold_sampled_s=tau(h[:,0],h[:,1]))
                    overlap=(h[:,0]>=t[0])&(h[:,0]<=t[-1]);ho=h[overlap]
                    if len(ho) and (np.diff(t)>0).all():
                        ti=np.interp(ho[:,0],t,T)*1e9;ri=np.interp(ho[:,0],t,rho)
                        r.update(old_overlapping_points=len(ho),temperature_points_agree_as_labelled=int(np.isclose(ti,ho[:,1],rtol=1e-4,atol=0).sum()),density_points_agree_as_labelled=int(np.isclose(ri,ho[:,2],rtol=1e-4,atol=0).sum()),max_relative_T_difference=relative_max(ti,ho[:,1]),max_relative_rho_difference=relative_max(ri,ho[:,2]),old_zero_T_denominators=int((ho[:,1]==0).sum()),old_zero_rho_denominators=int((ho[:,2]==0).sum()))
                records.append(r);member_records.append(dict(id=ident,member=member.filename,sha256=hashlib.sha256(data).hexdigest()))
                if len(records)%50000==0:log(f'Validated {len(records):,} trajectory files')
            except Exception as e:errors.append(dict(id=ident,member=member.filename,error=repr(e)))
    assert stamp(source)==before
    df=pd.DataFrame(records);df.to_csv(ROOT/'tables/trajectory_population.csv.gz',index=False,compression='gzip');csvsave('tables/trajectory_errors.csv',errors);csvsave('provenance/trajectory_member_hashes.csv',member_records)
    summary={k:int(df[k].fillna(0).sum()) for k in ['points','snapshot_index_matches','hydro_points_available','hydro_temperature_matches','hydro_native_density_matches','hydro_CGS_density_matches','hydro_physical_time_matches','old_overlapping_points','temperature_points_agree_as_labelled','density_points_agree_as_labelled'] if k in df}
    summary.update({k:int(df[k].fillna(0).sum()) for k in ['hydro_native_density_matches_export_precision','old_zero_T_denominators','old_zero_rho_denominators']})
    save('provenance/trajectory_summary.json',dict(source=before,sha256=sha(source),trajectories=len(df),errors=len(errors),unique_ids=int(df.id.nunique()),old_histories_matched=int(df.old_history_available.sum()),population_counts=summary,density_correction_factor=RHO_CF,physical_time_last_s=physical_time[-1],printed_time_last=percentiles(df.time_last),temperature_interpolation_test='Linear interpolation of archived corrected trajectory at legacy input times; diagnostic only, not an assertion about either solver interpolation.',historical_input_revision_caveat='The archived controls point at a mutable 50Am_output path. The corrected zip is an available candidate input, not a verified copy of every historical run input.'))
    log('Full candidate trajectory census complete')
if __name__=='__main__':main()
