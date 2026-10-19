#!/usr/bin/env python3
"""Validate legacy trajectories, numerical interpretation, and group source leads."""
import os,re,json,collections
from pathlib import Path
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def figsave(fig,name):
    dest=ROOT/'figures/validation';dest.mkdir(exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(dest/(name+'.'+ext),dpi=240,bbox_inches='tight')
    plt.close(fig)

def legacy_hydro():
    h=np.load(ROOT/'tables/hydro_histories.npy',mmap_mode='r');m=json.loads((ROOT/'provenance/hydro_history_sources.json').read_text());times=np.array([q['time_code']*100 for q in m['snapshots']]);rows=[];counts=collections.Counter();badexamples=[]
    for p in sorted((ROOT/'tables/raw_old').glob('out*/histories.npz')):
        d=np.load(p);v=d['time_temp_rho'];ids=np.repeat(d['id'],np.diff(d['offset']));loc=np.searchsorted(times,v[:,0]);left=np.maximum(loc-1,0);right=np.minimum(loc,len(times)-1);j=np.where(abs(v[:,0]-times[left])<=abs(v[:,0]-times[right]),left,right)
        raw=h[ids,j];tgood=np.isclose(v[:,0],times[j],rtol=1e-6,atol=1e-9);Tgood=np.isclose(v[:,1],raw[:,0],rtol=1e-10,atol=0);rgood=np.isclose(v[:,2],raw[:,1]*1.9889e27/6.955e10**3,rtol=1e-10,atol=0);finite=np.isfinite(v).all(axis=1)
        counts.update(points=len(v),physical_time_matches=int(tgood.sum()),temperature_matches=int((tgood&Tgood).sum()),density_cgs_matches=int((tgood&rgood).sum()),all_three_match=int((tgood&Tgood&rgood).sum()),nonfinite_points=int((~finite).sum()))
        for k,ident in enumerate(d['id']):
            a,b=d['offset'][k:k+2];rows.append(dict(id=int(ident),points=int(b-a),matched_times=int(tgood[a:b].sum()),matched_temperatures=int((tgood&Tgood)[a:b].sum()),matched_cgs_densities=int((tgood&rgood)[a:b].sum()),nonfinite=int((~finite[a:b]).sum())))
        for k in np.flatnonzero(~(tgood&Tgood&rgood))[:10]:badexamples.append(dict(id=int(ids[k]),legacy_time=v[k,0],hydro_time=times[j[k]],legacy_temperature=v[k,1],hydro_temperature=raw[k,0],legacy_density=v[k,2],hydro_density_cgs=raw[k,1]*1.9889e27/6.955e10**3))
    csvsave('tables/legacy_history_vs_hydro.csv',rows);save('provenance/legacy_history_vs_hydro.json',dict(counts=dict(counts),examples=badexamples,time_rtol=1e-6,time_atol_s=1e-9,temperature_and_density_rtol=1e-10,temperature_and_density_atol=0))

def supplement():
    nm,nt,ni=table('new');nu,nr=unique_rows(ni);cols=nm['columns'];ids,ro,rn=matched()
    initial=pd.read_csv(ROOT/'tables/initial_composition_population.csv.gz');trajectory=pd.read_csv(ROOT/'tables/trajectory_population.csv.gz');init=initial.set_index('id').reindex(nu);traj=trajectory.set_index('id').reindex(nu)
    e=np.load(ROOT/'tables/eoc_particle_validation.npz');common,ei,ii=np.intersect1d(e['id'],nu,return_indices=True);initial_ye=init.Ye_after_PRISM_initial_normalization.to_numpy()[ii];yechange=e['Ye'][ei]-initial_ye
    reported=nt[nr,cols.index('Y_{e}')];checks=dict(reported_Ye_zero=int((reported==0).sum()),reported_Ye_negative=int((reported<0).sum()),reported_Ye_over_one=int((reported>1).sum()),eoc_minus_archived_normalized_seed_Ye=percentiles(yechange),initial_vs_final_epoch_caveat='Archived seeds can precede restarts; this is an archive comparison, not proof of per-run input identity')
    save('provenance/ye_interpretation.json',checks)
    # All-species mass-weighted production by independently reconstructed physical cooling time.
    rows=[];mass=nt[nr,cols.index('mass')]/MSUN;tau=traj.hydro_physical_cooling_efold_sampled_s.to_numpy();edges=np.r_[0,np.logspace(-3,6,19),np.inf]
    for j,s in enumerate(cols):
        if species(s)[0]!='isotope':continue
        x=nt[nr,j]
        for lo,hi in zip(edges[:-1],edges[1:]):
            mask=(tau>=lo)&(tau<hi)&np.isfinite(x);rows.append(dict(species=s,cooling_time_lower_s=lo,cooling_time_upper_s=hi,particles=int(mask.sum()),yield_msun=float(np.sum(mass[mask]*x[mask]))))
        mask=~np.isfinite(tau)&np.isfinite(x);rows.append(dict(species=s,cooling_time_lower_s=np.nan,cooling_time_upper_s=np.nan,particles=int(mask.sum()),yield_msun=float(np.sum(mass[mask]*x[mask]))))
    csvsave('tables/all_species_physical_cooling_bins.csv',rows)
    fig,ax=plt.subplots(1,2,figsize=(10,4));ax[0].hist(e['sumX']-1,bins=150,histtype='step');ax[0].set_yscale('log');ax[0].set_xlabel('Archived final EOC sum X − 1');ax[0].set_ylabel('Particles');ok=np.isfinite(initial_ye)&np.isfinite(e['Ye'][ei]);ax[1].hexbin(initial_ye[ok],e['Ye'][ei][ok],gridsize=80,bins='log',mincnt=1);ax[1].set_xlabel('Archived normalized seed Ye');ax[1].set_ylabel('Archived final composition Ye');figsave(fig,'full_eoc_closure_and_seed_ye')
    fig,ax=plt.subplots(figsize=(6,4));ok=np.isfinite(reported[ii]);ax.hexbin(reported[ii][ok],e['Ye'][ei][ok],gridsize=90,bins='log',mincnt=1);ax.set(xlabel='Copied progenitor Ye (zero can mean missing)',ylabel='Archived final composition Ye');figsave(fig,'full_eoc_vs_reported_ye')
    q=pd.read_csv(ROOT/'tables/yield_coverage_decomposition.csv');q=q[(q.raw_covered_new_yield_msun+q.raw_absent_new_yield_msun)>0];fig,ax=plt.subplots(figsize=(max(10,len(q)*.2),4));frac=q.raw_absent_new_yield_msun/(q.raw_absent_new_yield_msun+q.raw_covered_new_yield_msun);ax.bar(np.arange(len(q)),frac);ax.set_xticks(np.arange(len(q)),q.species,rotation=90,fontsize=7);ax.set(ylabel='Fraction of new yield on IDs lacking legacy raw output',ylim=(0,1));figsave(fig,'legacy_coverage_yield_attribution')
    q=pd.read_csv(ROOT/'tables/eoc_mismatches_with_chunk_epoch.csv');fig,ax=plt.subplots(1,2,figsize=(10,4));ax[0].scatter(q.time_last,q.sum_last,s=10);ax[0].set(xscale='log',xlabel='Last chunk time (s, as written)',ylabel='Selected isotope sum X');ax[1].scatter(q.chunk,q.mismatched_species,s=10);ax[1].set(xlabel='Chunk number',ylabel='Discrepant species vs final EOC');figsave(fig,'eoc_mismatches_by_time_and_chunk')
    runtime=pd.read_csv(ROOT/'tables/archive_runtime_validation.csv');fig,ax=plt.subplots(1,2,figsize=(10,4));ax[0].scatter(runtime.hydro_initial_rho_cgs,runtime.first_rho_cgs_as_written);lo=runtime.hydro_initial_rho_cgs.min();hi=runtime.first_rho_cgs_as_written.max();ax[0].plot([lo,hi],[lo,hi],'k--',lw=.7);ax[0].set(xscale='log',yscale='log',xlabel='Hydrodynamic initial density (g cm⁻³)',ylabel='PRISM runtime initial density label (g cm⁻³)');hm=json.loads((ROOT/'provenance/hydro_history_sources.json').read_text());ax[1].plot([s['time_code']*100 for s in hm['snapshots']],[s['iteration'] for s in hm['snapshots']],'.-');ax[1].set(xscale='symlog',yscale='symlog',xlabel='Hydrodynamic physical time (s)',ylabel='Archived trajectory time label (s)');figsave(fig,'trajectory_unit_diagnostics')
    # Real decay chains include the 44Sc intermediate. Missing chains are not tested as closed.
    d=pd.read_csv(ROOT/'tables/decay_chain_diagnostics.csv',dtype={'old_missing':'str','new_missing':'str'});mask=d.chain=='44Ti -> 44Ca';d.loc[mask,'chain']='44Ti -> 44Sc -> 44Ca';d.loc[mask,'old_missing']='44Sc';d.loc[mask,'new_missing']='44Sc';d.loc[mask,'valid_closed_chain_pairs']=0;d.loc[mask,'inconsistent_with_isolated_decay_1e5']=0
    for c in d.columns:
        if c.startswith('chain_delta_'):d.loc[mask,c]=np.nan
    d.to_csv(ROOT/'tables/decay_chain_diagnostics.csv',index=False)

def source_leads():
    leads=[];direct=[];alternates=[];symlinks=[];code_candidates=0
    # Stream the completed metadata inventory; no second recursive filesystem walk.
    for df in pd.read_csv(ROOT/'provenance/all_group_files.csv.gz',chunksize=200000):
        for r in df[df.type.eq('symlink')].itertuples(index=False):symlinks.append(dict(path=r.path,target=os.readlink(r.path)))
        p=df.path.astype(str);mask=df.type.eq('file')&p.str.endswith(('.py','.sh','.ipynb','.f90','.f'))&(df.bytes<2000000)&~p.str.contains(r'/ext/|site-packages|/\.git/')
        for row in df[mask].itertuples(index=False):
            code_candidates+=1;path=Path(row.path)
            try:
                txt=path.read_text(errors='replace');hits=[dict(line=i+1,text=line[:350]) for i,line in enumerate(txt.splitlines()) if re.search(r'50Am_abund|50Am_output|50AmPlotting|50AmP|looperror.*1337',line,re.I)]
                if hits:leads.append(dict(path=str(path),sha256=sha(path),matches=hits[:25]))
            except OSError as ex:leads.append(dict(path=str(path),error=str(ex)))
        mask=p.str.contains(r'/(?:7753|11000)/(?:xtime\.txt|control\.json|x\.txt)$',regex=True)
        direct.extend(df.loc[mask,['path','bytes']].to_dict('records'))
    save('provenance/full_group_source_leads.json',dict(source_files_examined=code_candidates,leads=leads,direct_representative_outputs=direct,symlinks=symlinks))
    alt=Path('/data/pyoung3/gsvance/supernova_data/50Am/analysis/50Am_plotting.out');old=json.loads((ROOT/'provenance/old_scan.json').read_text());digest=sha(alt);save('provenance/alternate_legacy_copy.json',dict(source=stamp(alt),sha256=digest,byte_identical_sha256_to_main=digest==old['sha256']))

if __name__=='__main__':
    if not (ROOT/'provenance/legacy_history_vs_hydro.json').exists():legacy_hydro()
    supplement();log('Final physical consistency checks and figures complete');source_leads();log('Complete group source lead search finished')
