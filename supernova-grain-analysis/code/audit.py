#!/usr/bin/env python3
"""Reproducible bounded data audit. Never runs collaborator scripts or writes sources."""
from pathlib import Path
import argparse, csv, io, json, os, re, tarfile, hashlib, sys
from collections import Counter
from datetime import datetime, timezone

REPO = Path(__file__).resolve().parents[1]
OUT = Path('/scratch/rpraka15/supernova-grain-analysis')
os.environ.setdefault('MPLCONFIGDIR', str(OUT/'cache/matplotlib'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from readers import open_sdf, sdf_header, units, table, time_blocks, request_kind, read_extent

GROUP = Path('/data/pyoung3/groupdata/svierlin')
PRISM = Path('/scratch/svierlin/prism-v02')
LIST = Path('/scratch/rpraka15/List of isotopes_elements_variables.txt')
HYDRO = {
 '50Am': Path('/data/pyoung3/supernova_data/50Am/run3g_50Am6'),
 'jet3b': Path('/data/pyoung3/supernova_data/jet3b/jet3b'),
 'cco2': Path('/data/pyoung3/gsvance/supernova_data/cco2/cco2'),
 'vconv': Path('/data/pyoung3/gsvance/supernova_data/vconv/vconv_snsph'),
}

def setup():
    # Resolve the write boundary, including existing symlinks, before creating anything.
    if not OUT.resolve().is_relative_to(Path('/scratch/rpraka15').resolve()):
        raise ValueError('Output outside user scratch')
    for d in ['tables','plots','debug','cache','reports','code']:
        (OUT/d).mkdir(parents=True,exist_ok=True)

def dump(name, value):
    (OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def csvout(name, rows, fields=None):
    if not rows and fields is None:
        return
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader();w.writerows(rows)

def savefig(name):
    plt.savefig(OUT/'plots'/f'{name}.png',dpi=160,bbox_inches='tight')
    plt.savefig(OUT/'plots'/f'{name}.pdf',bbox_inches='tight')
    plt.close()

def catalog():
    rows,errors=[],[]
    for model,root in HYDRO.items():
        paths=sorted((p for p in root.iterdir() if re.search(r'_sph\.\d+$',p.name)),key=lambda p:int(p.name.rsplit('.',1)[1]))
        if model=='vconv':
            paths += [root.parent/'run3g_50_1milvc2_sph.4200']
        print('Reading snapshot headers',model,len(paths),flush=True)
        for path in paths:
            try:
                h=sdf_header(path); u=units(h['meta'])
                rows.append({'model':model,'path':str(path),'iteration':int(h['meta']['iter']),
                 'time_code':h['meta']['tpos'],'time_s':h['meta']['tpos']*u['time_s'],
                 'units_provenance':u['provenance'],'npart':h['npart'],
                 'record_bytes':h['record_bytes'],'header_bytes':h['offset'],
                 'file_bytes':h['file_bytes'],'field_count':len(h['dtype'].names),
                 'density_factor':u['density_g_cm3'],'time_factor':u['time_s'],
                 'header_sha256':hashlib.sha256(h['header'].encode()).hexdigest()})
            except (OSError,ValueError) as e:
                errors.append({'path':str(path),'error':str(e)})
    csvout('tables/snapshot_catalog.csv',rows);dump('debug/snapshot_errors.json',errors)
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for ax,model in zip(axes.flat,HYDRO):
        r=[r for r in rows if r['model']==model]
        ax.plot([r['iteration'] for r in r],[r['time_s'] for r in r],'.-',ms=3)
        ax.set(xlabel='Snapshot iteration (filename suffix)',ylabel='Time [s]',yscale='symlog',title=model+(' (inferred time unit)' if model in ['50Am','jet3b'] else ' (header time unit)'))
        ax.grid(alpha=.2)
    fig.suptitle('Iteration is not physical time; each header supplies tpos')
    savefig('01_snapshot_times')
    return rows

def hydro_samples(rows):
    diagnostics=[];selected={};positions={}
    for model in HYDRO:
        candidates=[r for r in rows if r['model']==model]
        r=max(candidates,key=lambda r:r['iteration'])
        h,a=open_sdf(r['path']);u=units(h['meta']); ids=np.asarray(a['ident']).copy()
        unique=np.unique(ids)
        idx=np.linspace(0,len(a)-1,min(20000,len(a)),dtype=int)
        s=a[idx].copy()
        selected[model]=(h,a,np.argsort(ids),ids)
        positions[model]=(s,r,u)
        sums=sum(s[f].astype(float) for f in a.dtype.names if re.fullmatch(r'f\d+',f))
        ye_from_x=np.zeros(len(s));bad_species=np.zeros(len(s),dtype=bool)
        species=[]
        for f in a.dtype.names:
            if not re.fullmatch(r'f\d+',f):continue
            i=f[1:]; p='p'+i; n='m'+i if 'm'+i in a.dtype.names else 'n'+i
            if p in a.dtype.names:
                combos=np.unique(np.column_stack([s[p],s[n]]),axis=0)
            else: combos=np.array([[h['meta'][p],h['meta'][n]]],dtype=int)
            for z,nn in combos:species.append({'field':f,'Z':int(z),'N':int(nn),'A':int(z+nn)})
            zz=s[p].astype(float) if p in a.dtype.names else np.full(len(s),h['meta'][p])
            nn=s[n].astype(float) if n in a.dtype.names else np.full(len(s),h['meta'][n])
            aa=zz+nn;bad_species|=(aa<=0)|(zz<0)|(nn<0)
            ye_from_x+=np.divide(zz,aa,out=np.zeros_like(zz),where=aa>0)*s[f]
        dump(f'debug/{model}_hydro_species_sample.json',species)
        (OUT/'debug'/f'{model}_selected_sdf_header.txt').write_text(h['header'])
        diag={'model':model,'path':r['path'],'npart':len(a),'unique_ids':len(unique),
              'duplicate_ids':len(a)-len(unique),'id_min':int(ids.min()),'id_max':int(ids.max()),
              'sample_size':len(s),'mass_nonpositive_sample':int(np.sum(s['mass']<=0)),
              'rho_nonpositive_sample':int(np.sum(s['rho']<=0)),
              'sample_sum_hydro_X_min':float(sums.min()),'sample_sum_hydro_X_max':float(sums.max()),
              'sample_Ye_min':float(s['Y_el'].min()),'sample_Ye_max':float(s['Y_el'].max()),
              'sample_Ye_zero_count':int(np.sum(s['Y_el']==0)),
              'sample_Ye_from_X_min':float(ye_from_x.min()),'sample_Ye_from_X_max':float(ye_from_x.max()),
              'sample_invalid_species_records':int(bad_species.sum()),
              'nonfinite_physical_sample':sum(int(np.sum(~np.isfinite(s[f]))) for f in ['x','y','z','rho','temp','mass','Y_el']),
              'nuclear_slots':len([f for f in a.dtype.names if re.fullmatch(r'f\d+',f)])}
        diagnostics.append(diag);print('Hydro validated',model,diag['npart'],'IDs',flush=True)
    csvout('tables/hydro_diagnostics.csv',diagnostics)
    fig,axes=plt.subplots(2,2,figsize=(11,9),layout='constrained')
    for ax,(model,(s,r,u)) in zip(axes.flat,positions.items()):
        c=np.log10(np.maximum(s['temp'],1.e-30))
        sc=ax.scatter(s['x'],s['z'],c=c,s=1,cmap='inferno',rasterized=True)
        ax.set(xlabel='x [code length]',ylabel='z [code length]',title=f'{model}: iteration {r["iteration"]:,}, n={len(s):,}')
        ax.set_aspect('equal',adjustable='datalim');fig.colorbar(sc,ax=ax,label='log10 temperature [K]')
    fig.suptitle('Hydro particle samples at different final available snapshots — not a common epoch')
    savefig('02_particle_geometry')
    fig,axes=plt.subplots(2,4,figsize=(15,7),layout='constrained')
    for col,(model,(s,r,u)) in enumerate(positions.items()):
        sums=sum(s[f].astype(float) for f in s.dtype.names if re.fullmatch(r'f\d+',f))
        axes[0,col].hist(sums,bins=50,color='#476eaa')
        axes[0,col].axvline(1,color='black',ls='--',lw=1)
        axes[0,col].set(title=model,xlabel='Sum hydro nuclear slots',ylabel='Sample particles',yscale='log')
        axes[1,col].hist(s['Y_el'],bins=40,color='#397e69')
        axes[1,col].set(xlabel='Stored Y_el',ylabel='Sample particles',yscale='log')
    fig.suptitle('Hydro composition diagnostics: zero Y_el and nonunit sums require interpretation')
    savefig('06_hydro_composition_flags')
    return selected

class LimitedReader:
    """Cap decompressed tar bytes: sampling must never traverse a multi-TB archive."""
    def __init__(self,f,limit):self.f=f;self.left=limit
    def read(self,n=-1):
        if self.left<=0:raise ValueError('Archive prefix byte budget reached')
        n=min(n if n>=0 else self.left,self.left)
        b=self.f.read(n);self.left-=len(b);return b

def archive_examples():
    import gzip
    examples=[];statuses=[]
    for batch in [2,3,4,5,6,7,10]:
        p=GROUP/'50Am_data'/f'eoc_batch_{batch:04d}.tar.gz'
        files={};first=None;count=0
        try:
            with gzip.open(p,'rb') as gz, tarfile.open(fileobj=LimitedReader(gz,64*1024**2),mode='r|') as t:
                for member in t:
                    count+=1
                    parts=Path(member.name).parts
                    if not parts or not parts[0].isdigit():continue
                    if first is None:first=parts[0]
                    if parts[0]!=first:break
                    if member.isfile() and parts[-1] in ['control.json','x.txt','xz.txt','y.txt','xtime.txt'] and member.size<24*1024**2:
                        files[parts[-1]]=t.extractfile(member).read().decode('ascii')
                    if {'control.json','x.txt','xz.txt','y.txt','xtime.txt'}<=files.keys():break
            status='prefix_only_not_full_integrity_check'
        except (OSError,EOFError,ValueError,tarfile.TarError) as e:status=str(e)
        statuses.append({'path':str(p),'first_id':first,'members_seen':count,'files_read':list(files),'status':status})
        if 'control.json' in files and 'x.txt' in files:
            examples.append(('50Am',int(first),str(p)+'::'+first,files))
        print('Archive prefix',batch,first,list(files),flush=True)
    dump('debug/archive_prefix_checks.json',statuses)
    return examples

def sample_dirs():
    # Deterministic direct probes across the ID range; no claim of random representativeness.
    probes=[0,23,1000,10076,50000,100000,102027,118326,191136,200000,251189,
            300000,400000,500000,600000,700000,800000,900000,950000,999999]
    samples=[];errors=[]
    census_samples={}
    sample_path=OUT/'inventory/sample_ids.json'
    if sample_path.exists():census_samples=json.loads(sample_path.read_text())
    for model in ['jet3b','cco2','vconv']:
        root=GROUP/f'{model}_data'/f'{model}_outputs'
        model_probes=sorted(set(probes+[int(p) for p in census_samples.get(model,[])]))
        for pid in model_probes:
            p=root/str(pid)
            if not p.is_dir():continue
            files={}
            for name in ['control.json','x.txt','xz.txt','y.txt','xtime.txt']:
                try:
                    f=p/name
                    if f.stat().st_size>32*1024**2:
                        errors.append({'path':str(f),'error':'over 32 MiB sample cap'});continue
                    files[name]=f.read_text()
                except OSError as e:errors.append({'path':str(p/name),'error':str(e)})
            if files:samples.append((model,pid,str(p),files))
        print('PRISM samples',model,sum(s[0]==model for s in samples),flush=True)
    dump('debug/sample_read_errors.json',errors)
    return samples

def prism_samples(samples,hydro):
    labels=[s.strip() for s in LIST.read_text().splitlines() if s.strip()]
    (OUT/'debug/requested_fields.txt').write_text('\n'.join(labels)+'\n')
    requested=[];diagnostics=[];trajectories={};headers=[];joins=[];errors=[]
    for model,pid,path,files in samples:
        try:
            config=json.loads(files['control.json']);x=table(files.get('x.txt',''),2)
            xz=table(files.get('xz.txt',''),1);y=table(files.get('y.txt',''),2)
            extent=read_extent(config['network']['extent']['path'])
            stop=config['network']['stop'].get('time')
            history=[];last={}
            if 'xtime.txt' in files:
                for head,last in time_blocks(io.StringIO(files['xtime.txt'])):history.append(head)
            if history:
                for head in history:headers.append({'model':model,'particle_id':pid,**head})
            end=history[-1]['time_s'] if history else None
            endpoint_match=bool(end is not None and stop is not None and np.isclose(end,stop,rtol=1e-5,atol=1e-10))
            maxdelta=max((abs(v-k[1]*y[k]) for k,v in x.items() if k in y),default=None)
            zsum={}
            for (z,a),v in x.items():zsum[(z,)]=zsum.get((z,),0)+v
            zdelta=max((abs(xz[z]-v) for z,v in zsum.items() if z in xz),default=None)
            trajectory=config['conditions']['trajectory']['path']
            traj=None
            try:
                traj=np.loadtxt(trajectory,skiprows=1,ndmin=2)
                if traj.shape[1]!=3 or not np.all(np.isfinite(traj)) or np.any(traj[:,1:]<=0):
                    raise ValueError('Invalid trajectory shape or values')
                trajectories.setdefault(model,[]).append((pid,traj,trajectory))
            except (OSError,ValueError) as e:
                traj=None
                errors.append({'path':trajectory,'error':str(e)})
            diag={'model':model,'particle_id':pid,'source':path,'x_present':'x.txt' in files,
                'xz_present':'xz.txt' in files,'x_rows':len(x),'sum_X':sum(x.values()) if x else None,
                'sum_Xz':sum(xz.values()) if xz else None,'X_minus_A_Y_max':maxdelta,'element_sum_delta_max':zdelta,
                'Y_e_from_saved_X':sum(z/a*v for (z,a),v in x.items()) if x else None,
                'configured_stop_s':stop,'last_saved_time_s':end,'endpoint_matches_config':endpoint_match,
                'endpoint_status':('not_checked_history_missing_or_capped' if end is None else
                                   'matches_config' if endpoint_match else 'does_not_match_config'),
                'time_blocks':len(history),'time_nondecreasing':bool(all(history[i]['time_s']>=history[i-1]['time_s'] for i in range(1,len(history)))),
                'max_final_last_X_delta':max((abs(x[k]-last[k]) for k in x.keys() & last.keys()),default=None),
                'trajectory_source':trajectory,'trajectory_points':len(traj) if traj is not None else None,
                'peak_T_K_in_input':float(traj[:,1].max()*1e9) if traj is not None else None,
                'peak_rho_in_input_as_written':float(traj[:,2].max()) if traj is not None else None,
                'trajectory_units_status':'unvalidated_see_provenance_report'}
            # Demonstrate an ID join, keeping hydro epoch separate from abundance epoch.
            h,a,order,ids=hydro[model]
            sorted_ids=ids[order]
            j=int(np.searchsorted(sorted_ids,pid));right=int(np.searchsorted(sorted_ids,pid,side='right'))
            index=int(order[j]) if right-j==1 else None
            joined={'model':model,'particle_id':pid,'hydro_found':index is not None,'hydro_id_matches':right-j,'hydro_source':h['path'],
                    'hydro_time_code':h['meta']['tpos'],'PRISM_stop_as_written_s':stop,
                    'epoch_status':'different_epochs_not_a_synchronized_science_table','PRISM_source':path}
            if index is not None:
                for f in ['x','y','z','vx','vy','vz','ax','ay','az','mass','h','rho','Y_el']:
                    joined[f+'_hydro_code' if f!='Y_el' else 'Y_el_hydro']=float(a[f][index])
            joins.append(joined)
            for label in labels:
                kind,key=request_kind(label);value=None
                if kind=='isotope':
                    value=x.get(key)
                    status=('missing_source' if 'x.txt' not in files else 'empty_output' if not x else
                            'stored' if key in x else 'outside_network' if key not in extent else 'omitted_sparse_or_unverified')
                    source=path+'/x.txt';meaning='final_mass_fraction';threshold=config.get('output',{}).get('x',{}).get('lower_bound',1e-15)
                elif kind=='element':
                    value=xz.get(key)
                    status=('missing_source' if 'xz.txt' not in files else 'empty_output' if not xz else
                            'stored' if key in xz else 'omitted_sparse_or_unverified')
                    source=path+'/xz.txt';meaning='final_element_mass_fraction';threshold=config.get('output',{}).get('xZ',{}).get('lower_bound',1e-15)
                else:
                    threshold=None
                    if label in ['Peak temp','Peak density']:
                        value=diag['peak_T_K_in_input' if label=='Peak temp' else 'peak_rho_in_input_as_written']
                        status='derived_input_peak_units_unvalidated' if value is not None else 'missing_input'
                        source=trajectory;meaning='maximum_of_saved_input_trajectory_not_continuous_peak'
                    else:
                        status='hydro_source_located_choose_epoch';source=h['path'];meaning='snapshot_field'
                        if label=='Id':value=pid;status='stored';meaning='particle_identifier'
                requested.append({'model':model,'particle_id':pid,'request':label,'kind':kind,
                    'Z':key[0] if key else None,'A':key[1] if key and len(key)==2 else None,
                    'value':value,'status':status,'source':source,'meaning':meaning,
                    'PRISM_stop_as_written_s':stop,'configured_output_threshold':threshold})
            diagnostics.append(diag)
        except (OSError,ValueError,KeyError) as e:
            errors.append({'path':path,'error':str(e)})
    csvout('tables/requested_values_sample.csv',requested)
    csvout('tables/prism_diagnostics.csv',diagnostics)
    csvout('tables/prism_time_headers_sample.csv',headers)
    csvout('tables/particle_id_join_demo.csv',joins)
    dump('debug/prism_errors.json',errors)
    mappings=[]
    for model in HYDRO:
        for label in labels:
            subset=[r for r in requested if r['model']==model and r['request']==label]
            kind,key=request_kind(label)
            counts=Counter(r['status'] for r in subset)
            source={'isotope':'<particle>/x.txt: Z,A,X','element':'<particle>/xz.txt: Z,X(Z)',
                    'physical':'SDF chosen snapshot; peak quantities from trajectory'}[kind]
            field={'Id':'ident','Density':'rho','Y_{e}':'Y_el','Peak temp':'max trajectory T * 1e9',
                   'Peak density':'max trajectory rho (conversion disputed)'}.get(label,label if kind=='physical' else '')
            mappings.append({'model':model,'request':label,'kind':kind,'Z':key[0] if key else None,
                'A':key[1] if key and len(key)==2 else None,'source_schema':source,'hydro_field_or_derivation':field,
                'sample_particles':len(subset),'stored_samples':counts.get('stored',0),
                'sparse_missing_samples':counts.get('omitted_sparse_or_unverified',0),
                'outside_network_samples':counts.get('outside_network',0),
                'status_counts':json.dumps(counts)})
    csvout('tables/requested_field_map.csv',mappings)
    return diagnostics,requested,trajectories

def plot_prism(diagnostics,requested,trajectories):
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for ax,model in zip(axes.flat,HYDRO):
        rows=[r for r in requested if r['model']==model and r['kind'] in ['isotope','element']]
        labels=['Explicitly stored','Sparse omission / unverified','Outside network','Missing / empty source']
        n=[sum(r['status']==s for r in rows) for s in ['stored','omitted_sparse_or_unverified','outside_network']]
        n.append(sum(r['status'] in ['missing_source','empty_output'] for r in rows))
        ax.barh(labels,n,color=['#24876b','#e6a23c','#b65050','#777777'])
        ax.set(title=f'{model}: sampled requested composition cells',xlabel='Particle × requested field pairs')
    fig.suptitle('Blank isotope cells must not silently become zeros; sample sizes differ')
    savefig('03_requested_composition_coverage')
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for i,model in enumerate(HYDRO):
        ds=[r for r in diagnostics if r['model']==model]
        if not ds:continue
        valid=[r for r in ds if r['sum_X'] is not None]
        axes[0].scatter([i]*len(valid),[r['sum_X']-1 for r in valid],s=20,label=model)
        axes[1].scatter([i]*len(ds),[r['last_saved_time_s'] or np.nan for r in ds],s=20)
    for ax in axes:
        ax.set_xticks(range(4),HYDRO);ax.grid(alpha=.2)
    axes[0].set(yscale='symlog',ylabel='Sum of saved final X minus 1',title='Mass-fraction normalization (unweighted)')
    axes[0].set_yscale('symlog',linthresh=1e-13)
    axes[1].set(yscale='log',ylabel='Last PRISM time [s, as written]',title='Final does not mean one shared epoch')
    savefig('04_prism_sanity')
    fig,axes=plt.subplots(2,4,figsize=(15,7),layout='constrained')
    for col,model in enumerate(HYDRO):
        for pid,t,p in trajectories.get(model,[])[:4]:
            axes[0,col].plot(t[:,0],t[:,1],label=str(pid));axes[1,col].plot(t[:,0],t[:,2])
        axes[0,col].set(title=model,yscale='log',ylabel='Temperature [GK]')
        axes[1,col].set(yscale='log',ylabel='Density [as-written g/cm³]',xlabel='Input time [as-written s]')
        for row in range(2):axes[row,col].set_xscale('symlog',linthresh=1e-3);axes[row,col].grid(alpha=.2)
        if trajectories.get(model):axes[0,col].legend(title='Particle ID',fontsize=7)
    fig.suptitle('Existing PRISM inputs: time/density provenance is disputed (see report)')
    savefig('05_input_trajectories')

def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    setup();rows=catalog();hydro=hydro_samples(rows)
    samples=archive_examples()+sample_dirs()
    diag,requested,trajectories=prism_samples(samples,hydro)
    plot_prism(diag,requested,trajectories)
    error_counts={p.name:len(json.loads(p.read_text())) for p in (OUT/'debug').glob('*errors.json')}
    dump('debug/run_summary.json',{'completed_utc':datetime.now(timezone.utc).isoformat(),
         'status':'completed_with_read_errors' if any(error_counts.values()) else 'completed_bounded_audit',
         'error_counts':error_counts,
         'snapshot_headers_validated':len(rows),'PRISM_particles_sampled':dict(Counter(r['model'] for r in diag)),
         'requested_fields':len(LIST.read_text().splitlines()),'sample_values_rows':len(requested),
         'scope':'All named snapshot headers; all IDs in four selected snapshots; bounded PRISM/archives samples. Not full archive integrity or all-particle abundance census.',
         'python':sys.version,'numpy':np.__version__,'matplotlib':matplotlib.__version__})
    print('Audit outputs:',OUT,flush=True)

if __name__=='__main__':main()
