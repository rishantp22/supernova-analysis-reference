#!/usr/bin/env python3
"""All-particle initial-composition closure and Ye census; preserve input precision."""
import zipfile,re,hashlib
from pathlib import Path
from collections import Counter,defaultdict
from common import *
def main():
    p=Path('/scratch/svierlin/prism-v02/input/conditions/50Am_abund.zip');before=stamp(p);records=[];counts=defaultdict(Counter);errors=[];hashes=[]
    nm,nt,ni=table('new');nu,nr=unique_rows(ni);idrow={int(i):int(r) for i,r in zip(nu,nr)};cols=nm['columns']
    with zipfile.ZipFile(p) as z:
        for member in z.infolist():
            if '__MACOSX' in member.filename:continue
            m=re.fullmatch(r'(\d+)_abund\.txt',Path(member.filename).name)
            if not m:continue
            ident=int(m[1])
            try:
                data=z.read(member);v=np.fromstring(data.replace(b'D',b'E'),sep=' ').reshape(-1,3);Z=v[:,0].astype('i4');A=v[:,1].astype('i4');x=v[:,2];code=Z*512+A
                if not np.isfinite(v).all() or (x<0).any() or len(np.unique(code))!=len(code):raise ValueError('Nonfinite/negative abundance or duplicate isotope')
                sx=x.sum();ye=np.sum(x*Z/A);r=dict(id=ident,species_rows=len(x),sumX=sx,Ye_unnormalized=ye,Ye_after_PRISM_initial_normalization=ye/sx if sx>0 else np.nan)
                if ident in idrow:
                    row=idrow[ident];r.update(reported_progenitor_Ye=float(nt[row,cols.index('Y_{e}')]),legacy_peak_temperature_K=float(nt[row,cols.index('peak temp')]))
                records.append(r)
                for c,val in zip(code,x):counts[int(c)]['rows']+=1;counts[int(c)]['positive']+=int(val>0)
                hashes.append(dict(id=ident,member=member.filename,sha256=hashlib.sha256(data).hexdigest()))
            except Exception as e:errors.append(dict(id=ident,member=member.filename,error=repr(e)))
            if len(records)%50000==0:log(f'Initial compositions: {len(records):,}')
    assert stamp(p)==before
    df=pd.DataFrame(records);df.to_csv(ROOT/'tables/initial_composition_population.csv.gz',index=False,compression='gzip');csvsave('tables/initial_composition_errors.csv',errors);csvsave('provenance/initial_composition_member_hashes.csv',hashes);csvsave('tables/initial_species_inventory.csv',[dict(Z=c//512,A=c%512,**v) for c,v in sorted(counts.items())])
    diff=df.Ye_after_PRISM_initial_normalization-df.reported_progenitor_Ye
    save('provenance/initial_composition_summary.json',dict(source=before,sha256=sha(p),particles=len(df),unique_ids=int(df.id.nunique()),errors=len(errors),sumX=percentiles(df.sumX),Ye_normalized=percentiles(df.Ye_after_PRISM_initial_normalization),normalized_initial_minus_reported_Ye=percentiles(diff),normalization='Computed only as a diagnostic of the normalization explicitly implemented by PRISM progression.f90; no source file changed.'))
    log('Initial composition census complete')
if __name__=='__main__':main()
