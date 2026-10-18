#!/usr/bin/env python3
"""Build an offline visual guide, interactive field lookup, and scratch code copy."""
from pathlib import Path
import csv,json,html,shutil,hashlib
from collections import Counter

REPO=Path(__file__).resolve().parents[1]
OUT=Path('/scratch/rpraka15/supernova-grain-analysis')

def readcsv(name):
    with (OUT/name).open() as f:return list(csv.DictReader(f))

def main():
    if not OUT.resolve().is_relative_to(Path('/scratch/rpraka15').resolve()):
        raise ValueError('Output boundary violated')
    summary=json.loads((OUT/'debug/run_summary.json').read_text())
    fields=readcsv('tables/requested_field_map.csv')
    values=readcsv('tables/requested_values_sample.csv')
    hydro=readcsv('tables/hydro_diagnostics.csv')
    inventory=json.loads((OUT/'inventory/inventory.json').read_text())
    modelrows=''
    for row in hydro:
        model=row['model'];n=inventory['models'].get(model,{}).get('numeric_entry_count','in archives; not enumerated')
        folder=f'/data/pyoung3/groupdata/svierlin/{model}_data'
        if model!='50Am':folder+=f'/{model}_outputs/&lt;ID&gt;/'
        else:folder+='/eoc_batch_NNNN.tar.gz::&lt;ID&gt;/'
        modelrows+=f'<tr><td><b>{model}</b></td><td>{int(row["npart"]):,}</td><td>{n:,}</td><td><code>{folder}</code></td></tr>' if isinstance(n,int) else f'<tr><td><b>{model}</b></td><td>{int(row["npart"]):,}</td><td>{n}</td><td><code>{folder}</code></td></tr>'
    plots={
      '01_snapshot_times':'Snapshot suffix versus physical time. 50Am and jet3b time factors are inferred from shared code; cco2 and vconv factors come from their headers.',
      '02_particle_geometry':'Where the fluid particles are. Each panel samples 20,000 rows at its own last available snapshot; epochs and color scales differ.',
      '03_requested_composition_coverage':'Which requested composition cells are explicitly stored in the inspected examples. Sparse omissions stay blank, not zero.',
      '04_prism_sanity':'Mass-fraction normalization and PRISM endpoint times. Reaching the configured endpoint does not validate input physics.',
      '05_input_trajectories':'Inputs as written, with disputed time/density units. The blank 50Am panel reflects missing configured input paths.',
      '06_hydro_composition_flags':'Raw hydro composition checks. Stored zero Y_el and isotope-slot sums above one need interpretation.',
      '08_unit_crosscheck':'Existing inputs compared with identified raw particles. The 50Am ZIP variants are candidates, not proven historical inputs.'}
    gallery=''.join(f'<figure><a href="../plots/{name}.png"><img loading="lazy" src="../plots/{name}.png" alt="{html.escape(caption)}"></a><figcaption>{html.escape(caption)} <a href="../plots/{name}.pdf">PDF</a></figcaption></figure>' for name,caption in plots.items() if (OUT/'plots'/f'{name}.png').exists())
    data=json.dumps({'fields':fields,'values':values},ensure_ascii=True).replace('</','<\\/')
    document='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Supernova particle data map</title>
<style>
:root{color-scheme:light;font:16px/1.6 system-ui,sans-serif;color:#193047;background:#f4f6f8}body{margin:0}main{max-width:1200px;margin:auto;padding:30px}h1{font-size:2.4rem;line-height:1.15}h2{margin-top:2.3rem}h3{margin-bottom:.4rem}a{color:#12668b}p{max-width:90ch}.meta{color:#526479}.lead{font-size:1.16rem}.card,.issue{background:white;border:1px solid #d5dce4;border-radius:10px;padding:18px;margin:14px 0}.issue{border-left:5px solid #b37421}.flow{display:grid;grid-template-columns:repeat(3,1fr);gap:18px}.flow .card{border-top:5px solid #207d75}.flow b{font-size:1.12rem}.tablewrap{overflow:auto;border:1px solid #d5dce4;background:white}table{border-collapse:collapse;width:100%;font-size:.9rem}th,td{text-align:left;padding:9px 12px;border-bottom:1px solid #e3e7eb;vertical-align:top}th{background:#eaf0f4}code{font-size:.84em;overflow-wrap:anywhere}input,select{padding:10px;border:1px solid #92a4b5;border-radius:5px;font:inherit;max-width:90%}.controls{display:flex;flex-wrap:wrap;gap:12px;margin:14px 0}figure{background:white;padding:16px;margin:20px 0;border-radius:10px;border:1px solid #d5dce4}img{width:100%;height:auto}figcaption{color:#44596e}nav{display:flex;gap:18px;flex-wrap:wrap}.pill{display:inline-block;padding:5px 10px;background:#e3f0ed;border-radius:20px;margin:5px 6px 5px 0}button{border:0;background:transparent;color:#12668b;text-decoration:underline;cursor:pointer;font:inherit;text-align:left}#detail{overflow:auto}.small{font-size:.9rem}@media(max-width:750px){main{padding:18px}.flow{grid-template-columns:1fr}h1{font-size:1.9rem}}
</style></head><body><main>
<p class="meta">ASU SOL · Read-only source audit · 25 September 2026</p>
<h1>Find a particle. Follow its history.<br>Know what each number means.</h1>
<p class="lead">A guide to 50Am, jet3b, cco2 and the vconv/vconvL association, built from actual files, headers, configurations and source code.</p>
<nav><a href="#map">Data map</a><a href="#findings">Audit findings</a><a href="#lookup">298-field lookup</a><a href="#plots">Debug plots</a><a href="#files">Files & scope</a></nav>
<p><span class="pill">483 snapshot headers validated</span><span class="pill">4 selected snapshot ID sets checked</span><span class="pill">SAMPLE_COUNT PRISM particles inspected</span></p>
<h2 id="map">Three linked kinds of data</h2>
<div class="flow"><div class="card"><b>1. Hydro snapshots</b><p>One file = one simulation time.<br>One record = one moving fluid particle.</p><code>ident, x, y, z, vx, mass, rho, temp …</code><p>Read the SDF header for the layout and time. A filename suffix is an iteration.</p></div>
<div class="card"><b>2. Particle trajectories</b><p>One file = one particle’s temperature and density history.</p><code>time(s), temp(GK), density(g/cm3)</code><p>PRISM uses this history as input. Its unit labels need to agree with the original simulation.</p></div>
<div class="card"><b>3. PRISM composition</b><p>One numbered folder = one particle’s nuclear calculation.</p><code>x.txt: Z, A, X<br>xz.txt: Z, X(element)<br>xtime.txt: repeated time blocks</code><p>Join by model and particle ID. Retain the abundance epoch separately.</p></div></div>
<div class="tablewrap"><table><thead><tr><th>Observed model</th><th>Particles in selected hydro snapshot</th><th>Numeric output entries</th><th>Composition location</th></tr></thead><tbody>MODEL_ROWS</tbody></table></div>
<p class="small">Numeric names are candidate particle entries, not proof of completed calculations. Hydro counts and output counts describe different populations. The user’s “grain” means SPH particle; smoothing length h is not dust radius. The observed vconv products are kept labeled vconv; complete equivalence with vconvL is unresolved.</p>
<h2 id="findings">What the adversarial checks found</h2>
<div class="issue"><b>Time and density discrepancies exist in the inputs.</b><p>Sampled cco2/vconv input times miss approximately a ×100 conversion; jet3b uses iteration numbers as seconds. All three sampled inputs have density approximately 0.240196 times the header-supported/reference value. Correcting output labels cannot repair nucleosynthesis calculated with incorrect inputs.</p><a href="unit_crosscheck.md">Numerical crosscheck</a> · <a href="provenance.md">Source-code evidence</a></div>
<div class="issue"><b>50Am has ambiguous input history and suspect archives.</b><p>The configured example trajectory paths are missing. Two ZIP variants exist: the original labels Kelvin as GK; the temperature-corrected variant still uses iteration times and raw density under physical-unit labels. We have not established which variant a historical PRISM job consumed. Archive creation logs show “No space left on device” for batches 1, 8 and 9; current archives have not had full integrity checks.</p><a href="coverage.md">Coverage and archive evidence</a></div>
<div class="issue"><b>Missing isotope rows are not measured zeros.</b><p>PRISM writes sparse outputs with a default cutoff of 10⁻¹⁵. Use blank values plus status, distinguish outside-network species, and sum all relevant isotopes for an element. Some earlier CSVs inserted zeros or lost particle/time identifiers.</p><a href="prism-formats.md">Formats and old-CSV audit</a></div>
<div class="issue"><b>Field presence does not prove a valid quantity.</b><p>The sampled late jet3b/vconv records store zero Y_el; 50Am/jet3b hydro nuclear-slot sums reach approximately 1.03. The plots preserve these values for diagnosis. Earlier “50AM_test” reports described a bundled Noh code test, not these explosion models.</p><a href="analysis-guide.md">Interpretation and extraction contract</a></div>
<h2 id="lookup">Search all requested fields</h2>
<p>Search a symbol such as <b>26Al</b>, <b>Y_{e}</b>, or <b>Peak density</b>. Click a request to inspect sampled values and source paths. Y is yttrium; Y_{e} is electron fraction. Counts refer only to inspected examples.</p>
<div class="controls"><label>Model <select id="model"><option value="">All models</option><option>50Am</option><option>jet3b</option><option>cco2</option><option>vconv</option></select></label><label>Field <input id="query" placeholder="e.g. 26Al" type="search"></label></div><p id="count" class="meta"></p>
<div class="tablewrap"><table><thead><tr><th>Model</th><th>Request</th><th>Kind</th><th>Source / field</th><th>Stored / sampled</th><th>Sparse omission</th><th>Outside network</th></tr></thead><tbody id="fieldrows"></tbody></table></div><div id="detail" class="card" hidden></div>
<h2 id="plots">Debug plots</h2><p>These plots explain structure and expose checks. They are not scientific model comparisons at a common age. Select a plot for its full resolution.</p>GALLERY
<h2 id="files">Files, reproducibility and limits</h2>
<p><b>Home Git repository:</b> <code>/home/rpraka15/supernova-grain-analysis</code><br><b>Analysis outputs:</b> <code>/scratch/rpraka15/supernova-grain-analysis</code></p>
<p>The audit validated all 483 named snapshot headers and their expected file lengths; checked all particle IDs in four selected snapshots; sampled 20,000 hydro records per model; listed direct numeric PRISM entries; and inspected bounded particle/archive examples. It did not decompress all 2.17 TB of 50Am archives or certify every particle calculation.</p>
<ul><li><a href="../tables/requested_field_map.csv">Complete requested-field map (CSV)</a></li><li><a href="../tables/requested_values_sample.csv">Sampled requested values, statuses and provenance (CSV)</a></li><li><a href="../tables/snapshot_catalog.csv">483-snapshot catalog (CSV)</a></li><li><a href="../tables/particle_id_join_demo.csv">ID join demonstration: differing epochs retained (CSV)</a></li><li><a href="../tables/hydro_diagnostics.csv">Hydro diagnostics (CSV)</a></li><li><a href="../tables/prism_diagnostics.csv">PRISM diagnostics (CSV)</a></li><li><a href="../inventory/inventory.json">Inventory, archive metadata and log findings (JSON)</a></li><li><a href="../debug/run_summary.json">Run scope, versions and error counts (JSON)</a></li><li><a href="../debug/prism_errors.json">Missing input paths / read errors (JSON)</a></li></ul>
<p>Code, reports and tests are tracked locally in Git. GitHub publication awaits authentication and a destination repository. No collaborator files were modified or deleted. A frozen copy of the analysis scripts lives in <code>../code/</code>.</p>
<p class="meta small">Generated from audit tables. RUN_DATE</p></main>
<script id="auditdata" type="application/json">DATA_JSON</script><script>
const data=JSON.parse(document.getElementById('auditdata').textContent);
const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function render(){const model=document.getElementById('model').value,q=document.getElementById('query').value.trim().toLowerCase();const rows=data.fields.filter(r=>(!model||r.model===model)&&(!q||r.request.toLowerCase().includes(q)));document.getElementById('count').textContent=`${rows.length} matching model–field pairs; showing ${Math.min(rows.length,100)}. Narrow the search to see more.`;document.getElementById('fieldrows').innerHTML=rows.slice(0,100).map(r=>`<tr><td>${esc(r.model)}</td><td><button data-model="${esc(r.model)}" data-request="${esc(r.request)}">${esc(r.request)}</button></td><td>${esc(r.kind)}</td><td>${esc(r.source_schema)}<br><code>${esc(r.hydro_field_or_derivation)}</code></td><td>${esc(r.stored_samples)} / ${esc(r.sample_particles)}</td><td>${esc(r.sparse_missing_samples)}</td><td>${esc(r.outside_network_samples)}</td></tr>`).join('');}
document.getElementById('fieldrows').addEventListener('click',e=>{const b=e.target.closest('button');if(!b)return;const rows=data.values.filter(r=>r.model===b.dataset.model&&r.request===b.dataset.request);const d=document.getElementById('detail');d.hidden=false;d.innerHTML=`<h3>${esc(b.dataset.model)} · ${esc(b.dataset.request)}</h3><p>Stored composition is final PRISM X; physical snapshot values are intentionally separate until an epoch is chosen. A blank is not zero.</p><table><thead><tr><th>Particle</th><th>Value</th><th>Status</th><th>PRISM stop as written [s]</th><th>Source</th></tr></thead><tbody>${rows.map(r=>`<tr><td>${esc(r.particle_id)}</td><td>${r.value===''?'—':esc(r.value)}</td><td>${esc(r.status)}</td><td>${esc(r.PRISM_stop_as_written_s)}</td><td><code>${esc(r.source)}</code></td></tr>`).join('')}</tbody></table>`;});
document.getElementById('model').addEventListener('change',render);document.getElementById('query').addEventListener('input',render);render();
</script></body></html>'''
    document=document.replace('SAMPLE_COUNT',str(sum(summary['PRISM_particles_sampled'].values()))).replace('MODEL_ROWS',modelrows).replace('GALLERY',gallery).replace('RUN_DATE',html.escape(summary['completed_utc'])).replace('DATA_JSON',data)
    (OUT/'reports/index.html').write_text(document)
    for p in (REPO/'docs').glob('*.md'):shutil.copy2(p,OUT/'reports'/p.name)
    for p in (REPO/'scripts').glob('*.py'):shutil.copy2(p,OUT/'code'/p.name)
    for p in ['README.md','requirements.txt']:shutil.copy2(REPO/p,OUT/'code'/p)
    manifest={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
              for directory in ['scripts','docs','tests'] for p in (REPO/directory).glob('*') if p.is_file()}
    manifest['requested_fields_input']=hashlib.sha256(Path('/scratch/rpraka15/List of isotopes_elements_variables.txt').read_bytes()).hexdigest()
    (OUT/'debug/code_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(OUT/'reports/index.html')

if __name__=='__main__':main()
