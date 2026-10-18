#!/usr/bin/env python3
"""Compare identified SPH particles with existing PRISM input trajectories.

Uses the main audit's provenance tables, read-only SDF memmaps, and bounded
particle/epoch samples. No collaborator scripts execute and no sources change.
"""
from pathlib import Path
import csv
import json
import os
import io
import zipfile

OUT = Path('/scratch/rpraka15/supernova-grain-analysis')
os.environ.setdefault('MPLCONFIGDIR', str(OUT / 'cache/matplotlib'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from readers import open_sdf, units

PARTICLES = {'50Am': 191136, 'jet3b': 10076, 'cco2': 118326, 'vconv': 251189}


def main():
    if not OUT.resolve().is_relative_to(Path('/scratch/rpraka15').resolve()):
        raise ValueError('Output escaped user scratch')
    with (OUT / 'tables/prism_diagnostics.csv').open() as f:
        diagnostics = list(csv.DictReader(f))
    with (OUT / 'tables/snapshot_catalog.csv').open() as f:
        snapshots = list(csv.DictReader(f))
    results, notes = [], []
    cases = [(m, p, 'existing_file') for m, p in PARTICLES.items() if m != '50Am']
    cases += [('50Am', PARTICLES['50Am'], v) for v in ['archived_original', 'archived_tempcorr']]
    for model, pid, variant in cases:
        config = next((r for r in diagnostics if r['model'] == model
                       and int(r['particle_id']) == pid), None)
        if config is None:
            notes.append({'model': model, 'particle_id': pid, 'status': 'no_config_in_audit_sample'})
            continue
        path = Path(config['trajectory_source'])
        trajectory_source = str(path)
        if model == '50Am' and not path.is_file():
            stem = '50Am_output' if variant == 'archived_original' else '50Am_output_tempcorr'
            archive = path.parent.parent / (stem + '.zip')
            member = f'{stem}/{pid}_trajectory.txt'
            with zipfile.ZipFile(archive) as z:
                info = z.getinfo(member)
                if info.file_size > 1024**2:
                    raise ValueError('Trajectory member exceeds 1 MiB bound')
                trajectory = np.loadtxt(io.StringIO(z.read(info).decode('ascii')), skiprows=1, ndmin=2)
            trajectory_source = str(archive) + '::' + member
            notes.append({'model': model, 'particle_id': pid, 'variant': variant,
                          'status': 'configured_path_missing_read_archived_candidate',
                          'trajectory_path': str(path), 'archive_member': trajectory_source,
                          'caveat': 'No proof that this archived variant was consumed by the sampled PRISM run'})
        elif not path.is_file():
            notes.append({'model': model, 'particle_id': pid, 'status': 'configured_trajectory_missing',
                          'trajectory_path': str(path), 'config_source': config['source']})
            continue
        else:
            trajectory = np.loadtxt(path, skiprows=1, ndmin=2)
        if trajectory.shape[1] != 3 or not np.isfinite(trajectory).all():
            raise ValueError(f'Invalid trajectory: {path}')
        candidates = sorted((r for r in snapshots if r['model'] == model),
                            key=lambda r: int(r['iteration']))
        # vconv trajectory generation ends at 3380; the separate 4200 binary
        # belongs to a later epoch and must not be paired with the final row.
        if model == 'vconv':
            candidates = [r for r in candidates if int(r['iteration']) <= 3380]
        selections = [candidates[i] for i in sorted({0, 1, len(candidates) - 1})]
        for snapshot in selections:
            h, records = open_sdf(snapshot['path'])
            factor = units(h['meta'])
            indices = np.flatnonzero(records['ident'] == pid)
            if len(indices) != 1:
                notes.append({'model': model, 'particle_id': pid, 'snapshot': snapshot['path'],
                              'status': 'particle_missing_or_duplicated', 'matches': len(indices)})
                del records
                continue
            particle = records[indices[0]].copy()
            del records
            temp_gk = float(particle['temp']) / 1e9
            rho_code = float(particle['rho'])
            time_code = h['meta']['tpos']
            iteration = int(h['meta']['iter'])
            # Independently compare temperature, allowing only text rounding.
            expected_written_temp = temp_gk * 1e9 if variant == 'archived_original' else temp_gk
            temperature_matches = np.flatnonzero(np.isclose(trajectory[:, 1], expected_written_temp,
                                                            rtol=8e-6, atol=0))
            expected_written_time = iteration if model in ['jet3b', '50Am'] else time_code
            convention_matches = np.flatnonzero(np.isclose(trajectory[:, 0], expected_written_time,
                                                           rtol=.006, atol=1e-16))
            matches = np.intersect1d(temperature_matches, convention_matches)
            if len(matches) != 1:
                notes.append({'model': model, 'particle_id': pid, 'snapshot': snapshot['path'],
                              'status': 'no_unique_temperature_and_time_convention_match',
                              'temperature_matches': temperature_matches.tolist(),
                              'convention_matches': convention_matches.tolist()})
                continue
            i = int(matches[0])
            written_time, written_temp, written_density = map(float, trajectory[i])
            physical_time = time_code * factor['time_s']
            physical_density = rho_code * factor['density_g_cm3']
            results.append({
                'model': model, 'variant': variant, 'particle_id': pid, 'iteration': iteration,
                'snapshot_path': snapshot['path'], 'trajectory_path': trajectory_source,
                'config_source': config['source'], 'trajectory_data_row_1based': i + 1,
                'epoch_match': 'particle_ID_plus_temperature_plus_script_time_convention',
                'temperature_matching_rows': len(temperature_matches),
                'raw_temperature_GK': temp_gk, 'written_temperature_GK': written_temp,
                'temperature_ratio_written_over_raw': written_temp / temp_gk,
                'time_code': time_code, 'time_s_reference': physical_time,
                'time_s_written': written_time,
                'time_ratio_written_over_reference': written_time / physical_time if physical_time else '',
                'rho_code': rho_code, 'rho_g_cm3_reference': physical_density,
                'rho_g_cm3_written': written_density,
                'density_ratio_written_over_reference': written_density / physical_density,
                'observed_density_conversion': written_density / rho_code,
                'reference_density_conversion': factor['density_g_cm3'],
                'reference_time_conversion': factor['time_s'],
                'reference_units_provenance': factor['provenance'],
            })
            print(model, pid, iteration, 't_written=', written_time, 't_reference=', physical_time,
                  'rho_ratio=', written_density / physical_density, flush=True)
    if results:
        with (OUT / 'tables/unit_crosscheck.csv').open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(results[0]))
            writer.writeheader()
            writer.writerows(results)
    (OUT / 'debug/unit_crosscheck_notes.json').write_text(json.dumps(notes, indent=2) + '\n')
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout='constrained')
    colors = {'50Am': '#805699', 'jet3b': '#bb4430', 'cco2': '#176d93', 'vconv': '#517b39'}
    for model in colors:
        rows = [r for r in results if r['model'] == model]
        times = [r for r in rows if r['time_s_reference'] > 0 and r['time_s_written'] > 0]
        axes[0].scatter([r['time_s_reference'] for r in times],
                        [r['time_s_written'] for r in times], label=model, c=colors[model], s=45)
        axes[1].scatter([model] * len(rows),
                        [r['density_ratio_written_over_reference'] for r in rows],
                        c=colors[model], s=45)
    axes[0].plot([1e-4, 1e9], [1e-4, 1e9], '--', color='0.5', label='Correct seconds')
    axes[0].set(xscale='log', yscale='log', xlabel='Reference physical time [s]',
                ylabel='Existing trajectory column labeled time(s)', title='Saved trajectory times disagree')
    axes[0].legend(fontsize=9)
    axes[1].axhline(1, ls='--', color='0.5')
    axes[1].axhline(1.42e-6 / (1.9889e27 / 6.955e10**3), ls=':', color='0.3')
    axes[1].set(ylabel='Written density / reference CGS density', yscale='log',
                title='Density labels disagree with reference units')
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('Empirical check: one identified particle per model, three sampled snapshots\n'
                 '50Am/jet3b factors inferred; cco2/vconv factors explicit in SDF headers', fontsize=11)
    for suffix in ['png', 'pdf']:
        fig.savefig(OUT / 'plots' / f'08_unit_crosscheck.{suffix}', dpi=180, bbox_inches='tight')
    plt.close(fig)
    report = ['# Existing trajectory unit crosscheck', '',
              'This check reads existing trajectory files referenced by audited particle configurations. '
              'Particle IDs are selected explicitly from the raw SDF records. Epoch matches require '
              'temperature agreement within text-rounding tolerance and agreement with the inspected '
              'generation script time convention. Rows are never joined by position. Some early '
              'vconv temperatures repeat, so temperature alone would be insufficient.', '',
              '| Model / variant | Particle | Iteration | Written time (s label) | Reference seconds | Density ratio |',
              '|---|---:|---:|---:|---:|---:|']
    for r in results:
        report.append(f"| {r['model']} / {r['variant']} | {r['particle_id']} | {r['iteration']} | "
                      f"{r['time_s_written']:.7g} | {r['time_s_reference']:.7g} | "
                      f"{r['density_ratio_written_over_reference']:.7g} |")
    report += ['', 'The sampled existing files corroborate the script-level discrepancies: jet3b '
               'writes iteration as time; cco2 and vconv write rounded code time, missing the factor '
               '100; all three use approximately 1.42e-6 for density instead of '
               '5.911824660348532e-6. The density ratio is approximately 0.2401965 '
               '(reference is 4.16326 times larger). These are observed output discrepancies, '
               'not only potential script hazards. They do not by themselves establish the '
               'downstream effect on nuclear abundances.', '',
               '50Am particle 191136 has a missing configured trajectory path; two archived candidate '
               'variants were read directly from ZIP members without extraction. The original variant '
               'writes Kelvin values under a GK label (factor 1e9); the tempcorr variant fixes that '
               'temperature conversion. Both store iteration as seconds and raw code density under '
               'the CGS label. Archived variants have not been established as the actual input consumed '
               'by the sampled PRISM run. '
               'The jet3b reference factors are inferred from shared source code because the sampled '
               'headers omit conversion factors; the same limitation applies to 50Am. '
               'cco2/vconv factors are explicit in their headers. '
               'The vconv samples here do not resolve vconvL model identity.', '',
               'Detailed numerical provenance: `tables/unit_crosscheck.csv`. Missing or ambiguous '
               'matches: `debug/unit_crosscheck_notes.json`.']
    (OUT / 'reports/unit_crosscheck.md').write_text('\n'.join(report) + '\n')


if __name__ == '__main__':
    main()
