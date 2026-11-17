"""Collect existing verification artifacts; never infer success from output presence."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / 'decks/subcmos2014/results'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    names = [
        'astra_final_ild', 'astra_final_ild_sanitized', 'astra_final_spacer_sanitized', 'astra_final_pmos_finish/pmos',
        'astra_regularized_spacer', 'astra_regularized_spacer_sanitized',
        'astra_domain_ild_float_grid', 'astra_domain_ild_sanitized',
        'astra_domain_spacer_float_grid', 'astra_domain_spacer_sanitized',
        'astra_domain_pmos_finish/pmos',
        'astra_domain_fresh/nmos', 'astra_domain_fresh/pmos',
        'astra_domain_stable_angles/pmos', 'astra_ordered_projection_nmos/nmos',
        'astra_be_activation/pmos', 'astra_be_error_diagnosis/pmos',
        'astra_p1_fresh/pmos', 'astra_p1_fresh/nmos',
        'astra_velocity_history_continue/pmos',
        'astra_finite_clip_continue/pmos', 'astra_activation_geometry/pmos',
        'astra_finite_spacer_sanitized', 'astra_reflecting_control',
        'astra_reflecting_sanitized', 'astra_tangent_boundary_continue/pmos',
        'astra_tangent_from_wells/pmos', 'astra_tangent_nmos/nmos',
        'sol_side_boundary_pad_replay', 'sol_side_boundary_pad_sanitized',
        'sol_side_boundary_full_activation/pmos',
        'sol_side_boundary_fresh/nmos', 'sol_side_boundary_fresh/pmos',
        'sol_side_boundary_threshold_control/pmos',
        'sol_kooi_back_pedal_control', 'sol_kooi_back_pedal_fixed',
        'sol_kooi_back_pedal_sanitized',
        'sol_back_pedal_double_continue/pmos',
        'sol_back_pedal_spacer_diag/pmos',
        'sol_pmos_spacer_clip_control',
        'sol_stable_weights_activation', 'sol_stable_weights_activation_sanitized',
        'sol_stable_weights_spacer_nmos', 'sol_stable_weights_spacer_pmos',
        'sol_stable_weights_locos', 'sol_stable_weights_full_activation/pmos',
        'sol_stable_weights_full_sanitized/pmos', 'sol_stable_weights_pmos_finish/pmos',
        'sol_stable_weights_fresh/nmos', 'sol_stable_weights_fresh/pmos',
        'sol_domain_final_full_activation/pmos', 'sol_domain_valid_estimator_full/pmos',
        'sol_domain_residual_trace/pmos', 'sol_ild_boundary_fixture_control',
        'sol_pad_motion_capture',
    ]
    records = []
    for name in names:
        p = RESULTS / name / 'result.json'
        record = {'case': name, 'result': str(p.relative_to(ROOT))}
        if p.exists():
            data = json.loads(p.read_text())
            record.update({k: data[k] for k in (
                'passed', 'returncode', 'timed_out', 'binary', 'binary_sha256',
                'elapsed_seconds', 'provenance', 'errors') if k in data})
            record['result_sha256'] = digest(p)
        else:
            record['status'] = 'INCOMPLETE'
        records.append(record)
    builds = []
    latest = ('native64_l3ozbabw', 'native64_92cvoyp7')
    for name in ('native64_m23qw5jt', 'native64_zhq5f5_p',
                 'native64_n0paar46', 'native64_k8ofsow3',
                 'native64_6jztdvdl', 'native64_32q9_n2u',
                 'native64_yuj0z1w7', 'native64_70li8k7h',
                 'native64_o54fuavs', 'native64_99j1uzl7', *latest):
        p = ROOT / 'suprem4/build' / name
        manifest = json.loads((p / 'build.json').read_text())
        mismatches = []
        for filename, expected in manifest['source_sha256'].items():
            # The builder intentionally generates Makefile from Makefile.proto.
            actual = ROOT / 'suprem4/src' / (
                'Makefile.proto' if filename == 'Makefile' else filename)
            if digest(actual) != expected:
                mismatches.append(filename)
        builds.append({'build': name, 'latest_candidate': name in latest,
                       'binary_sha256': digest(p / 'src/suprem4'),
                       'build_manifest_sha256': digest(p / 'build.json'),
                       'current_source_mismatches': mismatches})
    test_files = ['geometry_weights.c', 'run_geometry_weights.sh', 'time_stage.c',
                  'test_time_stage.py', 'material_domain.c', 'material_domain_db.c',
                  'run_material_domain.sh', 'moving_flux.c', 'run_moving_flux.sh',
                  'mesh_trial.c', 'run_mesh_trial.sh', 'diffusion_retry.c',
                  'test_diffusion_retry.py', 'mesh_patch_commit.c',
                  'test_mesh_patch_c.py', 'run_mesh_patch_commit.sh',
                  'run_rit_locos_strip.py', 'reflecting_normal.c',
                  'run_reflecting_normal.sh', 'etch_corner.c', 'run_etch_corner.sh',
                  'interface_measure.c', 'run_interface_measure.sh',
                  'be_time.c', 'run_be_time.sh', 'be_mass_term.c', 'run_be_mass_term.sh',
                  'p1_mass.c', 'run_p1_mass.sh', 'boundary_angle.c', 'run_boundary_angle.sh',
                  'rect_projection.c', 'run_rect_projection.sh',
                  'domain_offset.cpp', 'run_domain_offset.sh']
    report = {
        'date': '2026-11-15',
        'scope': 'Candidate numerical/domain repairs; full CMOS and physics acceptance remain open',
        'candidate_options': {'SUPREM4_MATERIAL_DOMAINS': 'oxide',
                              'SUPREM4_MESH_PATCH_REMOVAL': '1',
                              'SUPREM4_DEPOSITION_GEOMETRY': 'domain'},
        'experimental_build_flags': ['--ordering', 'amd', '--domain-offset'],
        'sanitizer_options': {'ASAN_OPTIONS': 'detect_leaks=0',
                              'UBSAN_OPTIONS': 'halt_on_error=1'},
        'validated_launcher_target': str((ROOT / 'suprem4/build/validated-native64').readlink()),
        'builds': builds, 'runs': records,
        'test_source_sha256': {name: digest(ROOT / 'suprem4/tests' / name) for name in test_files},
        'focused_check_report': {
            'path': 'suprem4/experiments/reflecting_boundary_20261115/focused_checks.json',
            'sha256': digest(ROOT / 'suprem4/experiments/reflecting_boundary_20261115/focused_checks.json'),
        },
        'domain_backend_sanitizer_log': {
            'path': 'suprem4/experiments/robust_offset_20261115/domain_offset_final.log',
            'sha256': digest(ROOT / 'suprem4/experiments/robust_offset_20261115/domain_offset_final.log'),
        },
        'negative_controls': {
            'geometry_weights': 'geometry_weights_negative_control.log',
            'time_stage': 'time_stage_negative_control.log',
            'directory': str(Path(__file__).parent.relative_to(ROOT)),
        },
        'independent_audits': [
            'decks/subcmos2014/results/astra_final_pmos_finish/pmos/contacts_partition_audit.json',
            'decks/subcmos2014/results/astra_final_pmos_finish/pmos/metal1_partition_audit.json',
            'decks/subcmos2014/results/astra_domain_pmos_finish/pmos/contacts_partition_audit.json',
            'decks/subcmos2014/results/astra_domain_pmos_finish/pmos/metal1_partition_audit.json',
            'decks/subcmos2014/results/sol_side_boundary_fresh/nmos/annealed_partition_audit.json',
            'decks/subcmos2014/results/sol_side_boundary_full_activation/pmos/annealed_partition_audit.json',
            'decks/subcmos2014/results/sol_stable_weights_full_activation/pmos/annealed_partition_audit.json',
            'decks/subcmos2014/results/sol_stable_weights_full_sanitized/pmos/annealed_partition_audit.json',
            'decks/subcmos2014/results/sol_stable_weights_ild_boundary/boundary_audit.json',
        ],
    }
    target = ROOT / 'suprem4/doc/material_domain_reliability_results.json'
    target.write_text(json.dumps(report, indent=2) + '\n')
    print(target)


if __name__ == '__main__':
    main()
