"""Frozen held-out comparison. Missing, adverse and nonfinite evidence stays visible."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
from summarize_public_pilot import calibration, finite, interval, paired_summary, safe_json

ROOT = Path(__file__).parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(a, b):
    return np.allclose(a, b, atol=1e-8, rtol=1e-8, equal_nan=False)


def validate_record(r, p, digest, arm):
    expected_algorithm = p['algorithm_commit'] if arm == 'baseline' else r['executionCommit']
    require(r['algorithmCommit'] == expected_algorithm, 'Algorithm revision mismatch')
    require(r['measureSha256'] == p[f'{arm}_measure_sha256'], 'Measure digest mismatch')
    require(r['protocolSha256'] == digest, 'Protocol mismatch')
    require(r['searchRngSeed'] == p['search_rng_seed'], 'Search RNG mismatch')
    require(r['executionCommit'] and r['workflowRun'] and r['runAttempt'], 'Missing execution identity')
    require('2026a' in r['matlab'] and r['eeglab'] == p['toolchain']['eeglab'] and r['eeglabSha'], 'Toolchain mismatch')
    entry = next(x for x in p['source']['files'] if x['subject'] == r['subject'])
    source = r['source']; excerpt = p['excerpt']
    require(source['subject'] == r['subject'] and source['source_sha256'] == entry['sha256'], 'Source identity mismatch')
    require(source['source_url'] == p['source']['base_url'] + entry['path'], 'Source URL mismatch')
    require(source['protocol_sha256'] == digest and source['channels'] == excerpt['channels'], 'Source protocol/channels mismatch')
    require(source['sampling_rate'] == excerpt['sampling_rate'] and source['samples'] == excerpt['duration_seconds'] * excerpt['sampling_rate'], 'Excerpt samples mismatch')
    require(source['offset_seconds'] == excerpt['offset_seconds'] and source['duration_seconds'] == excerpt['duration_seconds'], 'Excerpt duration mismatch')
    require(r['completed'] and not r['errors'], 'Incomplete execution: ' + repr(r.get('errors')))
    require(r.get('notRun') == 0, 'Incomplete Standard search')
    opts = r['searchOptions']
    require(opts['icaSeeds'] == 2 and opts['refineCutoffs'] and opts['eogCheck'] and opts['signalTemplate'] == 'data', 'Standard options changed')
    for name in ['erp', 'minimal', 'common']:
        x = r[name]
        require(x['available'], name + ': no estimate/recommendation')
        estimate = np.asarray(x['estimate'], float)
        pop = np.asarray(r['populationTruth'], float); real = np.asarray(r['realizedTruth'], float)
        require(estimate.shape == pop.shape == real.shape == (2,), 'Unexpected condition shape')
        require(finite(x['gain']) and x['gain'] > 0 and np.isfinite(estimate).all(), 'Invalid gain/estimate')
        corrected = estimate / x['gain']
        require(close(corrected, x['correctedEstimate']), 'Corrected estimate mismatch')
        for key, expected in [('actualError', corrected-pop), ('rawError', estimate-pop), ('conditionalError', corrected-real)]:
            require(close(expected, x[key]), 'Error field mismatch: ' + key)
        for key, field in [('mse', 'actualError'), ('rawMse', 'rawError'), ('conditionalMse', 'conditionalError')]:
            require(finite(x[key]) and close(float(np.mean(np.asarray(x[field], float)**2)), x[key]), 'MSE mismatch: ' + key)
        reported = np.asarray(x['reportedError'], float)
        require(reported.shape == (2,), 'Unexpected reported error shape')
        if np.isfinite(reported).all():
            require((reported >= 0).all() and finite(x['reportedMse']) and close(float(np.mean(reported**2)), x['reportedMse']), 'Reported MSE mismatch')
        else:
            require(not finite(x['reportedMse']), 'Nonfinite reported errors hidden by finite MSE')
    for key in ['replayMaxError', 'exportMaxError']:
        require(finite(r[key]) and 0 <= r[key] < p['evaluation']['replay_tolerance_uv'], key + ' failed')
    support = r['selectedSupport']
    require(len(support) == 2, 'Missing selected condition support')
    for c in support:
        ids = c['referenceIds']; scores = np.asarray(c['scores'], float)
        weights = np.asarray(c['weights'], float); flags = np.asarray(c['flagged'], bool)
        require(len(set(ids)) == len(ids) == c['referenceN'], 'Invalid trial identities')
        require(scores.shape == weights.shape == flags.shape == (len(ids),), 'Trial audit alignment mismatch')
        require(np.isfinite(weights).all() and (weights >= 0).all(), 'Invalid diagnostic weights')
        good = np.isfinite(scores) & (weights > 0)
        require((weights[~np.isfinite(scores)] == 0).all(), 'Absent trial has positive weight')
        require(c['kept'] == int(np.isfinite(scores).sum()) and close(c['retention'], c['kept']/c['referenceN']), 'Retention mismatch')
        require(c['positiveN'] == int(good.sum()) and c['flaggedN'] == int((good & flags).sum()) and c['unflaggedN'] == int((good & ~flags).sum()), 'Trial counts mismatch')
        for key, w in [('effectiveN', weights), ('effectiveFlaggedN', weights[flags]), ('effectiveUnflaggedN', weights[~flags])]:
            want = float(w.sum()**2/max(float((w*w).sum()), np.finfo(float).eps))
            require(close(c[key], want), 'Effective support mismatch')
    boot = r['bootstrapAudit']
    # MATLAB encodes one-element numeric vectors as JSON scalars.
    ids = np.atleast_1d(boot['feasibleIds']).tolist()
    counts = np.atleast_1d(boot['nonfiniteByCandidate']).tolist()
    require(boot['draws'] == 1999 and len(ids) == r['feasible'], 'Bootstrap audit coverage mismatch')
    require(len(counts) == r['feasible'], 'Bootstrap audit length mismatch')
    require(all(isinstance(x, int) and 0 <= x <= boot['draws'] for x in counts), 'Invalid nonfinite counts')
    require(0 <= boot['nonfiniteColumns'] <= boot['draws'], 'Invalid nonfinite column count')


def uncertainty_screen(r, p):
    g = p['decision_gates']; reasons = []
    for c in r['selectedSupport']:
        minimum = g['minimum_effective_unflagged_when_flagged'] if c['flaggedN'] > 0 else g['minimum_effective_unflagged_without_flagged']
        if not c['eyeBiasAssessed'] or not c['eyeBiasFinite']:
            reasons.append(c['condition'] + ': bias unassessed/nonfinite')
        if c['effectiveUnflaggedN'] < minimum - 1e-10:
            reasons.append(c['condition'] + ': insufficient effective unflagged support')
    if not finite(r['erp']['reportedMse']) or not all(finite(x) for x in r['erp']['reportedError']):
        reasons.append('Selected reported error nonfinite')
    if r['bootstrapAudit']['nonfiniteColumns'] > 0:
        reasons.append('Feasible-family bootstrap contains nonfinite scores')
    return dict(passed=not reasons, reasons=reasons)


def decide(rows, p):
    g = p['decision_gates']; n = p['evaluation']['expected_records']
    if len(rows) != n:
        return dict(decision='incomplete', pairs=len(rows), expected=n, gates={}, meanComparison=None)
    b = np.array([x['baselineMSE'] for x in rows]); v = np.array([x['variantMSE'] for x in rows])
    br = np.array([x['baselineRawMSE'] for x in rows]); vr = np.array([x['variantRawMSE'] for x in rows])
    d = v-b; dr = vr-br
    idx = np.random.RandomState(p['evaluation']['bootstrap_seed']).randint(0, n, (p['evaluation']['bootstrap_draws'], n))
    ci = interval(d[idx].mean(1)); tol = g['tie_tolerance_uv2']
    wins = int((d < -tol).sum()); losses = int((d > tol).sum())
    gates = dict(meanReductionAtLeastFivePercent=bool(b.mean() > 0 and v.mean() <= b.mean()*(1-g['minimum_relative_mean_mse_reduction'])), pairedCIEntirelyBelowZero=ci[1] < 0,
                 rawMeanNotWorse=bool(vr.mean() <= br.mean()+tol), winsExceedLosses=wins > losses,
                 noMaterialIndividualRegression=all(x['individualGuardrailPassed'] for x in rows),
                 variantUncertaintySupported=all(x['variantUncertainty']['passed'] for x in rows))
    worse = v.mean() > b.mean()+tol or vr.mean() > br.mean()+tol
    decision = 'reject-variant-retain-original' if worse else ('retain-heldout-supported' if all(gates.values()) else 'inconclusive-no-promotion')
    means = dict(baselineMSE=float(b.mean()), variantMSE=float(v.mean()), difference=float(d.mean()), ci95=ci,
                 relativeReduction=None if b.mean() == 0 else float(1-v.mean()/b.mean()), baselineRawMSE=float(br.mean()), variantRawMSE=float(vr.mean()), rawDifferenceCI95=interval(dr[idx].mean(1)),
                 wins=wins, losses=losses, ties=n-wins-losses)
    return dict(decision=decision, pairs=n, expected=n, gates=gates, meanComparison=means)


def compare(folder):
    p = json.loads((ROOT/'heldout_protocol.json').read_text())
    digest = hashlib.sha256((ROOT/'heldout_protocol.json').read_bytes()).hexdigest()
    require(hashlib.sha256((ROOT/'within_variant/Measure.m').read_bytes()).hexdigest() == p['variant_measure_sha256'], 'Frozen variant changed')
    coverage = []; rows = []; records = []; identities = set()
    for subject in [f['subject'] for f in p['source']['files']]:
        paths = list(folder.glob(f'heldout-pair-{subject}-attempt-*'))
        paths.sort(key=lambda x: int(x.name.rsplit('-', 1)[1]))
        cov = dict(subject=subject, present=bool(paths), validPair=False, attempts=[x.name for x in paths], error=None)
        try:
            require(bool(paths), 'Missing pair artifact')
            # Highest attempt is authoritative, even when incomplete. Never select the best-looking attempt.
            pair = paths[-1]; attempt = pair.name.rsplit('-', 1)[1]
            b = json.loads((pair/'baseline/result.json').read_text()); v = json.loads((pair/'variant/result.json').read_text())
            records.append(dict(subject=subject, baseline=b, variant=v, artifact=pair.name))
            require(b['subject'] == v['subject'] == subject, 'Subject mismatch')
            for arm, r in [('baseline', b), ('variant', v)]:
                validate_record(r, p, digest, arm)
                require(r['runAttempt'] == attempt, 'Attempt mismatch')
                if os.getenv('GITHUB_SHA'):
                    require(r['executionCommit'] == os.environ['GITHUB_SHA'], 'Collector revision mismatch')
                if os.getenv('GITHUB_RUN_ID'):
                    require(r['workflowRun'] == os.environ['GITHUB_RUN_ID'], 'Collector run mismatch')
            for key in ['executionCommit', 'workflowRun', 'runAttempt', 'source', 'searchOptions', 'populationTruth', 'realizedTruth', 'eeglabSha', 'searchRngSeed']:
                require(b[key] == v[key], 'Pair identity mismatch: ' + key)
            for arm in ['minimal', 'common']:
                require(close(b[arm]['estimate'], v[arm]['estimate']) and close(b[arm]['gain'], v[arm]['gain']), 'Fixed control changed: ' + arm)
            for bc, vc in zip(b['selectedSupport'], v['selectedSupport']):
                require(bc['referenceIds'] == vc['referenceIds'] and bc['flagged'] == vc['flagged'], 'Reference labels/IDs changed')
            identities.add((b['executionCommit'], b['workflowRun'], b['eeglabSha']))
            delta = v['erp']['mse']-b['erp']['mse']
            allowance = max(p['decision_gates']['individual_regression_absolute_allowance_uv2'], p['decision_gates']['individual_regression_relative_allowance']*b['erp']['mse'])
            rows.append(dict(subject=subject, baselineMSE=b['erp']['mse'], variantMSE=v['erp']['mse'], difference=delta,
                             baselineRawMSE=b['erp']['rawMse'], variantRawMSE=v['erp']['rawMse'], individualRegressionAllowance=allowance,
                             individualGuardrailPassed=delta <= allowance+p['decision_gates']['tie_tolerance_uv2'],
                             baselinePipeline=b['erp']['pipeline'], variantPipeline=v['erp']['pipeline'],
                             baselineSupport=b['selectedSupport'], variantSupport=v['selectedSupport'],
                             baselineUncertainty=uncertainty_screen(b, p), variantUncertainty=uncertainty_screen(v, p),
                             baselineICAStable=b['icaSeedCheck']['same'], variantICAStable=v['icaSeedCheck']['same']))
            cov['validPair'] = True
        except (ValueError, KeyError, TypeError, OSError, StopIteration) as ex:
            cov['error'] = str(ex)
        coverage.append(cov)
    result = decide(rows, p)
    if len(identities) > 1:
        result.update(decision='incomplete', identityError='Mixed revisions, workflow runs or EEGLAB revisions')
    valid_subjects = {r['subject'] for r in rows}
    valid_records = [r for r in records if r['subject'] in valid_subjects]
    result.update(scope=p['scope'], protocolSha256=digest, productionChanged=False, mergePerformed=False,
                  coverage=coverage, rows=rows, executionIdentities=sorted(identities), uncertainty={})
    for arm in ['baseline', 'variant']:
        rr = [r[arm] for r in valid_records]
        c = calibration(rr, p)
        c['missingOrNonfiniteSubjects'] = [s['subject'] for s in p['source']['files'] if not any(x['subject'] == s['subject'] and finite(x['erp']['reportedMse']) for x in rr)]
        c['completeCohort'] = c['assessable'] == p['evaluation']['expected_records']
        result['uncertainty'][arm] = c
    result['fixedComparators'] = {arm: {control: paired_summary([r[arm] for r in valid_records], control, p['evaluation']['bootstrap_draws'], p['evaluation']['bootstrap_seed']) for control in ['minimal', 'common']} for arm in ['baseline', 'variant']}
    return result, records


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('folder', type=Path); parser.add_argument('out', type=Path)
    args = parser.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    result, records = compare(args.folder)
    for name, value in [('summary', result), ('records', records)]:
        (args.out/(name+'.json')).write_text(json.dumps(safe_json(value), indent=2, allow_nan=False)+'\n')
    lines = ['# Frozen held-out scoring comparison', '', result['scope'], '', f"Decision: **{result['decision']}**. Valid pairs: {result['pairs']}/10.", '',
             'Mean comparison and prespecified gates:', '```json', json.dumps(safe_json({k: result[k] for k in ['meanComparison', 'gates']}), indent=2), '```', '',
             '| Participant | Baseline MSE | Variant MSE | Difference | Individual guardrail | Variant uncertainty |', '|---|---:|---:|---:|---|---|']
    for row in result['rows']:
        lines.append(f"| {row['subject']} | {row['baselineMSE']:.6g} | {row['variantMSE']:.6g} | {row['difference']:+.6g} | {row['individualGuardrailPassed']} | {row['variantUncertainty']['passed']} |")
    lines += ['', 'All ten planned subjects, including missing/invalid pairs:', '```json', json.dumps(result['coverage'], indent=2), '```', '',
              'Point accuracy, uncertainty support and engineering completion are separate outcomes. Ten same-source recordings cannot establish stable calibration or real-task superiority. No production change or merge. See summary.json for retention, support, calibration denominators, fixed comparators and reasons.']
    (args.out/'summary.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k: result[k] for k in ['decision', 'pairs', 'expected', 'gates', 'meanComparison']}, indent=2))
    if result['decision'] == 'incomplete':
        raise SystemExit('Incomplete or inconsistent paired evidence; inspect coverage. Scientific worsening alone is not a workflow failure.')


if __name__ == '__main__':
    main()
