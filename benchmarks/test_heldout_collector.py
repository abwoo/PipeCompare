"""Synthetic collector fixtures; these are not participant performance results."""
import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import compare_heldout as c


class HeldoutCollectorTests(unittest.TestCase):
    def setUp(self):
        self.p = json.loads((c.ROOT/'heldout_protocol.json').read_text())
        self.digest = hashlib.sha256((c.ROOT/'heldout_protocol.json').read_bytes()).hexdigest()

    def rows(self, value=.8):
        return [dict(subject=s, baselineMSE=1., variantMSE=value, baselineRawMSE=1., variantRawMSE=value,
                     individualGuardrailPassed=True, variantUncertainty=dict(passed=True)) for s in range(11, 21)]

    def test_improvement(self):
        result = c.decide(self.rows(), self.p)
        self.assertEqual(result['decision'], 'retain-heldout-supported')
        self.assertAlmostEqual(result['meanComparison']['relativeReduction'], .2)
        self.assertTrue(all(result['gates'].values()))

    def test_worsening_and_raw_guardrail(self):
        self.assertEqual(c.decide(self.rows(1.2), self.p)['decision'], 'reject-variant-retain-original')
        rows = self.rows()
        for r in rows:
            r['variantRawMSE'] = 1.2
        self.assertEqual(c.decide(rows, self.p)['decision'], 'reject-variant-retain-original')

    def test_inconclusive(self):
        self.assertEqual(c.decide(self.rows(1.), self.p)['decision'], 'inconclusive-no-promotion')
        rows = self.rows(); rows[0]['variantUncertainty']['passed'] = False
        self.assertEqual(c.decide(rows, self.p)['decision'], 'inconclusive-no-promotion')
        rows = self.rows(); rows[0]['individualGuardrailPassed'] = False
        self.assertEqual(c.decide(rows, self.p)['decision'], 'inconclusive-no-promotion')

    def test_missing_subject_is_incomplete(self):
        self.assertEqual(c.decide(self.rows()[:-1], self.p)['decision'], 'incomplete')
        with tempfile.TemporaryDirectory() as d:
            result, _ = c.compare(Path(d))
            self.assertEqual(result['decision'], 'incomplete')
            self.assertEqual(len(result['coverage']), 10)
            self.assertEqual(result['uncertainty']['variant']['assessable'], 0)

    def record(self, subject, arm):
        entry = next(f for f in self.p['source']['files'] if f['subject'] == subject)
        e = self.p['excerpt']
        x = dict(available=True, estimate=[2., 2.], gain=1., correctedEstimate=[2., 2.], actualError=[1., 1.], rawError=[1., 1.], conditionalError=[1., 1.], mse=1., rawMse=1., conditionalMse=1., reportedError=[1., 1.], reportedMse=1., pipeline='fixture')
        support = [dict(condition=name, referenceN=10, referenceIds=list(range(10)), scores=[2.]*10, weights=[1.]*10, flagged=[False]*10, kept=10, retention=1., positiveN=10, flaggedN=0, unflaggedN=10, effectiveN=10., effectiveFlaggedN=0., effectiveUnflaggedN=10., eyeBiasAssessed=True, eyeBiasFinite=True) for name in ['target', 'standard']]
        source = dict(subject=subject, source_sha256=entry['sha256'], source_url=self.p['source']['base_url']+entry['path'], protocol_sha256=self.digest, channels=e['channels'], sampling_rate=e['sampling_rate'], samples=9600, offset_seconds=0, duration_seconds=60)
        return dict(subject=subject, algorithmCommit=self.p['algorithm_commit'] if arm == 'baseline' else 'fixture-commit', executionCommit='fixture-commit', workflowRun='fixture-run', runAttempt='1', measureSha256=self.p[f'{arm}_measure_sha256'], protocolSha256=self.digest, searchRngSeed=self.p['search_rng_seed'], matlab='R2026a', eeglab='2026.0.0', eeglabSha='fixture-eeglab', source=source, completed=True, errors=[], notRun=0, feasible=2, searchOptions=dict(icaSeeds=2, refineCutoffs=True, eogCheck=True, signalTemplate='data'), populationTruth=[1., 1.], realizedTruth=[1., 1.], erp=copy.deepcopy(x), minimal=copy.deepcopy(x), common=copy.deepcopy(x), replayMaxError=0., exportMaxError=0., selectedSupport=support, bootstrapAudit=dict(draws=1999, feasibleIds=[1, 2], nonfiniteByCandidate=[0, 0], nonfiniteColumns=0), icaSeedCheck=dict(same=True))

    def test_record_recomputation_and_tampering(self):
        r = self.record(11, 'baseline'); c.validate_record(r, self.p, self.digest, 'baseline')
        for field, value in [('protocolSha256', 'wrong'), ('measureSha256', 'wrong')]:
            q = copy.deepcopy(r); q[field] = value
            with self.assertRaises(ValueError):
                c.validate_record(q, self.p, self.digest, 'baseline')
        r['erp']['mse'] = .1
        with self.assertRaises(ValueError):
            c.validate_record(r, self.p, self.digest, 'baseline')

    def test_nonfinite_uncertainty_stays_visible(self):
        r = self.record(11, 'variant'); r['erp']['reportedError'] = [None, 1.]; r['erp']['reportedMse'] = None
        c.validate_record(r, self.p, self.digest, 'variant')
        self.assertFalse(c.uncertainty_screen(r, self.p)['passed'])

    def test_single_feasible_matlab_scalar(self):
        r = self.record(11, 'baseline'); r['feasible'] = 1
        r['bootstrapAudit'].update(feasibleIds=1, nonfiniteByCandidate=0)
        c.validate_record(r, self.p, self.digest, 'baseline')
        r = self.record(11, 'variant'); r['bootstrapAudit']['nonfiniteColumns'] = 1
        self.assertFalse(c.uncertainty_screen(r, self.p)['passed'])

    def test_latest_incomplete_attempt_is_not_hidden(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for subject in range(11, 21):
                for arm in ['baseline', 'variant']:
                    path = root/f'heldout-pair-{subject}-attempt-1'/arm
                    path.mkdir(parents=True)
                    (path/'result.json').write_text(json.dumps(self.record(subject, arm)))
            result, _ = c.compare(root)
            self.assertEqual(result['pairs'], 10)
            self.assertEqual(result['decision'], 'inconclusive-no-promotion')
            (root/'heldout-pair-11-attempt-2').mkdir()
            result, _ = c.compare(root)
            self.assertEqual(result['pairs'], 9)
            self.assertEqual(result['decision'], 'incomplete')
            self.assertFalse(result['coverage'][0]['validPair'])


if __name__ == '__main__':
    unittest.main()
