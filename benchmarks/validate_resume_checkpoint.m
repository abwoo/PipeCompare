function out = validate_resume_checkpoint(manifest, currentMatlabVersion, currentEeglabVersion, artifactVariance)
%VALIDATE_RESUME_CHECKPOINT Strict resume guard (diagnostic only).
%
% This guard is stricter than Executor.resume itself: Executor.resume
% only WARNS on a version mismatch and resumes old checkpoints under
% the current artifactVariance default (see Executor.resume and the
% Round 6 search-adaptive notes). This guard ERRORS instead, so a
% sliced run can never silently mix versions or variance models.
%
% It does NOT resume anything and does NOT claim an exact resume:
% even when it passes, adaptive leaves beyond the manifest are
% re-derived by Executor.resume, so out.exactResume is false and
% out.adaptiveLeavesRederived is true. Passing this guard is a
% precondition for a resume attempt, not proof that a resume is
% correct; the resume itself still needs MATLAB verification in CI.
if ~isstruct(manifest)
    error('PipeCompare:Slice', 'resume checkpoint manifest must be a struct.');
end
for f = {'checkpointDir', 'artifactVariance', 'matlabVersion', 'eeglabVersion'}
    if ~isfield(manifest, f{1}) || isempty(manifest.(f{1}))
        error('PipeCompare:Slice', 'resume checkpoint manifest is missing field %s.', f{1});
    end
end
if ~isfolder(manifest.checkpointDir)
    error('PipeCompare:Slice', 'resume checkpoint directory %s does not exist.', manifest.checkpointDir);
end
if ~strcmp(manifest.matlabVersion, currentMatlabVersion)
    error('PipeCompare:Slice', 'resume version mismatch: manifest MATLAB %s, current %s; versions are pinned, not warned.', manifest.matlabVersion, currentMatlabVersion);
end
if ~strcmp(manifest.eeglabVersion, currentEeglabVersion)
    error('PipeCompare:Slice', 'resume version mismatch: manifest EEGLAB %s, current %s; versions are pinned, not warned.', manifest.eeglabVersion, currentEeglabVersion);
end
if ~ischar(artifactVariance) || isempty(artifactVariance) || ~strcmp(manifest.artifactVariance, artifactVariance)
    error('PipeCompare:Slice', 'resume artifactVariance mismatch or missing: manifest %s, requested %s; variance must be recorded explicitly, never defaulted.', manifest.artifactVariance, artifactVariance);
end
out = struct('ok', true, 'exactResume', false, 'adaptiveLeavesRederived', true, ...
    'note', 'Guard passed: versions and artifactVariance match the manifest. Executor.resume still re-derives adaptive leaves beyond the manifest; this is not an exact resume.');
end
