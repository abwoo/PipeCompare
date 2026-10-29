function m = write_slice_manifest(outDir, seed, sliceIndex, checkpointDir, artifactVariance, matlabVersion, eeglabVersion, eeglabSha)
%WRITE_SLICE_MANIFEST Write a slice manifest for a sliced search.
%
% Report-only: this writes a manifest, it does not resume anything
% and it does not change any pick. A manifest is not a completed
% search, and a sliced run is not a full run.
%
% The manifest records exactly what a later resume must check:
% the checkpoint directory, the artifactVariance in force, and the
% toolchain versions. It also records, in plain text, that an
% exact resume is NOT claimed: Executor.resume re-derives adaptive
% leaves beyond the manifest (see Executor.resume comments), so
% adaptiveLeavesRederived is true and exactResume is false.
if ~ischar(checkpointDir) || isempty(checkpointDir)
    error('PipeCompare:Slice', 'slice manifest requires a non-empty checkpointDir.');
end
if ~ischar(artifactVariance) || isempty(artifactVariance)
    error('PipeCompare:Slice', 'slice manifest requires an explicit artifactVariance; a default must never be inferred on resume.');
end
if ~ischar(outDir) || isempty(outDir)
    error('PipeCompare:Slice', 'slice manifest requires a non-empty outDir.');
end
if ~isfolder(outDir), mkdir(outDir); end
m = struct('seed', seed, 'sliceIndex', sliceIndex, ...
    'checkpointDir', checkpointDir, 'artifactVariance', artifactVariance, ...
    'matlabVersion', matlabVersion, 'eeglabVersion', eeglabVersion, ...
    'eeglabSha', eeglabSha, 'adaptiveLeavesRederived', true, ...
    'exactResume', false, 'mode', 'slice-diagnostic', ...
    'note', 'Sliced validation run only: never merged into a previous paired mean.');
fid = fopen(fullfile(outDir, 'slice_manifest.json'), 'w');
if fid < 0
    error('PipeCompare:Slice', 'Cannot write slice manifest in %s.', outDir);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(m));
clear cleanup;
end
