classdef IcaCheckpoint
    % Exact, search-local ICA memoization. No approximation or early stopping.
    methods (Static)
        function [E2,coms,info] = run(E,in,ctx,folder,identity,compute)
            assert(~isempty(identity),'PipeCompare:Checkpoint','Missing search identity.');
            code={};
            for name={'pipecompare.run.Steps','pipecompare.run.IcaCheckpoint','pop_runica','runica'}
                file=which(name{1});
                assert(isfile(file),'PipeCompare:Checkpoint','Missing ICA dependency %s.',name{1});
                code{end+1}=fileread(file); %#ok<AGROW>
            end
            token={1,identity,E,in,pipecompare.utils.fieldOr(ctx,'highpass',0), ...
                pipecompare.utils.fieldOr(ctx,'icaSeed',0),version,code};
            md=java.security.MessageDigest.getInstance('SHA-256');
            md.update(getByteStreamFromArray(token));
            key=lower(reshape(dec2hex(typecast(md.digest(),'uint8'))',1,[]));
            cacheDir=fullfile(folder,'ica_cache');if ~isfolder(cacheDir),mkdir(cacheDir);end
            file=fullfile(cacheDir,[key '.mat']);
            if isfile(file)
                try,d=load(file);catch ME
                    error('PipeCompare:Checkpoint','Cannot load ICA checkpoint: %s',ME.message);
                end
                assert(all(isfield(d,{'key','identity','E2','coms','info','rngState'})), ...
                    'PipeCompare:Checkpoint','Incomplete ICA checkpoint.');
                assert(isfield(d,'key')&&strcmp(d.key,key)&&isfield(d,'identity')&&strcmp(d.identity,identity), ...
                    'PipeCompare:Checkpoint','ICA checkpoint identity mismatch.');
                E2=d.E2;coms=d.coms;info=d.info;rng(d.rngState);
                pipecompare.utils.log('Reused exact ICA checkpoint (seed %d).',pipecompare.utils.fieldOr(ctx,'icaSeed',0));
                return;
            end
            [E2,coms,info]=compute();rngState=rng; %#ok<NASGU>
            tmp=[tempname(cacheDir) '.mat'];
            clean=onCleanup(@()pipecompare.run.IcaCheckpoint.removeTemp(tmp)); %#ok<NASGU>
            save(tmp,'key','identity','E2','coms','info','rngState','-v7.3');
            movefile(tmp,file,'f');
        end
        function removeTemp(file)
            if isfile(file),delete(file);end
        end
    end
end
