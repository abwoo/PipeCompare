function run_reference(out)
D=load(fullfile(out,'input.mat'));fs=D.fs;
[h,m]=pipecompare.eval.Hmm.events(D.X,fs);[f,mf]=pipecompare.eval.Hmm.events(D.F,fs);
r=D.X-movmedian(D.X,round(fs),2);a=abs(r-median(r,2));rob=any(a>5*max(1.4826*median(a,2),eps),1);
w=round(.05*fs);pad=@(x)conv(double(x),ones(1,2*w+1),'same')>0;
hp=pad(h);rp=pad(rob);fp=pad(f);
writetable(table(D.time(:),h(:),hp(:),rp(:),hp(:)|rp(:),fp(:),'VariableNames',{'seconds','hmmRaw','hmm','robust','union','fp1Hmm'}),fullfile(out,'masks.csv'));
meta=struct('hmmValid',m.valid,'fp1HmmValid',mf.valid,'hmmModel',m,'fp1Model',mf);
fid=fopen(fullfile(out,'models.json'),'w');fprintf(fid,'%s\n',jsonencode(meta));fclose(fid);
assert(m.valid&&mf.valid,'Invalid HMM fit; results retained but validation incomplete');
fprintf('EYE_REFERENCE_NATIVE_EXECUTION_PASSED\n');
end
