function check_replacement_gate()
b=repmat([2;1;3],1,1999)+sin(1:1999);
g=replacement_gate([2;1;3],b,1,2,.02);assert(g.switched&&g.selected==2);
h=replacement_gate([12;11;13],b+10,1,2,.02);assert(max(abs(g.upper-h.upper))<1e-12);
b(2,:)=1+4*sin(1:1999);g=replacement_gate([2;1;3],b,1,2,.02);assert(~g.switched);
b(2,1)=NaN;g=replacement_gate([2;1;3],b,1,2,.02);assert(~g.switched&&strcmp(g.reason,'nonfinite or insufficient bootstrap evidence'));
g=replacement_gate([2;1;3],zeros(3,99),1,2,.02);assert(~g.switched);
% A noisy extra challenger increases simultaneous uncertainty for the chosen one.
b=repmat([2;1;3],1,1999);b(3,:)=3+10*sin(1:1999);
g=replacement_gate([2;1;3],b,1,2,.02);assert(~g.switched);
fprintf('REPLACEMENT_GATE_CHECKS_PASSED\n');
end
