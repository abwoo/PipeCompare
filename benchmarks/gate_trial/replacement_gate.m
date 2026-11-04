function g=replacement_gate(points,boot,incumbent,challenger,alpha)
% One-sided basic-bootstrap simultaneous upper bounds for score differences.
% Inputs contain NO ground truth. Lower scores are better. The incumbent and
% challenger are selected upstream; all feasible alternatives enter correction.
g=struct('incumbent',incumbent,'challenger',challenger,'selected',incumbent, ...
 'switched',false,'alpha',alpha,'reason','unassessed','upper',[],'difference',[]);
points=points(:);n=numel(points);
assert(size(boot,1)==n&&ismember(incumbent,1:n)&&ismember(challenger,1:n));
assert(isscalar(alpha)&&alpha>0&&alpha<1);
if challenger==incumbent,g.reason='same recommendation';return;end
if size(boot,2)<100||any(~isfinite(points))||any(~isfinite(boot(:)))
 g.reason='nonfinite or insufficient bootstrap evidence';return;
end
alternatives=setdiff(1:n,incumbent);
d=points-points(incumbent);db=boot-boot(incumbent,:);
% Bootstrap analogue of (true difference - estimated difference).
err=d-db;mx=sort(max(err(alternatives,:),[],1));
at=min(numel(mx),ceil((numel(mx)+1)*(1-alpha)));
q=max(0,mx(at));upper=d+q;
g.difference=d;g.upper=upper;
if upper(challenger)<0
 g.selected=challenger;g.switched=true;g.reason='simultaneous score upper bound below zero';
else
 g.reason='improvement not supported; retain incumbent';
end
end
