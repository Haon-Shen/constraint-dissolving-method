function p = proj_budget(x, lo, up)
% PROJ_BUDGET  Projection of x onto the convex set {x >= 0 : lo <= e'*x <= up}.
p = max(x, 0);
sp = sum(p);
if sp > up
	p = proj_simplex(x, up);
elseif sp < lo
	p = proj_simplex(x, lo);
end
end

function p = proj_simplex(x, r)
% projection onto {x >= 0 : e'*x = r}
u = sort(x, 'descend');
css = cumsum(u) - r;
k = find(u - css./(1:numel(u))' > 0, 1, 'last');
p = max(x - css(k)/k, 0);
end
