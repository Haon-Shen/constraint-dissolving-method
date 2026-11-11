function z = proj_box(x, l, u)
% PROJ_BOX  Projection onto the box [l,u].
z = max(l, min(x, u));
end
