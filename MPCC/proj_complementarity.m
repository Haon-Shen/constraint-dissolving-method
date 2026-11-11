function p = proj_complementarity(x)
% PROJ_COMPLEMENTARITY  One element of the projection of x = [y; z] onto
%
%       {(y,z) : y >= 0, z >= 0, y_i*z_i = 0 for all i}.
%
% For every pair, the candidates (max(y_i,0), 0) and (0, max(z_i,0)) are
% compared; ties are resolved in favor of the first one.
l = numel(x)/2;
y = x(1:l);  z = x(l+1:2*l);
yp = max(y, 0);  zp = max(z, 0);

keep_y = (y - yp).^2 + z.^2 <= y.^2 + (z - zp).^2;
p = [yp.*keep_y; zp.*(~keep_y)];
end
