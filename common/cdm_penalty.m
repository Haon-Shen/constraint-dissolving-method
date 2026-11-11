function [h, gh] = cdm_penalty(z, beta, alpha, prob)
% CDM_PENALTY  Constraint dissolving penalty function and its exact gradient,
%
%       h_beta(x,s) = f(A_x(x,s)) + beta/2*||G(x)-s||^2,     z = [x; s],
%
% with H(x,s) = G(x)-s and the constraint dissolving mapping
%
%       A_x(x,s) = x - Q_D(x)*JG(x)'*( M(x,s) \ H(x,s) ),
%       M(x,s)   = JG(x)*Q_D(x)*JG(x)' + Q_C(s) + alpha*||H(x,s)||^2*I.
%
% Required fields of prob (n = dim x, m = dim s):
%   f(x), gradf(x)     objective and its gradient
%   G(x), JG(x)        constraint mapping and its (m x n) Jacobian
%   hessG(x,lam,u)     (sum_i lam_i*Hess G_i(x))*u;  set to [] if G is affine
%   QD(x,U)            Q_D(x)*U for U with n rows
%   dQD(x,u1,u2)       gradient w.r.t. x of u1'*Q_D(x)*u2
%   QC(s,W)            Q_C(s)*W for W with m rows
%   dQC(s,w1,w2)       gradient w.r.t. s of w1'*Q_C(s)*w2

n = prob.n;  m = prob.m;
x = z(1:n);  s = z(n+1:n+m);

JG  = prob.JG(x);
H   = prob.G(x) - s;
nH2 = H'*H;

QJt = prob.QD(x, JG');                          % Q_D(x)*JG(x)'
M   = JG*QJt + prob.QC(s, eye(m)) + alpha*nH2*eye(m);
M   = (M + M')/2;
w   = M\H;
Ax  = x - QJt*w;

h = prob.f(Ax) + 0.5*beta*nH2;

if nargout > 1
	% With v = gradf(A_x) fixed, v'*A_x = v'*x - psi, psi := a'*(M\H), a := JG*Q_D*v.
	% Below, y = M\a and the gradient of psi is assembled term by term.
	v    = prob.gradf(Ax);
	y    = M\(QJt'*v);
	yw   = y'*w;
	JGtw = JG'*w;
	JGty = JG'*y;
	JGtH = JG'*H;

	gpsi_x = prob.dQD(x, JGtw, v) + JGty ...
	       - prob.dQD(x, JGty, JGtw) - 2*alpha*yw*JGtH;
	if ~isempty(prob.hessG)
		gpsi_x = gpsi_x + prob.hessG(x, w, prob.QD(x, v)) ...
		       - prob.hessG(x, y, QJt*w) - prob.hessG(x, w, QJt*y);
	end
	gpsi_s = -y - prob.dQC(s, y, w) + 2*alpha*yw*H;

	gh = [v - gpsi_x + beta*JGtH;
	      -gpsi_s - beta*H];
end
end
