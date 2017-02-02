function model = PharmacokineticsPoly()
	syms x1 x2 x3 x4 xi1
	syms a1 a2 bb1 bb2 ka kc vm c0 g
	model.Name = 'PharmacokineticsPoly';
	model.Nder = 7;
	model.X = [x1,x2,x3,x4,xi1];
	model.Neq = 5;
	model.G = [0,0,0,0,0];
	model.Noc = 0;
	model.P = [];
	model.Par = [a1,a2,bb1,bb2,ka,kc,vm,c0,g];
	model.Npar = 9;
	model.IC = [c0,0,c0*g,0,1/(c0*ka + ka*kc + c0*g*kc)];
	model.H = [x1 x4];
	model.Nobs = 2;
	model.F = [- a1*(x1 - x2) - ka*vm*x1*xi1,a2*(x1 - x2),- bb1*(x3 - x4) - kc*vm*x3*xi1,bb2*(x3 - x4),ka*xi1^2*(a1*(x1 - x2) + ka*vm*x1*xi1) + kc*xi1^2*(bb1*(x3 - x4) + kc*vm*x3*xi1)];
end