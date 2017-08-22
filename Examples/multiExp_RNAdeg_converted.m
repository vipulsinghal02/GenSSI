function model = multiExp_RNAdeg_converted()
	syms mExp1 CExp1 mExp2 CExp2 mExp3 CExp3
	syms kr kc Et m0Exp1
	model.Name = 'multiExp_RNAdeg_converted';
	model.Nder = 4;
	model.X = [mExp1,CExp1,mExp2,CExp2,mExp3,CExp3];
	model.Neq = 6;
	model.G = [0,0,0,0,0,0];
	model.Noc = 0;
	model.P = [kr,kc,Et,m0Exp1];
	model.Par = [kr,kc,Et,m0Exp1];
	model.Npar = 4;
	model.IC = [m0Exp1,0,m0Exp1/2,0,2*m0Exp1,0];
	model.H = [mExp1,mExp2,mExp3];
	model.Nobs = 3;
	model.F = [CExp1*kr + 100*mExp1*(CExp1 - Et),- 100*mExp1*(CExp1 - Et) - CExp1*(kc + kr),CExp2*kr + 100*mExp2*(CExp2 - Et),- 100*mExp2*(CExp2 - Et) - CExp2*(kc + kr),CExp3*kr + 100*mExp3*(CExp3 - Et),- 100*mExp3*(CExp3 - Et) - CExp3*(kc + kr)];
end