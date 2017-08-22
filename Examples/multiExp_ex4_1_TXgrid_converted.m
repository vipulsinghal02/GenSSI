function model = multiExp_ex4_1_TXgrid_converted()
	syms m1Exp1 C1Exp1 D1Exp1 PD1Exp1 PD1NExp1 PD1tExp1 m2Exp1 C2Exp1 D2Exp1 PD2Exp1 PD2NExp1 PD2tExp1 PExp1 EExp1 NExp1 m1Exp2 C1Exp2 D1Exp2 PD1Exp2 PD1NExp2 PD1tExp2 m2Exp2 C2Exp2 D2Exp2 PD2Exp2 PD2NExp2 PD2tExp2 PExp2 EExp2 NExp2 m1Exp3 C1Exp3 D1Exp3 PD1Exp3 PD1NExp3 PD1tExp3 m2Exp3 C2Exp3 D2Exp3 PD2Exp3 PD2NExp3 PD2tExp3 PExp3 EExp3 NExp3 m1Exp4 C1Exp4 D1Exp4 PD1Exp4 PD1NExp4 PD1tExp4 m2Exp4 C2Exp4 D2Exp4 PD2Exp4 PD2NExp4 PD2tExp4 PExp4 EExp4 NExp4
	syms k_fm krm ktx kcm kfp krp kterm krn kfn E0 P0 D10ex4 D20ex4 N0 alph1 alph2
	model.Name = 'multiExp_ex4_1_TXgrid_converted';
	model.Nder = 8;
	model.X = [m1Exp1,C1Exp1,D1Exp1,PD1Exp1,PD1NExp1,PD1tExp1,m2Exp1,C2Exp1,D2Exp1,PD2Exp1,PD2NExp1,PD2tExp1,PExp1,EExp1,NExp1,m1Exp2,C1Exp2,D1Exp2,PD1Exp2,PD1NExp2,PD1tExp2,m2Exp2,C2Exp2,D2Exp2,PD2Exp2,PD2NExp2,PD2tExp2,PExp2,EExp2,NExp2,m1Exp3,C1Exp3,D1Exp3,PD1Exp3,PD1NExp3,PD1tExp3,m2Exp3,C2Exp3,D2Exp3,PD2Exp3,PD2NExp3,PD2tExp3,PExp3,EExp3,NExp3,m1Exp4,C1Exp4,D1Exp4,PD1Exp4,PD1NExp4,PD1tExp4,m2Exp4,C2Exp4,D2Exp4,PD2Exp4,PD2NExp4,PD2tExp4,PExp4,EExp4,NExp4];
	model.Neq = 60;
	model.G = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
	model.Noc = 0;
	model.P = [ktx,krp,kterm,krn,P0];
	model.Par = [ktx,krp,kterm,krn,P0];
	model.Npar = 5;
	model.IC = [0,0,D10ex4/2,0,0,0,0,0,D20ex4/2,0,0,0,P0,E0,N0,0,0,D10ex4/2,0,0,0,0,0,D20ex4,0,0,0,P0,E0,N0,0,0,D10ex4,0,0,0,0,0,D20ex4/2,0,0,0,P0,E0,N0,0,0,D10ex4,0,0,0,0,0,D20ex4,0,0,0,P0,E0,N0];
	model.H = [m1Exp1,m2Exp1,m1Exp2,m2Exp2,m1Exp3,m2Exp3,m1Exp4,m2Exp4];
	model.Nobs = 8;
	model.F = [C1Exp1*krm + PD1NExp1*ktx - EExp1*k_fm*m1Exp1,EExp1*k_fm*m1Exp1 - C1Exp1*krm - C1Exp1*kcm,PD1Exp1*krp + PD1tExp1*kterm - D1Exp1*PExp1*kfp,PD1NExp1*krn - PD1Exp1*krp + D1Exp1*PExp1*kfp - NExp1*PD1Exp1*kfn,NExp1*PD1Exp1*kfn - PD1NExp1*krn - PD1NExp1*alph1*ktx,PD1NExp1*ktx - PD1tExp1*kterm,C2Exp1*krm + PD2NExp1*ktx - EExp1*k_fm*m2Exp1,EExp1*k_fm*m2Exp1 - C2Exp1*krm - C2Exp1*kcm,PD2Exp1*krp + PD2tExp1*kterm - D2Exp1*PExp1*kfp,PD2NExp1*krn - PD2Exp1*krp + D2Exp1*PExp1*kfp - NExp1*PD2Exp1*kfn,NExp1*PD2Exp1*kfn - PD2NExp1*krn - PD2NExp1*alph2*ktx,PD2NExp1*ktx - PD2tExp1*kterm,PD1Exp1*krp + PD2Exp1*krp + PD1tExp1*kterm + PD2tExp1*kterm - D1Exp1*PExp1*kfp - D2Exp1*PExp1*kfp,C1Exp1*krm - C2Exp1*kcm - C1Exp1*kcm + C2Exp1*krm - EExp1*k_fm*m1Exp1 - EExp1*k_fm*m2Exp1,PD1NExp1*krn + PD2NExp1*krn - NExp1*PD1Exp1*kfn - NExp1*PD2Exp1*kfn,C1Exp2*krm + PD1NExp2*ktx - EExp2*k_fm*m1Exp2,EExp2*k_fm*m1Exp2 - C1Exp2*krm - C1Exp2*kcm,PD1Exp2*krp + PD1tExp2*kterm - D1Exp2*PExp2*kfp,PD1NExp2*krn - PD1Exp2*krp + D1Exp2*PExp2*kfp - NExp2*PD1Exp2*kfn,NExp2*PD1Exp2*kfn - PD1NExp2*krn - PD1NExp2*alph1*ktx,PD1NExp2*ktx - PD1tExp2*kterm,C2Exp2*krm + PD2NExp2*ktx - EExp2*k_fm*m2Exp2,EExp2*k_fm*m2Exp2 - C2Exp2*krm - C2Exp2*kcm,PD2Exp2*krp + PD2tExp2*kterm - D2Exp2*PExp2*kfp,PD2NExp2*krn - PD2Exp2*krp + D2Exp2*PExp2*kfp - NExp2*PD2Exp2*kfn,NExp2*PD2Exp2*kfn - PD2NExp2*krn - PD2NExp2*alph2*ktx,PD2NExp2*ktx - PD2tExp2*kterm,PD1Exp2*krp + PD2Exp2*krp + PD1tExp2*kterm + PD2tExp2*kterm - D1Exp2*PExp2*kfp - D2Exp2*PExp2*kfp,C1Exp2*krm - C2Exp2*kcm - C1Exp2*kcm + C2Exp2*krm - EExp2*k_fm*m1Exp2 - EExp2*k_fm*m2Exp2,PD1NExp2*krn + PD2NExp2*krn - NExp2*PD1Exp2*kfn - NExp2*PD2Exp2*kfn,C1Exp3*krm + PD1NExp3*ktx - EExp3*k_fm*m1Exp3,EExp3*k_fm*m1Exp3 - C1Exp3*krm - C1Exp3*kcm,PD1Exp3*krp + PD1tExp3*kterm - D1Exp3*PExp3*kfp,PD1NExp3*krn - PD1Exp3*krp + D1Exp3*PExp3*kfp - NExp3*PD1Exp3*kfn,NExp3*PD1Exp3*kfn - PD1NExp3*krn - PD1NExp3*alph1*ktx,PD1NExp3*ktx - PD1tExp3*kterm,C2Exp3*krm + PD2NExp3*ktx - EExp3*k_fm*m2Exp3,EExp3*k_fm*m2Exp3 - C2Exp3*krm - C2Exp3*kcm,PD2Exp3*krp + PD2tExp3*kterm - D2Exp3*PExp3*kfp,PD2NExp3*krn - PD2Exp3*krp + D2Exp3*PExp3*kfp - NExp3*PD2Exp3*kfn,NExp3*PD2Exp3*kfn - PD2NExp3*krn - PD2NExp3*alph2*ktx,PD2NExp3*ktx - PD2tExp3*kterm,PD1Exp3*krp + PD2Exp3*krp + PD1tExp3*kterm + PD2tExp3*kterm - D1Exp3*PExp3*kfp - D2Exp3*PExp3*kfp,C1Exp3*krm - C2Exp3*kcm - C1Exp3*kcm + C2Exp3*krm - EExp3*k_fm*m1Exp3 - EExp3*k_fm*m2Exp3,PD1NExp3*krn + PD2NExp3*krn - NExp3*PD1Exp3*kfn - NExp3*PD2Exp3*kfn,C1Exp4*krm + PD1NExp4*ktx - EExp4*k_fm*m1Exp4,EExp4*k_fm*m1Exp4 - C1Exp4*krm - C1Exp4*kcm,PD1Exp4*krp + PD1tExp4*kterm - D1Exp4*PExp4*kfp,PD1NExp4*krn - PD1Exp4*krp + D1Exp4*PExp4*kfp - NExp4*PD1Exp4*kfn,NExp4*PD1Exp4*kfn - PD1NExp4*krn - PD1NExp4*alph1*ktx,PD1NExp4*ktx - PD1tExp4*kterm,C2Exp4*krm + PD2NExp4*ktx - EExp4*k_fm*m2Exp4,EExp4*k_fm*m2Exp4 - C2Exp4*krm - C2Exp4*kcm,PD2Exp4*krp + PD2tExp4*kterm - D2Exp4*PExp4*kfp,PD2NExp4*krn - PD2Exp4*krp + D2Exp4*PExp4*kfp - NExp4*PD2Exp4*kfn,NExp4*PD2Exp4*kfn - PD2NExp4*krn - PD2NExp4*alph2*ktx,PD2NExp4*ktx - PD2tExp4*kterm,PD1Exp4*krp + PD2Exp4*krp + PD1tExp4*kterm + PD2tExp4*kterm - D1Exp4*PExp4*kfp - D2Exp4*PExp4*kfp,C1Exp4*krm - C2Exp4*kcm - C1Exp4*kcm + C2Exp4*krm - EExp4*k_fm*m1Exp4 - EExp4*k_fm*m2Exp4,PD1NExp4*krn + PD2NExp4*krn - NExp4*PD1Exp4*kfn - NExp4*PD2Exp4*kfn];
end