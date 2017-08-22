function model = multiExp_ex2_3_TXonly_converted()
	syms mExp1 EExp1 CExp1 PExp1 DExp1 PDExp1 PDNExp1 PDtExp1 NExp1 mExp2 EExp2 CExp2 PExp2 DExp2 PDExp2 PDNExp2 PDtExp2 NExp2 mExp3 EExp3 CExp3 PExp3 DExp3 PDExp3 PDNExp3 PDtExp3 NExp3 mExp4 EExp4 CExp4 PExp4 DExp4 PDExp4 PDNExp4 PDtExp4 NExp4
	syms k_fm krm ktx kcm kfp krp kterm krn kfn E0 P0 DtEx2 N0ex2 alph
	model.Name = 'multiExp_ex2_3_TXonly_converted';
	model.Nder = 7;
	model.X = [mExp1,EExp1,CExp1,PExp1,DExp1,PDExp1,PDNExp1,PDtExp1,NExp1,mExp2,EExp2,CExp2,PExp2,DExp2,PDExp2,PDNExp2,PDtExp2,NExp2,mExp3,EExp3,CExp3,PExp3,DExp3,PDExp3,PDNExp3,PDtExp3,NExp3,mExp4,EExp4,CExp4,PExp4,DExp4,PDExp4,PDNExp4,PDtExp4,NExp4];
	model.Neq = 36;
	model.G = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
	model.Noc = 0;
	model.P = [ktx,krp,kterm,krn,P0];
	model.Par = [ktx,krp,kterm,krn,P0];
	model.Npar = 5;
	model.IC = [0,E0,0,P0,DtEx2/2,0,0,0,N0ex2/2,0,E0,0,P0,DtEx2/2,0,0,0,N0ex2,0,E0,0,P0,DtEx2,0,0,0,N0ex2/2,0,E0,0,P0,DtEx2,0,0,0,N0ex2];
	model.H = [mExp1,mExp2,mExp3,mExp4];
	model.Nobs = 4;
	model.F = [CExp1*krm + PDNExp1*ktx - EExp1*k_fm*mExp1,CExp1*kcm + CExp1*krm - EExp1*k_fm*mExp1,EExp1*k_fm*mExp1 - CExp1*krm - CExp1*kcm,PDExp1*krp + PDtExp1*kterm - DExp1*PExp1*kfp,PDExp1*krp + PDtExp1*kterm - DExp1*PExp1*kfp,PDNExp1*krn - PDExp1*krp + DExp1*PExp1*kfp - NExp1*PDExp1*kfn,NExp1*PDExp1*kfn - PDNExp1*krn - PDNExp1*alph*ktx,PDNExp1*ktx - PDtExp1*kterm,PDNExp1*krn - NExp1*PDExp1*kfn,CExp2*krm + PDNExp2*ktx - EExp2*k_fm*mExp2,CExp2*kcm + CExp2*krm - EExp2*k_fm*mExp2,EExp2*k_fm*mExp2 - CExp2*krm - CExp2*kcm,PDExp2*krp + PDtExp2*kterm - DExp2*PExp2*kfp,PDExp2*krp + PDtExp2*kterm - DExp2*PExp2*kfp,PDNExp2*krn - PDExp2*krp + DExp2*PExp2*kfp - NExp2*PDExp2*kfn,NExp2*PDExp2*kfn - PDNExp2*krn - PDNExp2*alph*ktx,PDNExp2*ktx - PDtExp2*kterm,PDNExp2*krn - NExp2*PDExp2*kfn,CExp3*krm + PDNExp3*ktx - EExp3*k_fm*mExp3,CExp3*kcm + CExp3*krm - EExp3*k_fm*mExp3,EExp3*k_fm*mExp3 - CExp3*krm - CExp3*kcm,PDExp3*krp + PDtExp3*kterm - DExp3*PExp3*kfp,PDExp3*krp + PDtExp3*kterm - DExp3*PExp3*kfp,PDNExp3*krn - PDExp3*krp + DExp3*PExp3*kfp - NExp3*PDExp3*kfn,NExp3*PDExp3*kfn - PDNExp3*krn - PDNExp3*alph*ktx,PDNExp3*ktx - PDtExp3*kterm,PDNExp3*krn - NExp3*PDExp3*kfn,CExp4*krm + PDNExp4*ktx - EExp4*k_fm*mExp4,CExp4*kcm + CExp4*krm - EExp4*k_fm*mExp4,EExp4*k_fm*mExp4 - CExp4*krm - CExp4*kcm,PDExp4*krp + PDtExp4*kterm - DExp4*PExp4*kfp,PDExp4*krp + PDtExp4*kterm - DExp4*PExp4*kfp,PDNExp4*krn - PDExp4*krp + DExp4*PExp4*kfp - NExp4*PDExp4*kfn,NExp4*PDExp4*kfn - PDNExp4*krn - PDNExp4*alph*ktx,PDNExp4*ktx - PDtExp4*kterm,PDNExp4*krn - NExp4*PDExp4*kfn];
end