function model = multiExp_ex2_1_TXonly_converted()
	syms mExp1 CExp1 DExp1 PDExp1 PDNExp1 PDtExp1 NExp1 mExp2 CExp2 DExp2 PDExp2 PDNExp2 PDtExp2 NExp2 mExp3 CExp3 DExp3 PDExp3 PDNExp3 PDtExp3 NExp3 mExp4 CExp4 DExp4 PDExp4 PDNExp4 PDtExp4 NExp4
	syms k_fm krm ktx kcm kfp krp kterm krn kfn Et Pt DtEx2 N0ex2 alph
	model.Name = 'multiExp_ex2_1_TXonly_converted';
	model.Nder = 10;
	model.X = [mExp1,CExp1,DExp1,PDExp1,PDNExp1,PDtExp1,NExp1,mExp2,CExp2,DExp2,PDExp2,PDNExp2,PDtExp2,NExp2,mExp3,CExp3,DExp3,PDExp3,PDNExp3,PDtExp3,NExp3,mExp4,CExp4,DExp4,PDExp4,PDNExp4,PDtExp4,NExp4];
	model.Neq = 28;
	model.G = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
	model.Noc = 0;
	model.P = [ktx,krp,kterm,krn,Pt];
	model.Par = [ktx,krp,kterm,krn,Pt];
	model.Npar = 5;
	model.IC = [0,0,DtEx2/2,0,0,0,N0ex2/2,0,0,DtEx2/2,0,0,0,N0ex2,0,0,DtEx2,0,0,0,N0ex2/2,0,0,DtEx2,0,0,0,N0ex2];
	model.H = [mExp1,mExp2,mExp3,mExp4];
	model.Nobs = 4;
	model.F = [CExp1*krm + PDNExp1*ktx + k_fm*mExp1*(CExp1 - Et),- CExp1*kcm - CExp1*krm - k_fm*mExp1*(CExp1 - Et),PDExp1*krp + PDtExp1*kterm - DExp1*kfp*(DExp1 - DtEx2 + Pt),PDNExp1*krn - PDExp1*krp - NExp1*PDExp1*kfn + DExp1*kfp*(DExp1 - DtEx2 + Pt),NExp1*PDExp1*kfn - PDNExp1*krn - PDNExp1*alph*ktx,PDNExp1*ktx - PDtExp1*kterm,PDNExp1*krn - NExp1*PDExp1*kfn,CExp2*krm + PDNExp2*ktx + k_fm*mExp2*(CExp2 - Et),- CExp2*kcm - CExp2*krm - k_fm*mExp2*(CExp2 - Et),PDExp2*krp + PDtExp2*kterm - DExp2*kfp*(DExp2 - DtEx2 + Pt),PDNExp2*krn - PDExp2*krp - NExp2*PDExp2*kfn + DExp2*kfp*(DExp2 - DtEx2 + Pt),NExp2*PDExp2*kfn - PDNExp2*krn - PDNExp2*alph*ktx,PDNExp2*ktx - PDtExp2*kterm,PDNExp2*krn - NExp2*PDExp2*kfn,CExp3*krm + PDNExp3*ktx + k_fm*mExp3*(CExp3 - Et),- CExp3*kcm - CExp3*krm - k_fm*mExp3*(CExp3 - Et),PDExp3*krp + PDtExp3*kterm - DExp3*kfp*(DExp3 -  DtEx2 + Pt),PDNExp3*krn - PDExp3*krp - NExp3*PDExp3*kfn + DExp3*kfp*(DExp3 -  DtEx2 + Pt),NExp3*PDExp3*kfn - PDNExp3*krn - PDNExp3*alph*ktx,PDNExp3*ktx - PDtExp3*kterm,PDNExp3*krn - NExp3*PDExp3*kfn,CExp4*krm + PDNExp4*ktx + k_fm*mExp4*(CExp4 - Et),- CExp4*kcm - CExp4*krm - k_fm*mExp4*(CExp4 - Et),PDExp4*krp + PDtExp4*kterm - DExp4*kfp*(DExp4 -  DtEx2 + Pt),PDNExp4*krn - PDExp4*krp - NExp4*PDExp4*kfn + DExp4*kfp*(DExp4 -  DtEx2 + Pt),NExp4*PDExp4*kfn - PDNExp4*krn - PDNExp4*alph*ktx,PDNExp4*ktx - PDtExp4*kterm,PDNExp4*krn - NExp4*PDExp4*kfn];
end