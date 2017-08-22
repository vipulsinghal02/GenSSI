 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 TL IDENTIFIABILITY IN TXTL (v1)                         %%%
%%%                           2017                                          %%% 
%%%   Vipul Singhal, California Institute of Technology                     %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = ex3_1_TL_only()
    model.Name='ex3_1_TL_only';

    % Here we fix the RNA degradation numbers and the TX numbers
    % (assuming they are known from a previous experiment), and check the 
    % identifiability of the translation parameters. 
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     DECLARE SYMBOLIC VARIABLES:                          %
    %               - state variables                          %
    %               - parameters of the model                  %
    %               - initial state, if not known              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % SPECIES 
    % m = mRNA
    % E = RNAse enzyme
    % P = RNA Pol
    % D = free DNA
    % N = unused nucleotides
    % PDt = tx termination complex
    % PD = P and D complex
    % PDN = P, D and N complex
    % R = ribosome
    % AA = amino acids
    % mR, mRAA = bound complexes
    % mRt = translation termination complex. 
    % mE, mRE, mRAAE, mRtE = rna degradation bound complexes
    % 
    % BIOCHEMICAL EQUATIONS
    % 
    % mrna degradation (previously characterized)
    % (m || mR || mRAA || mRt) + E <-> (mE || mRE || mRAAE || mRtE)  -> E (+ R || + R + AA || + R) (k_fm, krm, kcm)
    % E + mE + mRE + mRAAE + mRtE = Et (Et known) we can choose to use this
    % law or not... since non of these are being manipulated, it will not
    % matter if it is used. can try both cases and see what happens. 
    % 
    % transcription (previously characterized)
    % P + D <-> PD (krp, kfp; kfp known)
    % PD + N <-> PDN (krn, kfn; kfn known)
    % PDN -> PD (ktx*(alph - 1))
    % PDN -> PDt + m (ktx)
    % PDt -> P + D (kterm)
    % alph base pairs, known
    % P + PD + PDN + PDt = Pt (unknown, parameter of interest)
    % D + PD + PDN + PDt = Dt (known, initial condition control)
    % these two conservation laws are actually the same. 
    % 
    % translation
    % m + R <-> mR
    % mR + AA <-> mRAA
    % mRAA -> mR
    % mRAA -> mRt + G
    % mRt -> m + R
    % 
    % in ODE form (we do not use conservation laws to reduce the ODEs.) 
    %
    % dP = -kfp * P * D + krp * PD + kterm * PDt   
    % dD = -kfp * P * D + krp * PD + kterm * PDt
    % dPD = kfp * P * D - krp * PD + krn * PDN - kfn * PD * N
    % dPDN = - krn * PDN + kfn * PD * N - alph * ktx * PDN
    % dPDt = ktx * PDN - kterm * PDt
    % dN = - kfn * PD * N + krn * PDN
    % dm = -k_fm * m * E + krm * mE + ktx * PDN - kfR * m * R + krR * mR + ktermTL * mRt
    % dR = - kfR * m * R + krR * mR + ktermTL * mRt + kcm * (mRE + mRAAE + mRtE)
    % dmR = kfR * m * R - krR * mR - kfAA * mR * AA + krAA * mRAA + (betaparam-1) * kTL * mRAA -k_fm * mR * E + krm * mRE
    % dmRAA = kfAA * mR * AA - krAA * mRAA - betaparam * kTL * mRAA - k_fm * mRAA * E + krm * mRAAE
    % dmRt = kTL * mRAA - ktermTL * mRt - k_fm * mRt * E + krm * mRtE
    % dAA = -kfAA * mR * AA + krAA * mRAA + kcm * mRAAE 
    % dG = kTL * mRAA 
    % dmE = k_fm * m * E - (krm + kcm) * mE
    % dmRE = k_fm * mR * E - (krm + kcm) * mRE
    % dmRAAE = k_fm * mRAA * E - (krm + kcm) * mRAAE
    % dmRtE = k_fm * mRt * E - (krm + kcm) * mRtE
    % dE = -k_fm * E * (m + mR + mRAA + mRt) + (krm + kcm) * (mE + mRE + mRAAE + mRtE)
    
    
    
    
    % dC = k_fm * m * (Et - C) - krm * C - kcm * C

    
    
    
    
    
    
    
    %
    % UNKNOWN PARAMETERS 
    % krp, krn, ktx, kterm 
    % 
    % KNOWN PARAMETERS
    % k_fm krm kcm kfp kfn 
    %
    % UNKNOWN IC: Pt = P(0)
    % 
    % KNOWN IC: Dt (injection), Nt (injection), all complexes at 0. Et from
    % previous characterization
    % 
    % ODEs (reduced via conservation laws)
    % 
    % 

    syms P D PD PDN PDt N m R mR mRAA mRt AA G mE mRE mRAAE mRtE E P0 D0 N0 R0 AA0 E0 
    syms kfp krp kterm krn kfn alph ktx k_fm krm kfR krR kcm ktermTL betaparam kfAA krAA kTL 

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %   MODEL RELATED DATA  %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    model.Nder=6;            % Number of derivatives
    model.Neq=18;             % Number of states 
    
    model.X=[P D PD PDN PDt N m R mR mRAA mRt AA G mE mRE mRAAE mRtE E];
    
    model.Noc=0;             % Number of controls
    
    % common species
    dP = -kfp * P * D + krp * PD + kterm * PDt ;      
    dN = - kfn * PD * N + krn * PDN ;
    dR = - kfR * m * R + krR * mR + ktermTL * mRt + kcm * (mRE + mRAAE + mRtE) ;
    dAA = -kfAA * mR * AA + krAA * mRAA + kcm * mRAAE ; 
    dE =  -k_fm * E * (m + mR + mRAA + mRt) + (krm + kcm) * (mE + mRE + mRAAE + mRtE);
    
    % DNA 1
    dD = -kfp * P * D + krp * PD + kterm * PDt;
    
    dPD = kfp * P * D - krp * PD + (krn + (alph2-1)*ktx) * PDN - kfn * PD * N; %MISTAKE!
    dPDN = - krn * PDN + kfn * PD * N - alph * ktx * PDN;
    dPDt = ktx * PDN - kterm * PDt;
    
    dm = -k_fm * m * E + krm * mE + ktx * PDN - kfR * m * R + krR * mR + ktermTL * mRt;
    
    dmR = kfR * m * R - krR * mR - kfAA * mR * AA + krAA * mRAA + (betaparam-1) * kTL * mRAA -k_fm * mR * E + krm * mRE;
    dmRAA = kfAA * mR * AA - krAA * mRAA - betaparam * kTL * mRAA - k_fm * mRAA * E + krm * mRAAE;
    dmRt = kTL * mRAA - ktermTL * mRt - k_fm * mRt * E + krm * mRtE;
    
    dG = kTL * mRAA; 
    
    dmE = k_fm * m * E - (krm + kcm) * mE;
    dmRE = k_fm * mR * E - (krm + kcm) * mRE;
    dmRAAE = k_fm * mRAA * E - (krm + kcm) * mRAAE;
    dmRtE = k_fm * mRt * E - (krm + kcm) * mRtE;
    
    model.F=[dP dD dPD dPDN dPDt dN dm dR dmR dmRAA dmRt dAA dG dmE dmRE dmRAAE dmRtE dE];
    g1=0; g2=0 ; g3 = 0; g4 = 0; g5 = 0; g6 = 0; g7=0;          % Controls
    g8=0; g9=0 ; g10 = 0; g11 = 0; g12 = 0; g13 = 0; g14=0;
    g15=0; g16=0 ; g17 = 0; g18 = 0;
    model.G=[g1, g2, g3, g4 , g5, g6 , g7, g8, g9 , g10, g11, g12, g13, g14 , g15, g16 , g17 , g18];
    h1=m;
    h2 =G;                    % Observables
    model.Nobs=2;                % Number of observables
    model.H=[h1, h2];
    model.IC=[P0 D0 0 0 0 N0 0 R0 0 0 0 AA0 0 0 0 0 0 E0];      % Initial conditions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  PARAMETERS CONSIDERED FOR IDENTIFIABILITY   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model.P=[kfp krp kterm krn kfn alph ktx k_fm krm kfR krR kcm ktermTL betaparam kfAA krAA kTL P0 D0 N0 R0 AA0 E0 ];
    model.Par=[krR ktermTL krAA kTL R0];% 
    model.Npar=5;            % Number of model parameters
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   GENERATING SERIES FUNCTION     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
