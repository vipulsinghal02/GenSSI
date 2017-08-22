%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Experiment 4: TX IDENTIFIABILITY IN TXTL (v1: TX GRID)                %%%
%%%                           2017                                          %%% 
%%%   Vipul Singhal, California Institute of Technology                     %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = ex4_1_TX_grid()
    model.Name='ex4_1_TX_grid';
    % tx grid, two DNAs, measuring two different mrnas
    % 
    % Here we fix the RNA degradation numbers (assuming they are known from
    % a previous experiment, and check the identifiability of the
    % transcription parameters. 
    % 
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     DECLARE SYMBOLIC VARIABLES:                          %
    %               - state variables                          %
    %               - parameters of the model                  %
    %               - initial state, if not known              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPECIES 
    % m1, m2 = mRNA 1 and 2
    % E = RNAse enzyme
    % P = RNA Pol
    % D1 = free DNA 1
    % D2 = free DNA 2 
    % N = unused nucleotides
    % PD1t, PD2t = tx termination complexes
    % PDi = P and Di complex
    % PDiN = P, Di and N complex
    % 
    % BIOCHEMICAL EQUATIONS
    % mrna degradation (previously characterized)
    % mi + E <-> Ci -> E(k_fm, krm, kcm, all known)
    % 
    % transcription
    % P + D1 <-> PD1 (krp, kfp; kfp known)
    % PD1 + N <-> PD1N (krn, kfn; kfn known)
    % PD1N -> PD1 (ktx*(alph1 - 1))
    % PD1N -> PD1t + m1 (ktx)
    % PD1t -> P + D1 (kterm)
    % P + D2 <-> PD2 (krp, kfp; kfp known)
    % PD2 + N <-> PD2N (krn, kfn; kfn known)
    % PD2N -> PD2 (ktx*(alph2 - 1))
    % PD2N -> PD2t + m2 (ktx)
    % PD2t -> P + D2 (kterm)
    % alph base pairs, known
    % 
    % 
    % the ODE model:
    % dm1 = -k_fm * m1 * E + krm * C1 + ktx * PD1N
    % dC1 = k_fm * m1 * E - krm * C1 - kcm * C1
    % dD1 = -kfp * P * D1 + krp * PD1 + kterm * PD1t
    % dPD1 = kfp * P * D1 - krp * PD1 + krn * PD1N - kfn * PD1 * N
    % dPD1N = - krn * PD1N + kfn * PD1 * N - alph1 * ktx * PD1N
    % dPD1t = ktx * PD1N - kterm * PD1t
    %
    % dm2 = -k_fm * m2 * E + krm * C2 + ktx * PD2N
    % dC2 = k_fm * m2 * E - krm * C2 - kcm * C2
    % dD2 = -kfp * P * D2 + krp * PD2 + kterm * PD2t
    % dPD2 = kfp * P * D2 - krp * PD2 + krn * PD2N - kfn * PD2 * N
    % dPD2N = - krn * PD2N + kfn * PD2 * N - alph2 * ktx * PD2N
    % dPD2t = ktx * PD2N - kterm * PD2t
    
    % dP = -kfp * P * D1 + krp * PD1 + kterm * PD1t -kfp * P * D2 + krp * PD2 + kterm * PD2t
    % dE = -k_fm * m1 * E + krm * C1 - kcm * C1 -k_fm * m2 * E + krm * C2 - kcm * C2
    % dN = - kfn * PD1 * N + krn * PD1N - kfn * PD2 * N + krn * PD2N
    %
    %
    %
    % UNKNOWN PARAMETERS 
    % krp, krn, ktx, kterm 
    % 
    % KNOWN PARAMETERS
    % k_fm krm kcm kfp kfn 
    %
    % UNKNOWN IC: Pt = P(0)
    % 
    % KNOWN IC: D10, D20 (injection), N0 (injection), all complexes at 0. E0 from
    % previous characterization
    % 

    syms m1 m2 C1 C2 E0 PD1N PD1 PD1t D1 D10 alph1 PD2N PD2 PD2t D2 D20 alph2 N P0 N0 E P 
    syms k_fm krm ktx kcm kfp krp kterm krn kfn 

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %   MODEL RELATED DATA  %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    model.Nder=8;            % Number of derivatives
    model.Neq=15;             % Number of states 
    model.X=[m1 C1 D1 PD1 PD1N PD1t m2 C2 D2 PD2 PD2N PD2t P E N];
    model.Noc=0;             % Number of controls
    dm1 = -k_fm * m1 * E + krm * C1 + ktx * PD1N;
    dC1 = k_fm * m1 * E - krm * C1 - kcm * C1;
    dD1 = -kfp * P * D1 + krp * PD1 + kterm * PD1t;
    dPD1 = kfp * P * D1 - krp * PD1 + krn * PD1N - kfn * PD1 * N;
    dPD1N = - krn * PD1N + kfn * PD1 * N - alph1 * ktx * PD1N;
    dPD1t = ktx * PD1N - kterm * PD1t;
    
    dm2 = -k_fm * m2 * E + krm * C2 + ktx * PD2N;
    dC2 = k_fm * m2 * E - krm * C2 - kcm * C2;
    dD2 = -kfp * P * D2 + krp * PD2 + kterm * PD2t;
    dPD2 = kfp * P * D2 - krp * PD2 + krn * PD2N - kfn * PD2 * N;
    dPD2N = - krn * PD2N + kfn * PD2 * N - alph2 * ktx * PD2N;
    dPD2t = ktx * PD2N - kterm * PD2t;
    
    dP = -kfp * P * D1 + krp * PD1 + kterm * PD1t -kfp * P * D2 + krp * PD2 + kterm * PD2t;
    dE = -k_fm * m1 * E + krm * C1 - kcm * C1 -k_fm * m2 * E + krm * C2 - kcm * C2;
    dN = - kfn * PD1 * N + krn * PD1N - kfn * PD2 * N + krn * PD2N;
    

    model.F=[dm1 dC1 dD1 dPD1 dPD1N dPD1t dm2 dC2 dD2 dPD2 dPD2N dPD2t dP dE dN];
    g1=0; g2=0 ; g3 = 0; g4 = 0; g5 = 0; g6 = 0; g7=0;  g8 = 0; g9=0; g10 = 0; g11 = 0; g12 = 0; g13=0;  g14 = 0; g15=0;         % Controls
    model.G=[g1, g2, g3, g4 , g5, g6 , g7, g8 , g9, g10, g11, g12, g13, g14, g15];
    h1=m1; h2 = m2 ;                  % Observables
    model.Nobs=2;                % Number of observables
    model.H=[h1, h2];
    model.IC=[0 0 D10 0 0 0 0 0 D20 0 0 0 P0 E0 N0];      % Initial conditions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  PARAMETERS CONSIDERED FOR IDENTIFIABILITY   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model.P=[k_fm krm ktx kcm kfp krp kterm krn kfn E0 P0 D10 D20 N0 alph1 alph2];
    model.Par=[ktx krp kterm krn P0];% 
    model.Npar=5;            % Number of model parameters
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   GENERATING SERIES FUNCTION     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
