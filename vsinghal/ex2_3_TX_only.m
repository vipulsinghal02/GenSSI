%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 TX IDENTIFIABILITY IN TXTL (v3: no conservation law)    %%%
%%%                           2017                                          %%% 
%%%   Vipul Singhal, California Institute of Technology                     %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = ex2_3_TX_only()
    model.Name='ex2_3_TX_only';
    % difference from version 1: Removed the conservation laws
    % 
    % overall:
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
    % m = mRNA
    % E = RNAse enzyme
    % uN = used Nucleotides
    % P = RNA Pol
    % D = free DNA
    % N = unused nucleotides
    % PDt = tx termination complex
    % PD = P and D complex
    % PDN = P, D and N complex
    % 
    % BIOCHEMICAL EQUATIONS
    % mrna degradation (previously characterized)
    % m + E <-> C -> E + 600 * uN (k_fm, krm, kcm, all known)
    % 
    % transcription
    % P + D <-> PD (krp, kfp; kfp known)
    % PD + N <-> PDN (krn, kfn; kfn known)
    % PDN -> PD (ktx*(alph - 1)/alph)
    % PDN -> PDt + m (ktx/alph)
    % PDt -> P + D (kterm)
    % alph base pairs, known
    % 
    % just rescaling ktx:
    % m + E <-> C -> E (k_fm, krm, kcm, all known)
    % E + C = Et (Et known)
    % P + D <-> PD (krp, kfp; kfp known)
    % PD + N <-> PDN (krn, kfn; kfn known)
    % PDN -> PD ((alph - 1)*ktx)
    % PDN -> PDt + m (ktx)
    % PDt -> P + D (kterm)
    % n = 600 base pairs, known
    % 
    % the ODE model:
    % dm = -k_fm * m * E + krm * C + ktx * PDN
    % dE = -k_fm * m * E + krm * C - kcm * C
    % dC = k_fm * m * E - krm * C - kcm * C
    % dP = -kfp * P * D + krp * PD + kterm * PDt
    % dD = -kfp * P * D + krp * PD + kterm * PDt
    % dPD = kfp * P * D - krp * PD + krn * PDN - kfn * PD * N
    % dPDN = - krn * PDN + kfn * PD * N - alph * ktx * PDN
    % dPDt = ktx * PDN - kterm * PDt
    % dN = - kfn * PD * N + krn * PDN
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
    % KNOWN IC: Dt (injection), Nt (injection), all complexes at 0. Et from
    % previous characterization
    % 
    % 
    % 

    syms m C E0 PDN PD PDt N D P0 D0 alph N0 E P 
    syms k_fm krm ktx kcm kfp krp kterm krn kfn 

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %   MODEL RELATED DATA  %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    model.Nder=7;            % Number of derivatives
    model.Neq=9;             % Number of states 
    model.X=[m E C P D PD PDN PDt N];
    model.Noc=0;             % Number of controls
    dm = -k_fm * m * E + krm * C + ktx * PDN;
    dE = -k_fm * m * E + krm * C + kcm * C;
    dC = k_fm * m * E - krm * C - kcm * C;
    dP = -kfp * P * D + krp * PD + kterm * PDt;
    dD = -kfp * P * D + krp * PD + kterm * PDt;
    dPD = kfp * P * D - krp * PD + krn * PDN - kfn * PD * N;
    dPDN = - krn * PDN + kfn * PD * N - alph * ktx * PDN;
    dPDt = ktx * PDN - kterm * PDt;
    dN = - kfn * PD * N + krn * PDN;
    

    model.F=[dm dE dC dP dD dPD dPDN dPDt dN];
    g1=0; g2=0 ; g3 = 0; g4 = 0; g5 = 0; g6 = 0; g7=0;  g8 = 0; g9=0;         % Controls
    model.G=[g1, g2, g3, g4 , g5, g6 , g7, g8 , g9];
    h1=m;                    % Observables
    model.Nobs=1;                % Number of observables
    model.H=[h1];
    model.IC=[0 E0 0 P0 D0 0 0 0 N0];      % Initial conditions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  PARAMETERS CONSIDERED FOR IDENTIFIABILITY   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model.P=[k_fm krm ktx kcm kfp krp kterm krn kfn E0 P0 D0 N0 alph];
    model.Par=[ktx krp kterm krn P0];% 
    model.Npar=5;            % Number of model parameters
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   GENERATING SERIES FUNCTION     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
