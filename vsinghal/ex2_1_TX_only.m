%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 TX IDENTIFIABILITY IN TXTL                              %%%
%%%                           2017                                          %%% 
%%%   Vipul Singhal, California Institute of Technology                     %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = ex2_1_TX_only()
    model.Name='ex2_1_TX_only';
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
    % E + C = Et (Et known)
    % 
    % transcription
    % P + D <-> PD (krp, kfp; kfp known)
    % PD + N <-> PDN (krn, kfn; kfn known)
    % PDN -> PD (ktx*(alph - 1)/alph)
    % PDN -> PDt + m (ktx/alph)
    % PDt -> P + D (kterm)
    % alph base pairs, known
    % P + PD + PDN + PDt = Pt (unknown, parameter of interest)
    % D + PD + PDN + PDt = Dt (known, initial condition control)
    % N + PDN + uN + 600*(m + C) = Nt (known, initial condition control)
    % 
    % We can rewrite this by removing uN and the corresponding conseravion
    % law, and just rescaling ktx:
    % m + E <-> C -> E (k_fm, krm, kcm, all known)
    % E + C = Et (Et known)
    % P + D <-> PD (krp, kfp; kfp known)
    % PD + N <-> PDN (krn, kfn; kfn known)
    % PDN -> PD ((alph - 1)*ktx)
    % PDN -> PDt + m (ktx)
    % PDt -> P + D (kterm)
    % n = 600 base pairs, known
    % P + PD + PDN + PDt = Pt (unknown, parameter of interest)
    % D + PD + PDN + PDt = Dt (known, initial condition control)
    % 
    % And reducing the model with conservation laws we get the ODE model:
    % dm = -k_fm * m * (Et - C) + krm * C + ktx * PDN
    % dC = k_fm * m * (Et - C) - krm * C - kcm * C
    % dP = -kfp * P * (P + Dt - Pt) + krp * PD + kterm * PDt
    % dPD = kfp * P * (P + Dt - Pt) - krp * PD + krn * PDN - kfn * PD * N
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
    % ODEs (reduced via conservation laws)
    % 
    % 

    syms m C Et PDN PD PDt N D Pt Dt alph N0
    syms k_fm krm ktx kcm kfp krp kterm krn kfn 

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %   MODEL RELATED DATA  %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    model.Nder=4;            % Number of derivatives
    model.Neq=7;             % Number of states 
    model.X=[m C D PD PDN PDt N];
    model.Noc=0;             % Number of controls
    dm = -k_fm * m * (Et - C) + krm * C + ktx * PDN;
    dC = k_fm * m * (Et - C) - krm * C - kcm * C;
    dD = -kfp * (Pt - Dt +D) * D + krp * PD + kterm * PDt;
    dPD = kfp * (Pt - Dt +D) * D - krp * PD + krn * PDN - kfn * PD * N;
    dPDN = - krn * PDN + kfn * PD * N - alph * ktx * PDN;
    dPDt = ktx * PDN - kterm * PDt;
    dN = - kfn * PD * N + krn * PDN;
    model.F=[dm dC dD dPD dPDN dPDt dN];
    g1=0; g2=0 ; g3 = 0; g4 = 0; g5 = 0; g6 = 0; g7=0;          % Controls
    model.G=[g1, g2, g3, g4 , g5, g6 , g7];
    h1=m;                    % Observables
    model.Nobs=1;                % Number of observables
    model.H=[h1];
    model.IC=[0 0 Dt 0 0 0 N0];      % Initial conditions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  PARAMETERS CONSIDERED FOR IDENTIFIABILITY   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model.P=[k_fm krm ktx kcm kfp krp kterm krn kfn Et Pt Dt N0 alph];
    model.Par=[ktx krp kterm krn Pt Dt N0];% 
    model.Npar=7;            % Number of model parameters
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   GENERATING SERIES FUNCTION     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
