%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 TL IDENTIFIABILITY IN TXTL (v1)                         %%%
%%%                           2017                                          %%% 
%%%   Vipul Singhal, California Institute of Technology                     %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = ex5_1_TL_grid()
    model.Name='ex5_1_TL_grid';

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
    % SPECIES (with 1 or 2 appended where it is specific to DNA 1 or 2. 
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
    % BIOCHEMICAL EQUATIONS (doe DNA 1, DNA 2 similar)
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
    syms P0 N0 R0 AA0 E0 D10 D20 
    syms P N R AA E 
    syms D1 PD1 PD1N PD1t m1  m1R m1RAA m1Rt G1 m1E m1RE m1RAAE m1RtE 
    syms D2 PD2 PD2N PD2t m2  m2R m2RAA m2Rt G2 m2E m2RE m2RAAE m2RtE 
    syms kfp krp kterm krn kfn alph1 alph2 ktx k_fm krm kfR krR kcm ktermTL betaparam1 betaparam2 kfAA krAA kTL 

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %   MODEL RELATED DATA  %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    model.Nder=9;            % Number of derivatives
    model.Neq=31;             % Number of states 
    
    model.X=[P N R AA E D1 PD1 PD1N PD1t m1 m1R m1RAA m1Rt G1 m1E m1RE ...
        m1RAAE m1RtE D2 PD2 PD2N PD2t m2 m2R m2RAA m2Rt G2 m2E m2RE m2RAAE m2RtE];
    
    model.Noc=0;             % Number of controls
    
    % common species
    dP = -kfp * P * D1 + krp * PD1 + kterm * PD1t -kfp * P * D2 + krp * PD2 + kterm * PD2t;      
    dN = - kfn * PD1 * N + krn * PD1N - kfn * PD2 * N + krn * PD2N;
    dR = - kfR * m1 * R + krR * m1R + ktermTL * m1Rt + kcm * (m1RE + m1RAAE ...
        + m1RtE) - kfR * m2 * R + krR * m2R + ktermTL * m2Rt + kcm * (m2RE + m2RAAE + m2RtE);
    dAA = -kfAA * m1R * AA + krAA * m1RAA + kcm * m1RAAE -kfAA * m2R * AA + ...
        krAA * m2RAA + kcm * m2RAAE; 
    dE = -k_fm * E * (m2 + m2R + m2RAA + m2Rt) + (krm + kcm) * (m2E + m2RE + ...
        m2RAAE + m2RtE) -k_fm * E * (m1 + m1R + m1RAA + m1Rt) + (krm + kcm) * ...
        (m1E + m1RE + m1RAAE + m1RtE);
    
    % DNA 1
    dD1 = -kfp * P * D1 + krp * PD1 + kterm * PD1t;
    
    dPD1 = kfp * P * D1 - krp * PD1 + (krn + (alph1-1)*ktx) * PD1N - kfn * PD1 * N;
    dPD1N = - krn * PD1N + kfn * PD1 * N - alph1 * ktx * PD1N;
    dPD1t = ktx * PD1N - kterm * PD1t;
    
    dm1 = -k_fm * m1 * E + krm * m1E + ktx * PD1N - kfR * m1 * R + krR * m1R ...
        + ktermTL * m1Rt;
    
    dm1R = kfR * m1 * R - krR * m1R - kfAA * m1R * AA + krAA * m1RAA + ...
        (betaparam1-1) * kTL * m1RAA -k_fm * m1R * E + krm * m1RE;
    dm1RAA = kfAA * m1R * AA - krAA * m1RAA - betaparam1 * kTL * m1RAA - ...
        k_fm * m1RAA * E + krm * m1RAAE;
    dm1Rt = kTL * m1RAA - ktermTL * m1Rt - k_fm * m1Rt * E + krm * m1RtE;
    
    dG1 = kTL * m1RAA; 
    
    dm1E = k_fm * m1 * E - (krm + kcm) * m1E;
    dm1RE = k_fm * m1R * E - (krm + kcm) * m1RE;
    dm1RAAE = k_fm * m1RAA * E - (krm + kcm) * m1RAAE;
    dm1RtE = k_fm * m1Rt * E - (krm + kcm) * m1RtE;
    
    % DNA 2
    dD2 = -kfp * P * D2 + krp * PD2 + kterm * PD2t;
    
    dPD2 = kfp * P * D2 - krp * PD2 + (krn + (alph2-1)*ktx) * PD2N - kfn * PD2 * N;
    dPD2N = - krn * PD2N + kfn * PD2 * N - alph2 * ktx * PD2N;
    dPD2t = ktx * PD2N - kterm * PD2t;
    
    dm2 = -k_fm * m2 * E + krm * m2E + ktx * PD2N - kfR * m2 * R + krR * ...
        m2R + ktermTL * m2Rt;
    
    dm2R = kfR * m2 * R - krR * m2R - kfAA * m2R * AA + krAA * m2RAA + ...
        (betaparam2-1) * kTL * m2RAA -k_fm * m2R * E + krm * m2RE;
    dm2RAA = kfAA * m2R * AA - krAA * m2RAA - betaparam2 * kTL * m2RAA - ...
        k_fm * m2RAA * E + krm * m2RAAE;
    dm2Rt = kTL * m2RAA - ktermTL * m2Rt - k_fm * m2Rt * E + krm * m2RtE;
    
    dG2 = kTL * m2RAA; 
    
    dm2E = k_fm * m2 * E - (krm + kcm) * m2E;
    dm2RE = k_fm * m2R * E - (krm + kcm) * m2RE;
    dm2RAAE = k_fm * m2RAA * E - (krm + kcm) * m2RAAE;
    dm2RtE = k_fm * m2Rt * E - (krm + kcm) * m2RtE;
    
    model.F=[dP dN dR dAA dE dD1 dPD1 dPD1N dPD1t dm1 dm1R dm1RAA dm1Rt dG1 ...
        dm1E dm1RE dm1RAAE dm1RtE dD2 dPD2 dPD2N dPD2t dm2 dm2R dm2RAA ...
        dm2Rt dG2 dm2E dm2RE dm2RAAE dm2RtE];
    g1=0; g2=0 ; g3 = 0; g4 = 0; g5 = 0; g6 = 0; g7=0;          % Controls
    g8=0; g9=0 ; g10 = 0; g11 = 0; g12 = 0; g13 = 0; g14=0;
    g15=0; g16=0 ; g17 = 0; g18 = 0;
    g19=0; g20=0 ; g21 = 0 ; g22 = 0; g23 = 0; g24 = 0; g25 = 0; g26 = 0; g27=0;          % Controls
    g28=0; g29=0 ; g30 = 0; g31 = 0;
    model.G=[g1, g2, g3, g4 , g5, g6 , g7, g8, g9 , g10, g11, g12, g13,...
        g14 , g15, g16 , g17 , g18, g19, g20, g21, g22, g23, g24 , g25, g26 , g27, g28, g29 , g30, g31 ];
    h1=m1;
    h2 =G1;
    h3 = m2; 
    h4 = G2; % Observables
    model.Nobs=4;                % Number of observables
    model.H=[h1, h2 h3 h4];
    model.IC=[P0 N0 R0 AA0 E0 D10 0 0 0 0 0 0 0 0 0 0 ...
        0 0 D20 0 0 0 0 0 0 0 0 0 0 0 0];      % Initial conditions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  PARAMETERS CONSIDERED FOR IDENTIFIABILITY   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model.P=[kfp krp kterm krn kfn ktx k_fm krm kfR krR kcm ktermTL ...
        betaparam2 betaparam1 kfAA krAA kTL P0 D10 D20 N0 R0 AA0 E0 alph1 alph2];
    model.Par=[krR ktermTL krAA kTL R0];% 
    model.Npar=5;            % Number of model parameters
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   GENERATING SERIES FUNCTION     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
