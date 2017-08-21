%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 THE CONSTITUTIVE EXPRESSION IN TXTL                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = genSSI_geneexpr()
    model.Name='geneexpr';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     DECLARE SYMBOLIC VARIABLES:                          %
    %               - state variables                          %
    %               - parameters of the model                  %
    %               -initial state, if not known               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %   STATE VARIABLES     %
    %%%%%%%%%%%%%%%%%%%%%%%%%    

RNAP
[protein sigma70] 
[protein sigma28] = 20
Ribo = 50
RNAP70 = 0
RecBCD = 0
[protein gamS] = 0
[RecBCD:gamS] = 0
RNase = 600
AGTP = 3.18e+06
CUTP = 1.908e+06
AA = 3.18e+07
[protein deGFP] = 0
[protein deGFP*] = 0
[RNA rbs--deGFP] = 0
[Ribo:RNA rbs--deGFP] = 0
[DNA p70--rbs--deGFP] = 0
[RNAP70:DNA p70--rbs--deGFP] = 0
[CUTP:AGTP:RNAP70:DNA p70--rbs--deGFP] = 0
[term_RNAP70:DNA p70--rbs--deGFP] = 0
[AA:AGTP:Ribo:RNA rbs--deGFP] = 0
[AGTP:RNAP70:DNA p70--rbs--deGFP] = 0
[CUTP:RNAP70:DNA p70--rbs--deGFP] = 0
[RNA rbs--deGFP:RNase] = 0
[AA:AGTP:Ribo:RNA rbs--deGFP:RNase] = 0
[Ribo:RNA rbs--deGFP:RNase] = 0
RNAP28 = 0
AGTP_UNUSE = 0
AGTPdeg_F = 0
AGTPdegrate = 0.000147
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %       PARAMETERS      %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    sym   TXTL_PROT_deGFP_MATURATION
    sym TXTL_RNAPBOUND_TERMINATION_RATE
    sym TXTL_RNAdeg_catalysis
    sym TXTL_RBS_Ribo_deGFP_Kd
    sym TXTL_RBS_Ribo_deGFP_F
    sym TXTL_RBS_Ribo_deGFP_R
    sym TXTL_P70_RNAPbound_deGFP_Kd
    sym TXTL_P70_RNAPbound_deGFP_F
    sym TXTL_P70_RNAPbound_deGFP_R
    sym TXTL_RNAP_S70_Kd
    sym TXTL_RNAP_S70_F
    sym TXTL_RNAP_S70_R
    sym TXTL_NTP_RNAP_1_Kd
    sym TXTL_NTP_RNAP_1_F
    sym TXTL_NTP_RNAP_1_R
    sym TXTL_NTP_RNAP_2_Kd
    sym TXTL_NTP_RNAP_2_F
    sym TXTL_NTP_RNAP_2_R
    sym TXTL_AA_Kd
    sym TXTL_AA_F
    sym TXTL_AA_R
    sym TXTL_RNAdeg_Kd
    sym TXTL_RNAdeg_F
    sym TXTL_RNAdeg_R
    sym TXTL_transcription_rate1
    sym k_elon
    sym TXTL_k_con_TX
    sym TXTL_TL_rate
    sym k_elon_prot
    sym TXTL_k_con_TL
    sym GamS_RecBCD_f
    sym GamS_RecBCD_r
    sym RNAPdeg_F
    sym RNAP70deg_F
    sym RNAP28deg_F
    sym contents
    sym TotalRNA
    syms x1 x2 x3 p1 p2 p3 p4 p5 p6 p7 p8

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %   MODEL RELATED DATA  %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    model.Nder=8;            % Number of derivatives

    model.Neq=3;             % Number of states 
    model.X=[x1 x2 x3];

    model.Npar=8;            % Number of model parameters

    model.Noc=0;             % Number of controls

    model.Nobs=1;                % Number of observables

    A1 = -p4*x1+p1/(p2+x3^p3);   % Equations of the model 
    A2 = p5*x1-p6*x2;
    A3 = p7*x2-p8*x3;
    model.F=[A1 A2 A3];

    g1=0;g2=0;g3=0;              % Controls
    model.G=[g1 g2 g3];

    h1=x1;h2=x2;h3=x3;           % Observables
    model.H=[h1 h2 h3];

    model.IC=[0.3 0.9 1.3];      % Initial conditions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  PARAMETERS CONSIDERED FOR IDENTIFIABILITY   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    model.Par=[p1 p2 p4 p5 p6 p7 p8];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   GENERATING SERIES FUNCTION     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
