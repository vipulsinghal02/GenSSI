%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 RNA DEGRADATION IN TXTL                                 %%%
%%%                           2017                                          %%% 
%%%   Vipul Singhal, California Institute of Technology                     %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = ex1_genSSI_RNAdeg()
    model.Name='TXTL_RNAdeg';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     DECLARE SYMBOLIC VARIABLES:                          %
    %               - state variables                          %
    %               - parameters of the model                  %
    %               - initial state, if not known              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPECIES 
    % m = mRNA
    % E = RNAse enzyme
    % 
    % BIOCHEMICAL EQUATIONS
    % m + E <-> C -> E
    % E + C = Et
    % 
    % Known IC:
    % m = dosing at known conc: 50 100 500 1000 2000 5000
    % C(0) = 0
    %
    % UNKNOWN PARAMETERS 
    % kf kr kc 
    % 
    % UNKNOWN IC: Et = E(0)
    % 
    % ODEs
    % dm = -kf m E + kr C
    % dC = +kf m E - (kr+kc) C 
    % dE = -dC

    syms m C m0
    syms kr kc Et

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %   MODEL RELATED DATA  %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    model.Nder=4;            % Number of derivatives
    model.Neq=2;             % Number of states 
    model.X=[m C];
    model.Noc=0;             % Number of controls
    dmdt = -100 * m * (Et - C) + kr * C;
    dCdt =  100 * m * (Et - C) - (kr+kc) * C;   % Equations of the model 
    model.F=[dmdt dCdt];
    g1=0; g2=0 ;           % Controls
    model.G=[g1, g2];
    h1=m;                    % Observables
    model.Nobs=1;                % Number of observables
    model.H=[h1];
    model.IC=[m0 0];      % Initial conditions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  PARAMETERS CONSIDERED FOR IDENTIFIABILITY   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model.P=[kr kc Et m0];
    model.Par=[kr kc Et m0];% 
    model.Npar=4;            % Number of model parameters
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   GENERATING SERIES FUNCTION     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
