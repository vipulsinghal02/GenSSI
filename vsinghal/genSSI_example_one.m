%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 THE GOODWIN OSCILLATOR                                   %%%
%%%  Bibliography: Goodwin, B.C., Oscillatory behavior in enzymatic control  %%% 
%%%                processes, Adv. Enzyme Regul. (3), 425–428, 1965.        %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = genSSI_example_one()
    model.Name='genSSI_example_one';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     DECLARE SYMBOLIC VARIABLES:                          %
    %               - state variables                          %
    %               - parameters of the model                  %
    %               - initial state, if not known              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    syms mrna prot am ap dm dp km kp dna

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %   MODEL RELATED DATA  %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    model.Nder=8;            % Number of derivatives

    model.Neq=3;             % Number of states 
    model.X=[dna mrna prot];

    model.Npar=6;            % Number of model parameters

    model.Noc=0;             % Number of controls

    model.Nobs=3;                % Number of observables
    A1 = 0;
    A2 = am*dna/(km+dna)-dm*mrna;   % Equations of the model 
    A3 = ap*mrna/(kp+mrna)-dp*prot;
    model.F=[A1 A2 A3];

    g1=0;g2=0;g3=0;            % Controls
    model.G=[g1 g2 g3];

    h1=dna; h2=mrna; h3=prot;           % Observables
    model.H=[h1 h2 h3];

    model.IC=[5 0 0];      % Initial conditions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  PARAMETERS CONSIDERED FOR IDENTIFIABILITY   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    model.Par=[  am ap dm dp ];% km   kp

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   GENERATING SERIES FUNCTION     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
