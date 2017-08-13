function model = JAK_STAT()
    % JAK_STAT provides the GenSSI implementation of the JAK-STAT signaling
    % pathway model introduced by
    % 
    %    Raue et al. (2009). Structural and practical identifiability 
    %    analysis of partially observed dynamical models by exploiting 
    %    the profile likelihood. Bioinformatics 25(15): 1923-1929.

    % Model name
    model.Name = 'JAK_STAT';

    % Symbolic variables
    syms STAT pSTAT pSTAT_pSTAT npSTAT_npSTAT nSTAT1 nSTAT2 nSTAT3 nSTAT4 nSTAT5
    syms p1 p2 p3 p4 init_STAT Omega_cyt Omega_nuc

    % State variables
    model.X = [STAT,...
               pSTAT,...
               pSTAT_pSTAT,...
               npSTAT_npSTAT,...
               nSTAT1,...
               nSTAT2,...
               nSTAT3,...
               nSTAT4,...
               nSTAT5];

    % Right-hand side of differential equation
    model.F = [Omega_nuc/Omega_cyt*p4*nSTAT5,...
               -2*p2*pSTAT^2,...
               p2*pSTAT^2-p3*pSTAT_pSTAT,...
               Omega_cyt/Omega_nuc*p3*pSTAT_pSTAT-p4*npSTAT_npSTAT,...
               p4*(2*npSTAT_npSTAT-nSTAT1),...
               p4*(nSTAT1-nSTAT2),...
               p4*(nSTAT2-nSTAT3),...
               p4*(nSTAT3-nSTAT4),...
               p4*(nSTAT4-nSTAT5)];

    % Initial conditions
    model.IC = [init_STAT,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0];
    
    % Controls
    model.G = [-p1*STAT,...
                p1*STAT,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0,...
                0];

    % Observables
    model.H = [(pSTAT+2*pSTAT_pSTAT)/init_STAT,...
               (STAT+pSTAT+2*pSTAT_pSTAT)/init_STAT];

    % Parameters (for which structural identifiability is analyzed)
    model.Par = [p1,p2,p3,p4,init_STAT,Omega_cyt,Omega_nuc];
    
    % Number of Lie derivatives
    model.Nder = 7;
end

