function multiExp = multiExp_RNAdeg()
% /Users/vipulsinghal/Dropbox/Documents/vipul_repo/Identifiability/Code_Identifiability/multiExp_RNAdeg.m
    %   original configuration
    %   change in transfection (m0) via IC, 2 additional concentrations
    multiExp.Name='multiExp_RNAdeg';
    syms kr kc Et m0Exp1
    multiExp.Nexp=3;
    multiExp.U = [0;...
                  0;...
                  0];
    multiExp.IC = [m0Exp1,0,...
                   0.5*m0Exp1,0,...
                   2*m0Exp1,0];
    multiExp.P=[kr kc Et m0Exp1];
    multiExp.Par=[kr kc Et m0Exp1];
end