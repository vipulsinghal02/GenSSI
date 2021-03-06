function multiExp = M1_1_eDef4()
    % mRNA model M1 (trivial), 1 observable, original model, 2 controls
    % 4 experiments: 
    %   original configuration
    %   change in transfection (m0) via IC
    %   change in translation via control u1 = uInh
    %   change in GFP degradation via control u2 = uDeg
    multiExp.Name='M1_1_eDef4';
    syms d b kTL m0Exp1
    multiExp.Nexp=4;
    multiExp.U = [1,1;...
                  1,1;...
                  .5,1;...
                  1,.75];
    multiExp.IC = [m0Exp1,0,...
                   0.5*m0Exp1,0,...
                   m0Exp1,0,...
                   m0Exp1,0];
    multiExp.P=[d,b,kTL,m0Exp1];
    multiExp.Par=[d,b,kTL,m0Exp1];
end