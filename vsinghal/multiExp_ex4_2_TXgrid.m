function multiExp = multiExp_ex4_2_TXgrid()

    multiExp.Name='multiExp_ex4_2_TXgrid';
    syms k_fm krm ktx kcm kfp krp kterm krn kfn E0 P0 D10ex4 D20ex4 N0ex4 alph1 alph2
    
    multiExp.Nexp=27;
    
    multiExp.U = [0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0;...
                  0];
              
%     multiExp.IC = [0,0, 0.5*Dt0Exp2, 0 0 0,  0.5*N0exp2,...
%                    0,0, 0.5*Dt0Exp2, 0 0 0,  N0exp2,...
%                    0,0, 0.5*Dt0Exp2, 0 0 0,  2*N0exp2,...
%                    0,0, Dt0Exp2, 0 0 0,  0.5*N0exp2,...
%                    0,0, Dt0Exp2, 0 0 0,  N0exp2,...
%                    0,0, Dt0Exp2, 0 0 0,  2*N0exp2,...
%                    0,0, 2*Dt0Exp2, 0 0 0,  0.5*N0exp2,...
%                    0,0, 2*Dt0Exp2, 0 0 0,  N0exp2,...
%                    0,0, 2*Dt0Exp2, 0 0 0,  2*N0exp2];
               
    multiExp.IC = [0 0 0.5*D10ex4 0 0 0 0 0 0.5*D20ex4 0 0 0 P0 E0 0.5*N0ex4,...
        0 0 0.5*D10ex4 0 0 0 0 0 D20ex4 0 0 0 P0 E0 0.5*N0ex4,...
        0 0 0.5*D10ex4 0 0 0 0 0 2*D20ex4 0 0 0 P0 E0 0.5*N0ex4,...
        0 0 D10ex4 0 0 0 0 0 0.5*D20ex4 0 0 0 P0 E0 0.5*N0ex4,...
        0 0 D10ex4 0 0 0 0 0 D20ex4 0 0 0 P0 E0 0.5*N0ex4,...
        0 0 D10ex4 0 0 0 0 0 2*D20ex4 0 0 0 P0 E0 0.5*N0ex4,...
        0 0 2*D10ex4 0 0 0 0 0 0.5*D20ex4 0 0 0 P0 E0 0.5*N0ex4,...
        0 0 2*D10ex4 0 0 0 0 0 D20ex4 0 0 0 P0 E0 0.5*N0ex4,...
        0 0 2*D10ex4 0 0 0 0 0 2*D20ex4 0 0 0 P0 E0 0.5*N0ex4,...
        0 0 0.5*D10ex4 0 0 0 0 0 0.5*D20ex4 0 0 0 P0 E0 N0ex4,...
        0 0 0.5*D10ex4 0 0 0 0 0 D20ex4 0 0 0 P0 E0 N0ex4,...
        0 0 0.5*D10ex4 0 0 0 0 0 2*D20ex4 0 0 0 P0 E0 N0ex4,...
        0 0 D10ex4 0 0 0 0 0 0.5*D20ex4 0 0 0 P0 E0 N0ex4,...
        0 0 D10ex4 0 0 0 0 0 D20ex4 0 0 0 P0 E0 N0ex4,...
        0 0 D10ex4 0 0 0 0 0 2*D20ex4 0 0 0 P0 E0 N0ex4,...
        0 0 2*D10ex4 0 0 0 0 0 0.5*D20ex4 0 0 0 P0 E0 N0ex4,...
        0 0 2*D10ex4 0 0 0 0 0 D20ex4 0 0 0 P0 E0 N0ex4,...
        0 0 2*D10ex4 0 0 0 0 0 2*D20ex4 0 0 0 P0 E0 N0ex4,...
        0 0 0.5*D10ex4 0 0 0 0 0 0.5*D20ex4 0 0 0 P0 E0 2*N0ex4,...
        0 0 0.5*D10ex4 0 0 0 0 0 D20ex4 0 0 0 P0 E0 2*N0ex4,...
        0 0 0.5*D10ex4 0 0 0 0 0 2*D20ex4 0 0 0 P0 E0 2*N0ex4,...
        0 0 D10ex4 0 0 0 0 0 0.5*D20ex4 0 0 0 P0 E0 2*N0ex4,...
        0 0 D10ex4 0 0 0 0 0 D20ex4 0 0 0 P0 E0 2*N0ex4,...
        0 0 D10ex4 0 0 0 0 0 2*D20ex4 0 0 0 P0 E0 2*N0ex4,...
        0 0 2*D10ex4 0 0 0 0 0 0.5*D20ex4 0 0 0 P0 E0 2*N0ex4,...
        0 0 2*D10ex4 0 0 0 0 0 D20ex4 0 0 0 P0 E0 2*N0ex4,...
        0 0 2*D10ex4 0 0 0 0 0 2*D20ex4 0 0 0 P0 E0 2*N0ex4];
               
    multiExp.P=[k_fm krm ktx kcm kfp krp kterm krn kfn E0 P0 D10ex4 D20ex4 N0ex4 alph1 alph2];
    multiExp.Par=[ktx krp kterm krn P0];
end