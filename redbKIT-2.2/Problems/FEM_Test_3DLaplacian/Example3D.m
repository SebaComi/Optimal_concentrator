%EXAMPLE3D shows how to solve a 3D problem with the FEM library. Also performs
%convergence analysis.
%
% GENERATE the meshes before running (see readme in gmsh subfolder)

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

clc
clear all

[~,~,~] = mkdir('Figures');

fem      =  'P2';
dim      =  3;

h = [2/4 2/9];% 2/19 2/29];

%% solve for different level of refinement
for i = 1 : length(h)
    
    %% load P1 mesh
    [vertices, boundaries, elements] = msh_to_Mmesh(strcat('gmsh/Cubeh',num2str(i)), dim);
%     return
    %% Solve
    [U, FE_SPACE, MESH, DATA, errorL2(i), errorH1(i)]  = Elliptic_Solver(dim, elements, vertices, boundaries, fem, 'Dirichlet_data');
%     X = MESH.vertices(1,:)';
%     Y = MESH.vertices(2,:)';
%     Z = MESH.vertices(3,:)';
%     [Xq,Yq,Zq] = meshgrid(linspace(min(X),max(X),10),linspace(min(Y),max(Y),10),linspace(min(Z),max(Z),10));
%     Uq = interp3(X,Y,Z,U,Xq,Yq,Zq);
%     slice(Xq,Yq,Zq,Uq,0,0,0)
%     pdeplot3D(MESH.nodes,MESH.elements,'ColorMapData',U,'Mesh','on','FaceAlpha',1)
% return

    %% export solution
    ADR_export_solution(dim, U(1:MESH.numVertices), MESH.vertices, MESH.elements, ['Figures/SOL_',num2str(i)]);
    
end

%% Plot convergence order
r_FE = str2double(fem(2));

figure
loglog(h,errorL2,'-ob','LineWidth',2)
hold on
loglog(h,errorH1,'-or','LineWidth',2)
loglog(h, errorH1(1)/2/h(1)*h.^r_FE, 'r--', h, errorL2(1)/2/(h(1)^2)*h.^(r_FE+1), 'b--')
legend(['errL2',fem],['errH1',fem], ['h^',num2str(r_FE)],  ['h^',num2str(r_FE+1)])
title(['Convergence ',fem,' Finite elements'])
xlim([min(h) max(h)])
grid on

