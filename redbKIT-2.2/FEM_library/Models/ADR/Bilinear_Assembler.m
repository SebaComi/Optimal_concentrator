function [X] = Bilinear_Assembler(MESH, FE_SPACE, Operator)

% Last MOD: Carlo Sinigaglia, 2 Feb 2022
%TODO add subdomain assembly for third dimension

% Operator can be "diffusion", "transport_x" , "transport_y" , "reaction"
% Tensors are defined A_{ijk} =  grad.phi_i grad.phi_j      phi_k
%                       Bx_{ijk} = phi_i      gradx.phi_j     phi_k
%                       By_{ijk} = phi_i      grady.phi_j     phi_k
%                       C_{ijk} =  phi_i      phi_j           phi_k
    
% C assembly, returns matrices in sparse vector format
index_subd = [1:MESH.numElem];
[Arows, Acols,Atens,Mtens] = bilinear_assembler_C_comp(MESH.dim,MESH.elements, FE_SPACE.numElemDof,FE_SPACE.quad_weights, ...
                                                       MESH.invjac(index_subd,:,:), MESH.jac(index_subd), FE_SPACE.phi,FE_SPACE.dphi_ref,Operator);

subs = [Arows , Acols , Atens]; %<-- Subscripts of the nonzeros.
vals = Mtens; %<-- The values of the nonzeros.
X = sptensor(subs,vals);



return