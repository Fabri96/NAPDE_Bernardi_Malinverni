function [M,A,Adv_glob,Ms_glob,f,tauu]=C_matrix2D_theory(Dati,femregion)
%% [M,A,f] = C_matrix2D(Dati,femregion)
%==========================================================================
% Assembly of the mass matrix M, stiffness matrix A and rhs f
%==========================================================================
%    called in C_main2D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%
%    OUTPUT:
%          M           : (sparse(ndof,ndof) real) mass matrix
%          A           : (sparse(ndof,ndof) real) stiffnes matrix
%          f           : (sparse(ndof,1) real) rhs vector


addpath FESpace
addpath Assembly

% fprintf('============================================================\n')
% fprintf('Assembling matrices and right hand side ... \n');
% fprintf('============================================================\n')


% connectivity infos
ndof         = femregion.ndof; % degrees of freedom
nln          = femregion.nln;  % local degrees of freedom
ne           = femregion.ne;   % number of elements
connectivity = femregion.connectivity; % connectivity matrix


% shape functions
[basis] = C_shape_basis(Dati);

% quadrature nodes and weights for integrals
[nodes_2D, w_2D] = C_quadrature(Dati);

% evaluation of shape bases 
[dphiq,Grad] = C_evalshape(basis,nodes_2D);


% Assembly begin ...
M = sparse(ndof,ndof);  % Global Mass matrix
A = sparse(ndof,ndof);  % Global Stiffness matrix
f = sparse(ndof,1);     % Global Load vector
Adv_glob = sparse(ndof,ndof);
Ms_glob = sparse(ndof,ndof);
% Beta Functions
b_x = inline(Dati.beta_x);
b_y = inline(Dati.beta_y);
%div_b = inline(Dati.div_beta);


for ie = 1 : ne
     
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie);
  
    
    [BJ, pphys_2D] = C_get_Jacobian(femregion.coord(iglo,:), nodes_2D, Dati.MeshType);
    % BJ        = Jacobian of the elemental map 
    % pphys_2D = vertex coordinates in the physical domain 
   
    %=============================================================%
    % STIFFNESS MATRIX
    %=============================================================%
    
    % Local stiffness matrix 
    [A_loc] = C_lap_loc(Grad,w_2D,nln,BJ);
    [M_loc] = C_mass_loc(dphiq,w_2D,nln,BJ);
    [Adv_loc]=C_adv_loc(Grad,dphiq,b_x,b_y,w_2D,nln,BJ,pphys_2D);
    
    % Assembly phase for stiffness matrix
    A(iglo,iglo) = A(iglo,iglo) + Dati.mu*A_loc + Dati.sigma*M_loc;
    Adv_glob(iglo,iglo)= Adv_glob(iglo,iglo) + Adv_loc ;
    % Assembly phase for mass matrix
    M(iglo,iglo) = M(iglo,iglo) + M_loc;
    
    %==============================================
    % FORCING TERM --RHS
    %==============================================

    % Local load vector
    [load] = C_loc_rhs2D(Dati.force,dphiq,BJ,w_2D,pphys_2D,nln,Dati.mu);   

    % Assembly phase for the load vector
    f(iglo) = f(iglo) + load;

        switch Dati.name_method
        case{'SUPG'}
            S_loc = zeros(nln,nln);


            for i = 1:nln
                for j = 1:nln
                    for q = 1:length(w_2D)
                        B = BJ(:,:,q); % 2x2
                        x = pphys_2D(q,1);                         % x quadrature node
                        y = pphys_2D(q,2);
                        grad_phi_i = Grad(q,:,i)'; % 2x1
                        grad_phi_j = Grad(q,:,j)'; % 2x1
                        beta = [b_x(x,y) b_y(x,y)];
                        S_loc(i,j) = S_loc(i,j) + ...
                              (beta * (B'\grad_phi_i)) ...
                            * (beta * (B'\grad_phi_j) + Dati.sigma .* dphiq(1,q,j)) ...
                            * det(B) * w_2D(q);
                    end
                end
            end

             tau_k = 1 / sqrt(1/(Dati.dt)^2 + 30*Dati.mu^2 /femregion.h^4 + norm(beta)^2 / femregion.h^2);
             %tau_k = 1 / sqrt(30*Dati.mu^2 /femregion.h^4 + norm(beta)^2 / femregion.h^2);
             %tau_k = tau;
             %tau_k = 1 / sqrt(30*Dati.mu^2 /femregion.h^4 + norm(beta)^2 / femregion.h^2);
             tauu = tau_k;
            
             Ms_loc=zeros(nln,nln);

            for i=1:nln
                for j=1:nln
                    for k=1:length(w_2D)
                        Binv = inv(BJ(:,:,k));      % inverse
                        Jdet = det(BJ(:,:,k));      % determinant 
                        Ms_loc(i,j) = Ms_loc(i,j) + (Jdet.*w_2D(k)) .* dphiq(1,k,j).* beta * ((Grad(k,:,i) * Binv)');
                    end
                end
            end

           % h_k = polygon_diameter(femregion.dof(iglo,:));
            A(iglo,iglo) = A(iglo,iglo) + tau_k*S_loc;
            Ms_glob(iglo,iglo) = Ms_glob(iglo,iglo) + tau_k * Ms_loc;


            g = zeros(nln,1);
            x = pphys_2D(:,1);                         
            y = pphys_2D(:,2);
            mu = Dati.mu;
            F = eval(Dati.force);

            for i = 1:nln
                for k = 1:length(w_2D)
                    B = BJ(:,:,q); % 2x2
                    Binv = inv(B);
                    x = pphys_2D(k,1);                         
                    y = pphys_2D(k,2);
                    Jdet = det(BJ(:,:,k));  % determinant 
                    beta = [b_x(x,y) b_y(x,y)];
                    g(i) = g(i) + w_2D(k) * Jdet * F(k) * beta * ((Grad(k,:,i) * Binv)');
                end    
            end

            f(iglo) = f(iglo) + tau_k * g;
    end
end
end
