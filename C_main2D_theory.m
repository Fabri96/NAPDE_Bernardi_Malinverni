function [errors,solutions,femregion,Dati,tauu]=C_main2D_theory(TestName,nRef,mu)
%==========================================================================
% Solution of the Poisson's problem with linear finite elements
% (non homogeneous Dirichlet boundary conditions)
%==========================================================================
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          nRef        : (int)     refinement level
%          mu          : (int)     diffusion coefficient
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) infos about finite elements
%                        discretization
%          Dati        : (struct)  see C_dati.m
%          
% Usage: 
%    [errors,solutions,femregion,Dati] = C_main2D('Test1',3,0.001)
 


addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


%==========================================================================
% LOAD DATA FOR TEST CASE
%==========================================================================

Dati = C_dati(TestName);
Dati.nRefinement = nRef;
Dati.mu=mu;
mu = Dati.mu;

switch TestName
    case 'Test1'
        f = @(x,y,mu) - sin(2*pi*x) .* (y.^3-y) - mu .* (-4*pi^2 * sin(2*pi*x).*(y.^3-y) + 6.*y.*sin(2*pi*x))  + 100*2*pi*cos(2*pi*x).*(y.^3-y) + 100*sin(2*pi*x).*(3*y.^2-1);
        
    case 'Test2'
        f = @(x,y,mu) -(sin(2*pi*x)).^2.*sin(2*pi*y) - mu*(8*pi^2*sin(2*pi*y).*((cos(2*pi*x)).^2-(sin(2*pi*x)).^2)-4*pi^2*(sin(2*pi*x)).^2.*sin(2*pi*y)) ...
            +100*(4*pi.*sin(2*pi*y).*sin(2*pi*x).*cos(2*pi*x)+2*pi.*(sin(2*pi*x)).^2.*cos(2*pi*y));
        
    case 'Test3'
        f = @(x,y,mu) -10*exp(4*y).*(y.^3-7/4*y.^2+3/4*y).*(x.^3-3/2*x.^2+1/2*x) ...
            -mu*(5.*exp(4.*y).*(2.*x - 1).*(32.*x.^2.*y.^3 - 8.*x.^2.*y.^2 - 20.*x.^2.*y + 5.*x.^2 - 32.*x.*y.^3 + 8.*x.*y.^2 + 20.*x.*y ...
            - 5.*x + 12.*y.^3 - 21.*y.^2 + 9.*y))/2 ...
            +100*( 5.*x.*exp(4.*y).*(2.*x.^2 - 3.*x + 1).*(3.*y.^2 - (7.*y)./2 + 3/4) + (5.*y.*exp(4.*y).*(3.*x.^2 - 3.*x + 1/2).* ...
            (4.*y.^2 - 7.*y + 3))/2 + 5.*x.*y.*exp(4.*y).*(2.*x.^2 - 3.*x + 1).*(4.*y.^2 - 7.*y + 3) );
end

Dati.force = strrep(char(f),'@(x,y,mu)','');

%==========================================================================
% MESH GENERATION
%==========================================================================

[Region] = C_create_mesh(Dati);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[femregion] = C_create_femregion(Dati,Region); 

%==========================================================================
% BUILD FINITE ELEMENT MATRICES (M and A) and RIGHT-HAND SIDE (b)
%==========================================================================

[M_no_bc,A_no_bc,Adv_glob,Ms_glob,b_no_bc,tauu] = C_matrix2D_theory(Dati,femregion);
% time integration parameters
t = 0;              %initial time
T = Dati.T;         %final time
dt = Dati.dt;       %time-step
theta = Dati.theta; %theta-method

err_time = zeros(0.01 * T/dt,1);
err_timeH = zeros(0.01 * T/dt,1);
ii = 1; %counter

x = femregion.coord(:,1);
y = femregion.coord(:,2);

% Evaluation of initial condition
u0 = eval(Dati.initialcond);

% Snapshot of the initial condition
% % % C_snapshot_sol(femregion, u0, Dati, t)

% evaluation of force at time tk
f0 = eval(Dati.forcetime);

% Definition of the matrix
K_no_bc = (M_no_bc + Ms_glob) + theta * dt * (A_no_bc + Adv_glob);

for t = dt : dt : T 

    % evaluation of force at time tk+1
    f1 = eval(Dati.forcetime);
    
    % evaluation of rhs
    rhs_no_bc = dt * b_no_bc * (theta*f1 + (1-theta)*f0) + (M_no_bc - (1 - theta) * dt * A_no_bc) * u0;
    
    %========================================================
    % COMPUTE BOUNDARY CONDITIONS -- MODIFICATION OF K and rhs
    %========================================================

    [K,rhs,u_g] = C_bound_cond2D(K_no_bc,rhs_no_bc,femregion,Dati,t);
    
    %==========================================
    % SOLVE THE LINEAR SYSTEM
    %==========================================

    KKKK = full(K);
    rhssss = full(rhs);
    u1_no_bc = pinv(KKKK)*rhssss;
    
    %==================================================================
    % ASSIGN DIRICHLET BOUNDARY CONDITIONS -- through the lifting ug
    %==================================================================

    u1 = u1_no_bc + u_g;    
    
    %==================================================================
    % Snapshot of the solution
    %==================================================================
 
%     if mod(round(t/dt), 100) == 0
%       C_snapshot_sol(femregion, u1, Dati, t)
%     end

    %==================================================================
    % UPDATE
    %==================================================================
    
    f0 = f1;
    u0 = u1;

%   to obtain error on the cylinder

    if rem(t/dt,100) == 0
        [sols] = C_postprocessing(Dati,femregion,u1,t);
        [E_L2, E_SEMI_H1] = C_error_L2_H1(femregion, sols.uh, Dati,t);
        err_time(ii) = E_L2;
        err_timeH(ii) = sqrt(E_L2.^2 + E_SEMI_H1.^2);
        ii = ii + 1;
    end
        
        
end
%C_snapshot_sol(femregion, u1, Dati, t)

%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

[solutions] = C_postprocessing(Dati,femregion,u1,t);

for i=1:length(solutions.uh)
    if solutions.u_ex(i) == 0
        solutions.uh(i) = 0;
    end
end


%==========================================================================
% ERROR ANALYSIS
%==========================================================================
errors = [];
if (Dati.plot_errors)
    [errors] = C_compute_errors(Dati,femregion,solutions,t);
end

err_in_time = (100*dt) * trapz(err_time);
err_in_timeH = (100*dt) * trapz(err_timeH);

errors.Time_error = err_in_time;
errors.Time_errorH = err_in_timeH;



