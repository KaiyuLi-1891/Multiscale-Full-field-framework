function [out, sol] = MatSolver(ca1, cx)
%MATSOLVER Quasi-1D Poisson-Nernst-Planck boundary-layer solver.
%
% This is the cleaned, drop-in version used with the accompanying COMSOL
% model. The numerical formulation and hard-coded parameters are retained
% from the original implementation.
%
% Inputs
%   ca1 : cation concentration at z = h (mol/m^3)
%   cx  : anion concentration at z = h (mol/m^3)
%
% Outputs
%   out : electric displacement epsilon*E at z = h (C/m^2)
%   sol : bvp4c solution structure from the final local solve

h = 20*10^(-9);              % Quasi-1D domain thickness (m)
i_s = 2000;                  % Current density (A/m^2)
J = i_s/96500;               % Prescribed anion molar flux (mol/m^2/s)
epsilon = 6.9505E-10;        % Permittivity (F/m)
DA = 1.96e-9;                % Cation diffusion coefficient (m^2/s)
DX = 5.270e-9;               % Anion diffusion coefficient (m^2/s)
R = 8.3145;                  % Gas constant (J/mol/K)
F = 96500;                   % Faraday constant (C/mol)
T = 298.15;                  % Temperature (K)

% Electric field imposed at the Stern-layer boundary (V/m).
E_0 = -log10(i_s/0.1)*0.12/(0.5*10^(-9));

% Initial mesh for the boundary-value problem.
xmesh = 0:0.1*10^(-9):h;

CA = 0;
CX = 0;
solinit = bvpinit(xmesh, @guess);

CAmatrix = ca1;
CXmatrix = cx;
sizeofmatrix = size(CAmatrix);
sizeofmatrix = sizeofmatrix(1);
D_hmatrix = [];

% Solve one independent quasi-1D boundary layer for each COMSOL input row.
for i = 1:sizeofmatrix
    CA = CAmatrix(i);
    CX = CXmatrix(i);
    sol = bvp4c(@bvpfcn, @bcfcn, solinit);
    E_h = sol.y(1,end)*10^6;
    D_h = E_h*epsilon;
    D_hmatrix = [D_hmatrix; D_h]; %#ok<AGROW>
end

out = D_hmatrix;

    function dydx = bvpfcn(~, y)
        % State scaling:
        %   y(1)*1e6 = electric field E (V/m)
        %   y(2)*1e3 = cation concentration (mol/m^3)
        %   y(3)     = anion concentration (mol/m^3)
        dydx = zeros(3,1);
        dydx = [F*(y(2)*10^3-y(3))/epsilon/10^6; ...
                (DA*y(2)*10^3*F/(R*T)*y(1)*10^6)/DA/10^3; ...
                (-J+DX*y(3)*(-1)*F/(R*T)*y(1)*10^6)/DX];
    end

    function res = bcfcn(ya, yb)
        % Boundary conditions:
        %   E(0) = E_0
        %   cation(h) = CA
        %   anion(h) = CX
        res = [ya(1)-E_0/10^6;
               yb(2)-CA/10^3;
               yb(3)-CX];
    end

    function g = guess(x)
        % Original initial guess retained for numerical consistency.
        g = [-1/(x+1);
              1/(x+1);
              CX/h*x];
    end
end
