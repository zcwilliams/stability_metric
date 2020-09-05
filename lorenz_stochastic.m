function [XYZ]=lorenz_stochastic(r,s,b,sigma,tf,dt)

t0      = 0;
tspan   = t0:dt:tf;
options = sdeset('RandSeed',1,'SDEType','ito','DiagonalNoise','yes');

% specify multiplicative noise form
g = @(t,X)[sigma*X(1);...
           sigma*X(2);...
           sigma*X(3)];  

% specify model equations
f = @(t,X)[-s*X(1) + s*X(2);...         % x
           r*X(1) - X(2) - X(1)*X(3);...% y
           -b*X(3) + X(1)*X(2)];        % z

% sde solve
[XYZ] = sde_euler(f,g,tspan,randn(1,3),options); 