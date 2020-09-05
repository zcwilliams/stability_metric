% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Accompanying code for manuscript:
%   
%   "Variations in stability revealed by temporal asymmetries in contraction of phase space flow" 
%   
%   Authors: Zachary C. Williams and Dylan E. McNamara
%   
%   Date: September 1, 2020. 
%
%   This script demonstrates the calculation of the stability metric, 'delta lambda', 
%   described in the manuscript, for the reconstructed Lorenz attractor. 
% 
%   Either the standard Lorenz-83  system or the Lorenz system with mulitiplicative 
%   noise can be evaluated. 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Contents
%   
%   1. stability_metric_demo.m  % Main demo file
%
%   2. lorenz_RK.m              % 4th order Runge-Kutta scheme solving
%                                 Lorenz system
%   3. lorenz_stochastic.m      % Stochastic Lorenz solver
%
%   4. delay_embed.m            % Attractor reconstruction delay embedding 
%   5. getlambda.m              % Calculate the forward and backward time
%                                 lambda (as described in manuscript)
%   
%   Note on choice of embedding time lag 
%   To determine a reasonable value of embedding time lag (tau), we use the average mutual
%   information function. A matlab function for this calculation can be found here: 
%          
%   https://www.mathworks.com/matlabcentral/fileexchange/10040-average-mutual-information
%
%   Tau is determined by taking the first minimum of the average mutual
%   information function. Tau is kept fixed when comparing stability
%   across systems. Example using output from the Lorenz model in
%   "stability_metric_demo.m":
%
%          [v,lag] = ami(XYZ(:,1),XYZ(:,1),[1:1000]);
%          findtau=find(diff(v)>0);
%          tau=lag(findtau(1));          

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Overview of stability_metric_demo.m
%
%    (0) Enter path to folder containing files.
%    (1) Generate time series output from Lorenz model, this is used to create the 
%            library set and test set. 
%    (2) Reconstruct a library and test attractor by delay embedding 
%            the solution output from the Lorenz system. Then use 
%            knnsearch to find the nearest neighbor (residing in the
%            library) for each point in the test trajectory.
%    (3) Calculate lambda^- and lambda^+ as a function of local time (L) for the
%            forward and backward trajectories. This is done in the function 'getlambda.m'.
%    (4) Plot (lambda^-)-(lambda^+) as a function of L. The metric 
%            'delta lambda', is the maximum of the curve (lambda^-)-(lambda^+),
%            displayed in figure(1). 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Required Toolboxes
%
%   1) SDETools toolbox (available from https://github.com/horchler/SDETools)
%       SDETools is used to generate time series from the Lorenz system with noise 
%   2) Statistics and Machine Learning toolbox (available from Mathworks)
%       for the functions knnsearch.m and createns.m.
%