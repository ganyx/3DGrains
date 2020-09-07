%% This function is used for calculating the spherical harmonic coefficiencs

function  fvec = SH_expansion(vertices,maxDeg,Pid)
resultdir = fullfile(cd, 'Ini_coeff');
vertnum = size(vertices,1);
max_d = maxDeg;
% Note that degree 'd' we want to use depends on the vertnum 
% The total number of unknowns is (d+1)*(d+1)
% The total number of equations is vertnum
% We want equ_num >= unk_num
deg = max(1, floor(sqrt(vertnum)*1/2));
deg = min(deg, max_d);
disp(sprintf('Use spharm up to %d degree (vec_len=%d).',deg,(deg+1)^2));
center = mean(vertices,1);
for i = 1:size(vertices,1)
    vertices(i,:) = vertices(i,:)-center; % the polar radius from the center
    rads(i,:) = norm( vertices(i,:));
end
% visulization of the initial image
% subplot(1,2,1);
% patch_lighta(vertices, faces); hold on;

% calculate the spherical basis
Z = calculate_SPHARM_basis(vertices, deg); %% vertices is the direct mapping
[x,y] = size(Z);
disp(sprintf('Least square for %d equations and %d unknowns',x,y));

% Least square fitting
 fvec = Z\rads;   %This does not work as it is expected to work in certain environment
%  if ~exist(resultdir,'dir')
%      mkdir(resultdir);
%  end
%  new_name = [resultdir '\' [Pid, '_LSM_COEF.mat']];
%  save(new_name, 'vertices', 'faces', 'fvec');

return;

