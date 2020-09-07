%% The main funtion to implement the sphercial harmonic analysis
clc; clear;
%% The input of particle Number and result dictionary
for n=10
Pid = ['Surf10000-' num2str(n)]; 
% Input the Data after spherical paremetrazation
inputname = [Pid, '_CALD_ini.mat'];
inputdir = fullfile(cd, 'Surf_input\ParameterizePart_Origin');

load([inputdir '\' inputname]);
center = mean(vertices,1);
vertices = vertices -center(ones(size(vertices,1),1),:);
figure(2*n-1)
patch_lighta(vertices, faces);
axis([-1.5*10^3,1.5*10^3,-1.5*10^3,1.5*10^3,-1.5*10^3,1.5*10^3]);

%% spherical expansion and reconstructon
clear sph_verts faces Z; load('L3_icosa.mat');
maxDeg = 15 ; % the degree for spherical basis calculation
reconDeg = 15  ; % the degree used for spherical harmonic reconstrution
fvec = SH_expansion(vertices,maxDeg,Pid); % Rotational spherical expansion
spharm_verts=SH_reconstruction(faces,sph_verts,reconDeg,fvec);  
figure (2*n);
patch_lightmesh(spharm_verts, faces);    
axis([-1.5*10^3,1.5*10^3,-1.5*10^3,1.5*10^3,-1.5*10^3,1.5*10^3]);
Fvec(:,n)=fvec;
end


