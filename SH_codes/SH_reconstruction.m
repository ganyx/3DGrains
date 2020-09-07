function spharm_verts=SH_reconstruction(faces,sph_verts,reconDeg,fvec)
resultdir = fullfile(cd, 'Norm_coeff');

[PHI,THETA] = cart2sph(sph_verts(:,1),sph_verts(:,2),sph_verts(:,3));
degree = reconDeg;
lb = 1;
ub = (degree+1)^2;
Z = calculate_SPHARM_basis(sph_verts, degree); 
spharm_rads = real(Z(:,lb:ub)*fvec(lb:ub,:));
min_rad = min(spharm_rads);
[rx, ry, rz] = sph2cart(PHI,THETA,spharm_rads);
spharm_verts = [rx,ry,rz];


% patch_lighta(spharm_verts, faces); hold on; 
% patch_lightmesh(spharm_verts, faces); hold on; 
% plot3([0,1.8*RR(1,1)],[0,1.8*RR(2,1)],[0,1.8*RR(3,1)],'linewidth', 2 ); hold on;
% plot3([0,-1.8*RR(1,1)],[0,-1.8*RR(2,1)],[0,-1.8*RR(3,1)],'linewidth', 2  ); hold on;
% plot3([0,1.5*RR(1,2)],[0,1.5*RR(2,2)],[0,1.5*RR(3,2)],'linewidth', 2  ); hold on;
% plot3([0,-1.5*RR(1,2)],[0,-1.5*RR(2,2)],[0,-1.5*RR(3,2) ],'linewidth', 2 ); hold on;
% plot3([0,1.2*RR(1,3)],[0,1.2*RR(2,3)],[0,1.2*RR(3,3)] ,'linewidth', 2 ); hold on;
% plot3([0,-1.2*RR(1,3)],[0,-1.2*RR(2,3)],[0,-1.2*RR(3,3)],'linewidth', 2  ); hold on;
% axis([-1.5,1.5,-1.5,1.5,-1.5,1.5]); 
return;

