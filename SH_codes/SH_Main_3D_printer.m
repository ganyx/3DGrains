%% The main funtion to implement the sphercial harmonic analysis to generate virtual particle shapes
% If you use this code for publication, please cite our paper in it.
% by Wei Deheng, Zhongzheng Wang, Jean-Michel Pereira and Yixiang Gan.
clc; 
clear;
% spherical expansion and reconstructon
clear sph_verts faces Z; load('L3_icosa.mat');
maxDeg = 15 ; % the degree for spherical basis calculation
reconDeg = 15  ; % the degree used for spherical harmonic reconstrution
% fvec = SH_expansion(vertices,maxDeg,Pid); % Rotational spherical expansion
for cc=1:12 %:11
    fractal= 2.0+cc*0.05; %2+(cc-1)*0.1;
for ccc=1:10
    D_2= 0.03*ccc; %0.01*(ccc-1);% 0.1+0.01*(ccc-1);
for bb=1:1000
    FD(bb,1)=normrnd(fractal*2-6,0);% 0.0909074*2, noise
end
for c=1:1000
n = 1;
while n<16
    J = rand(n+1,1)-rand(n+1,1);
    K = flipud(J(1:n,1));
    A(n,1) = [(-1)^n];
    B = K.*A;
    L(n^2:(n+1)^2-1,1) = [J;B];
    M = rand(n,1)-rand(n,1);
    C(n,1) = [(-1)^(n+1)];
    N(n^2:(n+1)^2-1,1) = [M;[0];flipud(M).*C];
    O = L+N*i;
    P = [1;O];
    n = n+1;
end
Q = conj(P);
b=1;
while b<16
    R(1,1) = sqrt(P(1,1).*Q(1,1));
    R(b+1,1) = sqrt(sum(P(b^2+1:(b+1)^2,1).*Q(b^2+1:(b+1)^2,1)));
    b = b+1;
end
Nor(1,1:2)=[1 0];
Nor(2,1:2)=[0 0];
Nor(3,1:2)=[D_2 0];% 0.0428
% Nor(3,2)=Nor(3,2)*2;
for aa=4:16
    Nor(aa,1)=Nor(3,1)*((aa-1)/2)^(FD(c,1));
    Nor(aa,2)=0.*Nor(aa,1); % 0.3.*Nor(aa,1);
end
for f=1:16
    S(f,1) = normrnd(Nor(f,1),Nor(f,2));
end
T = S./R;
d = 0;
while d<16
fvec(d^2+1:(d+1)^2,1) = P(d^2+1:(d+1)^2,1).*T(d+1,1);
d=d+1;
end
fvec=fvec*3.5;
spharm_verts=SH_reconstruction(faces,sph_verts,reconDeg,fvec);  
vertices = spharm_verts;
number=c+(ccc-1)*1000+(cc-1)*10000;
strcat(number);
filename = strcat(num2str(number))
save(filename,'vertices','faces','fvec');
J=[];K=[];A=[];B=[];L=[];M=[];C=[];N=[];O=[];P=[];Q=[];R=[];S=[];T=[];U=[];V=[];
end
end
end