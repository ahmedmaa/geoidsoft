% This function computes legendre polynomial
% 
%
%
% Inputs:
%         t : cos(theta)
%         L : the maximum degree of expansion
% 
%
% Output: 
%         Pn: The Legendre Polynomial
%        wgf: factorial for WG Kernel
%        lsf: factorial for VK and LSM knernels
%                        
%                      
%
%                            Ahmed Abdalla
%                     Louisiana State University
%                          Jan 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pnn,wgf, lsf] = lgpoly(t,L)
% t is the cos(theta)
% m is the maximum degree
%  is the coefficients from degree 0 up to m
% load c.mat
% L=200;
Pnn=zeros(L+1,length(t));
wgf=zeros(L+1,1);
lsf=zeros(L+1,1);
Pnn(1,1:length(t)) = 1;    
Pnn(2,1:length(t)) = t; 
for n = 2:L 
    Pnn(n+1,:) = -(n-1)./n.*Pnn(n-1)+(2.*n-1)./n.*t.*Pnn(n);
    wgf(n+1,1) = (2.*n+1)./(n-1);
    lsf(n+1,1) = (2.*n+1)./2;
end