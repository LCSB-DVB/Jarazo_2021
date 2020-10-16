%{
SCRIPT:         distmatrix.m
SUPPLEMENT TO:  Chapter 13
DESCRIPTION:    This function generates a matrix of size M-by-N where all 
                values represent the distance from that pixel to the center
                of the image.
%}

%   Copyright 2011 O. Marques
%   Practical Image and Video Processing Using MATLAB, Wiley-IEEE, 2011.
%   $Revision: 1.0 Date: 2011/06/21 13:12:00 

function y = distmatrix(M,N)

u = 0:(M - 1);
v = 0:(N - 1);

ind_u = find(u > M/2);
u(ind_u) = u(ind_u) - M;
ind_v = find(v > N/2);
v(ind_v) = v(ind_v) - N;

[V, U] = meshgrid(v, u);

%calculate distance matrix
y = sqrt((U .^ 2) + (V .^ 2));