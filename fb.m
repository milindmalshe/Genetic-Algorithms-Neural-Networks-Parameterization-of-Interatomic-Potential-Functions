%%%returns fC
function [fB] = fb(i,j,r,B,lambda2)

rij=r(i,j);
fB = B*exp(-1.*lambda2.*rij);