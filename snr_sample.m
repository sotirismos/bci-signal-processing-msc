function [SP,NP]=snr_sample(X)
%
% [sp,np]=mine_snr(X)
% snr of the curves / wavelets  tabulated  in a sample X 


[N,p]=size(X); ave=mean(X);

d=dmatrix(X);

v=ones(1,N);
   
NP= ( 1/(2 * p*N*(N-1)) ) * ( v*d*v');

SP=(1/p)* (ave*ave') - (1/N)*NP ; 
 
