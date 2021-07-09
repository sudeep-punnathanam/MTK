% This is to test arbitrary trial distributions in cbmc for bond stretch
% and angle bending

clc; clear all; close all;
rng('shuffle','twister');
k=10;
s=1/sqrt(k);
r=1.54;
u=k/2*(r-1.54)^2;

npts=1000000;
trj=zeros(npts,1);

s_trial=2*s;
for l=1:npts
    rnew=randn(1)*s_trial+1.54;
    utrial_old=((r-1.54)/s_trial)^2/2;
    utrial_new=((rnew-1.54)/s_trial)^2/2;
    unew=k/2*(rnew-1.54)^2;
    arg=(rnew/r)^2*exp(-((unew-utrial_new)-(u-utrial_old)));
    if rand(1) < arg
        u=unew;
        r=rnew;
    end
    trj(l)=r;
end

bin_low=0.5;
bin_high=3.0;
[n,edges]=histcounts(trj,'BinWidth',0.01,'BinLimits',[bin_low, bin_high],'Normalization','pdf');
plot(edges(1:end-1)+0.005,n,'k.');
x=linspace(bin_low,bin_high,1001);
y=1/sqrt(2*pi)/s*exp(-0.5*((x-1.54)/s).^2);
hold on;
plot(x,y,'b');
y=(x.^2).*exp(-0.5*((x-1.54)/s).^2);
C=trapz(y)*(bin_high-bin_low)/1000;
plot(x,y/C,'r');
hold off;
