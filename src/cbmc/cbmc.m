clc; clear all; close all;
k=50.0;
theta0=120;
theta0=theta0*pi/180;

N=10000;
theta1=zeros(N,1);
theta2=zeros(N,1);
theta3=zeros(N,1);

theta1(1)=theta0;
theta2(1)=theta0;
theta3(1)=theta0;
for i=2:N
    theta1(i)=theta0+randn(1)/sqrt(k);
    r2=[-cos(theta1(i)) sin(theta1(i)) 0.0];
    
    while 1
        theta2_trial=rand(1)*pi;
        phi_trial=rand(1)*2*pi;
        r3=[-cos(theta2_trial)  sin(theta2_trial)*cos(phi_trial)  ...
            sin(theta2_trial)*sin(phi_trial)];
        
        theta3_trial=acos(dot(r2,r3));
        arg=exp(-k/2*(theta2_trial-theta0)^2)*exp(-k/2*(theta3_trial-theta0)^2);
        if rand(1) < arg
            theta2(i)=theta2_trial;
            theta3(i)=theta3_trial;
            break;
        end
    end
end

edges=linspace(0,pi,101);
plot(edges(1:100),histcounts(theta1,edges),'k-', ...
    edges(1:100),histcounts(theta2,edges),'r-', ...
    edges(1:100),histcounts(theta3,edges),'b-')
