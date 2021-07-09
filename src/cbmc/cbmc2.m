clc; clear all; close all;
k=50.0;
theta0=120;
theta0=theta0*pi/180;

N=100000;
theta1=zeros(N,1);
theta2=zeros(N,1);
theta3=zeros(N,1);
phi=zeros(N,1);

theta1(1)=theta0;
theta2(1)=theta0;
theta3(1)=theta0;
phi(1)=theta0;
NT=100;
for i=2:N
    theta1_trial=theta0+randn(1)/sqrt(k);
    theta2_trial=theta0+randn(1)/sqrt(k);
    phi_trial=rand(NT,1)*2*pi;
    r2=[-cos(theta1_trial) sin(theta1_trial) 0.0];
    r3=[-cos(theta2_trial)*ones(NT,1)  sin(theta2_trial)*cos(phi_trial)  ...
            sin(theta2_trial)*sin(phi_trial)];
    theta3_test=acos(dot(repmat(r2,NT,1),r3,2));
    w=exp(-k/2*(theta3_test-theta0).^2);
    t=rand(1)*sum(w);
    ws=0;
    for j=1:NT
        ws=ws+w(j);
        if t < ws
            theta3_trial=theta3_test(j);
            break;
        end
    end
    rosen_new=sum(w);
    
    phi_trial=[phi(i-1); rand(NT-1,1)*2*pi];
    r2=[-cos(theta1(i-1)) sin(theta1(i-1)) 0.0];
    r3=[-cos(theta2(i-1))*ones(NT,1)  sin(theta2(i-1))*cos(phi_trial)  ...
            sin(theta2(i-1))*sin(phi_trial)];
    theta3_test=acos(dot(repmat(r2,NT,1),r3,2));
    w=exp(-k/2*(theta3_test-theta0).^2);
    rosen_old=sum(w);
    
    arg=rosen_new/rosen_old;
    if rand(1) < arg
        theta1(i)=theta1_trial;
        theta2(i)=theta2_trial;
        theta3(i)=theta3_trial;
    else
        theta1(i)=theta1(i-1);
        theta2(i)=theta2(i-1);
        theta3(i)=theta3(i-1);
    end

end

edges=linspace(0,pi,101);
plot(edges(1:100)*180/pi,histcounts(theta1,edges),'k-', ...
    edges(1:100)*180/pi,histcounts(theta2,edges),'r-', ...
    edges(1:100)*180/pi,histcounts(theta3,edges),'b-')
