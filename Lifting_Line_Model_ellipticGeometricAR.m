%Original: Tianshu Liu, WMU; Modified by Keaton Warn to write f(AoA) data to files 12/11/2020
%CHECK: 1. Elliptical (r=1e-15) or bird (r=1e-6) planform? Variable AoA or AR? alpha...(N,1) or (ydb)? 
clear all;
close all;

% define the cuts and their positions in the y-coordinate y/b
% (span coordinate)
e1=0.0001;
y_left=-(0.5-e1);   % left wing tip
y_right=(0.5-e1);   % right wing tip
N_segment=50;       % number of cuts
dy=(y_right-y_left)/N_segment;
ydb=[y_left:dy:y_right];  % span locations of the cuts

% define the wing planform or chord distribution 
b=8;    % span
c0=1;   % wing root chord
V_inf=10; % freestream velocity

% elliptic planform
c=c0*(1-(2*ydb).^2).^0.5;

% bird planform
% for i=1:length(ydb)
%     c(i)=c0*Bird_Planform(ydb(i));
% end

% calculate the wing area and aspect ration (AR)
%S=c0*b*trapz(ydb,c);
%AR=b^2/S;%---------------------------------------------------------------AR

for i=1:length(ydb)
    x_LE(i)=c(i)/4/b;
    x_TE(i)=-c(i)*(3/4)/b;
end

% ansformation from y/b to the angular variable thetaT
theta=acos(-2*ydb);

% calculate the right-hand vector B based the AoAs
    %Set up loop parameters
%AR Range
AR=0; %starting aspect ratio value
AR_end=35;%b^2/S;
AR_range=AR_end-AR;
res=100; %Discrete resolution (evaluation points per AR_range)
step=(AR_range)/res; %how much the AR value is increased by per loop iter
for k=1:(res+1)
    AR=(step*(k-1));   %AR increased from 0 by (b^2/S)/res each iter

    b=AR*c0*trapz(ydb,c); %changing b (span) for downstream computations
    
    alpha_inf=3;                            %frestream AoA in degrees°
    N=length(ydb);                          %number of cuts
    %alpha=(alpha_inf*pi/180)*ones(N,1);     %define the local geometrical AoA
    alpha=(alpha_inf*pi/180)*(1+1*abs(ydb));
    alpha_in=(0*pi/180)*ones(N,1);          %define the zero-lift AoA
    B=alpha-alpha_in;

    % calculate the elements of the matrix C
    for m=1:N
        for n=1:N
            C(m,n)=2*b*sin(n*theta(m))/(pi*c(m))+n*sin(n*theta(m))/sin(theta(m));
        end
    end

    % solve the system A*C=B for the coefficients in the Fourier series
    r=1e-15;
    if(rcond(C)>r)
       A=inv(C)*B;
       else
       [U,S,V]=svd(C); %singular value decomposition. [A]=[U]*[S]*[V]'
       sv=diag(S);
       for i=1:length(sv)
          if sv(i)>max(max(sv))*r
             rsv(i)=1/sv(i);
          else rsv(i)=0;
          end
       end
       S1=diag(rsv);
       A=V*S1*(U'*B);
    end

    % calculate the vortex strength distribution, gamma
    for m=1:N
        gamma_add=0;
        n=1;
        while n<=N
            gamma_add=gamma_add+A(n)*sin(n*theta(m));
            n=n+1;
        end
        gamma(m)=gamma_add;
    end
    gamma=2*b*V_inf*gamma;

    % calculate the lift coefficient
    CL=pi*A(1)*AR;

    % calculate the span efficiency--------------------------------------e
    delta_add=0;
    for i=2:N
        delta_add=delta_add+i*(A(i)/A(1))^2;
    end

    delta=delta_add;
    e=1/(1+delta); % span efficiency

    % calculate the induced drag
    CDi=CL^2/(pi*AR*e);
    
    % output data
    data_out=[k; alpha_inf; AR; CL; CDi; e];
    
    for m=1:N
        arrGamma(k,m)=gamma(m);
    end
    
    for j=1:6
        arrData(j,k)=data_out(j);
    end
end
%Write to .txt file (read by Windows Notepad)
%fileID=fopen('LLME_AR4_Data.txt','w');
%fprintf(fileID,'%0s %12s %16s %20s %25s %30s\r\n','k','alpha_inf','AR', 'CL','CDi','e');
%fprintf(fileID, '%0f %12f %16f %20f %25f %30f\r\n',arrData);
%fclose(fileID);

%Write to Excel
%if loop # of iterations changed, remember to manually clear cell data
%2D
AD=arrData.';
T=array2table(AD);
fileName = 'FTD7.xlsx';
header = {'k','alpha_inf','AR','CL','CDi','e'};
xlswrite(fileName,header,1,'A1');
xlswrite(fileName,AD,1,'A2');
%EXCEL PLOT: 3D plot data across changing AoA
for i=1:(res+1)
    arrAR(i)=arrData(3,i);
end
AT=arrAR.';
fileName='FTD7 - 3DPlot.xlsx';
xlswrite(fileName,ydb,1,'B1');
xlswrite(fileName,arrGamma,1,'B2');
xlswrite(fileName,AT,1,'A2');

%MATLAB PLOTS
% plot the chord distribution
figure(1);
plot(ydb,c/c0,'-o');
grid;
xlabel('y/b');
ylabel('Chord Length Distribution, c/c_0');

% plot the vortex strength distribution
figure(2);
plot(ydb,gamma,'-o');
grid;
xlabel('y/b');
ylabel('Circulation Distribution, Gamma');

%PLOT 3 plot of the wing planform
figure(3);
plot(ydb,x_LE,'-o',ydb,x_TE,'-o');
grid;
xlabel('y/b');
ylabel('x/b');
axis image;
title('Wing Planform');




