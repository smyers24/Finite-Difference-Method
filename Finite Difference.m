clear all;
clc;

%Open data files and parse them into variables
stress=fopen('stress.dat','r'); 
strain=fopen('strain.dat','r');
x=csvread('stress.dat');
y=csvread('strain.dat');

%take the mean of the respective data points
stressmean = mean(x);
strainmean = mean(y);

%Bisection Method
maxiter=100;
a = 30;
b = 40;

%for modularity, the user could be allowed to enter their desired tolerance
%delta = input('Please enter root tolerance.\n');
%eps = input('Please enter residual tolerance.\n');

delta = 0.000001; eps = 0.000001;
E0=1; 
%equation with left end point plugged in
u = (E0 / a)*(exp(a*strainmean) - 1) - stressmean;

%function value at left point
e = b - a; % interval length
for k = 1:maxiter
e = e*0.5; 
% shrink interval by half
c = a + e; 
%update middle point
w = (E0 / c)*(exp(c*strainmean) - 1) - stressmean;
%function value at middle point
    if (abs(e) < delta || abs(w) < eps)
    break;
    end
    if ((u > 0) && (w < 0)) || ((u < 0) && (w > 0)) 
    b = c; 
    else
a = c; u = w;
    end
    fprintf('Nr. of iterations it took the bisection method to converge to root: %f \n',k);
    fprintf('Approximated root found by the bisection method is: %f \n',c);
end

n=12;

%Forward and Backward difference schemes
for i=1:n
    FDiff(i) = (x(i+1) - x(i)) / (y(i+1) - y(i)); 
end

for i=n+1:-1:2
    BDiff(i) = (x(i) - x(i-1)) / (y(i) - y(i-1)); 
end

%Plot Forward Difference scheme
plot(FDiff(5:12),x(5:12),'ro');
hold on
%Plot Backward Difference scheme
plot(BDiff(6:13),x(6:13),'gd');
hold off
title('Forwards and Backwards Diff.');
xlabel('Forward and Backward Difference')
ylabel('Stress')

% Least Square Method with Forward Diff.
Ex = 0; Exy = 0; st=0; Ey = 0; Ex2=0; sr=0; 
	
for i = 5:12
    Ex = Ex + x(i);
		Ey = Ey + FDiff(i);
		Exy = Exy + (x(i) * FDiff(i));
		Ex2 = Ex2 + x(i)^2;
end
     
xm = Ex/8;
ym = Ey/8;
a1 = (8*Exy - Ex*Ey) / (8*Ex2 - Ex^2);
a0 = ym - (a1*xm);

for i = 5:12
		st = st + (FDiff(i) - ym)^2;
		sr = sr + (FDiff(i) - a0 - a1*x(i))^2;
end

r2 = (st - sr) / st;
% Outputting solution to  system
fprintf('The line of regression is \n y= % fx + % f \n',a1,a0);
fprintf('r^2 = % f \n',r2);

%Plot Line of regression for Forward Difference scheme
figure, plot(FDiff(5:12), x(5:12), 'ro', a0+a1*x(5:12), x(5:12), 'b-');
title('Line of regression from Forward Difference');
xlabel('Forward Difference Data')
ylabel('Stress')

% Least Square Method with Backward Diff.
Ex = 0; Exy = 0; st=0; Ey = 0; Ex2=0; sr=0; 
	
for i = 6:13
    Ex = Ex + x(i);
		Ey = Ey + BDiff(i);
		Exy = Exy + (x(i) * BDiff(i));
		Ex2 = Ex2 + x(i)^2;
end
     
xm = Ex/8;
ym = Ey/8;
a1B = (8*Exy - Ex*Ey) / (8*Ex2 - Ex^2);
a0B = ym - (a1B*xm);

for i = 6:13
		st = st + (BDiff(i) - ym)^2;
		sr = sr + (BDiff(i) - a0B - a1B*x(i))^2;
end

r2 = (st - sr) / st;
% Outputting solution to system
fprintf('The line of regression is \n y= % fx + % f \n',a1B,a0B);
fprintf('r^2 = % f \n',r2);

%Plot Line of regression for Backward Difference scheme
figure, plot(BDiff(6:13),x(6:13), 'ro', a0B+a1B*x(6:13), x(6:13), 'b-');
title('Line of regression for Backward Difference');
xlabel('Backwards Difference')
ylabel('Stress')

FNew = (a0 / a1)* (exp(a1*y)-1);
BNew = (a0B / a1B)* (exp(a1B*y)-1);

%Plotting analytic curves
figure, plot(y,FNew, 'r', y,x, 'b');
title('Forward Analytic Curve');
xlabel('Strain')
ylabel('New Forward Difference')
figure, plot(y,BNew, 'r', y,x, 'b');
title('Analytic Backwards Diff');
xlabel('Strain')
ylabel('New Forward Difference')

%Plot new analytic curves
EA=(stressmean/(exp(a1*strainmean)-1))*(exp(a1*y)-1);
figure, plot(y,EA, 'r', y,x, 'b');
title('New Forward Difference Analytic Curve');
xlabel('Strain')
ylabel('New Forward Difference Analytic')
EB=(stressmean/(exp(a1B*strainmean)-1))*(exp(a1B*y)-1);
figure, plot(y,EB, 'r', y,x, 'b');
title('New Backward Difference Analytic Curve');
xlabel('Strain')
ylabel('New Forward Difference Analytic')
