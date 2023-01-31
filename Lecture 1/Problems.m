

%% -------------- PROBLEM 1.1 -------------- 

x1 = -5:0.005:5;
x2 = -5:0.005:5;
[X1,X2] = meshgrid(x1,x2);
F = X1.^2 - 2*X1 + 3*X1.*X2 + 4*X2.^2;
v = [0:2:10 10:10:100 100:20:200];
[c,h]=contour(X,Y,F,v,"LineWidth",2);
xlabel("x1", "FontSize", 14);
xlabel("x2", "FontSize", 14);
colorbar;

%yc1 = (x+2).^2;
%yc2 = (4*x)/10;
%hold on
%    fill(x,yc1,[0.7 0.7 0.7],"FaceAlpha",0.2)
%    fill([x x(end) x(1)],[yc2 -5 -5],[0.7 0.7 0.7],"FaceAlpha",0.2)
%hold off

%%


%% -------------- PROBLEM 1.2 -------------- 

% Then we need to derive an analytical expression for the gradient and the
% hessian.

% DONE

%%

%% -------------- PROBLEM 1.3 -------------- 

% Then we need to derive an analytical expression for the gradient and the
% hessian.

% DONE

%%

%% -------------- PROBLEM 1.4 -------------- 

x0 = [2 3];
[fVal, grad, H] = FunEx1(x0);

%%

%% -------------- PROBLEM 1.5 -------------- 

f = @(x) x(1).^2 - 2*x(1) + 3*x(1).*x(2) + 4*x(2).^2;
x0 = [2 3];
[fValApprox, gradApprox, HApprox] = FiniteDifference(f, x0, 0.1);
fprintf("The function value f(x) is: \n");
fValApprox
fprintf("The function value grad(f(x)) is: \n");
gradApprox
fprintf("The function value H(f(x)) is: \n");
HApprox


%%



%% -------------- PROBLEM 3.1: Derivatives of a Multivariate Vector Function -------------- 

x0 = [1.7, 2.1];
f = @(x) [exp(x(1)) - x(2); x(1).^2 - 2.*x(2)];
[c, dc, d2c] = FunEx2(f, x0);

%%
