
clc;clear all;
close all;

N=50; %
Function_name=2;
%[lb,ub,dim,fobj]=Get_Functions_details(Function_name); % benchmark  Function 
[dim,lb,ub, Vio, GloMin, fobj]=ProbInfo(Function_name); %engineering problems
T=2000;
[xmin,fmin,CNVG]=IHBA(fobj,dim,lb,ub,T,N,Vio);

disp('--------Best Fitness--------');
display(['The average fitness of IHBA is: ', num2str(fmin)]);


figure,
semilogy(CNVG,'r')
xlim([0 T]);
title('Convergence curve')
xlabel('Iteration');
ylabel('Best fitness obtained so far');
legend('IHBA')
display(['The best location= ', num2str(xmin)]);
display(['The best fitness score = ', num2str(fmin)]);

        

