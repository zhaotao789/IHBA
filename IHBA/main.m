
clc;clear all;
close all;

N=50;T=2000;
% N=50; % Number of Squirrel
Function_name=5; % Name of the test function that can be from F1 to F23 

[xmin,fmin,CNVG]=IHBA(Function_name,T,N);

disp('--------Best Fitness--------');
display(['The best fitness of IHBA is: ', num2str(min(CNVG))]);

figure,
semilogy(CNVG,'r')
xlim([0 T]);
title('Convergence curve')
xlabel('Iteration');
ylabel('Best fitness obtained so far');
legend('HBA')

display(['The best location= ', num2str(xmin)]);
display(['The best fitness score = ', num2str(fmin)]);

        

