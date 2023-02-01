
function [Xprey, Food_Score,CNVG] = IHBA(objfunc, dim,lb,ub,tmax,N,Vio)
beta       =6;     % the ability of HB to get the food  Eq.(4)
gamma=3;
w       = 0.8;     %constant in Eq. (3)
vec_flag=[1,-1];
CNVG=zeros(tmax,1);
%initialization
X=initialization(N,dim,ub,lb);
%Evaluation
Xnew=X;
fitness = fun_calcobjfunc(objfunc, X,Vio);
[GYbest, gbest] = min(fitness);
Xprey = X(gbest,:);
for t = 1:tmax
    %alpha=C*exp(-t/tmax);   %density factor in Eq. (3)
    a=(cos(2*rand)+1)*(1-t/tmax);
    alpha=a*(2*rand-1); 
    C=alpha;
    I=Intensity(N,Xprey,X); %intensity in Eq. (2)
    for i=1:N
        F=vec_flag(floor(2*rand()+1));        
            di=((Xprey-X(i,:)));
            if rand<.5
                Z = unifrnd(-a,a,1,dim);    %生成(连续)均匀分布的随机数   其中a是约束
                H=Z.*X(i,:);   %定义是放大还是缩小
                r3=rand;                r4=rand;                r5=rand;                
                Xnew(i,:)=C.*H+(rand-a)*X(randi([1,N]),:)+F*beta*I(i)* Xprey+F*r3*alpha*(di)*abs(cos(2*pi*r4)*(1-cos(2*pi*r5)));
            else
                r7=rand;
                Xnew(i,:)=X(i,:)-C.*(C*(X(i,:)- Xnew(randi([1,N]),:))+rand*(X(i,:)-Xnew(randi([1,N]),:)))+F*r7*alpha*di;
            end     
        
    end
    for i=1:N
        FU=Xnew(i,:)>ub;FL=Xnew(i,:)<lb;Xnew(i,:)=(Xnew(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        tempFitness = fun_calcobjfunc(objfunc, Xnew(i,:),Vio);
        if tempFitness<fitness(i)
            fitness(i)=tempFitness;
            X(i,:)= Xnew(i,:);
        end
        if tempFitness<GYbest
            GYbest=tempFitness;
            Xprey=Xnew(i,:);
        end
    end

    %开发
    for i=1:N
        if a>=w  
            delta= mean(Xnew,1);
            Xnew(i,:)=C*delta.*(X(i,:)-Xprey)+X(i,:);  %探索
       else
           
           if rand>=0.5
              h=randn(1,dim);%理论上最大值是正无穷，它产生的数均值是0，标准差是1，但实际上你不会看到很大的数产生（概率很小，几乎为0） 
           else
              h=randn(1,1);%理论上最大值是正无穷，它产生的数均值是0，标准差是1，但实际上你不会看到很大的数产生（概率很小，几乎为0） 
              % randn(m,n) 或 Y = randn([m n])     randn（random normal distribution）是一种产生标准正态分布的随机数或矩阵的函数
           end
           r1=rand; 
           Xnew(i,:)= Xprey-(Xprey*(2*r1-1)-X(i,:)*(2*r1-1)).*(gamma*h);  %2*r1-1=(-1,1)    Beta*h加强正态分布           
        end
      
    end
    for i=1:N
        FU=Xnew(i,:)>ub;FL=Xnew(i,:)<lb;Xnew(i,:)=(Xnew(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        tempFitness = fun_calcobjfunc(objfunc, Xnew(i,:),Vio);
        if tempFitness<fitness(i)
            fitness(i)=tempFitness;
            X(i,:)= Xnew(i,:);
        end 
        if tempFitness<GYbest
            GYbest=tempFitness;
            Xprey=Xnew(i,:);
        end
    end
    
%     [~, SortOrder]=sort(fitness);
%     fitness=fitness(SortOrder);
%     PopPos_=X;
%     for i=1:N
%         X(i,:)=PopPos_(SortOrder(i));
%     end

    CNVG(t)=GYbest;
   
end
Food_Score = GYbest;
end


function Y = fun_calcobjfunc(func, X,Vio)
N = size(X,1);
for i = 1:N
    Y(i) = CostFunction(X(i,:),Vio,func); 
end
end
function I=Intensity(N,Xprey,X)
for i=1:N-1
    di(i) =( norm((X(i,:)-Xprey+eps))).^2;
    S(i)=( norm((X(i,:)-X(i+1,:)+eps))).^2;
end
di(N)=( norm((X(N,:)-Xprey+eps))).^2;
S(N)=( norm((X(N,:)-X(1,:)+eps))).^2;
for i=1:N
    r2=rand;
    I(i)=r2*S(i)/(4*pi*di(i));
end
end
function [X]=initialization(N,dim,up,down)
if size(up,2)==1
    X=rand(N,dim).*(up-down)+down;
end
if size(up,2)>1
    for i=1:dim
        high=up(i);low=down(i);
        X(:,i)=rand(N,1).*(high-low)+low;
    end
end
end
