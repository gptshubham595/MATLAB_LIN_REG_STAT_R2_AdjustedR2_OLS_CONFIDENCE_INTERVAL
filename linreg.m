%Model:
% y=a0+a1*x

%Least Square Method

%Sr= Σ e(i)^2 = Σ (y(i)-a0-a1*x(i))^2

%| n 	   Σx(i)   | |a0| = |   Σy(i) |
%| Σx(i)   Σx(i)^2 | |a1|   |Σx(i)y(i)|
%Aa=B
% a1= (n Σ(xi*yi) - Σxi*Σyi)/(n Σ(xi*xi) - Σxi*Σxi)
% a0 = mean(y) - a1*mean(x)

X=[0,2,4,6,9,11,12,15,17,19];
Y=[5,6,7,6,9,8,8,10,12,12];
n=length(X);

% or just 
a=polyfit(X,Y,1)
Ym=a(2)+a(1)*X;

%ERROR
e=Y-Ym
Sr=sum(e.^2)
%St=Deviation of y wrt mean y
St=sum((Y-mean(Y)).^2)

R2=(St-Sr)/St
K=1;
Adjusted_R2=1-((Sr/(n-K))/(St/(n-1)))

%PLOT
X1=min(X)-5:max(X)+5;
Ym=a(2)+a(1)*X1;
plot(X,Y,'bo') 
hold on,
plot(X1,Ym)
legend('data','linear')

lower_95=round(a.*1.05);
upper_95=round(a.*0.95);    
degree_fre=N-K-1;
alpha=0.05; % for 95%
Prob=alpha/2;
%Int    Std_Err  tStat Pvalue Low95 Up95
%a(2)
%a(1)
std_err=std(Y)/sqrt(N);
%[(100-95)/200 (95+(100-95/2))/100]
ts = tinv([0.025  0.975],N-1);      % T-Score
CI = mean(Y) + ts*std_err;