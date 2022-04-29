%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Khanh Dao Duc, Jalal Khouhak
%         Department of Mathematics 
%         The University of British Columbia (UBC)
% 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
% close all
%%%tolerance 
tol = 0.95;


%%%%%%%%%%%%%%%%%%DATA%%%%%%%%%%%%%%%%%%%

%Simulate some data%%%%%%%%%%%%%%%%
% xdata=[5,10,15,20,30,45,60];
% xdata=[15,20,30,45,60];
% xdata = [30,45 60 90 120 150];
% xdata = [45 60 90 120 150];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
xdata = [5,10,15,20,30,45,60];
ydataAve1 =[134814.2,985829.5,2096233.333,2372001.66,2586933.33,2641990,2803356.667];
ydataAve2 =[535553.5,3611828.333,6740431.667,7418885,8796735,7746235,8146876.667];
ydataAve3 =[1117699.167,5887583.333,10085211.67,10977265,1.17E+07,1.07E+07,10604476.67];
ydataAve4 =[2765650,11820630,1.62E+07,16983333.33,1.59E+07,1.42E+07,1.60E+07];
%%%%%%%%%%%%
ydata1 = [122274	48941.2	198857	118714	188295	131804
881883	687200	1003220	800384	1364390	1177900
1633750	1715720	2302490	1849480	2617000	2458960
1861590	2014680	2224520	1766310	3234960	3129950
1795750	2114770	2072660	2433020	3644590	3460810
2176520	2086970	2717050	1729960	3552270	3589170
2591690	2440720	2123590	2557590	3531330	3575220]';

%%%%%%%%%%%%%%%%
ydata2 = [332140	237628	898097	612235	547600	585621
2400170	2546680	4899220	3702760	4069870	4052270
5589020	5744740	7655200	6992130	7473480	6988020
5967130	7040710	4940980	9270310	8614570	8679610
6720410	6519710	9697930	1.18E+07	9288030	8754330
5238980	6953060	7195720	1.08E+07	8114260	8175390
5034900	6899380	8533690	8542230	9963730	9907330]';

%%%%%%%%%%%%%%%%
ydata3 = [702709	729238	1255970	575898	1596410	1845970
5075330	5439450	4501700	4397790	7882040	8029190
9607690	1.02E+07	9652550	7866430	1.18E+07	1.14E+07
9717650	1.05E+07	1.08E+07	8645940	1.41E+07	1.21E+07
1.16E+07	1.09E+07	1.22E+07	1.01E+07	1.27E+07	1.30E+07
1.27E+07	1.04E+07	9714260	9532970	1.12E+07	1.05E+07
9455460	1.08E+07	6433200	1.10E+07	1.39E+07	1.20E+07]';
%%%%%%%%%%%%%%%%%%%%
ydata4 = [1990540	2670260	3261810	2758850	3023790	2888650
9803880	1.04E+07	1.34E+07	1.13E+07	1.28E+07	1.32E+07
1.86E+07	1.67E+07	1.41E+07	1.67E+07	1.62E+07	1.50E+07
1.52E+07	1.58E+07	1.42E+07	1.63E+07	1.71E+07	2.33E+07
1.66E+07	1.64E+07	1.23E+07	1.37E+07	1.71E+07	1.92E+07
1.72E+07	1.39E+07	1.01E+07	1.55E+07	1.35E+07	1.51E+07
1.83E+07	2.10E+07	1.31E+07	1.27E+07	1.60E+07	1.51E+07]';
err10ng=[];
err25ng=[];
err50ng=[];
err100ng=[];
for i=1:7
    err10ng(i) = max(ydata1(:,i) - ydataAve1(i));
    err25ng(i) = max(ydata2(:,i) - ydataAve2(i));
    err50ng(i) = max(ydata3(:,i) - ydataAve3(i));
    err100ng(i) = max(ydata4(:,i) - ydataAve4(i));
end
LW=1.5;
hold on
errorbar(xdata,ydataAve1,err10ng,'--o','color',[0.40,0.51,0.94],'LineWidth',LW)
errorbar(xdata,ydataAve2,err25ng,'--o','color',[0.06,0.14,0.92],'LineWidth',LW)
errorbar(xdata,ydataAve3,err50ng,'--o','color',[0.00,0.09,0.61],'LineWidth',LW)
errorbar(xdata,ydataAve4,err100ng,'--o','color',[0.04,0.03,0.12],'LineWidth',LW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
a=-2; b=0.5; c=3;

% color1=[0.40,0.51,0.94];
% color2=[0.06,0.14,0.92];
% color3 = [0.00,0.09,0.61];
% color4=[0.04,0.03,0.12];
% 
ydata4 = [2765650,11820630,1.62E+07,16983333.33,1.59E+07,1.42E+07,1.60E+07];
ydata3 = [1117699.167,5887583.333,10085211.67,10977265,1.17E+07,1.07E+07,10604476.67];
ydata2=[535553.5,3611828.333,6740431.667,7418885,8796735,7746235,8146876.667];
ydata1 = [134814.2,985829.5,2096233.333,2372001.66,2586933.33,2641990,2803356.667];
% 
% ydata = ydata4;
% color = color4;
% plot(xdata,ydata,'--o','color',color,'MarkerEdgeColor',color,'MarkerFaceColor',color)
% legend('10ng-6min','10ng-12min','50ng-6min','50-12min')
% hold on 
% legend('10ng-6min','10ng-12min','50ng-6min','59ng-12min')
% y=[36.91671667,35.4445,36.80561667,39.38896667]
% x=[10 10 10 10]
% plot(x,y,'ko','MarkerEdgeColor','k','MarkerFaceColor','k')


%%%%%Calculating heads and tails
mxdata = breakdownvec(xdata);
mydata1 = breakdownvec(ydata1);
mydata2 = breakdownvec(ydata2);
mydata3 = breakdownvec(ydata3);
mydata4 = breakdownvec(ydata4);
%%%%Finding where the delay ends - Calculating tails
for i = 1:length(mxdata)
    [Rvalues(i),a1,b1,c1,delta] = fitfun(mxdata{i},mydata1{i});
    if Rvalues(i) >= tol
        break
    end
end
for i = 1:length(mxdata)
    [Rvalues(i),a2,b2,c2,delta] = fitfun(mxdata{i},mydata2{i});
    if Rvalues(i) >= tol
        break
    end
end
for i = 1:length(mxdata)
    [Rvalues(i),a3,b3,c3,delta] = fitfun(mxdata{i},mydata3{i});
    if Rvalues(i) >= tol
        break
    end
end
for i = 1:length(mxdata)
    [Rvalues(i),a4,b4,c4,delta] = fitfun(mxdata{i},mydata4{i});
    if Rvalues(i) >= tol
        break
    end
end
[a1,b1,c1];
[a2,b2,c2];
[a3,b3,c3];
[a4,b4,c4];

figure(1);

syms x
xdata = 0:60;
plot(xdata, a1*b1.^xdata + c1,'color',[0.40,0.51,0.94],'linewidth',1.5);
hold on
plot(xdata, a2*b2.^xdata + c2,'color',[0.06,0.14,0.92],'linewidth',1.5);
plot(xdata, a3*b3.^xdata + c3,'color',[0.00,0.09,0.61],'linewidth',1.5);
plot(xdata, a4*b4.^xdata + c4,'color',[0.04,0.03,0.12],'linewidth',1.5);


% axes('position',[0.25 0.25 0.5 0.5]);
% box on
% plot(xdata,-(6.2982455e+04)*0.908020.^xdata + 1.86719517e+04)
% % plot(xdata,(-2.12471e+04)*0.97424.^xdata + 1.5193e+04,'linewidth',1.5)
% axis tight
% hold on		
xlabel('time (min)')
ylabel('Protein Level (a.u.)')
title('Rabbit CFPSS')
% legend('50ng-6min','50ng-12min')
% legend('10ng,15min-fit','10ng,30min-fit','50ng,15min-fit','50ng,30min-fit')
% legend('10ng-fit','25ng-fit','50ng-fit','100ng-fit')
legend('10ng','25ng','50ng','100ng')
ylim([0 inf])
xlim([0 60])
% legend('9.2ng,15min-sim','54.8ng,15min-sim','9.2ng,30min-sim','54.8ng,30min-sim','9.2ng,60min-sim','54.8ng,60min-sim')
 hold on
%plot(xdata,ydata,'r+')
% tauu = -1/log(b);
% deltaT = tauu*log((-a/c));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ydata1 = [42.77783333;40.00006667;3376.965;10645.95667;14418.48333;17159.78333;19304.88333]';
% ydata2 = [40.80563333;37.77786667;547.2648333;3146.668333;5371.008333;8050.38;11157.77833]';
% ydata3 = [42.8334;45.77788333;78194.06667;413853.5;633987.1667;696276.5;718321.1667]';
% ydata4 = [39.30561667;70.91698333;48271.33333;170006.3333;274315.1667;324035.5;343845.3333]';

% syms x
% xdata = 0:60;
% axes('position',[0.25 0.25 0.5 0.5]);
% box on
% plot(xdata(xdata>=15),-(6.2982455e+04)*0.908020.^xdata(xdata>=15) + 1.86719517e+04,xdata(xdata>=15),(-2.12471e+04)*0.97424.^xdata(xdata>=15) + 1.5193e+04,'linewidth',1.5)
% xlim([0 60])
% ylim([0 inf])
% hold on
% xdata=[5,10,15,20,30,45,60];
% plot(xdata,ydata1,'--o',xdata,ydata2,'--o')
% xlim([0 60])
% ylim([0 inf])
% axis tight
		
		


%plot(xl, alpha + beta.*xl); %this is plotting the linearfit for the delay - ignore if you prefer. ​
%%%Functions%%%%
function mdata = breakdownvec(wdata)
    for i = 1:length(wdata)-1
        mdata{i} = (wdata(i:end));
    end
end
function [Rsq,a,b,c,delta] = fitfun(xdata,ydata)
    %Finding relevant values:
    n = length(xdata);
    xmax = max(xdata);
    xmin = min(xdata);
    %Compute Si
    %S = zeros(1,n);
    %S(2:end) = S(1:end-1)+(((ydata(1:end-1) + ydata(2:end))/2).*(xdata(2:end) - xdata(1:end-1)));
    
    %%%%% For comparison, I rewrote the formula for S from https://math.stackexchange.com/questions/1163618/exponential-curve-fit
    S=zeros(1,n);S(1)=0;
    for i=2:length(S)
    S(i)=S(i-1)+0.5*(ydata(i)+ydata(i-1))*(xdata(i)-xdata(i-1));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %S
    %S2  %%% this displays different values from S 
    %Solve the linear system to find b 
    M = zeros(2,2);
    M(1,1) = sum((xdata(1:end) - xdata(1)).^2);
    M(1,2) = sum(((xdata(1:end) - xdata(1))).*S); 
    M(2,1) = M(1,2);
    M(2,2) = sum(S.^2);
   
    V(1) = sum((xdata(1:end) - xdata(1)).*(ydata(1:end) - ydata(1)));
    V(2) = sum((ydata(1:end) - ydata(1)).*S);
    B = M\V';
    b = exp(B(2)); %% on the stackoverflow they say to do b = exp(B(2)), but it seems to fit better without the exp
    theta = b.^xdata;
    %Solve the linear system to find a,c
    CA = [n sum(theta); sum(theta) sum(theta.^2)]\[sum(ydata) ; sum(ydata.*theta)];
    c = CA(1);
    a = CA(2);
    delta = log(b);
    
    
    %Calculating R^2
    sqtot = sum((ydata+mean(ydata)).^2);
    yest = a*b.^(xdata) + c;
    %yest = a*exp(delta*(xdata)) + c; 
    sqres = sum((ydata-yest).^2);
    Rsq = 1 - (sqres/sqtot);
end
