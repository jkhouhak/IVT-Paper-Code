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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
xdata = [15 30 45 60 90 120 150];
ydataAve1 =[23.0278,17815.31667,61178.38333,88255.15,98001.66667,102727.8333,92647.28333];
ydataAve2 =[20.0278,11197.47667,45311.63333,69455.95,81043.03333,82493.28333,74704.38333];
ydataAve3 =[20.944466,24768.45,165426.3333,307826.6667,455530.5,479963.6667,467688.2];
ydataAve4 =[20.8889,30144.08333,185069,312446.8333,436292.1667,467082.6667,449602.3333];
%%%%%%%%%%%%
%10ng-15min
ydata1 = [23.5	15.8333	28.8334	22.8334	26.5	20.6667
16574.5	12027.6	20973.5	19148.6	21578	16589.7
61548.9	54463.3	64487.7	60215.9	64961.9	61392.6
93896.8	91571.3	87061.7	88091.1	83214.2	85695.8
97224.1	102795	89951.4	95275.4	97571.1	105193
107669	102970	102022	101254	100313	102139
89512.8	91885.1	90544.7	92894.4	92485.7	98561]';

%%%%%%%%%%%%%%%%
%10ng-30min
ydata2 = [19	19.8334	18.3333	19.1667	19.8334	24
8637.13	7606.93	14453.9	14391	11455.3	10640.6
41315.9	40946.5	50220.2	49417.9	44789.1	45180.2
68664.5	66751.3	75101.1	70613	67099.4	68506.4
75237.8	76095.3	83009.9	89275.8	79021.6	83617.8
75011.4	84466.1	84414.1	86206.3	81099.8	83762
66562.6	68450.6	78669.9	79845.4	75101.1	79596.7]';

%%%%%%%%%%%%%%%%
%50ng-15min
ydata3 = [19	18.8333	20	23.8334	17.6667	26.3334
15970.9	18499.5	33408	32231.1	22802.3	25698.9
142599	152275	182372	183903	162337	169072
301707	303277	321027	314281	299850	306818
429097	405042	483029	489026	466595	460394
492790	423193	499328	495215	479304	489952
439205	450000	473616	495712	450809	479099]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%N/A%%%%%%%%%%%%% : (7,2)
%%%%%%%%%%%%%%%%%%%%
%50ng-30min
ydata4 = [15	18.5	19.6667	32	24.5	15.6667
20981.6	28133.7	34126	42306.9	23743.8	31572.5
160201	170257	198016	212073	175334	194533
307728	309891	309099	320918	300063	326982
375655	386793	475744	476374	443155	460032
452199	449587	470616	499573	467201	463320
421391	410809	477746	473151	458331	456186]';

err10ng15min=[];
err10ng30min=[];
err50ng15min=[];
err50ng30min=[];
for i=1:7
    err10ng15min(i) = max(ydata1(:,i) - ydataAve1(i));
    err10ng30min(i) = max(ydata2(:,i) - ydataAve2(i));
    err50ng15min(i) = max(ydata3(:,i) - ydataAve3(i));
    err50ng30min(i) = max(ydata4(:,i) - ydataAve4(i));
end
LW=1.5;
hold on

errorbar(xdata,ydataAve1,err10ng15min,'--o','color',[0.06,0.14,0.92],'LineWidth',LW)
errorbar(xdata,ydataAve2,err10ng30min,'--o','color',[0.00,0.09,0.61],'LineWidth',LW)
errorbar(xdata,ydataAve3,err50ng15min,'--o','color',[0.97,0.35,0.35],'LineWidth',LW)
errorbar(xdata,ydataAve4,err50ng30min,'--o','color',[0.66,0.06,0.06],'LineWidth',LW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
a=-2; b=0.5; c=3;

% 
ydata1 =[17815.31667,61178.38333,88255.15,98001.66667,102727.8333,92647.28333];
ydata2 =[11197.47667,45311.63333,69455.95,81043.03333,82493.28333,74704.38333];
ydata3 =[24768.45,165426.3333,307826.6667,455530.5,479963.6667,467688.2];
ydata4 =[30144.08333,185069,312446.8333,436292.1667,467082.6667,449602.3333];
%%%%%Calcuating heads and tails
xdata = [30 45 60 90 120 150];
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
xdata = 0:150;
plot(xdata, a1*b1.^xdata + c1,'color',[0.06,0.14,0.92],'linewidth',1.5);
hold on
plot(xdata, a2*b2.^xdata + c2,'color',[0.00,0.09,0.61],'linewidth',1.5);
plot(xdata, a3*b3.^xdata + c3,'color',[0.97,0.35,0.35],'linewidth',1.5);
plot(xdata, a4*b4.^xdata + c4,'color',[0.66,0.06,0.06],'linewidth',1.5);


xlabel('time (min)')
ylabel('Protein Level (a.u.)')
title('P.Falciparum CFPSS')

legend('10ng-15min','10ng-30min','50ng-15min','50ng-30min')
ylim([0 inf])
xlim([0 150])
% legend('9.2ng,15min-sim','54.8ng,15min-sim','9.2ng,30min-sim','54.8ng,30min-sim','9.2ng,60min-sim','54.8ng,60min-sim')
 hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
