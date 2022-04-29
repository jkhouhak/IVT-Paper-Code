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
xdata = [5 10 15 20 30 45 60];
ydataAve1 =[195026.8333,1038039.167,1648073.333,1953148.333,2027258.333,2154196.667,2106806.667];
ydataAve2 =[128436.5,674758.8333,964648.1667,1170175.333,1418876.667,1539318.333,1543413.333];
ydataAve3 =[1930060,7076921.667,8311923.333,9172461.667,8611016.667,7409325,6635273.333];
ydataAve4 =[1866635,5255980,6081458.333,6040480,6075821.667,5770991.667,5019161.667];
%%%%%%%%%%%%
%10ng-15min
ydata1 = [269143	218472	158067	140628	210639	173212
1231140	1147190	996803	991782	966192	895128
1811880	1893710	1796770	1704960	1354710	1326410
2209100	2317630	1979180	1975380	1663110	1574490
2334780	2401560	2073520	2157740	1624110	1571840
2566650	2640430	2114760	2275990	1718550	1608800
2535890	2538170	2286080	2180860	1564650	1535190]';

%%%%%%%%%%%%%%%%
%10ng-30min
ydata2 = [132493	121747	125531	130068	120125	140655
710799	716451	731696	763984	558370	567253
1083600	1031080	944108	1068320	837214	823567
1270910	1331700	1271770	1228830	981263	936579
1507200	1599760	1661720	1508060	1132320	1104200
1775190	1762770	1718350	1578160	1215370	1186070
1742770	1713240	1740570	1644030	1240550	1179320]';

%%%%%%%%%%%%%%%%
%50ng-15min
ydata3 = [2115730	2057720	1859600	1811210	1931490	1804610
7634680	7397830	7533660	7612220	6346600	5936540
9042940	9111220	9342410	9794160	6458650	6122160
9920580	9529310	1.01E+07	9973780	8126730	7357170
8854940	8887740	1.16E+07	1.19E+07	5491960	4939760
8987130	8852960	8347890	8227860	4607240	5432870
6638620	6321450	9472490	9329800	3906630	4142650]';
%%%%%%%%%%%%%%%%%%%%
%50ng-30min
ydata4 = [1863610	1887470	2003150	1779490	1912130	1753960
5810520	5782400	5848660	5611400	4385100	4097800
6832160	6811090	6696810	6868480	4533650	4746560
6948760	6864490	6504080	6371900	4727130	4826520
6634820	6643290	7166010	7404100	4279850	4326860
6712140	7023020	6657940	6439110	3991070	3802670
6079240	5914540	5474020	5232040	4085500	3329630]';
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
ydata1 =[195026.8333,1038039.167,1648073.333,1953148.333,2027258.333,2154196.667,2106806.667];
ydata2 =[128436.5,674758.8333,964648.1667,1170175.333,1418876.667,1539318.333,1543413.333];
ydata3 =[1930060,7076921.667,8311923.333,9172461.667,8611016.667,7409325,6635273.333];
ydata4 =[1866635,5255980,6081458.333,6040480,6075821.667,5770991.667,5019161.667];

%%%%%Calcuating heads and tails
% xdata = [45 60 90 120 150];
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
xlim([0 60])
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
