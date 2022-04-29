%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Khanh Dao Duc, Jalal Khouhak
%         Department of Mathematics 
%         The University of British Columbia (UBC)
% 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
%%%toleran

ce 
tol = 0.95;


%%%%%%%%%%%%%%%%%%DATA%%%%%%%%%%%%%%%%%%%

%Simulate the data%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
xdata = [15,30,45,60,90,120,150];
ydataAve1 = [18.58333333,6092.51,31081.61667,47467.83333,58311.1,61394.15,57774.4];
ydataAve2 = [37.30566667,11866.83333,80876.33333,142514.1667,193246.5,194404.5,192915];
ydataAve3 = [19.54631113,24289.32777,155640.111,281197.889,402180.1113,429591.5557,415283.122];
ydataAve4 = [19.11111667,25277.35,199506.8333,424508.8333,709548.6667,780483,789057.5];

%%%%%%%%%%%%
ydata1 = [19.5	20.5	17.5	14.8333	20.5	18.6667
6685.04	7188.84	4434.16	4558.04	7641.95	6047.03
38861.8	35559.2	28522.6	26317.6	28554	28674.5
53362.7	51896.2	45632.1	44991.9	43754.9	45169.2
62811.6	62486.1	53241.6	55322.7	57867	58137.6
63888.3	65765.2	64098.3	63776	54923.1	55914
59556.9	61613	58131.9	57541.7	54367	55435.9]';

%%%%%%%%%%%%%%%%
ydata2 = [21	17.6667	15.3333	20.6667	129.834	19.3333
12826.4	11924.5	13693.2	12149	10360.2	10247.7
85822.6	89571	92232.7	88177.2	64321.8	65132.7
145388	145026	162408	154715	122454	125094
197010	208596	219808	215480	156045	162540
204807	208412	219180	229713	150477	153838
204091	205536	214681	231244	150978	150960]';

%%%%%%%%%%%%%%%%
ydata3 = [16.6667	17.1667	12	18	17.6667	19.3333
17817.2	26899.9	15167.2	18545.8	15583	13719.6
135231	135742	106897	103787	108732	108161
243411	238922	198748	216575	224383	217882
304228	316230	311968	329596	314397	311887
321735	329151	365910	366008	324218	343348
324726	324345	351367	362198	288361	320356]';
%%%%%%%%%%%%%%%%%%%%
ydata4 = [19.6667	19.8334	17.3333	16.3333	22.1667	19.3333
22318.4	26434.5	23221.1	27580.6	24979.8	27129.7
181148	211463	193536	220206	192070	198618
390575	419012	413939	468424	418685	436418
678472	748700	679609	748913	706577	695021
688717	798826	752829	842288	765290	834948
718506	791679	765440	862356	770266	826098]';
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

xdata = [45 60 90 120 150];
ydata1 = [31081.61667,47467.83333,58311.1,61394.15,57774.4];
ydata2 = [80876.33333,142514.1667,193246.5,194404.5,192915];
ydata3 = [155640.111,281197.889,402180.1113,429591.5557,415283.122];
ydata4 = [199506.8333,424508.8333,709548.6667,780483,789057.5];

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
xdata = 0:150;
plot(xdata, a1*b1.^xdata + c1,'color',[0.40,0.51,0.94],'linewidth',1.5);
hold on
plot(xdata, a2*b2.^xdata + c2,'color',[0.06,0.14,0.92],'linewidth',1.5);
plot(xdata, a3*b3.^xdata + c3,'color',[0.00,0.09,0.61],'linewidth',1.5);
plot(xdata, a4*b4.^xdata + c4,'color',[0.04,0.03,0.12],'linewidth',1.5);

xlabel('time (min)')
ylabel('Protein Level (a.u.)')
title('P.Falciparum CFPSS')

legend('10ng','25ng','50ng','100ng')
ylim([0 inf])
xlim([0 150])
% legend('9.2ng,15min-sim','54.8ng,15min-sim','9.2ng,30min-sim','54.8ng,30min-sim','9.2ng,60min-sim','54.8ng,60min-sim')
 hold on
%plot(xdata,ydata,'r+')
% tauu = -1/log(b);
% deltaT = tauu*log((-a/c));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
