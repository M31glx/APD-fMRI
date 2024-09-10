
data_log = spm_select(1,'mat4GATd*','Select mat4GATd .mat ');


ff = load(data_log,'mat4GATd');
Group1 = ff.mat4GATd.g1;Group2 = ff.mat4GATd.g2;

MinMesPlot=ff.mat4GATd.MinThr;
MaxMesPlot=ff.mat4GATd.MaxThr;
MesStepPlot=ff.mat4GATd.MesStep;

D1 = spm_select(1,['NetMesBin_D_' Group1 '_FinalThrRange'],['Select ' Group1 ' Unadjusted NetMesBin_D_*.mat file that you want to examine']);
D2 = spm_select(1,['NetMesBin_D_' Group2 '_FinalThrRange'],['Select ' Group2 ' Unadjusted NetMesBin_D_*.mat file that you want to examine']);

if isempty(D1) || isempty(D2)

    error('no input identified')
    
end

f1=load(D1,'NetMes_Bin');data1=f1.NetMes_Bin;
f2=load(D2,'NetMes_Bin');data2=f2.NetMes_Bin;


%% 

xxx = [MinMesPlot:MesStepPlot:MaxMesPlot];
MinThr=input('minimum threshold:');MinIdx=find(single(xxx)==single(MinThr));
MaxThr=input('maximum threshold:');MaxIdx=find(single(xxx)==single(MaxThr));

if MaxIdx > size(data1,2)
    
    MaxIdx = size(data1,2);

end


NetMes1=data1(:,MinIdx:MaxIdx);
NetMes2=data2(:,MinIdx:MaxIdx);


deg1=[];

for i=1:size(NetMes1,1)
       
    temp_deg1=[];
    
    for j=1:size(NetMes1,2)
        
        temp_deg1=[temp_deg1;NetMes1{i,j}{1,1}];
        
    end
    
    deg1=[deg1;mean(temp_deg1)];
    

end

deg1 = mean(deg1);



deg2=[];

for i=1:size(NetMes2,1)
       
    temp_deg2=[];
    
    for j=1:size(NetMes2,2)
        
        temp_deg2=[temp_deg2;NetMes2{i,j}{1,1}];
        
    end
    
    deg2=[deg2;mean(temp_deg2)];
    

end

deg2 = mean(deg2);
    
    
%% 

Step=1;

%%
deg=deg1;
MIN=min(deg);MAX=max(deg);X=MIN:Step:MAX;
[n,xout] = hist(deg,X);
p=n/sum(n);
p_cum=[];
count=0;

for i=xout

    count=count+1;
    p_cum=[p_cum 1-sum(p(1:count-1))];
    
end

figure;HistStep=round((MAX-MIN)./5);
hist(deg,MIN:HistStep:MAX);
xlabel('degree','fontsize',12,'fontweight','b')
ylabel('number of nodes','fontsize',12,'fontweight','b')
title('Degree Distribution','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1)
grid on
hgsave(['DegHist_' Group1 '.fig'])

figure;
plot(xout,p_cum,'+');grid on
xlabel('degree','fontsize',12,'fontweight','b')

title('Cumulative Degree Distribution','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group1)
hgsave(['DegDistCum_NormalPlot_' Group1 '.fig'])

figure
plot(log(xout),log(p_cum),'+')
hold


f=fittype('x.^(-a).*exp(-x./b)','independent','x');

[fitfun,GOF]=fit(xout',p_cum',f);
params=coeffvalues(fitfun);
a=params(1)
b=params(2)
b_log=log(b)
r2=GOF.rsquare
fprintf('%s\n',['fitted function:  x.^(' num2str(a) '-1).*exp(-x./' num2str(b) ')']);

fitted_f=xout.^(-a).*exp(-xout./b);

plot(log(xout),log(fitted_f),'r.')

hold off
xlabel('log(degree)','fontsize',12,'fontweight','b')
ylabel('log(cumulative distribution)','fontsize',12)
title('Log Plot of Cumulative Degree Distribution','fontsize',14,'fontweight','b','fontangle','italic')

legend(Group1,'fit')
grid on

for i=1:10

    Finish=input('press 1 if done;2 for try again: ');
    
    if isequal(Finish,1)
    
        break
    
    end
    
    bb=input('type your choice of b: ');
    figure
    plot(log(xout),log(p_cum),'+')
    hold
    
    fitted_f=xout.^(a-1).*exp(-xout./bb);
    
    plot(log(xout),log(fitted_f),'r.')
    xlabel('log(degree)','fontsize',12,'fontweight','b')
    ylabel('log(cumulative distribution)','fontsize',12)
    title('Log Plot of Cumulative Degree Distribution','fontsize',14,'fontweight','b','fontangle','italic')
    
    legend(Group1,'fit')
    grid on

end

hgsave(['DegDist_LogPlotFit' Group1 '.fig'])
hold off 

%% 

deg=deg2;
MIN=min(deg);MAX=max(deg);X=MIN:Step:MAX;
[n,xout] = hist(deg,X);
p=n/sum(n);
p_cum=[];
count=0;

for i=xout
    
    count=count+1;
    p_cum=[p_cum 1-sum(p(1:count-1))];
    
end

figure;HistStep=round((MAX-MIN)./5);
hist(deg,MIN:HistStep:MAX);
xlabel('degree','fontsize',12,'fontweight','b')
ylabel('number of nodes','fontsize',12,'fontweight','b')
title('Degree Distribution','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group2)
grid on
hgsave(['DegHist_' Group2 '.fig'])

figure;
plot(xout,p_cum,'+');grid on
xlabel('degree','fontsize',12,'fontweight','b')
title('Cumulative Degree Distribution','fontsize',14,'fontweight','b','fontangle','italic')
legend(Group2)
hgsave(['DegDistCum_NormalPlot_' Group2 '.fig'])

figure
plot(log(xout),log(p_cum),'+')
hold

f=fittype('x.^(-a).*exp(-x./b)','independent','x');

[fitfun,GOF]=fit(xout',p_cum',f);
params=coeffvalues(fitfun);
a=params(1)
b=params(2)
b_log=log(b)
r2=GOF.rsquare
fprintf('%s\n',['fitted function:  x.^(' num2str(a) '-1).*exp(-x./' num2str(b) ')']);

fitted_f=xout.^(-a).*exp(-xout./b);

plot(log(xout),log(fitted_f),'r.')
hold off
xlabel('log(degree)','fontsize',12,'fontweight','b')
ylabel('log(cumulative distribution)','fontsize',12,'fontweight','b')
title('Log Plot of Cumulative Degree Distribution','fontsize',14,'fontweight','b','fontangle','italic')

legend(Group2,'fit')
grid on

for i=1:10

    Finish=input('type 1 if done, 2 for try again: ');
    
    if isequal(Finish,1)
    
        break

    end
    
    bb=input('type your choice of b: ');
    figure
    plot(log(xout),log(p_cum),'+')
    hold
    
    fitted_f=xout.^(-a).*exp(-xout./bb);
    
    plot(log(xout),log(fitted_f),'r.');
    hold off
    xlabel('log(degree)','fontsize',12,'fontweight','b')
    ylabel('log(cumulative distribution)','fontsize',12,'fontweight','b')
    title('Log Plot of Cumulative Degree Distribution','fontsize',14,'fontweight','b','fontangle','italic')
    
    legend(Group2,'fit')

end

hgsave(['DegDist_LogPlotFit' Group2 '.fig'])
hold off

fprintf('%-4s\n','.... done ....');


