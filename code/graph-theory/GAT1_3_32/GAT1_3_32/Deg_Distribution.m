function Deg_Distribution

input_type=input('type 1 to use previous GAT output or 2 for manual degree input (a mat file named "deg": ');
data_log = spm_select(1,'mat4GAT*','Select mat4GAT .mat (output from Load Data)');

ff = load(data_log,'mat4GAT');
Group1 = ff.mat4GAT.g1;Group2 = ff.mat4GAT.g2;
switch input_type

    case 1
    
        DminOrAuc=input('degree distribution at Dmin(type 1) or based on AUC (type 2): ');
        
        if isequal(DminOrAuc,1)
        
            D = spm_select(1,'KAN_FixedDens*','Select the KAN_FixedDens_Results.mat  ');
            
            if isempty(D)
            
                error('no input identified')
            
            end
            f1=load(D,'NetMes_Bin1');data1=f1.NetMes_Bin1;deg1=data1{1,1};
            f2=load(D,'NetMes_Bin2');data2=f2.NetMes_Bin2;deg2=data2{1,1};
            
        elseif isequal(DminOrAuc,2)%AUC
        
            D1 = spm_select(1,'AUC_NetMesReg_*','Select Group1 AUC_NetMesReg_Group1.mat  ');
            D2 = spm_select(1,'AUC_NetMesReg_*','Select Group2 AUC_NetMesReg_Group2.mat  ');
            
            if isempty(D1) || isempty(D2)
            
                error('no input identified')
            
            end
            
            f1=load(D1,'auc_MDeg1o');data1=f1.auc_MDeg1o;deg1=data1;
            f2=load(D2,'auc_MDeg2o');data2=f2.auc_MDeg2o;deg2=data2;
            
        
        end
        
    case 2
    
        D1 = spm_select(1,'.mat','Select 1st group "deg" .mat ');
        D2 = spm_select(1,'.mat','Select 2nd group "deg" .mat ');
        
        if isempty(D1) || isempty(D2)
        
            error('no input identified')
        
        end
        
        f1=load(D1,'deg');data1=f1.deg;deg1=data1;
        f2=load(D2,'deg');data2=f2.deg;deg2=data2;
        
end
    
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

f=fittype('x.^(a-1).*exp(-x./b)','independent','x');


[fitfun,GOF]=fit(xout',p_cum',f);
params=coeffvalues(fitfun);
a=params(1)
b=params(2)
b_log=log(b)
r2=GOF.rsquare

fprintf('%s\n',['fitted function:  x.^(' num2str(a) '-1).*exp(-x./' num2str(b) ')']);
fitted_f=xout.^(a-1).*exp(-xout./b);


plot(log(xout),log(fitted_f),'r.')

hold off
xlabel('log(degree)','fontsize',12,'fontweight','b')
ylabel('log(cumulative distribution)','fontsize',12)
title('Log Plot of Cumulative Degree Distribution','fontsize',14,'fontweight','b','fontangle','italic')

legend(Group1,'fit')
grid on

for i=1:10

    Finish=input('type 1 if done;2 for try again: ');
    
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

f=fittype('x.^(a-1).*exp(-x./b)','independent','x');


[fitfun,GOF]=fit(xout',p_cum',f);
params=coeffvalues(fitfun);
a=params(1)
b=params(2)
b_log=log(b)
r2=GOF.rsquare

fprintf('%s\n',['fitted function:  x.^(' num2str(a) '-1).*exp(-x./' num2str(b) ')']);
fitted_f=xout.^(a-1).*exp(-xout./b);


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
    fitted_f=xout.^(a-1).*exp(-xout./bb);
    
    
    plot(log(xout),log(fitted_f),'r.');
    hold off
    xlabel('log(degree)','fontsize',12,'fontweight','b')
    ylabel('log(cumulative distribution)','fontsize',12,'fontweight','b')
    title('Log Plot of Cumulative Degree Distribution','fontsize',14,'fontweight','b','fontangle','italic')
    
    legend(Group2,'fit')

end

hgsave(['DegDist_LogPlotFit' Group2 '.fig'])
hold off

