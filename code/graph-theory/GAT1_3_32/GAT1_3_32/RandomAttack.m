function MeanSize=RandomAttack(G,nsim)
%This function do random attack analysis on the input binary graph
%The inputs are:
%G(binary graph),
%nsim is the number of simuations 

%% Hadi Hosseini, June 2011

%%
z=size(G,1);
sz_full=size(largest_component(G),2);%largest connected component
msz=[];%mean size for each i
for i=1:z
    sz=[];
    for j=1:nsim%number of simulations for each iteration 
        GG=G;%a copy of the original graph
        RandInx=randperm(z);%we randomize say 90 and select the first i variables
        Randmat=RandInx(1:i);
        GG(:,Randmat)=0;GG(Randmat,:)=0;
        %calculate the size of the largest remained connected component
        sz=[sz;size(largest_component(GG),2)];
    end
    msz=[msz;mean(sz)];
end
MeanSize=msz./sz_full;      
    
    
    
    
