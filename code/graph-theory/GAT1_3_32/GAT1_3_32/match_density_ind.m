function G2thr = match_density_ind(G2,TD)%
    
% G2: correlation matrix with zero diagonal and positive correlation values
% TD: Target density of ineterest (number of edges in the original thresholded network)
% TD can be calculated using "nnz(triu(Graph))"

% G2thr : G2 thresholded at the target density TD

%% -- Hadi Hosseini, May 2012


%%
Step = 0.0001;%0.00001;%step for changing the threshold

%finding the maximum threshold

G2max = max(max(G2));
ss = -1;
ii = 0;

while ii <= G2max 

    ss = ss + 1;
    ii = Step + ss.*Step; %thr

    G2(G2<=ii) = 0;


    NN = nnz(triu(G2));%number of edges    

    %old delta
    if ss ~= 0

        delta_o = delta;

    else

        delta_o = NN - TD;

    end

    %current delta
    delta = NN - TD;


    %check the transition
    if delta <= 0

        if abs(delta) <= delta_o

            Thresh = ii;

            %exit
            ii = G2max + .1;

        elseif abs(delta) > delta_o

            Thresh = ii- Step;

            %exit
            ii = G2max + .1;

        end

    end

end

G2thr = G2; G2thr(G2thr <= Thresh) = 0;

clear G2 Thresh_Dens Ind

end