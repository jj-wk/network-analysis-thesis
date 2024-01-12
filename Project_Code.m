%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choose which altered edgelists to use %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PFa = 0.1;

extraEdges = false;

alteredWeights = false;

shuffledWeights = false;

airports = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating synthetic adjancency matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 1000; % Number of nodes

m0 = 6; % Number of nodes in intial core
m = 4; % Number of links formed by new nodes

% Initialising adjacency matrix of initial core (fully connected)
A = ones(m0) - eye(m0);
A = sparse(A);

for i = 1:N

    k = full(sum(A)); % Degree sequence
    p = k/sum(k); % Probability of linking to a node

    p = cumsum(p);

    r = rand(m,1);

    %%% Adding new nodes and links
    ind = [];

    for j = 1:m

        aux = p - r(j);
        aux = find(aux > 0);
        ind = [ind; aux(1)];

    end

    ind = unique(ind);

    A = [A; zeros(1,size(A,2))];
    A = [A zeros(size(A,1),1)];

    A(end,ind) = 1; %sets edge weights
    A(ind,end) = 1;

end


%finds length of A
[i,j,s] = find(A);

%produces a power law distribution from the normal distribution for both
%lists
pd = makedist('Normal','mu',250,'sigma',80); %heterogeneous
%pd = makedist('Normal','mu',250,'sigma',10); %homogeneous

%iterates through A's nonzero values to draw strengths from the normal
%distributions
%change the pd to change the power distribution
for k = 1:length(i)
    A(i(k),j(k)) = round(abs(random(pd)));
end


%saves network to txt file
[i,j,s] = find(A);
syntheticNetwork = [i j s];
syntheticNetwork2 = [i j];
writematrix(syntheticNetwork);
writematrix(syntheticNetwork2);


%testing the filter on the initial synthetic dataset
title("Filtering Techniques on Unchanged Synthetic Network");
%height(nonzeros(A))
filter(A);
global initialHyGe;
initialHyGe = 0;

%testing PF filter on initial synthetic dataset
PFResults = PF(A,PFa,10,0);
validatePF(PFResults);
global initialPF;
initialPF = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating altered edgelists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if extraEdges
    global validatedArrayHyGe;
    validatedArrayHyGe = [];
    global validatedArrayPF;
    validatedArrayPF = [];
    global intervals;
    intervals = [];
    
    for x=0:0.005:0.2 %remove normally
        for y=1:3 %remove normally, does 120 loops
            x
            alteredA1 = A;
            %for loop adding extra fake links to certain nodes
            selectedNodeProbability = 0.05;  %probability that new edges are made for a node 0.05
            edgeProbability = x;  %probability that new edges are created between two nodes 0.02
            for i = 1:height(alteredA1)
                x1 = rand;
                x2 = rand;
                if x1 < selectedNodeProbability
                    for j = 1:height(alteredA1)
                        if x2 < edgeProbability
                            alteredA1(i,j) = round(abs(random(pd)));
                            alteredA1(j,i) = round(abs(random(pd)));
                        end
                    end
                end
            end
            %testing the filter on this new altered matrix with fake links
            %figure();
            title("Filtering Techniques on Altered Synthetic Network");
            height(nonzeros(alteredA1))
            filter(alteredA1);
            PFResults = PF(alteredA1,PFa,10,0);
            validatePF(PFResults);
        
            %saves adjacency matrix as txt file
            [i,j,s] = find(alteredA1);
            alteredAdjacency = [i j];
            writematrix(alteredAdjacency);
    
            global intervals;
            intervals(end + 1) = x;
        end
    end
    
    figure();
    hold on
    plot(intervals,validatedArrayHyGe,'LineWidth',3);
    plot(intervals,validatedArrayPF,'LineWidth',3);
    title("Filtering Techniques on Unchanged Synthetic Network with Increasing Number of Edges");
    xlabel('Percentage used to add edges');
    ylabel('Percentage of edges validated');
end


if alteredWeights
    global validatedArrayHyGe;
    validatedArrayHyGe = [];
    global validatedArrayPF;
    validatedArrayPF = [];
    global intervals;
    intervals = [];

    for x=0:0.05:5 %remove normally
        for y=1:1 %remove normally, does 120 loops
            x
            alteredA2 = A;
            %for loop changing edge weights for certain nodes
            selectedNodeProbability = 0.05;  %probability that weights are altered for a node
            gamma = x; %multiple used to alter edge weights
            for i = 1:height(alteredA2)
                x1 = rand;
                if x1 < selectedNodeProbability
                    %alters out weights by 1+gamma
                    for j = 1:height(alteredA2)
                        if alteredA2(i,j) ~= 0
                            alteredA2(i,j) = round(alteredA2(i,j) * gamma);
                        end
                    end
                    %alters in weights by 1+gamma
                    for j = 1:height(alteredA2)
                        if alteredA2(j,i) ~= 0
                            alteredA2(j,i) = round(alteredA2(j,i) * gamma);
                        end
                    end
                end
            end
            %testing the filter on this new altered matrix changed edge weights
            %figure();
            title("Filtering Techniques on Altered Synthetic Network");
            height(nonzeros(alteredA2))
            filter(alteredA2);
            PFResults = PF(alteredA2,PFa,10,0);
            validatePF(PFResults);
        
            %saves adjacency matrix as txt file
            [i,j,s] = find(alteredA2);
            alteredAdjacency = [i j s];
            writematrix(alteredAdjacency);

            global intervals;
            intervals(end + 1) = x;
        end
    end

    figure();
    hold on
    plot(intervals,validatedArrayHyGe,'LineWidth',3);
    plot(intervals,validatedArrayPF,'LineWidth',3);
    title("Filtering Techniques on Unchanged Synthetic Network with Increasing Number of Edges");
    xlabel('Percentage used to add edges');
    ylabel('Percentage of edges validated');
end


if shuffledWeights
    global validatedArrayHyGe;
    validatedArrayHyGe = [];
    global validatedArrayPF;
    validatedArrayPF = [];

    for x=1:100
        alteredA3 = A;
        %finds length of A
        [i,j,s] = find(A);
        shuffleArray = [];
        %adds all current edge weights to an array
        for k = 1:length(i)
            shuffleArray(end + 1) = A(i(k),j(k));
        end
        %shuffles edge weight array and sets them as new weights
        shuffleArray = shuffleArray(randperm(length(shuffleArray)));
        [i,j,s] = find(alteredA3);
        for k = 1:length(i)
            alteredA3(i(k),j(k)) = shuffleArray(k);
        end
        %testing the filter on this new altered matrix with shuffled edge weights
        %figure();
        title("Filtering Techniques on Altered Synthetic Network");
        height(nonzeros(alteredA3))
        filter(alteredA3);
        PFResults = PF(alteredA3,PFa,10,0);
        validatePF(PFResults);
    
        %saves adjacency matrix as txt file
        [i,j,s] = find(alteredA3);
        alteredAdjacency = [i j];
        writematrix(alteredAdjacency);
    end

    validatedArrayHyGe = (validatedArrayHyGe - initialHyGe);
    validatedArrayPF = (validatedArrayPF - initialPF);

    figure();
    histogram(validatedArrayHyGe);
    title("Percentage differences of validated links after edge weight shuffling");
    xlabel("Difference in percentage of validated links");
    ylabel("Frequency");
    figure();
    histogram(validatedArrayPF);
    title("Percentage differences of validated links after edge weight shuffling");
    xlabel("Difference in percentage of validated links");
    ylabel("Frequency");
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Using airport network%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airportsPFtable = table();

if airports
    %chooses what file to open which is the Aiport Network and comes from the Bureau of
    %Transportation Statisitcs
    filename = 'AIRPORT_NETWORK.csv';
    
    %reads csv file data into matlab table
    data = readtable(filename);
    
    
    %produces table with unique originIDs and DestIDs
    uniqueOriginIDs = unique(data(:,"ORIGIN_AIRPORT_ID"));
    uniqueDestIDs = unique(data(:,"DEST_AIRPORT_ID"));
    %concatenates the two tables
    uniqueOriginIDs = renamevars(uniqueOriginIDs,"ORIGIN_AIRPORT_ID","ID");
    uniqueDestIDs = renamevars(uniqueDestIDs,"DEST_AIRPORT_ID","ID");
    uniqueIDs = cat(1,uniqueOriginIDs,uniqueDestIDs);
    %removes duplicates
    uniqueIDs = unique(uniqueIDs);
    
    %number of unique IDs
    numOfNodes = height(uniqueIDs);
    
    
    %number of unique flight routes in the data
    uniqueFlights = unique(data(:,["ORIGIN_AIRPORT_ID", "DEST_AIRPORT_ID"]));
    %array holding passenger data for unique flights
    passengerSum = [1;1];
    
    %variable for total network strength
    totalStrength = 0;
    
    for flight = 1:height(uniqueFlights)
        %creates a logical to find passenger numbers for unique flights
        flightPass = (data.ORIGIN_AIRPORT_ID == uniqueFlights(flight,"ORIGIN_AIRPORT_ID").Variables) & (data.DEST_AIRPORT_ID == uniqueFlights(flight,"DEST_AIRPORT_ID").Variables);
        %sums up passenger data for these flights
        allPass = sum(data.PASSENGERS(flightPass));
        %adds this data to average passengers array
        passengerSum(end+1) = allPass;
    
        %adds to the total network strength
        totalStrength = totalStrength + allPass;
    end
    
    
    %a table is created using the arrays containing the data on passengers
    passengerTable = table(passengerSum);
    
    %the first two rows are deleted as these were only for formatting the data
    passengerTable([1,2],:) = [];
    
    
    %produces a table with unique flight routes and their average passengers
    %per journey
    uniqueFlights = [uniqueFlights passengerTable];
    
    
    %converting to adjacency matrix for Polya Urn filtering
    temp1 = table2array(uniqueFlights(:,1));
    temp2 = table2array(uniqueFlights(:,2));
    tempWeights = table2array(uniqueFlights(:,3));
    
    tempG = digraph(temp1,temp2,tempWeights);
    
    adjacencyA = adjacency(tempG,"weighted")

    [i,j,s] = find(adjacencyA);

    airportsPFtable = [i j s];
    
    filter(adjacencyA);
    
    flightPF = PF(adjacencyA,0.6,10,0);
    validatePF(flightPF)

    %table containing location IDs, passengers and p-values
    airportsPFtable = [airportsPFtable flightPF];

    %removes elements not on backbone
    alphaSignificance = 0.05;
    airportsPFtable(airportsPFtable(:, 4) >= alphaSignificance, :) = [];
    
    writematrix(airportsPFtable);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Running hygecdf on input adjancency matrix%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function to run the hyper geometric filter on a given edgelist matrix
function filter(AMatrix)
    [i,j,s] = find(AMatrix);

    alphaSignificance = 0.05/length(i);
    numValidated = 0;
    numChecked = 0;

    totalStrength = sum(sum(AMatrix));
    inStrengths = sum(AMatrix,2);
    outStrengths = sum(AMatrix,1);
    
    %creates a new matrix and array where calculated p-values will be stored
    pMatrix = AMatrix;
    pArray = [];
    percentageArray = [];

    %for loops iterating through rows and columns of adjancency matrix
    for i = 1:height(AMatrix)
        for j = 1:width(AMatrix)
            %disp(num2str(i) + "   " + num2str(j))
            if AMatrix(i,j) == 0
            else 
                %weight of A->B
                x = AMatrix(i,j) - 1;
                %in strength of node A
                M = inStrengths(i,1);
                %out strength of node B
                K = outStrengths(1,j);
                %total network strength
                N = totalStrength;
    
                p = hygecdf(x,N,K,M,"upper");

                %adds calculated values to matrix and array
                pMatrix(i,j) = p;
                pArray(end + 1) = p;

                numChecked = numChecked + 1;
                if p <= alphaSignificance
                    numValidated = numValidated + 1;
                end

                percentageArray(end + 1) = numValidated/numChecked;
            end
        end
    end

    numValidated


    %%%%%
    global initialHyGe;
    if initialHyGe == 0
        initialHyGe = numValidated/numChecked;
    end

    global validatedArrayHyGe;
    validatedArrayHyGe(end + 1) = numValidated/numChecked;
    %%%%%


    %creates a histogram illustrating the frequency of p-values where
    %validated values are in the left most column
    hold on
    edges = [-0.1 alphaSignificance:0.1:1.1];
    %figure1 = histogram(pArray,edges);
    
    %line graph of percentage checked
    checkedArray = 1:1:numChecked;
    plot(checkedArray,percentageArray,'LineWidth',3); %normally this is plot not semilogy
    xlabel('Number of edges checked');
    ylabel('Percentage of edges validated');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Running Polya Urn on input adjancency matrix%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validatePF(AMatrix)
    [i,j,s] = find(AMatrix);

    alphaSignificance = 0.05;
    numValidatedPF = 0;
    numChecked = 0;

    percentageArray = [];

    for i = 1:height(AMatrix)
        for j = 1:width(AMatrix)
            if AMatrix(i,j) == 0
            else 
                numChecked = numChecked + 1;
                if AMatrix(i,j) <= alphaSignificance %< is normal
                        numValidatedPF = numValidatedPF + 1;
                end
                percentageArray(end + 1) = numValidatedPF/numChecked;
            end
        end
    end
    numValidatedPF
    %numChecked

    %%%%%
    global initialPF;
    if initialPF == 0
        initialPF = numValidatedPF/numChecked;
    end

    global validatedArrayPF;
    validatedArrayPF(end + 1) = numValidatedPF/numChecked;
    %%%%%

    %line graph of percentage checked
    checkedArray = 1:1:numChecked;
    plot(checkedArray,percentageArray,'LineWidth',3); %normally this is plot not semilogy
    xlabel('Number of edges checked');
    ylabel('Percentage of edges validated');
end


function [ P ] = PF( W,a,apr_lvl,parallel )
%PF gives the P-values prescribed by the Polya Filter for each link of the network 
%   W: is the adjecency matrix
%   a: is the free parameter of the filter
%      note1: if a<0 estimates the free parameter using Maximum Likelihood extrimation
%      note2: when a=1 the Polya Filter concides with the disparity filter
%   apr_lvl: is a constant defining the regime under which the approximate form of the p-value (Eq. (6) of the paper below) can be used.
%     For example, setting apr_lvl = 10 will trigger the use of the approximate form for every link such that s > 10*k/a, w > 10, and s-w>10*k/a. 
%     The approximate form of the p-value is much faster to compute. If not given it is 
%     automatically set to apr_lvl=10, if apr_lvl = 0 the approximate form is always used.
%     The option apr_lvl=0 is prescribed (and automatically done) in networks with non-intger
%     weights.
%   parallel: 0 => no parallel computing, ~=0 => parallel computing    

% The Polya Filter is described in the academic paper:
% "A Pólya urn approach to information filtering in complex networks" by R. Marcaccioli and G. Livan available at https://rdcu.be/bmIjz

% NOTE: The precision of the calculation is set to the precision of the machine (usually 1e-12).
%  to reach the precision reached in the paper we used the Matlab toolbox "Multiple Precision Toolbox"

%check if the network is symmetric (i.e. undirected) and get the edge list for both cases
if issymmetric(W)==1
    U = triu(W);
    [i,j,w] = find(U); 
else
    [i,j,w] = find(W); 
end

%get the degrees and strengths
k_in(:,1) = full(sum(W~=0,1));
s_in(:,1) = full(sum(W,1));
k_out(:,1) = full(sum(W~=0,2));
s_out(:,1) = full(sum(W,2));

k_in = k_in(j);
k_out = k_out(i);
s_in = s_in(j);
s_out = s_out(i);

%if a<0 get the ML estimates
if a<0
   disp('Starting the ML estimation')
   [a,err] = get_ML_estimate(W);
   disp(['Estimation ended, a = ', num2str(a)])
end

%use the asymptotic form if non integer wights are present
if sum(mod(w,1))~=0
    apr_lvl = 0;
end

%calculate p-values
p(:,1) = polya_cdf(w,s_in,k_in,a,apr_lvl,parallel);
p(:,2) = polya_cdf(w,s_out,k_out,a,apr_lvl,parallel);

P = min(p(:,1),p(:,2));

%handle the case k=1
P(k_in==1) = p(k_in==1,2);
P(k_out==1) = p(k_out==1,1);
P((double(k_in==1)+double(k_out==1))==2) = 1;

end

function [p] = polya_cdf(w,s,k,a,L,parallel)

p = nan(length(w),1);

if a==0 %binomial case
    p = binocdf(w,s,1./k);
else
    %check where approximation can be performed
    idx1 = s-w>=L*((k-1)./a+1);
    idx2 = w>=L*max(1./a,1);
    idx3 = s>=L*max(k./a,1);
    idx4 = k>=L*(a-1+1e-20);
    idx = (double(idx1)+double(idx2)+double(idx3)+double(idx4))==4;

    %calculate the p-values that can be approximated
    p(idx) = 1/gamma(1/a)*(1-w(idx)./s(idx)).^((k(idx)-1)/a).*(w(idx).*k(idx)./(s(idx)*a)).^(1/a-1);

    %calculate the p-values that cannot be aprroximated
    idx = find(double(idx)==0);
    if isempty(idx)==0
        if parallel==0
            for ii=1:length(idx)
                n = s(idx(ii));
                A = 1/a;
                B = (k(idx(ii))-1)./a;
                x = (0:1:w(idx(ii))-1)';
                p(idx(ii)) = 1 - sum(exp(gammaln(n+1)+betaln(x+A,n-x+B)-gammaln(x+1)-gammaln(n-x+1)-betaln(A,B)));
            end
        else
            aux = nan(length(idx),1);
            parfor ii=1:length(idx)
                n = s(idx(ii));
                A = 1/a;
                B = (k(idx(ii))-1)./a;
                x = (0:1:w(idx(ii))-1)';
                aux(ii,1) = 1 - sum(exp(gammaln(n+1)+betaln(x+A,n-x+B)-gammaln(x+1)-gammaln(n-x+1)-betaln(A,B)));
            end
            p(idx) = aux;
        end
    end
end

%catch rounding errors (use the mp-toolbox for higher prcision)
p(p<0) = 0;

end

function [a_best,err] = get_ML_estimate(W) 
%This function can be used to get the maximum likelihood estimation of the
%free parameter "a". It will set the filter on the network's own
%heterogeneity and it will therefore produce ultra-sparse backbones
%   W: is the adjecency matrix

%check if the network is symmetric (i.e. undirected) and get the edge list for both cases
if issymmetric(W)==1
    U = triu(W);
    [i,j,w] = find(U); 
else
    [i,j,w] = find(W); 
end

%get the degrees and strengths
k_in(:,1) = full(sum(W~=0,1));
s_in(:,1) = full(sum(W,1));
k_out(:,1) = full(sum(W~=0,2));
s_out(:,1) = full(sum(W,2));

w = [w;w];
k = [k_out(i);k_in(j)];
s = [s_out(i);s_in(j)];

%get rid of the links with degree 1
w(k==1) = [];
s(k==1) = [];
k(k==1) = [];

f = @(a)eq_from_lhood([w k s],a);
g = @(a)lhood([w k s],a);
options_lsq = optimoptions(@lsqnonlin,'Display','off');

%first find the maximum value of the likelihood
x0 = lsqnonlin(g,0.5,0,15,options_lsq);
%use this value to calculate the value that put the derivative to 0
[a_best,err] = lsqnonlin(f,x0,0,15,options_lsq);

%try higher precision if feval is not close to 0
if err>1e-6
    options_lsq = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',2000,'FunctionTolerance',1e-15,'StepTolerance',1e-25,'OptimalityTolerance',1e-25,'Display','off');
    [a_best,err] = lsqnonlin(f,x0,0,15,options_lsq);
    if err>1e-6
        disp('Try stricter minimization options')
    end
end

end

function [ out ] = eq_from_lhood( X,a )
%derivative of the likelihood that must be put equal to 0
    w = X(:,1);
    k = X(:,2);
    s = X(:,3);
    DL = a.^(-2).*(psi(a.^(-1))+((-1)+k).*psi(a.^(-1).*((-1)+k))+(-1).*k.*psi(a.^(-1).*k)+k.*psi(a.^(-1).*k+s)+(-1).*psi(a.^(-1)+w)+(-1).*((-1)+k).*psi(a.^(-1).*((-1)+k+a.*s+(-1).*a.*w)));
    DL(isnan(DL)) = 0;
    out = double(sum(DL));
    if isinf(out)==1
        out=1e100;
    end
end

function [ out ] = lhood( X,a )
%likelihood that needs to be maximised
    w = X(:,1);
    k = X(:,2);
    s = X(:,3);
    P = polya_pdf(w,s,k,a);
    P(P==0)=1e-20;
    out = sum(-log(P));
    if isinf(out)==1
        out = sign(out)*1e100;
    end
end

function [p] = polya_pdf(w,s,k,a)
%pdf of the distribution
if a==0 %binomial case
    p = binocdf(w,s,1./k);
else       
    n = s;
    A = 1/a;
    B = (k-1)./a;
    x = w;
    p = exp(gammaln(n+1)+betaln(x+A,n-x+B)-gammaln(x+1)-gammaln(n-x+1)-betaln(A,B));
end
end



