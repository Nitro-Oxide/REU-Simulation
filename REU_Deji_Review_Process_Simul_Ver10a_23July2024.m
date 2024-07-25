% REU_Deji_MATLAB_Simul_Ver10_15Jul2024.m

% NSF REU 2024: Mathematical modeling in the sciences'
% A Probablistic Study of the Student Admission Process at a University
% By undergraduate Researcher Adedeji Kuforiji, 
%      in collaboration with Nathanael E Johnson.
% Mentored by Dr. Aliakbar M. Haghighi


%The purpose of this program is to simulate the admission process of a 
...university's admission process by modeling the university as a U/U/c 
...queue. 

%This is the simulation with a varying service rate.
clear all;

%These variables represent the outcomes of the simulation

admitted = 0; %number of admitted
rejected = 0; %number of rejected
reneged = 0; %number of renegers
balked = 0; %number of balkers
QL = 0; %Current queue length
TNA = 0; %Total number of applications submitted

% The following arrays are containers of the above variables, 
%      for analysis across iterations. 
% They are initialized as an array of zeros for the sake of computational 
%      efficiency.

a_pop = zeros(1000,1);
r_pop = zeros(1000,1);

rng_pop = zeros(1000,1);
b_pop = zeros(1000,1);
AQL = zeros(180,1);
MAQL = zeros(1000,1);
ANA = zeros(1000,1);
% The following arrays are containers for the mean values of the above 
%      containers over iterations with different arrival/service rates


MA = zeros(5,1); %Mean Admitted
MR = zeros(5,1); %Mean Rejected
MRNG = zeros(5,1); %Mean Reneged
MB = zeros(5,1); %Mean Balked
MP = zeros(5,1); %Mean Probability of Admittance
Mr = zeros(5,1); %Mean Rate of Service
VA = zeros(5,1); %Variance of Admission
MNA = zeros(5,1); %Mean number of applications submitted

% The following loop allows us to vary the upper bound of the intervals by 
% which arrival and service rates are distributed.

% m is the number of arrivals in a day, 
% it is assumed to be uniformly distributed discretely over interval 
% [1,39], so that the mean is 20 arrivals per day

% n_a is the number of departures in a day, 
% it is uniformly distributed, so that the mean is 30. 
% Service rate (review rate) is continuously over interval [1,2r-1], where 
% r is our desired mean rate. 
% This also allows us to run the simulation with a varying mean rate of
% service e.g: when m = 1, the interval is [1,39] so that the mean  is 20,
% when m=2, the interval is [1,59]..

% For balking, we need a loop so that every arrival gets an individual 
% chance to balk. 
% This will allows for cases where we have either no balk or everyone balks


for m = 1:10 
   n_a = randi(39); 
   cmu = unifrnd(1,2*(10+5*m)-1); %This loop is for the number of seasons 
                                  % that are being simulated
   for iteration = 1:1000  %This represents the number of days being simulated
       for day = 1:180 
           TNA = TNA + n_a;       
           while(n_a ~= 0) %Applicants arrive from infinite source
               n_a = n_a - 1;
               if randi(100) > 28
                   QL = QL + 1;
               else
                   balked = balked +1;
                   continue;
               end
           end
            n_a = randi(39); %The number of arrivals for the next day is generated
           if (QL-cmu >= 0) %Reviewing Applications
               for i = 1:cmu
                   if randi(100) > 40
                      admitted = admitted +1;
                   else
                       rejected = rejected +1;
                   end
               end
               QL = QL - cmu;
           else 
               for i = 1:QL
                   if randi(100) > 40
                      admitted = admitted +1;
                   else
                       rejected = rejected +1;
                   end
               end
               QL = 0;
               cmu = unifrnd(1,2*(10+10*m)-1);
               continue;
           end
           cmu = unifrnd(1,2*(10+10*m)-1);
           
           
           renegers = 1/unifrnd(0.01,0.09);%Inverse of time willing to 
                                            ...wait determines rate of
                                            ...reneging
           if(QL - renegers < 0) %Reneging
               reneged = reneged + QL;
               QL = 0;
               continue;
           else
               reneged = reneged + renegers;
           end
            AQL(day) = QL;
       end
       %data collection
       ANA = TNA;
       MAQL(iteration) = mean(AQL);
       a_pop(iteration) = admitted;
       r_pop(iteration) = rejected;
       rng_pop(iteration) = reneged;
       b_pop(iteration) = balked;
       admitted = 0;
       rejected = 0;
       reneged = 0;
       balked = 0;
       QL = 0;
       TNA = 0;

   end
   %data collection
   MA(m) = mean(a_pop);
   VA(m) = var(a_pop);
   MR(m) = mean(r_pop);
   MP(m) = mean(a_pop) / (mean(a_pop + r_pop) + mean(MAQL));
   MRNG(m) = mean(rng_pop);
   MB(m) = mean(b_pop);
   MNA(m)= mean(ANA);
   Mr(m) = 10+5*m;


end
tbl = table(Mr,MA,MR,MRNG,MB,MP,MNA);
tbl = renamevars(tbl,["Mr","MA","MR","MRNG","MB","MP","MNA"],[ ...
    "Mean Rate of Service","Mean Admitted","Mean Rejected", ...
    "Mean Reneged","Mean Balked","Probability of Admittance", ...
    "Mean Number of Arrivals"])
figure(1)
plot(tbl,1,2,'Color',"black")
xlabel("Mean Rate of Serivce")
ylabel("Mean Number of Applicants Admitted")

figure(2)
plot(tbl,1,4,'Color',"black")
xlabel("Mean Rate of Service")
ylabel("Mean Number of Reneged Application")
tbl2 = table(Mr,MA,MRNG,VA);
tbl2 = renamevars(tbl2,["Mr","MA","MRNG","VA"],["Mean Rate of Service", ...
    "Mean Admitted","Mean Reneged", "Variance of Admitted"])
