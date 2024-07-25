% REU_Deji_MATLAB_Simul_Ver10_15Jul2024.m

% NSF REU 2024: Mathematical modeling in the sciences'
% A Probablistic Study of the Student Admission Process at a University
% By undergraduate Researcher Adedeji Kuforiji, 
%      in collaboration with Nathanael E Johnson.
% Mentored by Dr. Aliakbar M. Haghighi


%The purpose of this program is to simulate the admission process of a 
...university's admission process by modeling the university as a U/U/c 
...queue.  

%This is the special case of the simulation
clear all;

%These variables represent the outcomes of the simulation

admitted = 0; %number of admitted
rejected = 0; %number of rejected
reneged = 0; %number of renegers
balked = 0; %number of balkers
QL = 0; %Current queue length
attended = 0; %number of applications attended

% The following arrays are containers of the above variables, 
%      for analysis across iterations. 
% They are initialized as an array of zeros for the sake of computational 
%      efficiency.

a_pop = zeros(1000,1);
r_pop = zeros(1000,1);
rng_pop = zeros(1000,1);
b_pop = zeros(1000,1);

% The following arrays are containers for the mean values of the above 
%      containers over iterations with different arrival/service rates

MA = zeros(5,1);
MR = zeros(5,1);
MRNG = zeros(5,1);
MB = zeros(5,1);
Mr = zeros(5,1);

%This loop allows us to vary the upper bound of the intervals by ...
%which arrival and service rates are distributed.
for m = 1:5 
%This is the number of arrivals in a day, it is uniformly distributed ...
% discretely over interval [1,39], so that the mean is 20 arrivals a day
   n_a = randi(39); 
%This is the number of departures in a day. it is uniformly distributed ...
%continuously over interval [1,2r-1],where r is our desired mean rate. 
%This also allows us to run the simulation with a varying mean rate of...
%service e.g: when m = 1, the interval is [1,39] so that the mean  is 20,..
% when m=2, the interval is [1,59], so that the mean is 30 etc.
   cmu = unifrnd(1,2*(10+10*m)-1); 
   %This loop is for the number of seasons that are being simulated
   for iteration = 1:20000 
       %This represents the number of days being simulated
       for day = 1:180 
           
           %balking does not occur in the special case
            QL = n_a;
           
            n_a = randi(39); %The number of arrivals for the next 
                             %day is generated
           if(QL - cmu < 0) %Reviewing Applications
               admitted = admitted + QL;
               QL = 0;
               continue;
           else
               admitted = admitted+cmu;
           end
           cmu = unifrnd(1,2*(10+10*m)-1); 
        
       end
       %data collection
       a_pop(iteration) = admitted;
       r_pop(iteration) = rejected;
       rng_pop(iteration) = reneged;
       b_pop(iteration) = balked;
       admitted = 0;
       rejected = 0;
       reneged = 0;
       balked = 0;
       QL = 0;

   end
   %data collection
   MA(m) = mean(a_pop);
   MR(m) = mean(r_pop);
   MRNG(m) = mean(rng_pop);
   MB(m) = mean(b_pop);
   Mr(m) = 10+10*m;
end
P_Error = abs((3600-MA) / 100);
tbl = table(Mr,MA,P_Error)