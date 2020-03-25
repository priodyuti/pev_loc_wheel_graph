% Run the SIS model on the given networks
% email for communication: prioyutipradhan@gmail.com

%% SIS part has taken from GIAN slides of "Inmaculada Leyva", Complex Systems Group & GISC, 
%Universidad Rey Juan Carlos, 28933 Móstoles, Madrid, Spain held at IIT Indore
%% 

clear all
%%------------------------------------------
%fd = fopen('wheel_infect-dublin_comb_loc.txt','rt');
%fd = fopen('ER_random.txt','rt');
%fd = fopen('wheel_email_comb_deloc.txt','rt');
%fd = fopen('SF_greedy_1.txt','rt');
% fd = fopen('wheel_random_comb_loc.txt','rt');
fd = fopen('regular.txt','rt');


% f_out = fopen('wheel_random_comb_loc_infection.txt','wt');
%f_out = fopen('wheel_infect-dublin_comb_loc_infection.txt','wt')
f_out = fopen('regular_infection.txt','wt');


formatSpec = '%d %d';
sizeA = [2 Inf];

Y = fscanf(fd,formatSpec,sizeA);
Y = Y';
fclose(fd);

N = max(max(Y(:,1)),max(Y(:,2))); % Size of the network
A = zeros(N,N);
t = size(Y,1);
fprintf('Number of nodes:%d\n',N);
fprintf('Number of edges:%d\n',t);
for i = 1:t
    A(Y(i,1),Y(i,2)) = 1; 
    A(Y(i,2),Y(i,1)) = 1;
end    
clear Y;

[max_deg,index] = max(sum(A,2));
max_eigval = max(eig(A));

              % States vector: 0=healthy/susceptible, 1=Infected, 2=Recovered
V = zeros(1,N); % Initially all the nodes are susceptible  
%n = index;
n = randi(N)   % Choosing patient 0
V(n) = 1;       % Patient 0 becomes ill % SIR Parameters
beta = (1/max_eigval) + 0.1     % Infection rate
mu =  0.1;       % Recovering rate
T = 100;        % Time
t=1;

I(t) = length(find(V==1));
S(t) = length(find(V==0));

for t=2:T  
  V0 = V;
  infected = [];
  [x infected] = find(V0 == 1);
  for n = 1:length(infected)
     [neighs,y] = find(A(:,infected(n))==1);
     badluck = 1.*(rand(1,length(neighs))< beta );
     contacts = find(badluck>0);
     healthy_neighs = find(V(neighs(contacts)) == 0);
     V(neighs(contacts(healthy_neighs))) = 1;
     clear neighs badluck contacts healthy_neighs
  end
  %%% Healings
  for i=1:length(infected)
     if(rand < mu)
        V(infected(i)) = 0;
     end
  end 
  % Populations at time t
  I(t) = length(find(V==1));
  S(t) = length(find(V==0));
end
% The plot
for i=1:T
  fprintf(f_out,'%f \n',I(i)/N);
end
fclose(f_out);

plot(I./N,'ro-','MarkerSize',3)
%hold on
%plot(S./N,'k*-','MarkerSize',3)
%hold off
xlabel('Time')
ylabel('I')
xlim([1 T])
ylim([0 1])
%legend('I','S')
