%% Generating a random Conservation graph based on different algorithms

function [Inc_Con,Cc_Con,Ab,n_c,e_c,Sptree_branch,Sptree_chord,b,c] = Network_Generation(n,network_flag)
del_node = [];

if network_flag==1
    % Generating a random network based on Erdos Renyi algorithm
    G_ER = erdosRenyi(n,0.5,1);         
    Adj = full(G_ER.Adj);             % Adjacent matrix 
elseif network_flag==2
    % Generating a random small world network
    G_SW = smallw(n,3,0.05);
    Adj = full(adjacency(G_SW));   % Adjacent matrix 
else
    % Generating a random scale free network
    Adj = BAgraph_dir(n,5,3);      % Adjacent matrix 
end
    

%% Converting the generated graph into a conservation graph
while isempty(del_node)
    Adj_dir = Adj;
    % Randomly assigning directions to the edges
    for i=1:n
        for j =i:n
            if Adj_dir(i,j)==1
                if i==j
                    Adj_dir(i,j)=0;
                else
                    if(rand()<=0.6)
                        Adj_dir(i,j)=0;
                        Adj_dir(j,i)=1;
                    else
                        Adj_dir(i,j)=1;
                        Adj_dir(j,i)=0;
                    end
                end
            end
        end
    end      

    % Getting Incidence Matrix of the generated directedgraph
    e = sum(Adj_dir(:));
    Inc_dir = zeros(n,e);
    m=1;
    for i=1:n
        for j=1:n
            if Adj_dir(i,j)==1
                Inc_dir(i,m)=-1;
                Inc_dir(j,m)=1;
                m=m+1;
            end
        end
    end

    % Merging all the source and sink nodes into a single environment node
    for i=1:n
        if ~(any(Inc_dir(i,:)==1)&&any(Inc_dir(i,:)==-1))
            del_node = [del_node i];
        end
    end
end
temp = sum(Inc_dir(del_node,:),1);
Inc_Con = Inc_dir;
Inc_Con(del_node,:) = [];
Inc_Con = [temp;Inc_Con];
Inc_Con(:,~any(Inc_Con,1)) = [];  % Deleting empty columns -  edges which form self loops at environment node
n_c = size(Inc_Con,1);
e_c = size(Inc_Con,2);

%% Getting Adjancy matrix of conservation graph
Adj_Con = zeros(size(n_c,n_c)); 
for i = 1:e_c
    s(i) = find(Inc_Con(:,i)==-1);
    t(i)= find(Inc_Con(:,i)==1);
    Adj_Con(s(i),t(i)) = i;
end

% %% Drawing the conservation graph
% G_con = digraph(Adj_Con);
% plot(G_con,'Layout','force','EdgeLabel',G_con.Edges.Weight)

%% Finding a spanning tree of the graph
Adj_temp = Adj_Con;
for i=1:n_c
    for j=1:n_c
        if Adj_Con(i,j)==0
            Adj_temp(i,j) = Adj_Con(j,i);
        end
    end
end
G_temp = graph(Adj_temp,'upper');
T_con = minspantree(G_temp);
Sptree_branch = T_con.Edges.Weight';
Sptree_chord = setdiff(1:e_c,Sptree_branch);
b = n_c-1; % No of branches
c = e_c-b; % No of chords        
     
%% Getting the f-cutset Matrix from Incidence Matrix 
Cc_Con = (inv(Inc_Con(1:end-1,Sptree_branch))*Inc_Con(1:end-1,Sptree_chord)); % Non-identity part

%% Partial Network Information
Ab = Inc_Con(1:end-1,Sptree_branch);






