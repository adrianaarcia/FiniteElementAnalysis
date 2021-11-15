% --- Initialization
% load data
load elements.txt
load fixed.txt
load force.txt
load mesh.txt

[nelements,~] = size(elements);
[nnodes,~] = size(mesh);

% define Given Vars
Ar = 1;
E = 1e7; 

EA = E*Ar; % calculate commonly used value

% initialize Master stiffness matrix
M = zeros(nnodes*2);

% initialize plot
newplot
hold on

Li = zeros(nelements,1); Ks = Li; % initial lengths and spring constant
for n=1:nelements
    i = elements(n,1); j = elements(n,2);
    nodes = mesh(elements(n,:),:); % [xi, yi; xj, yj]
    
    % --- Disconnection & Localization
    % calculate phi, the angle between global and local coordinate systems
    phi = atan((nodes(2,2)-nodes(1,2))/(nodes(2,1)-nodes(1,1)));%nodes(2,:)-nodes(1,:)
    L = pdist(nodes, "euclidean"); % calculate length
    Li(n) = L;
    % --- Member/element formation
    ks = EA/L; Ks(n) = ks; % spring constant
    Kb = zeros([4,4]); % element stiffness matrix in local coords
    Kb(1,1) = 1; Kb(3,3) = 1; Kb(1,3) = -1; Kb(3,1) = -1;
    Kb = ks * Kb; 
    
    % --- Globalization
    Td = eye(4); Td=Td*cos(phi); % displacement transformation matrix
    Td(1,2)=sin(phi);Td(3,4)=sin(phi); Td(2,1)=-sin(phi); Td(4,3)=-sin(phi);
    Tf = transpose(Td); % force transformation matrix
    
    K = Tf*Kb*Td; %element stiffness matrix in global coords
    
    % --- Merge
    A = zeros(2*nnodes); % augment the stiffness matrix
    yi = 2*i; yj = 2*j; xi = yi-1; xj = yj-1;  
    A(xi,xi)=K(1,1); A(xi,yi)=K(1,2); A(xi,xj)=K(1,3); A(xi,yj)=K(1,4);
    A(yi,xi)=K(2,1); A(yi,yi)=K(2,2); A(yi,xj)=K(2,3); A(yi,yj)=K(2,4);
    A(xj,xi)=K(3,1); A(xj,yi)=K(3,2); A(xj,xj)=K(3,3); A(xj,yj)=K(3,4);
    A(yj,xi)=K(4,1); A(yj,yi)=K(4,2); A(yj,xj)=K(4,3); A(yj,yj)=K(4,4);
    
    M = M + A; % add to master stiffness matrix
end

% --- Apply Boundary Conditions
for n=1:size(fixed,2)
    p1 = fixed(1,n)*2;  % zero out fixed rows/cols and leave 1s on diagonal
    p2 = p1-1;
    M(p1,:) = 0; M(:,p1) = 0; M(p1,p1) = 1;
    M(p2,:) = 0; M(:,p2) = 0; M(p2,p2) = 1;
end

% initialize forces with values from forces.txt
forces = zeros([nnodes*2,1]);
for n=1:size(force,1)
    y = force(n,1)*2; x = y - 1;
    forces(x,:) = force(n,2);
    forces(y,:) = force(n,3);
end

% --- Solve for Displacements & Recover Derived Quantities
U = inv(M)*forces;

D = zeros([size(U,1)/2, 3]); % find max displacement magnitude
for n=1:(size(U,1)/2)
    c = n*2; v = U(c-1:c,:);
    D(n,:)= [transpose(v) norm(v)];
end
[m, index] = max(D(:,3));
fprintf("Max displacement of %f occurs at node %i\n", m, index) % print out max displacement

% --- Show magnified displacement & find max stress
mg = 10 ; % set factor by which to magnify
D = D(:,1:2); Mag = mg*D;

S = zeros([nelements,1]);
for n=1:nelements
    i = elements(n,1); j = elements(n,2);
    
    nodes = mesh(elements(n,:),:); % [xi, yi; xj, yj]
    plot(nodes(:,1),nodes(:,2),'b-') % plot original
    
    show(1,:) = nodes(1,:) + Mag(i,:);
    show(2,:) = nodes(2,:) + Mag(j,:);
    plot(show(:,1),show(:,2),'r-') % plot magnified displacement
    
    nodes(1,:) = nodes(1,:) + D(i,:);
    nodes(2,:) = nodes(2,:) + D(j,:);
    L = pdist(nodes, "euclidean");
    Ax = Ks(n)*(L-Li(n));
    S(n,1) = Ax/Ar; % store stress
end


[m, index] = max(S(:,1));
fprintf("Max stress of %f occurs at element %i", m, index) % print out max displacement
hold off
