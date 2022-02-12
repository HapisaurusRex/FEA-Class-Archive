%% Part A (done)
clear
clc
format longg
Cross_A = readmatrix('crossSectionData.csv');
plot(Cross_A(:,1),Cross_A(:,2))
xlabel('Distance from tip of sample(mm)')
ylabel('Cross Sectional Area(mm^2)')
avrg_cross_area = mean(Cross_A(:,2))


%% Part B (Done)
format longg
%import and figure out the slope via polyfit()
force_disp = readmatrix('ForceDisplacement.csv');
mx = polyfit(force_disp(:,1),force_disp(:,2),1) % x = displacments, y = force, 1st degree poly is a linear

%plot the imported data and generate a line with the slope found
plot (force_disp(:,1),force_disp(:,2))
ylabel('Force(N)')
xlabel('Displacement(mm)')
title('A plot of Force versus Displacement')
hold on
plot (force_disp(:,1), (mx(1,1).*force_disp(:,1))) %print line with slope found in polyfit

%% Part C (done)
clc
clear
format longg
L = 205; %length in mm
Cross_A = readmatrix('crossSectionData.csv'); %import data
A = mean(Cross_A(:,2)); % find the average area (mm^2)
force_disp = readmatrix('ForceDisplacement.csv'); %import data
mx = polyfit(force_disp(:,1),force_disp(:,2),1); % x = displacments, y = force, 1st degree poly is a linear
S_a = mx(1,1); %N/mm
E = (S_a * L)/ A %N/mm^2


%% Part D (done)_
clear
clc 
format longg
element_no = [1:1:100]; 
node_no = element_no + 1;
e = 10e3; % youngs in N/mm^2
area_data = readmatrix('crossSectionData.csv');
answer = zeros(100,2);


for i = element_no(1,:)
    %create alot of matricies that will be used
    element_relation = zeros(element_no(1,i),2);
    k_matrix = zeros((i+1),(i+1));
    force_matrix = zeros(i+1,1); 
    disp_matrix = zeros(i+1,1); 
    split = linspace(1,205,i+1);
    element_length = 205/i; %determine element length based on desired no of elements
    
    for j = 1:1:i %create element relations
        element_relation(j,1) = split(1,j);
        element_relation(j,2)= split(1,(j+1)) ;
    end
    for k = 1:1:element_no(1,i) %figure out average area of each element
        min = round(element_relation(k,1));
        max = round(element_relation(k,2));
        element_area = mean(area_data(min:max,2));
        %determine element stiffness constant and add to global stiffness
        %matrix
        constant = (e*element_area)/element_length;
        n1 = k;
        n2 = k+1;  
        k_matrix(n1,n1) = k_matrix(n1,n1) + constant ;
        k_matrix(n1,n2) = k_matrix(n1,n2) - constant ;
        k_matrix(n2,n1) = k_matrix(n2,n1) - constant;
        k_matrix(n2,n2) = k_matrix(n2,n2) + constant ;
    end
    %apply boundary conditions
    force_matrix(i+1,1) = -1200; %force always will be at the last node
    K_modified = k_matrix;
    fixed_nodes = [1]; %fixed at first node
    for fixed = fixed_nodes
        K_modified (fixed,:) = 0;
        K_modified(fixed,fixed) = 1;
    end
    disp_matrix = K_modified\force_matrix; %find displacement
    answer(i,1) = i;
    S_fea = -1200/disp_matrix(i+1,1); 
    %Sfea = force over displacement at x = 205 aka the last node always
    answer(i,2) = S_fea;
end

plot(answer(:,1),answer(:,2))

for test = 2:1:100 %we can't technically check the first value against the 0th value
    past = answer((test-1),2);
    current = answer(test,2);
    error = abs((past-current)/current) * 100;
    if error < 0.2
        disp('the number of elements required to get 0.2% error BETWEEN 2 simulataneous sims')
        disp(test)
        break
    end
end
%% Part E (done)

%Sa does not match Sfea, Sa is larger than Sfea.
%S is related to Youngs modulus and since we used a lower one than the one
%found in part C, that would explain the lower Sfea

%S = F/disp = AE/L

%Also note that the displacements found for the free end using FEA appear to be greater than
%the free end displacement gathered via the experiment. Which also used
%youngs modulus for calculations

%1200, 0.03 (more accurate one) - experimental data
%1200, 0.057 (from 100 elements) - computational data

%% Part F (elements fixed at 10) (done)
clear
clc
format longg
elements = 10;
nodes = elements + 1;
area_data = readmatrix('crossSectionData.csv');
k_matrix = zeros(nodes,nodes);
disp_matrix = zeros(nodes,1);
force_matrix = zeros(nodes,1);
force_matrix(11,1) = -1200;

constant_matrix = zeros(10,1);
constant_matrix(:,1) = 1/(205/elements);
split = round(linspace(1,205,nodes));

for j = 1:1:10
        element_relation(j,1) = split(1,j);
        element_relation(j,2)= split(1,(j+1)) ;
        constant_matrix(j,1) = constant_matrix(j,1) * (mean(area_data(element_relation(j,1):element_relation(j,2),2)));
end

e = 10e3:0.1e3:25e3; %in N/mm^2 cause everything else is in mm
S_fea = zeros(length(e),2);
formula= zeros(length(e),1);

for i = 1:length(e)
    youngs = e(1,i);
    S_fea(i,1) = youngs;
    for k =1:1:10
        n1 = k;
        n2 = k+1;
        c = constant_matrix(k,1);
        c = youngs * c;
        k_matrix(n1,n1) = k_matrix(n1,n1) + c ;
        k_matrix(n1,n2) = k_matrix(n1,n2) - c ;
        k_matrix(n2,n1) = k_matrix(n2,n1) - c;
        k_matrix(n2,n2) = k_matrix(n2,n2) + c ;
    end
    
    K_modified = k_matrix;
    fixed_nodes = [1];
    for fixed = fixed_nodes
        K_modified (fixed,:) = 0;
        K_modified(fixed,fixed) = 1;
    end
    disp_matrix = K_modified\force_matrix;
    S_fea(i,2) = -1200/disp_matrix(11,1);
end

%import Sa
force_disp = readmatrix('ForceDisplacement.csv'); %import data
mx = polyfit(force_disp(:,1),force_disp(:,2),1); % x = displacments, y = force, 1st degree poly is a linear
S_a = mx(1,1); %N/mm

for rar = 1:length(e)
    formula(rar,1) = sqrt((S_fea(rar,2)-S_a)^2);
end

plot(e,formula(:,1))
ylabel('error')
xlabel('youngs modulus')

%based on my graph, which I think is wrong, the best E to minimize error is
%the lowest value. E = 10000 MPa, this is different from the value
%calcuated in part c which was E = 18027 MPa

%part g

e = 10e6; %based on previous comment
element_stress = zeros(10,1);

for element_counter = 1:1:10
    element_length = 205/10;
    n1 = element_counter;
    n2 = element_counter+1;
    stress_matrix = [-1/element_length, 1/element_length];
    local_disp_matrix = [disp_matrix(n1,1);disp_matrix(n2,1)]; %take values from the displacement matrix
    element_stress(element_counter,1) = e * stress_matrix *local_disp_matrix; 
end

%largest magnitude value of stress occurs at the end of the bone under
%compressive load (x = 205) it makes sense as the fixed end has a larger
%cross sectional area to spread for force around, the end has the lowest
%refering to figure 2 meaning the force is more concentrated over a smaller
%area

%% Part H

%I learnt that FEA allows you to simplify complex geometries into more
%manageable segments that you can analyze and play with

%however it is hard to get correct and if done wrong can present errors

%things that could improve this project
%checks e.g. Sfea should be in the same magnitude of Sa as a hint
%a worked through error calculation so one can know the answer is correct
%these two are important as the later questions stem from getting the
%correct answers before





