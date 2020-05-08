% Import forcechain files and conduct analysis
% Abhishek Bihani, March 2020

clear

%% iterate through all VL's...

for iter=1:11
perc = (iter-1)*10

orig = 'D:\liggghts\LIGGGHTS-PUBLIC\examples\LIGGGHTS\Tutorials_public\tut\compaction\post\';

NewFileName = sprintf('VL%g',perc);% create output directory
folder = [orig NewFileName '\'];

num_files = dir([folder '*.dump']); %number of force chain files

[~,index] = sortrows({num_files.date}.');
num_files = num_files(index);
clear index %sort files in correct order

num_files1 = dir([folder '*.']); %number of grain files
num_files1 = num_files1(3:end);
[~,index] = sortrows({num_files1.date}.');
num_files1 = num_files1(index);
clear index %sort files in correct order

total_files=length(num_files);
real_count=1;
for count = 9 : total_files %start looping from 16000 step (compression begins)
    tic
    clearvars -except count total_files num_files num_files1 folder phi VL z_top strain_axial real_count meanF coord_big_all_avg coord_small_all_avg BB_strong_frac BS_strong_frac SS_strong_frac coord_big_strong_avg coord_small_strong_avg coord_big_weak_avg coord_small_weak_avg loop_number zmax_border zmin_border BB_frac_str BS_frac_str SS_frac_str BB_frac_weak BS_frac_weak SS_frac_weak BB_frac_all BS_frac_all SS_frac_all
    count
    %import force chain files
    filename=num_files(count).name;
    path_=[folder filename];
    
    opts = detectImportOptions(path_, 'NumHeaderLines', 9, 'FileType','text');
    opts.VariableNames = {'x1','y1','z1','x2','y2','z2','f1','f2','f3','fn1','fn2','fn3','ft1','ft2','ft3','cx','cy','cz'}; % variable names x1 y1 z1 x2 y2 z2 f1 f2 f3 fn1 fn2 fn3 ft1 ft2 ft3 cx cy cz
    force_chain_data = readtable(path_,opts);
    force_chain_data = force_chain_data(:,1:18); %remove extra column
    force_chain = table2array(force_chain_data);
    
    %import grain files
    filename=num_files1(count).name;
    path_=[folder filename];
    
    opts = detectImportOptions(path_, 'NumHeaderLines', 9, 'FileType','text');
    grain_data = readtable(path_,opts);
    grain_data = grain_data(:,1:10); %remove extra column
    grains = table2array(grain_data);
    
    for g_count = 1: length(grains)
        grains(g_count,11) = grains(g_count,5) + grains(g_count,6);  % z_center + radius for zmax
        grains(g_count,12) = grains(g_count,5) - grains(g_count,6);  % z_center - radius for zmin
        
    end
    
       
    %% select grains/forcechains within specified boundary range
    
    [xmin,ix1] = min(grains(:,3));
    [ymin,iy1] = min(grains(:,4));
    [zmin,iz1] = min(grains(:,12));
    [xmax,ix2] = max(grains(:,3));
    [ymax,iy2] = max(grains(:,4));
    [zmax,iz2] = max(grains(:,11));
    
    %xmax = 0; %new limit
    %ymax = 0; %new limit
    
    x_len = (xmax+grains(ix2,6))-(xmin-grains(ix1,6));
    y_len = (ymax+grains(iy2,6))-(ymin-grains(iy1,6));
    
    if real_count == 1  % find zmax and zmin before compression begins for reference
        zmax_border = zmax;
        zmin_border = zmin;
        zmax_border=(zmax_border+0.1); % give extra limits
        zmin_border=(zmin_border-0.1);
    end
    

    
    g_count_new = 1; 
    for g_count = 1: length(grains)
        zmax1 = grains(g_count,11);
        zmin1 = grains(g_count,12);

        if zmax1 <= (zmax_border) && zmin1 >= (zmin_border) %if only z > zmin_ global or z < zmax_ global, keep it 
            g_new(g_count_new,:) = grains(g_count,:);
            g_new(g_count_new,1) = g_count_new;
            g_count_new = g_count_new + 1; 
        else
            passs=1;
        end
    end
    
   
    grains=[];
    grains=g_new;
    
    [zmin,iz1] = min(grains(:,12)); %recalculate zmax and zmin after checking if within bounds
    [zmax,iz2] = max(grains(:,11));
    

  % z_len = (zmax+grains(iz2,6))-(zmin+grains(iz1,6));
    z_len = (zmax)-(zmin);
    
    
    fc_new = force_chain;
    
%     fc_count_new = 1;
%     g_count_new = 1;
%     
%     for fc_count = 1 : length(force_chain)
%         if force_chain(fc_count,1) <= xmax && force_chain(fc_count,1) >= xmin && force_chain(fc_count,2) <= ymax && force_chain(fc_count,2) >= ymin && force_chain(fc_count,3) <= zmax && force_chain(fc_count,3) >= zmin
%             fc_new(fc_count_new,:) = force_chain(fc_count,:);
%             fc_count_new = fc_count_new + 1;
%         end
%     end
%     for g_count = 1: length(grains) 
%         if grains(g_count,3) <= xmax && grains(g_count,3) >= xmin && grains(g_count,4) <= ymax && grains(g_count,4) >= ymin && grains(g_count,5) <= zmax && grains(g_count,5) >= zmin
%             g_new(g_count_new,:) = grains(g_count,:);
%             g_count_new = g_count_new + 1; 
%         end
%     end
%     
%     grains=[];
%     grains=g_new;
    
    % assign id's to force chain adjacent grains
    adjGrains1=[fc_new(:,1),fc_new(:,2),fc_new(:,3)];
    adjGrains2=[fc_new(:,4),fc_new(:,5),fc_new(:,6)];
  
    [~,locG1] = ismember(adjGrains1,grains(:,3:5),'rows');
    [~,locG2] = ismember(adjGrains2,grains(:,3:5),'rows');
    
    Index_adjGrains1= grains(locG1,1);
    Index_adjGrains2= grains(locG2,1);
    fc_new(:,24) = Index_adjGrains1;
    fc_new(:,25) = Index_adjGrains2;
    fc_count= (1:1:length(fc_new))'; 
    
    adjGrains_all=[adjGrains1,Index_adjGrains1,adjGrains2,Index_adjGrains2,fc_count];
    
    % assign id's to force chains in grains
    for g_count = 1 : length(grains)
        g_center = grains(g_count,3:5);
               
        locg1 = find(ismember(adjGrains1, g_center,'rows'));
        locg2 = find(ismember(adjGrains2, g_center,'rows'));
        
        g_coord{g_count,:}= union(locg1,locg2);
        
    end
    
    %% calculate grain radii of all pairs, normal force chain magnitude, bond length
    
    for fc_count = 1 : length(fc_new)
        fc_new(fc_count,19)= ((fc_new(fc_count,1)-fc_new(fc_count,16))^2+(fc_new(fc_count,2)-fc_new(fc_count,17))^2+(fc_new(fc_count,3)-fc_new(fc_count,18))^2)^0.5;%radius 1
        fc_new(fc_count,20)= ((fc_new(fc_count,4)-fc_new(fc_count,16))^2+(fc_new(fc_count,5)-fc_new(fc_count,17))^2+(fc_new(fc_count,6)-fc_new(fc_count,18))^2)^0.5;%radius 2
        %fc_new(fc_count,21)= ((fc_new(fc_count,10))^2+(fc_new(fc_count,11))^2+(fc_new(fc_count,12))^2)^0.5; %magnitude normal force
        fc_new(fc_count,21)= abs(fc_new(fc_count,12)); %magnitude normal force
        
        fc_new(fc_count,26)=fc_count;
        if (fc_new(fc_count,19) + fc_new(fc_count,20)) < 0.03
            fc_new(fc_count,22) = 0; % Small-Small (0.01 +0.01)
        elseif (fc_new(fc_count,19) + fc_new(fc_count,20)) > 0.09
            fc_new(fc_count,22) = 2; % Large-Large (0.07 + 0.07)
        else
            fc_new(fc_count,22) = 1; % Large-Small (0.07 + 0.01)
        end
    end
    
    %% Calculate avg force- allot strong (1)/weak force chains (0) for every chain pair 
    
    %%%%%%%%% --------->>>  at least 3 grains & low angle between grains &
    %%%%%%%%% going downward? (not implemented)
    
    meanF(real_count) = mean(fc_new(:,21));
    fc_strong_count = 1;
    for fc_count = 1 : length(fc_new)
        if fc_new(fc_count,21) > meanF(real_count) && grains(Index_adjGrains1(fc_count),10) > 1 && grains(Index_adjGrains2(fc_count),10) > 1 % force chain larger than mean with no floating grains
            fc_new(fc_count,23) = 1; %strong force chain
            fc_strong(fc_strong_count,:) = fc_new(fc_count,:); %save all strong force chains separately for vtk
            fc_strong_count = fc_strong_count+1;
        else
            fc_new(fc_count,23) = 0; %weak force chain
        end
    end
    
       
    %% Calculate fraction of big-big pairs in strong chains etc. (% grain size
    % in strong force chain). Similar for big-small and small-small
    
    BB_strong = nnz(fc_new(:,22) == 2 & fc_new(:,23) == 1);
    BB_weak = nnz(fc_new(:,22) == 2 & fc_new(:,23) == 0);

    BS_strong = nnz(fc_new(:,22) == 1 & fc_new(:,23) == 1);
    BS_weak = nnz(fc_new(:,22) == 1 & fc_new(:,23) == 0);

    SS_strong = nnz(fc_new(:,22) == 0 & fc_new(:,23) == 1);
    SS_weak = nnz(fc_new(:,22) == 0 & fc_new(:,23) == 0);
    
    
    % strong fraction among big-big, big-small and small-small chains
    BB_strong_frac(real_count) = BB_strong /(BB_strong + BB_weak);    
    BS_strong_frac(real_count) = BS_strong /(BS_strong + BS_weak);
    SS_strong_frac(real_count) = SS_strong /(SS_strong + SS_weak);
    
    % BB, BS, SS fraction among strong chains
    
    BB_frac_str(real_count) = BB_strong /(BB_strong + BS_strong + SS_strong);
    BS_frac_str(real_count) = BS_strong /(BB_strong + BS_strong + SS_strong); 
    SS_frac_str(real_count) = SS_strong /(BB_strong + BS_strong + SS_strong); 
    
    % BB, BS, SS fraction among weak chains
    
    BB_frac_weak(real_count) = BB_weak /(BB_weak + BS_weak + SS_weak);
    BS_frac_weak(real_count) = BS_weak /(BB_weak + BS_weak + SS_weak); 
    SS_frac_weak(real_count) = SS_weak /(BB_weak + BS_weak + SS_weak);
    
    % BB, BS, SS fraction among all chains
    
    BB_all = nnz(fc_new(:,22) == 2);
    BS_all = nnz(fc_new(:,22) == 1);
    SS_all = nnz(fc_new(:,22) == 0);
    
    BB_frac_all(real_count) = BB_all /(BB_all + BS_all + SS_all);
    BS_frac_all(real_count) = BS_all /(BB_all + BS_all + SS_all); 
    SS_frac_all(real_count) = SS_all /(BB_all + BS_all + SS_all);  
        
    %% calculate number of unique grains and their radii.

    grains = sortrows(grains); 
    
    r_unique = grains(:,6);
    grains_unique=grains(:,3:5);
    
    %% coordinate number of each grain (all, strong, weak), average of big and small grains
    
    %a_counts = accumarray(i_new,1);
    value_counts = [grains_unique,r_unique,grains(:,10)];
    
    coord_big_all = 0;
    coord_big_strong = 0;
    coord_big_weak = 0;
    coord_small_all = 0;
    coord_small_strong = 0;
    coord_small_weak = 0;
    
    for i = 1: length(value_counts)
        if value_counts(i,4) > 0.05 % big grain
            
            coord_big_all = value_counts(i,5) + coord_big_all; % all neighboring coordinates
            
            strong = 0;
            weak = 0;
            
            all_fc_neighbor = (g_coord{i,:}); % fc id's of all chains connected to the grain
            
            for j = 1 : length(all_fc_neighbor)
                if fc_new(all_fc_neighbor(j),23) == 1 %strong fc
                    strong=strong+1;    
                end
            end
            
            coord_big_strong = strong + coord_big_strong; % all strong neighboring coordinates
            weak = value_counts(i,5) - strong; 
            coord_big_weak = weak + coord_big_weak; % all weak neighboring coordinates
            
        else                        % small grain
            coord_small_all = value_counts(i,5) + coord_small_all; % all neighboring coordinates
            
            strong = 0;
            weak = 0;
            
            all_fc_neighbor = (g_coord{i,:}); % fc id's of all chains connected to the grain
            
            for j = 1 : length(all_fc_neighbor)
                if fc_new(all_fc_neighbor(j),23) == 1 %strong fc
                    strong=strong+1;    
                end
            end
            
            coord_small_strong = strong + coord_small_strong; % all strong neighboring coordinates
            weak = value_counts(i,5) - strong; 
            coord_small_weak = weak + coord_small_weak; % all weak neighboring coordinates           
              
        end
    end
    
    Nsmall = nnz(r_unique < 0.05);
    Nbig = nnz(r_unique > 0.05);
    
    coord_big_all_avg(real_count) =  coord_big_all/ Nbig; %AVG COORD BIG GRAINS- ALL CHAINS
    coord_small_all_avg(real_count) =  coord_small_all/ Nsmall; %AVG COORD SMALL GRAINS- ALL CHAINS
    
    coord_big_strong_avg(real_count) =  coord_big_strong/ Nbig; %AVG COORD BIG GRAINS- STRONG CHAINS
    coord_small_strong_avg(real_count) =  coord_small_strong/ Nsmall; %AVG COORD SMALL GRAINS- STRONG CHAINS
    
    coord_big_weak_avg(real_count) =  coord_big_weak/ Nbig; % AVG COORD BIG GRAINS- WEAK CHAINS
    coord_small_weak_avg(real_count) =  coord_small_weak/ Nsmall; % AVG COORD BIG GRAINS- WEAK CHAINS
    
    
    %% Calculate phi and VL, and strain (decreasing Zmax from 1st in series),compression index
    Rsmall = round(min(r_unique),2);
    Rbig = round(max(r_unique),2);
    
    Vl= Nbig*4/3*pi*Rbig^3;
    Vs= Nsmall*4/3*pi*Rsmall^3;
    
    phi(real_count) = 1 - (Vl+Vs)/(x_len*y_len*z_len);
    VL(real_count) = Vl/(Vl+Vs);
    z_top(real_count)= zmax;   %%%%%%% add radius to ztop
    strain_axial(real_count)=abs((zmax-z_top(1))/(zmax-zmin));
    
       
    %% Output vtk files- grains (pore filter), forcechain (throat filter)- all
    
    all = 1; %all force chains
    forcechain_VTK_ascii(VL(real_count),real_count,folder,value_counts,fc_new,Index_adjGrains1,Index_adjGrains2,all)
    
    %% Output vtk files- grains (pore filter), forcechain (throat filter)- strong
    
    % assign id's to force chain adjacent grains
    adjGrains1=[fc_strong(:,1),fc_strong(:,2),fc_strong(:,3)];
    adjGrains2=[fc_strong(:,4),fc_strong(:,5),fc_strong(:,6)];
  
    [~,locG1] = ismember(adjGrains1,grains_unique,'rows');
    [~,locG2] = ismember(adjGrains2,grains_unique,'rows');
    
    Index_adjGrains1= grains(locG1,1);
    Index_adjGrains2= grains(locG2,1);
    
    adjGrains_strong=[adjGrains1,Index_adjGrains1,adjGrains2,Index_adjGrains2];
    all = 0; %all force chains
    forcechain_VTK_ascii(VL(real_count),real_count,folder,value_counts,fc_strong,Index_adjGrains1,Index_adjGrains2,all)
    
    %% Create grainpack input for calculating permeability/ percolation tortuosity, Pc using output xyzr at evert step
    
    xyzr=value_counts(:,1:4);
    basename = sprintf('xyzr_%g.csv',real_count);
    csvwrite([folder basename],xyzr);
    
    %%
    loop_number(real_count) = real_count;
    real_count = real_count+1;
    toc
end
    %% Save values with compaction csv
    finale =[loop_number; phi; VL; z_top; strain_axial; meanF; coord_big_all_avg; coord_small_all_avg; coord_big_strong_avg; coord_small_strong_avg; coord_big_weak_avg; coord_small_weak_avg; BB_strong_frac; BS_strong_frac; SS_strong_frac; BB_frac_str; BS_frac_str; SS_frac_str; BB_frac_weak; BS_frac_weak; SS_frac_weak; BB_frac_all; BS_frac_all; SS_frac_all];
    finale=finale';
    csvwrite([folder 'post-processing.csv'],finale);
end