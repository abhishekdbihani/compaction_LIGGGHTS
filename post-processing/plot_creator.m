% Import post-processing data and conduct analysis
% Abhishek Bihani, March 2020

clear
close all

%% Import csv files and organize as matrix
finale=zeros(1,25);
finale1= zeros(11,25,11);
for iter=1:11
    perc = (iter-1)*10;
    
    orig = 'D:\liggghts\LIGGGHTS-PUBLIC\examples\LIGGGHTS\Tutorials_public\tut\compaction\post\';
    
    NewFileName = sprintf('VL%g',perc);% create output directory
    folder = [orig NewFileName '\'];
    
    files = dir([folder 'post-processing.csv']); %number of force chain files
    
    thisdir = files(1).name;
    Data=load([folder thisdir]);
    
    c2 = perc/100*ones(size(Data,1),1) ;  % last column added
    Data = [Data c2] ;
    
    [r,c]=size(Data);
    finale=[finale;Data];
    finale1(1:r,1:c,iter) = Data;
end


% create a default color map ranging from red to green
length = 11;
red = [1, 0, 0];
green = [0, 1, 0];
colors_p = [linspace(red(1),green(1),length)', linspace(red(2),green(2),length)', linspace(red(3),green(3),length)'];

%% plot all graphs for different V%L

show_poro= 1;
if show_poro == 1
    % Porosity vs V%L
    
    figure()
    hold on
    
    for count =1:11
        y=finale1(1:11,2,count);  % Porosity
        x=finale1(1:11,25,count); % V%L
        x=x*100;
        plot(x,y,'o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:));
    end
    grid on
  %  title('Porosity vs V%L')
    xlabel('V%L')
    xlim([0 100])
    ylabel('Porosity')
    ylim([0 1])
    
    Legend=cell(11,1);
    for iter=1:11
        perc1 = (iter-1)*10;
        Legend{iter}=strcat('V%L ', num2str(perc1));
    end
    legend(Legend)
    hold off
    legend show
end

% Strain vs count
figure()
hold on
for count =1:11
    x=count;  % coordination number: Large grains - All force chains
    y=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Strain vs count')
ylabel('Axial Strain')
%xlim([0 0.7])
xlabel('count')
%ylim([0 200])
legend(Legend)
hold off
legend show



% coordination number: Large grains - All force chains

figure()
hold on
for count =1:11
    y=finale1(1:11,7,count);  % coordination number: Large grains - All force chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Mean coordination number: Large grains - All force chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('coordination number')
ylim([0 200])
legend(Legend)
hold off
legend show


% coordination number: Small grains - All force chains
figure()
hold on
for count =1:11
    y=finale1(1:11,8,count);  % coordination number: Small grains - All force chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Mean coordination number: Small grains - All force chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('coordination number')
ylim([0 10])
legend(Legend)
hold off
legend show

% coordination number: Large grains - Strong force chains

figure()
hold on
for count =1:11
    y=finale1(1:11,9,count);  % Large grains - Strong force chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Mean coordination number: Large grains - Strong force chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('coordination number')
ylim([0 200])
legend(Legend)
hold off
legend show

% coordination number: Small grains - Strong force chains

figure()
hold on
for count =1:11
    y=finale1(1:11,10,count);  % coordination number: Small grains - Strong force chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Mean coordination number: Small grains - Strong force chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('coordination number')
ylim([0 10])
legend(Legend)
hold off
legend show

% Strong fraction of Large-Large grain chains

figure()
hold on
for count =1:11
    y=finale1(1:11,13,count);  % Strong fraction of Large-Large grain chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Strong fraction of Large-Large grain chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('Fraction')
ylim([0 1])
set(gca, 'YScale', 'log')
legend(Legend)
hold off
legend show

% Strong fraction of Large-Small grain chains
figure()
hold on
for count =1:11
    y=finale1(1:11,14,count);  % Strong fraction of Large-Small grain chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Strong fraction of Large-Small grain chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('Fraction')
ylim([0 1])
set(gca, 'YScale', 'log')
legend(Legend)
hold off
legend show

% Strong fraction of Small-Small grain chains
figure()
hold on
for count =1:11
    y=finale1(1:11,15,count);  % Strong fraction of Small-Small grain chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Strong fraction of Small-Small grain chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('Fraction')
set(gca, 'YScale', 'log')
ylim([0 1])
legend(Legend)
hold off
legend show

% Large-Large fraction of Strong chains
figure()
hold on
for count =1:11
    y=finale1(1:11,16,count);  % Large-Large fraction of Strong chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Large-Large fraction of Strong chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('Fraction')
ylim([0 1])
set(gca, 'YScale', 'log')
legend(Legend)
hold off
legend show

% Large-Small fraction of Strong chains
figure()
hold on
for count =1:11
    y=finale1(1:11,17,count);  % Large-Small fraction of Strong chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Large-Small fraction of Strong chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('Fraction')
ylim([0.01 0.5])
%set(gca, 'YScale', 'log')
legend(Legend)
hold off
legend show

% Small-Small fraction of Strong chains
figure()
hold on
for count =1:11
    y=finale1(1:11,18,count);  % Small-Small fraction of Strong chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Small-Small fraction of Strong chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('Fraction')
ylim([0.5 1])
%set(gca, 'YScale', 'log')
legend(Legend)
hold off
legend show

% Large-Large fraction of All chains
figure()
hold on
for count =1:11
    y=finale1(1:11,22,count);  % Large-Large fraction of All chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Large-Large fraction of All chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('Fraction')
ylim([0 1])
set(gca, 'YScale', 'log')
legend(Legend)
hold off
legend show

% Large-Small fraction of All chains
figure()
hold on
for count =1:11
    y=finale1(1:11,23,count);  % Large-Small fraction of All chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Large-Small fraction of All chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('Fraction')
ylim([0 1])
set(gca, 'YScale', 'log')
legend(Legend)
hold off
legend show

% Small-Small fraction of All chains
figure()
hold on
for count =1:11
    y=finale1(1:11,24,count);  % Small-Small fraction of All chains
    x=finale1(1:11,5,count); % axial strain
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(count,:),'MarkerEdgeColor',colors_p(count,:),'color',colors_p(count,:)); %
end
grid on
title('Small-Small fraction of All chains')
xlabel('Axial Strain')
xlim([0 0.4])
ylabel('Fraction')
ylim([0 1])
set(gca, 'YScale', 'log')
legend(Legend)
hold off
legend show

% Contribution of each contact type in strong chains (before & after compaction)
figure()
hold on
    last = 5;
    
    y=squeeze(finale1(1,16,1:11));  % Large-Large fraction of Strong chains- Initial
    x=squeeze(finale1(1,25,1:11)); % V%L
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(1,:),'MarkerEdgeColor',colors_p(1,:),'color',colors_p(1,:)); %
    
    y=squeeze(finale1(last,16,1:11));  % Large-Large fraction of Strong chains- Final
    x=squeeze(finale1(last,25,1:11)); % V%L
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(3,:),'MarkerEdgeColor',colors_p(3,:),'color',colors_p(3,:)); %
    

    y=squeeze(finale1(1,17,1:11));  % Large-Small fraction of Strong chains- Initial
    x=squeeze(finale1(1,25,1:11)); % V%L
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(5,:),'MarkerEdgeColor',colors_p(5,:),'color',colors_p(5,:)); %
    
    y=squeeze(finale1(last,17,1:11));  % Large-Small fraction of Strong chains- Final
    x=squeeze(finale1(last,25,1:11)); % V%L
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(7,:),'MarkerEdgeColor',colors_p(7,:),'color',colors_p(7,:)); %
    

    y=squeeze(finale1(1,18,1:11));  % Small-Small fraction of Strong chains- Initial
    x=squeeze(finale1(1,25,1:11)); % V%L
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(9,:),'MarkerEdgeColor',colors_p(9,:),'color',colors_p(9,:)); %
    
    y=squeeze(finale1(last,18,1:11)); % Small-Small fraction of Strong chains- Final
    x=squeeze(finale1(last,25,1:11)); % V%L
    plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(11,:),'MarkerEdgeColor',colors_p(11,:),'color',colors_p(11,:)); %

     
grid on
title('Contact type contribution: Before Compaction (BC), Limited Compaction (LC)')
xlabel('V%L')
xlim([0 1])
ylabel('Fraction')
ylim([0 1])
set(gca, 'YScale', 'log')
legend('Large-Large: BC','Large-Large: LC','Large-Small: BC','Large-Small: LC','Small-Small: BC', 'Small-Small: LC')
hold off
legend show

% % Contribution of each contact type in strong chains (after compaction)
% figure()
% hold on
% 
%     y=squeeze(finale1(9,16,1:11));  % Large-Large fraction of Strong chains
%     x=squeeze(finale1(9,25,1:11)); % V%L
%     plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(1,:),'MarkerEdgeColor',colors_p(1,:),'color',colors_p(1,:)); %
% 
%     y=squeeze(finale1(9,17,1:11));  % Large-Small fraction of Strong chains
%     x=squeeze(finale1(9,25,1:11)); % V%L
%     plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(5,:),'MarkerEdgeColor',colors_p(5,:),'color',colors_p(5,:)); %
% 
%     y=squeeze(finale1(9,18,1:11));  % Small-Small fraction of Strong chains
%     x=squeeze(finale1(9,25,1:11)); % V%L
%     plot(x,y,'-o','MarkerSize',6,'MarkerFaceColor',colors_p(11,:),'MarkerEdgeColor',colors_p(11,:),'color',colors_p(11,:)); %
% 
%      
% grid on
% title('Contribution of each contact type in strong chains (after compaction)')
% xlabel('% V%L')
% xlim([0 1])
% ylabel('Fraction')
% ylim([0 1])
% set(gca, 'YScale', 'log')
% legend('Large-Large','Large-Small','Small-small')
% hold off
% legend show


