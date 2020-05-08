    function [] = forcechain_VTK_ascii(VL,real_count,folder,value_counts,fc_new,Index_adjGrains1,Index_adjGrains2,all)
    
    % Creates vtk files with grains and force-chains with magnitude of normal force
    % Abhishek Bihani, Mar 2020
    
    NewFileName = sprintf('Paraview_grains_VL-%g_count-%g_all-%g',round(VL,2),real_count,all);% create output directory
    NewFileName = [folder NewFileName '.vtk'];
    
    fid = fopen(NewFileName, 'w');
    fprintf(fid, '# vtk DataFile Version 3.0\nNetwork Flow\nASCII\nDATASET POLYDATA\n');
    
    % writing points position for grains
    ngrains = length(value_counts);
    
    fprintf(fid, '\nPOINTS %u float\n', ngrains);
    for i=1:ngrains
        fprintf(fid, '%f %f %f\n',value_counts(i,1),value_counts(i,2),value_counts(i,3));
    end
        
%     % writing points position for force chains
%     nchains = length(fc_new);
%     fprintf(fid, '\nPOINTS %u float\n', nchains);
%     for i=1:nchains
%         x_fc = (fc_new(i,1)+fc_new(i,4))/2;
%         y_fc = (fc_new(i,2)+fc_new(i,5))/2;
%         z_fc = (fc_new(i,3)+fc_new(i,6))/2;
%         fprintf(fid, '%f %f %f\n',x_fc,y_fc,z_fc);
%     end

    % writing lines connection
    
    nchains = length(fc_new);
    fprintf(fid, '\nLINES %u %u\n', nchains, nchains*3);
    for i=1:nchains
        fprintf(fid, '2 %u %u\n',Index_adjGrains1(i)-1,Index_adjGrains2(i)-1);
    end
    
    % writing grain' scalars (radius)
    
    fprintf(fid, '\nPOINT_DATA %u\nSCALARS pores_radius float\nLOOKUP_TABLE default\n', ngrains);
    for i=1:ngrains
        fprintf(fid, '%f\n', value_counts(i,4));
    end
    
    % writing lines' scalars (throat radius - normal force magnitude of chains)
    
    fprintf(fid, '\nCELL_DATA %u\nSCALARS force_norm_mag float\nLOOKUP_TABLE default\n', nchains);
    for i=1:nchains
        fprintf(fid, '%f\n', fc_new(i,21));
    end
    
    fclose(fid);
    end