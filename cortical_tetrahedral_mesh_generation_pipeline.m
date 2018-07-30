clear;close all;

% Matlab 2013b and later version is needed.

addpath('E:\lyx\tetgen\iso2mesh'); % iso2mesh lib
addpath('E:\lyx\tetgen\FreesurferMatlab'); % FreesurferMatlab lib

exe_dir = 'E:\lyx\tetgen'; % path to exe
start_path = 'E:\lyx\tetgen'; % Path to the data
toplevel = uigetdir(start_path);

% generate paths to each subject
dire = dir(toplevel);
dire(~[dire.isdir]) = [];
dire = dire(3:end); % remove . and ..(the first two directories)

subFolder = genpath(toplevel);
remain = subFolder;
listOfFolderNames = {};

while(true)
    [singleSubFolder, remain] = strtok(remain, ';');
    if isempty(singleSubFolder)
        break;
    end
    listOfFolderNames = [listOfFolderNames singleSubFolder];
end
j = 1;
for i = 2:1:length(listOfFolderNames) %%%
    SubFolderName{j} = listOfFolderNames{i};
    j = j+1;
end

numberOfFolders = length(SubFolderName);
meshfix_dir = [exe_dir,'\meshfix.exe'];
FilterMesh_dir = [exe_dir,'\FilterMesh.exe'];
MeshSimplify_dir = [exe_dir,'\MeshSimplify.exe'];
tetgen_dir = [exe_dir,'\tetgen.exe'];

% processing each subject
for n = 1:numberOfFolders
    
    tic
    
    dir = SubFolderName{n};
    % input data: pial surface and white surface in the format of original
    % Freesurfer's result
    [pv, pf] = freesurfer_read_surf([dir,'\lh.pial']);
    [wv, wf] = freesurfer_read_surf([dir,'\lh.white']);
    % Format transform
    saveoff(pv,pf,[dir,'\lh_pial.off']);
    saveoff(wv,wf,[dir,'\lh_white.off']);
    % Initially fixing some geometric errors by meshfix.exe
    pcommand = [meshfix_dir,' ',[dir,'\lh_pial.off -q -a 0.01 -o '],[dir,'\lh_pial_fixed']];
    system(pcommand);
    wcommand = [meshfix_dir,' ',[dir,'\lh_white.off -q -a 0.01 -o '],[dir,'\lh_white_fixed']];
    system(wcommand);
    
    [pv,pf]=readoff([dir,'\lh_pial_fixed.off']);
    [wv,wf]=readoff([dir,'\lh_white_fixed.off']);
    
    write_surf_m2( [dir,'\lh_pial'], pv, pf);
    write_surf_m2( [dir,'\lh_white'], wv, wf);
    
    % genus-0 preserving downsampling
    pcommand = [FilterMesh_dir,' ',dir,'\lh_pial.m -trisubdiv > ',dir,'\ptemp1.m'];
    system(pcommand);
    pcommand = [FilterMesh_dir,' ',dir,'\ptemp1.m -removekey wid > ',dir,'\ptemp2.m'];
    system(pcommand);
    pcommand = [MeshSimplify_dir,' ',dir,'\ptemp2.m -minedgelength -nfaces 120000 -simplify -nooutput -outmesh ',dir,'\lh_pial_sub.m'];
    system(pcommand);
    wcommand = [FilterMesh_dir,' ',dir,'\lh_white.m -trisubdiv > ',dir,'\wtemp1.m'];
    system(wcommand);
    wcommand = [FilterMesh_dir,' ',dir,'\wtemp1.m -removekey wid > ',dir,'\wtemp2.m'];
    system(wcommand);
    wcommand = [MeshSimplify_dir,' ',dir,'\wtemp2.m -minedgelength -nfaces 120000 -simplify -nooutput -outmesh ',dir,'\lh_white_sub.m'];
    system(wcommand);
    
    % Organising downsampled data
    [pv, pf, pnormal, p_vid, p_fid] = readDownSampleMesh([dir,'\lh_pial_sub.m']);
    [wv, wf, wnormal, w_vid, w_fid] = readDownSampleMesh([dir,'\lh_white_sub.m']);
    pv = pv';
    wv = wv';
    pf = pf';
    wf = wf';
    for i = 1:size(p_vid,2)
        [row,col] = find(pf == p_vid(i));
        for m = 1:size(row,1)
            pf(row(m),col(m)) = i;
        end
    end
    for i = 1:size(w_vid,2)
        [row,col] = find(wf == w_vid(i));
        for m = 1:size(row,1)
            wf(row(m),col(m)) = i;
        end
    end
    
    saveoff(pv,pf,[dir,'\tpial.off']);
    saveoff(wv,wf,[dir,'\twhite.off']);
    
    % Another round of error check to guarantee the mesh is error free
    pcommand = [meshfix_dir,' ',dir,'\tpial.off -q -a 0.01 -o ',dir,'\tpial_fixed'];
    system(pcommand);
    wcommand = [meshfix_dir,' ',dir,'\twhite.off -q -a 0.01 -o ',dir,'\twhite_fixed'];
    system(wcommand);
    [pv,pf]=readoff([dir,'\tpial_fixed.off']);
    [wv,wf]=readoff([dir,'\twhite_fixed.off']);
    
    % intersection detection and mark
    [P, w_index, p_index] = intersect(wv,pv,'rows');
    
    [pn,connnum,count]=meshconn(pf,size(pv,1));
    [wn,connnum,count]=meshconn(wf,size(wv,1));
    
    count = 1;
    remove_pv = p_index';
    remove_wv = w_index';
    while count < 3
        for i = 1:size(p_index,1)
            remove_pv = [remove_pv, unique(pn{remove_pv(1,i)})];
        end
        remove_pv = unique(remove_pv);
        count = count + 1;
    end
    
    while count < 6
        for i = 1:size(w_index,1)
            remove_wv = [remove_wv, unique(wn{remove_wv(1,i)})];
        end
        remove_wv = unique(remove_wv);
        count = count + 1;
    end
    
    pial = createMeshFacesOfVertex(pv, pf);
    white = createMeshFacesOfVertex(wv, wf);
    
    %pial = removeVerticesPatch(pial,pv(remove_pv,:));
    %white = removeVerticesPatch(white,wv(remove_wv,:));
    pv = pial.vertices;
    pf = pial.faces;
    
    wv = white.vertices;
    wf = white.faces;
    
    %saveoff(pv,pf,'./lh_pial.off');
    %saveoff(wv,wf,'./lh_white.off');
    write_surf_m2( [dir,'\lh_pial'], pv, pf);
    write_surf_m2( [dir,'\lh_white'], wv, wf);
    
    % Now, merge pial and white surfaces into one whole cortical surface
    command = [FilterMesh_dir,' ',dir,'\lh_pial.m -merge ',dir,'\lh_white.m > ',dir,'\lh_combine.m'];
    system(command);
    command = [MeshSimplify_dir,' ',dir,'\lh_combine.m -outmesh ',dir,'\lh_combine.m'];
    system(command);
    
    % Organizing combined surface
    [v, f, normal, v_id, f_id] = readDownSampleMesh([dir,'\lh_combine.m']);
    v = v';
    f = f';
    for i = 1:size(v_id,2)
        [row,col] = find(f == v_id(i));
        for m = 1:size(row,1)
            f(row(m),col(m)) = i;
        end
    end
    saveoff(v,f,[dir,'\lh_combine.off']);
    
    % Self-intersection detection and mark
    command = [tetgen_dir,' -d ',dir,'\lh_combine.off > ',dir,'\self_intersection.txt'];
    system(command);
    
    [neigh,connnum,count]=meshconn(f,size(v,1));
    
    intersect_face_pair = readIntersectionFile([dir,'\self_intersection.txt']);
    intersect_face = unique(intersect_face_pair);
    inter_vertex = unique(f(intersect_face(:),:));
    
    count = 1;
    remove_v = inter_vertex';
    while count < 2
        for i = 1:size(inter_vertex,1)
            remove_v = [remove_v, unique(neigh{remove_v(1,i)})];
        end
        remove_v = unique(remove_v);
        count = count + 1;
    end
    remove_f = [];
    parfor i = 1:size(f,1)
        f(i,:) = sort(f(i,:));
        if ismember(f(i,:)',remove_v','rows')
            remove_f = [remove_f,i];
        end
    end
   
    remove_pv = remove_v(remove_v <= size(pv,1));
    remove_wv = remove_v(remove_v > size(pv,1))-size(pv,1);
    a = setdiff(remove_v',inter_vertex,'rows');
    
    % this part is only for visualization
    rgbp = zeros(size(v,1),3);
    rgbp(a,1) = 200;
    rgbp(inter_vertex,1) = 0;
    other = setdiff(1:size(v,1),a);
    other = setdiff(other,inter_vertex);
    %rgbp(ia,2) = 113;
    %rgbp(ia,3) = 158;
    rgbp(other,1) = 200;
    rgbp(other,2) = 200;
    rgbp(other,3) = 235;
    writeOFF([dir,'\overlap_p.off '], v,f,[],rgbp,[]);
    
    % Solving the self-intersection
    [pn,connnum,count]=meshconn(pf,size(pv,1));
    [wn,connnum,count]=meshconn(wf,size(wv,1));
    
    WTR = triangulation(wf,wv);
    WN = vertexNormal(WTR);
    %
    PTR = triangulation(pf,pv);
    PN = vertexNormal(PTR);
    % write
    wv(remove_wv,:) = wv(remove_wv,:) - 0.5 * WN(remove_wv,:);
    pv(remove_pv,:) = pv(remove_pv,:) + 0.5 * PN(remove_pv,:);
    
    %[pf,I] = orient_outward(pv,pf,(1:size(pf,1))');
    %[wf,I] = orient_outward(wv,wf,(1:size(wf,1))');
    
    saveoff(pv,pf,[dir,'\tpial.off']);
    saveoff(wv,wf,[dir,'\twhite.off']);
    
    % Another checking to confirm that no more errors are introduced
    pcommand = [meshfix_dir,' ',dir,'\tpial.off -q -a 0.01 -o ',dir,'\tpial_fixed'];
    system(pcommand);
    wcommand = [meshfix_dir,' ',dir,'\twhite.off -q -a 0.01 -o ',dir,'\twhite_fixed'];
    system(wcommand);
    
    [pv2,pf2]=readoff([dir,'\tpial_fixed.off']);
    [wv2,wf2]=readoff([dir,'\twhite_fixed.off']);
    
    % Mark vertices that are modified for solving the self-intersection
    [C,ia] = setdiff(wv2,wv,'rows');
    [C,ib] = setdiff(pv2,pv,'rows');
    % wC = 255 .* ones(size(wv2,1),3);
    % wC(ia,1) = 120;
    % writeOFF('lh_white_second_remove.off',wv2,wf2,[],wC,[]);
    remove_wv = ia;
    remove_pv = ib;
    remove_id = [remove_pv;remove_wv + size(pv2,1)];
    writeRemoveID(remove_id,[dir,'\remove_id.txt']);
    
    % write out final surface file
    write_surf_m2( [dir,'\lh_pial_final'], pv2, pf2);
    write_surf_m2( [dir,'\lh_white_final'], wv2, wf2);
    
    
    command = [FilterMesh_dir,' ',dir,'\lh_pial_final.m -merge ',dir,'\lh_white_final.m > ',dir,'\lh_combine_temp.m'];
    system(command);
    command = [MeshSimplify_dir,' ',dir,'\lh_combine_temp.m -outmesh ',dir,'\lh_combine_final.m'];
    system(command);
    [v, f, normal, v_id, f_id] = readDownSampleMesh([dir,'\lh_combine_final.m']);
    v = v';
    f = f';
    for i = 1:size(v_id,2)
        [row,col] = find(f == v_id(i));
        for m = 1:size(row,1)
            f(row(m),col(m)) = i;
        end
    end
    %saveoff(v,f,'lh_combine_final.off');
    
    % Tetrahedral mesh generation
    [seeds_region,p0,v0,t,idx]=surfinterior(v,f);
    [seeds_hole,p0,v0,t,idx]=surfinterior(wv2,wf2);
    %seeds=surfseeds(v,f);
    p0=min(v(:,1:3));
    p1=max(v(:,1:3));
    seeds = [seeds_region;seeds_hole];
    savesurfpoly(v,f,[],seeds,p0,p1,[dir,'\lh_combine_final.poly'],[]);
    
    command = [tetgen_dir,' -pq1.2a0.8AY ',dir,'\lh_combine_final.poly'];
    system(command);
    savesurfpoly(v,f,seeds_hole,seeds_region,p0,p1,[dir,'\lh_combine_final2.poly'],[]);
    command = [tetgen_dir,' -pq1.2a0.8AY ',dir,'\lh_combine_final2.poly'];
    system(command);
    %cd(userpath);
    %system('cd ..');
    %system(['cp -R ',SubFolderDir{n},' /home/hohokam/YHFan/Data/ADNI2/AD']);
    disp(['Subject ', SubFolderName{n},' is treated!!']);
    
    t = toc;
    
    % recording completed id
    record = fopen('E:\lyx\tetgen\record422.txt', 'At');
    fprintf(record, [dir(end-6:end), ' ']);
    fprintf(record, [num2str(t), '\n']);
    fclose(record);
    
end
%toc

