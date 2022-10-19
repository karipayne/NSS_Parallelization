

%% automated data reshaping, posterior calculations and NSS for each video
% get data and load them into the appropriate dataframes
% note that 'context' is eye-tracking DIEM data and
% 'nocont' (no context) is the mouse-blur MCBRD data
%subjects = (1:1:7); %{11};

%diem_subjects = (1:1:42);
mcbrd_subjects = (1:1:54);
videos = (1:1:29);

mean_x_binocular_context = []
mean_y_binocular_context = []
mean_x_binocular_nocont = []
mean_y_binocular_nocont = []

% loop over each video (or the specified range of videos)
%for v = 1:length(videos)
for v = 6 % this was used for (length(videos)), but now iterates through only one video so I can submit as separate slurm jobs. A lazy fix.
    
    video = ['vid' num2str(v)]; % don't forget to convert the video number to a string.
    path = ['NSS_ready_data/', video, '/'];
    disp("set path correctly")
    disp(video)
    disp(path)
    
    % get the number of DIEM subjects for the current video
    %fileList = dir(fullfile('vid3', '*.txt'));
    fileList = dir(fullfile(path, '*.txt'));
    fileList2 = string(transpose({fileList.name}))
    num_DIEM = sum(contains(fileList2, "DIEM") )
    
    diem_subjects = (1:1:num_DIEM);
    disp("got the number of DIEM subjects")
    disp(num_DIEM)

    % separate looping due to different subject numbers
    % Start looping through every subject
    for s = 1:length(diem_subjects)
        % load the subject's DIEM data
        %diem = readtable([path, 'DIEM_nile_', num2str(subjects(s)), '.txt']);
        %diem = readtable([path, 'DIEM_', num2str(diem_subjects(s)), '_vid4_NSS.txt'], 'ReadVariableNames', false, 'TreatAsEmpty', 'NA');
        diem = readtable([path, 'DIEM_', num2str(diem_subjects(s)),  '_', video, '_NSS.txt'], 'ReadVariableNames', false, 'TreatAsEmpty', 'NA');
        %disp(size(diem))
        mean_x_binocular_context = [mean_x_binocular_context, diem{:,1}];
        mean_y_binocular_context = [mean_y_binocular_context, diem{:,2}];
    end
    disp("done with diem. starting mcbrd")

    for s = 1:length(mcbrd_subjects)
        % load the subjects MCBRD data
        %mcbrd = readtable([path, 'MCBRD_nile_', num2str(subjects(s)), '.txt']);
        %mcbrd = readtable([path, num2str(mcbrd_subjects(s)), '_vid4_NSS.txt'], 'ReadVariableNames', false, 'TreatAsEmpty', 'NA');
        mcbrd = readtable([path, num2str(mcbrd_subjects(s)), '_', video, '_NSS.txt'], 'ReadVariableNames', false, 'TreatAsEmpty', 'NA');
        %disp(size(mcbrd))
        mean_x_binocular_nocont = [mean_x_binocular_nocont, mcbrd{:,1}];
        mean_y_binocular_nocont = [mean_y_binocular_nocont, mcbrd{:,2}];
    end
    disp("done with mcbrd")


    % now do the posterior computations for the 
    % full mouse blur data and full DIEM dataset
    [px_nocont] = getPosteriorsOnAdultDistribution(mean_x_binocular_nocont, mean_y_binocular_nocont, mean_x_binocular_context, mean_y_binocular_context)
    %csvwrite('Output/vid4_GSimil_MCBRD.csv',px_nocont); %This compares context to no context. When you have more than 2 groups and more than 1 comparison, then copy this multiple times and change parameters
    csvwrite(['Output/', video, '_GSimil_MCBRD.csv'],px_nocont); %This compares context to no context. When you have more than 2 groups and more than 1 comparison, then copy this multiple times and change parameters
    plot(mean(px_nocont,2),'b'); %Context group was the comparison group
    hold on;
    
    disp("done with first posteriors")

    [px_context] = getPosteriorsOnAdultDistributionLOO(mean_x_binocular_context, mean_y_binocular_context)
    %csvwrite('Output/vid4_GSimil_DIEM.csv',px_context); %Leave 1 out for context
    csvwrite(['Output/', video, '_GSimil_DIEM.csv'],px_nocont); %This compares context to no context. When you have more than 2 groups and more than 1 comparison, then copy this multiple times and change parameters
    plot(mean(px_context,2),'r');
    
    disp("done with second posteriors")

    % My normalization across all  frames
    px_contextLOO_normalized = (px_context - mean(mean(px_context,1),2)) ./ std(px_context(:));
    px_nocont_normalized = (px_nocont - mean(mean(px_context,1),2)) ./ std(px_context(:));
    disp("done with normalizaltion")

    figure(2);
    plot(mean(px_nocont_normalized,2),'r');
    hold on;
    plot(mean(px_contextLOO_normalized,2),'b');
    legend({'MCBRD', 'DIEM'});

    save('results_AFET_GSimil', 'px_nocont', 'px_context', 'px_nocont_normalized', 'px_contextLOO_normalized');

    px_nocont_normalized_trans = px_nocont_normalized.'
    px_contextLOO_normalized_trans = px_contextLOO_normalized.'
    px_all = [px_nocont_normalized_trans; px_contextLOO_normalized_trans];
    px_all(1,1:10)
    %px_trans(1,1:10)
    %csvwrite('Output/vid4_NSS_all_trans.csv',px_all);
    csvwrite(['Output/', video, '_NSS_all_trans.csv'],px_all);
    %diary off

end
