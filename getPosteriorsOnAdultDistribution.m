function [px] = getPosteriorsOnAdultDistribution(baby_x, baby_y, adult_x, adult_y, spatial_sigma, temporal_sigma, save_file)

if nargin < 5
    spatial_sigma = 96; %Karinote: updated for 2 deg at 48 pixels per degree %1 deg=59.76pixels --> 2 degs on 120 at 60cm
    temporal_sigma = 8; %Karinote: updated to 30fps % number of frames in moving window --> 280ms at 25fps
    save_file = 'results_nocont';
end

[num_frames num_babies] = size(baby_x);
[num_frames num_adults] = size(adult_x);
mu_all = [];

% find all good fixations and aggregate into huge matrix
for frame_i = 1 : num_frames
    good_data = logical(adult_x(frame_i,:) & adult_y(frame_i,:));
    %disp(good_data)
    %disp("NEW")
    mu_all = [mu_all; [adult_x(frame_i,good_data)' adult_y(frame_i,good_data)' ones(sum(good_data), 1).*frame_i]];
end
%disp(good_data)
%disp(mu_all)
%size(mu_all)
% weighting of all fixations
p_all = ones(1,length(mu_all))/length(mu_all);

% start with empty posteriors
px = zeros(num_frames, num_babies);

% for every frame, every subject, compute the probability they are explained
% by the distribution adult fixations
for frame_i = 1 : num_frames
    fprintf('%d/%d\n', frame_i, num_frames);
    
    
    % only check mu's within a temporal window, as we know the sigma is too
    % low to have mu's more than 3 sigmas away contribute to
    % probabilitiy.
    mu = [];
    for frame_j = max(1,frame_i-temporal_sigma*6) : min(num_frames,frame_i+temporal_sigma*6)
        good_data = logical(adult_x(frame_j,:) & adult_y(frame_j,:));
        mu = [mu; [adult_x(frame_j,good_data)' adult_y(frame_j,good_data)' ones(sum(good_data), 1).*frame_j]];
    end
    
    % compose the covariances
    single_sigma = eye(3) .* spatial_sigma;
    single_sigma(3,3) = temporal_sigma;
    
    % for every fixation, same sigma
    sigmas = repmat(single_sigma, [1 1 length(mu)]);
    
    % get the correct number of weights
    p = p_all(1:length(mu));
    
    % compose final distribution
    obj = gmdistribution(mu,sigmas,p);
    
    % check every subject's fixation, calculate probability
    for baby_i = 1 : num_babies
        cur_baby_x = baby_x(frame_i,baby_i);
        cur_baby_y = baby_y(frame_i,baby_i);
        
        px(frame_i,baby_i) = pdf(obj, [cur_baby_x cur_baby_y frame_i]);
    end
    
    fprintf('%d\n', mean(px(frame_i,:)) ./ std(px(frame_i,:)));
    plot(mean(px(1:frame_i,:),2) ./ std(px(1:frame_i,:),0,2));
    drawnow;
    
    %save('results, 'px', 'frame_i', 'spatial_sigma', 'temporal_sigma', 'num_frames', 'num_babies', 'num_adults');
end