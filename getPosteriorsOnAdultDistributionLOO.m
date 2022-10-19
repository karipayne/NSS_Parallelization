function [px] = getPosteriorsOnAdultDistributionLOO(adult_x, adult_y, spatial_sigma, temporal_sigma)
 
if nargin < 3
    spatial_sigma = 96; %Karinote: updated for 2 deg at 48 pixels per degree %1 deg=59.76pixels --> 2 degs on 120 at 60cm
    temporal_sigma = 8; %Karinote: updated to 30fps % number of frames in moving window --> 280ms at 25fps
end
 
[num_frames num_adults] = size(adult_x);
 
px = zeros(num_frames, num_adults);
for frame_i = 1 : num_frames
    fprintf('%d/%d\n', frame_i, num_frames);
    
    for adult_i = 1 : num_adults
    
        cur_adult_x = adult_x(frame_i,adult_i); %get current baby's X/Y coordinates
        cur_adult_y = adult_y(frame_i,adult_i);
        
        if cur_adult_x == 0 || cur_adult_y == 0
            fprintf('Skip participant %d\n', adult_i)
        else
            
           adults_vec = 1 : num_adults;
           adults_vec(adult_i) = [];
        
           mu = [];
           for frame_j = max(1,frame_i-temporal_sigma*6) : min(num_frames,frame_i+temporal_sigma*6)
                this_frame_adults_x = adult_x(frame_j, adults_vec);
                this_frame_adults_y = adult_y(frame_j, adults_vec);
                good_data = logical(this_frame_adults_x & this_frame_adults_y);
                mu = [mu; [this_frame_adults_x(good_data)' this_frame_adults_y(good_data)' ones(sum(good_data), 1).*frame_j]];
            end
 
            single_sigma = eye(3) .* spatial_sigma;
            single_sigma(3,3) = temporal_sigma;
            sigmas = repmat(single_sigma, [1 1 length(mu)]);
            p = ones(1,length(mu)) ./ length(mu);
           obj = gmdistribution(mu,sigmas,p); %create the GMM distribution
        
          px(frame_i,adult_i) = pdf(obj, [cur_adult_x cur_adult_y frame_i]); %sample the adult pdf at the current XY coords
        end
    end
    fprintf('%d\n', mean(px(frame_i,:)) ./ std(px(frame_i,:)));
    plot(mean(px(1:frame_i,:),2) ./ std(px(1:frame_i,:),0,2));
    drawnow;
   
end
 
save('results_adultsLOO', 'px', 'spatial_sigma', 'temporal_sigma', 'num_frames', 'num_adults');



