%% EE415 Term Project Final Code
% _by Tuba Gök 2575306_

clc; close all; clearvars;
%% 
% 

forward_backprojection_gui
%% I. Defining GUI Function

function forward_backprojection_gui
    % the main figure window
    hFig = figure('Name', 'Forward & Backprojection GUI', 'NumberTitle', 'off', ...
                  'Position', [100, 100, 1200, 700]);

    % Forward Projection Section
    uicontrol('Style', 'text', 'Position', [20, 650, 150, 20], 'String', 'Forward Projection', 'FontWeight', 'bold');

    uicontrol('Style', 'text', 'Position', [20, 620, 100, 20], 'String', 'Image File:');
    imageFileInput = uicontrol('Style', 'edit', 'Position', [120, 620, 150, 20]);

    uicontrol('Style', 'text', 'Position', [20, 590, 100, 20], 'String', 'Number of Beams:');
    numBeamsInput = uicontrol('Style', 'edit', 'Position', [120, 590, 150, 20]);

    uicontrol('Style', 'text', 'Position', [20, 560, 100, 20], 'String', 'Step Size:');
    stepSizeInput = uicontrol('Style', 'edit', 'Position', [120, 560, 150, 20]);

    uicontrol('Style', 'text', 'Position', [20, 530, 100, 20], 'String', 'Specific Angle:');
    specificAngleInput = uicontrol('Style', 'edit', 'Position', [120, 530, 150, 20]);

    forwardButton = uicontrol('Style', 'pushbutton', 'Position', [20, 500, 100, 30], ...
                              'String', 'Run Forward', 'Callback', @runForwardProjection);

    % Backprojection Section
    uicontrol('Style', 'text', 'Position', [320, 650, 150, 20], 'String', 'Backprojection', 'FontWeight', 'bold');

    uicontrol('Style', 'text', 'Position', [320, 620, 120, 20], 'String', 'Reconstructed Size:');
    recSizeInput = uicontrol('Style', 'edit', 'Position', [450, 620, 120, 20]);

    uicontrol('Style', 'text', 'Position', [320, 590, 120, 20], 'String', 'Filter Option:');
    filterOptionInput = uicontrol('Style', 'edit', 'Position', [450, 590, 120, 20]);

    uicontrol('Style', 'text', 'Position', [320, 560, 150, 20], 'String', 'Select Projection Source:');
    projectionSource = uicontrol('Style', 'popupmenu', 'Position', [450, 560, 150, 20], ...
                                 'String', {'Use Forward Result', 'Upload .txt File'});

    backButton = uicontrol('Style', 'pushbutton', 'Position', [320, 500, 100, 30], ...
                           'String', 'Run Back', 'Callback', @runBackprojection);

    % Axes for Plotting outputs
    forwardAxes = axes('Parent', hFig, 'Position', [0.1, 0.1, 0.35, 0.4]);
    backAxesOriginal = axes('Parent', hFig, 'Position', [0.55, 0.55, 0.35, 0.4]);
    backAxesReconstructed = axes('Parent', hFig, 'Position', [0.55, 0.1, 0.35, 0.4]);


    forwardResult = []; % for result matrix 
    originalImage = []; % the original image for comparison

    % Forward Projection Callback
    function runForwardProjection(~, ~)
        imageFile = get(imageFileInput, 'String');
        numBeams = str2double(get(numBeamsInput, 'String'));
        stepSize = str2double(get(stepSizeInput, 'String'));
        specificAngle = str2double(get(specificAngleInput, 'String'));

        forwardResult = forward_projection(imageFile, numBeams, stepSize);
        data = load(imageFile);
        originalImage = cell2mat(struct2cell(data));
        

        % for specific angle
        angleIndex = round(specificAngle / stepSize) + 1;
        specificProjection = forwardResult(angleIndex, :);

        axes(forwardAxes);
        plot(specificProjection);
        title(['Projection at theta =', num2str(specificAngle)]);
        xlabel('Beam Index');
        ylabel('p-values');
    end

    % Backprojection Callback
    function runBackprojection(~, ~)
        recSizeStr = get(recSizeInput, 'String');
        filterOption = str2double(get(filterOptionInput, 'String'));
        selectedSource = get(projectionSource, 'Value');

        if selectedSource == 1
            projectionData = forwardResult;
        else
            % loading from file
            [fileName, filePath] = uigetfile('*.txt', 'Select a Projection Data File');
            if isequal(fileName, 0)
                return;
            end
            projectionData = load(fullfile(filePath, fileName));
        end

        recSize = str2num(recSizeStr); %#ok<ST2NM>
        reconstructed_image = backprojection(projectionData, recSize, filterOption);
        

        % for comparison
        axes(backAxesOriginal);
        imagesc(originalImage);
        colormap('gray');
        title('Original Image');
        colorbar;

        axes(backAxesReconstructed);
        imagesc(reconstructed_image);
        colormap('gray');
        if filterOption == 1
            title('Reconstructed Image without Filtering');
        elseif filterOption == 2
            title('Reconstructed Image with Filtering');
        else
            title('Reconstructed Image with Filtering, Hamming Window');
        end
        colorbar;

    
    end
end

%% II. Defining Forward Projection Function

function result = forward_projection(image, number_of_beams, step_size)
    data = load(image);
    example_image = cell2mat(struct2cell(data));
    
    % getting the size of the image for assigning M value
    [rows, cols] = size(example_image);
    M = rows;
    N = cols;
    
    % Generate t and theta values
    t_values = linspace(-sqrt((M*M)/4  + (N*N)/4 ), sqrt((M*M)/4  + (N*N)/4 ), number_of_beams);
    theta_values = linspace(0, 180 - step_size, (180 / step_size));
    
    theta_rad = deg2rad(theta_values);
    
    % Defining result vector
    result = zeros(length(theta_values), number_of_beams);
    
    % Loop through each theta value (angle)
    for i = 1:length(theta_values)
        theta = theta_rad(i);  
    
        % Loop through each t value (beam position)
        for m = 1:length(t_values)
            t = t_values(m);  
    
            %to make sure that they are integer values
            x_points = round(linspace(-cols/2, cols/2, cols)); 
            y_points = round(linspace(-rows/2, rows/2, rows));
    
            % Calculate y-values from x (line equation)
            y_determined_x = (t - x_points * cos(theta)) / sin(theta);
    
            % Calculate x-values from y (line equation)
            x_determined_y = (t - y_points * sin(theta)) / cos(theta);
    
            % Find valid intersection points within image bounds, excluding
            % the points outside the image
            valid_x = (y_determined_x >= -rows/2) & (y_determined_x <= rows/2);
            valid_y =  (x_determined_y >= -cols/2) & (x_determined_y <= cols/2);
    
            % Extract the valid intersection points
            intersection_points_x = [x_points(valid_x)', y_determined_x(valid_x)'];
            intersection_points_y = [x_determined_y(valid_y)', y_points(valid_y)'];
    
            % Combine the intersection points from both x and y intersections
            intersection_points = [intersection_points_x; intersection_points_y];
            
    
            % Sort the points, excluding the repeated sorted_points if
            % exist any
            sorted_points = unique(sortrows(intersection_points),'rows');
            
    
            % Extract the x and y coordinates of the sorted points
            x_coords = sorted_points(:, 1);
            y_coords = sorted_points(:, 2);
            
    
            % Calculate distances between consecutive points (vectorized)
            dx = diff(x_coords);
            dy = diff(y_coords);
            distances = sqrt(dx.^2 + dy.^2);
            
    
            % Calculate midpoints between consecutive points (vectorized)
            midpoints_x = (x_coords(1:end-1) + x_coords(2:end)) / 2;
            midpoints_y = (y_coords(1:end-1) + y_coords(2:end)) / 2;
            
    
            % Calculate row and column indices for each midpoint
            row_data = (M / 2) - floor(midpoints_y);  
            column_data = (N / 2) + ceil(midpoints_x);  
    
            % Remove NaN values from distances, row_data, and column_data
            valid_idx = ~isnan(distances) & ~isnan(row_data) & ~isnan(column_data);
            distances = distances(valid_idx);
            row_data = row_data(valid_idx);
            column_data = column_data(valid_idx);
    
            % Ensure that row_data and column_data are within the valid bounds
            row_data = max(1, min(rows, round(row_data))); 
            column_data = max(1, min(cols, round(column_data))); 
    
            % Calculate the weighted sum of pixel values and distances
            for t = 1:length(row_data)
                result(i, m) = result(i, m) + distances(t) * example_image(row_data(t), column_data(t));
            end
        end
    end
    
    
    %Create the Sonogram
    figure;
    imagesc(theta_values, 1:number_of_beams, result'); % transposing result for correct orientation
    colormap('gray'); 
    colorbar; 
    title('Sonogram (Projection Data)');
    xlabel('Projection Angle (Degrees)');
    ylabel('Beam Number');
    
   


end 
%% 
% 
%% III. Defining Backward Projection Function

function reconstructed_image = backprojection(input_data, reconstructed_image_size, fltr_no)
    
    if ischar(input_data) || isstring(input_data)
        [~, ~, ext] = fileparts(input_data);
        if strcmp(ext, '.mat')
            loaded_data = load(input_data);  
            projection_data = loaded_data.projection_data; %if the loaded data as in case 1, projection data will directly equall to result matrix
        elseif strcmp(ext, '.txt')
            fileID = fopen(input_data, 'r');
            num_projections = fscanf(fileID, '%d', 1);
            num_samples = fscanf(fileID, '%d', 1);
            projection_data = zeros(num_projections, num_samples);
            for i = 1:num_projections
                fscanf(fileID, '%d', 1); % skipping projection index number
                projection_data(i, :) = fscanf(fileID, '%f', num_samples); %creating the result matrix that the function can handle from .txt
            end
            fclose(fileID);
        end
    else
        projection_data = input_data;
    end

    
    [num_angles, num_beams] = size(projection_data); % number of angles and beams must be found from the result matrix since it is agnostic 
    rows = reconstructed_image_size(1);
    cols = reconstructed_image_size(2);
    M = rows;
    N = cols; % our image size may not be square    
    
    t_values = linspace(-sqrt((M*M)/4  + (N*N)/4 ), sqrt((M*M)/4  + (N*N)/4 ), num_beams);
    theta_values = linspace(0, pi, num_angles);
    
    % This part is different than the forward projection, when linspace()
    % is used the x and y axis values are distorted and shown as lines in
    % the reconstructed image
    x = -cols/2: cols/2;
    y = -rows/2: rows/2;
    

    reconstructed_image = zeros(rows, cols);

    % Apply filtering (if needed)
    if fltr_no > 1 % for 1, there's no filtering and for 2, there is filtering without window
        % Ram-lak filter
        freqs = linspace(-1, 1,  num_beams);
        ram_Lak_filter = abs(freqs); 
        
        if fltr_no == 3
            window = hamming( num_beams)'; % hanning() or blackman() can also be used for comparison, lowpass filter
            ram_Lak_filter = ram_Lak_filter .* window; 
        end
        
        % Filtering each slice
        for i = 1:num_angles
            proj_data_fft = fftshift(fft(projection_data(i, :))); % fftshift() to obtain high pass filter
            proj_data_fft = proj_data_fft .* ram_Lak_filter; 
            projection_data(i, :) = real(ifft(ifftshift(proj_data_fft))); % Inverse FFT
        end
    end


    for i = 1:length(theta_values)
        theta_value = theta_values(i);  
    
        % beam position
        for m = 1:length(t_values)
            t = t_values(m);  
    
  
            % calculating x and y values from line eqn
            y_determined_x = (t - x * cos(theta_value)) / sin(theta_value);
            x_determined_y = (t - y * sin(theta_value)) / cos(theta_value);
    
            valid_x = (y_determined_x >= -rows/2) & (y_determined_x <= rows/2);
            valid_y =  (x_determined_y >= -cols/2) & (x_determined_y <= cols/2);
    
            intersection_points_x = [x(valid_x)', y_determined_x(valid_x)'];
            intersection_points_y = [x_determined_y(valid_y)', y(valid_y)'];
    
            intersection_points = [intersection_points_x; intersection_points_y];
    
            sorted_points = unique(sortrows(intersection_points),'rows');
    
            x_coords = sorted_points(:, 1);
            y_coords = sorted_points(:, 2);
    
            %distances between consecutive points (vectorized)
            dx = diff(x_coords);
            dy = diff(y_coords);
            distances = sqrt(dx.^2 + dy.^2);
    
            %midpoints between consecutive points (vectorized)
            midpoints_x = (x_coords(1:end-1) + x_coords(2:end)) / 2;
            midpoints_y = (y_coords(1:end-1) + y_coords(2:end)) / 2;
    
            row_data = (M / 2) - floor(midpoints_y);  % Row indices from y-coordinates
            column_data = (N / 2) + ceil(midpoints_x);  % Column indices from x-coordinates
    
            % Remove NaN values from distances, row_data, and column_data
            valid_idx = ~isnan(distances) & ~isnan(row_data) & ~isnan(column_data);
            distances = distances(valid_idx);
            row_data = row_data(valid_idx);
            column_data = column_data(valid_idx);
    
            row_data = max(1, min(rows, round(row_data))); % Ensure rows are within bounds
            column_data = max(1, min(cols, round(column_data))); % Ensure columns are within bounds
  

            for k = 1:length(row_data)
                reconstructed_image(row_data(k),column_data(k)) = reconstructed_image(row_data(k),column_data(k)) + projection_data(i,m)*distances(k);
            end
        end
    end
    
end