function vesselanalysis()
    % Vessel Analysis
    % Semi-automatic blood vessel diameter, length, and area, from 2D images.
    %
    % Workflow:
    %   1. Select folder with images.
    %   2. Pick one example image to set a global threshold (slider preview).
    %   3. Enter scale (px / µm). For ECHO 20X it is 5.211 px / µm. 
    %   4. For each image, it will:
    %       - Segment vessels
    %       - Fill holes / bridge gaps
    %       - Skeletonize
    %       - Split skeleton into segments
    %       - For each segment, sample multiple cross-sections perpendicular
    %         to the skeleton and measure vessel diameter.
    %   5. Save CSV + QC overlays (with segment IDs).
    %
    % Output CSV columns:
    %   ImageName, SegmentID,
    %   MeanDiameter_um, StdDiameter_um, MinDiameter_um, MaxDiameter_um,
    %   SegmentLength_um, ApproxArea_um2, TotalArea_um2, NumSamples
    %
    % ApproxArea_um2 = MeanDiameter_um * SegmentLength_um
    %   (rectangular approximation of 2D outlined vessel area for each segment)
    %
    % TotalArea_um2 = sum(ApproxArea_um2) for all segments in that image.
    %
    % Requires Image Processing Toolbox (and Computer Vision Toolbox for insertText).

    %% --- Select folder with images ---
    imgDir = uigetdir('', 'Select folder with vessel images');
    if isequal(imgDir, 0)
        disp('No folder selected. Exiting.');
        return;
    end

    % Get list of images
    exts = {'*.jpg','*.jpeg','*.png','*.tif','*.tiff','*.bmp'};
    fileList = [];
    for k = 1:numel(exts)
        fileList = [fileList; dir(fullfile(imgDir, exts{k}))]; %#ok<AGROW>
    end

    if isempty(fileList)
        error('No image files found in the selected folder.');
    end

    %% --- Select an example image to set threshold ---
    exampleIdx = listdlg('PromptString', 'Select an example image to preview threshold:', ...
                         'ListString', {fileList.name}, ...
                         'SelectionMode', 'single');

    if isempty(exampleIdx)
        disp('No example image selected. Exiting.');
        return;
    end

    examplePath = fullfile(imgDir, fileList(exampleIdx).name);
    Iexample = imread(examplePath);

    % Convert to grayscale 
    if ndims(Iexample) == 3
        IgrayEx = rgb2gray(Iexample);
    else
        IgrayEx = Iexample;
    end
    IgrayEx = im2double(IgrayEx);

    % Otsu as initial suggestion
    initT = graythresh(IgrayEx);

    %% --- Threshold preview GUI ---
    manualThreshold = thresholdPreviewGUI(IgrayEx, initT);

    if isempty(manualThreshold)
        disp('Threshold preview cancelled. Exiting.');
        return;
    end
    fprintf('Chosen threshold = %.4f\n', manualThreshold);

    %% --- Ask for scale: pixels per micron ---
    prompt = {'Enter scale: pixels per micron (px / µm):'};
    dlgTitle = 'Pixel Scale';
    defAns = {'5.21'};
    answer = inputdlg(prompt, dlgTitle, 1, defAns);

    if isempty(answer)
        disp('No scale entered. Exiting.');
        return;
    end

    pxPerMicron = str2double(answer{1});
    if isnan(pxPerMicron) || pxPerMicron <= 0
        error('Invalid scale value. Must be positive numeric.');
    end

    umPerPixel = 1 / pxPerMicron;

    %% --- Parameters (can adjust if needed) ---
    minObjectAreaPixels     = 1000;   % removes tiny specks, original was 50. 1000 seems to exclude background
    minSegmentLengthPixels  = 100;    % ignores very short segments
    gapFillRadiusPixels     = 5;      % closing radius for weak signal gaps

    % Cross-section measurement parameters
    sampleStep              = 5;      % use ~every 'n'th skeleton pixel
    neighborRadiusPixels    = 5;      % for local orientation (PCA neighborhood)
    maxProfileLengthPixels  = 40;     % max half-width to probe along normal
    profileStepPixels       = 0.5;    % step size along normal (subpixel OK)

    %% --- Prepare results storage ---
    allImageNames      = {};
    allSegmentIDs      = [];
    allMeanDiam_um     = [];
    allStdDiam_um      = [];
    allMinDiam_um      = [];
    allMaxDiam_um      = [];
    allLength_um       = [];
    allArea_um2        = [];   % per-segment approximate area
    allTotalArea_um2   = [];   % per-image total area, repeated for each segment
    allNSamples        = [];

    %% --- QC folder for overlays ---
    qcDir = fullfile(imgDir, 'QC_overlays');
    if ~exist(qcDir, 'dir')
        mkdir(qcDir);
    end

    %% --- Process each image ---
    for i = 1:numel(fileList)
        fname = fileList(i).name;
        fpath = fullfile(imgDir, fname);

        fprintf('Processing image %d/%d: %s\n', i, numel(fileList), fname);

        I = imread(fpath);

        [segmentStats, overlayRGB] = analyze_single_vessel_image( ...
            I, umPerPixel, minObjectAreaPixels, ...
            minSegmentLengthPixels, manualThreshold, ...
            gapFillRadiusPixels, sampleStep, neighborRadiusPixels, ...
            maxProfileLengthPixels, profileStepPixels);

        % --- Total area per image (sum of all segment areas, in µm^2) ---
        if ~isempty(segmentStats)
            imageTotalArea_um2 = sum([segmentStats.area]);
        else
            imageTotalArea_um2 = 0;
        end

        % Append to global table arrays
        numSegs = numel(segmentStats);
        for s = 1:numSegs
            allImageNames{end+1,1}    = fname;                        %#ok<AGROW>
            allSegmentIDs(end+1,1)    = s;                             %#ok<AGROW>
            allMeanDiam_um(end+1,1)   = segmentStats(s).meanDiam;      %#ok<AGROW>
            allStdDiam_um(end+1,1)    = segmentStats(s).stdDiam;       %#ok<AGROW>
            allMinDiam_um(end+1,1)    = segmentStats(s).minDiam;       %#ok<AGROW>
            allMaxDiam_um(end+1,1)    = segmentStats(s).maxDiam;       %#ok<AGROW>
            allLength_um(end+1,1)     = segmentStats(s).length;        %#ok<AGROW>
            allArea_um2(end+1,1)      = segmentStats(s).area;          %#ok<AGROW>
            allTotalArea_um2(end+1,1) = imageTotalArea_um2;            %#ok<AGROW>
            allNSamples(end+1,1)      = segmentStats(s).nSamples;      %#ok<AGROW>
        end

        % Save QC overlay
        if ~isempty(overlayRGB)
            [~, baseName, ~] = fileparts(fname);
            imwrite(overlayRGB, fullfile(qcDir, [baseName '_QC.png']));
        end
    end

    %% --- Build and save table ---
    resultsTable = table( ...
        allImageNames, allSegmentIDs, ...
        allMeanDiam_um, allStdDiam_um, ...
        allMinDiam_um, allMaxDiam_um, ...
        allLength_um, allArea_um2, allTotalArea_um2, allNSamples, ...
        'VariableNames', { ...
            'ImageName','SegmentID', ...
            'MeanDiameter_um','StdDiameter_um', ...
            'MinDiameter_um','MaxDiameter_um', ...
            'SegmentLength_um','ApproxArea_um2','TotalArea_um2','NumSamples'});

    outCSV = fullfile(imgDir, 'vessel_diameter_results_crossSections.csv');
    writetable(resultsTable, outCSV);

    fprintf('Done. Results saved to:\n%s\n', outCSV);
    fprintf('QC overlay images saved to:\n%s\n', qcDir);
end


function [segmentStats, overlayRGB] = analyze_single_vessel_image( ...
    I, umPerPixel, minObjectAreaPixels, minSegmentLengthPixels, ...
    manualThreshold, gapFillRadiusPixels, sampleStep, neighborRadiusPixels, ...
    maxProfileLengthPixels, profileStepPixels)

    %% --- Convert to grayscale double [0,1] ---
    if ndims(I) == 3
        Igray = rgb2gray(I);
    else
        Igray = I;
    end
    Igray = im2double(Igray);

    %% --- Threshold using user-chosen threshold ---
    T = manualThreshold;
    if isempty(T) || isnan(T)
        T = graythresh(Igray); % fallback
    end
    BW = imbinarize(Igray, T);

    % Ensure vessels are white
    if nnz(BW) > numel(BW)/2
        BW = ~BW;
    end

    % Remove tiny specks
    BW = bwareaopen(BW, minObjectAreaPixels);

    % Fill internal holes
    BW = imfill(BW, 'holes');

    % Bridge small gaps
    BW = bwmorph(BW, 'bridge', 5);

    % Morphological closing to fill weak gaps
    if gapFillRadiusPixels > 0
        se = strel('disk', gapFillRadiusPixels);
        BW = imclose(BW, se);
    end

    if nnz(BW) == 0
        segmentStats = struct('meanDiam',{},'stdDiam',{}, ...
                              'minDiam',{},'maxDiam',{}, ...
                              'length',{},'nSamples',{}, ...
                              'area',{});
        overlayRGB = [];
        return;
    end

    %% --- Skeletonize ---
    if exist('bwskel','file')  % newer MATLAB
        skel = bwskel(BW);
    else
        skel = bwmorph(BW, 'skel', Inf);
    end
    skel = bwmorph(skel, 'spur', 5);

    if nnz(skel) == 0
        segmentStats = struct('meanDiam',{},'stdDiam',{}, ...
                              'minDiam',{},'maxDiam',{}, ...
                              'length',{},'nSamples',{}, ...
                              'area',{});
        overlayRGB = [];
        return;
    end

    %% --- Segment skeleton (remove branchpoints) ---
    branchPoints = bwmorph(skel, 'branchpoints');
    skelNoBranch = skel & ~branchPoints;
    [Lseg, numSegs] = bwlabel(skelNoBranch);

    [nRows, nCols] = size(BW);

    segmentStats = struct('meanDiam',{},'stdDiam',{}, ...
                          'minDiam',{},'maxDiam',{}, ...
                          'length',{},'nSamples',{}, ...
                          'area',{});

    %% --- Analyze each segment with cross-section measurements ---
    for s = 1:numSegs
        segMask = (Lseg == s);
        nPix = nnz(segMask);
        if nPix < minSegmentLengthPixels
            continue;
        end

        [rAll, cAll] = find(segMask);
        nPoints = numel(rAll);
        if nPoints == 0
            continue;
        end

        % pick subset of skeleton points to sample
        if nPoints <= sampleStep
            idxSample = 1:nPoints;
        else
            idxSample = 1:sampleStep:nPoints;
        end

        diamList_px = [];

        for k = idxSample
            r0 = rAll(k);
            c0 = cAll(k);

            % local neighborhood on skeleton for orientation
            d2 = (double(rAll) - double(r0)).^2 + ...
                 (double(cAll) - double(c0)).^2;
            neighMask = d2 <= neighborRadiusPixels^2;
            rn = double(rAll(neighMask));
            cn = double(cAll(neighMask));

            if numel(rn) < 2
                continue; % not enough points for PCA
            end

            % PCA: columns = [x (=col), y (=row)]
            X = [cn rn];
            mu = mean(X,1);
            Xc = X - mu;
            C = (Xc' * Xc) / size(Xc,1);
            [V,D] = eig(C);
            [~, idxMax] = max(diag(D));
            tangent = V(:, idxMax);
            tangent = tangent / norm(tangent + eps);

            % normal vector perpendicular to vessel
            normal = [-tangent(2); tangent(1)]; % [nx; ny] in (x,y) order

            % measure diameter along normal in binary mask
            diam_px = measureDiameterAlongNormal( ...
                BW, r0, c0, normal, ...
                maxProfileLengthPixels, profileStepPixels, nRows, nCols);

            if ~isnan(diam_px) && diam_px > 0
                diamList_px(end+1) = diam_px; %#ok<AGROW>
            end
        end

        if isempty(diamList_px)
            continue;
        end

        % convert to microns
        diamList_um = diamList_px * umPerPixel;

        % segment length ~ #skeleton pixels * pixel size
        segLength_um = nPix * umPerPixel;

        st.meanDiam = mean(diamList_um);
        st.stdDiam  = std(diamList_um);
        st.minDiam  = min(diamList_um);
        st.maxDiam  = max(diamList_um);
        st.length   = segLength_um;
        st.nSamples = numel(diamList_um);

        % Approximate 2D area of outlined vessel segment (µm^2)
        % Treat segment as a long rectangle: area ≈ mean diameter * length.
        st.area     = st.meanDiam * st.length;

        segmentStats(end+1) = st; %#ok<AGROW>
    end

    %% --- QC overlay with segment IDs ---
    Idisp = mat2gray(Igray);
    IdispRGB = repmat(Idisp, [1 1 3]);

    vesselEdges = bwperim(BW);
    IdispRGB(:,:,1) = max(IdispRGB(:,:,1), vesselEdges); % red edges
    IdispRGB(:,:,2) = max(IdispRGB(:,:,2), skel);        % green skeleton

    overlayRGB = im2uint8(IdispRGB);

    % relabel segments (same as above) to put IDs on overlay
    branchPoints = bwmorph(skel, 'branchpoints');
    skelNoBranch = skel & ~branchPoints;
    [Lseg, numSegs] = bwlabel(skelNoBranch);

    for s = 1:numSegs
        segMask = (Lseg == s);
        if nnz(segMask) == 0
            continue;
        end
        stats = regionprops(segMask, 'Centroid');
        c = stats.Centroid;  % [x,y]
        overlayRGB = insertText(overlayRGB, c, num2str(s), ...
            'FontSize', 18, 'TextColor', 'yellow', ...
            'BoxOpacity', 0, 'AnchorPoint', 'Center');
    end
end


function diam_px = measureDiameterAlongNormal( ...
        BW, r0, c0, normal, maxLen, step, nRows, nCols)
    % BW: binary vessel mask
    % (r0,c0): skeleton point (row, col)
    % normal: [nx; ny] (x,y) unit vector
    % returns diameter in pixels along the normal

    plusDist  = marchToBoundary(BW, r0, c0, normal, +1, maxLen, step, nRows, nCols);
    minusDist = marchToBoundary(BW, r0, c0, normal, -1, maxLen, step, nRows, nCols);

    if isnan(plusDist) || isnan(minusDist)
        diam_px = NaN;
    else
        diam_px = plusDist + minusDist;
    end
end


function d = marchToBoundary(BW, r0, c0, normal, directionSign, maxLen, step, nRows, nCols)
    % walk from (r0,c0) along ± normal until we leave the vessel

    d = NaN;
    nx = normal(1); % in x (col)
    ny = normal(2); % in y (row)

    for t = step:step:maxLen
        r = r0 + directionSign * ny * t;
        c = c0 + directionSign * nx * t;

        rRound = round(r);
        cRound = round(c);

        if rRound < 1 || rRound > nRows || cRound < 1 || cRound > nCols
            d = t;
            return;
        end

        if ~BW(rRound, cRound)
            d = t;
            return;
        end
    end
end


function selectedT = thresholdPreviewGUI(Igray, initT)
    % THRESHOLDPREVIEWGUI
    % GUI for selecting threshold manually with live preview.

    selectedT = [];
    fig = figure('Name','Threshold Preview', ...
                 'NumberTitle','off', ...
                 'MenuBar','none', ...
                 'ToolBar','none');

    % Original
    subplot(1,2,1);
    imshow(Igray, []);
    title('Original');

    % Binary view
    ax2 = subplot(1,2,2);
    BW = imbinarize(Igray, initT);
    hBW = imshow(BW, []);
    title(sprintf('Threshold = %.3f', initT));

    % Slider label
    uicontrol('Style', 'text', ...
              'String', 'Adjust Threshold', ...
              'Units', 'normalized', ...
              'Position', [0.25 0.05 0.5 0.05], ...
              'FontSize', 10);

    % Slider
    slider = uicontrol('Style', 'slider', ...
                       'Min', 0, 'Max', 1, 'Value', initT, ...
                       'Units', 'normalized', ...
                       'Position', [0.25 0.01 0.5 0.04], ...
                       'Callback', @(src,~) updatePreview(src, Igray, hBW, ax2));

    % OK button
    uicontrol('Style', 'pushbutton', ...
              'String', 'Use This Threshold', ...
              'Units', 'normalized', ...
              'Position', [0.80 0.01 0.18 0.06], ...
              'Callback', @(~,~) assignAndClose());

    % Cancel button
    uicontrol('Style', 'pushbutton', ...
              'String', 'Cancel', ...
              'Units', 'normalized', ...
              'Position', [0.02 0.01 0.18 0.06], ...
              'Callback', @(~,~) cancelAndClose());

    uiwait(fig);

    function updatePreview(src, Igray, hBW, ax2)
        T = src.Value;
        BWnew = imbinarize(Igray, T);
        set(hBW, 'CData', BWnew);
        title(ax2, sprintf('Threshold = %.3f', T));
    end

    function assignAndClose()
        selectedT = slider.Value;
        if ishghandle(fig)
            uiresume(fig);
            close(fig);
        end
    end

    function cancelAndClose()
        selectedT = [];
        if ishghandle(fig)
            uiresume(fig);
            close(fig);
        end
    end
end
