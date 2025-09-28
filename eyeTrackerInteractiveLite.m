function eyeTrackerInteractiveLite(pathIn)
% eyeTrackerInteractiveLite  Interactive pupil/iris circle detection.
% Works with a video file OR a folder of images.
%
% AUTO-SESSION MODE (the part you asked for):
% If pathIn is a VIDEO file, this program will FIRST create a session folder
% (under ./eye_sessions/<videoName>_<timestamp>/), then EXTRACT FRAMES at a
% stride of 15 (i.e., keep 1 of every 15 frames) up to a MAX of 25 frames,
% and finally run the tracker on that folder of PNGs.  <-- You can change these two numbers.

    % ---- pick a source if not given ----
    if nargin==0 || isempty(pathIn)
        choice = questdlg('Use a video file or a folder of images?', 'Input', ...
                          'Video file','Folder of images','Video file');
        switch choice
            case 'Video file'
                [f,p] = uigetfile({'*.mp4;*.mov;*.avi','Video Files'},'Select video');
                if isequal(f,0), return; end
                pathIn = fullfile(p,f);
            case 'Folder of images'
                p = uigetdir(pwd,'Select image folder');
                if isequal(p,0), return; end
                pathIn = p;
            otherwise
                return;
        end
    end

    % ================== AUTO-SESSION (your request) ==================
    % If the input is a video file, convert it to a small image-sequence session
    % using: stride = 15 (keep 1 of every 15 frames), maxFrames = 25.
    if ~isfolder(pathIn) && exist(pathIn,'file') ~= 0
        stride    = 15;   % <-- keep 1 of every 15 frames
        maxFrames = 25;   % <-- cap total frames to 25
        pathIn = autoMakeSessionFromVideo(pathIn, 'eye_sessions', stride, maxFrames);
        if isempty(pathIn), return; end  % failed to build session
    end
    % =================================================================

    % ---- reader (abstracts folder vs video) ----
    R = makeReader(pathIn);
    if R.nFrames==0, errordlg('No frames found.'); return; end

    % ---- state ----
    S.R = R; S.pathIn = pathIn;
    S.frameIdx = 1; S.playing = true;
    S.sens = 0.85; S.edgeThr = 0.10;
    S.rmin = 15; S.rmax = max(30, round(min(R.height,R.width)/3));
    S.polarity = 'dark'; S.metricThr = 0.25;
    S.minCenterDist = 15;

    % ---- figure / UI ----
    S.hFig = figure('Name','Eye Tracker (Lite)','NumberTitle','off', ...
        'Color','k','Units','normalized','Position',[0.1 0.1 0.8 0.8], ...
        'KeyPressFcn',@(src,evt) onKey(src,evt), ...
        'CloseRequestFcn',@(src,evt) onClose());
    S.hAx  = axes('Parent',S.hFig,'Position',[0.03 0.12 0.74 0.85]); 
    S.hTxt = uicontrol('Style','text','Units','normalized','Position',[0.03 0.02 0.94 0.06], ...
                       'BackgroundColor','k','ForegroundColor','w','HorizontalAlignment','left');

    uicontrol('Style','text','String','Sensitivity','Units','normalized','Position',[0.80 0.87 0.10 0.03], 'ForegroundColor','w','BackgroundColor','k');
    S.hSens = uicontrol('Style','slider','Min',0.5,'Max',0.99,'Value',S.sens,'Units','normalized','Position',[0.90 0.87 0.07 0.03], 'Callback',@(s,e) refreshParams());
    uicontrol('Style','text','String','EdgeThr','Units','normalized','Position',[0.80 0.82 0.10 0.03], 'ForegroundColor','w','BackgroundColor','k');
    S.hEdg  = uicontrol('Style','slider','Min',0.01,'Max',0.5,'Value',S.edgeThr,'Units','normalized','Position',[0.90 0.82 0.07 0.03], 'Callback',@(s,e) refreshParams());
    uicontrol('Style','text','String','Rmin','Units','normalized','Position',[0.80 0.77 0.10 0.03], 'ForegroundColor','w','BackgroundColor','k');
    S.hRmn  = uicontrol('Style','slider','Min',5,'Max',400,'Value',S.rmin,'Units','normalized','Position',[0.90 0.77 0.07 0.03], 'Callback',@(s,e) refreshParams());
    uicontrol('Style','text','String','Rmax','Units','normalized','Position',[0.80 0.72 0.10 0.03], 'ForegroundColor','w','BackgroundColor','k');
    S.hRmx  = uicontrol('Style','slider','Min',10,'Max',500,'Value',S.rmax,'Units','normalized','Position',[0.90 0.72 0.07 0.03], 'Callback',@(s,e) refreshParams());

    S.hPlay = uicontrol('Style','togglebutton','String','Play/Pause (space)','Units','normalized','Position',[0.80 0.64 0.17 0.05], 'Callback',@(s,e) setPlay());
    uicontrol('Style','pushbutton','String','<< Prev','Units','normalized','Position',[0.80 0.58 0.08 0.05], 'Callback',@(s,e) step(-1));
    uicontrol('Style','pushbutton','String','Next >>','Units','normalized','Position',[0.89 0.58 0.08 0.05], 'Callback',@(s,e) step(+1));
    S.hPol  = uicontrol('Style','pushbutton','String','Polarity: DARK','Units','normalized','Position',[0.80 0.50 0.17 0.05], 'Callback',@(s,e) togglePol());
    S.hBoth = uicontrol('Style','togglebutton','String','Try both','Units','normalized','Position',[0.80 0.44 0.17 0.05], 'Callback',@(s,e) redraw());
    uicontrol('Style','pushbutton','String','Save annotated frame','Units','normalized','Position',[0.80 0.36 0.17 0.05], 'Callback',@(s,e) saveFrame());
    uicontrol('Style','pushbutton','String','Batch save all frames','Units','normalized','Position',[0.80 0.30 0.17 0.05], 'Callback',@(s,e) saveAll());

    % ---- timer ----
    S.t = timer('ExecutionMode','fixedSpacing','Period',0.05,'TimerFcn',@(s,e) onTimer());
    guidata(S.hFig,S);
    redraw(); start(S.t);

    % ===== inline callbacks (use guidata; no additional function defs) =====
    function refreshParams()
        S = guidata(S.hFig);
        S.sens    = get(S.hSens,'Value');
        S.edgeThr = get(S.hEdg,'Value');
        S.rmin    = round(get(S.hRmn,'Value'));
        S.rmax    = round(get(S.hRmx,'Value'));
        if S.rmax<=S.rmin, S.rmax = S.rmin + 1; set(S.hRmx,'Value',S.rmax); end
        guidata(S.hFig,S); redraw();
    end
    function setPlay()
        S = guidata(S.hFig);
        S.playing = ~S.playing; set(S.hPlay,'Value', ~S.playing);
        guidata(S.hFig,S);
    end
    function step(d)
        S = guidata(S.hFig);
        S.frameIdx = max(1, min(S.R.nFrames, S.frameIdx + d));
        guidata(S.hFig,S); redraw();
    end
    function togglePol()
        S = guidata(S.hFig);
        if strcmpi(S.polarity,'dark')
            S.polarity='bright'; set(S.hPol,'String','Polarity: BRIGHT');
        else
            S.polarity='dark'; set(S.hPol,'String','Polarity: DARK');
        end
        guidata(S.hFig,S); redraw();
    end
    function onTimer(~,~)
        if ~ishandle(S.hFig), try stop(S.t); delete(S.t); end; return; end
        S = guidata(S.hFig);
        if S.playing && S.frameIdx < S.R.nFrames
            S.frameIdx = S.frameIdx + 1; guidata(S.hFig,S); redraw();
        end
    end
    function onKey(~,evt)
        switch lower(evt.Key)
            case 'space', setPlay();
            case 'rightarrow', step(+1);
            case 'leftarrow',  step(-1);
            case 's', saveFrame();
            case 'q', onClose();
        end
    end
    function onClose(~,~)
        if isvalid(S.hFig), try stop(S.t); delete(S.t); end; end
        if isvalid(S.hFig), delete(S.hFig); end
    end
    function redraw()
        S = guidata(S.hFig);
        I = S.R.read(S.frameIdx);
        imshow(I,'Parent',S.hAx); hold(S.hAx,'on');

        G = imgaussfilt(im2gray(I), 1.0);
        rRange = [max(5,S.rmin) max(6,S.rmax)];
        [C,Rr,M] = detectCircles(G, rRange, S.sens, S.edgeThr, S.polarity, logical(get(S.hBoth,'Value')));

        n = min([size(C,1), numel(Rr), numel(M)]); C = C(1:n,:); Rr = Rr(1:n); M = M(1:n);
        keep = true(n,1);
        if ~isempty(M), keep = keep & (M>=S.metricThr); end
        for i=1:n
            if ~keep(i), continue; end
            d = sqrt(sum((C - C(i,:)).^2,2));
            dup = d < S.minCenterDist; dup(i)=false;
            keep(dup) = false;
        end
        C = C(keep,:); Rr = Rr(keep);

        if ~isempty(C)
            th = linspace(0,2*pi,120);
            for k=1:numel(Rr)
                xs = C(k,1) + Rr(k)*cos(th);
                ys = C(k,2) + Rr(k)*sin(th);
                plot(S.hAx, xs, ys, 'g-', 'LineWidth',1.5);
                plot(S.hAx, C(k,1), C(k,2), 'r+', 'LineWidth',1.5);
            end
        end
        title(S.hAx, sprintf('Frame %d / %d', S.frameIdx, S.R.nFrames), 'Color','w');
        set(S.hTxt,'String', sprintf('Sens=%.2f  EdgeThr=%.2f  R=[%d %d]  Pol=%s  TryBoth=%d', ...
             S.sens, S.edgeThr, rRange(1), rRange(2), S.polarity, logical(get(S.hBoth,'Value'))));
        drawnow limitrate;
    end
    function saveFrame(~,~)
        S = guidata(S.hFig);
        [outDir, base] = defaultOut(S.pathIn);
        if ~exist(outDir,'dir'), mkdir(outDir); end
        fpath = fullfile(outDir, sprintf('%s_%05d.png', base, S.frameIdx));
        F = getframe(S.hAx);
        imwrite(F.cdata, fpath);
        fprintf('Saved %s\n', fpath);
    end
    function saveAll(~,~)
        S = guidata(S.hFig);
        wasPlaying = S.playing; S.playing = false; guidata(S.hFig,S);
        answer = inputdlg({'Stride (every Nth):','Max frames:'},'Batch save',1,{'15','25'});
        if isempty(answer), S.playing = wasPlaying; guidata(S.hFig,S); return; end
        stride = max(1, round(str2double(answer{1})));
        maxFrames = max(1, round(str2double(answer{2})));
        [outDir, base] = defaultOut(S.pathIn);
        if ~exist(outDir,'dir'), mkdir(outDir); end
        saved = 0;
        for k = 1:stride:S.R.nFrames
            S.frameIdx = k; guidata(S.hFig,S); redraw();
            F = getframe(S.hAx);
            imwrite(F.cdata, fullfile(outDir, sprintf('%s_%05d.png', base, k)));
            saved = saved + 1;
            if saved >= maxFrames, break; end
            S = guidata(S.hFig); %#ok<NASGU>
        end
        S = guidata(S.hFig); S.playing = wasPlaying; guidata(S.hFig,S);
        msgbox(sprintf('Saved %d annotated frames to:\n%s', saved, outDir));
    end
end

% ---------------- local functions (NOT nested) ----------------
function [C,R,M] = detectCircles(G, rRange, sens, edgeThr, polarity, tryBoth)
    if tryBoth
        [c1,r1,m1] = imfindcircles(G, rRange, 'ObjectPolarity','bright', ...
            'Sensitivity',sens,'EdgeThreshold',edgeThr,'Method','TwoStage');
        [c2,r2,m2] = imfindcircles(G, rRange, 'ObjectPolarity','dark', ...
            'Sensitivity',sens,'EdgeThreshold',edgeThr,'Method','TwoStage');
        C = [c1;c2]; R = [r1;r2]; M = [m1;m2];
    else
        [C,R,M] = imfindcircles(G, rRange, 'ObjectPolarity',polarity, ...
            'Sensitivity',sens,'EdgeThreshold',edgeThr,'Method','TwoStage');
    end
end

function R = makeReader(pathIn)
    R.isFolder = isfolder(pathIn);
    if R.isFolder
        exts = {'*.png','*.jpg','*.jpeg','*.bmp','*.tif','*.tiff'};
        L = []; for i=1:numel(exts), L = [L; dir(fullfile(pathIn, exts{i}))]; end %#ok<AGROW>
        if isempty(L), R.nFrames=0; R.read=@(~) []; R.height=0; R.width=0; return; end
        [~,ord] = sort({L.name}); L = L(ord);
        R.nFrames = numel(L);
        I1 = imread(fullfile(L(1).folder, L(1).name));
        R.height = size(I1,1); R.width = size(I1,2);
        R.L = L;
        R.read = @(k) imread(fullfile(R.L(k).folder, R.L(k).name));
    else
        % For safety, if a video sneaks through here, read by timestamp
        [nF,h,w] = getVideoInfo(pathIn);
        R.nFrames = nF; R.height = h; R.width = w; R.path = pathIn;
        R.read = @(k) readVideoFrameAt(R.path, k);
    end
end

function [nFrames,h,w] = getVideoInfo(videoPath)
    v = VideoReader(videoPath);
    nFrames = max(1, floor(v.FrameRate * v.Duration));
    I1 = readVideoFrameAt(videoPath, 1);
    h = size(I1,1); w = size(I1,2);
end

function I = readVideoFrameAt(videoPath, k)
    v = VideoReader(videoPath);
    t = max(0, (k-1)/v.FrameRate);
    v.CurrentTime = t;
    I = readFrame(v);
end

function [outDir, base] = defaultOut(pathIn)
    if isfolder(pathIn)
        outDir = fullfile(pathIn, '_annotated'); base = 'frame';
    else
        [p,b,~] = fileparts(pathIn);
        outDir = fullfile(p, b + "_annotated_frames"); base = b;
    end
end

function sessionDir = autoMakeSessionFromVideo(videoPath, baseOutDir, stride, maxFrames)
% autoMakeSessionFromVideo  Build a lightweight image-sequence session
% from a video, saving PNG frames at the given stride and cap.

    if nargin<2 || isempty(baseOutDir), baseOutDir = fullfile(pwd,'eye_sessions'); end
    if nargin<3 || isempty(stride),     stride = 15;  end
    if nargin<4 || isempty(maxFrames),  maxFrames = 25; end

    ts = datestr(now,'yyyy-mm-dd_HHMMSS');
    [~,b,~] = fileparts(videoPath);
    sessionName = sprintf('%s_%s', b, ts);
    sessionDir  = fullfile(baseOutDir, sessionName);
    if ~exist(sessionDir,'dir'), mkdir(sessionDir); end

    try
        v = VideoReader(videoPath);
        kOut = 0; k = 0;
        while hasFrame(v)
            img = readFrame(v); k = k + 1;
            if mod(k-1, stride) ~= 0, continue; end
            kOut = kOut + 1;
            imwrite(img, fullfile(sessionDir, sprintf('frame_%05d.png', kOut)));
            if kOut >= maxFrames, break; end
        end
        fprintf('Auto-session: saved %d frames to %s\n', kOut, sessionDir);
    catch ME
    % Use identifier + formatted message (satisfies Code Analyzer)
    warning(ME.identifier, 'Auto-session failed: %s', ME.message);
    sessionDir = '';
end

end
