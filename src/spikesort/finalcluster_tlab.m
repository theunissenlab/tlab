function varargout = finalcluster_tlab(varargin)
% FINALCLUSTER_TLAB MATLAB code for finalcluster_tlab.fig
%      FINALCLUSTER_TLAB, by itself, creates a new FINALCLUSTER_TLAB or raises the existing
%      singleton*.
%
%      H = FINALCLUSTER_TLAB returns the handle to a new FINALCLUSTER_TLAB or the handle to
%      the existing singleton*.
%
%      FINALCLUSTER_TLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINALCLUSTER_TLAB.M with the given input arguments.
%
%      FINALCLUSTER_TLAB('Property','Value',...) creates a new FINALCLUSTER_TLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before finalcluster_tlab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to finalcluster_tlab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help finalcluster_tlab

% Last Modified by GUIDE v2.5 28-Jun-2012 17:51:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @finalcluster_tlab_OpeningFcn, ...
                   'gui_OutputFcn',  @finalcluster_tlab_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before finalcluster_tlab is made visible.
function finalcluster_tlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to finalcluster_tlab (see VARARGIN)
global PF vv snips currsort uvsorts vtrials vsorts labels sort validtrials tmpv2 tmpclusters sortparams kms kmslast plotmin plotmax tmpunsortedidx tmpsortidx vs tridx tmps2
global trashlabels
% Choose default command line output for finalcluster_tlab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% keyboard
% UIWAIT makes finalcluster_tlab wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% index of presorted
idxPresorted = tmpsortidx;

kms = max(1, length(unique(tmpsortidx(tmpsortidx>0))));

set(handles.numcluster, 'String', num2str(kms));

% Sort the rest of snips
labels = cluster(snips(tmpunsortedidx,:), kms, handles);

% Fill in the sort index for all
uniqueLabels = unique(labels);
for ll = uniqueLabels'
    tmpsortidx(tmpunsortedidx(labels==ll)) = ll;
end

% Create the vector that indicates which snipets can't be sorted correctly
trashlabels = trashsnips(tmpsortidx, snips, idxPresorted);

% Print the number of spikes in each sort.
nSpikesTot = length(tmpsortidx);
nSorts = length(uniqueLabels);
nSpikesSort = zeros(1, nSorts);

fprintf(1, 'Number of spikes sorted %d:\n', nSpikesTot);
is = 0;
for ll = uniqueLabels'
    is = is + 1;
    nSpikesSort(is) = length(find(tmpsortidx == ll & trashlabels == 0));
    fprintf(1,'\t Spikes in sort %d = %d (%.1f%%)\n', ll, nSpikesSort(is), 100*nSpikesSort(is)/nSpikesTot);
end
fprintf(1,'Number of Trash spikes %d (%.1f%%)\n', sum(trashlabels), 100*sum(trashlabels)/nSpikesTot);   

kmslast = kms;

% Make nice plots
namp = 100; 
nt = size(snips,2);
nsnips = size(snips,1);
histval = zeros(nt,namp);

% Calculate the mean and sd for all snippets class
meanSnip = mean(snips,1);
sdSnip = std(snips);

% Find bounds for histograms
minAmp = min(meanSnip - 3*sdSnip);
maxAmp = max(meanSnip + 3*sdSnip);

% x and y labels for the snip plots
t=1:nt;
a=linspace(minAmp, maxAmp, namp);


for iso = 1:min(6, nSorts)
    hh = ['handles.axes' num2str(iso)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset');
    axes(eval(hh));
    
    % Calculate the mean and sd for each snippet class
    meanSnip = mean(snips(tmpsortidx==iso & trashlabels==0,:),1);
    sdSnip = std(snips(tmpsortidx==iso & trashlabels==0,:));
    
   
    % make a histogram
    histval = zeros(namp,nt);
    nkeep = 0;
    for isn=1:nsnips
        if (tmpsortidx(isn) == iso && trashlabels(isn)==0)            
            for it=1:nt
                ampind = floor(namp*(snips(isn,it)-minAmp)/(maxAmp-minAmp))+1;
                if (ampind <= namp && ampind >= 1)
                    histval(ampind,it) = histval(ampind,it) + 1; 
                end
            end
            nkeep = nkeep+1;
        end
        
    end
    histval = histval./nkeep;
    
    pcolor(t, a, (histval));
    shading interp;
    caxis([0 0.001*namp]);
    colormap('gray');
    axis xy;
            
    hold on;
    plot(meanSnip,'r', 'linewidth', 2);
    plot(meanSnip + sdSnip, 'r--');
    plot(meanSnip - sdSnip, 'r--');
    hold 'off';
        
    set(eval(hh), 'ButtonDownFcn', bdf);
end

if (nSorts < 6)
    hh = ['handles.axes' num2str(nSorts+1)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), snips(trashlabels==1,:)'); set(eval(hh), 'YLim', [minAmp maxAmp]);
    hold(eval(hh)); plot(eval(hh), mean(snips(trashlabels==1,:)',2), 'k', 'linewidth', 4)
    set(eval(hh), 'ButtonDownFcn', bdf)
end    


        
% --- Outputs from this function are returned to the command line.
function varargout = finalcluster_tlab_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function numcluster_Callback(hObject, eventdata, handles)
% hObject    handle to numcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numcluster as text
%        str2double(get(hObject,'String')) returns contents of numcluster as a double


% --- Executes during object creation, after setting all properties.
function numcluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numpcs_Callback(hObject, eventdata, handles)
% hObject    handle to numpcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numpcs as text
%        str2double(get(hObject,'String')) returns contents of numpcs as a double


% --- Executes during object creation, after setting all properties.
function numpcs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numpcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in kmeans.
function kmeans_Callback(hObject, eventdata, handles)
% hObject    handle to kmeans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global PF vv ch currsort uvsorts vtrials vsorts labels sort validtrials tmpv tmpv2 subspace plotmin plotmax kmslast tmpunsortedidx tmpsortidx snips
set(handles.status, 'String', 'Wait'); drawnow
kms = str2num(get(handles.numcluster, 'String'));

if get(handles.randomforrest, 'Value') == 1
    set(handles.kmeanscluster, 'Value', 1)
    set(handles.spectralcluster, 'Value',0)
    set(handles.kmedoidscluster, 'Value',0)
    set(handles.randomforrest, 'Value',0)
end

subspace = 0;
% tmpv = cat(2,vv{uvsorts, ch, :});

if kms ~= kmslast
    tmpsortidx = zeros(1, size(snips,1));
    tmpunsortedidx = find(tmpsortidx==0);
    kmslast = kms;
end

if kms == 1
    labels = ones(size(snips,1),1);
else

    labels = cluster(snips(tmpunsortedidx,:), kms, handles);
end

for ll = unique(labels)'
    tmpsortidx(tmpunsortedidx(labels==ll)) = ll;
end

% Create the vector that indicates which snipets can't be sorted correctly
trashlabels = trashsnips(tmpsortidx, snips);

for kk = 1:6
    hh = ['handles.axes' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), snips((labels==kk),:)'); hold(eval(hh)); plot(eval(hh),mean(snips((labels==kk),:)',2), 'k', 'linewidth', 4);
    set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end

 set(handles.status, 'String', 'Ready'); drawnow
% --- Executes on button press in pca.
function pca_Callback(hObject, eventdata, handles)
% hObject    handle to pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PF vv snips currsort uvsorts vtrials vsorts labels sort validtrials tmpv tmpv2 subspace r5pc plotmin plotmax kmslast tmpunsortedidx tmpsortidx
set(handles.status, 'String', 'Wait'); drawnow
kms = str2num(get(handles.numcluster, 'String'));

if kms ~= kmslast
    tmpsortidx = zeros(1, size(snips,1));
    tmpunsortedidx = find(tmpsortidx==0);
    kmslast = kms;
end

if get(handles.randomforrest, 'Value') == 1
    set(handles.kmeanscluster, 'Value', 1)
    set(handles.spectralcluster, 'Value',0)
    set(handles.kmedoidscluster, 'Value',0)
    set(handles.randomforrest, 'Value',0)
end

subspace = 2;
% tmpv = cat(2,vv{uvsorts, ch, :});
if kms == 1
    labels = ones(size(snips(tmpunsortedidx,:),1),1);
%     labels = ones(size(tmpv2,2),1);
else

%     covtv = tmpv2*tmpv2';
%     [U,S,V]= svds(covtv,str2num(get(handles.numpcs, 'String')));
% 
%     ktmpv = (U'*tmpv2);
    dsnips = bsxfun(@minus, snips(tmpunsortedidx,:), msnips);
    ktmpv = (r5pc(:,1:str2num(get(handles.numpcs, 'String')))'*dsnips');
%     ktmpv = (r5pc(:,1:str2num(get(handles.numpcs, 'String')))'*tmpv2);
    
    labels = cluster(ktmpv', kms, handles);
end

for ll = unique(labels)'
    tmpsortidx(tmpunsortedidx(labels==ll)) = ll;
end

% Create the vector that indicates which snipets can't be sorted correctly
trashlabels = trashsnips(tmpsortidx, snips);

for kk = 1:6
    hh = ['handles.axes' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), snips(tmpunsortedidx(labels==kk),:)'); hold(eval(hh)); plot(eval(hh),mean(snips(tmpunsortedidx(labels==kk),:)',2), 'k', 'linewidth', 4);
%     plot(eval(hh), tmpv2(:,(labels==kk))); hold(eval(hh)); plot(eval(hh),mean(tmpv2(:,(labels==kk)),2), 'k', 'linewidth', 4);
    set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end


set(handles.status, 'String', 'Ready'); drawnow
% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 subspace tmpclusters sortparams ens
global r5pc s5pc vs tridx tmps2 tmpsortidx tmpunsortedidx retry sc unit nclasses nfiles fname
global trashlabels

set(handles.status, 'String', 'Wait'); drawnow

kms = str2num(get(handles.numcluster, 'String'));


for ll = unique(labels)'
    tmpsortidx(tmpunsortedidx(labels==ll)) = ll;
end


% Loop through the number of sorts and save each one to a separate hdf5
% file
nsorts = max(tmpsortidx);

for isort=1:nsorts
    
    % First make a copy of the unit structure
    unitSorted = unit;
    
    % Then add new structure fields describing unit type
    unitSorted.sortid = isort;
    if get(eval(['handles.singleunit' num2str(isort)]), 'Value')
        unitSorted.sortType  = 'single';
    elseif get(eval(['handles.multiunit' num2str(isort)]), 'Value')
        unitSorted.sortType  = 'multi';
    else
        unitSorted.sortType  = 'noise';
    end
    
    % Now select relevant spikes and save only those times.
    cnt=0;
    for iclass=1:nclasses
        
        for nfi=1:nfiles(iclass)
            response = unit.class_responses.(unit.classes{iclass}){nfi};              
            ntrials = length(response.trials);
            
            for it=1:ntrials
                trial = response.trials{it};
                ns = length(trial.spikeIds);
                idxsort = [];
                
                for is=1:ns
                    cnt = cnt+1;
                    if (tmpsortidx(cnt) == isort  && trashlabels(cnt)==0)
                        idxsort = [idxsort is];
                    end
                    
                    if mod(cnt, 5000)==0
                        fprintf('.')
                    end
                end
                if isempty(idxsort)
                    unitSorted.class_responses.(unit.classes{iclass}){nfi}.trials{it}.spikeIds = [];
                    unitSorted.class_responses.(unit.classes{iclass}){nfi}.trials{it}.spikeTimes = [];
                    unitSorted.class_responses.(unit.classes{iclass}){nfi}.trials{it}.spikeWav = [];
                else
                    
                    unitSorted.class_responses.(unit.classes{iclass}){nfi}.trials{it}.spikeIds = trial.spikeIds(idxsort);
                    unitSorted.class_responses.(unit.classes{iclass}){nfi}.trials{it}.spikeTimes = trial.spikeTimes(idxsort);
                    unitSorted.class_responses.(unit.classes{iclass}){nfi}.trials{it}.spikeWav = trial.spikeWav(idxsort,:);
                end
            end
        end
    end
    
    % Update the z-score calculation
    unitSorted = calc_zscore_unit(unitSorted, []);   % The empty array forces to use all data available
    
    % Write the new h5 file.
    write_unit_h5file(sprintf('%s_ss%d.h5', fname(1:end-3), isort), unitSorted);
    
end

set(handles.status, 'String', 'Saved'); drawnow

retry=0;
delete(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in peaks.
function peaks_Callback(hObject, eventdata, handles)
% hObject    handle to peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global PF vv ch currsort uvsorts vtrials vsorts labels sort validtrials tmpv tmpv2 subspace r5pc plotmin plotmax kmslast tmpsortidx tmpunsortedidx snips
set(handles.status, 'String', 'Wait'); drawnow
kms = str2num(get(handles.numcluster, 'String'));

if get(handles.randomforrest, 'Value') == 1
    set(handles.kmeanscluster, 'Value', 1)
    set(handles.spectralcluster, 'Value',0)
    set(handles.kmedoidscluster, 'Value',0)
    set(handles.randomforrest, 'Value',0)
end


if kms ~= kmslast
    tmpsortidx = zeros(1, size(snips,1));
    tmpunsortedidx = find(tmpsortidx==0);
    kmslast = kms;
end


subspace = 1;


if kms == 1
    labels = ones(size(snips,1),1);
else

    heights = max(snips(tmpunsortedidx,:)') - min(snips(tmpunsortedidx,:)');
    heights = heights - mean(heights);
    heights = heights ./ std(heights);
    [minr,minc] = find(bsxfun(@(x,y) (x == y), snips(tmpunsortedidx,:)', min(snips(tmpunsortedidx,:)'))==1);
    [maxr,maxc] = find(bsxfun(@(x,y) (x == y), snips(tmpunsortedidx,:)', max(snips(tmpunsortedidx,:)'))==1);
    minr = minr(unique(minc));
    maxr = maxr(unique(maxc));
    widths = minr - maxr;
    widths = widths - mean(widths);
    widths = widths ./ std(widths);
    dsnips = bsxfun(@minus, snips(tmpunsortedidx,:), msnips);
    ktmpv2 = r5pc(:,1:str2num(get(handles.numpcs, 'String')))'*dsnips';

    ktmpv = [widths'; heights; ktmpv2];


    labels = cluster(ktmpv', kms, handles);
end

for ll = unique(labels)'
    tmpsortidx(tmpunsortedidx(labels==ll)) = ll;
end

% Create the vector that indicates which snipets can't be sorted correctly
trashlabels = trashsnips(tmpsortidx, snips);

for kk = 1:6
    hh = ['handles.axes' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), snips((labels==kk),:)'); hold(eval(hh)); plot(eval(hh),mean(snips((labels==kk),:)',2), 'k', 'linewidth', 4);
    set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end

set(handles.status, 'String', 'Ready'); drawnow


function labels = cluster(data, k, handles)
global ch tmpv2 validtrials subspace sortparams
global heightsmean heightsstd widthsmean widthsstd ens r5pc s5pc kms tmpsortidx tmpunsortedidx clusize vs tridx tmps2 snips
global msnips mdsnips

if k ==1 | k ==0
    labels = ones(size(data,1),1);
    ens = 1;
else

    if get(handles.spectralcluster, 'Value') == 1
    % keyboard
        ens = [];
        W = SimGraph(data', 1);
        A = SpectralClustering(W, k, 3);
        labels = full(sum(bsxfun(@times, A, 1:kms),2));
    elseif get(handles.kmeanscluster, 'Value') == 1
        ens = [];
        if size(data,1)<size(data,2)
            labels = fkmeans(data,k);
        else
            try

                centers = zeros(k,size(data,2));
                for ii = 1:k
                    dsnips = bsxfun(@minus, snips(tmpsortidx==ii,:), msnips);
                    switch subspace
                        case 1
                            trn = snips(tmpsortidx==ii,:)';
                            heights = max(trn) - min(trn);
                            heights = heights - heightsmean;
                            heights = heights ./ heightsstd;
                            [minr,minc] = find(bsxfun(@(x,y) (x == y), trn, min(trn))==1);
                            [maxr,maxc] = find(bsxfun(@(x,y) (x == y), trn, max(trn))==1);
                            minr = minr(unique(minc));
                            maxr = maxr(unique(maxc));
                            widths = minr - maxr;
                            widths = widths - widthsmean;
                            widths = widths ./ widthsstd;
                            trn2 = r5pc(:,1:str2num(get(handles.numpcs, 'String')))'*dsnips';

                            trn = [widths'; heights; trn2];
                        case 2
                            trn = r5pc(:,1:str2num(get(handles.numpcs, 'String')))'*dsnips';
                    end

                    centers(ii,:)= mean(trn,2)';
                end

                labels = kmeans(data, k, 'start', centers);
            catch
                labels = kmeans(data, k);
            end
            
        end
    elseif get(handles.kmedoidscluster, 'Value') == 1
        ens = [];
        if size(data,1)==1
            labels=1;
        else
            labels = kmedoids(data', k)';
        end
    elseif get(handles.randomforrest, 'Value') == 1
        

            try
                
%                 tcv = cat(3, tmpclusters(validtrials).v);
%                 tcv = tcv(:,ch,:);
                cnt=0;
                n = sum(tmpsortidx>0);
                trnData = zeros(18+17+18+10+2,n);
                Y = zeros(1,n);
                for ii = 1:k

%                     trn = cat(2,tcv{ii,:});
                    trn = snips(tmpsortidx==ii,:)';
                    idx=cnt+1:cnt+size(trn,2);
                    trnData(1:18,idx) = trn;
                    trnData(19:35,idx) = diff(trn);
                    trnData(36:53,idx) = abs(fft(trn));
                    dsnips = bsxfun(@minus, snips(tmpsortidx==ii,:), msnips);
                    trnData(54:58,idx) = r5pc(:,1:5)'*dsnips';
                    dsnips = bsxfun(@minus, diff(snips(tmpsortidx==ii,:)')', mdsnips);
                    trnData(59:63,idx) = s5pc(:,1:5)'*dsnips';
                    heights = max(trn) - min(trn);
                    trnData(64,idx) = heights;
                    [minr,minc] = find(bsxfun(@(x,y) (x == y), trn, min(trn))==1);
                    [maxr,maxc] = find(bsxfun(@(x,y) (x == y), trn, max(trn))==1);
                    minr = minr(unique(minc));
                    maxr = maxr(unique(maxc));
                    widths = minr - maxr;
                    trnData(65,idx) = widths;
                    Y(idx) = ii;
                    cnt = cnt+size(trn,2);
                end

                ens = TreeBagger(300, trnData', Y');
                if ~isempty(data)
      
                    data = data';

                    X = zeros(18+17+18+10+2, size(data,2));
                    X(1:18,:) = data;
                    X(19:35,:) = diff(data);
                    X(36:53,:) = abs(fft(data));
                    dsnips = bsxfun(@minus, data', msnips);
                    X(54:58,:) = r5pc(:,1:5)'*dsnips';
                    dsnips = bsxfun(@minus, X(19:35,:)', mdsnips);
                    X(59:63,:) = s5pc(:,1:5)'*dsnips';
                    heights = max(data) - min(data);
                    X(64,:) = heights;
                    [minr,minc] = find(bsxfun(@(x,y) (x == y), data, min(data))==1);
                    [maxr,maxc] = find(bsxfun(@(x,y) (x == y), data, max(data))==1);
                    minr = minr(unique(minc));
                    maxr = maxr(unique(maxc));
                    widths = minr - maxr;
                    X(65,:) = widths;

                    [labels pobs]= predict(ens, X');
                    labels = str2num(cell2mat(labels));
                else
                    labels = [];
                end
            catch
                
                error('There is no training data for the RF')
                
            end
    end
end



% --- Executes on button press in kmeanscluster.
function kmeanscluster_Callback(hObject, eventdata, handles)
% hObject    handle to kmeanscluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kmeanscluster
set(hObject, 'Value', 1)
set(handles.spectralcluster, 'Value',0)
set(handles.kmedoidscluster, 'Value',0)
set(handles.randomforrest, 'Value',0)

% --- Executes on button press in spectralcluster.
function spectralcluster_Callback(hObject, eventdata, handles)
% hObject    handle to spectralcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spectralcluster
set(hObject, 'Value', 1)
set(handles.kmeanscluster, 'Value',0)
set(handles.kmedoidscluster, 'Value',0)
set(handles.randomforrest, 'Value',0)

% --- Executes on button press in kmedoidscluster.
function kmedoidscluster_Callback(hObject, eventdata, handles)
% hObject    handle to kmedoidscluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kmedoidscluster
set(hObject, 'Value', 1)
set(handles.spectralcluster, 'Value',0)
set(handles.kmeanscluster, 'Value',0)
set(handles.randomforrest, 'Value',0)


% --- Executes on button press in randomforrest.
function randomforrest_Callback(hObject, eventdata, handles)
% hObject    handle to randomforrest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of randomforrest
set(hObject, 'Value', 1)
set(handles.spectralcluster, 'Value',0)
set(handles.kmedoidscluster, 'Value',0)
set(handles.kmeanscluster, 'Value',0)


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in discard.
function discard_Callback(hObject, eventdata, handles)
% hObject    handle to discard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ens ch retry
ens = [];
retry=0;
delete(handles.figure1);


% --- Executes on button press in pcplot.
function pcplot_Callback(hObject, eventdata, handles)
% hObject    handle to pcplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ch tmpclusters validtrials subspace sortparams labels
global heightsmean heightsstd widthsmean widthsstd ens r5pc s5pc kms tmpv2 snips
figure;
colors = {'r', 'g', 'b', 'y', 'c', 'k'}

for kk = 1:kms
    dsnips = bsxfun(@minus, snips(labels==kk,:), msnips);
    tmp = r5pc(:,1:2)'*dsnips';
    hold on; scatter(tmp(1,:), tmp(2,:), colors{kk})
end


% --- Executes on button press in retry.
function retry_Callback(hObject, eventdata, handles)
% hObject    handle to retry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global retry ens
ens = [];
retry=1;
delete(handles.figure1);

function [plotmin, plotmax] = findplotminmax(tmpv2)

plotokflag = 0;
cnt = 0;

plotmin = min(tmpv2(:));
plotmax = max(tmpv2(:));

while plotokflag == 0
    

    [~,maxidx] = find(tmpv2 == max(tmpv2(:)));
    [~,minidx] = find(tmpv2 == min(tmpv2(:)));

    if plotmax >  max(4*std(tmpv2'))
        tmpv2(:,maxidx) = [];
        plotmin = min(tmpv2(:));
        plotmax = max(tmpv2(:));
    else
        plotokflag=1;
    end
    
    cnt = cnt+1;
    if cnt>10
        plotokfflag=1;
    end
end


% --- Executes on button press in multiunit1.
function multiunit1_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit1
set(handles.singleunit1, 'Value',0)
set(handles.noise1, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit1.
function singleunit1_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit1
set(handles.multiunit1, 'Value',0)
set(handles.noise1, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in multiunit2.
function multiunit2_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit2
set(handles.singleunit2, 'Value',0)
set(handles.noise2, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit2.
function singleunit2_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit2
set(handles.multiunit2, 'Value',0)
set(handles.noise2, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in multiunit3.
function multiunit3_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit3
set(handles.singleunit3, 'Value',0)
set(handles.noise3, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit3.
function singleunit3_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit3
set(handles.multiunit3, 'Value',0)
set(handles.noise3, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in multiunit4.
function multiunit4_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit4
set(handles.singleunit4, 'Value',0)
set(handles.noise4, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit4.
function singleunit4_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit4
set(handles.multiunit4, 'Value',0)
set(handles.noise4, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in multiunit5.
function multiunit5_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit5
set(handles.singleunit5, 'Value',0)
set(handles.noise5, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit5.
function singleunit5_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit5
set(handles.multiunit5, 'Value',0)
set(handles.noise5, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in multiunit6.
function multiunit6_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit6
set(handles.singleunit6, 'Value',0)
set(handles.noise6, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit6.
function singleunit6_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit6
set(handles.multiunit6, 'Value',0)
set(handles.noise6, 'Value',0)
set(hObject, 'Value',1)


% --- Executes on button press in noise1.
function noise1_Callback(hObject, eventdata, handles)
% hObject    handle to noise1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise1
set(handles.singleunit1, 'Value',0)
set(handles.multiunit1, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in noise2.
function noise2_Callback(hObject, eventdata, handles)
% hObject    handle to noise2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise2
set(handles.singleunit2, 'Value',0)
set(handles.multiunit2, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in noise3.
function noise3_Callback(hObject, eventdata, handles)
% hObject    handle to noise3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise3
set(handles.singleunit3, 'Value',0)
set(handles.multiunit3, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in noise4.
function noise4_Callback(hObject, eventdata, handles)
% hObject    handle to noise4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise4
set(handles.singleunit4, 'Value',0)
set(handles.multiunit4, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in noise5.
function noise5_Callback(hObject, eventdata, handles)
% hObject    handle to noise5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise5
set(handles.singleunit5, 'Value',0)
set(handles.multiunit5, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in noise6.
function noise6_Callback(hObject, eventdata, handles)
% hObject    handle to noise6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise6
set(handles.singleunit6, 'Value',0)
set(handles.multiunit6, 'Value',0)
set(hObject, 'Value',1)

% ---
function trashlabels = trashsnips(SnipIdx, snips, idxPresorted)
% Performx a chi-square test on the PCs of the snips to decide whether
% the sorthed snips belong to a cluster. Returns an index of the same size
% as SnipIds with one for trash spikes.


pThresh = 0.0001;        % Cut-off value for chi-square distribution.

% calclate the mean pf snips and sd of snips for each cluster
uniqueLabels = unique(SnipIdx);
trashlabels= zeros(size(SnipIdx));

if (nargin == 2)
    idxPresorted = SnipIdx;
end 

 
for kk=uniqueLabels

    LocalIdx = find(SnipIdx==kk);

    % Calculate PCs for presorted
    preSortedSnips = snips(idxPresorted == kk, :);
    [pcs coeff] = princomp(preSortedSnips);
    meanSnip = mean(preSortedSnips);
    meanCoeff = mean(coeff);
    sdCoeff= std(coeff);

      
    for ss=1:length(LocalIdx)
        Idx = LocalIdx(ss);
        coeffSS = (snips(Idx,:)-meanSnip)*pcs;

        chival = sum(((coeffSS - meanCoeff)./sdCoeff).^2);
        pval = 1- chi2cdf(chival, size(snips, 2));  
  
        if pval < pThresh
               trashlabels(Idx)=1;
        end
    end
end
