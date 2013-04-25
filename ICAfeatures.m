% ICAfeatures.m
% 
% Jessica Thompaon
% April 15, 2013

% Calculate EEG features based on ICA stability, dipolarity, and
% correlation with musical audio feature
% 

clear
sessions = {'08feb13cd', '15feb13cd', '22feb13cd', '26feb13cd', '27feb13cd', '01mar13cd'};
%sessions = {'15feb13cd', '22feb13cd', '26feb13cd', '27feb13cd', '01mar13cd'};
srate=200;

%load ICASSOresults % contains sR, Iq, in_avg, ext_avg, in_min, ext_max
%% Stability Analysis
disp('Beginning stability analysis...')

%[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % start EEGLAB
stable_comps = zeros(64,6); % indicies of stable components
rdv = ones(64,6); % Residual dipole variance
% 
%load('stabilityresults.mat')
for i=1:6
    % Load EEGLAB datasets, ICA and dipole fitting must have been performed
    EEG{i} = pop_loadset('filename',[sessions{i},'_filt+resamp200+music+ICA.set'],'filepath','/Users/jthompson/data/EEG-Thesis/');
end
for i=1:6
    disp(sprintf('session %d',i))
    % Load EEGLAB datasets, ICA and dipole fitting must have been performed
    EEG{i} = pop_loadset('filename',[sessions{i},'_filt+resamp200+music+ICA.set'],'filepath','/Users/jthompson/data/EEG-Thesis/');
    % do not reduce dimensionality
    % Resample 30 times using random initial conditions and bootstrapping
    % (use FastICA parameters: kurtosis as non-linearity,
    % symmetric estimation approach
    sR{i}=icassoEst('both', EEG{i}.data, 30, 'lastEig', 64, 'g', 'pow3', ...
                 'approach', 'symm','maxNumIterations', 150); %perform ICA 30 times
    sR{i} = icassoCluster(sR{i}); %cluster
    [Iq{i},in_avg{i},ext_avg{i},in_min{i},ext_max{i}]=icassoStability(sR{i},64); % calculate stability index

    % Reject components with IQ less than .97
    stable_comps(:,i) = Iq{i}>=.96;
    %rdv(:,i) = getfield(EEG{i}, {'dipfit' 'model' 'rv'}, ':');
    for j=1:64
        rdv(j,i) = EEG{i}.dipfit.model(1,j).rv;      
    end      
 end

rdv20= rdv<.2;% indicies of dipolar components (redidual dipole variance less than 20%)
accepted_comps = stable_comps & rdv20; % accept only stable, dipolar components
accepted_spatial_maps = zeros(sum(accepted_comps(:)),64);% each row is a spatial map (components by dimensions)

m=1;
for k=1:6
    numcomps = sum(accepted_comps(:,k));
    accepted_spatial_maps(m:m+numcomps-1,:) = EEG{k}.icawinv(:,accepted_comps(:,k))';
    m=m+numcomps;
end
save('stabilityresults.mat','accepted_comps', 'accepted_spatial_maps', 'sR')
remaining = sum(accepted_comps(:));
rejected = (64*4)-remaining;
disp(sprintf('%d components remain. (%d rejected by stability and dipolarity analysis.)', remaining, rejected))
% OR
%load('stabilityresults.mat')

%% cluster spatial maps of accepted components
%[CENTRES, OPTIONS, POST, ERRLOG] = kmeans(CENTRES, DATA);
disp('Beginning clustering of components...')
maxclust = 7;
X = accepted_spatial_maps;
C = clusterdata(X,'maxclust',maxclust);
%scatter3(X(:,1),X(:,2),X(:,3),100,T,'filled')

sessionidx = [];
compidx=[];
for j=1:6
    sessionidx= [sessionidx j*ones(1, sum(accepted_comps(:,j)))];
    compidx = [compidx find(accepted_comps(:,j))'];
end
idx = cat(1, sessionidx, compidx)'; %transposed 
%run(mainClusterValidationNC) % 3 or 7 clusters
% reject clusters that contain less than half of the sessions
accepted_clusters = [];
for i=1:maxclust
    %clusters{cluster}(session,comp)
    clustidx = find(C==i);
    if sum(ismember([1 2 3 4 5 6], idx(clustidx,1)))>2 % if half or more of the session are present in cluster i, accept
        clusters{i} = idx(clustidx,:);
        accepted_clusters = [accepted_clusters i];
    end   
end
save('acceptedclusters.mat','clusters', 'accepted_clusters')
numclusters = length(accepted_clusters);
disp(sprintf('%d out of %d clusters accepted.',numclusters, maxclust));
% OR
%load('acceptedclusters.mat')
%% extract component oscilation bands
% get start and stop latencies, need to look at type to not include
% boundary events
disp('Beginning time frequency analysis of accepted cluster components...')
for i=1:6
    %EEG{i} = pop_loadset('filename',[sessions{i},'_filt+resamp200+music+ICA.set'],'filepath','/Users/jthompson/data/EEG-Thesis/');%load EEG data
    start{i}= [];
    stop{i} = [];
    for j=1:size(EEG{i}.event, 2)
        if str2num(EEG{i}.event(1,j).type) < 150
            if mod(str2num(EEG{i}.event(1,j).type),2) == 1
                start{i} = [start{i} EEG{i}.event(1,j).latency];
            elseif mod(str2num(EEG{i}.event(1,j).type),2) == 0
                stop{i} = [stop{i} EEG{i}.event(1,j).latency];
            end
        else
            if mod(str2num(EEG{i}.event(1,j).type),2) == 1
                stop{i} = [stop{i} EEG{i}.event(1,j).latency];
            elseif mod(str2num(EEG{i}.event(1,j).type),2) == 0
                start{i} = [start{i} EEG{i}.event(1,j).latency];
            end
        end
    end
end
stop{1} = [stop{1} size(EEG{1}.data,2)]; 
% extract alpha (8--12Hz), theta (4--7Hz), and beta (12--30) time courses
% for each component in each accepted cluster for each stimulus
for i=1:length(accepted_clusters) %for each cluster
    clustidx = clusters{accepted_clusters(i)};
    for j=1:size(clustidx,1) %for each component
        session= clustidx(j,1);
        comp= clustidx(j,2);
        for k=1:length(start{session}) %for each stimulus
            tmpsig = eeg_getdatact(EEG{session}, 'component', comp, 'samples', [int32(start{session}(k)):int32(stop{session}(k))]);
            %S=spectrogram(x,window,noverlap,nfft,fs)
            [S,F,T,P] = spectrogram(double(tmpsig),40,20,256,200);% 20th of sec window, half overlap (like the audio features)
            %surf(T,F,10*log10(P),'edgecolor','none'); axis tight;
            %view(0,90);
            thetaidx1=find(diff(F>4));
            thetaidx2=find(diff(F>7))+1;
            alphaidx1=find(diff(F>8));
            alphaidx2=find(diff(F>12))+1;
            betaidx1=find(diff(F>12));
            betaidx2=find(diff(F>30))+1; 
            theta = mean(P(thetaidx1:thetaidx2,:),1);
            alpha = mean(P(alphaidx1:alphaidx2,:),1);
            beta = mean(P(betaidx1:betaidx2,:),1);
            oscillations{i}{j}{k}.alpha = alpha;
            oscillations{i}{j}{k}.beta = beta;
            oscillations{i}{j}{k}.theta = theta;
            fname = sprintf('oscillations/%02d-%02d_clust%02d_comp%02d.mat',session-1, k-1,accepted_clusters(i),comp);
            save(fname, 'alpha','beta','theta')
            
        end
    end
end
save('alloscillations.mat','oscillations','start','stop')
load('alloscillations.mat')

%% Correlation Analysis
% Correlate 45s of the three eeg freq bands from each session with 5 audiofeatures from Cong paper, extracted with
% MIRToolobox
%load('audiofeatures.mat')
numtracks=[83,10,8,9,12,14];
featuresrate=10;
%

len = (30*featuresrate); %30 seconds
for session=1:6 % for each session
    track = random('unid',numtracks(session)); % choose a random track
    tracklen = length(features{session}(track).centroid);
    %tracklen = (stop{session}(track)-start{session}(track));
    while tracklen<len%length of track less than 30 seconds
        track = random('unid',numtracks(session));
        tracklen = length(features{session}(track).centroid);
    end
    start = random('unid',int32(tracklen-len));% randomly select starting point 45 seconds or
    stop = start+len;
    feature_selection_set_idx{session} = [track,start,stop];
end

% end
for i=1:length(accepted_clusters)%for each cluster
    %for each component
    clustidx = accepted_clusters(i);
    numcomps(i) = size(clusters{clustidx},1);
    %corrcoefs{i} = zeros(numcomps,3,5);
    for j=1:numcomps(i)
       session = clusters{clustidx}(j,1);
       comp = clusters{clustidx}(j,2);
       trck = feature_selection_set_idx{session}(1);
       strt = feature_selection_set_idx{session}(2);
       stp = feature_selection_set_idx{session}(3);
       a = oscillations{i}{j}{trck}.alpha(strt:stp);
       b = oscillations{i}{j}{trck}.beta(strt:stp);
       t = oscillations{i}{j}{trck}.theta(strt:stp);
       osc = {a, b, t};
       c = features{session}(trck).centroid(strt:stp);
       e = features{session}(trck).envresamp(strt:stp);
       p = features{session}(trck).pulseclar(strt:stp);
       m = features{session}(trck).mode(strt:stp);
       k = features{session}(trck).keyclar(strt:stp);
       audio = {c,e,p,m,k};
       for audiofeat=1:5
           for eegosc=1:3
               [R,P] = corrcoef(osc{eegosc}',audio{audiofeat}');
               corrcoefs{i}(j,audiofeat,eegosc,:)= [R(1,2),P(1,2)];
           end
       end
    end
end
% accept only clusters where half of components are significantly
% correlated with at least one of the musical audio features
for i=1:length(accepted_clusters)
    if sum(sum(sum(sum(corrcoefs{i}(:,:,:,2)<.05,4),3),2)>0) > numcomps(i)/2
        disp(sprintf('cluster %s accepted!'),accepted_clusters(i))
    else
        accepted_clusters(i) = [];
    end
end

%% Calcuate average spatial maps

for clust=1:length(accepted_clusters)
    num_comps = length(clusters{accepted_clusters(clust)});
    spatial_maps_2avg{clust} = zeros(num_comps,64);
    m=1;
    for session=1:6
       compindx = find(clusters{accepted_clusters(clust)}(:,1)==session);
       comps = clusters{accepted_clusters(clust)}(compindx,2);
       spatial_maps_2avg{clust}(m:m+length(comps)-1,:) = EEG{session}.icawinv(:,comps)'; %(channel x components)
       m=m+length(comps);
    end
    spatialmap(clust,:) = mean(spatial_maps_2avg{clust},1);
end
save('spatialmaps.mat','spatialmap')

%% calculate frequency band activity and save as final features
featvecspertrack = zeros(136,1);
trackidx = 0;
for session=1:6
    clear theta alpha beta
    for clust=1:2
        eegfeatsig{session}(clust,:) = spatialmap(clust,:)*EEG{session}.data;
        m=1;
        for k=1:length(start{session}) % for each track
            
            [S,F,T,P] = spectrogram(double(eegfeatsig{session}(clust,int32(start{session}(k)):int32(stop{session}(k)))),40,20,256,200);
            featlen=length(T);
            if clust==1
                trackidx = trackidx+1;
                featvecspertrack(trackidx) = featlen;
            end
            thetaidx1=find(diff(F>4));
            thetaidx2=find(diff(F>7))+1;
            alphaidx1=find(diff(F>8));
            alphaidx2=find(diff(F>12))+1;
            betaidx1=find(diff(F>12));
            betaidx2=find(diff(F>30))+1;
            theta(clust,m:m+featlen-1) = mean(P(thetaidx1:thetaidx2,:),1);
            alpha(clust,m:m+featlen-1) = mean(P(alphaidx1:alphaidx2,:),1);
            beta(clust,m:m+featlen-1) = mean(P(betaidx1:betaidx2,:),1);
            m=m+featlen;
        end
    end
    eegfeatures{session} = [theta' alpha' beta'];
end
save('eegfeatures.mat','eegfeatures')
csvwrite('featvecspertrack.csv',featvecspertrack)
feats2save = [];
for i=1:6
    feats2save = cat(1,feats2save, eegfeatures{i});
end
csvwrite('eegfeatures.csv',feats2save)

% resample env feature to be same length as everything else
% for i=1:6
%     for j=1:numtracks(i)
%         features{i}(j).envresamp = resample(features{i}(j).env, length(features{i}(j).centroid),length(features{i}(j).env));
%     end
% end

% Plot accepted averaged spatial maps
subplot(121)
topoplot(spatialmap(1,:), EEG{1}.chanlocs);
title('Cluster 1','FontSize',18)
subplot(122)
topoplot(spatialmap(2,:), EEG{1}.chanlocs);
title('Cluster 2','FontSize',18)

