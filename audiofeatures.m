% audiofeatures.m

% Jessica Thompson
% April 15, 2013

% use MIRToolbox to extract several audio features from musical stimuli
% 1. Spectral centroid
% 2. Amplitude envelope
% 3. Pulse Clarity
% 4. Mode
% 5. Key clarity

%cd stimuli
files = dir('stimuli/*.wav');
numfiles = length(files);
%test = 'stimuli/05-12-GraviesAndTheMainDishSauce-BeefBallsPoachedInSpinachSoup.wav';
numtracks=[83,10,8,9,12,14];
trackidx=zeros(sum(numtracks),2);
m=1;
for i=1:6
    trackidx(m:m+numtracks(i)-1,1) = i;
    trackidx(m:m+numtracks(i)-1,2) = 1:numtracks(i);
    m=m+numtracks(i);
end

for i=84:numfiles
    fname = ['stimuli/' files(i).name];
    %a = miraudio(fname, 'Frame',8192,'sp', 4410,'sp');
    % Spectral centroid
    s = mirspectrum(fname, 'Frame',8192,'sp', 4410,'sp', 'Window','hamming','Max',16000,'dB'); % 100 ms windows, half overlap, dB scaling, 
    c= mircentroid(s); 
    % Amplitude envelope
    e = mirenvelope(fname,'Halfwave','PostDecim', 44100/10, 'Gauss', 2);
    % Pulse Clarity
    ac=mirautocor(fname, 'Frame',8192,'sp', 4410,'sp', 'Window','hamming','Max',16000);
    p = mirpulseclarity(ac);
    % Mode
    m = mirmode(s);
    % Key clarity
    [k,kc] = mirkey(s);
    % store
    session = trackidx(i,1);
    track = trackidx(i,2);
    features(session,track).centroid = mirgetdata(c);
    features(session,track).env = mirgetdata(e);
    features(session,track).pulseclar = mirgetdata(p);
    features(session,track).mode = mirgetdata(m);
    features(session,track).keyclar = mirgetdata(kc);
end
save('audiofeatures.mat','features')