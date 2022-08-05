function offlineBRAutoSort(filename,multiplier,inputDataDirectory,outputPpnevDirectory)

% Spikes are detected using an envelope detector that is based on first
% calculating the amount of noise in the signal and applying a multiplier
% to that level. Then an envelope signal is calculated and compared to that
% calculated threshold.

if nargin < 2
     multiplier = 1.5;
end
swin  = -9:9;
win   = -30:30;

% [pathstr,name,~] = fileparts(filename); 
% filename = [pathstr filesep name]; 
% NS6_header = openNSx([filename '.ns6'],'noread');
cd(inputDataDirectory)
NS6_header = openNSx(filename,'noread');
ConnectorBank = {NS6_header.ElectrodesInfo.ConnectorBank};


cd(outputPpnevDirectory)
saveName = strcat(filename(1:end-4),'.mat');
empty = struct;
save(saveName,'empty','-v7.3')
writeObject = matfile(saveName,'Writable',true);

for e = 1:length(ConnectorBank)
    
    displayText = strcat('starting Electrode Number',string(e));
    disp(displayText)

    if strcmp(ConnectorBank{e},'E')
        continue
    end
    
    clear electrode NS x env noise
    
    cd(inputDataDirectory)
    electrode = sprintf('c:%u',e);
    NS = openNSx(filename, electrode,'read');
    if iscell(NS.Data)
        NS.Data =  cell2mat(NS.Data);
    end
    electrodeID = NS.ElectrodesInfo.ElectrodeID;
    x = double(NS.Data)';
    Fs = double(NS.MetaTags.SamplingFreq);
    nyq = Fs/2;
    
    clear NS

    if e == 1
        % Calculating the time vector
        clear TM
        TM = 1:length(x);
        TM = downsample(TM,2);
        TM = downsample(TM,3);
    end
    
    
    % - Second Order Butterworth Lowpass at 5kHz applied to the 30kHz raw signal
    lpc = 5000;
    lWn = lpc/nyq;
    [bwb,bwa] = butter(2,lWn,'low');
    env = filtfilt(bwb,bwa,x);
    clear bwb bwa
    
    % - Downsampled by 2
    env = downsample(env,2); 
    Fs  = Fs/2;
    nyq = Fs/2;
    
    % - Second Order Butterworth Highpass at 1kHz applied to the 15 kHz filtered signal
    hpc = 1000;
    hWn = hpc/nyq;
    [bwb,bwa] = butter(2,hWn,'high');
    env = filtfilt(bwb,bwa,env);
    clear bwb bwa

    % - Absolute value taken
    env = abs(env); 
    noise = env;
    
    
    % Calculating the Envelope:
    % - A 1.5kHz low pass second order butterworth is applied to the absolute bandpassed signal
    lpc = 1500;
    lWn = lpc/nyq;
    [bwb,bwa] = butter(2,lWn,'low');
    env = filtfilt(bwb,bwa,env);
    clear bwb bwa

    % - The signal is downsampled by 3 and is used as the envelope
    env = downsample(env,3);
    env = single(env); % Variable size is reduced due to memory issues (64GB ram on Lenovo ThinkPad memory limit reached by electrode 15)

    % Calculating the Noise:
    % - A boxcar smoother is applied and the signal [N/A is downsampled by 60, This gives a noise level which is only calculated every second]
    v = ones(Fs,1); % 1 second
    v = v ./ sum(v);
    noise = convnfft(noise,v,'same'); %can also use 'doConv' or 'conv' but slower
    
    % - The signal is downsampled by 3 to match envelope
    noise = downsample(noise,3);
    noise = single(noise); % Variable size is reduced due to memory issues (64GB ram on Lenovo ThinkPad memory limit reached by electrode 15)
    
    % - The noise level is multiplied by the noise multiplier (default 2.2)
    noise = noise .* multiplier;
        
    
    % Calculating spikes and extract spike waves;
    clear spk
    spk  = TM(env>noise);
    clear env noise 
    for i = 1:length(spk)
        w = spk(i) + swin;
        if any(w<1) || any(w>length(x))
            continue
        end
        dat     = x(spk(i) + swin);
        if all(dat == 0)
            continue
        end
        % align wave on max slope
            v       = diff(dat); 
            [~,mi]  = max(v);
            mi      = swin(mi);
            spk(i)  = spk(i) + mi;
            w       = spk(i) + win; 
            if any(w<1) || any(w>length(x))
                continue
            end
    end
    % get rid of redundant spikes (new 7/25/17)
    [spk, ui, ~] = unique(spk); 
    elec = zeros(1,length(spk)); 
    elec(:) = electrodeID;

    % store
    % store
    cd(outputPpnevDirectory)
    spkVarName = strcat('spk',string(e)); %formely TimeStamp
    writeObject.(spkVarName) = spk;
    elecVarName = strcat('elec',string(e)); %formerly Electrode
    writeObject.(elecVarName)  = elec; 

    


clear spk elec 
 
end
end

% % 
% % Unit = ones(size(TimeStamp));
% % [TimeStamp,idx] = sort(TimeStamp);
% % Electrode = Electrode(idx);
% % WAVES = WAVES(:,idx); 
% % 
% % % load NEV
% % clear NEV ppNEV
% % NEV = openNEV(strcat(filename,'.nev'),'noread','nomat','nosave');
% % 
% % % organize ppNEV structure for outpu
% % ppNEV.MetaTags             = NEV.MetaTags; 
% % ppNEV.ElectrodesInfo       = NEV.ElectrodesInfo; 
% % ppNEV.Data.SerialDigitalIO = NEV.Data.SerialDigitalIO; 
% % clear NEV
% % 
% % ppNEV.Data.Spikes.TimeStamp = TimeStamp; 
% % ppNEV.Data.Spikes.Electrode = Electrode; 
% % ppNEV.Data.Spikes.Unit = Unit; 
% % 
% % ppNEV.filename = filename; 
% % ppNEV.multiplier = multiplier; 
% % ppNEV.timestamp = now; 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
