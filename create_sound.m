function sound_sample = create_sound(instrument,notes,constants)
% load constants
fs = constants.fs;
durationScale = constants.durationScale;
durationChord = constants.durationChord;

% find frequency of notes
freq_vec = zeros();
switch instrument.temperament
    case {'Equal'}
        if length(notes) > 1
            for ii = 1:length(notes)
                freq_vec(ii) = note_to_freq(notes{ii}.note);
%                 root = notes{ii}.note;
%                 if length(root) == 1
%                     n_keys = 48 + root - 'A';
%                 elseif length(root) == 2
%                     oct = str2double(root(2));
%                     n_keys = 12*oct + (root(1) - 'A') + 1;
%                 else
%                     error('Inproper root specified');
%                 end
%                 freq_vec(ii) = 2^((n_keys - 49)/12) * 440;
            end
        else
            root = notes.note;
            if length(root) == 1
                n_keys = 48 + root - 'A';
            elseif length(root) == 2
                oct = str2double(root(2));
                n_keys = 12*oct + (root(1) - 'A') + 1;
            else
                error('Inproper root specified');
            end
            freq_vec = 2^((n_keys - 49)/12) * 440;
        end
    case {'Just'}
%         if length(root) == 1
%             n_keys = 48 + root - 'A';
%         elseif length(root) == 2
%             oct = str2double(root(2));
%             n_keys = 12*oct + (root(1) - 'A') + 1;
%         else
%             error('Inproper root specified');
%         end
%         root_freq = 2^((n_keys - 49)/12) * 440;
        if length(notes) > 1
            root = notes{1}.note;
            root_freq = note_to_freq(root);
            freq_vec(1) = root_freq;
            switch instrument.mode
                case {'Major'}
                    ratio_prev = cumprod([9/8*10/9 16/15*9/8]);
                case {'Minor'}
                    ratio_prev = cumprod([9/8*16/15 10/9*9/8]);
            end
            ratios = [1 ratio_prev];
            freq_vec = root_freq*ratios;
        else
            root = notes.note;
            root_freq = note_to_freq(root);
            freq_vec = root_freq;
        end
end

t = 0:1/fs:(instrument.totalTime-1/fs);
switch instrument.sound
    case {'Additive'} % bell from Jerse 4.28
        amp_mult = [1 .67 1 1.8 2.67 1.67 1.46 1.33 1.33 1 1.33];
        dur_mult = [1 .9 .65 .55 .325 .35 .25 .2 .15 .1 .075];
        freq_mult = [.56 .56 .92 .92 1.19 1.7 2 2.74 3 3.76 4.07];
        freq_add = [0 1 0 1.7 0 0 0 0 0 0 0];
        
        tones = zeros(length(notes),length(t));
        for ii = 1:length(notes)
            if length(notes) > 1
                dur_vec = notes{ii}.duration*dur_mult.';
            else
                dur_vec = notes.duration*dur_mult.';
            end
            t_mat = repmat(t,length(dur_vec),1);
            envelope = zeros(length(dur_vec),length(t));
            for jj = 1:length(dur_vec)
                
                if dur_vec(jj) < length(t)
                    t_mat(jj,int64(dur_vec(jj))+1:end) = 0;
                    envelope(jj,:) = [ones(1,int64(dur_vec(jj)/2)) ...
                        linspace(1,0,int64(dur_vec(jj)/2)) zeros(1,length(t)-int64(dur_vec(jj)))];
                else
                    envelope(jj,:) = [ones(1,int64(dur_vec(jj)/2)) ...
                        linspace(1,0,int64(dur_vec(jj)/2))];
                end
%                 t_comp = 0:1/fs:(dur_vec(jj)-1/fs);
%                 t_note = [t_comp zeros(1,length(t)-length(t_comp))];
            end
            freq_mat = repmat((freq_vec(ii)*freq_mult + freq_add).',1,length(t));
            amp_mat = repmat(amp_mult.',1,length(t));
            tones(ii,:) = sum(amp_mat.*sin(2*pi*freq_mat.*t_mat).*envelope);
        end
        sound_sample = sum(tones,1);
%         fund = sum(sin(2*pi*freq_vec.*repmat(t,length(freq_vec),1))/2.5);
%         inharm = sum(sin(2*pi*(freq_vec/10).*repmat(t,length(freq_vec),1))/6);
%         noise = sin(2*pi*500*normrnd(1)*repmat(t,length(freq_vec),1))
    case {'Subtractive'}
        vbw = dsp.VariableBandwidthFIRFilter('FilterType','Bandpass',...
            'FilterOrder',500,'SampleRate',fs,'CenterFrequency',440,...
            'Bandwidth',440/4,'Window','Chebyshev');
        
        tones = zeros(length(notes),length(t));
        for ii = 1:length(notes)
            vbw.CenterFrequency = freq_vec(ii);
            vbw.Bandwidth = vbw.CenterFrequency/16;
            freq_start = 0.55*freq_vec(ii);
            freq_end = 1.05*freq_vec(ii);
            source = square(2*pi*100.*t);
            
            p = 200;
            cent_vec = linspace(freq_start,freq_end,p);
            for jj = 1:p
                vbw.CenterFrequency = cent_vec(jj);
                tones(ii,((jj-1)*length(t)/p+1):jj*length(t)/p) = ...
                    vbw(source(((jj-1)*length(t)/p+1):jj*length(t)/p).');
            end
        end
        attack = linspace(0,1,length(t)/2);
        envelope = [attack ones(1,length(t)-length(t)/16-length(attack)) linspace(1,0,length(t)/16)];
        tones = tones.*repmat(envelope,length(notes),1);
        sound_sample = sum(tones,1);
        
%         bw = 100;
%         b2 = exp(-2*pi*bw/fs);
%         
%         attack = linspace(0,1,length(t)/2);
%         envelope = [attack ones(1,length(t)-length(t)/8-length(attack)) linspace(1,0,length(t)/8)];
%         
%         tones = zeros(length(notes),length(t)+2);
%         tp2 = 0:1/fs:(instrument.totalTime+1/fs);
%         for ii = 1:length(notes)
%             b1 = -4*b2/(1+b2)*cos(2*pi*(linspace(freq_vec(ii)*0.7,freq_vec(ii)*1.1,length(t)))/fs);
%             a0 = (1-b2)*sqrt(1-b1.^2/(4*b2));
%             source = sawtooth(2*pi*freq_vec(ii).*tp2);
%             for jj = 1:length(t)
%                 tones(ii,jj+2) = [a0(jj) -b1(jj) -b2]*[source(jj+2);tones(ii,jj+1);tones(ii,jj)];
%             end
%         end
%         tones = tones(:,3:end).*repmat(envelope,length(notes),1);
%         sound_sample = sum(tones,1);

     case {'FM'} % clarinet from Jerse 5.10
         IMAX = 0.001;
         f1_env = [((1:length(t)/4)/(length(t)/4)).^2 ones(1,(7*length(t)/8)-(length(t)/4)) ...
             fliplr(((1:length(t)/8)/(length(t)/8)).^2)];
         f2_env = [fliplr(((1:length(t)/4)/(length(t)/4)).^2) zeros(1,3*length(t)/4)];
         
         tones = zeros(length(notes),length(t));
         for ii = 1:length(freq_vec)
             fc = freq_vec(ii);
             fm = 2/3*fc;
             IMAX = IMAX;
             IMIN = IMAX/2;
             
             mod_amp = fm*(IMAX-IMIN)*f2_env + fm*IMIN;
             mod_sig = mod_amp.*sin(2*pi*fm*t);
             tones(ii,:) = f1_env.*sin(2*pi*(fc.*t + mod_sig));
             if length(notes) > 1
                 tones(ii,:) = [tones(ii,1:notes{ii}.duration) zeros(1,(notes{ii}.duration+1):length(t))];
             else
                 tones(ii,:) = [tones(ii,1:notes.duration) zeros(1,(notes.duration+1):length(t))];
             end
         end
         sound_sample = sum(tones,1);
    case {'Waveshaper'} % clarinet from Jerse 5.28
        attack = 0.085;
        decay = 0.025;
        release = 0.64;
        envelope = [linspace(0,1,attack*length(t)) linspace(1,0.75,decay*length(t)) ...
            0.75*ones(1,(1-attack-decay-release)*length(t)) linspace(0.75,0,release*length(t))];
        
        tones = zeros(length(notes),length(t));
        for ii = 1:length(notes)
            source_sig = 255*sin(2*pi*freq_vec(ii).*t); % generate sinewave
            env_sig = source_sig.*envelope + 256; % apply envelope
            tones(ii,:) = wave_shape(env_sig); % wave-shaping
        end
        sound_sample = sum(tones,1);
end

end

function output = wave_shape(input)
    output = zeros(size(input));
    tf_left = @(x) (0.5/(511-312))*x - 1; % transfer function
    tf_center = @(x) (0.5/(312-255))*(x-255);
    tf_right = @(x) (0.5/200)*(x-511) + 1;
    
    intermed = [tf_left(input);tf_center(input);tf_right(input)]; % apply transfer function
    output(input <= 200) = intermed(1,(input <= 200));
    output(input > 200 & input < 312) = intermed(2,(input > 200 & input < 312));
    output(input >= 312) = intermed(3,(input >= 312));
end

function freq = note_to_freq(note) % look up frequency
    note_table = readtable('note_freq.csv','ReadRowNames',1);
    freq = note_table(note,1).Hz;
end