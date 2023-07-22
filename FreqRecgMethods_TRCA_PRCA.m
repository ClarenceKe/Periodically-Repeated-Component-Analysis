function [acc, itrr]=FreqRecgMethods_TRCA_PRCA(eeg, labels, num_targs, ...
    dataLength, train_data_l, sti_f, sti_phase, num_cv, cv_label, delay_s, time, fs, num_fbs, is_ensemble, num_of_harmonics)

for loocv_i = 1:1:num_cv
    % PRCA
    segment_data_fixDL=time>=(0) & time<=(delay_s+train_data_l*0.1);
    traindata_fixDL = eeg(:, :, segment_data_fixDL, cv_label(loocv_i, :));
    model.prca_fixDL = train_prca_trca(traindata_fixDL, fs, delay_s, num_fbs, sti_f, num_of_harmonics); % fixed training data length for PRCA 

    for n_l=1:length(dataLength)%1:1:30
        data_l=dataLength(n_l);
    
        segment_data=time>=(0) & time<=(delay_s+data_l*0.1);
        len_sel_s=0.5+data_l*0.1;

        % Training stage 
        traindata = eeg(:, :, segment_data, cv_label(loocv_i, :));

        model.prca_trca = train_prca_trca(traindata, fs, delay_s, num_fbs, sti_f, num_of_harmonics); % trca and prca: training data length = training data length 

        test = 1:1:size(eeg, 4); test(cv_label(loocv_i, :)) = [];
        testdata = squeeze(eeg(:, :, segment_data, test));

        for nblock_test=1:size(testdata, 4)
            estimated.prca_trca= test_prca_trca(testdata(:, :, :, nblock_test), model.prca_trca, delay_s, is_ensemble); % trca and prca: training data length = training data length 
            estimated.prca_fixDL= test_prca(testdata(:, :, :, nblock_test), model.prca_fixDL, delay_s, is_ensemble);% fixed training data length for PRCA 

            tmp_fd=fields(estimated);
            for nfd=1:length(tmp_fd)
                tmp_fd2=fields(estimated.(tmp_fd{nfd}));
                for nfd2=1:length(tmp_fd2)
                    tmp_is_correct.(tmp_fd{nfd}).(tmp_fd2{nfd2})(:, nblock_test)=estimated.(tmp_fd{nfd}).(tmp_fd2{nfd2})==labels;
                end
            end
        end

        tmp_fd=fields(estimated);
        for nfd=1:length(tmp_fd)
            tmp_fd2=fields(estimated.(tmp_fd{nfd}));
            for nfd2=1:length(tmp_fd2)
                accs.(tmp_fd{nfd}).(tmp_fd2{nfd2})(n_l, loocv_i) = ...
                    mean(tmp_is_correct.(tmp_fd{nfd}).(tmp_fd2{nfd2}), 'all')*100;
                itrs.(tmp_fd{nfd}).(tmp_fd2{nfd2})(n_l, loocv_i) = ...
                    itr(num_targs, mean(tmp_is_correct.(tmp_fd{nfd}).(tmp_fd2{nfd2}), 'all'), len_sel_s);
            end
        end

        fprintf('\r PRCA: # CV: %d, DL = %2.2f s: Acc = %2.2f%%, ITR = %2.2f bpm',...
            loocv_i, data_l*0.1, accs.prca_trca.prca(n_l, loocv_i), itrs.prca_trca.prca(n_l, loocv_i));
        fprintf('\r PRCA_fixDL: # CV: %d, DL = %2.2f s: Acc = %2.2f%%, ITR = %2.2f bpm',...
            loocv_i, data_l*0.1, accs.prca_fixDL.prca(n_l, loocv_i), itrs.prca_fixDL.prca(n_l, loocv_i));
        fprintf('\r TRCA: # CV: %d, DL = %2.2f s: Acc = %2.2f%%, ITR = %2.2f bpm',...
            loocv_i, data_l*0.1, accs.prca_trca.trca(n_l, loocv_i), itrs.prca_trca.trca(n_l, loocv_i));
    end
end% loocv_i

tmp_fd=fields(estimated);
for nfd=1:length(tmp_fd)
    tmp_fd2=fields(estimated.(tmp_fd{nfd}));
    for nfd2=1:length(tmp_fd2)
        acc.(tmp_fd{nfd}).(tmp_fd2{nfd2})=mean(accs.(tmp_fd{nfd}).(tmp_fd2{nfd2}), [2]);
        itrr.(tmp_fd{nfd}).(tmp_fd2{nfd2})=mean(itrs.(tmp_fd{nfd}).(tmp_fd2{nfd2}), [2]);
    end
end

%%
function model = train_prca_trca(eeg, fs, delay_s, num_fbs, sti_f, num_of_harmonics)
if nargin < 6
    error('stats:train_trca:LackOfInput', 'Not enough input arguments.'); 
end

if ~exist('num_fbs', 'var') || isempty(num_fbs), num_fbs = 3; end
[num_targs, num_chans, num_smpls, ~] = size(eeg);
%%
for targ_i = 1:1:num_targs
    for fb_i = 1:1:num_fbs
        eeg_tmp0 = squeeze(eeg(targ_i, :, :, :));
        eeg_tmp0 = filterbank(eeg_tmp0(:,:,:), fs, fb_i);
        eeg_tmp = eeg_tmp0(:,round(delay_s*fs):end,:);

        trains_spatial.Tt(targ_i,fb_i,:,:) = squeeze(mean(eeg_tmp,3));

        [w_tmp, trains_spatial.Pt{targ_i}(fb_i,:,:)]=prca_trca(eeg_tmp, sti_f(targ_i), fs);

        tmp_fd=fields(w_tmp);
        for nfd=1:length(tmp_fd)
            W.(tmp_fd{nfd})(fb_i, targ_i, :)=w_tmp.(tmp_fd{nfd})(:,1);
        end
    end
end % fb_i

model = struct('trains_spatial', trains_spatial, 'W', W, ...
    'num_fbs', num_fbs, 'fs', fs, 'num_targs', num_targs, 'sti_f', sti_f, ...
    'num_of_harmonics', num_of_harmonics);

%%
function [W, tmplate_P]=  prca_trca(eeg, sti_f, fs)
%%
[num_chans, num_smpls, num_trials]  = size(eeg);
Nsamp_T=round(fs/sti_f);
N_T=floor(num_smpls/Nsamp_T);
N_T4tmplate=1;

for nt=1:num_trials
    for n=1:N_T-(N_T4tmplate-1)
        tmp=eeg(:, (1:Nsamp_T*N_T4tmplate)+Nsamp_T*(n-1), nt); 
        tmp_eeg(:, :, n, nt)=bsxfun(@minus, tmp, mean(tmp,2));
    end
    eeg(:, :, nt)=bsxfun(@minus, eeg(:, :, nt), mean(eeg(:, :, nt),2));
end

%%  PRCA
X1 = tmp_eeg(:,:);
X2 = sum(tmp_eeg, [3 4]);
Q = X1*X1'/size(X1,2);
S = X2*X2'/size(X2,2)-Q;
% PRCA eigenvalue algorithm
[eig_vec,eig_val] = eig(Q\S);
[V,sort_idx]=sort(diag(eig_val),'descend');
W.prca=eig_vec(:,sort_idx);

%% TRCA
X1 = eeg(:,:);
X2 = sum(eeg, 3);
Q = X1*X1'/size(X1,2);
S = X2*X2'/size(X2,2)-Q;
% TRCA eigenvalue algorithm
[eig_vec,eig_val] = eig(Q\S);
[V,sort_idx]=sort(diag(eig_val),'descend');
W.trca=eig_vec(:,sort_idx);
%%
tmplate_P=squeeze(mean(tmp_eeg, [3 4]));

%%
function [results]= test_prca_trca(eeg, model, delay_s, is_ensemble)

global a; global b;

if ~exist('is_ensemble', 'var') || isempty(is_ensemble)
    is_ensemble = 1; end

if ~exist('model', 'var')
    error('Training model based on TRCA is required. See train_trca().'); 
end

%fb_coefs = [1:model.num_fbs].^(-1.25)+0.25;

fb_coefs = [1:model.num_fbs].^(-a)+b;

% tmp_fd=fields(model.W);

for targ_i = 1:1:model.num_targs
    test_tmp = squeeze(eeg(targ_i, :, :));
    tmp_fd=fields(model.W);
    for fb_i = 1:1:model.num_fbs
        testdata = filterbank(test_tmp, model.fs, fb_i);
        testdata_tmp=testdata(:, round(delay_s*model.fs):end);
        testdata_tmp=bsxfun(@minus, testdata_tmp, mean(testdata_tmp,2));

        for class_i = 1:1:model.num_targs
            if ~is_ensemble
                for nfd=1:length(tmp_fd)
                    w.(tmp_fd{nfd})=real(squeeze(model.W.(tmp_fd{nfd})(fb_i, class_i, :)));
                end
            else
                for nfd=1:length(tmp_fd)
                    w.(tmp_fd{nfd})=real(squeeze(model.W.(tmp_fd{nfd})(fb_i, :, :))');
                end
            end

            traindata_P =  squeeze(model.trains_spatial.Pt{class_i}( fb_i, :, :));
            traindata_T =  squeeze(model.trains_spatial.Tt(class_i, fb_i, :, 1:size(testdata_tmp, 2)));

            r.prca(fb_i,class_i) =corr2(testdata_tmp'*w.prca, ...
                templatemaker(testdata_tmp, traindata_P'*w.prca));

            r.trca(fb_i,class_i) =corr2(testdata_tmp'*w.trca, ...
                traindata_T'*w.trca);
        end % class_i
    end % fb_i

    [~, results.prca(targ_i)] = max(fb_coefs*r.prca);
    [~, results.trca(targ_i)] = max(fb_coefs*r.trca);
end % targ_i

%%
function [results]= test_prca(eeg, model, delay_s, is_ensemble)

global a; global b;

if ~exist('is_ensemble', 'var') || isempty(is_ensemble)
    is_ensemble = 1; end

if ~exist('model', 'var')
    error('Training model based on TRCA is required. See train_trca().'); 
end

%fb_coefs = [1:model.num_fbs].^(-1.25)+0.25;

fb_coefs = [1:model.num_fbs].^(-a)+b;

% tmp_fd=fields(model.W);

for targ_i = 1:1:model.num_targs
    test_tmp = squeeze(eeg(targ_i, :, :));
    tmp_fd=fields(model.W);
    for fb_i = 1:1:model.num_fbs
        testdata = filterbank(test_tmp, model.fs, fb_i);
        testdata_tmp=testdata(:, round(delay_s*model.fs):end);
        testdata_tmp=bsxfun(@minus, testdata_tmp, mean(testdata_tmp,2));

        for class_i = 1:1:model.num_targs
            if ~is_ensemble
                for nfd=1:length(tmp_fd)
                    w.(tmp_fd{nfd})=real(squeeze(model.W.(tmp_fd{nfd})(fb_i, class_i, :)));
                end
            else
                for nfd=1:length(tmp_fd)
                    w.(tmp_fd{nfd})=real(squeeze(model.W.(tmp_fd{nfd})(fb_i, :, :))');
                end
            end

            traindata_P =  squeeze(model.trains_spatial.Pt{class_i}( fb_i, :, :));
%             traindata_T =  squeeze(model.trains_spatial.Tt(class_i, fb_i, :, 1:size(testdata_tmp, 2)));

            r.prca(fb_i,class_i) =corr2(testdata_tmp'*w.prca, ...
                templatemaker(testdata_tmp, traindata_P'*w.prca));

%             r.trca(fb_i,class_i) =corr2(testdata_tmp'*w.trca, ...
%                 traindata_T'*w.trca);
        end % class_i
    end % fb_i

    [~, results.prca(targ_i)] = max(fb_coefs*r.prca);
%     [~, results.trca(targ_i)] = max(fb_coefs*r.trca);
end % targ_i

%%
function y = filterbank(eeg, fs, idx_fb)
% Filter bank design for decomposing EEG data into sub-band components [1].
%
% function y = filterbank(eeg, fs, idx_fb)
%
% Input:
%   eeg             : Input eeg data
%                     (# of channels, Data length [sample], # of trials)
%   fs              : Sampling rate
%   idx_fb          : Index of filters in filter bank analysis
%
% Output:
%   y               : Sub-band components decomposed by a filter bank.
%
% Reference:
%   [1] X. Chen, Y. Wang, S. Gao, T. -P. Jung and X. Gao,
%       "Filter bank canonical correlation analysis for implementing a 
%       high-speed SSVEP-based brain-computer interface",
%       J. Neural Eng., vol.12, 046008, 2015.
%
% Masaki Nakanishi, 22-Dec-2017
% Swartz Center for Computational Neuroscience, Institute for Neural
% Computation, University of California San Diego
% E-mail: masaki@sccn.ucsd.edu

if nargin < 2
    error('stats:test_fbcca:LackOfInput', 'Not enough input arguments.'); 
end

if nargin < 3 || isempty(idx_fb)
    warning('stats:filterbank:MissingInput',...
        'Missing filter index. Default value (idx_fb = 1) will be used.'); 
    idx_fb = 1;
elseif idx_fb < 1 || 10 < idx_fb
    error('stats:filterbank:InvalidInput',...
        'The number of sub-bands must be 0 < idx_fb <= 10.'); 
end

[num_chans, ~, num_trials] = size(eeg);
fs=fs/2;

passband = [6, 14, 22, 30, 38, 46, 54, 62, 70, 78];
stopband = [4, 10, 16, 24, 32, 40, 48, 56, 64, 72];
pass_end = 80; stop_end = 90;

Wp = [passband(idx_fb)/fs, stop_end/fs];
Ws = [stopband(idx_fb)/fs, pass_end/fs];

[B, A] = cheby1(6, 0.5, Wp);

filterorder = 8;
filtercutoff = [95 105]/(fs);
[f_b, f_a] = butter(filterorder,filtercutoff,'stop');
filtercutoff = [49 51]/(fs);
[f_b1, f_a1] = butter(5,filtercutoff,'stop');

y = zeros(size(eeg));
if num_trials == 1
    for ch_i = 1:num_chans
        tmpy=filtfilt(f_b,f_a, eeg(ch_i, :));
        tmpy=filtfilt(f_b1,f_a1, tmpy);
        y(ch_i, :) = filtfilt(B, A, tmpy);
        
    end % ch_i
else
    for trial_i = 1:num_trials
        for ch_i = 1:1:num_chans
            tmpy=filtfilt(f_b,f_a, eeg(ch_i, :, trial_i));
            tmpy=filtfilt(f_b1,f_a1, tmpy);
            y(ch_i, :, trial_i) = filtfilt(B, A, tmpy);
        end % trial_i
    end % ch_i
end % if num_trials == 1

%%
function tmplate=templatemaker(testdata_tmp, tmplate_filtered)
if size(testdata_tmp, 2)<size(tmplate_filtered, 1)
    tmplate=tmplate_filtered(1:size(testdata_tmp, 2), :);
elseif size(testdata_tmp, 2)>size(tmplate_filtered, 1)
    n=floor(size(testdata_tmp, 2)/size(tmplate_filtered, 1));
    tmplate=repmat(tmplate_filtered, n, 1);
    tmplate=[tmplate; tmplate_filtered( 1:(size(testdata_tmp, 2)-size(tmplate, 1)), :)];
else
    tmplate=tmplate_filtered;
end

