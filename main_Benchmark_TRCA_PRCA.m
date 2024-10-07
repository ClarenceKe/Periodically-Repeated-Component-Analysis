% This code shows how to use the Periodically Repeated Component Analysis for SSVEP recognition.

% Please refer the following papers for more details:
% Ke, Yufeng; Liu, Shuang; Ming, Dong (2024). Enhancing SSVEP Identification with Less Individual Calibration Data 
% Using Periodically Repeated Component Analysis.  IEEE Transactions on Biomedical Engineering. 71(4): 1319 - 1331. https://doi.org/10.1109/TBME.2023.3333435

% This code is prepared by Yufeng Ke (clarenceke AT tju.edu.cn)
% Date: 
% 22 July 2023 (v1.0)

% The TRCA Algorithm was adapted from https://github.com/edwin465/Comparison-between-two-implementation-approaches-of-TRCA-in-Matlab

% if you use this code for a publication, please cite the following paper:
% Ke, Yufeng; Liu, Shuang; Ming, Dong (2024). Enhancing SSVEP Identification with Less Individual Calibration Data 
% Using Periodically Repeated Component Analysis.  IEEE Transactions on Biomedical Engineering. 71(4): 1319 - 1331. https://doi.org/10.1109/TBME.2023.3333435
%%

clear

% Please download the SSVEP benchmark dataset for this code
% Wang, Y., et al. (2016). A benchmark dataset for SSVEP-based brain-computer interfaces. IEEE Transactions on Neural Systems and Rehabilitation Engineering, 25(10), 1746-1752.

filepath='Benchmark Dataset\';
files=findfilenames(filepath, 'S', 'mat');
load('Benchmark Dataset\Freq_Phase.mat')

num_targs=40;
labels = [1:1:num_targs];    
is_ensemble=1;
base_time=0.5;

fs=250;
num_fbs=5;
delay_s = 0.14;  

global passband stopband pass_end stop_end;
passband = [6, 14, 22, 30, 38, 46, 54, 62, 70, 78];
stopband = [4, 10, 16, 24, 32, 40, 48, 56, 64, 72];
pass_end = 80; stop_end = 90;

num_of_harmonics=5;

%%      single trial training,  effects of testing data length
for nsub=1:35
%% load data
        load([filepath files{nsub}]);
        eeg_tmp=permute(data(:, :, :, :), [3 1 2 4]);
        clear data
        for ntar=1:size(eeg_tmp, 1)
            for ntri=1:size(eeg_tmp, 4)
                eeg(ntar, :, :, ntri)=resample(squeeze(eeg_tmp(ntar, :, :, ntri))', 1000, 250)'; % upsampling to 1000 Hz
            end
        end
        clear eeg_tmp

        fs=1000;
        N=size(eeg,3); 
        time    = linspace(0, N, N) / fs - base_time;
        fprintf('\n Sub  %d', nsub);  
        %% Estimate classification performance
        global a; global b; a=1.25; b=0.25;
        sti_f=freqs;
        sti_phase=phases;
        channels={[61:63], [55:57 61:63], [54:58 61:63], [48 53:59, 61:63]};
        chan_num=4;
        ChannelUse=channels{chan_num};
        dataLength=[2:15];  %  testing data length (×0.1 s)

        datL_train=10;% fixed training data length = 1 s
        channels={[61:63], [55:57 61:63], [54:58 61:63], [48 53:59, 61:63]};
        train_num=[1:5]; %  number of training trials

        for train_n=1:length(train_num)
            fprintf('\n Sub  %d:  %d trials \n ', nsub, train_num(train_n));
            cv_label=nchoosek(1:size(eeg, 4),train_num(train_n));  
            num_cv = size(cv_label,1);
            ChannelUse=channels{4};

            [acc, itrr]=FreqRecgMethods_TRCA_PRCA(eeg(:, ChannelUse, :, :), labels, num_targs, ...
                dataLength, datL_train, sti_f, sti_phase, num_cv, cv_label, ...
                delay_s, time, fs, num_fbs, is_ensemble, num_of_harmonics);
         
            tmp_fd=fields(acc);
            for nfd=1:length(tmp_fd)
                tmp_fd2=fields(acc.(tmp_fd{nfd}));
                for nfd2=1:length(tmp_fd2)
                    ssvep_Recg.acc.(tmp_fd{nfd}).(tmp_fd2{nfd2})(:, train_n, nsub)=acc.(tmp_fd{nfd}).(tmp_fd2{nfd2});
                    ssvep_Recg.itr.(tmp_fd{nfd}).(tmp_fd2{nfd2})(:, train_n, nsub)=itrr.(tmp_fd{nfd}).(tmp_fd2{nfd2});
                    ssvep_Recg.data_length=0.1*dataLength;
                    ssvep_Recg.train_num_trials=train_num;
                end
            end
        end
    save('Benchmark_TRCA_PRCA', 'ssvep_Recg')
end

