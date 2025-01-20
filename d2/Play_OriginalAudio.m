% This script is used for TELECOM201b Lab
% It only plays the audio file out.wav

% Rev: March 2023, Germain

% Load file
[InputAudio,fsAudio]    = audioread('lab_data/out.wav');

% Prepare signal
truncate   = 2^19;  % 2^19 = 524288
achannel   = 1;     % 1 = left channel, 2 = right channel
InputAudio = InputAudio(1:truncate,achannel)/max(InputAudio(:,achannel));

% Play audio
soundsc(InputAudio,fsAudio);