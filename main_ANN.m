clear all
setGlobalx(0);

mbs = 5;
sz = 50;

layer = Output1RegressionLayer('mse');

ins = linspace(1e-4, 1e-2,sz);
ins = normalize(ins,'range',[0.1,1]);
targs = 10.^-4*ones(1,sz);

for i = 1:sz
    targs(i) = targs(i) + i*10^-6;
end

A = (randperm(length(1:sz)));

for i = 1:sz
    inputs(i) = ins(A(i));
    targets(i) = targs(A(i));
end


layers = [
    sequenceInputLayer(1)
    lstmLayer(200,'OutputMode',"sequence")  
    fullyConnectedLayer(1)
    Output1RegressionLayer('mse')];

N = 30; %number of sequences

cellArrTrain = cell(N,1);
cellArrTarget = cell(N,1);
for i = 1:N
    seqq = inputs(i);
    sepp = targets(i);
    seqq = num2cell(seqq,1);
    sepp = num2cell(sepp,1);
    cellArrTrain(i) = seqq;
    cellArrTarget(i) = sepp;
end

XValidation = cell((sz-N),1);
YValidation = cell((sz-N),1);
for i = N+1:sz
    seq = inputs(i);
    sep = targets(i);
    seq = num2cell(seq,1);
    sep = num2cell(sep,1);
    XValidation(i-N) = seq;
    YValidation(i-N) = sep;
end

options = trainingOptions('adam', ...
    'MaxEpochs',2, ...
    'GradientThreshold',3, ...
    'InitialLearnRate',0.003, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',1, ... 
    'LearnRateDropFactor',0.4, ... 
    'ValidationData',{XValidation,YValidation}, ...
    'ValidationFrequency',1, ...
    'LearnRateDropFactor',1,...
    'Verbose',1, ...
    'MiniBatchSize',mbs, ...
    'Plots','training-progress');


net = trainNetwork(cellArrTrain,cellArrTarget,layers,options);

YPred = predict(net,cellArrTrain);

