%function mnist_svd
clear; close all; clc;
tic;
[ImageTrain,LabelTrain] = digitTrain4DArrayData();
[ImageValidate,LabelValidate] = digitTest4DArrayData();

ImageValidate = squeeze(ImageValidate);
ImageTrain = squeeze(ImageTrain);

nTrain = 5000;
nValidate = 5000;
nrows = 28;
ncols = 28;
totalPixels = nrows*ncols;

% vectorize the images
ImageTrainTmp = NaN(totalPixels,nTrain);
ImageValidateTmp = ImageTrainTmp;
toc;

for i = 1:nTrain
    im_ = squeeze(ImageTrain(:,:,i));
    im_ = im_(:);
    ImageTrainTmp(:,i) = im_;
end

for i = 1:nValidate
    im_ = squeeze(ImageValidate(:,:,i));
    im_ = im_(:);
    ImageValidateTmp(:,i) = im_;
end
ImageTrain = normalizeWaveforms(ImageTrainTmp);
ImageValidate = normalizeWaveforms(ImageValidateTmp);
clear ImageTrainTmp ImageValidateTmp

%
maxDim = 10;
nCategories = 10; %10 digits therefore 10 categories
wjk1 = NaN(maxDim*nCategories,totalPixels);
n = 0;
for i = 0:nCategories-1
    catI = LabelTrain == categorical(i);
    n = i;
    theseImages = ImageTrain(:,catI);
    [U,~,~] = svd(wiener2(theseImages,[3,3]));
    %[U,~,~] = svd(theseImages);
    wjk1(maxDim*n+1:maxDim*n+maxDim,:) = U(:,1:maxDim)';
end
toc;

%
first_row = [ones(1,maxDim) zeros(1,maxDim*(nCategories-1))];
wjk2 = zeros(nCategories,maxDim*nCategories);
wjk2(1,:) = first_row;
for i = 2:nCategories
    wjk2(i,:) = circshift(wjk2(i-1,:),maxDim);
end
toc;

%% training section -- incomplete, finish
z1 = wjk1*ImageTrain;
hidden = z1.^2;
out = wjk2*hidden;
y = zeros(size(out));
LabelTrain = double(string(LabelTrain));
for i = 1:nTrain
    index = LabelTrain(i);
    y(index+1,i) = 1;
end
C = sum((out-y).^2);
toc;

[~,maxI] = max(out,[],1);
labelGuess = maxI' - 1;
labelError = labelGuess-LabelTrain;
sum(labelError == 0)/nTrain
toc;

%% validatation section
z1 = wjk1*ImageValidate;
hidden = z1.^2;
out = wjk2*hidden;
y = zeros(size(out));
LabelValidate = double(string(LabelValidate));
for i = 1:nValidate
    index = LabelValidate(i);
    y(index+1,i) = 1;
end
C = sum((out-y).^2);
toc;

[~,maxI] = max(out,[],1);
labelGuess = maxI' - 1;
labelError = labelGuess-LabelValidate;
sum(labelError == 0)/nValidate
toc;


%% from here on out, the code works to apply CNN per a matlab online example
% clear; close all;
% dataFolder = fullfile(toolboxdir('nnet'),'nndemos','nndatasets','DigitDataset');
% imds = imageDatastore(dataFolder, ...
%     'IncludeSubfolders',true,'LabelSource','foldernames');
% figure
% tiledlayout("flow");
% perm = randperm(10000,20);
% for i = 1:20
%     nexttile
%     imshow(imds.Files{perm(i)});
% end
%
% %%
% classNames = categories(imds.Labels);
% labelCount = countEachLabel(imds);
%
% %%
% img = readimage(imds,1);
% size(img)
%
% %%
% numTrainFiles = 750;
% [imdsTrain,imdsValidation] = splitEachLabel(imds,numTrainFiles,"randomize");
%
% %%
% layers = [
%     imageInputLayer([28 28 1])
%
%     convolution2dLayer(3,8,Padding="same")
%     batchNormalizationLayer
%     reluLayer
%
%     maxPooling2dLayer(2,Stride=2)
%
%     convolution2dLayer(3,16,Padding="same")
%     batchNormalizationLayer
%     reluLayer
%
%     maxPooling2dLayer(2,Stride=2)
%
%     convolution2dLayer(3,32,Padding="same")
%     batchNormalizationLayer
%     reluLayer
%
%     fullyConnectedLayer(10)
%     softmaxLayer];
%
% %%
% options = trainingOptions("sgdm", ...
%     InitialLearnRate=0.01, ...
%     MaxEpochs=4, ...
%     Shuffle="every-epoch", ...
%     ValidationData=imdsValidation, ...
%     ValidationFrequency=30, ...
%     Plots="training-progress", ...
%     Metrics="accuracy", ...
%     Verbose=false);
%
% %%
% net = trainnet(imdsTrain,layers,"crossentropy",options);
%
% %%
% scores = minibatchpredict(net,imdsValidation);
% YValidation = scores2label(scores,classNames);
% TValidation = imdsValidation.Labels;
%
% %%
% accuracy = testnet(net,imdsValidation,"accuracy")