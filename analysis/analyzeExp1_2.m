%analyze exp1_2 data

data_dC= data; %de-cluttered
data_dC(1:100) = 0;
% data_dC = data_dC - median(data_dC);
data_err = sqrt(data_dC);

gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewedgaussian = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x);
fitSkewGauss = @(x,a,b,c,d) a*skewedgaussian((x-b)/c,d);
weights = 1./data_err;
weights(data_err == 0) = 0;
% std = 4;
% alpha = -4;
sumOfFSG = @(a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8,...
    std,alpha,x)...% d1,d2,d3,d4,d5,d6,d7,d8,x)...
    ...%     c1,c2,c3,c4,c5,c6,c7,c8,d1,d2,d3,d4,d5,d6,d7,d8,x)...
    fitSkewGauss(x,a1,b1,std,alpha)+...
    fitSkewGauss(x,a2,b2,std,alpha)+...
    fitSkewGauss(x,a3,b3,std,alpha)+...
    fitSkewGauss(x,a4,b4,std,alpha)+...
    fitSkewGauss(x,a5,b5,std,alpha)+...
    fitSkewGauss(x,a6,b6,std,alpha)+...
    fitSkewGauss(x,a7,b7,std,alpha)+...
    fitSkewGauss(x,a8,b8,std,alpha);

[fitResult,gof] = fit([0:2047]', data_dC', sumOfFSG,...
    'StartPoint',[1e3*ones(1,8),[800,810,850,900,910,940,1010,1310],7,-7]...
    ,'Lower',[10*ones(1,8),zeros(1,8),0.5,-10]...
    ,'Upper',[1e4*ones(1,8),2e3*ones(1,8),10,-0.1]...
    ,'MaxFunEvals',2e3,'MaxIter',1e4,'weights',weights);

figure
hold on
h = plot(fitResult,0:2047,data_dC);
set(h(2),'linewidth',2)
plot(0:2047,data_dC-data_err,':b')
plot(0:2047,data_dC+data_err,':b')

errors = diff(confint (fitResult));

energy_dataByChannel = [fitResult.b1,fitResult.b2,fitResult.b3,...
    fitResult.b4,fitResult.b5,fitResult.b6,fitResult.b7,fitResult.b8];
[energy_fromLit,energyIX] = sort([5423,5340,5685,6288,6778,6051,6090,8784]); % keV

energy_dataError = errors(9:end-2); 

eFrequency_dataByChannel = [fitResult.a1,fitResult.a2,fitResult.a3,...
    fitResult.a4,fitResult.a5,fitResult.a6,fitResult.a7,fitResult.a8];
eFrequency_dataError = errors(1:8)/sum(eFrequency_dataByChannel(1:2));
eFrequency_dataByChannel = eFrequency_dataByChannel/sum(eFrequency_dataByChannel(1:2));
eFrequency_fromLit = 1e-2*[71,28,95,100,100,25,10,64];
eFrequency_fromLit = eFrequency_fromLit(energyIX);

weightEBin = energy_dataByChannel./energy_dataError;
[fitEnergyChannel,gofEnCh] = fit(energy_dataByChannel',energy_fromLit','poly1','weights',weightEBin);
% [fitEnergyChannelLin,gofEnChLin] = fit(energy_dataByChannel',energy_fromLit','poly1','Upper',[inf,0],'Lower',[-inf,0]);
figure
subplot(2,1,1)
hold on
plot(fitEnergyChannel,energy_dataByChannel,energy_fromLit)
subplot(2,1,2)
title('residuals')
plot(energy_dataByChannel,feval(fitEnergyChannel,energy_dataByChannel) - energy_fromLit','.')

figure
hold on
plot(eFrequency_fromLit,'o')
errorbar(eFrequency_dataByChannel,eFrequency_dataError,'.')
plot(eFrequency_dataByChannel,'x')