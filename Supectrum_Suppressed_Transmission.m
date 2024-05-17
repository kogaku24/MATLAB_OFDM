%% OFDM
BW = 20e6; %帯域幅[Hz]
SCS = 15e3; %サブキャリア間隔[Hz]
FFTsize = 2048; %FFT size
cp = 4.76; %CP[μs]  
M=4;

rate = 602/1024;
I = 25;
DecisionType='hard'; 
Algorithm = 'Normalized min-sum';
OutputType = 'llr';

Nsub = fix(BW/SCS); %サブキャリア数
ceil_GBlen = ceil((FFTsize-Nsub)/2);%ガードバンド
floor_GBlen = floor((FFTsize-Nsub)/2);%ガードバンド
Nsym = 1024;
Ndata = Nsub*(Nsym*2)*rate;%データ数
symlen = (1/SCS)*1e6;%シンボル長[μs]
CPlen = fix(Nsub/(symlen/cp));%CP
samplerate = (FFTsize*SCS);%サンプルレート[Hz]


% 送信側
%% PN生成器
pnSequence = comm.PNSequence('Polynomial',[21 19 0],'SamplesPerFrame',Ndata,'InitialConditions',1);
pnSignal= pnSequence();

%% LDPCEncode
A = length(pnSignal); %transport block length
rv = 0; 
modulation = 'QPSK'; %変調
nlayers = 1; 
%UL-SCH coding parameters
%UL-SCH符号化パラメータを決定
cbsInfo = nrULSCHInfo(A,rate);
%LDPC
in = pnSignal;
tbIn = nrCRCEncode(in,cbsInfo.CRC); %transport block CRC attachment　CRC:Code Block Concatenation
cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN); %Code block segmentation and CRC attachment
enc = nrLDPCEncode(cbsIn,cbsInfo.BGN); %LDPC encoding
outlen = ceil(A/rate);
chIn = nrRateMatchLDPC(enc,outlen,rv,modulation,nlayers);

%% S-P変換
s_pSignal = reshape(chIn,FFTsize-(ceil_GBlen + floor_GBlen),[]); %reshape

%% サブキャリア変調
subcarriermod = zeros(Nsub,Nsym);
for i = 1:Nsub
    subcarriermod(i,:) = qammod(s_pSignal(i,:).',M,'InputType','bit','UnitAveragePower',true).';
end

%% スペクトラム圧縮
s_sub=1;
f_sub=134;

for j = s_sub:f_sub
    subcarriermod(j,:) = zeros(1,Nsym);
end

%% OFDM変調
hMod = comm.OFDMModulator("FFTLength",FFTsize,"NumGuardBandCarriers",[ceil_GBlen;floor_GBlen],"CyclicPrefixLength",CPlen,"NumSymbols",Nsym);
hmodSignal = hMod(subcarriermod);
paprSignal = hmodSignal;

snr=0;
[AWGNSignal,noiseVariance] = awgn(hmodSignal,snr,-10*log10(FFTsize));

% 受信側
%% OFDMdemod
hdemod = comm.OFDMDemodulator("FFTLength",FFTsize,"NumGuardBandCarriers",[ceil_GBlen;floor_GBlen],"CyclicPrefixLength",CPlen,"NumSymbols",Nsym);
ofdmdemod = hdemod(AWGNSignal);

%% サブキャリア復調
subcarrierdemod = zeros(Nsub,Nsym*2);
for ii = 1:Nsub
    subcarrierdemod(ii,:) = qamdemod(ofdmdemod(ii,:).',M,'OutputType',OutputType,'UnitAveragePower',true).';
end

%% P-S変換
p_sSignal = reshape(subcarrierdemod,[],1);

%% LDPCDecode
raterec = nrRateRecoverLDPC(p_sSignal,A,rate,rv,modulation,nlayers);
decBits = nrLDPCDecode(raterec,cbsInfo.BGN,I,'DecisionType',DecisionType,'Algorithm',Algorithm); %LDPC decoding
blk = nrCodeBlockDesegmentLDPC(decBits,cbsInfo.BGN,A+cbsInfo.L);
out = nrCRCDecode(blk,cbsInfo.CRC); %transport block CRC decoding


% 特性
%% PAPR
valPAPR = zeros(Nsym,3);
paprSignal2 = reshape(paprSignal,FFTsize+CPlen,[]);%FFTsize+CPlen

for index = 1:Nsym
    valPAPR(index,1) = max(abs(paprSignal2(:,index)).^2);
    valPAPR(index,2) = mean(abs(paprSignal2(:,index)).^2);
    valPAPR(index,3) = 10*log10(valPAPR(index,1)./valPAPR(index,2));
end

data = valPAPR(1:end,3);
data = sort(data);

count = 0;
datanow = 0;
mn = zeros(length(data),1);


for index1 = 1:length(data)
    if datanow == data(index1,1)
        count2 = 1;
    end
    if datanow < data(index1,1)
        count = count + 1;
        datanow = data(index1,1);
        mn(count) = data(index1,1);
    end
end

% Graph
CDFdata = zeros(count,5);

for index2 = 1:count
    CDFdata(index2,1) = mn(index2);
end
for index3 =  1:count
    CDFdata(index3,2) = length(find(data==CDFdata(index3,1)));
end
sum = 0;
for index4 = 1:count
    sum = sum + CDFdata(index4,2);
    CDFdata(index4,3) = sum;
end
for index5 = 1:count
    CDFdata(index5,4) = CDFdata(index5,3)/length(data);
end
for index6 = 1:count
    CDFdata(index6,5) = 1 - CDFdata(index6,4);
end
%
x = CDFdata(1:end,1);
y = CDFdata(1:end,5);
G=semilogy(x,y);
xlabel('PAPR(dB)')
ylabel('CCDF')
grid on
xlim([0 20]);
ylim([1e-5 1])


%% BER
err = comm.ErrorRate('ReceiveDelay',0);
ber = err(pnSignal, double(out));