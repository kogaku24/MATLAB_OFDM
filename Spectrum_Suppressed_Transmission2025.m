% SST
%profile on;
clear
tStart = tic;
%% Parameters
%% input
%rng(123)
%% EXCEL BER
filename = 'file.xlsx';
%% shift
min_shift = 0;
max_shift = 100;
%% SNR
min_snrdB = -5;
max_snrdB = 5;
snr = 1;
%% S Ratio 圧縮率
min_SR = 0; % [%]
max_SR = 10; % [%]
SRS = 10; % [%]

%% MCS
MCS_sw = 1;  %1:MCS_on, 0:custom
load("MCS256QAM.mat")
%0~4:QPSK, 5~10:16QAM, 11~19:64QAM, 20~27:256QAM
min_MCSindex = 30;
max_MCSindex = 30;
space = 1;

% switch
BER_short = 1; %0:off, :onフロアを引いたら終了
inter_sw = 1; %0:off, 1:on



%% Base

%60kHz,100MHz
BW = 100e6; %bandwidth[Hz]
SCS = 60e3; %sub-carrier spacing[Hz]
Nrb = 132;  % Number of resource blocks

%{
%30kHz,100MHz
BW = 100e6; %bandwidth[Hz]
SCS = 30e3; %sub-carrier spacing[Hz]
Nrb = 273;
%}
%{
%15kHz,20MHz
BW = 20e6;
SCS = 15e3;
Nrb = 79;
%}
rbSize = 12;% Number of subcarriers per resource block
Nslot=100; %Number of slot
%Nsym = 14*Nslot;

%Nsym = 1500;
Nsym = 14;
Nsymbol = 14;
cp = 1.17; %cp[μs]


Nexcel = 5;

%% p-generate data

%% p-enc
%% p-ldpc
%% para ldpc enc
rv = 0;
nlayers = 1;
%% para ldpc dec
%I = 8;
I = 10;
DecisionType = 'hard';
%DecisionType='soft';
Algorithm = 'Normalized min-sum';

%% para SC
%% para scMod

%% para scDeMod
OutputType ='approxllr';
%OutputType ='llr';
%OutputType ='bit';
%% para OFDM
%% para ofdmMod
FFTsize = 2048;   % Number of FFT points
%FFTsize = 4096;   % Number of FFT points

%% calc p-Generate data
Nsub = Nrb*rbSize; % number of data subcarriers in sub-band
%BPSC = log2(M); %bitsPerSubCarrier 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM
%Npn = floor(Nsub*Nsym*rate);

%% calc para-ldpc
%tblen = Npn; %transport block length  

%% calc para-ofdm mod
ceil_GBlen = ceil((FFTsize-Nsub)/2);%ガードバンド
floor_GBlen = floor((FFTsize-Nsub)/2);%ガードバンド
symlen = (1/SCS)*1e6;%シンボル長[μs]
CPlen = fix(Nsub/(symlen/cp)); %CP % Cyclic prefix length in samples ??????????

%% calc-para SA
samplerate = (FFTsize*SCS);%サンプルレート[Hz]

%% para BER
BER_txrx = zeros(50,1);
snr_txrx = zeros(50,1);
BER_pre = zeros(50,1);
snr_pre = zeros(50,1);



    for MCSindex = min_MCSindex:space:max_MCSindex %MCS

        M = MCS(MCSindex+1,2);
        rate = MCS(MCSindex+1,3); 
 
        switch M
            case 4
                mod = 'QPSK';
            case 16
                mod = '16QAM';
            case  64
                mod = '64QAM';
            case  256
                mod = '256QAM';
        end

        %disp("MSC:"+MCS(MCSindex+1,1)+", mod:"+mod)
        %% calc p-Generate data
        BPSC = log2(M); %bitsPerSubCarrier 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM

        %% calc p-Generate data
        %Nsym = Nsymbol/BPSC;
        Npn = floor(Nsub*Nsym*rate*BPSC);
        %Npn = 20000;
        pnSequence = comm.PNSequence('Polynomial',[21 19 0], ...
                        'SamplesPerFrame',Npn,'InitialConditions',1,'Mask',0);

        for SR = min_SR:SRS:max_SR %圧縮率
            nsc = ceil(Nsub*(SR/100));
            %disp("圧縮率="+SR+"%, Null="+nsc)
           
            
            Nber = 0;
            for snrdB = min_snrdB:snr:max_snrdB % snrdB
                Nber = Nber+1;
                tSnr = tic;
                for shift = min_shift:max_shift %pn shift
                    %% TX
                    
                    
                    
                    %% Generate data
                    %pnSequence = comm.PNSequence('Polynomial',[21 19 0], ...
                    %    'SamplesPerFrame',Npn,'InitialConditions',1,'Mask',(Npn*(shift)-1));
                    release(pnSequence);
                    pnSequence.Mask = Npn*(shift);
                    pnSig= pnSequence();
                    
                    %% Enc
                    %% LDPC enc
                    cbsInfo = nrULSCHInfo(Npn,rate); %UL-SCH coding parameters %UL-SCH符号化パラメータを決定
                    in = pnSig;
                    tbIn = nrCRCEncode(in,cbsInfo.CRC); %transport block CRC attachment　CRC:Code Block Concatenation
                    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN); %Code block segmentation and CRC attachment
                    enc = nrLDPCEncode(cbsIn,cbsInfo.BGN); %LDPC encoding
                    %outlen = ceil(tblen/rate);
                    outlen = ceil(Nsub*Nsym*BPSC);
                    chIn = nrRateMatchLDPC(enc,outlen,rv,mod,nlayers);
                    encSig=chIn;
            
                    if inter_sw == 1
                        % パターン1の処理
                        b_inter = encSig;
                        a_inter = reshape(reshape(b_inter, Nsub, []).', [], 1);
                        encSig = a_inter;
                        disp("interleave:" + inter_sw);
                    end
            
                    %% S/P
                    tx_spSig = reshape(encSig,Nsub,[]);
                    
                    %% SC mod
                    %scMod = zeros(Nsub, (outlen/Nsub)/BPSC);
                    
                    scMod = qammod(tx_spSig.',M, 'InputType', 'bit', ...
                                'UnitAveragePower', true).';
                    scSig = scMod;

                    %%電力補正用scSig
                    scSig0 = scMod; 
                    
                    %% SST low (after SC mod)
                    s_sub=1;
                    f_sub=nsc;

                    scSig(s_sub:f_sub,:) = zeros(f_sub-s_sub+1,(outlen/Nsub)/BPSC);

                    %% OFDM mod
                    ofdmMod = comm.OFDMModulator("FFTLength",FFTsize, ...
                        "NumGuardBandCarriers",[ceil_GBlen;floor_GBlen], ...
                        "CyclicPrefixLength",CPlen, ...
                        "NumSymbols",(outlen/Nsub)/BPSC); %?????????????
                    
                    ofdmSig = ofdmMod(scSig);

                    %電力補正用ofdmSig
                    ofdmMod0 = comm.OFDMModulator("FFTLength",FFTsize, ...
                        "NumGuardBandCarriers",[ceil_GBlen;floor_GBlen], ...
                        "CyclicPrefixLength",CPlen, ...
                        "NumSymbols",(outlen/Nsub)/BPSC); %?????????????
                    ofdmSig0 = ofdmMod0(scSig0);
                    
                    
                    
                    
                    
                    %% Channel
                    
                    rng(123) % awgn rng
                    %% AWGN
                    %sigPower = sum(abs(ofdmSig(:)).^2)/numel(ofdmSig); %電力補正なし
                    sigPower = sum(abs(ofdmSig0(:)).^2)/numel(ofdmSig0); %電力補正あり
                    %sigPower =3.776224407607219e-04;
                    %awgnSig = sAWGN(ofdmSig,sigPower,snrdB+10*log10(Nsub/FFTsize));
                    awgnSig = sAWGN(ofdmSig,sigPower,snrdB);
                    %awgnSig = awgn(ofdmSig,snrdB,"measured");
                    
                    
                    
                    
                    
                    %% RX
                    
                    
                    %% OFDM demod
                    ofdmDemod = comm.OFDMDemodulator("FFTLength",FFTsize, ...
                        "NumGuardBandCarriers",[ceil_GBlen;floor_GBlen], ...
                        "CyclicPrefixLength",CPlen, ...
                        "NumSymbols",(outlen/Nsub)/BPSC); %?????????????
                    
                    ofdmDeSig = ofdmDemod(awgnSig);
                    
                    %% SC demod
                    %scDeMod = zeros(Nsub,outlen/Nsub);
                    
                    scDeMod = qamdemod(ofdmDeSig.',M,'OutputType',OutputType,'UnitAveragePower',true).';
   
                    scDeSig = scDeMod;
                    

                    %% P/S
                    rx_psSig = reshape(scDeSig,[],1);

                    if inter_sw == 1
                        % パターン1の処理
                        b_deinter = rx_psSig;
                        a_deinter = reshape(reshape(b_deinter, outlen/Nsub, []).', [], 1);
                        rx_psSig = a_deinter;
                    end
                    

                    %% LDPC dec
                    raterec = nrRateRecoverLDPC(rx_psSig,Npn,rate,rv,mod,nlayers);
                    decBits = nrLDPCDecode(raterec,cbsInfo.BGN,I,'DecisionType',DecisionType,'Algorithm',Algorithm); %LDPC decoding
                    blk = nrCodeBlockDesegmentLDPC(decBits,cbsInfo.BGN,Npn+cbsInfo.L);
                    decSig = nrCRCDecode(blk,cbsInfo.CRC); %transport block CRC decoding
                    %decSig = decSig < 0;  %% soft
                    
                    out = decSig;
                    
                    if shift == min_shift
                        pnSigData = pnSig;
                        outData = out;
                        encSigData = encSig;
                        rx_psSigData = rx_psSig;
                        
                    else
                        pnSigData = cat(1,pnSigData,pnSig);
                        outData = cat(1,outData,out);
                        encSigData = cat(1,encSigData,encSig);
                        rx_psSigData = cat(1,rx_psSigData,rx_psSig);
                    end

                    if shift == max_shift
                        pnSig = pnSigData;
                        out = outData;
                        encSig = encSigData;
                        rx_psSig = rx_psSigData;
                        %% Analyzer
                        %% BER
                        err = comm.ErrorRate('ReceiveDelay',0,'ResetInputPort',true);
                        ber = err(pnSig, double(out),1);
                        
                        %% preBER
                        preBits = rx_psSig <= 0;
                        err_pre = comm.ErrorRate('ReceiveDelay',0,'ResetInputPort',true);
                        ber_pre = err_pre(encSig, double(preBits),1);
                    end
                    disp("---------------------------------")
                    disp("shift:"+shift)
                    disp("データ数:"+sprintf('%e',Npn)+", Power:"+sigPower)
                    disp(" ")
                    disp("MSC:"+MCS(MCSindex+1,1)+", mod:"+mod+", rate:"+string(rate))
                    disp("圧縮率:"+SR+"%, Null:"+nsc)
                    disp("snrdB:"+snrdB)
                    disp("Nsymbol:"+Nsymbol)
                    disp("---------------------------------")
    
                    disp(toc(tSnr)+"[sec], Total:"+ceil(toc(tStart))+"[sec]("+toc(tStart)/60+"[min])")
                
                end % shift
            
                BER_txrx(Nber,1) = ber(1,:);
                snr_txrx(Nber,1) =snrdB;
                BER_pre(Nber,1) = ber_pre(1,:);
                snr_pre(Nber,1) =snrdB;
                
                disp("データ数:"+sprintf('%e',Npn)+", Power:"+sigPower)
                disp(" ")
                disp("MSC:"+MCS(MCSindex+1,1)+", mod:"+mod+", rate:"+string(rate))
                disp("圧縮率:"+SR+"%, Null:"+nsc)
                disp("snrdB:"+snrdB)
                disp("BER:"+ber(1,:)+", preBER:"+ber_pre(1,:))
                disp("Nsymbol:"+Nsymbol)
                disp("---------------------------------")

                disp(toc(tSnr)+"[sec], Total:"+ceil(toc(tStart))+"[sec]("+toc(tStart)/60+"[min])")
                %% 処理時間短縮
               %if BER_short == 0
               if BER_short == 1 && Nber >= 5
                    % 必要な部分だけを一度に丸める
                    BER_samples = round(BER_txrx(Nber-2:Nber, 1), 5);
                
                    % 比較結果を一度に計算
                    BER_case01 = BER_samples(3) == BER_samples(2); % BER_sample0 == BER_sample1
                    BER_case02 = BER_samples(3) == BER_samples(1); % BER_sample0 == BER_sample2
                
                    % 条件分岐を簡略化
                    if BER_case01 && BER_case02
                        break;
                    end
                end

            end % snrdB
        
        %% Constellation
        constDiagram = comm.ConstellationDiagram;
        %constDiagram(scSig(1,:).');
        %constDiagram(ofdmDeSig(1,:).');
        
        %% Spectrum Analyzer
        %% matlab spectrum analyzer
        scope = spectrumAnalyzer("SampleRate",samplerate);
        %scope(ofdmSig);
        %scope(reshape(ifftSig,[],1));
        
        %% custom spectrum analyzer
        %{
        % fft spectrum
        up=4;
        %up = 1;
        point = up*up;
        Pfft= FFTsize*point;
        speFFT = fftshift(fft((ofdmSig),Pfft));
        %speFFT = (fft((ofdmSig),Pfft));
        speLog = 10*log10(abs(speFFT).^2);
        speLog = speLog - max(speLog);
        x = 1:1:Pfft;
        xd = -61.44:(60/(point))*1e-3:61.44-(60/(point))*1e-3;
        y = speLog(x,1);
        plot(xd,y);
        xlim([-61.44 61.44-15e-3])
        ylim([-90,0])
        %}
        
        
        %% EXCEL_BER
        excel26 = excel_column(Nexcel);
        time = datetime('now',Format='MMdd_HHmmss_yyyy');
        t = string(time);
        %name = 'BER';
        %extension = '.xlsx';
        %filename = [name, t, extension];
        %filename ='BER2024.xlsx';

        %sheet = t;
        sheet = string(MCS(MCSindex+1,1));

        writematrix("MCS:"+MCS(MCSindex+1,1),filename,'sheet',sheet,'Range','B2');
        writematrix(['mod:',mod],filename,'sheet',sheet,'Range','B3');
        writematrix("rate:" + string(rate),filename,'sheet',sheet,'Range','B4');
        writematrix("I:" + string(I),filename,'sheet',sheet,'Range','B5');
        writematrix(['DecType:',DecisionType],filename,'sheet',sheet,'Range','B6');
        writematrix(['Algorithm :',Algorithm],filename,'sheet',sheet,'Range','B7');
        writematrix("BW:"+ BW,filename,'sheet',sheet,'Range','B8');
        writematrix("SCS:" + SCS,filename,'sheet',sheet,'Range','B9');
        writematrix("DecInput:" + OutputType,filename,'sheet',sheet,'Range','B10');
        writematrix("圧縮率:" + SR,filename,'sheet',sheet,'Range','B11');
        writematrix("Ndata:" + Npn,filename,'sheet',sheet,'Range','B12');
        writematrix("fsub" + f_sub,filename,'sheet',sheet,'Range','B13');
        writematrix("ssub" + s_sub,filename,'sheet',sheet,'Range','B14');
        writematrix("i_off:0,i_on:1 = " + inter_sw,filename,'sheet',sheet,'Range','B15');
        writematrix("信号パワー:" + sigPower,filename,'sheet',sheet,'Range','B16');
        writematrix("インターリーブのシンボル数:" + Nsymbol,filename,'sheet',sheet,'Range','B16');
        
        writematrix('SNR',filename,'sheet',sheet,'Range','D2');
        %writematrix('BER',filename,'sheet',sheet,'Range','E2');
        %writematrix('preBER',filename,'sheet',sheet,'Range','G2');
        writematrix(snr_txrx,filename,'sheet',sheet,'Range','D3');
        writematrix("圧縮率:" + SR,filename,'sheet',sheet,'Range',excel26+'1');
        writematrix(Nsymbol,filename,'sheet',sheet,'Range',excel26+'2');
        writematrix(BER_txrx,filename,'sheet',sheet,'Range',excel26+'3');
        %writematrix(snr_pre,filename,'sheet',sheet,'Range','D3');
        %writematrix(BER_pre,filename,'Sheet',sheet,'Range','G3');

        BER_txrx = zeros(50,1);
        snr_txrx = zeros(50,1);
        BER_pre = zeros(50,1);
        snr_pre = zeros(50,1);
        
        Nexcel = Nexcel +1;
        
        end %圧縮率
    
    end % MCS
    Nexcel = 5;
%end % inter
tEnd = toc(tStart);
disp("合計処理時間:"+tEnd/60+"分")

%end %symbol

%% function AWGN 
%多値変調に変更
%Ebの換算処理をする

function [r] = sAWGN(gouseisignal,TP,snrdB)
%入力パラメータ
%gouseisignal : 入力信号
%TP:送信信号の電力
%EbN0dB :Eb/N0(dB)
%k:1シンボルあたりのビット数
%sps:1シンボルあたりのサンプル数
%rate:符号化率

%多値変調の時 EbN0ではなくSNで決める必要があるためノイズ調整の必要がある
%% SNRを求める
     snr = 10.^(snrdB/10);
           
%% -----------ノイズ電力------------------------
NP = TP /snr;
%% ----------分散と標準偏差------------------------------
sigma2 = NP /2;                            %分散
sigma =sqrt(sigma2);                     %標準偏差

%% -----------------ノイズベクトルの実数と虚数----------------

x = sigma.* randn(length(gouseisignal),1);     %ノイズベクトルの実数

y = sigma.* randn(length(gouseisignal),1);     %ノイズベクトルの虚数

%% ノイズベクトル
noise = x + 1i*y;                           %ノイズベクトル

%% ----------------受信信号(出力パラメータ)-------------------
r = gouseisignal + noise;

end

%%　EXCEL変換
function e26 = excel_column(column)

e26 = dec2base(round(abs(column)) - 1, 26);
ascii = double(e26);
ascii(1:end-1) = ascii(1:end-1) - 1;
ascii(ascii > 63) = ascii(ascii > 63) + 10;
ascii(ascii > 47 & ascii <= 63) = ascii(ascii > 47 & ascii <= 63) + 17;
% 文字列に変換
e26 = string(char(ascii));

end

%profile viewer