% WorkingP='\\fmri-t9\users\Moran\DCE\DCE_Duration\001_HAVLIN_HAIM_MORDECHAY\Study20140831_120511_day224_T4\DCE_38min\HaMo_20140831\';
%% Load basic stuff
load('Export.mat');
handles=Export;
Idxs=Export.Idxs;
HConvIdxM=CreateConvIdxMFromSampleTs(numel(handles.HSampleTs));
HTriB=HConvIdxM>0;
HConvIdxMTriB=HConvIdxM(HTriB);

% BATfinal VpFinal KtransFinal Kepfinal VeFinal
BATIdx=1;
VpIdx=2;
KtransIdx=3;
KepIdx=4;
CurKeps=handles.CurPKs(:,KepIdx);
CurKeps(isnan(CurKeps))=0;

CurVps=handles.CurPKs(:,VpIdx);
CurVps(isnan(CurVps))=0;

Hdt=diff(handles.HSampleTs);
Hdt=Hdt(1);

HHSampleTs=0:Hdt:handles.HSampleTs(end);
HHAIF=interp1(handles.HSampleTs,handles.HAIF',HHSampleTs);
HHConvIdxM=CreateConvIdxMFromSampleTs(numel(HHSampleTs));
HHTriB=HHConvIdxM>0;
HHConvIdxMTriB=HHConvIdxM(HHTriB);

HHConvd2=DCECostFuncgrT1ForConv(HHAIF',CurKeps,HHSampleTs,HHConvIdxMTriB,HHTriB);
for i=1:size(HHConvd2,1)
    HConvd2(i,:)=interp1(HHSampleTs,HHConvd2(i,:),handles.HSampleTs,[],'extrap');
end
% HConvd2=DCECostFuncgrT1ForConv(handles.HAIF',CurKeps,handles.HSampleTs,HConvIdxMTriB,HTriB);

Hdt=handles.HSampleTs(2)-handles.HSampleTs(1);
dt=diff(handles.SampleTs(1:2));
% ThreeSec=ceil(3/(Hdt*60));
% TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;

CurBATs=handles.CurPKs(:,BATIdx)/-60;
CurBATs(isnan(CurBATs))=1;
for i=1:numel(Idxs)
%     SHConvd2(i,:)=interp1(handles.HSampleTs,HConvd2(i,:)',handles.HSampleTs+handles.TDif(CurBATs(i)),[],'extrap')';
%     SHSAIF(i,:)=interp1(handles.HSampleTs,handles.HAIF,handles.HSampleTs+handles.TDif(CurBATs(i)),[],'extrap');
    SHConvd2(i,:)=interp1(handles.HSampleTs,HConvd2(i,:)',handles.HSampleTs+CurBATs(i),[],'extrap')';
    SHSAIF(i,:)=interp1(handles.HSampleTs,handles.HAIF,handles.HSampleTs+CurBATs(i),[],'extrap');
end
%
CurKtranses=handles.CurPKs(:,KtransIdx);
CurKtranses(isnan(CurKtranses))=0;
for i=1:numel(Idxs)
%     Regressors=[SAIF(PKs(c,1),:); squeeze(SHConvd2(PKs(c,1),c,:))'];
    Regressors=[SHSAIF(i,:); squeeze(SHConvd2(i,:))];
        % BATfinal VpFinal KtransFinal Kepfinal VeFinal
    Sims(i,:)=((Regressors')*([handles.CurPKs(i,[VpIdx]) CurKtranses(i)]'));        
%     figure;plot(handles.SampleTs,handles.CTC2D(Idxs,:),'k*',handles.SampleTs,Sims,'r',handles.SampleTs,SHSAIF*handles.PKs(Idxs(i),3),'g',handles.SampleTs,SHConvd2*handles.PKs(Idxs(i),4)/dt,'m')

%     Regressors=[handles.HSAIF(handles.MIdxs(Idxs(i),1),:); squeeze(handles.HSHConvd(handles.MIdxs(Idxs(i),1),handles.Keps1I(handles.MIdxs(Idxs(i),2)),:))'];
%     Sims(i,:)=((Regressors')*handles.CXs(:,Idxs(i)));
%     for j=1:numel(handles.Titles)
%         Strs{j}=[Strs{j} ' ' num2str(handles.Vols{j}(ChosenVoxels(i,1),ChosenVoxels(i,2),ChosenVoxels(i,3)))];
%     end
end
%%
SAIF=zeros([handles.nTDif numel(handles.SampleTs)]);
CSAIF=zeros([handles.nTDif numel(handles.SampleTs)]);
CHAIF = cumtrapz(handles.HSampleTs,handles.HAIF);
for i=1:handles.nTDif
    SAIF(i,:)=interp1(handles.HSampleTs,handles.HAIF,handles.SampleTs+handles.TDif(i),[],'extrap');
    CSAIF(i,:)=interp1(handles.HSampleTs,CHAIF,handles.SampleTs+handles.TDif(i),[],'extrap');
end
%%
for i=1:size(Sims,1)
    SSims(i,:)=interp1(handles.HSampleTs,Sims(i,:),handles.SampleTs,[],'extrap');
end
%%
DTs=diff(handles.SampleTs);
nRegTimePoints=find(DTs>DTs(1)*2,1);
nExtraTPs=size(handles.SampleTs,2)-nRegTimePoints;
%% Extract for simulated

%    1        2         3        4          5     6        7    8   9    10
% BATfinal VpFinal KtransFinal Kepfinal VeFinal RSSFinal RSS0 RSS1 RSS2 RSS3
for k=1:nExtraTPs+1
    PKs_check(k,:,:) = FindPKBATgAIFMuraseF4Models_TProb(SSims(:,1:end-nExtraTPs-1+k),SAIF(:,1:end-nExtraTPs-1+k),handles.SampleTs(:,1:end-nExtraTPs-1+k),CSAIF(:,1:end-nExtraTPs-1+k));
end

for k=1:nExtraTPs+1
    SVps(:,k)=PKs_check(k,:,VpIdx);
    SKtranss(:,k)=PKs_check(k,:,KtransIdx);
    SKeps(:,k)=PKs_check(k,:,KepIdx);
end
figure;plot(SKeps','*');hold on;plot(repmat(CurKeps,[1 nExtraTPs+1])','r-')
title('Real values - red, Blue - Estimated values')