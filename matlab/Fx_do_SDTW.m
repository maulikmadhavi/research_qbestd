%% Code Information
% This MATLAB code is used to compute segmental DTW between qrycoef and
% refcoef. SDTW is executed by several execution of basic dtw module. Here,
% we can impose the constraints on warping window adjustment.
%% Input Output parameters
% % % % % % % Input % % % % % % %
% refcoef := feature vectors from reference waveform (ND X Nframes)
% qrycoef := feature vectors from query waveform (ND X Nframes)
% cost_DTW := costs for accumulated distance and warping path calculation.
% (deletion, substitution, insertion) respectively.
% Type_localdist := Type of local distance.
% Case: 's' --> Euclidean, 'i'--> Inner product based distance.
% Warping_adjustment := Percentage permissible warping window adjustment
% value.
% Maskflag := flag associated with masking of distance matrix. 1-->
% Masking, 0 --> No Masking
% inc := frame increment duration (in ms or in samples)

% % % % % Output % % % % % %
% startpos := Hypothetical starting position  (in ms or in samples)
% endpos := Hypothetical starting position  (in ms or in samples)
% dist := Effective DTW distance
% Type_localdist='s';% R=Warping_adjustment;% R=Warping_adjustment;



function [dist, startpos,endpos,DistMtrx]= Fx_do_SDTW(refcoef,qrycoef,Type_localdist,Warping_adjustment)
%#codegen
% coder.inline('never');

%% Local distance
D=0;
switch Type_localdist
    case('i'); % If "Inner Product" is used as local distance computation metric
        D=-log(inDistance(refcoef,qrycoef)); % Inner product based distance
    case('in'); % If "Inner Product" is used as local distance computation metric
        refcoef=refcoef./repmat(sum(refcoef.*refcoef),size(refcoef,1),1);
        qrycoef=qrycoef./repmat(sum(qrycoef.*qrycoef),size(qrycoef,1),1);
        D=-log(inDistance(refcoef,qrycoef)); % Inner product based distance
    case('s'); % If "Euclidean" is used as local distance computation metric
        D=sqDistance(refcoef,qrycoef);       % Euclidean distance
    case('k'); % If "KL distance" (symmetric version) as a local distance
        %         D=slmetric_pw(refcoef,qrycoef,'kldiv')+slmetric_pw(qrycoef,refcoef,'kldiv')';
        D = KL_symdistance(refcoef,qrycoef);
    case('b'); % If "Bhattacharya Distance"  as a local distance        
        D = bhattDistance(refcoef,qrycoef);
end

% D=(D-repmat(min(D),size(D,1),1))./(repmat(max(D),size(D,1),1)-repmat(min(D),size(D,1),1));


%% Preallocation
[N1,N2]=size(D);                     % Num ref X Num query
R=round(N2*Warping_adjustment/100);
Mask=masking_distMtrx_segDTW(N2,R);
DistMtrx=zeros(1,max(1,N1-N2-R));   % Preallocate to store the dtw distances
ep=zeros(1,max(1,N1-N2-R));   % Preallocate to store the end pos


%% Compute SDTW for all possible segments
if length(DistMtrx)==1
    startpos=1;
    endpos=size(D,1);
    dist=DTW_c_basic_skel_nobt(D);
    DistMtrx=dist;
else
    for k1=1:1:max(N1-N2-R,1)     % Run FOR loop for every segments to do segmental DTW
        SF=k1;                  % Starting frame of the segmental DTW
        EF=min(k1+N2+R-1,size(D,1));    % Ending frame of the segmental DTW
        Dtemp=D(SF:EF,:);  % Temporary distance matrix formed by local distance computation
        [DistMtrx(1,k1),ep(1,k1)]=DTW_c_skel_nobt(Dtemp+Mask);  % Perform basic skeleton operation on 1 segment
    end
    % Optimum index
    [dist,optind]=min(DistMtrx);   % Get the optimum index of segmental DTW
    
    startpos=optind;
    endpos=optind+ep(optind);
    
end


