
% get data and store them in some matrices
function data = read_raw_data(filename_data, NCol,NMeas,inplanLine, mb, NSlc,NRefLine)
% filename_data = 'meas_MID28_cmrr_Grp2_PF7_8_mb4_disf0_tr1_mtx74_te19_ref32_FID3152.dat';
os = 2;
szChannelHeader = 128;
szScanHeader = 0;
NCha = 32;
NCol = NCol*os;
nByte = NCha*(szChannelHeader + 8*NCol);
sz = [2 nByte/8];
shape = [NCol+szChannelHeader/8, NCha];
cut   = szChannelHeader/8 + (1 : NCol);
fileID_data = fopen(filename_data,'r','ieee-le');
fseek(fileID_data,0,'eof');
fileSize_data = ftell(fileID_data);
measLength_data = fileSize_data;
measOffset = 0;
version = 'vb';
isInplacePATref = false;
fseek(fileID_data,0,'bof');
hdr_len_data = fread(fileID_data,1,'uint32');
fseek(fileID_data,hdr_len_data,'bof');

[mdh_blob_data,filePos_data,isEOF] = loop_mdh_read(fileID_data, version, 1, 1, measOffset, measLength_data);
[mdh_data, mask_data] = evalMDH( mdh_blob_data, version, isInplacePATref);
%%
clear image_data
clear other_data
clear phasecor_data
clear refscan_data
line_pc = 1;
line_img = 1;
line_others = 1;
line_rs = 1;
for k = 1:length(filePos_data)-1
    fseek(fileID_data,filePos_data(k),'bof');
    raw = fread(fileID_data,sz,'float=>single');
    new_raw = single(zeros([sz(1) sz(2)]));
    temp_real = raw(1,:);
    temp_imag = raw(2,:);
    temp_raw = complex(temp_real,temp_imag)';
    temp_raw = reshape(temp_raw,shape);
    temp_raw = temp_raw(cut,:);
    %temp_new_real = zeros(1, size(new_raw,2));
    %temp_new_imag = zeros(1, size(new_raw,2));
%     temp_new_real = temp_real(1:2:end);
%     temp_new_imag = temp_imag(1:2:end);
%     if(ismember(k, find(mask_data.MDH_NOISEADJSCAN == 1)))
%         noise_scan = temp_raw(1:2:end,:);
%         %line_noise = line_noise + 1;
%     end
%     if(ismember(k, find(mask_data.MDH_PATREFSCAN == 1)))
%         refscan(:,:,line_ref) = temp_raw(1:2:end,:);
%         line_ref = line_ref + 1;
%     end
%     if(ismember(k, find(mask_data.MDH_PHASCOR == 1)))
%         phase_cor(:,:,line_PC) = temp_raw(1:2:end,:);
%         line_PC = line_PC + 1;
%     end
    if(ismember(k, find(mask_data.MDH_IMASCAN == 1)))
%         new_temp_raw = single(zeros(size(temp_raw,1),size(temp_raw,2)));
%         for ch = 1:size(temp_raw,2)
%             input = temp_raw(:,ch);
%             %new_temp_raw(:,ch) = reweight(input,F1d_rounded_cropped,os);
%             new_temp_raw(:,ch) = input;
%         end
        image_data(:,:,line_img) = temp_raw;
        line_img = line_img + 1;
        % need to do something
    elseif(ismember(k,find(mask_data.MDH_PATREFSCAN ==1)))
       refscan_data(:,:,line_rs) = temp_raw;
       line_rs = line_rs +1;
    elseif(ismember(k,find(mask_data.MDH_PHASCOR ==1)))
        phasecor_data(:,:,line_pc) = temp_raw;
        line_pc = line_pc + 1;
    else
        other_data(:,:,line_others) = temp_raw;
        line_others = line_others + 1;
    end
    

end
image_data_p1 = image_data(:,:,1:inplanLine*NSlc);
image_data_p2 = image_data(:,:,inplanLine*NSlc + 1:end);
data.image_data_p1 = reshape(image_data_p1, [size(image_data_p1,1),size(image_data_p1,2),inplanLine,NSlc]);
data.image_data_p2 = reshape(image_data_p2, [size(image_data_p2,1),size(image_data_p2,2),inplanLine,NSlc/mb,NMeas]);
data.other_data = other_data;
% data.phasecor_data = reshape;
data.phasecor_data_p1 = reshape(phasecor_data(:,:,1:3*NSlc),[size(phasecor_data,1),size(phasecor_data,2),3,NSlc]);
data.phasecor_data_p2 = reshape(phasecor_data(:,:,3*NSlc+1:end),[size(phasecor_data,1),size(phasecor_data,2),3,NSlc/mb,NMeas]);
% data.image_data = image_data;
data.refscan_data = reshape(refscan_data,[size(refscan_data,1),size(refscan_data,2),NRefLine,NSlc]);
end
%%
function [mdh_blob, filePos, isEOF] = loop_mdh_read( fid, version, Nscans, scan, measOffset, measLength)
% Goal of this function is to gather all mdhs in the dat file and store them
% in binary form, first. This enables us to evaluate and parse the stuff in
% a MATLAB-friendly (vectorized) way. We also yield a clear separation between
% a lengthy loop and other expressions that are evaluated very few times.
%
% The main challenge is that we never know a priori, where the next mdh is
% and how many there are. So we have to actually evaluate some mdh fields to
% find the next one.
%
% All slow things of the parsing step are found in the while loop.
% => It is the (only) place where micro-optimizations are worthwhile.
%
% The current state is that we are close to sequential disk I/O times.
% More fancy improvements may be possible by using workers through parfeval()
% or threads using a java class (probably faster + no toolbox):
% http://undocumentedmatlab.com/blog/explicit-multi-threading-in-matlab-part1

    switch version
        case 'vb'
            isVD    = false;
            byteMDH = 128;
        case 'vd'
            isVD    = true;
            byteMDH = 184;
            szScanHeader    = 192; % [bytes]
            szChannelHeader =  32; % [bytes]
        otherwise
            % arbitrary assumptions:
            isVD    = false;
            byteMDH = 128;
            warning( [mfilename() ':UnknownVer'], 'Software version "%s" is not supported.', version );
    end

    cPos          = ftell(fid);
    n_acq         = 0;
    allocSize     = 4096;
    ulDMALength   = byteMDH;
    isEOF         = false;
    last_progress = 0;

    mdh_blob = zeros( byteMDH, 0, 'uint8' );
    szBlob   = size( mdh_blob, 2 );
    filePos  = zeros(0, 1, class(cPos));  % avoid bug in Matlab 2013b: https://scivision.co/matlab-fseek-bug-with-uint64-offset/

    fseek(fid,cPos,'bof');

    % ======================================
    %   constants and conditional variables
    % ======================================
        bit_0 = uint8(2^0);
        bit_5 = uint8(2^5);
        mdhStart = 1-byteMDH;
        
        u8_000 = zeros( 3, 1, 'uint8'); % for comparison with data_u8(1:3)

        % 20 fill bytes in VD (21:40)
        evIdx   = uint8(    21  + 20*isVD); % 1st byte of evalInfoMask
        dmaIdx  = uint8((29:32) + 20*isVD); % to correct DMA length using NCol and NCha
        if isVD
            dmaOff  = szScanHeader;
            dmaSkip = szChannelHeader;
        else
            dmaOff  = 0;
            dmaSkip = byteMDH;
        end
    % ======================================

    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave % octave does not support a cancel button
        h = waitbar(0,'','Name', sprintf('Reading Scan ID %d/%d', scan, Nscans));
    else
        h = waitbar(0,'','Name', sprintf('Reading Scan ID %d/%d', scan, Nscans),...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0)
    end

    t0 = tic;
    while true
        % Read mdh as binary (uint8) and evaluate as little as possible to know...
        %   ... where the next mdh is (ulDMALength / ushSamplesInScan & ushUsedChannels)
        %   ... whether it is only for sync (MDH_SYNCDATA)
        %   ... whether it is the last one (MDH_ACQEND)
        % evalMDH() contains the correct and readable code for all mdh entries.
         
        try
            % read everything and cut out the mdh
            data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
            data_u8 = data_u8( mdhStart+end :  end );
        catch exc
            warning( [mfilename() ':UnxpctdEOF'],  ...
                      [ '\nAn unexpected read error occurred at this byte offset: %d (%g GiB)\n'...
                        'Will stop reading now.\n'                                             ...
                        '=== MATLABs error message ================\n'                         ...
                        exc.message                                                            ...
                        '\n=== end of error =========================\n'                       ...
                       ], cPos, cPos/1024^3 )
            isEOF = true;
            break
        end

        if ~isOctave && getappdata(h,'canceling') 
            break;
        end
        
        bitMask = data_u8(evIdx);   % the initial 8 bit from evalInfoMask are enough

        if   isequal( data_u8(1:3), u8_000 )    ... % probably ulDMALength == 0
          || bitand(bitMask, bit_0)                 % MDH_ACQEND

            % ok, look closer if really all *4* bytes are 0:
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );

            if ulDMALength == 0 || bitand(bitMask, bit_0)
                cPos = cPos + ulDMALength;
                % jump to next full 512 bytes
                if mod(cPos,512)
                    cPos = cPos + 512 - mod(cPos,512);
                end
                break;
            end
        end
        if bitand(bitMask, bit_5)   % MDH_SYNCDATA
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
            cPos = cPos + ulDMALength;
            continue
        end

        % pehses: the pack bit indicates that multiple ADC are packed into one
        % DMA, often in EPI scans (controlled by fRTSetReadoutPackaging in IDEA)
        % since this code assumes one adc (x NCha) per DMA, we have to correct
        % the "DMA length"
        %     if mdh.ulPackBit
        % it seems that the packbit is not always set correctly
        NCol_NCha = double( typecast( data_u8(dmaIdx), 'uint16' ) );  % [ushSamplesInScan  ushUsedChannels]
        ulDMALength = dmaOff + (8*NCol_NCha(1) + dmaSkip) * NCol_NCha(2);

        n_acq = n_acq + 1;

        % grow arrays in batches
        if n_acq > szBlob
            mdh_blob( :, end + allocSize ) = 0;
            filePos( end + allocSize ) = 0;
            szBlob = size( mdh_blob, 2 );
        end
        mdh_blob(:,n_acq) = data_u8;
        filePos( n_acq )  = cPos;

        progress = (cPos-measOffset)/measLength;
        
        if progress > last_progress  + 0.01
            last_progress = progress;
            elapsed_time  = toc(t0);
            time_left     = elapsed_time * (1/progress-1);
            progress_str  = sprintf('%3.0f %% read in %4.0f s;\nestimated time left: %4.0f s', round(100*progress), elapsed_time, time_left);
            waitbar(progress, h, progress_str);
        end

        cPos = cPos + ulDMALength;
    end % while true
    delete(h);
    
    if isEOF
        n_acq = n_acq-1;    % ignore the last attempt
    end

    filePos( n_acq+1 ) = cPos;  % save pointer to the next scan

    % discard overallocation:
    mdh_blob = mdh_blob(:,1:n_acq);
    filePos  = reshape( filePos(1:n_acq+1), 1, [] ); % row vector

    fprintf('%8.1f MB read in %4.0f s\n', measLength/1024^2, round(toc(t0)));

end % of loop_mdh_read()

function [mdh,mask] = evalMDH( mdh_blob, version, isInplacePATref )
% see pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h
% and pkg/MrServers/MrMeasSrv/SeqIF/MDH/MdhProxy.h

if ~isa( mdh_blob, 'uint8' )
    error([mfilename() ':NoInt8'], 'mdh data must be a uint8 array!')
end

if version(end) == 'd'
    isVD = true;
    mdh_blob = mdh_blob([1:20 41:end], :);  % remove 20 unnecessary bytes
else
    isVD = false;
end

Nmeas   = size( mdh_blob, 2 );

mdh.ulPackBit   = bitget( mdh_blob(4,:), 2).';
mdh.ulPCI_rx    = bitset(bitset(mdh_blob(4,:), 7, 0), 8, 0).'; % keep 6 relevant bits
mdh_blob(4,:)   = bitget( mdh_blob(4,:),1);  % ubit24: keep only 1 bit from the 4th byte

% unfortunately, typecast works on vectors, only
data_uint32     = typecast( reshape(mdh_blob(1:76,:),  [],1), 'uint32' );
data_uint16     = typecast( reshape(mdh_blob(29:end,:),[],1), 'uint16' );
data_single     = typecast( reshape(mdh_blob(69:end,:),[],1), 'single' );

data_uint32 = reshape( data_uint32, [], Nmeas ).';
data_uint16 = reshape( data_uint16, [], Nmeas ).';
data_single = reshape( data_single, [], Nmeas ).';
                                                        %  byte pos.
%mdh.ulDMALength               = data_uint32(:,1);      %   1 :   4
mdh.lMeasUID                   = data_uint32(:,2);      %   5 :   8
mdh.ulScanCounter              = data_uint32(:,3);      %   9 :  12
mdh.ulTimeStamp                = data_uint32(:,4);      %  13 :  16
mdh.ulPMUTimeStamp             = data_uint32(:,5);      %  17 :  20
mdh.aulEvalInfoMask            = data_uint32(:,6:7);    %  21 :  28
mdh.ushSamplesInScan           = data_uint16(:,1);      %  29 :  30
mdh.ushUsedChannels            = data_uint16(:,2);      %  31 :  32
mdh.sLC                        = data_uint16(:,3:16);   %  33 :  60
mdh.sCutOff                    = data_uint16(:,17:18);  %  61 :  64
mdh.ushKSpaceCentreColumn      = data_uint16(:,19);     %  66 :  66
mdh.ushCoilSelect              = data_uint16(:,20);     %  67 :  68
mdh.fReadOutOffcentre          = data_single(:, 1);     %  69 :  72
mdh.ulTimeSinceLastRF          = data_uint32(:,19);     %  73 :  76
mdh.ushKSpaceCentreLineNo      = data_uint16(:,25);     %  77 :  78
mdh.ushKSpaceCentrePartitionNo = data_uint16(:,26);     %  79 :  80

if isVD
    mdh.SlicePos                    = data_single(:, 4:10); %  81 : 108
    mdh.aushIceProgramPara          = data_uint16(:,41:64); % 109 : 156
    mdh.aushFreePara                = data_uint16(:,65:68); % 157 : 164
else
    mdh.aushIceProgramPara          = data_uint16(:,27:30); %  81 :  88
    mdh.aushFreePara                = data_uint16(:,31:34); %  89 :  96
    mdh.SlicePos                    = data_single(:, 8:14); %  97 : 124
end

% inlining of evalInfoMask
evalInfoMask1 = mdh.aulEvalInfoMask(:,1);
mask.MDH_ACQEND            = min(bitand(evalInfoMask1, 2^0), 1);
mask.MDH_RTFEEDBACK        = min(bitand(evalInfoMask1, 2^1), 1);
mask.MDH_HPFEEDBACK        = min(bitand(evalInfoMask1, 2^2), 1);
mask.MDH_SYNCDATA          = min(bitand(evalInfoMask1, 2^5), 1);
mask.MDH_RAWDATACORRECTION = min(bitand(evalInfoMask1, 2^10),1);
mask.MDH_REFPHASESTABSCAN  = min(bitand(evalInfoMask1, 2^14),1);
mask.MDH_PHASESTABSCAN     = min(bitand(evalInfoMask1, 2^15),1);
mask.MDH_SIGNREV           = min(bitand(evalInfoMask1, 2^17),1);
mask.MDH_PHASCOR           = min(bitand(evalInfoMask1, 2^21),1);
mask.MDH_PATREFSCAN        = min(bitand(evalInfoMask1, 2^22),1);
mask.MDH_PATREFANDIMASCAN  = min(bitand(evalInfoMask1, 2^23),1);
mask.MDH_REFLECT           = min(bitand(evalInfoMask1, 2^24),1);
mask.MDH_NOISEADJSCAN      = min(bitand(evalInfoMask1, 2^25),1);
mask.MDH_VOP               = min(bitand(mdh.aulEvalInfoMask(2), 2^(53-32)),1); % was 0 in VD
mask.MDH_IMASCAN           = ones( Nmeas, 1, 'uint32' );


noImaScan = ( mask.MDH_ACQEND           | mask.MDH_RTFEEDBACK   | mask.MDH_HPFEEDBACK    ...
            | mask.MDH_PHASCOR          | mask.MDH_NOISEADJSCAN | mask.MDH_PHASESTABSCAN ...
            | mask.MDH_REFPHASESTABSCAN | mask.MDH_SYNCDATA );

if ~isInplacePATref
    noImaScan = (noImaScan | (mask.MDH_PATREFSCAN & ~mask.MDH_PATREFANDIMASCAN));
end
            
mask.MDH_IMASCAN( noImaScan ) = 0;

end 