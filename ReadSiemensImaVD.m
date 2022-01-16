function [data, Mag, Phase, Real, Imag, rows, cols] = ReadSiemensImaVD(ima_file, dispopt,ifShow)
% [data, Mag, Phase, Real, Imag, rows, cols] = ReadSiemensIma(ima_file, dispopt)
% Read .ima and .IceHeader files obtained by ICE simulation.The names of the .ima and .IceHead
% have either the same base, or
%    WriteToFile_0001.ima <---> MiniHead_0001.IceHead
%    WriteToFile_0002.ima <---> MiniHead_0002.IceHead
%    WriteToFile_0003.ima <---> MiniHead_0003.IceHead
%                  ... ...
%
% EXAMPLE:
%        [data Mag Phase Real Imag rows cols] = ReadSiemensIma(ima_file);
%
% Maolin Qiu YALE 12-02-2007 maolin.qiu(at)yale(dot)edu
%
if nargin < 1
    help ReadSiemensIma;
    [filename pathname] = uigetfile( ...
       {'*.ima';'*.IceHead';'*.*'}, ...
        'Pick a Siemens IMA(IceHead) file');
    if ~filename & ~pathname
        disp(['You selected no file.']);
        return;
    else
        ima_file = fullfile(pathname, filename);
    end
end

if nargin < 2
    dispopt = 'on';
end

[pa na] = fileparts(ima_file);
header_file = fullfile(pa, [na '.IceHead']);
if exist(header_file, 'file')
    disp(['IceHead file: ' header_file]);
else
    header_file = fullfile(pa, ['MiniHead_ima_' na(end-4:end) '.IceHead']);
    if exist(header_file, 'file')
        disp(['IceHead file: ' header_file]);
    else
        disp(['Measurement data file does not exist: ' header_file]);
        return;
    end
end

ima_file = fullfile(pa, [na '.ima']);
if exist(ima_file, 'file')
    disp(['Image data file: ' ima_file]);
else
    ima_file = fullfile(pa, ['WriteToFile_' na(end-4:end) '.ima']);
    if exist(ima_file, 'file')
        disp(['Image data file: ' ima_file]);
    else
        disp(['Image data file does not exist: ' ima_file]);
        return;
    end
end

Mag = []; Phase = []; Real = []; Imag = [];

[cols, rows] = parseIceHead(header_file);

fid = fopen(ima_file, 'r');
[data N] = fread(fid, inf, 'int16');
if N == cols*rows
  data = reshape(data, cols, rows);
  data = data';
  if(ifShow)
    figure, imagesc(data), axis image, colormap(gray);
  end
else
  frewind(fid);
  [data N] = fread(fid, inf, 'float');
  if N == cols*rows
    data = reshape(data, cols, rows);
    data = data';
    figure, imagesc(data), axis image, colormap(gray);
  elseif N == 2*cols*rows
    Real = reshape(data(1:2:end), cols, rows)';
    Imag = reshape(data(2:2:end), cols, rows)';
    subplot(2,3,1), imagesc(Real), axis image, colormap(gray);
    subplot(2,3,4), imagesc(Imag), axis image, colormap(gray);    
    data = complex(Real, Imag);
    subplot(2,3,2), imagesc(abs(data)), axis image, colormap(gray);
    subplot(2,3,5), imagesc(angle(data)), axis image, colormap(gray);    
    I = fftshift(fft2(fftshift(data)));
    Mag = abs(I);
    Phase = angle(I);
    subplot(2,3,3), imagesc(Mag), axis image, colormap(gray);    
    subplot(2,3,6), imagesc(Phase), axis image, colormap(gray);
  end
end
fclose(fid);

function [nCols, nRows] = parseIceHead(header_file)
    fid = fopen(header_file, 'r');
    
    
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(strfind(line, '<ParamLong."NoOfCols">'))
            discard = fgetl(fid);
            discard = fgetl(fid);
            nCols = str2num(strtrim(fgetl(fid)));
        end
        if ~isempty(strfind(line, '<ParamLong."NoOfRows">'))
            discard = fgetl(fid);
            discard = fgetl(fid);
            nRows = str2num(strtrim(fgetl(fid)));
        end
    end
    
    fclose(fid);
    
