function [ image, meta ] = readflexcube(fname, varargin)
% READFLEXCUBE reads XYZC data from a Flex file
% [ cube, meta ] = READFLEXCUBE(fname)
% * fname is the filename of the flex file
% * channels is a numeric array with the indexes of the channels to read
%   Indices start with 1, if not specified read all channels
% * meta is the BioFormats metadata store object

parser = inputParser;
addParamValue(parser, 'Info', []);
addParamValue(parser, 'PlaneCount', []);
addParamValue(parser, 'ChannelCount', []);
addParamValue(parser, 'PixelSizeZ', []);
parse(parser, varargin{:});
args = parser.Results;

[~, f] = fileattrib(fname);
fname = f.Name;

if not(isempty(args.Info))
    info = args.Info;
else
    info = imfinfo(fname);
end

sizeX = info(1).Width;
sizeY = info(1).Height;
sizeZ = 1;
sizeC = 1;

if not(isempty(args.PlaneCount))
    sizeZ = args.PlaneCount;
end

if not(isempty(args.ChannelCount))
    sizeC = args.ChannelCount;
end

totalcount = numel(info);

if totalcount ~= sizeZ*sizeC
    if isempty(args.PlaneCount)
       sizeZ = floor(totalcount / sizeC);
    elseif isempty(args.ChannelCount)
       sizeC = floor(totalcount / sizeZ);
    end
end

cube = zeros([sizeY(), sizeX(), sizeZ, sizeC], 'uint16');

% getZCTCoords = @(i) [mod(i, sizeZ), ...
%                      mod(floor(i / sizeZ), sizeC), ...
%                      floor(floor(i / sizeZ) / sizeC)];

getCZTCoords = @(i) [mod(floor(i / sizeC), sizeZ), ...
                     mod(i, sizeC), ...
                     floor(floor(i / sizeC) / sizeZ)];

w = warning('off','all');
tf = Tiff(fname);
for i=1:totalcount
   zct = getCZTCoords(i - 1) + 1;
   cube(:, :, zct(1), zct(2)) = tf.read();
   if i < totalcount
       tf.nextDirectory();
   end
end
warning(w);

image.data = cube;
image.filename = fname;

xmlstr = info(1).UnknownTags(1).Value;
xml = xmlreadstring(xmlstr);

image.PixelSizeX = str2double(char(xml.getElementsByTagName('ImageResolutionX').item(0).getTextContent()));
image.PixelSizeY = str2double(char(xml.getElementsByTagName('ImageResolutionY').item(0).getTextContent()));

if not(isempty(args.PixelSizeZ))
    pixelZ = args.PixelSizeZ;
else
    pixelZ = nan;
end

image.PixelSizeZ = pixelZ;

meta = xml;

end

