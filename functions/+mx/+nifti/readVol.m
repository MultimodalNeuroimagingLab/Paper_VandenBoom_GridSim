% 
%   Read a nifti file by filepath (uses SPM to read)
%
%   [volumeData, volume, xyzData] = readVol(filepath)
%
%   filepath  		= the filepath to the nifti to be read
%
%   Returns: 
%       volumeData  = 3D or 4D matrix of image data
%       volume      = a structure array containing image volume information
%       xyzData     = 3xn matrix of XYZ locations returned (in mm)
%
%   Copyright (c) 2016, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [volumeData, volume, xyzData] = readVol(filepath)
    
    % set 0 to return no data
    volumeData = [];
    volume = [];
    
    % check if the input file exist
    if exist(filepath, 'file') == 0
        fprintf(2, ['Error: invalid filepath to nifti volume (', filepath, ')\n']);
        return;
    end
    
    % read the file
    volume = spm_vol(filepath);
    if nargout > 2
        [volumeData, xyzData] = spm_read_vols(volume);
    else
        volumeData = spm_read_vols(volume);
    end
    
    % message
    disp(['Volume ', filepath, ' was read succesfully']);

end