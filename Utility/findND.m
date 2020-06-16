function varargout=findND(X,varargin)
%Find non-zero elements in ND-arrays. Replicates all behavior from find.
% The syntax is equivalent to the built-in find, but extended to
% multi-dimensional input.
%
% [...] = findND(X,K) returns at most the first K indices. K must be a
% positive  scalar of any type.
%
% [...] = findND(X,K,side) returns either the first K or the last K
% inidices. The input side  must be a char, either 'first' or 'last'. The
% default behavior is 'first'.
%
% [I1,I2,I3,...,In] = findND(X,...) returns indices along all the
% dimensions of X.
%
% [I1,I2,I3,...,In,V] = findND(X,...) returns indices along all the
% dimensions of X, and additionally returns a vector containg the values.
%
% Note for Matlab 6.5:
% The syntax with more than one input is present in the online doc for R14
% (Matlab 7.0), so this might be the latest release without support for
% this syntax.
%
% Compatibility:
% Matlab: should work on all releases (tested on R2018b, R2012b and 6.5)
% Octave: tested on 4.2.1
% OS:     should work cross-platform
%
% Version: 1.2
% Date:    2018-10-17
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 )
% Email=  'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})
%
% Changes from v1.0 to v1.1:
% - Added support for Matlab 6.5 (R13).
% - Fixed a minor bug where the orientation of the output vector did not
%   match the orientation of the built-in find function. It also always
%   returned 1 as value, instead of the true value.
% - Added support upper case third input.
%
% Changes from v1.1 to v1.2:
% - Fixed a minor bug where the findND(___,'last') syntax returns 1 element
%   too much on releases prior to version 7.
% - Fixed a minor issue where the version number would not be retrieved
%   correctly, which will become relevant when version 10.0 will be
%   released.
% - Some cosmetic changes.

%#dependencies{getversion}

% Copyright (c) 2018, Rik Wisselink
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of  nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%Parse inputs
if ~(isnumeric(X) || islogical(X)) || numel(X)==0
    error('HJW:findND:FirstInput',...
        ['Expected first input (X) to be a non-empty numeric ',...
        'or logical array.'])
end
switch nargin
    case 1
        %[...] = findND(X);
        side='first';
        K=inf;
    case 2
        %[...] = findND(X,K);
        side='first';
        K=varargin{1};
        if ~(isnumeric(K) || islogical(K)) || numel(K)~=1 || any(K<0)
            error('HJW:findND:SecondInput',['Expected second input (K)',...
                ' to be a positive numeric or logical scalar.'])
        end
    case 3
        %[...] = FIND(X,K,'first');
        K=varargin{1};
        if ~(isnumeric(K) || islogical(K)) || numel(K)~=1 || any(K<0)
            error('HJW:findND:SecondInput',['Expected second input (K)',...
                ' to be a positive numeric or logical scalar.'])
        end
        side=varargin{2};
        if ~isa(side,'char') ||...
                ~( strcmpi(side,'first') || strcmpi(side,'last'))
            error('HJW:findND:ThirdInput',...
                'Third input must be either ''first'' or ''last''.')
        end
        side=lower(side);
    otherwise
        error('HJW:findND:InputNumber','Incorrect number of inputs.')
end

%parse outputs
nDims=length(size(X));
%allowed outputs: 0, 1, nDims, nDims+1
if nargout>1 && nargout<nDims
    error('HJW:findND:Output','Incorrect number of output arguments.')
end

varargout=cell(nargout,1);
v=getversion;
if v<7
    %The find(X,k,side) syntax was introduced between 6.5 and 7
    if nargout>nDims
        [ind,col_index_equal_to_one,val]=find(X(:));%#ok no tilde pre-R2009b
        %X(:) converts X to a column vector. Treating X(:) as a matrix
        %forces val to be the actual value, instead of the column index.
        if length(ind)>K
            if strcmp(side,'first')
                %select first K outputs
                ind=ind(1:K);
                val=val(1:K);
            else
                %select last K outputs
                ind=ind((end-K+1):end);
                val=val((end-K+1):end);
            end
        end
        [varargout{1:(end-1)}] = ind2sub(size(X),ind);
        varargout{end}=val;
    else
        ind=find(X);
        if length(ind)>K
            if strcmp(side,'first')
                %select first K outputs
                ind=ind(1:K);
            else
                %select last K outputs
                ind=ind((end-K):end);
            end
        end
        [varargout{:}] = ind2sub(size(X),ind);
    end
else
    if nargout>nDims
        %Tilde (~) to ignore outputs was introduced in R2009b. It is
        %probably faster to ignore the extra output than to use an if.
        [ind,col_index_equal_to_one,val]=find(X(:),K,side);%#ok
        %X(:) converts X to a column vector. Treating X(:) as a matrix
        %forces val to be the actual value, instead of the column index.
        [varargout{1:(end-1)}] = ind2sub(size(X),ind);
        varargout{end}=val;
    else
        ind=find(X,K,side);
        [varargout{:}] = ind2sub(size(X),ind);
    end
end
end
function [v,isOctave]=getversion(Rxxxxab)
%Get current version or determine numbered version
%
%Without any input this function returns the version number for the current
%release. If run with a char input, this function determines the version
%number. The conversion is based on a manual list and therefore needs to be
%updated manually, so it might not be complete (although it is possible to
%load the list from Wikipedia).
%
%Don't use v=version;v=str2double(v(1:3)); as it is incomplete for several
%releases (like e.g. 7.14 or in the future 10.1).
%
% Version: 1.0
% Date:    2018-10-17
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 )
% Email=  'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

persistent v_num v_dict octave
if isempty(v_num)
    %test if Octave is used instead of Matlab
    octave=exist('OCTAVE_VERSION', 'builtin');
    
    %get current version number
    v_num=version;
    ind=strfind(v_num,'.');
    if numel(ind)==1,ind=numel(v_num);else,ind=ind(2)-1;end
    v_num=str2double(v_num(1:ind));
    
    %get dictionary to use for ismember
    v_dict={'R13' 6.5;'R13SP1' 6.5;'R13SP2' 6.5;'R14' 7.0;'R14SP1' 7.0;...
        'R14SP2' 7.0;'R14SP3' 7.1;'R2006a' 7.2;'R2006b' 7.3;...
        'R2007a' 7.4;'R2007b' 7.5;'R2008a' 7.6;'R2008b' 7.7;...
        'R2009a' 7.8;'R2009b' 7.9;'R2010a' 7.1;'R2010b' 7.11;...
        'R2011a' 7.12;'R2011b' 7.13;'R2012a' 7.14;'R2012b' 8.0;...
        'R2013a' 8.1;'R2013b' 8.2;'R2014a' 8.3;'R2014b' 8.4;...
        'R2015a' 8.5;'R2015b' 8.6;'R2016a' 9.0;'R2016b' 9.1;...
        'R2017a' 9.2;'R2017b' 9.3;'R2018a' 9.4;'R2018b' 9.5};
end
if nargin==0
    v=v_num;
	isOctave=octave;
else
    L=ismember(v_dict(:,1),Rxxxxab);
    if sum(L)~=1
        v=NaN;
        warning('HJW:getversion:NotInDict',...
            'The requested version is not in the hard-coded list.')
    else
        v=v_dict{L,2};
    end
end
end
