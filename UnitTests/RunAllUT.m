%{ 
mwave - A water wave and wave energy converter computation package 
Copyright (C) 2014  Cameron McNatt

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contributors:
    C. McNatt
%}
%% runs all the unit tests

folder = [mwavePath '\UnitTests'];

tests = matlab.unittest.TestSuite.fromFolder(folder);

profile on;
res = run(tests);
profile off;

%% Collect failed/incomplete tests

N = length(res);

ifail = 0;
iincomp = 0;
for n = 1:N
    if res(n).Failed
        ifail = ifail + 1;
        failedTests(ifail) = tests(n);
    end
    if res(n).Incomplete
        iincomp = iincomp + 1;
        incompTests(iincomp) = tests(n);
    end
end
    

%% run specfic test to debug

test = computeMooringKUT;

run(test, 'test2')

%res = run(test, 'testIncCirWave');

%% run specfic test to debug

test = WaveCompUT;
profile on;
res = run(test);
profile off;

%% - Check the code coverage of a unit test
testClassName = 'ModesOfMotion';

test = eval([testClassName 'UT']);

fpath = fileparts(which(testClassName));

curdir = cd(fpath);

profile on;
res = run(test);
profile off;

% now go to the 'Current Folder' widget, click on the little down array
% select 'Reports' then 'Coverage Report'

%%

cd(curdir);