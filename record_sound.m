function x=record_sound(n,fs,b,channel)
%record_sound Record sound from soundcard.
%	x=record_sound(n,fs,bits,channel) record n seconds of sound taken from
%	channel at sample frequency fs coding the samples in b bits. The output
%	x contains the samples recorded.
%
%	If any of the parameters fs, b, or channel are left empty then defaults
%	are used. The defaults are:
%		fs=8192Hz;
%		b=8;
%		channel=1;
%
%	Obs.: Calling the function with less than four inputs will result in an
%	error. If any input does not want to specified please leave it empty by
%	using [].
%
%	Example:
%		x=record_sound(3,[],[],[])
%	will record 3 second of the sound at default values of fs, b and
%	channel.
%
%	See also audiorecorder, sound.

% Created by Fernando Gama (fgama@seas.upenn.edu) for ESE-224 course. Spring
% 2016.

if nargin < 4
    channel=1;
    if nargin < 3
        b = 8;
        if nargin < 2
            fs = 8192;
        end
    end
end

if ~isnumeric(n) || ~isscalar(n)
	error_msg='Wrong input duration n.';
	error(error_msg);
end
if n<0
	error_msg='Time cannot be negative.';
	error(error_msg);
end
if ~n
	x=[];
	warning_msg='Zero duration selected. Returning empty vector.';
	warning(warning_msg);
	return;
end
if ~isnumeric(fs)
	error_msg='Wrong input sampling frequency fs.';
	error(error_msg);
end
if isempty(fs)
	fs=8192;
elseif ~isscalar(fs)
	error_msg='Wrong input sampling frequency fs.';
	error(error_msg);
elseif ~fs
	x=[];
	warning_msg='Sampling frequency selected to be 0. No signal returned';
	warning(warning_msg);
	return;
elseif fs<0
	warning_msg='Absolute value of sampling frequency is used instead';
	warning(warning_msg);
	fs=-fs;
end
if ~isnumeric(b)
	error_msg='Wrong input bits b.';
	error(error_msg);
end
if isempty(b)
	b=8;
elseif ~isscalar(b)
	error_msg='Wrong input bits b.';
	error(error_msg);
elseif mod(b,1)
	error_msg='Number of bits to code must be an integer.';
	error(error_msg);
elseif b<0
	warning_msg='Negative number of bits detected. Taking absolute value.';
	warning(warning_msg);
	b=-b;
elseif ~b
	warning_msg=['Number of bits to code is set to zero. '...
		'Returning empty vector.'];
	warning(warning_msg);
	x=[];
	return;
end
if ~isnumeric(channel)
	error_msg='Wrong input of channel to use.';
	error(error_msg)
end
if isempty(channel)
	channel=1;
elseif ~isscalar(channel)
	error_msg='Wrong input of channel to use.';
	error(error_msg)
elseif channel<=0
	error_msg=['Wrong channel selected. Check your PC configuration '...
		'available channels.'];
	error(error_msg);
elseif mod(channel,1)
	error_msg='The selected channel should be a positive integer.';
	error(error_msg);
end

my_sound = audiorecorder(fs, b, channel);
fprintf('Start speaking.\n')
% record for 10 seconds before moving on    
recordblocking(my_sound, n);
fprintf('End of Recording.\n');
x=getaudiodata(my_sound)';
x=reshape(x,length(x),1);
end