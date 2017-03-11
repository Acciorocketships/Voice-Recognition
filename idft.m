function [x,t]=idft(X,fs,varargin)
% idft Inverse discrete Fourier transform.
%	[t,x]=idft(X,fs) is the inverse discrete Fourier transform of
%	frequency component X sampled at frequency sampling fs. The output
%	consists of vector t which has the time index (x axis) and vector x
%	which contains the values of x(t) (y axis). The length of the
%	vector x is equal to the number of frequency components in K.
%
%	[t,x]=idft(X,fs,N) computes the iDFT of frequency components X
%	sampled at frequency fs and outputs a vector x of length N which
%	contains the samples of x(t) sampled at sampling frequency fs and
%	of duration T=N/fs.
%
%	[t,x]=idft(X,fs,'centered') computes the iDFT of frequency
%	components X where it is assumed that these components have been
%	centered around k=0 and which contains the X(k) for k=0 (see output
%	X_c of dft function; help dft).
%
%	[t,x]=idft(X,fs,N,'centered') computes the iDFT of frequency
%	components X where it is assumed they are centered around k=0 (see
%	above) and outputs a vector x with N samples of x(t).
%
%	See also ifft, fft, dft.

% Created by Fernando Gama (fgama@seas.upenn.edu) for ESE-224 course. Spring
% 2016.

%%\\\\\ START: SANITY CHECK OF INPUTS.

% Check that the inputs are correct
if ~isnumeric(fs) % Check if fs is a number.
	error_msg='The sampling frequency should be a positive number.';
	error(error_msg);
end
if ~isnumeric(X) % Check if x is a number.
	error_msg='The input signal must be a vector';
	error(error_msg);
end
if sum(ismember(size(X),1))~=1 % Check if it is actually a vector 
	% (and not a matrix).
	error_msg='The input signal must be a vector';
	error(error_msg);
end
X=reshape(X,length(X),1); % Turn X into a column vector.
% First of all, we need to compute the length N of the signal x:
K=length(X);
% Check if an input K was used:
if ~isempty(varargin) % If the third argument is there
	if length(varargin)>2
		error_msg='Too many input arguments.';
		error(error_msg);
	end
	if length(varargin)==2
		N=varargin{1}; % Assign it to K because it is supposed to be that
		% parameter that is set on the third argument.
		situation=varargin{2}; % Saves if it has to be centered or not.
	end
	if length(varargin)==1
		aux=varargin{1};
		if ischar(aux)
			situation=aux;
			N=K;
		elseif isnumeric(aux)
			N=aux;
			situation=[];
		else
			N='error';
		end	
	end
	if ~isnumeric(N) % If the input is not a number, we get an error.
		error_msg=['Third argument should be a positive ' ...
			'integer corresponding to the number of DFT' ...
			' coefficients or the legend centered. See' ...
			' help idft.'];
		error(error_msg);
	end
	if N<=0
		error_msg=['The number of frequencies of the DFT has to be'...
			' a positive number'];
		error(error_msg);
	end
	if ~isempty(situation)
		if ~strcmp(situation,'centered')
			error_msg=['Unidentified request for computing iDFT.'...
				' See help idft'];
			error(error_msg);
		end
	end
else
	N=K; % If no specific number of frequencies was selected, use the same
		% as the length of the signal.
	situation=[];
end
% Frequency cannot be zero
if ~fs
	error('Frequency cannot be zero.');
end
% If the sampling frequency is a negative number, just use the positive.
if fs<0
	warning_msg='Absolute value of sampling frequency is used instead';
	warning(warning_msg);
	fs=-fs;
end

%%\\\\\ END: SANITY CHECK OF INPUTS.



%%\\\\\ START: COMPUTATION OF DFT.

% Several methods are used to output such DFT. Only the most efficient one is
% left uncommented. But the rest of them are included for educative purposes.

% The computation of times is rather straightforward. If we know that we
% took samples every Ts=1/fs seconds, and we know we have N samples, then,
% the duration of the signal is T=N*Ts=N/fs. Then, we need to fit, starting
% from t=0, N samples, one every time Ts. Then,
T=N/fs;
t=(0:1/fs:T-1/fs)';

% \\\\\ CENTERED FFT.
% If the DFT was centered, then just uncenter it before doing all the idft
% processing.
if strcmp(situation,'centered')
	X=circshift(X,-(ceil(K/2)-1));
end

% \\\\\ METHOD 1: For loops
% The for loops methos is pretty analogous to the one adopted in Method 1
% of dft() function. Try to adapt that one to see if you understood how the
% for loops method works. For that, you'll need to recall the iDFT
% definition (eq. 1, lab 3).
% x(n) = 1/sqrt(K) sum_{k=0 to K-1} X(k) exp(j * 2pi * k * n / N);

% \\\\ METHOD 2: Matrix form
% Using for loops in matlab is rather expensive. Matlab is already optimized to
% work with vectors and matrices, so it's a good idea to take advantage of this,
% together with the built-in elementwise operators.
% Following the mathematical development from Method 2 as explained in the dft.m
% file, adapt it to perform the matrix computation of the IDFT. Analyze if K and
% N should be treated any differently. Is there any relationship between the DFT
% matrix WKN and the iDFT matrix?

% \\\\\ METHOD 3: Built-in fft() function
% Even though the matrix form is fast, it is still not fast enough for large
% signals x. For that, it is better to use the built in fft() function which is
% the optimal way to compute a dft. Besides, it is really easy to code:

x=ifft(X,N)*sqrt(N);

%%\\\\\ END: COMPUTATION OF DFT.
end