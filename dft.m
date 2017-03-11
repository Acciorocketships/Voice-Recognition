function [X_c,f_c,X,f]=dft(x,fs,varargin)
% dft Discrete Fourier transform.
%	[f,X,f_c,X_c]=dft(x,fs) is the discrete Fourier transform (DFT) of 
%	vector x. X is a vector that contains the DFT coefficients and f is a
%	vector that contains the real frequencies determined by sampling 
%	frequency fs. f are the frequencies starting at f=0 and X are the
%	corresponding frequency components. f_c is a vector containing the
%	frequencies such that f_c=0 is at the center, and X_c contains the
%	frequency componentes corresponding to f_c. In essence, f_c and X_c
%	are the centered counterparts of f and X. Frequency 0 is always
%	present.
%
%	[f,X,f_c,X_c]=dft(x,fs,K) computes the DFT with only K coefficients. 
%	Recall that this periodizes signal x with period K. If the length of
%	the signal x is less than K, then the signal will be padded with zeros.
%	If the length of the signal x is greater than K, then there will be
%	aliasing occured from periodizing singal x with period K. f_c and X_c
%	are the centered counterparts of f and X
%
%	See also fft, idft, ifft.

% Created by Fernando Gama (fgama@seas.upenn.edu) for ESE-224 course. Spring
% 2016.

%%\\\\\ START: SANITY CHECK OF INPUTS.

% Check that the inputs are correct
if ~isnumeric(fs) % Check if fs is a number.
	error_msg='The sampling frequency should be a positive number.';
	error(error_msg);
end
if ~isnumeric(x) % Check if x is a number.
	error_msg='The input signal must be a vector';
	error(error_msg);
end
if sum(ismember(size(x),1))~=1 % Check if it is actually a vector 
	% (and not a matrix).
	error_msg='The input signal must be a vector';
	error(error_msg);
end
x=reshape(x,length(x),1); % Turn x into a column vector.
% First of all, we need to compute the length N of the signal x:
N=length(x);
% Check if an input K was used:
if ~isempty(varargin) % If the third argument is there
	if length(varargin)>1
		error_msg='Too many input arguments.';
		error(error_msg);
	end
	K=varargin{1}; % Assign it to K because it is supposed to be that
		% parameter that is set on the third argument.
	if ~isnumeric(K) % If the input is not a number, we get an error.
		error_msg=['Third argument should be a positive integer ' ...
			'corresponding to the number of DFT coefficients.'];
		error(error_msg);
	end
	if K<=0
		error_msg=['The number of frequencies of the DFT has to be'...
			'a positive number'];
		error(error_msg);
	end
else
	K=N; % If no specific number of frequencies was selected, use the same
		% as the length of the signal.
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

% The computation of the frequencies is rather straightforward and does not
% change. If we take a look at (eq. 6, lab 1) we have that
% f_k=k/N*fs
% for each k=0,1,...,K-1. That is, for each discrete frequency k we can convert it
% to real time frequency according to this equation.
% Now, we know the output variable f is actually a vector, and that each element
% of this vector corresponds to a real frequency that was obtained by the
% previous equation. This can be done with a for loop. But that is not necessary
% if we think in terms of vectors: if we have a vector with discrete frequencies
% k=0,1,...,K-1, then we just need to multiply each element of the vector by fs/K
% (again, N is changed by K because K is the total number of frequencies we will
% be considering). That is a valid vector operation, so the fastest way to
% achieve f is just by multipling a vector of discrete frequencies by fs/K. That
% is:
f=(0:K-1)'*fs/K; % (0:K-1) just creates a vector from 0 to K by steps of 1. That is
	% (0:K-1)=[0 1 2 ... K-1]. By doing ' we are transposing the vector so as to
	% work with column vectors which is the standard.
f_c=(-ceil(K/2-1):floor(K/2))'*fs/K; % This accounts for the frequencies 
	% centered at zero. I want to be guaranteed that k=0 is always a
	% possible k. Then, I also have to account for both even and odd choices
	% of K, and that's why the floor() function appears to round down the
	% numbers.

% \\\\\ METHOD 1: For loops
% By definition of DFT (eq. 1, lab 2) we have that
% X(k) = 1/sqrt(N) * sum_{n=0 to N-1} x(n) exp(-j 2 pi k n / N)
% where n is the discrete time index, k is the discrete time frequency and N is
% the length of the time signal x. Observe that this gives the DFT coefficient
% for a single coefficient k. We are trying to look for all coefficients
% k=0,1,...,K.
% This means that, for each value of k=0,1,...,K we will need to compute N
% multiplications x(n)*exp(-j 2 pi k n / N), for n=0,1,...,N-1, and sum the
% result.

% First thing we need is to create the variable X that we will output. This will
% be a vector of length K.

%!!!!!! Uncomment the following for Method 1 to work:
% (Rember to comment other methods)

% X=zeros(K,1); % This is a mere initialization to "ask for memory", although it 
% 	% is not strictly necessary, it speeds up the process. It can be
% 	% initialized in any value since all values will be overwritten during
% 	% the for loops.
% for k=0:K-1 % For each time index k=0,1,...,K;
% 	for n=0:N-1 % For each frequency n=0,1,...,N-1:
% 		X(k+1)=X(k+1)+1/sqrt(N)*x(n+1)*exp(-1j*2*pi*k*n/K);
% 		% Obs. 1: the iterative variable k starts at 0, but the first
% 		% index of a vector in Matlab is 1, so we need to start with k+1
% 		% so that when k=0, the result is saved on the first element of
% 		% the vector which is for index=1.
% 		% Obs. 2: in the case we have K different from N, then the
% 		% signal will be periodized with period K. That is why the
% 		% exponential is divided by K instead of N.
% 	end
% end

% \\\\ METHOD 2: Matrix form
% Using for loops in matlab is rather expensive. Matlab is already optimized to
% work with vectors and matrices, so it's a good idea to take advantage of this,
% together with the built-in elementwise operators.
% Observe that we want a vector X=[X(0) X(1) ... X(K-1)]'=[X(k), k=0,...,K-1]'
% and that each element of the vector is computed as before:
% X(k) = 1/sqrt(N) * sum_{n=0 to N-1} x(n) exp(-j 2 pi k n / N)
% But, as mentioned in (eq. 2, lab. 2) this is nothing more than the inner
% product of x=[x(0) ... x(N-1)]'=[x(n), n=0,...,N-1]' with the complex
% exponential of frequency K and length N ekN=[ekN(0) ... ekN(N-1)]'=[ekN(n),
% n=0,...,N-1]', so that X(0)=<x,e0N>=e0N'*x, X(1)=<x,e1N>=e1N*x', ...,
% X(K-1)=<x,e(K-1)N>=e(K-1)N'*x. So we see that we have a bunch of vectors
% {ekN, k=0,...,K-1} that always multiply the same vector x. We can achieve
% this by creating a matrix where each of the ekN is a row, and then get
% the full vector X by multiplying this matrix by x. We will denote the
% matrix as WKN for future reference:
% WKN=[e0N'; e1N'; ... ; e(K-1)N'];
% Finally, observe that as we move from row to row, k grows from 0 to K-1.
% And as we move from column to column, it is n that grows from 0 to N-1.
% Also, observe that ekN is a function of (k,n) given by exp(-j 2 pi k n /K).
% So if we create a matrices of indices (k,n) then we can directly apply
% the exponential by making use of the elementwise nature of this operation
% in Matlab.

%!!!!!! Uncomment the following for Method 2 to work:
% (Rember to comment other methods)

% matrix_k=repmat((0:K-1)',1,N);
% matrix_n=repmat((0:N-1),K,1);
% indices=matrix_k.*matrix_n;
% WKN=1/sqrt(N)*exp(-1j*2*pi*indices/K);
% X=WKN*x;

% \\\\\ METHOD 3: Built-in fft() function
% Even though the matrix form is fast, it is still not fast enough for large
% signals x. For that, it is better to use the built in fft() function which is
% the optimal way to compute a dft. Besides, it is really easy to code:

X=fft(x,K)/sqrt(N);

% \\\\\ CENTERED FFT.
X_c=circshift(X,ceil(K/2)-1); % Circularly shift X to get it centered in f_c==0

%%\\\\\ END: COMPUTATION OF DFT.
end