%% 1
clc
% Recording voice
norecording = false;
if norecording
    T = 1.5;
    fs = 8192;
    samples = 10;
    one = zeros(samples, T*fs);
    two = zeros(samples, T*fs);
    for i = 1:10
        one(i,:) = record_sound(T,fs);
    end
    disp('Starting two');
    pause(2);
    for i = 1:10
        two(i,:) = record_sound(T);
    end
else
    if ~(exist('one') && exist('two'))
        load('recordings.mat');
    end
end
    

%% 2.1

Fone = zeros(size(Fone));
Ftwo = zeros(size(Ftwo));
[samples,~] = size(Fone);
fs = 8192;

for i = 1:samples
    Fone(i,:) = dft(one(i,:),fs);
    Ftwo(i,:) = dft(two(i,:),fs);
    Fone(i,:) = Fone(i,:)/norm(Fone(i,:));
    Ftwo(i,:) = Ftwo(i,:)/norm(Ftwo(i,:));
end

Foneavg = mean(Fone);
Ftwoavg = mean(Ftwo);


%% 2.2

T = 1.5;
fs = 8192;
x = record_sound(T,fs)';
X = dft(x,fs)';
diffone = innerprod(abs(X),abs(Foneavg));
difftwo = innerprod(abs(X),abs(Ftwoavg));
fprintf('Average Correlation to ''one'': %0.2f\n',diffone);
fprintf('Average Correlation to ''two'': %0.2f\n',difftwo);
if diffone > difftwo
    word = 'one';
else
    word = 'two';
end
fprintf('Average Prediction: %s\n',word);


%% 3

diffoneall = zeros(samples,1);
difftwoall = zeros(samples,1);
for i = 1:samples
    diffoneall(i,1) =  innerprod(abs(X),abs(Fone(i,:)));
    difftwoall(i,1) =  innerprod(abs(X),abs(Ftwo(i,:)));
end
diffone = max(diffoneall);
difftwo = max(difftwoall);
if diffone > difftwo
    word = 'one';
else
    word = 'two';
end
fprintf('Max Correlation to ''one'': %0.2f\n',diffone);
fprintf('Max Correlation to ''two'': %0.2f\n',difftwo);
fprintf('Maximum Prediction: %s\n',word);

%% 4

% Fnums(digitnumber,data,samplenum)
% i: sample number
% k: digit number
norecording = false;
if norecording
    nums = zeros(samples,length(Fone),10);
    nums(:,:,1) = one;
    nums(:,:,2) = two;
    for k = 3:10
        for i = 1:samples
            disp(sprintf('starting recording %d',k));
            Fnums(i,:,k) = record_sound(T,fs);
        end
    end
end

nodft = false;
if nodft
    for k = 1:10
        for i = 1:samples
            Fnums(i,:,k) = dft(nums(i,:,k),fs);
            Fnums(i,:,k) = Fnums(i,:,k)/norm(Fnums(i,:,k));
        end
    end
end

newsample = false;
if newsample
    x = record_sound(T,fs)';
    X = dft(x,fs)';
end

diffnums = zeros(10,samples);
for k = 1:10
    for i = 1:samples
        diffnums(i,k) = innerprod(abs(X),abs(Fnums(i,:,k)));
    end
end
correlations = max(diffnums);
for k=1:10
    fprintf('Correlation to %d: %0.2f\n',k,correlations(k));
end
number = find(correlations==max(correlations));
number = mod(number,10);
fprintf('You Said: %d\n',number);