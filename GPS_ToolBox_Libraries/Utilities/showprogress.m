function showprogress(j,varargin)
%Show loop progress if varargin=1
if isempty(varargin) 
    PROG=0;
else
    PROG=varargin{1};
end
if PROG==1;
    if j==1
        fprintf('0');
    end
    if mod(j,10)==0
        n=0;
        while j+1-10*n>0
            n=n+1;
        end
        if mod(n,10)==1
            fprintf('\n')
        end
        fprintf('%i',mod(n-1,10));
    else
        fprintf('.')
    end
end