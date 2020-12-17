function [f]=Bird_Planform(ydb)

% K-O planform
if abs(ydb)<=0.25
    f=1;
elseif abs(ydb)>0.25
    f=8*abs(ydb)*(1-2*abs(ydb));
end











