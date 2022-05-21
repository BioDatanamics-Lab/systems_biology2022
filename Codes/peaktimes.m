function [peakt,cnt] = peaktimes(v,t,tmin,tmax)

cnt = 0;

peakt=zeros(1,1);
for j=1+1:length(t)-1
    if t(j) >= tmin && t(j) <= tmax
        if v(j-1)<v(j) && v(j+1)<v(j)
            cnt = cnt+1;
            peakt(cnt)=t(j);
            aux = 1;
        end
    end
end