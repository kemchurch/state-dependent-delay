function Sval = eval1spline(s,t)
global pushInf
OutOfDomain = 0;
M = get_segments(s);
a = get_a(s);
b = get_b(s);
dom = get_dom(s);
pieces = get_pieces(s);
y = [b;b(end,:)+a(end,:)/M];
x = dom(1,:) + (dom(2,:)-dom(1,:)).*((0:M).'/M);
for k=1:pieces-1
   x(end,k) = NaN;
   y(end,k) = NaN;
end
x = reshape(x,[(M+1)*pieces,1]); x(isnan(x))=[];
y = reshape(y,[(M+1)*pieces,1]); y(isnan(y))=[];
% Process t and x to avoid weird zero conflicts
x(abs(x)<1E-15) = 0;
t(abs(t)<1E-15) = 0;
Sval = interp1(x,y,t);
% fix NaN indices
NaNInd = isnan(Sval);
if sum(NaNInd)>0
    NaNInd_left = NaNInd.*(t<=dom(1,1)).*((1:length(t))');
    NaNInd_right = NaNInd.*(t>=dom(2,end)).*((1:length(t))');
    OutOfDomain = 1;
    if sum(NaNInd_left)>0
        NaNInd_left(NaNInd_left==0)=[];
        if NaNInd_left(end)+2>length(t)
           Sval(NaNInd_left) = b(1,1);
        else
            diff = (Sval(NaNInd_left(end)+2)-Sval(NaNInd_left(end)+1))*M/(dom(2,1)-dom(1,1));
            Sval(NaNInd_left) = Sval(NaNInd_left(end)+1) + diff*(t(NaNInd_left)-t(NaNInd_left(end)+1));
        end
    end
    if sum(NaNInd_right)>0
        NaNInd_right(NaNInd_right==0)=[];
        if NaNInd_right(1)-2<1
           Sval(NaNInd_right) = b(end,end);
        else
            diff = (Sval(NaNInd_right(1)-2)-Sval(NaNInd_right(1)-1))*M/(dom(2,end)-dom(1,end));
            Sval(NaNInd_right) = Sval(NaNInd_right(1)-1) + diff*(t(NaNInd_right)-t(NaNInd_right(1)-1));
        end
    end
end
if sum(isnan(Sval))>0
    error('There is a NaN.');
end
if OutOfDomain == 1 & pushInf == 1
    Sval = Inf*Sval;
end
end