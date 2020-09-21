function y=quantili(x,w,q)
% computes quantiles q of x with weights w
% w is discrete prob function
% need not be sorted or normalized

[xs,ix]=sort(x);
ws=w(ix);
ws=ws/sum(ws);
cums=cumsum(ws);
y=zeros(length(q),1);
for i=1:length(q)
	y(i)=interp1q(cums,xs,q(i));
end