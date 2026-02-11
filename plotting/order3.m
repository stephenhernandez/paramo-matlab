function [v,d]=order3(v1,d1)
%
% [v,d]=order(v1,d1)
%
% Order eigenvalues, d, into descending order
% Arrange eigenvectors accordingly
%

%%
%v = v1;
% d = d1;
[d,sI] = sort(diag(d1),'descend'); %sort([d1(1,1) d1(2,2) d1(3,3)]);

%%
%ld = length(d2);
v = v1(:,sI);

%%
% for k = 1:3
%     v(1,3-k+1)=v1(1,index(k));
%     v(2,3-k+1)=v1(2,index(k));
%     v(3,3-k+1)=v1(3,index(k));
%     d(3-k+1)=d2(k);
% end