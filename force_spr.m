function [f] = force_spr(X,links,l_rest,k)
%   Detailed explanation goes here
num_points=size(X,1);
num_links=size(links,1);

f=zeros(num_points,2);
dx=X(links(:,1),:)-X(links(:,2),:);
l_links=sqrt(sum(dx.^2,2));
l_spr=l_links-l_rest;
dir=dx./l_links;
ff=-k*l_spr.*dir;

for i=1:num_links
    f(links(i,1),:)=f(links(i,1),:)+ff(i,:);
    f(links(i,2),:)=f(links(i,2),:)-ff(i,:);
end

end

