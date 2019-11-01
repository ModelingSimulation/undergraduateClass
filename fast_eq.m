function [Nbar_a,Nbar_b,Nbar_ea,Nbar_eb] = fast_eq(N_a,N_b,N_e,ra,rb)
%to achieve fast equilibrium
p=zeros(N_a+1,N_b+1);
p(1,1)=1;

for i=1:N_a
    p(i+1,1)=p(i,1)*(N_a-i+1)*(N_e-i+1)/(ra*i);
end

i=1:N_a+1;

for j=1:N_b
    p(i,j+1)=(N_b-j+1)*p(i,j).*(N_e-i'-j+2)/(rb*j);
end
psum=sum(sum(p));
p=p/psum;

% row_max=min(N_e,N_a);
% col_max=min(N_e,N_b);
% p(row_max:end,col_max:end)=0;
Nbar_ea=sum((0:N_a).*sum(p,2)');
Nbar_eb=sum((0:N_b).*sum(p,1));
Nbar_a=N_a-Nbar_ea;
Nbar_b=N_b-Nbar_eb;

end


