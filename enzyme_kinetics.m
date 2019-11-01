%10/18
clear
for num_e=[20]
    
num_a=1000; %set number of different kinds of molecules
num_b=100;
 %number of enzymes should be much less than a or b
num_ea=1;
num_eb=1;
ra=5; %set different probability per unit time for 
rb=3; %ra&rb

alpha=1;
beta=1;

alphae=10;
betae=1;

Tmax=10; %maximum time for all the reaction
t=0; %time
result_list=zeros(1,6);
i=0; %number of events

%%An event-driven simulation
while t<Tmax

    i=i+1; 
    %%
    [Nbar_a,Nbar_b,Nbar_ea,Nbar_eb]=fast_eq(num_a,num_b,num_e,ra,rb);
    %%
    t1=-(log(rand))/(Nbar_a*alpha); %creating exponential distributed
    t2=-(log(rand))/(Nbar_b*beta); %random variables by different  
    t3=-(log(rand))/(Nbar_ea*alphae); %probability per unit time
    t4=-(log(rand))/(Nbar_eb*betae);
    tlist=[t1 t2 t3 t4];
    
    [t_event i_event]=min(tlist); %choose the first event that happens to happen
 
    if i_event ==1
        %first reaction: a->b
        num_a=num_a-1;
        num_b=num_b+1;        
        t=t+t1;
    elseif i_event==2     
        %second reaction: b->a
        num_a=num_a+1;
        num_b=num_b-1;
        t=t+t2;
    elseif i_event==3 
        %third reaction: ea->eb
        num_a=num_a-1;
        num_b=num_b+1;
        t=t+t3;  
    else 
        %fourth reaction: eb->ea
        num_a=num_a+1;
        num_b=num_b-1;
        t=t+t4;
    end
    result_list(i,1:4)=[t num_a num_b num_e];
    result_list(i,5)=(num_a-num_b)/t;
end
%%plot results
num_event=length(result_list);
% for j =4
%     plot(result_list(ceil(:,1),result_list(:,j+1));
%     hold on
% end
%hold on


end
%legend ('NA','NB','NE')