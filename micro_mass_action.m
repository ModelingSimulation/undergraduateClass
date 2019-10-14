%10/11
clear

num_a=100; %set number of different kinds of molecules
num_b=200;
num_c=100;
num_ab=1;
num_abc=10;
r1=1; %set different probability per unit time for 
r2=0.5; %different reactions, r1,2,3,4.
r3=0.3;
r4=0.1;
Tmax=20; %maximum time for all the reaction
t=0; %time
result_list=zeros(1,6);
i=0; %number of events

%%An event-driven simulation
while t<Tmax
%     if num_abc>50 %An alterlate option for a catalyst
%         r4=0.2;
%     end
    
%     if t>4.5  %An alterlate option for a catalyst
%         r2=1;
%     end
    i=i+1; 
    t1=-(log(rand))/(r1*num_a*num_b); %creating exponential distributed
    t2=-(log(rand))/(r2*num_ab); %random variables by different  
    t3=-(log(rand))/(r3*num_a*num_b*num_c); %probability per unit time
    t4=-(log(rand))/(r4*num_abc);
    tlist=[t1 t2 t3 t4];
    
    [t_event i_event]=min(tlist); %choose the first event that happens to happen
    
    if i_event ==1
        %first reaction: a+b->ab
        num_a=num_a-1;
        num_b=num_b-1;
        num_ab=num_ab+1;
        t=t+t1;
    elseif i_event==2     
        %second reaction: ab->a+b
        num_a=num_a+1;
        num_b=num_b+1;
        num_ab=num_ab-1;
        t=t+t2;
    elseif i_event==3 
        %third reaction: a+b+c->abc
        num_a=num_a-1;
        num_b=num_b-1;
        num_c=num_c-1;
        num_abc=num_abc+1;
        t=t+t3;  
    else 
        %fourth reaction: abc->ab+c
        num_abc=num_abc-1;
        num_ab=num_ab+1;
        num_c=num_c+1;
        t=t+t4;
    end
    result_list(i,:)=[t num_a num_b num_c num_ab num_abc];
end
%%plot results
for j =1:5
    plot(result_list(:,1),result_list(:,j+1));
    hold on
end
