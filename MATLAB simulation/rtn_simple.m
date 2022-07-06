function [ x ] = rtn_simple(t,dt,tc,impact)
% The time constants (tc) are [avg_time_in_high avg_time_in_low]
% The trap's impact is deltaId

tau=((1./tc(:,2))+(1./tc(:,1))).^-1;

trap_state=[round(rand())==1;]; % initialize trap state

for i=2:length(t)
     if trap_state(1,i-1)
        trap_state(1,i)=(rand()>=(1-exp(-dt/tau(1)))*(tc(1,2)/(tc(1,2)+tc(1,1))));
     else
        trap_state(1,i)=(rand()<=(1-exp(-dt/tau(1)))*(tc(1,1)/(tc(1,2)+tc(1,1))));
     end
end

avg=0;
for i=1:length(t)
    if (trap_state(1,i))
        x(i)=impact/2;
    else
        x(i)=-impact/2;
    end
    avg=x(i)+avg;
end
x=x-avg/length(t);

