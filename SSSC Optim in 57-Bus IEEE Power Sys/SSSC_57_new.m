function object=SSSC_57_new(x);

 x=[0.0	0.078125	0.02734375	0.0	0.03588280512247244	0.015625	0.03515625	0.0625	0.0625	0.109375	0.0546875	0.015625	0.06640625	0.0	-0.06640625	0.171875	0.01171875	0.01171875	0.31606605803402177	0.046875	0.046875	0.0	0.0078125	0.3702492614929396	0.03515625	0.0625	0.22265625	0.0	0.015625	0.0625	0.1484375	0.0	0.0	0.0	0.3828526226439568	0.0	0.4375	0.07421875	0.28125	0.0	0.078125	0.835576173704977	0.25	0.1328125	0.140625	0.09375	0.0	1.0	0.03515625	0.0	-0.9375	0.6625175089707293	0.0	0.0390625	0.0	0.0625	0.0	0.0	0.0390625	0.0	0.109375	0.0	0.046875	0.0	1.0	0.0	0.08203125	0.0	0.8564359266908602	0.0625	0.0	0.0	0.3629916366985626	0.3828125	0.01171875	0.0	0.12109375	0.625	0.015625	0.0625];

 global line voltage qwe voltagee fazz faz

%-----------------   Network Information   ---------------------------------------------------------

for n=1:57
    voltagee(n,1)=voltage(n,1);
end

admittance=zeros(57,57);

for i=1:80
    admittance(line(i,1),line(i,2))=(-1/line(i,3))*182.25;
    admittance(line(i,2),line(i,1))=(-1/line(i,3))*182.25;
end

for i=1:57
   for m=1:57 
       admittance(i,i)=admittance(i,i)-(admittance(i,m));
   end
end

dd=inv(admittance);
xt=.1*j;


%-----------------   This box makes inputs discrete   ---------------------------------------------------------

Vm=1;

for z=1:80
    if x(1,z)<.1*Vm && x(1,z)>=0*Vm
        x(1,z)=0*Vm;
    end
    if x(1,z)<0.2*Vm && x(1,z)>=0.1*Vm
        x(1,z)=.1*Vm;
    end
    if x(1,z)<0.3*Vm && x(1,z)>=0.2*Vm
        x(1,z)=0.2*Vm;
    end
    if x(1,z)<.4*Vm && x(1,z)>=0.3*Vm
        x(1,z)=0.3*Vm;
    end
    if x(1,z)<.5*Vm && x(1,z)>=.4*Vm
        x(1,z)=0.4*Vm;
    end
    if x(1,z)<.6*Vm && x(1,z)>=0.5*Vm
        x(1,z)=0.5*Vm;
    end
    if    x(1,z)<.7*Vm && x(1,z)>=.6*Vm
        x(1,z)=.6*Vm;
    end
    if     x(1,z)<.8*Vm && x(1,z)>=0.7*Vm
        x(1,z)=0.7*Vm;
    end
    if    x(1,z)<.9*Vm && x(1,z)>=0.8*Vm
        x(1,z)=0.8*Vm;
    end
    if    x(1,z)<1*Vm && x(1,z)>=0.9*Vm
        x(1,z)=0.9*Vm;
    end
    if   x(1,z)>=1*Vm
        x(1,z)=1*Vm;
    end
    
    
    
    if x(1,z)<0*Vm && x(1,z)>=-0.1*Vm
        x(1,z)=0*Vm;
    end
    if x(1,z)<-.1*Vm && x(1,z)>=-0.2*Vm
        x(1,z)=-0.1*Vm;
    end
    if x(1,z)<-.2*Vm && x(1,z)>=-0.3*Vm
        x(1,z)=-0.2*Vm;
    end
    if x(1,z)<-.3*Vm && x(1,z)>=-0.4*Vm
        x(1,z)=-0.3*Vm;
    end
    if x(1,z)<-.4*Vm && x(1,z)>=-.5*Vm
        x(1,z)=-0.4*Vm;
    end
    if x(1,z)<-.5*Vm && x(1,z)>=-0.6*Vm
        x(1,z)=-0.5*Vm;
    end
    if    x(1,z)<-.6*Vm && x(1,z)>=-.7*Vm
        x(1,z)=-.6*Vm;
    end
    if     x(1,z)<-.7*Vm && x(1,z)>=-0.8*Vm
        x(1,z)=-0.7*Vm;
    end
    if    x(1,z)<-.8*Vm && x(1,z)>=-0.9*Vm
        x(1,z)=-0.8*Vm;
    end
    if    x(1,z)<-.9*Vm && x(1,z)>=-1*Vm
        x(1,z)=-0.9*Vm;
    end
    if    x(1,z)<-1*Vm 
        x(1,z)=-1*Vm;
    end

end

%----------------------------       Sens. Analysis    ----------------------------------------------

comp=find(x~=0);
numb=length(comp);
deltay=zeros(57,57);

if numb>40
    objectt=100;
else
    powerM=0;
    VV=0;
    phase=0;

    if numb==0
        objectt=0;
    else
        for k=1:100
            deltay=zeros(57,57);
            for i=1:numb
                current(i,k)=(voltagee(line(comp(i),1),k)-voltagee(line(comp(i),2),k)-x(1,comp(i)))/(xt+line(comp(i),3));
                z(k,comp(i))=x(1,comp(i))/current(i,k);
                qwe(1,i)= (x(1,comp(i))/current(i,k))- line(i,3);
                        
                deltay(line(comp(i),1),line(comp(i),1))=1/(xt+z(k,comp(i))+line(comp(i),3))-1/(line(comp(i),3));
                deltay(line(comp(i),2),line(comp(i),2))=deltay(line(comp(i),1),line(comp(i),1));
                deltay(line(comp(i),1),line(comp(i),2))=-deltay(line(comp(i),1),line(comp(i),1));
                deltay(line(comp(i),2),line(comp(i),1))=-deltay(line(comp(i),1),line(comp(i),1));
        
              
                deltavoltage(:,k)=-dd*deltay*voltagee(:,1);   %------ Eq. (15)---------------
                voltagee(:,k+1)= voltagee(:,1)+ deltavoltage(:,k);
            end
        end
    
%---------------------------------   O.F.   -------------------------------------------    
    
        for p=1:80
                faz(p,1)=abs(angle(voltagee(line(p,1),1))-angle(voltagee(line(p,2),1)));
                faz(p,2)=abs(angle(voltagee(line(p,1),k))-angle(voltagee(line(p,2),k)));
                faz(p,3)=(faz(p,1)-faz(p,2))/faz(p,1);         %--------to be edited
                if faz(p,1)>pi
                    faz(p,1)=abs(faz(p,1)-2*pi);
                end
                if faz(p,2)>pi
                    faz(p,2)=abs(faz(p,2)-2*pi);
                end
        end
        
        for p=1:80
                power(p,1)=(abs(voltagee(line(p,1),1))*abs(voltagee(line(p,2),1)))/imag(line(p,3));
                power(p,2)=(abs(voltagee(line(p,1),k+1))*abs(voltagee(line(p,2),k+1)))/imag(line(p,3));
                power(p,3)=(power(p,2)-power(p,1))/power(p,1);
        end
        
%-------------------   Constraints   ---------------------------------------------------        
        
        objectt=0;
        for l=1:57
            if abs(voltagee(l,k))>1.2 || abs(voltagee(l,k))<0.8
               objectt=objectt+10000000000;
            end
        end
        for p=1:80
            if faz(p,2)>pi/3
               objectt=objectt+10000000000000;
            else
               fazz(p,1)=faz(p,3);
            end
        end
                              
%---------------------   Ideas for individual weight coefficients   -----------------------------------------------------
       
        aa=0;
        for t=1:80
            aa=aa+abs(sin(faz(t,1)))+abs(sin(faz(t,2)));
        end
        for e=1:80
            powerM=powerM+power(e,3)*((abs(sin(faz(e,1)))+abs(sin(faz(e,2))))/aa);
        end
        
        
        
        bb=0;
        for y=1:57
            bb=bb+abs(abs(voltagee(y,1))-1)+abs(abs(voltagee(y,k+1))-1);
        end
        for q=1:57
            V_0=abs(abs(voltagee(q,1))-1);
            V_comp=abs(abs(voltagee(q,k+1))-1);
            VV=VV+((V_0-V_comp)/V_0)*((V_0+V_comp)/bb);
        end
    
        
        cc=0;
        for u=1:80
            cc=cc+faz(u,1)/pi+faz(u,2)/pi;
        end
        for r=1:80
            phase=phase+(fazz(r,1))*(faz(r,1)+faz(r,2))/(pi*cc);
        end
        
        
        objectt=-(phase+powerM+VV)+objectt;

end
end

object=objectt;

%--------------------------------------------------------------------------
% 
x=x
numb 

V_0=abs(voltagee(:,1)) 
V_angle_0=(angle(voltagee(:,1)))*180/pi
  
V_comp=abs(voltagee(:,k+1))
V_angle_comp=(angle(voltagee(:,k+1)))*180/pi
  
bbb=(abs(abs(voltagee(:,1))-1)+abs(abs(voltagee(:,k+1))-1))/bb


phase_0=faz(:,1)*180/pi
phase_comp=faz(:,2)*180/pi

ccc=(faz(:,1)/pi+faz(:,2)/pi)/cc


power_0=power(:,1)
power_comp=power(:,2)
aaa=(abs(sin(faz(:,1)))+abs(sin(faz(:,2))))/aa


a_0=abs(sin(faz(:,1)))
a_comp=abs(sin(faz(:,2)))

object



%--------------------------------------------------------------------------

% A=0;
% B=0;
% C=0;
% for i=1:6
%     A=A+aaa(i);
%     B=B+bbb(i);
%     C=C+ccc(i);
% end
% A;
% B;
% C;

%--------------------------------------------------------------------------

%end
 