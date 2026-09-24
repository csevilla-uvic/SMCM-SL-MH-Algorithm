function  [r01,r02,r12,r23,r10,r20,r30]=  equilibriumdistributionSL(Timescales, C,Cl,D,Cn,Wn)
%CS Remove Global, set Timescales manually
tau01 = Timescales(1);
tau02 = Timescales(2);
tau12 = Timescales(3);
tau23 = Timescales(4);
tau10 = Timescales(5);
tau20 = Timescales(6);
tau30 = Timescales(7);

f = 0.5*(gammabb(Cn)+gammabb(Wn));
%the F function used to make the formulas a bit simpler.

r01 = gammabb(Cl)*gammabb(D)*(1.d0-f)/tau01;
r02 = gammabb(C)*(1.d0-gammabb(D))*(1.d0-f)/tau02;
r12 = gammabb(C)*(1.d0-gammabb(D))*(1.d0-f)/tau12;
r23 = 1.0/tau23;

r10 = 1.0/tau10;
r20 = 1.0/tau20;
r30 = 1.0/tau30;

%r01 = gammabb(Cl)*gammabb(D)/tau01;
%r02 = gammabb(C)*(1.d0-gammabb(D))/tau02;
%r10 = gammabb(D)/tau10;
%r12 = gammabb(C)*(1.d0-gammabb(D))/tau12;
%r20 = (1-gammabb(C))/tau20;
%r30 = 1.0/tau30;
%IG update rates
 
    
%IG not sure what's going on in last part here. Think about logic. This is done here to get pi_eq only. leaving for now but definitely should be changed.
%x(1) = 1;
%if(r01==0)
%    x(2) = 0;
%else
%    x(2) = r01/(r10+r12);
%end
%     
%if(r12*r01==0)
%    x(3) = 1/(r20+r23)*(r02);
%else
%    x(3) = 1/(r20+r23)*(r02+(r12*r01)/(r10+r12));
%end
%    
%x(4)=x(3)*r23/r30;
%
%pi_eq=x/sum(x); 
%    
