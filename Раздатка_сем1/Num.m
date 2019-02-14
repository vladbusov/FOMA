% листинг файла Num.m
function [f,psi]=Num(gamma,E,V,Xmin,Xmax,Ngreed)
% функци€, возвраща€ дл€ заданной энергии величину f, вычисл€емую
% в соответствие с (18), и значени€ волновой функции в узлах 
% координатной сетки
% gamma ? коэффициент, вход€щий в обезразмеренное уравнение 
% Ўредингера (22)
% E ? заданное значение энергии
% V ? вектор, содержащий значени€ потенциала в узлах координатной 
% сетки
% Xmin ? лева€ граница координатной сетки
% Xmax ? права€ граница координатной сетки
% Ngreed ? число узлов координатной сетки
dx=(Xmax-Xmin)/Ngreed;
c=dx.^2*gamma^2/12;
Imatch=1;
psi(1)=0;
psi(2)=9.99999*10^-10;
% интегрирование уравнени€ Ўредингера справа налево
Kim1=c*(E-V(1));
Ki=c*(E-V(2));
for i=2:Ngreed-1
   Kip1=c*(E-V(i+1));
   if and(Ki*Kip1<=0,Ki>0)
      Imatch=i;
      i=Ngreed;
   end;
   if i<=Ngreed-1
      psi(i+1)=(psi(i)*(2-10*Ki)-psi(i-1)*(1+Kim1))/(1+Kip1);
      if abs(psi(i+1))>=10^10
        for k=1:i+1
           psi(k)=psi(k)*9.99999*10^-6;    
        end;
     end;
     Kim1=Ki;
     Ki=Kip1;
  end;
end;
if Imatch==1
  Imatch=Ngreed-10;
end;
% интегрирование уравнени€ Ўредингера слева направо
psi_Left=psi(Imatch);
psi(Ngreed)=0;
psi(Ngreed-1)=9.99999*10^-10;
Kip1=c*(E-V(Ngreed));
Ki=c*(E-V(Ngreed-1));
for i=Ngreed-1:-1:Imatch+1
   Kim1=c*(E-V(i-1));
   psi(i-1)=(psi(i)*(2-10*Ki)-psi(i+1)*(1+Kip1))/(1+Kim1);
   if abs(psi(i))>10^10
     for k=Ngreed:-1:i
        psi(k)=psi(k)*9.99999*10^-6;
     end;
  end;
  Kip1=Ki;
  Ki=Kim1;
end;
if psi_Left<0
   for i=Imatch:Ngreed
      psi(i)=-psi(i);
   end;
end;
psi1=abs(psi);
Psimax=max(psi1);
% вычисление разности между волновыми функци€, полученными 
% интегрированием уравнени€ Ўредингера слева направо и справа налево, 
% в узле с номером Imatch
f=(psi_Left+psi(Imatch)-(psi(Imatch-1)+psi(Imatch+1)))/Psimax;

