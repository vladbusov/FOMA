% листинг файла Elevel.m
function [EL,Psi]=Elevel(gamma,dE,V,Emax,Xmin,Xmax,Ngreed)
% функция, возвращающая собственные значения и собственные функции 
% уравнения Шредингера
% gamma ? коэффициент, входящий в обезразмеренное уравнение 
% Шредингера (22)
% dE ? приращение энергии
% V ? вектор, содержащий значения потенциала в узлах координатной 
% сетки
% Emax ? максимальное значение энергии
% Xmin ? левая граница координатной сетки
% Xmax ? правая граница координатной сетки
% Ngreed ? число узлов координатной сетки
Tolf=10^-6;
Emin=-1;
E=Emin+dE;
m=1;
Start=0;
while E<Emax
  if Start==0
    E1=E;
    [f,psi]=Num(gamma,E,V,Xmin,Xmax,Ngreed);
    E=E+dE;
    F1=f;
    Start=1;
  end;
  E2=E;
  [f,psi]=Num(gamma,E,V,Xmin,Xmax,Ngreed);
  F2=f;
  if F1*F2>0
    E1=E;
    F1=F2;
    E=E+dE;
  end;
  if F1*F2<0
% уточнение собственного значения энергии  и вычисление значений 
% собственной функции в узлах координатной сетки
    a=(E2-E1)/(F2-F1);  
    E=E1-a*F1;
    [f,psi]=Num(gamma,E,V,Xmin,Xmax,Ngreed);
    F1=F2;
    E1=E2;
    EL(m,1)=m;
    EL(m,2)=E;
    if m==1
      Psi0=psi';
    else
      Psi0=cat(2,Psi0,psi');
    end;
    m=m+1;
    E=E+dE;
  end;
end;
dx=(Xmax-Xmin)/(Ngreed-1);
N1=size(Psi0,2);
% нормировка собственных волновых функций
for i=1:N1
  S=Psi0(:,i);
  S1=S.^2;
  Norm=0;
  for j=1:Ngreed-1
    Norm=Norm+0.5*(S1(j)+S1(j+1))*dx;    
  end;
  S=S./(Norm^0.5);
  if i==1
    Psi=S;
  else
    Psi=cat(2,Psi,S);
  end;
end;
