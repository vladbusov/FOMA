% ������� ����� Num.m
function [f,psi]=Num(gamma,E,V,Xmin,Xmax,Ngreed)
% �������, ��������� ��� �������� ������� �������� f, �����������
% � ������������ � (18), � �������� �������� ������� � ����� 
% ������������ �����
% gamma ? �����������, �������� � ��������������� ��������� 
% ���������� (22)
% E ? �������� �������� �������
% V ? ������, ���������� �������� ���������� � ����� ������������ 
% �����
% Xmin ? ����� ������� ������������ �����
% Xmax ? ������ ������� ������������ �����
% Ngreed ? ����� ����� ������������ �����
dx=(Xmax-Xmin)/Ngreed;
c=dx.^2*gamma^2/12;
Imatch=1;
psi(1)=0;
psi(2)=9.99999*10^-10;
% �������������� ��������� ���������� ������ ������
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
% �������������� ��������� ���������� ����� �������
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
% ���������� �������� ����� ��������� �������, ����������� 
% ��������������� ��������� ���������� ����� ������� � ������ ������, 
% � ���� � ������� Imatch
f=(psi_Left+psi(Imatch)-(psi(Imatch-1)+psi(Imatch+1)))/Psimax;

