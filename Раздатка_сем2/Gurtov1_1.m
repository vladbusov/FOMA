mz_eSi=1.08;
mz_pSi=0.56;
mz_eGe=0.56
mz_pGe=0.35
k=1.38*1e-23;
T1=77;
T2=300;
h=6.63*1e-34;
m0=9.1*1e-31;
%h=1.05*1e-34
Eg0Si=1.21; %��
Eg0Ge=0.8; %��
aSi=2.4*1e-4; %��/�
aGe=5.8*1e-4;
NcSi77=(2*(2*pi*mz_eSi*m0.*k.*T1./h.^2)^(3/2))/1000000 %����� �� ������� ����� �� ������ ������� � ����������
NvSi77=(2*(2*pi*mz_pSi*m0.*k.*T1./h.^2)^(3/2))/1000000
EgSi77=(Eg0Si-aSi*T1)
nSi77=(sqrt(NcSi77*NvSi77))*2.71^(-EgSi77*1.6*1e-19/2/k/T1)
%n1=sqrt(Nc*1e-6*Nv*1e-6) %� ������� ����
%pok=EgSi/2/k/T
%n2=2.71^(-pok)
%nn=n1*n2
%pok2=-(1.2/2/k/T)
%n3=2.71^pok2
%m=n1*n3
NcSi300=(2*(2*pi*mz_eSi*m0.*k.*T2./h.^2)^(3/2))/1000000 %����� �� ������� ����� �� ������ ������� � ����������
NvSi300=(2*(2*pi*mz_pSi*m0.*k.*T2./h.^2)^(3/2))/1000000
EgSi300=(Eg0Si-aSi*T2)
nSi300=(sqrt(NcSi300*NvSi300))*2.71^(-EgSi300*1.6*1e-19/2/k/T2)
%��� �������� ��� �=77�
NcGe77=(2*(2*pi*mz_eGe*m0.*k.*T1./h.^2)^(3/2))/1000000 %����� �� ������� ����� �� ������ ������� � ����������
NvGe77=(2*(2*pi*mz_pGe*m0.*k.*T1./h.^2)^(3/2))/1000000
EgGe77=(Eg0Ge-aGe*T1)
nGe77=(sqrt(NcGe77*NvGe77))*2.71^(-EgGe77*1.6*1e-19/2/k/T1)
%
%������ ������ �� �
T=70:10:300;
NcSi_T=(2*(2*pi*mz_eSi*m0.*k.*T./h.^2).^(3/2))/1000000 %����� �� ������� ����� �� ������ ������� � ����������
NvSi_T=(2*(2*pi*mz_pSi*m0.*k.*T./h.^2).^(3/2))/1000000
EgSi_T=(Eg0Si-aSi.*T)
nSi_T=(sqrt(NcSi_T.*NvSi_T)).*2.71.^(-EgSi_T.*1.6*1e-19/2/k/T)
plot(T,nSi_T);
grid on 