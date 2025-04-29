function K = elemKVEM(d,e,rmt,coor)
Nv = size(coor,1);

C = [d(1,1),d(1,2),d(1,3),e(1,1),e(1,2);
    d(2,1),d(2,2),d(2,3),e(2,1),e(2,2);
    d(3,1),d(3,2),d(3,3),e(3,1),e(3,2);
    e(1,1),e(2,1),e(3,1),rmt(1,1),rmt(1,2);
    e(1,2),e(2,2),e(3,2),rmt(2,1),rmt(2,2);];

[dM,Pis,Pi,area] = calculatePi(coor);

Pis1 = blkdiag(Pis,Pis,Pis);

A1 = [1,0;0,0;0,1;0,0;0,0];
A2 = [0,0;0,1;1,0;0,0;0,0];
A3 = [0,0;0,0;0,0;1,0;0,1];
A = [A1,A2,A3];

dM1 = blkdiag(dM',dM',dM');

G0 = area*dM1'*A'*C*A*dM1;

Kc = Pis1'*G0*Pis1;

Kcuu = Kc(1:Nv*2,1:Nv*2);
Kcdd = Kc(Nv*2+1:end,Nv*2+1:end);
Piuu = blkdiag(Pi,Pi);
IIuu = eye(Nv*2);
IIdd = eye(Nv);
Pidd = Pi;

alpha1 = 1/Nv*trace(Kcuu);
alpha2 = 1/Nv*trace(Kcdd);

alpha3 = 1/Nv*max(max(e))*0.1;
kes = [alpha3 * ((IIdd - Pidd)' * (IIdd - Pidd));alpha3 * ((IIdd - Pidd)' * (IIdd - Pidd))];

Ks = [alpha1 * ((IIuu - Piuu)' * (IIuu - Piuu)),kes;kes',alpha2 * ((IIdd - Pidd)' * (IIdd - Pidd))];

AK = Kc+Ks;
K = AK;
