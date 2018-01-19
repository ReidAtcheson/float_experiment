using Distributions
N=Normal(-0.940,7.38);
function f(x)
  return pdf(N,x);
end
(P1,err)=quadgk(f,20,100);
(P2,err)=quadgk(f,-100,-20);
P=P1+P2;
println(P);


(P3,err)=quadgk(f,-10,10);


println(P3);
