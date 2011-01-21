undef("setdepth")
procedure setdepth(zeta, h, theta_s, theta_b, hc, N,z_r,z_w,Hz)

begin

;  must allocate space for local arrays

sc_w = new(N, double)
Cs_w = new(N, double)
sc_r = new(N, double)
Cs_r = new(N, double)


cff1=1./sinh(theta_s)
cff2=0.5/tanh(0.5*theta_s)
sc_w0=-1.0
Cs_w0=-1.0
ds = 1.0/N

do k=0,N-1    
    ; S-coordinate stretching curves at RHO-points (C_r) and  at W-points (C_w)
    ; S-coordinate at RHO-points (sc_r) and at W-points (sc_w)
    sc_w(k)=ds*(k-N+1)
    Cs_w(k)=(1.-theta_b)*cff1*sinh(theta_s*sc_w(k)) +theta_b*(cff2*tanh(theta_s*(sc_w(k)+0.5))-0.5)
    
    sc_r(k)=ds*((k-N+1)-0.5);
    Cs_r(k)=(1.-theta_b)*cff1*sinh(theta_s*sc_r(k))  +theta_b*(cff2*tanh(theta_s*(sc_r(k)+0.5))-0.5)
end do

z_w(:,:,0)=-h
hinv=1./h
z_r(:,:,0)=-h

do k=0,N-1
    cff_w=hc*(sc_w(k)-Cs_w(k));
    cff1_w=Cs_w(k);
    cff2_w=sc_w(k)+1.;
    
    cff_r=hc*(sc_r(k)-Cs_r(k));
    cff1_r=Cs_r(k);
    cff2_r=sc_r(k)+1.
    ; Depth of sigma coordinate at W-points
    z_w0=cff_w+cff1_w*h
    z_w(:,:,k+1)=z_w0+zeta(:,:)*(1.+z_w0*hinv)
    
    ; Depth of sigma coordinate at RHO-points
    z_r0=cff_r+cff1_r*h;
    z_r(:,:,k)=z_r0+zeta(:,:)*(1.+z_r0*hinv)
end do

z_w(:,:,N)=zeta(:,:)

;
; Hz
;
kk=ispan(1,N,1);  ;Note:  can't use "k" again, unless deleted!!!!
Hz=(z_w(:,:,kk)-z_w(:,:,kk-1));


;if opt ==1
;    Hz=ods*(z_w(:,:,k,t)-z_w(:,:,k-1,t));
;end

end
