function [XYZ]=lorenz_RK(r, s, b,  total_time, tstep)

vec=randn(3,1);
f=zeros(3,1);
nsteps = total_time/tstep;
store_vec=zeros(3,nsteps);


for n=1:nsteps;
    
    h2=tstep/2;
    
    f(1)=s*(vec(2)-vec(1));
    f(2)=vec(1)*(r-vec(3))-vec(2);
    f(3)=vec(1)*vec(2)-b*vec(3);
    k1=tstep.* f;
    
    vec2=vec+k1/2;
    f(1)=s*(vec2(2)-vec2(1));
    f(2)=vec2(1)*(r-vec2(3))-vec2(2);
    f(3)=vec2(1)*vec2(2)-b*vec2(3);
    k2=tstep.* f;
    
    vec3=vec+k2/2;
    f(1)=s*(vec3(2)-vec3(1));
    f(2)=vec3(1)*(r-vec3(3))-vec3(2);
    f(3)=vec3(1)*vec3(2)-b*vec3(3);
    k3=tstep.* f;
    
    vec4=vec+k3;
    f(1)=s*(vec4(2)-vec4(1));
    f(2)=vec4(1)*(r-vec4(3))-vec4(2);
    f(3)=vec4(1)*vec4(2)-b*vec4(3);
    k4=tstep.*f;
    
    vec_nplus1=vec+(k1/6)+(k2/3)+(k3/3)+(k4/6);
    store_vec(:,n)=vec;
    vec=vec_nplus1;
    
end

XYZ=store_vec';

