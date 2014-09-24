function z=ZBfilter(f,fcenter,fwidth);

z= exp(-(f-fcenter).^2/(2*fwidth^2));
