function w3 = caltitweight(new_alpha1,new_alpha2,M)
fun3 = @(w) (new_alpha1)*(-(igamma(M + 1/2, w^2) - w*igamma(M, w^2))/gamma(M))+new_alpha2*w;
w=[0,10];
w3=fzero(fun3,w);
end

 