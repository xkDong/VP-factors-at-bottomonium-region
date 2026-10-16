function newenergy=Newenergy(energy,k)
o=length(energy);
newenergy=zeros(1,o*k-k+1);
newenergy(o*k-k+1)=energy(o);
for i=1:o-1
    dl=(energy(i+1)-energy(i))/k;
    newenergy(k*i-k+1)=energy(i);
    for j=1:k-1
        newenergy(k*i-k+j+1)=energy(i)+dl*j;
    end    
end

end