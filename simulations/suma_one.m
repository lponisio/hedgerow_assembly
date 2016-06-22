function nn_new=suma_one(nn)
         
         [ma im]=max(nn);
         
         nn1=nn;
         nn1(im)=[];
         s1=sum(nn1);
         new=1-s1;
         
         nn_new=nn;
         nn_new(im)=new;
end
         
