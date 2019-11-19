function image=imresize4D(image0,dim)
    if size(image0,1)~=dim(1) || size(image0,2)~=dim(2) || size(image0,3)~=dim(3)
        image=zeros([dim size(image0,4)]);
        for i4=1:size(image0,4)
            image(:,:,:,i4)=imresize3D(image0(:,:,:,i4),dim);
        end
    else
        image=image0;
    end
end