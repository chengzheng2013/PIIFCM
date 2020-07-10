function [ result ] = t( img )
    [~,sigma] = possibility(img,3e-6,10);
    if (sigma(1,1)>sigma(2,2))
        sigma_Max = round(sigma(1,1).^(0.5));
        sigma_Min = round(sigma(2,2).^(0.5));
    else
        sigma_Max = round(sigma(2,2).^(0.5));
        sigma_Min = round(sigma(1,1).^(0.5));
    end
    SE=strel('disk',round(2*sigma_Min));
    result=imtophat(img,SE)/255;
end

