function display_image(img)
    
    img = norm_image(img,1);
    img = img.^0.4545;
    img = uint8(norm_image(img,255));
    imshow(img);
end