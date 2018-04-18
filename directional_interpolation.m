function [mos, maski] = directional_interpolation(rgb, mask)
    maski = zeros(size(rgb));
    maski= maski(:,:,1:3);
    Rh = zeros(size(sum(maski,3)));
    Rv = zeros(size(Rh));
    Bh = zeros(size(Rh));
    Bv = zeros(size(Rh));
    G= rgb(:,:,2);
    G1=zeros(size(G));
    G1(1:2:end,2:2:end)=G(1:2:end,2:2:end);
    G2=zeros(size(G));
    G2(2:2:end,1:2:end)=G(2:2:end,1:2:end);
  
    % bggr color filter array (mask)
    maski(1:2:end,2:2:end,2)=1;
    maski(2:2:end,1:2:end,2)= 1;
    maski(1:2:end,1:2:end,3)=1;
    maski(2:2:end, 2:2:end,1)=1;
    % Bh
    for i=1:2:size(Bh,1)
        for j=1:2:size(Bh,2)
            
            if mask(i,j,3)~=1
                if j>1 && j<(size(Bh,2)-1)
                    Bh(i,j)= (rgb(i,j-2,3)+rgb(i,j+2,3))/2;
                elseif j==1
                    Bh(i,j)= rgb(i,j+2,3);
                elseif j==size(Bh,2) || j==(size(Bh,2)-1)
                    Bh(i,j)= rgb(i,j-2,3);
                end
            else
                Bh(i,j)=rgb(i,j,3);
            end
            
        end
    end
    % Bv
    for j=1:2:size(Bv,2)
        for i=1:2:size(Bv,1)
            
            if mask(i,j,3)~=1
                if i>1 && i<(size(Bv,1)-1)
                    Bv(i,j)= (rgb(i-2,j,3)+rgb(i+2,j,3))/2;
                elseif i==1
                    Bv(i,j)= rgb(i+2,j,3);
                elseif i==size(Bv,1) ||i==(size(Bv,1)-1)
                    Bv(i,j)= rgb(i-2,j,3);
                end
            else
                Bv(i,j)=rgb(i,j,3);
            end
            
        end
    end
    
    %Rh
    for i=2:2:size(Rh,1)
        for j=2:2:size(Rh,2)
            
            if mask(i,j,1)~=1
                if j>2 && j<(size(Rh,2)-1)
                    Rh(i,j)= (rgb(i,j-2,1)+rgb(i,j+2,1))/2;
                elseif j==2
                    Rh(i,j)= rgb(i,j+2,1);
                elseif j==size(Rh,2) || j==(size(Rh,2)-1)
                    Rh(i,j)= rgb(i,j-2,1);
                end
            else
                Rh(i,j)=rgb(i,j,1);
            end
            
        end
    end
    %Rv
    for j=2:2:size(Rv,2)
        for i=2:2:size(Rv,1)
            
            if mask(i,j,1)~=1
                if i>2 && i<(size(Rv,1)-1)
                    Rv(i,j)= (rgb(i-2,j,1)+rgb(i+2,j,1))/2;
                elseif i==2
                    Rv(i,j)= rgb(i+2,j,1);
                elseif i==size(Rv,1) ||i==(size(Bv,1)-1)
                    Rv(i,j)= rgb(i-2,j,1);
                end
            else
                Rv(i,j)=rgb(i,j,1);
            end
            
        end
    end
    
    %%%%%  for RH horizontal estimation *********************
    RGh = zeros(size(Rh));
    for i=2:2:size(Rh,1)
        for j=2:2:size(Rh,2)
            
            if j>2 && j<(size(Rh,2)-1)
                RGh(i,j)= (-1/4)*Rh(i,j-2) + (1/2)*G(i,j-1)+(1/2)*Rh(i,j)+(1/4)*G(i,j+1)+(-1/4)*Rh(i,j+2);
            elseif j==2
                RGh(i,j)= (-1/4)*Rh(i,j+2) + (1/2)*G(i,j-1)+(1/2)*Rh(i,j)+(1/4)*G(i,j+1)+(-1/4)*Rh(i,j+2);
            elseif j==size(Rh,2) || j==(size(Rh,2)-1)
                RGh(i,j)= (-1/4)*Rh(i,j-2) + (1/2)*G(i,j-1)+(1/2)*Rh(i,j)+(1/4)*G(i,j-1)+(-1/4)*Rh(i,j-2);
            end

            
        end
    end
    
    %GRh
    GRh = zeros(size(Rh));
    for i=2:2:size(GRh,1)
        for j=1:2:size(GRh,2)
            
            if j>1 && j<(size(Rh,2)-1)
                GRh(i,j)= (-1/4)*G(i,j-2) + (1/2)*Rh(i,j-1)+(1/2)*G(i,j)+(1/4)*Rh(i,j+1)+(-1/4)*G(i,j+2);
            elseif j==1
                GRh(i,j)= (-1/4)*G(i,j+2) + (1/2)*Rh(i,j+1)+(1/2)*G(i,j)+(1/4)*Rh(i,j+1)+(-1/4)*G(i,j+2);
            elseif j==size(Rh,2) || j==(size(Rh,2)-1)
                GRh(i,j)= (-1/4)*G(i,j-2) + (1/2)*Rh(i,j-1)+(1/2)*G(i,j)+(1/4)*Rh(i,j-1)+(-1/4)*G(i,j-2);
            end

            
        end
    end
    %BGh
    BGh = zeros(size(Bh));
    for i=1:2:size(BGh,1)
        for j=1:2:size(BGh,2)
            
            if j>1 && j<(size(BGh,2)-1)
                BGh(i,j)= (-1/4)*Bh(i,j-2) + (1/2)*G(i,j-1)+(1/2)*Bh(i,j)+(1/4)*G(i,j+1)-(1/4)*Bh(i,j+2);
            elseif j==1
                BGh(i,j)= (-1/4)*Bh(i,j+2) + (1/2)*G(i,j+1)+(1/2)*Bh(i,j)+(1/4)*G(i,j+1)-(1/4)*Bh(i,j+2);
            elseif j==size(BGh,2) || j==(size(BGh,2)-1)
                BGh(i,j)= (-1/4)*Bh(i,j-2) + (1/2)*G(i,j-1)+(1/2)*Bh(i,j)+(1/4)*G(i,j-1)-(1/4)*Bh(i,j-2);
            end

            
        end
    end
    % GBh
    GBh = zeros(size(Bh));
    for i=1:2:size(GBh,1)
        for j=2:2:size(GBh,2)
            
            if j>2 && j<(size(GBh,2)-1)
                GBh(i,j)= (-1/4)*G(i,j-2) + (1/2)*Bh(i,j-1)+(1/2)*G(i,j)+(1/4)*Bh(i,j+1)-(1/4)*G(i,j+2);
            elseif j==2
                GBh(i,j)= (-1/4)*G(i,j+2) + (1/2)*Bh(i,j-1)+(1/2)*G(i,j)+(1/4)*Bh(i,j+1)-(1/4)*G(i,j+2);
            elseif j==size(GBh,2) || j==(size(GBh,2)-1)
                GBh(i,j)= (-1/4)*G(i,j-2) + (1/2)*Bh(i,j-1)+(1/2)*G(i,j)+(1/4)*Bh(i,j-1)-(1/4)*G(i,j-2);
            end

            
        end
    end
    
    
    %%%  vertical estimation *************************
    BGv = zeros(size(Bv));
    for j=1:2:size(BGv,2)
        for i=1:2:size(BGv,1)
            
            if i>1 && i<(size(Rh,1)-1)
                BGv(i,j)= (-1/4)*Bv(i-2,j) + (1/2)*G(i-1,j)+(1/2)*Bv(i,j)+(1/4)*G(i+1,j)+(-1/4)*Bv(i+2,j);
            elseif i==1
                BGv(i,j)= (-1/4)*Bv(i+2,j) + (1/2)*G(i+1,j)+(1/2)*Bv(i,j)+(1/4)*G(i+1,j)+(-1/4)*Bv(i+2,j);
            elseif i==size(Rh,1) || i==(size(Rh,1)-1)
                BGv(i,j)= (-1/4)*Bv(i-2,j) + (1/2)*G(i-1,j)+(1/2)*Bv(i,j)+(1/4)*G(i-1,j)+(-1/4)*Bv(i-2,j);
            end
   
        end
    end
    % GBv
    GBv = zeros(size(Bv));
    for j=1:2:size(GBv,2)
        for i=2:2:size(GBv,1)
            
            if i>2 && i<(size(GBh,1)-1)
                GBv(i,j)= (-1/4)*G(i-2,j) + (1/2)*Bv(i-1,j)+(1/2)*G(i,j)+(1/4)*Bv(i+1,j)+(-1/4)*G(i+2,j);
            elseif i==2
                GBv(i,j)= (-1/4)*G(i+2,j) + (1/2)*Bv(i-1,j)+(1/2)*G(i,j)+(1/4)*Bv(i+1,j)+(-1/4)*G(i+2,j);
            elseif i==size(Rh,1) || i==(size(Rh,1)-1)
                GBv(i,j)= (-1/4)*G(i-2,j) + (1/2)*Bv(i-1,j)+(1/2)*G(i,j)+(1/4)*Bv(i-1,j)+(-1/4)*G(i-2,j);
            end
    
        end
    end
    % GRv
    GRv = zeros(size(Rv));
    for j=2:2:size(GRv,2)
        for i=1:2:size(GRv,1)
            
            if i>1 && i<(size(Rh,1)-1)
                GRv(i,j)= (-1/4)*G(i-2,j) + (1/2)*Rv(i-1,j)+(1/2)*G(i,j)+(1/4)*Rv(i+1,j)+(-1/4)*G(i+2,j);
            elseif i==1
                GRv(i,j)= (-1/4)*G(i+2,j) + (1/2)*Rv(i+1,j)+(1/2)*G(i,j)+(1/4)*Rv(i+1,j)+(-1/4)*G(i+2,j);
            elseif i==size(Rh,1) || i==(size(Rh,1)-1)
                GRv(i,j)= (-1/4)*G(i-2,j) + (1/2)*Rv(i-1,j)+(1/2)*G(i,j)+(1/4)*Rv(i-1,j)+(-1/4)*G(i-2,j);
            end
   
        end
    end
    % RGv
    RGv = zeros(size(Rv));
    for j=2:2:size(RGv,2)
        for i=2:2:size(RGv,1)
            
            if i>2 && i<(size(Rh,1)-1)
                RGv(i,j)= (-1/4)*Rv(i-2,j) + (1/2)*G(i-1,j)+(1/2)*Rv(i,j)+(1/4)*G(i+1,j)+(-1/4)*Rv(i+2,j);
            elseif i==2
                RGv(i,j)= (-1/4)*Rv(i+2,j) + (1/2)*G(i-1,j)+(1/2)*Rv(i,j)+(1/4)*G(i+1,j)+(-1/4)*Rv(i+2,j);
            elseif i==size(Rh,1) || i==(size(Rh,1)-1)
                RGv(i,j)= (-1/4)*Rv(i-2,j) + (1/2)*G(i-1,j)+(1/2)*Rv(i,j)+(1/4)*G(i-1,j)+(-1/4)*Rv(i-2,j);
            end
   
        end
    end
    % color difference **************************
    RGh_dif= RGh-Rh;
    RGh_dif = RGh_dif+(G2-GRh);
    RGv_dif = RGv - Rv;
    RGv_dif = RGv_dif -(G1-GRv);
    
    BGh_dif = BGh-Bh;
    BGh_dif =BGh_dif +(G1-GBh);
    BGv_dif = BGv -Bv;
    BGv_dif = BGv_dif-(G2-GBv);
    % step (ii) directional color difference *********
    % directional weight
    Bh_g = zeros(size(Bh));
    Bh_g(1:2:end,2:2:end) =Bh(1:2:end,1:2:end); 
    Rh_g = zeros(size(Rh));
    Rh_g(2:2:end,1:2:end) =Rh(2:2:end,2:2:end); 
    
    Bv_g = zeros(size(Bh));
    Bv_g(2:2:end,1:2:end) =Bv(1:2:end,1:2:end); 
    Rv_g = zeros(size(Rv));
    Rv_g(1:2:end,2:2:end) =Rv(2:2:end,2:2:end); 
    
    difh = G + RGh+ BGh -Rh-Bh- Rh_g - Bh_g;%GRh -GBh; % differential horizontal
    difv = G + RGv+ BGv -Rv-Bv- Rv_g -Bv_g; % GRv - GBv; differential vertical
    % color difference gradient
    Kh = [1,0,-1];
    Kv = Kh';
    difh2 = abs(imfilter(difh, Kh, 'replicate'));
    difv2 = abs(imfilter(difv, Kv, 'replicate'));
    % directional weight
    K = ones(5,5);
    wh = imfilter(difh2, K, 'replicate');
    wv = imfilter(difv2, K, 'replicate');
    Kw = [1,0,0,0,0]; 
    Ke = [0,0,0,0,1];
    Ks = Ke'; 
    Kn = Kw';
    Ww = imfilter(wh, Kw, 'replicate');
    We = imfilter(wh, Ke, 'replicate');
    Wn = imfilter(wv, Kn, 'replicate');
    Ws = imfilter(wv, Ks, 'replicate');
    Ww = 1 ./ (Ww.*Ww + 1E-32);
    We = 1 ./ (We.*We + 1E-32);
    Ws = 1 ./ (Ws.*Ws + 1E-32);
    Wn = 1 ./ (Wn.*Wn + 1E-32);
    
    h = [0.56,0.35,0.08,0.01,0]; %fspecial('gaussian', [1,9], 1);
    Ke = [0,0,1,1,1] .* h;  %[0,0,0,0,1,1,1,1,1] .* h; 
    Kw = [1,1,1,0,0] .* h; %[1,1,1,1,1,0,0,0,0] .* h;
    Ke = Ke / sum(Ke, 2);
    Kw = Kw / sum(Kw, 2);
    Ks = Ke'; 
    Kn = Kw';
    difn = imfilter(difv, Kn, 'replicate');
    difs = imfilter(difv, Ks, 'replicate');
    difw = imfilter(difh, Kw, 'replicate');
    dife = imfilter(difh, Ke, 'replicate');
    Wt = Ww + We + Wn + Ws;
    dif = (Wn.*difn + Ws.*difs + Ww.*difw + We.*dife) ./ Wt ;
    
    % Step (iii) ************************************
    kh=[0,0,0,0,1];
    kv=kh';
    WH = imfilter(wh, kh, 'replicate');
    WV = imfilter(wv, kv, 'replicate');
    WV = 1./(WV.*WV + 1E-32);
    WH = 1./(WH.*WH + 1E-32);
    R_hat = ((WH.*Rh)+(WV.*Rv))./(WH+WV);
    B_hat = (((WH.*Bh)+(WV.*Bv))./(WH+WV));
    

    
    mos(:,:,1) = R_hat;
    mos(:,:,2)=rgb(:,:,2);
    mos(:,:,3) = B_hat;

end