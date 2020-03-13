function se_ImageCut(filename,numbers)


A = imread(filename,'png');


for i=1:numel(numbers)
    switch numbers(i)
        case 'R'
            imwrite(A(630:1000,35:500,:),strrep(filename,'.png','_R.png'),'png')
        case 'L'
            imwrite(A(630:1000,[35:500]+525,:),strrep(filename,'.png','_L.png'),'png')
            
        case 'P'
            imwrite(A([615:1020]-480,[35:500]+525,:),strrep(filename,'.png','_P.png'),'png')
        case 'A'
            imwrite(A([615:1020]-480,[70:475]-10,:),strrep(filename,'.png','_A.png'),'png')
            
        case 'B'
            imwrite(flipdim(permute(A([615:990]+480,[55:505]-10,:),[2 1 3]),1),strrep(filename,'.png','_B.png'),'png')
        case 'T'
            imwrite(flipdim(permute(A([615:990]+485,[55:505]+520,:),[2 1 3]),1),strrep(filename,'.png','_T.png'),'png')
            
        case 'X'
            
            top = double(flipdim(permute(A([625:980]+485,[55:505]+520,:),[2 1 3]),1));
            
            r   = A(630:1000,40:495,:);
            l   = A(630:1000,[40:495]+525,:);
            
            Z1 = interp2(top(:,:,1),linspace(1,size(top,2),round(size(r,1)/size(top,1)*size(top,2))),linspace(1,size(top,1),size(r,1))');
            Z1 = cat(3,Z1,interp2(top(:,:,2),linspace(1,size(top,2),round(size(r,1)/size(top,1)*size(top,2))),linspace(1,size(top,1),size(r,1))'));
            Z1 = cat(3,Z1,interp2(top(:,:,3),linspace(1,size(top,2),round(size(r,1)/size(top,1)*size(top,2))),linspace(1,size(top,1),size(r,1))'));
            
            
            img =  cat(2, A(630:1000,[35:500]+525,:), uint8(Z1));
            img =  cat(2, img, A(630:1000,35:500,:));
            
            imwrite(img,strrep(filename,'.png','_LTR.png'),'png')
            
    end
end