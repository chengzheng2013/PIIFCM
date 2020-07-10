function [up, left, right, down, upleft, upright, downleft, downright] = eigNeig(img)  
% Return the eight neighbour images of the current images 
% Example: [up, left, right, down, upleft, upright, downleft, downright] = eigNeig(img)  

[h, w] = size(img);
% Eight neighbours
up = zeros(h, w);
left = up;
right = up;
down = up;

upleft = up;
upright = up;
downleft = up;
downright = up;

% Copy the corresponding regions
up(2 : end, : ) = img(1 : end-1, : );
down(2 : end - 1, : ) = img(3 : end, : );
left( : , 2 : end) = img( : , 1 : end - 1);
right( : , 2 : end - 1) = img( : , 3 : end);

upleft(2:end, 2:end) = img(1:end-1, 1:end-1);
upright(2 : end, 2 : end - 1) = img(1 : end - 1, 3 : end); 
downleft(2 : end - 1, 2 : end) = img(3 : end, 1 : end - 1);
downright(2 : end - 1, 2 : end - 1) = img(3 : end, 3 : end);

% Extract the central part
upleft = upleft(3 : end - 2, 3 : end - 2);
up = up(3 : end - 2, 3 : end - 2);
upright = upright(3 : end - 2, 3 : end - 2);
left =left(3 : end - 2, 3 : end - 2);
right = right(3 : end - 2, 3 : end - 2);
down = down(3 : end - 2, 3 : end - 2);
downleft = downleft(3 : end - 2, 3 : end - 2);
downright = downright(3 : end - 2, 3 : end - 2);

end