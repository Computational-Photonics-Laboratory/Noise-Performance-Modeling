function output = inner_product(v_left,v_right,dt)
% this innper product is only in computational sense, NOT the physical
% sense. 

output = v_left.'*v_right*dt;
end