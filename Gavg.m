function gavg = Gavg(x)
    gavg = zeros(size(x, 1), 32, 3);

    for i = 1:size(x, 1)
        for t =1:25
            gavg(i,:, 1) = gavg(i,:, 1) + x(i, :, t);
            gavg(i,:, 2) = gavg(i,:, 2) +  x(i, :, t+25);
            gavg(i,:, 3) = gavg(i,:, 3)+ x(i, :, t+50);
        end
    end
        gavg = gavg(:,:, :)./25;
end