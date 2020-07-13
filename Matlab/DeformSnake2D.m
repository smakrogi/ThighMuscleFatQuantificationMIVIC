function [xout, yout] = DeformSnake2D(x, y, px, py, options)
%syntax: [xout, yout] = DeformSnake2D(x, y, px, py, options);

% sampling and record number to N
    n=length(x);
    delta=2;
    x = [x;x(1,1)];
    y = [y;y(1,1)];
    t = 1:n+1;
    ts = [1:delta:n+1]';
    xi = interp1(t,x,ts);
    yi = interp1(t,y,ts);
    n = length(xi);
    x = xi(1:n-1);
    y = yi(1:n-1);
    
    if( options.debug )
        figure(1)
    end
    [x,y] = snakeinterp(x,y,3,1);
    if( options.debug )
        snakedisp(x,y,'r')
    end
    
    %     sprintf('\n');
    %     infoString = [];
    %     infoString = sprintf('Scale #%d', itScale);
    %     fprintf('%c',char(8*ones(1,length(infoString))));
    %     fprintf(infoString);
    
    for i=1:options.nIterations,
        % Using traditional evolution.
        %     [x,y] = snakedeform(x,y,alpha,beta,gamma,kappa,px,py,itStep);
        % Using a balloon force.
        [x,y] = snakedeform2(x,y,options.alpha,options.beta, ...
            options.gamma,options.kappa,options.kappap,...
            px,py,options.itStep);
        % "student" version.
%         [x,y] = snakeinterp(x,y,3,1);
        % "professional" version
        [x,y] = snakeinterp(x,y,2,0.5);
        if( options.debug && mod(i,5)==0)
            snakedisp(x,y,'r')
            title(['Deformation in progress,  iter = ' num2str(i*options.itStep)])
            saveas(gcf, ['snake_iteration_', num2str(i*options.itStep)], 'png');
            pause(0.005);
        end
    end
    
    if( options.debug )
        saveas(gcf, 'snake_evolution', 'png')
    end
    
    xout = x;
    yout = y;
    
end