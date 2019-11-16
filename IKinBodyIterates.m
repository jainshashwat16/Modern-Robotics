function [thetalist, success] = IKinBodyIterates(Blist, M, T, thetalist0, eomg, ev)
% *** CHAPTER 6: INVERSE KINEMATICS ***
% Takes Blist: The joint screw axes in the end-effector frame when the
%              manipulator is at the home position, in the format of a 
%              matrix with the screw axes as the columns,
%       M: The home configuration of the end-effector,
%       T: The desired end-effector configuration Tsd,
%       thetalist0: An initial guess of joint angles that are close to 
%                   satisfying Tsd,
%       eomg: A small positive tolerance on the end-effector orientation
%             error. The returned joint angles must give an end-effector 
%             orientation error less than eomg,
%       ev: A small positive tolerance on the end-effector linear position 
%           error. The returned joint angles must give an end-effector
%           position error less than ev.
% Returns thetalist: Joint angles that achieve T within the specified 
%                    tolerances,
%         success: A logical value where TRUE means that the function found
%                  a solution and FALSE means that it ran through the set 
%                  number of maximum iterations without finding a solution
%                  within the tolerances eomg and ev.
% Uses an iterative Newton-Raphson root-finding method.
% The maximum number of iterations before the algorithm is terminated has 
% been hardcoded in as a variable called maxiterations. It is set to 20 at 
% the start of the function, but can be changed if needed.  
% Example Inputs:
% 
% clear; clc;
% Blist = [[0; 0; -1; 2; 0; 0], [0; 0; 0; 0; 1; 0], [0; 0; 1; 0; 0; 0.1]];
% M = [[-1, 0, 0, 0]; [0, 1, 0, 6]; [0, 0, -1, 2]; [0, 0, 0, 1]];
% T = [[0, 1, 0, -5]; [1, 0, 0, 4]; [0, 0, -1, 1.6858]; [0, 0, 0, 1]];
% thetalist0 = [1.5; 2.5; 3];
% eomg = 0.01;
% ev = 0.001;
% [thetalist, success] = IKinBody(Blist, M, T, thetalist0, eomg, ev)
% 
% Output:
% thetalist =
%    1.5707
%    2.9997
%    3.1415
% success =
%     1

thetalist = thetalist0;
i = 0;
maxiterations = 20;
Vb = se3ToVec(MatrixLog6(TransInv(FKinBody(M, Blist, thetalist)) * T));
err = norm(Vb(1: 3)) > eomg || norm(Vb(4: 6)) > ev;

% fopen will open the txt file where will be printed all the iteration values 
fileId = fopen('log.txt','w');
    
while err && i < maxiterations
    % Check if it's the first iteration. This is needed for printing also the
    % 0 iteration (with initial guess) to the text file.
    if i ~= 0
        thetalist = thetalist + pinv(JacobianBody(Blist, thetalist)) * Vb;

        Vb = se3ToVec(MatrixLog6(TransInv(FKinBody(M, Blist, thetalist)) * T));
        err = norm(Vb(1: 3)) > eomg || norm(Vb(4: 6)) > ev;
    end
    
    % Save the thetalist at every iteration inside a matrix
    thetaMatrix(i+1,:) = thetalist;
    
    % Calculate the Tsb matrix
    Tsb = FKinBody(M,Blist,thetalist);
    
    % Now it prints the values from the iteration, so the result will be
    % like:
    % 
    %   Iteration 1:
    %
    %   joint vector:
    %   0.500 , 0.512 , 0.523 , 0.534 , 0.545
    %
    %   SE(3) end-effector config:
    %   1   0   0   5
    %   0   1   0   2
    %   0   0   1   0
    %   0   0   0   1
    %
    %   error twist V_b:
    %   0.000 , 0.523 , 0.234 , 0.253 , 0.522 , 0.222
    %
    %   angular error magnitude ||omega_b||: 0.002
    %   linear error magnitude ||v_b||: 0.0002
    %   __________________________________________
    % (the line under the iteration values is to make it more readable
    
    fprintf(fileId,'Iteration %d : \n \n',i);
    fprintf(fileId,'joint vector :\n %.3f , %.3f , %.3f , %.3f , %.3f , %.3f \n\n',thetalist);
    fprintf(fileId,'SE(3) end-effector config: \n');
    fprintf(fileId,'%.3f %.3f %.3f %.3f \n \n',Tsb'); 
    fprintf(fileId,'error twist V_b: \n %.3f , %.3f , %.3f , %.3f , %.3f , %.3f \n\n',Vb);
    fprintf(fileId,'angular error magnitude ||omega_b||: %.4f \n',norm(Vb(1:3)));
    fprintf(fileId,'linear error magnitude ||v_b||: %.4f \n',norm(Vb(4:6)));
    fprintf(fileId,'______________________________________________\n');
    
    i = i + 1;
end

% Close the file
fclose(fileId);

% Create a csv file for store the thetalist data
csvwrite('iterates.csv',thetaMatrix);

success = ~ err;
end