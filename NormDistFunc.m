function y=NormDistFunc(X1,X2)

    % Define points for the first curved line
    x1 = X1(:,1)'; % x-coordinates
    y1 = X1(:,2)'; % y-coordinates

    % Define points for the second curved line
    x2 = X2(:,1); % x-coordinates
    y2 = X2(:,2); % y-coordinates

    % Plot the two curved lines - Uncomment it for debugging
% % %     figure;
% % %     plot(x1, y1, 'b', 'LineWidth', 2);
% % %     hold on;
% % %     plot(x2, y2, 'r', 'LineWidth', 2);
% % %     xlabel('X');
% % %     ylabel('Y');
% % %     title('Curved Lines');

    % Calculate the normal vector at each point of the first curve
    normal_vectors = zeros(2, length(x1)); % Initialize array for normal vectors

    for i = 1:length(x1)
        % Calculate the tangent vector at each point
        if i == 1
            tangent_vector = [x1(i+1) - x1(i), y1(i+1) - y1(i)]; % Forward difference at the first point
        elseif i == length(x1)
            tangent_vector = [x1(i) - x1(i-1), y1(i) - y1(i-1)]; % Backward difference at the last point
        else
            tangent_vector = [x1(i+1) - x1(i-1), y1(i+1) - y1(i-1)] / 2; % Central difference at other points
        end

        % Normalize the tangent vector
        tangent_vector = tangent_vector / norm(tangent_vector);

        % Calculate the normal vector by rotating the tangent vector by 90 degrees counterclockwise
        normal_vector = [0, -1; 1, 0] * tangent_vector';

        % Store the normal vector
        normal_vectors(:, i) = normal_vector;
    end

    % Plot the normal vectors - Uncomment it for debugging
% % %     quiver(x1, y1, normal_vectors(1, :), normal_vectors(2, :), 'g');
% % %     legend('First Curved Line', 'Second Curved Line', 'Normal Vectors');

    % Calculate the distance between the two curved lines in the direction of the first curve's normal vector
    distances_along_normal = zeros(size(x1));
    for i = 1:length(x1)
        % Calculate the vector between corresponding points on the two curves
        vector_between_points = [x2(i) - x1(i), y2(i) - y1(i)];

        % Project the vector onto the normal vector of the first curve
        distance_along_normal = dot(vector_between_points, normal_vectors(:, i)');

        % Store the projected distance
        distances_along_normal(i) = distance_along_normal;
    end
    
    y=abs(distances_along_normal);

    % Plot the distances along the normal direction - Uncomment it for debugging
% % %     figure;
% % %     plot(x1, y, 'k', 'LineWidth', 2);
% % %     xlabel('X');
% % %     ylabel('Distance along Normal');
% % %     title('Distance between Curved Lines along Normal Direction');

end