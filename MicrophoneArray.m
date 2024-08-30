% Define coordinates of microphones in cm
M1 = [0, 0, 0];
M2 = [0, 10.2, 0];
M3 = [0, -4.8, -8.66];
M4 = [10, -4.8, 8.66];
M5 = [-10, -4.8, 8.66];
M6 = [-40, 0, 0];

% Combine coordinates into a matrix for easier plotting
mic_positions = [M1; M2; M3; M4; M5; M6];

% Plot microphones M1 to M5 forming the tetrahedron
figure;
hold on;
plot3(mic_positions(1:5,1), mic_positions(1:5,2), mic_positions(1:5,3), 'ks', 'MarkerSize',10,'MarkerFaceColor','k', 'DisplayName', 'Mics 1-5');
text(mic_positions(1:5,1), mic_positions(1:5,2), mic_positions(1:5,3), {'M1','M2','M3','M4','M5'}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Plot M6 with a different color
plot3(M6(1), M6(2), M6(3), 'ro', 'MarkerSize', 10, 'DisplayName', 'Mic M6');
text(M6(1), M6(2), M6(3), 'M6', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Connect the points to form the tetrahedron
tetra_edges = [1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4; 3 5; 4 5];

% for i = 1:size(tetra_edges, 1)
%    plot3(mic_positions(tetra_edges(i,:), 1), mic_positions(tetra_edges(i,:), 2), mic_positions(tetra_edges(i,:), 3), 'b-', 'LineWidth', 1.5);
% end

% Connect M6 with M1 with a different color line
plot3([M1(1), M6(1)], [M1(2), M6(2)], [M1(3), M6(3)], 'r-', 'LineWidth', 1.5);
for i = 2:5
    plot3([M1(1), mic_positions(i,1)], [M1(2), mic_positions(i,2)], [M1(3), mic_positions(i,3)], 'b-', 'LineWidth', 1.5); % Green lines to connect each mic to M1
end
% Set up the plot for better visualization
grid on;
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
title('3D Visualization of Microphones');
legend('show');
axis equal;
view(3);
hold off;
