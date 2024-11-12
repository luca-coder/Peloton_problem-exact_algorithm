import matplotlib.pyplot as plt
import math
import Main
import numpy as np

global m
global n
n = Main.peloton_speed / Main.ralph_speed


def plot_circuit(vertices, ralph_position, points):
    plt.scatter(ralph_position[0], ralph_position[1])
    plt.scatter(5.181817939458401, 9.909090363781402)
    for vertex in vertices:
        plt.scatter(vertex[0], vertex[1], color='green')
    for i in range(0, len(vertices) - 1):
        plt.plot([vertices[i][0], vertices[i + 1][0]], [vertices[i][1], vertices[i + 1][1]], color='blue',
                 linestyle='-', linewidth=2)
    for point in points:
        for element in point:
            try:
                if(element[0] != 6.393939313152802):
                    plt.scatter(element[0], element[1], color='red', label='Interception intervals')
                    plt.xlabel('X')
                    plt.ylabel('Y')
                    plt.title('Optimal solution')
            except Exception:
                print("empty")
    plt.show()


def solve_equation(x_ralph, y_ralph, n):
    solutions = []
    if(y_ralph == 0):
        if(x_ralph <= 0):
            return []
        else:
            solutions.append[((Main.peloton_speed / (Main.peloton_speed + Main.ralph_speed)) * x_ralph), y_ralph]
            solutions.append[((Main.peloton_speed / (Main.peloton_speed - Main.ralph_speed)) * x_ralph), y_ralph]
            return solutions
    else:
        m = x_ralph / y_ralph
        angles = []
        if (m ** 2 - n ** 2 + 1) < 0:
            return angles
        t1 = (m + math.sqrt(m ** 2 - n ** 2 + 1)) / (n + 1)
        t2 = (m - math.sqrt(m ** 2 - n ** 2 + 1)) / (n + 1)
        if (t1 == t2):
            angles.append(2 * math.atan(t1))
            return angles
        angles.append(2 * math.atan(t1))
        angles.append(2 * math.atan(t2))
        return find_extreme_interception_points(x_ralph, y_ralph, angles)


def find_extreme_interception_points(x_ralph, y_ralph, angles):
    points = []

    for angle in angles:
        points.append([x_ralph + y_ralph * (math.cos(angle) / math.sin(angle)), 0])
    return points


def calculate_slope(vertices):
    slopes = []

    for i in range(0, len(vertices) - 1):
        slopes.append((vertices[i + 1][1] - vertices[i][1]) / (vertices[i + 1][0] - vertices[i][0]))
    return slopes


def calculate_lengths(vertices):
    lengths = []

    for i in range(0, len(vertices) - 1):
        lengths.append(math.sqrt((vertices[i][0] - vertices[i + 1][0]) ** 2 +
                                 (vertices[i][1] - vertices[i + 1][1]) ** 2))
    return lengths


def calculate_distance_ralph_endingpoint(ralph_position, index, vertices):

    return math.sqrt((ralph_position[0] - vertices[index][0]) ** 2 +
                     (ralph_position[1] - vertices[index][1]) ** 2)


def resize_following_edge(ralph_position, starting_index, edge_index, vertices, lenghts):
    direction_vector = (vertices[edge_index + 1][0] - vertices[edge_index][0],
                        vertices[edge_index + 1][1] - vertices[edge_index][1])

    magnitude = (direction_vector[0] ** 2 + direction_vector[1] ** 2) ** 0.5

    unit_vector = (direction_vector[0] / magnitude, direction_vector[1] / magnitude)
    length_to_add = calculate_distance_ralph_endingpoint(ralph_position, starting_index + 1, vertices)
    for i in range(starting_index + 1, edge_index):

        length_to_add += lenghts[i]
    new_starting_vertex = (vertices[edge_index][0] - length_to_add * unit_vector[0],
                           vertices[edge_index][1] - length_to_add * unit_vector[1])

    return new_starting_vertex, length_to_add


def resize_previous_edge(ralph_position, starting_index, edge_index, vertices, lenghts):
    direction_vector = (vertices[edge_index + 1][0] - vertices[edge_index][0],
                        vertices[edge_index + 1][1] - vertices[edge_index][1])

    magnitude = (direction_vector[0] ** 2 + direction_vector[1] ** 2) ** 0.5

    unit_vector = (direction_vector[0] / magnitude, direction_vector[1] / magnitude)

    length_to_add = calculate_distance_ralph_endingpoint(ralph_position, starting_index, vertices)
    for i in range(starting_index - 1, edge_index, -1):
        length_to_add += lenghts[i]
    new_starting_vertex = (vertices[edge_index + 1][0] + length_to_add * unit_vector[0],
                           vertices[edge_index + 1][1] + length_to_add * unit_vector[1])
    return new_starting_vertex, length_to_add


def find_ralph_current_interval(slopes, vertices, ralph_position):
    index = 0
    for i in range(0, len(slopes)):
        intercept = vertices[i][1] - slopes[i] * vertices[i][0]
        if (abs(ralph_position[1] - (slopes[i] * ralph_position[0] + intercept)) < 0.0001):
            index = i
            break
    return index


def find_interception_interval_forward(vertices, ralph_position, ralph_index, edge_index, slopes, lenghts):
    new_ralph_position = []
    new_endpoints = []
    new_first_points = []
    interception_points = []
    verified_points = []
    final_points = []
    counter = 0
    distances = []

    new_starting_vertex = resize_following_edge(ralph_position, ralph_index, edge_index, vertices, lenghts)[0]
    distance = resize_following_edge(ralph_position, ralph_index, edge_index, vertices, lenghts)[1]
    new_ralph_position = rotation_and_traslation(slopes[edge_index], new_starting_vertex, vertices[edge_index + 1], ralph_position)[0]
    new_endpoint = rotation_and_traslation(slopes[edge_index], new_starting_vertex, vertices[edge_index + 1], ralph_position)[1]


    result = solve_equation(new_ralph_position[0], new_ralph_position[1], n)

    if(len(result) == 0):
        final_points = []
    if(len(result) == 1):
        print("one point was found")
    if(len(result) > 1):
        verified_points = verify_interception_points(distance, new_endpoint, result[0], result[1])

    if(len(verified_points) > 1):
        final_points = find_real_interception_points(verified_points[1], verified_points[0], new_starting_vertex, vertices[edge_index + 1])

    return final_points


def find_interception_interval_backward(vertices, ralph_position, ralph_index, edge_index, slopes, lenghts):
    new_ralph_position = []
    new_endpoints = []
    new_first_points = []
    interception_points = []
    verified_points = []
    final_points = []
    counter = 0
    distances = []


    new_starting_vertex = resize_previous_edge(ralph_position, ralph_index, edge_index, vertices, lenghts)[0]
    distance = resize_previous_edge(ralph_position, ralph_index, edge_index, vertices, lenghts)[1]
    new_ralph_position = rotation_and_traslation(slopes[edge_index], new_starting_vertex, vertices[edge_index + 1], ralph_position)[0]
    new_endpoint = rotation_and_traslation(slopes[edge_index], new_starting_vertex, vertices[edge_index], ralph_position)[1]


    result = solve_equation(new_ralph_position[0], new_ralph_position[1], n)

    if(len(result) == 0):
        verified_points = []
    if(len(result) == 1):
        print("one point was found")
    if(len(result) > 1):
        verified_points = verify_interception_points(distance, new_endpoint, result[0], result[1])


    if (len(verified_points) > 1):
        final_points = find_real_interception_points(verified_points[1], verified_points[0], new_starting_vertex, vertices[edge_index])


    return final_points


def rotation_and_traslation(slope, starting_point, end_point, ralph_starting_position):


    end_point_translated = [end_point[0] - starting_point[0], end_point[1] - starting_point[1]]


    ralph_position_traslated = [ralph_starting_position[0] - starting_point[0],
                                ralph_starting_position[1] - starting_point[1]]


    if (starting_point[0] < end_point[0] and starting_point[1] < end_point[1]):
        theta = math.atan(slope)
        rotation_matrix = [
            [math.cos(theta), math.sin(theta)],
            [-math.sin(theta), math.cos(theta)]
        ]

    if (starting_point[0] > end_point[0] and starting_point[1] < end_point[1]):
        theta = math.pi - abs(math.atan(slope))

        rotation_matrix = [
            [math.cos(theta), math.sin(theta)],
            [-math.sin(theta), math.cos(theta)]
        ]

    if (starting_point[0] < end_point[0] and starting_point[1] > end_point[1]):
        theta = abs(math.atan(slope))


        rotation_matrix = [
            [math.cos(theta), -math.sin(theta)],
            [math.sin(theta), math.cos(theta)]
        ]

    if (starting_point[0] > end_point[0] and starting_point[1] > end_point[1]):
        theta = (math.pi) - (math.atan(slope))

        rotation_matrix = [
            [math.cos(theta), -math.sin(theta)],
            [math.sin(theta), math.cos(theta)]
        ]

    ralph_position_rotated = [
        rotation_matrix[0][0] * ralph_position_traslated[0] + rotation_matrix[0][1] * ralph_position_traslated[1],
        rotation_matrix[1][0] * ralph_position_traslated[0] + rotation_matrix[1][1] * ralph_position_traslated[1]]

    end_point_rotated = [
        rotation_matrix[0][0] * end_point_translated[0] + rotation_matrix[0][1] * end_point_translated[1], 0]
    return ralph_position_rotated, end_point_rotated


def verify_interception_points(first_point, last_point, first_interception_point, last_interception_point):
    points = []


    if (first_interception_point[0] <= first_point and last_interception_point[0] <= first_point):
        points.append([])
        return points
    if (first_interception_point[0] <= first_point):
        points.append([first_point, 0])
    elif (first_interception_point[0] >= last_point[0]):
        points.append(last_point)
    else:
        points.append(first_interception_point)

    if (last_interception_point[0] <= first_point):
        points.append([first_point, 0])
    elif (last_interception_point[0] >= last_point[0]):
        points.append(last_point)
    else:
        points.append(last_interception_point)
    return points


def find_real_interception_points(first_rotated_point, last_rotated_point, first_real_point, last_real_point):
    first_real_point = np.array(first_real_point)
    last_real_point = np.array(last_real_point)

    points = []
    if (last_real_point[0] - first_real_point[0] > 0):
        direction_vector = last_real_point - first_real_point
        normalized_direction = direction_vector / np.linalg.norm(direction_vector)
        points.append((first_real_point + first_rotated_point[0] * normalized_direction).tolist())
        points.append((first_real_point + last_rotated_point[0] * normalized_direction).tolist())
    if (last_real_point[0] - first_real_point[0] < 0):
        direction_vector = first_real_point - last_real_point
        normalized_direction = direction_vector / np.linalg.norm(direction_vector)
        points.append((first_real_point - first_rotated_point[0] * normalized_direction).tolist())
        points.append((first_real_point - last_rotated_point[0] * normalized_direction).tolist())

    return points


def find_next_interception_interval(extreme_points, vertices, first_edge_index, edge_index, slopes, lengths):
    check_completed = False
    index = 0
    new_extreme_points = []
    new_extreme = []
    new_point = []
    new_forward_interval_temp = []
    new_backward_interval_temp = []
    new_points_list = []
    new_point_forward_temp = []
    extreme_points_list = []

    for extreme_point in extreme_points:

        if(abs(extreme_point[0] - vertices[edge_index][0]) < 0.0000001) or ():
            continue
        new_extreme_points = find_interception_interval_forward(vertices, extreme_point, first_edge_index, edge_index, slopes, lengths)
        extreme_points_list.append(new_extreme_points)
        if(new_extreme_points == []):
            continue
        if(is_comprehended_interval(new_extreme_points, [vertices[edge_index], vertices[edge_index + 1]])):
            return [vertices[edge_index], vertices[edge_index + 1]]
        new_forward_interval_temp = new_extreme_points
        for new_extreme in new_extreme_points:
            new_backward_interval = find_interception_interval_backward(vertices, new_extreme, edge_index, first_edge_index, slopes, lengths)
            new_backward_interval_temp = new_backward_interval
            for point in new_backward_interval:
                if (is_comprehended_point(extreme_points[0][0], extreme_points[1][0], point[0])):
                    new_point = point
                    break
            if(new_point == []):
                return new_point
            while(True):
                new_intermediate_extreme_points_forward = find_interception_interval_forward(vertices, new_point, first_edge_index, edge_index, slopes, lengths)
                if(len(new_intermediate_extreme_points_forward) > 0):
                    if(abs(new_intermediate_extreme_points_forward[0][0] - new_intermediate_extreme_points_forward[1][0]) < 0.1):
                        new_points_list.append(new_intermediate_extreme_points_forward[0])
                        break
                for new_intermediate_extreme_point_forward in new_intermediate_extreme_points_forward:
                    found1 = False
                    for new_forward_point_temp in new_forward_interval_temp:
                        if (abs(new_forward_point_temp[0] - new_intermediate_extreme_point_forward[0]) < 0.1):
                            found1 = True
                            break
                    if(not found1):
                        new_point = new_intermediate_extreme_point_forward
                if(is_comprehended_point(vertices[edge_index][0], vertices[edge_index + 1][0], new_point[0]) == False):
                    final_point1 = calculate_max_extreme(vertices[edge_index], vertices[edge_index + 1], new_point)
                    new_points_list.append(final_point1)
                    break
                new_forward_interval_temp = new_intermediate_extreme_points_forward
                new_intermediate_extreme_points_backward = find_interception_interval_backward(vertices, new_point, edge_index, first_edge_index, slopes, lengths)
                if (len(new_intermediate_extreme_points_backward) > 0):
                    if(abs(new_intermediate_extreme_points_backward[0][0] - new_intermediate_extreme_points_backward[1][0]) < 0.1):
                        new_points_list.append(new_intermediate_extreme_points_backward[0])
                        break
                for new_intermediate_extreme_point_backward in new_intermediate_extreme_points_backward:
                    found2 = False
                    for new_backward_point_temp in new_backward_interval_temp:
                        if (abs(new_backward_point_temp[0] - new_intermediate_extreme_point_backward[0]) < 0.1):
                            found2 = True
                            break
                    if (not found2):
                        new_point = new_intermediate_extreme_point_backward
                if (is_comprehended_point(extreme_points[0][0], extreme_points[1][0], new_point[0]) == False):
                    final_point2 = calculate_max_extreme(extreme_points[0], extreme_points[1], new_point)
                    new_points_list.append(final_point2)
                    break
                new_backward_interval_temp = new_intermediate_extreme_points_backward
    if(new_points_list == []):
        return []
    final_interval = calculate_final_interval(new_points_list, extreme_points_list)
    return final_interval


def calculate_final_interval(new_points_list, original_points_list):
    max_distance = 0
    max_distance_pair = []

    for new_point in new_points_list:
        for original_points in original_points_list:
            for original_point in original_points:
                distance = euclidean_distance(new_point, original_point)
                if distance > max_distance:
                    max_distance = distance
                    max_distance_pair.append([new_point, original_point])

    max_distance_pair = max_distance_pair[len(max_distance_pair) - 1]
    return max_distance_pair


def calculate_max_extreme(point1, point2, value):
    if(euclidean_distance(point1, value) >= euclidean_distance(point2, value)):
        return point1
    else:
        return point2


def euclidean_distance(point1, point2):
    return math.sqrt((point2[0] - point1[0]) ** 2 + (point2[1] - point1[1]) ** 2)


def is_comprehended_point(number1, number2, value):
    min_value = min(number1, number2)
    max_value = max(number1, number2)
    if(min_value <= value <=  max_value):
        return True

    return False


def is_comprehended_point_with_tolerance(number1, number2, value):
    min_value = min(number1, number2)
    max_value = max(number1, number2)
    min_value = [x - 0.00001 for x in min_value]
    max_value = [x + 0.00001 for x in max_value]
    if(min_value <= value <=  max_value):
        return True

    return False


def is_comprehended_interval(interval1, interval2):

    comprehended = False

    for point in interval2:
        if(is_comprehended_point_with_tolerance(interval1[0], interval1[1], point) == False):
            comprehended = False
            break
        comprehended = True
    return comprehended