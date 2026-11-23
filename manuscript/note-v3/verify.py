import argparse
import itertools
import json
from fractions import Fraction
from math import isqrt
from pathlib import Path


ROWS = tuple(itertools.product((0, 1), repeat=2))
STATES = tuple(itertools.product(ROWS, repeat=2))
OEIS_PREFIX = (2, 3, 6, 10, 21, 42, 86, 179, 370, 770, 1601, 3330,
               6930, 14419, 30006, 62442, 129941, 270410, 562726, 1171043,
               2436962, 5071362, 10553601, 21962242, 45703842, 95110563,
               197926886, 411889610, 857150101, 1783745642, 3712008566,
               7724760339)


def boundary(row, adjacent):
    return all((row[column] == row[1 - column]) +
               (row[column] == adjacent[column]) != 2 for column in range(2))


def interior(previous, row, following):
    return all((row[column] == row[1 - column]) +
               (row[column] == previous[column]) +
               (row[column] == following[column]) != 2 for column in range(2))


def step(vector):
    return {
        (previous, row): sum(vector[row, following] for following in ROWS
                             if interior(previous, row, following))
        for previous, row in STATES
    }


def total(vector):
    return sum(vector[first, second] for first, second in STATES
               if first[0] == 0 and boundary(first, second))


def direct_count(height, width=2):
    if height == 0:
        return 0
    result = 0
    for board in range(0, 1 << (width * height), 2):
        legal = True
        for row in range(height):
            for column in range(width):
                location = width * row + column
                value = (board >> location) & 1
                neighbors = []
                if column > 0:
                    neighbors.append(location - 1)
                if column + 1 < width:
                    neighbors.append(location + 1)
                if row > 0:
                    neighbors.append(location - width)
                if row + 1 < height:
                    neighbors.append(location + width)
                if sum(((board >> neighbor) & 1) == value
                       for neighbor in neighbors) == 2:
                    legal = False
                    break
            if not legal:
                break
        result += legal
    return result


def fibonacci(index):
    previous, current = 0, 1
    for unused in range(index):
        previous, current = current, previous + current
    return previous


def line_count(length):
    return sum(all(not (bits[index] == bits[index + 1] == bits[index + 2])
                   for index in range(length - 2))
               for bits in itertools.product((0, 1), repeat=length) if bits[0] == 0)


def rectangle_count(height, width):
    if height == 1:
        return line_count(width)
    rows = tuple(itertools.product((0, 1), repeat=width))
    states = tuple(itertools.product(rows, repeat=2))

    def horizontal(row, column):
        return sum(row[neighbor] == row[column]
                   for neighbor in (column - 1, column + 1) if 0 <= neighbor < width)

    def end(row, neighbor):
        return all(horizontal(row, column) + (row[column] == neighbor[column]) != 2
                   for column in range(width))

    def middle(previous, row, following):
        return all(horizontal(row, column) + (row[column] == previous[column])
                   + (row[column] == following[column]) != 2 for column in range(width))

    vector = {(first, second): int(end(second, first)) for first, second in states}
    for unused in range(height - 2):
        vector = {(first, second): sum(vector[second, third] for third in rows
                                        if middle(first, second, third))
                  for first, second in states}
    return sum(vector[first, second] for first, second in states
               if first[0] == 0 and end(first, second))


def rectangle_model(width):
    rows = tuple(itertools.product((0, 1), repeat=width))
    states = tuple(itertools.product(rows, repeat=2))
    row_indices = {row: index for index, row in enumerate(rows)}

    def horizontal(row, column):
        return sum(row[column] == row[neighbor] for neighbor in (column - 1, column + 1)
                   if 0 <= neighbor < width)

    def end(row, adjacent):
        return all(horizontal(row, column) + (row[column] == adjacent[column]) != 2
                   for column in range(width))

    transitions = []
    for first, second in states:
        transitions.append([row_indices[second] * len(rows) + row_indices[third]
                            for third in rows
                            if all(horizontal(second, column) + (second[column] == first[column])
                                   + (second[column] == third[column]) != 2
                                   for column in range(width))])
    initial = [int(end(second, first)) for first, second in states]
    functional = [int(first[0] == 0 and end(first, second)) for first, second in states]
    return transitions, initial, functional


def modular_determinant(matrix, prime):
    values = [[entry % prime for entry in row] for row in matrix]
    determinant = 1
    for column in range(len(values)):
        pivot = next((row for row in range(column, len(values)) if values[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            values[pivot], values[column] = values[column], values[pivot]
            determinant = -determinant
        value = values[column][column]
        determinant = determinant * value % prime
        inverse = pow(value, -1, prime)
        for row in range(column + 1, len(values)):
            factor = values[row][column] * inverse % prime
            for target in range(column, len(values)):
                values[row][target] = (values[row][target] - factor * values[column][target]) % prime
    return determinant


def compare_bfile(counts, maximum, path):
    if path is None:
        return None
    records = []
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            index, value = map(int, line.split())
            if 1 <= index <= maximum:
                assert counts[index] == value
                records.append(index)
    assert records == list(range(1, maximum + 1))
    return len(records)


def wider_checks(maximum, width_three_bfile=None, width_four_bfile=None):
    certificates = json.loads(Path(__file__).with_name('width-certificates.json').read_text())
    result = []
    for certificate in certificates['widths']:
        width, degree = certificate['width'], certificate['degree']
        coefficients = certificate['coefficients_in_increasing_degree']
        assert degree == (20 if width == 3 else 48)
        assert len(coefficients) == degree + 1 and coefficients[0] == 1
        transitions, initial, functional = rectangle_model(width)
        vectors = [initial]
        for unused in range(max(maximum, 2 * degree) - 2):
            vectors.append([sum(vectors[-1][target] for target in targets) for targets in transitions])
        counts = [0, fibonacci(width + 1)]
        counts += [sum(weight * value for weight, value in zip(functional, vector)) for vector in vectors]
        assert counts[1:degree + 1] == certificate['initial_n1_to_degree']
        annihilator = [sum(coefficient * vectors[degree - lag][state]
                           for lag, coefficient in enumerate(coefficients))
                       for state in range(len(initial))]
        assert not any(annihilator)
        for height in range(degree + 1, maximum + 1):
            assert sum(coefficient * counts[height - lag]
                       for lag, coefficient in enumerate(coefficients)) == 0
        hankel = [[counts[row + column + 2] for column in range(degree)] for row in range(degree)]
        prime = certificate['hankel_modulus']
        assert prime == 1000000007
        assert all(prime % divisor for divisor in range(2, isqrt(prime) + 1))
        determinant = modular_determinant(hankel, prime)
        assert determinant == certificate['hankel_determinant_mod_prime'] != 0
        published = ([3, 6, 10, 28, 73, 196, 515, 1376, 3686, 9914] if width == 3 else
                     [5, 10, 28, 99, 326, 1080, 3765, 13282, 46928, 166611])
        assert counts[1:11] == published
        bfile = width_three_bfile if width == 3 else width_four_bfile
        result.append({'width': width, 'states': len(initial), 'annihilator_degree': degree,
                       'annihilator_nonzero_coordinates': sum(value != 0 for value in annihilator),
                       'recurrence_comparison_heights': [degree + 1, maximum],
                       'hankel_size': degree, 'hankel_modulus': prime,
                       'hankel_determinant_mod_prime': determinant,
                       'published_prefix_terms': len(published),
                       'external_bfile_matching_terms': compare_bfile(counts, maximum, bfile)})
    return result


def series_coefficients(numerator, denominator, maximum):
    result = []
    for degree in range(maximum + 1):
        coefficient = numerator[degree] if degree < len(numerator) else 0
        coefficient -= sum(denominator[shift] * result[degree - shift]
                           for shift in range(1, min(degree + 1, len(denominator))))
        result.append(coefficient)
    return result


def norm_sibling_checks(norms, maximum):
    matrix_polynomials = [(Fraction(1), Fraction(0)), (Fraction(0), Fraction(1))]
    for unused in range(maximum):
        constant, linear = matrix_polynomials[-1]
        previous_constant, previous_linear = matrix_polynomials[-2]
        matrix_polynomials.append((2 * linear - previous_constant,
                                   2 * constant + linear - previous_linear))
    assert all(matrix_polynomials[index + 1][1] == norms[index]
               for index in range(maximum + 1))
    second_kind = [(Fraction(1), Fraction(0)), (Fraction(2), Fraction(0))]
    for index in range(2, maximum + 1):
        real, imaginary = second_kind[-1]
        previous_real, previous_imaginary = second_kind[-2]
        multiplier = (-imaginary / 2, real / 2) if index % 2 == 0 else (2 * real, 2 * imaginary)
        second_kind.append((multiplier[0] - previous_real, multiplier[1] - previous_imaginary))
    assert all((real * real + imaginary * imaginary) / (4 if index % 2 else 1) == norms[index]
               for index, (real, imaginary) in enumerate(second_kind))
    quartic = (1, -1, -2, -1, 1)
    differences = series_coefficients((0, 1), quartic, maximum + 1)
    assert all(norms[index] == differences[index + 1] - differences[index - 1]
               for index in range(1, maximum + 1))
    gaussian = [(1, 0), (0, 1)]
    for unused in range(2, maximum + 1):
        real = gaussian[-1][0] + gaussian[-2][0]
        imaginary = gaussian[-1][1] + gaussian[-2][1]
        gaussian.append((-imaginary, real))
    denominator = (1, 0, 1, 2, 1)
    real_parts = series_coefficients((1,), denominator, maximum)
    imaginary_parts = series_coefficients((0, 1, 1), denominator, maximum)
    assert gaussian == list(zip(real_parts, imaginary_parts))
    assert all(real * real + imaginary * imaginary == norms[index]
               for index, (real, imaginary) in enumerate(gaussian))
    extended = dict(enumerate(norms))
    for index in range(3, 1 - maximum, -1):
        extended[index - 4] = (extended[index - 1] + 2 * extended[index - 2]
                               + extended[index - 3] - extended[index])
    assert all(extended[-index - 2] == norms[index] for index in range(maximum + 1))
    divisibility_pairs = [(index, multiplier) for index in range(1, 33)
                          for multiplier in range(1, (maximum + 1) // index + 1)]
    assert all(norms[index * multiplier - 1] % norms[index - 1] == 0
               for index, multiplier in divisibility_pairs)
    transformed = series_coefficients((0, 1, 0, -1), quartic, maximum + 1)
    assert transformed == [0] + norms
    printed_transform = series_coefficients((0, 1, 0, -1), (1, -1, 6, -1, 1), maximum)
    assert printed_transform[0] == 0 != norms[0]
    return {'chebyshev_first_kind_and_matrix_indices': [0, maximum],
            'chebyshev_second_kind_indices': [0, maximum],
            'a092886_difference_indices': [1, maximum],
            'a201837_a201838_sum_of_squares_indices': [0, maximum],
            'two_sided_reflection_nonnegative_indices': [0, maximum],
            'shifted_divisibility_pairs': len(divisibility_pairs),
            'shifted_divisibility_first_index_range': [1, 32],
            'corrected_chebyshev_transform_degrees': [0, maximum + 1],
            'source_formula_counterexamples': {
                'printed_T_n_at_n2': {'printed_value': 1, 'c2': norms[2]},
                'printed_matrix_T_n_at_n2': {'printed_value': 1, 'c2': norms[2]},
                'printed_U_n_minus_1_at_n2': {'printed_value': 1, 'c2': norms[2]},
                'printed_transform_constant': {'printed_value': printed_transform[0], 'c0': norms[0]}}}


def check(maximum, bfile=None, width_three_bfile=None, width_four_bfile=None):
    assert 100 <= maximum <= 1000
    vectors = [{(previous, row): int(boundary(row, previous))
                for previous, row in STATES}]
    for index in range(maximum - 2):
        vectors.append(step(vectors[-1]))
    for state in STATES:
        assert vectors[5][state] == (2 * vectors[4][state] + vectors[3][state]
                                     - vectors[2][state] - 2 * vectors[1][state]
                                     + vectors[0][state])
    counts = [0, 2] + [total(vector) for vector in vectors]
    assert len(counts) == maximum + 1
    assert tuple(counts[1:33]) == OEIS_PREFIX
    direct = [direct_count(height) for height in range(9)]
    assert direct == counts[:9]
    for height in range(6, maximum + 1):
        assert counts[height] == (2 * counts[height - 1] + counts[height - 2]
                                   - counts[height - 3] - 2 * counts[height - 4]
                                   + counts[height - 5])
    denominator = (1, -2, -1, 1, 2, -1)
    numerator = (0, 2, -1, -2, -3, 2)
    for degree in range(maximum + 1):
        coefficient = sum(value * counts[degree - shift]
                          for shift, value in enumerate(denominator) if shift <= degree)
        assert coefficient == (numerator[degree] if degree < len(numerator) else 0)
    complex_terms = [(1, 0), (1, 0)]
    for index in range(2, maximum + 1):
        last_real, last_imaginary = complex_terms[-1]
        old_real, old_imaginary = complex_terms[-2]
        complex_terms.append((last_real - old_imaginary, last_imaginary + old_real))
    norms = [real * real + imaginary * imaginary for real, imaginary in complex_terms]
    assert all(counts[index] == 1 + norms[index] for index in range(1, maximum + 1))
    assert all(norms[index] == norms[index - 1] + 2 * norms[index - 2]
               + norms[index - 3] - norms[index - 4] for index in range(4, maximum + 1))
    assert all(line_count(length) == fibonacci(length + 1) for length in range(1, 15))
    rectangle_counts = [[rectangle_count(height, width) for width in range(1, 5)]
                        for height in range(1, 5)]
    for height in range(1, 5):
        for width in range(1, 5):
            assert rectangle_counts[height - 1][width - 1] == direct_count(height, width)
            assert rectangle_counts[height - 1][width - 1] == rectangle_counts[width - 1][height - 1]
    bfile_matches = None
    if bfile:
        records = []
        for line in Path(bfile).read_text().splitlines():
            if line.strip() and not line.lstrip().startswith('#'):
                index, value = map(int, line.split())
                if 1 <= index <= maximum:
                    assert counts[index] == value
                    records.append(index)
        assert records == list(range(1, maximum + 1))
        bfile_matches = len(records)
    return {
        'status': 'pass',
        'state_order': [''.join(map(str, first)) + '/' + ''.join(map(str, second))
                        for first, second in STATES],
        'continuation_vectors_v0_to_v5': [[vector[state] for vector in vectors[:6]]
                                         for state in STATES],
        'finite_certificate_states': 16,
        'literal_board_comparison_heights': [1, 8],
        'literal_counts_including_zero_convention': direct,
        'recurrence_comparison_heights': [6, maximum],
        'generating_function_coefficient_degrees': [0, maximum],
        'quaternion_identity_comparison_heights': [1, maximum],
        'width_one_comparison_heights': [1, 14],
        'general_width_literal_comparison': 'All rectangles with 1 <= height,width <= 4',
        'rectangle_counts_heights_1_to_4_widths_1_to_4': rectangle_counts,
        'published_prefix_terms': len(OEIS_PREFIX),
        'external_bfile_matching_terms': bfile_matches,
        'wider_recurrence_checks': wider_checks(maximum, width_three_bfile, width_four_bfile),
        'norm_sibling_checks': norm_sibling_checks(norms, maximum),
        'counts_1_to_32': counts[1:33],
        'coverage': 'Exact finite integer computations; all-height conclusions are proved in note.tex.'
    }


def main():
    parser = argparse.ArgumentParser(description='Finite checks for the binary ladder formulas.')
    parser.add_argument('--maximum', type=int, default=210)
    parser.add_argument('--bfile', type=Path)
    parser.add_argument('--width-three-bfile', type=Path)
    parser.add_argument('--width-four-bfile', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = json.dumps(check(args.maximum, args.bfile, args.width_three_bfile,
                              args.width_four_bfile), indent=2) + '\n'
    if args.output:
        args.output.write_text(result)
    else:
        print(result, end='')


if __name__ == '__main__':
    main()
