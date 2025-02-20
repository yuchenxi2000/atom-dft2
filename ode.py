import numpy as np


def runge_kutta_4(func, func_mid, h_step, num_steps, y):
    for i in range(num_steps-1):
        k1 = func(i, y[i])
        k2 = func_mid(i, y[i] + 0.5 * h_step * k1)
        k3 = func_mid(i, y[i] + 0.5 * h_step * k2)
        k4 = func(i+1, y[i] + h_step * k3)
        y[i+1] = y[i] + h_step / 6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4)


def adam_4(func, h_step, num_steps, y):
    if num_steps < 4:
        raise Exception(f'4th-order explicit adams for linear 2d ode requires at least 4 points, but {num_steps} were given!')
    for i in range(num_steps-4):
        dy0 = func(i, y[i])
        dy1 = func(i + 1, y[i + 1])
        dy2 = func(i + 2, y[i + 2])
        dy3 = func(i + 3, y[i + 3])
        y4p = y[i + 3] + h_step/24.0 * (55.0 * dy3 - 59.0 * dy2 + 37.0 * dy1 - 9.0 * dy0)
        dy4 = func(i + 4, y4p)
        y[i + 4] = y[i + 3] + h_step/24.0 * (9.0 * dy4 + 19.0 * dy3 - 5.0 * dy2 + dy1)


def adam_4_reversed(func, h_step, num_steps, y):
    if num_steps < 4:
        raise Exception(f'4th-order explicit adams for linear 2d ode requires at least 4 points, but {num_steps} were given!')
    for i in range(num_steps-1, 3, -1):
        dy0 = func(i, y[i])
        dy1 = func(i - 1, y[i - 1])
        dy2 = func(i - 2, y[i - 2])
        dy3 = func(i - 3, y[i - 3])
        y4p = y[i - 3] - h_step / 24.0 * (55.0 * dy3 - 59.0 * dy2 + 37.0 * dy1 - 9.0 * dy0)
        dy4 = func(i - 4, y4p)
        y[i - 4] = y[i - 3] - h_step / 24.0 * (9.0 * dy4 + 19.0 * dy3 - 5.0 * dy2 + dy1)


def implicit_adam_4(coeff_func, h_step, num_steps, y):
    if num_steps < 3:
        raise Exception(f'4th-order implicit adams for linear 2d ode requires at least 3 points, but {num_steps} were given!')
    for i in range(num_steps-3):
        ti0 = i
        ti1 = i + 1
        ti2 = i + 2
        ti3 = i + 3
        mat0, vec0 = coeff_func(ti0)
        mat1, vec1 = coeff_func(ti1)
        mat2, vec2 = coeff_func(ti2)
        mat3, vec3 = coeff_func(ti3)
        z = y[i+2] + h_step/24.0 * (
                19.0 * (mat2 @ y[i+2] + vec2)
                - 5.0 * (mat1 @ y[i+1] + vec1)
                + (mat0 @ y[i] + vec0)
                + 9.0 * vec3
        )
        y[i+3] = np.linalg.solve(np.eye(2) - 9.0/24.0 * h_step * mat3, z)


def implicit_adam_4_reversed(coeff_func, h_step, num_steps, y):
    if num_steps < 3:
        raise Exception(f'4th-order implicit adams for linear 2d ode requires at least 3 points, but {num_steps} were given!')
    for i in range(num_steps-1, 2, -1):
        ti0 = i
        ti1 = i - 1
        ti2 = i - 2
        ti3 = i - 3
        mat0, vec0 = coeff_func(ti0)
        mat1, vec1 = coeff_func(ti1)
        mat2, vec2 = coeff_func(ti2)
        mat3, vec3 = coeff_func(ti3)
        z = y[i-2] - h_step/24.0 * (
                19.0 * (mat2 @ y[i-2] + vec2)
                - 5.0 * (mat1 @ y[i-1] + vec1)
                + (mat0 @ y[i] + vec0)
                + 9.0 * vec3
        )
        y[i-3] = np.linalg.solve(np.eye(2) + 9.0/24.0 * h_step * mat3, z)
